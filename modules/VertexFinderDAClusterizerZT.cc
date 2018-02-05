/** \class VertexFinderDAClusterizerZT
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  \authors M. Selvaggi, L. Gray
 *
 */


#include "modules/VertexFinderDAClusterizerZT.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMatrixT.h"
#include "TVector3.h"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>

using namespace std;

static const Double_t c_light   = 2.99792458e+8;


//------------------------------------------------------------------------------

VertexFinderDAClusterizerZT::VertexFinderDAClusterizerZT() :
  fVerbose(0),
  fMaxIterations(0),
  fBetaMax(0),
  fBetaStop(0),
  fBetaPurge(0),
  fVertexZSize(0),
  fVertexTSize(0),
  fCoolingFactor(0),
  fDzCutOff(0),
  fD0CutOff(0),
  fDtCutOff(0),
  fD2Merge(0),
{
}

//------------------------------------------------------------------------------

VertexFinderDAClusterizerZT::~VertexFinderDAClusterizerZT()
{
}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::Init()
{
  fVerbose         = GetBool("Verbose", 1);

  // !!FIX defaul values
  fMaxIterations   = GetInt("MaxIterations", 100);

  fBetaMax         = GetDouble("BetaMax ", 0.1);
  fBetaPurge       = GetDouble("BetaPurge", 1.0);
  fBetaStop        = GetDouble("BetaStop", 1.0);

  fVertexZSize = GetDouble("VertexZSize", 0.5); //in mm
  fVertexTSize  = 1E12*GetDouble("VertexTimeSize", 10E-12); //Convert from [s] to [ps]

  fCoolingFactor   = GetDouble("CoolingFactor", 0.8); // Multiply T so to cooldown must be <1

  // !!FIX defaul values
  // DzCutOff also used as initializer to Z_init to do outlayers, as D0CutOff must be O(1) and find anew way to the flowwing line meaning (add new var)
  fDzCutOff        = GetDouble("DzCutOff", 40);  // don't know yet. For the moment 30*DzCutOff is hard cut off for the considered tracks
  fD0CutOff        = GetDouble("D0CutOff", 3);  // d0/sigma_d0, used to compute the pi (weight) of the track
  fDtCutOff        = GetDouble("DtCutOff", 3);  // Used as first line for outlayer rejection O(1). Put only a parameter maybe and use the other for hard cutoff

  fD2Merge         = GetDouble("DzMerge", 2); // (dz/ZSize)^2+(dt/TSize)^2 limit for merging vertices

  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));

  if (fBetaMax < fBetaPurge)
  {
    fBetaPurge = fBetaMax;
    if (fVerbose)
    {
      cout << "BetaPurge set to " << fBetaPurge << endl;
    }
  }

  if (fBetaPurge < fBetaStop)
  {
    fBetaStop = fBetaPurge;
    if (fVerbose)
    {
      cout << "BetaPurge set to " << fBetaPurge << endl;
    }
  }
}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::Process()
{
  fInputArray->Sort();

  TLorentzVector pos, mom;
  if (fVerbose)
  {
     cout<<" start processing vertices ..."<<endl;
     cout<<" Found "<<fInputArray->GetEntriesFast()<<" input tracks"<<endl;
     //loop over input tracks
     fItInputArray->Reset();
     Candidate *candidate;
     while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
     {
        pos = candidate->InitialPosition;
        mom = candidate->Momentum;

        cout<<"pt: "<<mom.Pt()<<", eta: "<<mom.Eta()<<", phi: "<<mom.Phi()<<", z: "<<candidate->DZ/10<<endl;
     }
  }

  // clusterize tracks
  TObjArray *ClusterArray = new TObjArray;
  clusterize(*ClusterArray);

  if (fVerbose){std::cout <<  " clustering returned  "<< ClusterArray->GetEntriesFast() << " clusters  from " << fInputArray->GetEntriesFast() << " selected tracks" <<std::endl;}

  //loop over vertex candidates
  TIterator * ItClusterArray = ClusterArray->MakeIterator();
  ItClusterArray->Reset();
  Int_t ivtx = 0;
  while((candidate = static_cast<Candidate*>(ItClusterArray->Next())))
  {


     if(fVerbose)cout<<"this vertex has: "<<candidate->GetCandidates()->GetEntriesFast()<<" tracks"<<endl;

     /*
     // Somehow fit the vertex from the tracks now
     // loop over tracks belonging to this vertex
     // TIter it1(candidate->GetCandidates());
     // it1.Reset();
     //
     // Candidate *track;
     // while((track = static_cast<Candidate*>(it1.Next())))
     // {
     //    itr++;
     //
     //    double t = track->InitialPosition.T()/c_light;
     //    double dt = track->ErrorT/c_light;
     //    const double time = t;
     //    const double inverr = 1.0/dt;
     //    meantime += time*inverr;
     //    expv_x2  += time*time*inverr;
     //    normw    += inverr;
     //
     //    // compute error position TBC
     //    const double pt = track->Momentum.Pt();
     //    const double z = track->DZ/10.0;
     //    const double err_pt = track->ErrorPT;
     //    const double err_z = track->ErrorDZ;
     //
     //    const double wi = (pt/(err_pt*err_z))*(pt/(err_pt*err_z));
     //    meanpos += z*wi;
     //
     //    meanerr2 += err_z*err_z*wi;
     //    normpos += wi;
     //    sumpt2 += pt*pt;
     //
     // }
     //
     // candidate->Position.SetXYZT(0.0, 0.0, meanpos*10.0 , meantime*c_light);
     // candidate->PositionError.SetXYZT(0.0, 0.0, errpos*10.0 , errtime*c_light);
     // candidate->SumPT2 = sumpt2;
     // candidate->ClusterNDF = itr;
     // candidate->ClusterIndex = ivtx;
     // ivtx++;
     */

     fVertexOutputArray->Add(candidate);



   }// end of cluster loop

  delete ClusterArray;

}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::clusterize(TObjArray &clusters)
{
  tracks_t tks;
  fill(tks);

  unsigned int nt=tks.getSize();
  double rho0=0.0;  // start with no outlier rejection

  vector< Candidate* > clusters;
  if (nt == 0) return clusters;

  vertex_t vtx; // the vertex prototypes
  // initialize:single vertex at infinite temperature
  vtx.addItem(0, 0, 1);

  // Fit the vertex at T=inf and return the starting temperature
  double beta=beta0(tks, vtx);

  // Cool down untill reaching the temperature to finish increasing the number of vertexes
  while(beta < fBetaPurge)
  {
    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0);
      niter++;
    }
    while (delta2 > 1. &&  niter < fMaxIterations)

    beta /= fCoolingFactor;
    if(beta < fBetaStop) split(beta, vtx);
  }

  // Merge vertexes which are too close
  while(merge(vtx, fD2Merge))
  {
    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0);
      niter++;
    }
    while (delta2 > 1. &&  niter < fMaxIterations)
  }

  // Cooldown untill the limit before assigning track to vertices
  while(beta < fBetaMax)
  {
    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0);
      niter++;
    }
    while (delta2 > 1. &&  niter < fMaxIterations)

    beta /= fCoolingFactor;
  }


  // Build the cluster candidates
  for(unsigned int k = 0; k < vtx.getSize(); k++)
  {
    DelphesFactory *factory = GetFactory();
    Candidate * candidate = factory->NewCandidate();

    candidate->ClusterIndex = k;
    candidate->Position.SetXYZT(0.0, 0.0, vtx.z[k] , vtx.t[k]*1E-9*c_light);

    clusters.Add(candidate);
  }


  // Assign each track to the closest vertex
  for(unsigned int i = 0; i< tks.getSize(); i++)
  {
    double d2_min = 0;
    unsigned int k_min;

    for(unsigned int k = 0; k < vtx.getSize(); k++)
    {
      double d2 = Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i])

      if (k == 0 || d2 < d2_min)
      {
        d2_min = d2;
        k_min = k;
      }
    }

    tks.tt[i]->ClusterIndex = k_min;
    ((Candidate *) clusters.At(k_min))->AddCandidate(tks.tt[i]);
  }
}

//------------------------------------------------------------------------------
// Definition of the distance metrci between track and vertex
double VertexFinderDAClusterizerZT::Energy(double t_z, double v_z, double dz2_o, double t_t, double v_t, double dt2_o)
{
  return (t_z - v_z)*(t_z - k_z)* dz2_o + (t_t - v_t)(t_t - k_t)*dt2_o;
}

//------------------------------------------------------------------------------
// Fill tks with the input candidates array
void VertexFinderDAClusterizerZT::fill(tracks_t &tks)
{
  tks.sum_w = 0;
  tks.sum_dz2_o = 0;
  tks.sum_dt2_o = 0;

  Candidate *candidate;

  double z, t, dz_o, dt_o, w;
  double p, e, bz;
  fItInputArray->Reset();
  while(candidate = static_cast<Candidate*>(fItInputArray->Next()))
  {
    z = candidate->DZ; // [mm]
    if(fabs(z) > 3*fDzCutOff) continue;

    //Temporary in v0 where the right mass is assumed
    p = candidate->Momentum.Pt() * sqrt(1 + candidate->CtgTheta*candidate->CtgTheta);
    e = sqrt(p*p + candidate->Mass*candidate->Mass);
    cout << "M:" << candidate->Mass << endl;
    bz = candidate->Momentum.Pt() * candidate->CtgTheta/e;
    t = candidate->Position.T()*1.E9/c_light; // from [mm] to [ps]
    t += (z - candidate->Position.Z())*1E9/(c_light*bz);

    dz2_o = candidate->ErrorDZ*candidate->ErrorDZ;
    dz2_o += fVertexZSize*fVertexZSize;
    // when needed add beam spot width (x-y)?? mha??
    dz2_o = 1/dz2_o; //Multipling is faster than dividing all the times
    tks.sum_dz2_o += dz2_o;

    dt2_o = candidate->ErrorT*1.E9/c_light; // [ps]
    dt2_o *= dt2_o;
    dt2_o += fVertexTSize*fVertexTSize; // [ps^2]
    dt2_o = 1/dt2_o;
    tks.sum_dt2_o += dt2_o;

    if(fD0CutOff > 0 && candidate->ErrorD0 > 0)
    {
      double d0_sig = candidate->D0/candidate->ErrorD0;
      w = exp(d0_sig*d0_sig - fD0CutOff*fD0CutOff);
      w = 1./(1. + w);
    }
    else
    {
      w = 1;
    }

    if(fVerbose)
    {
      cout << "tks add: " << z << "  " << t << "  " << dz2_o << "  " << dt2_o << "  " << &(*candidate) << "  " << w << endl;
    }

    tks.addItem(z, t, dz2_o, dt2_o, &(*candidate), w); //PROVA: rimuovi &(*---)
    tks.sum_w += w;
  }
  tks.NormalizeWeights();

  return tks;
}

//------------------------------------------------------------------------------
// Compute higher phase transition temperature
static double VertexFinderDAClusterizerZT::beta0(tracks_t &tks, vertex_t &vtx)
{
  if(vtx.getSize() != 1)
  {
    throw std::invalid_argument( "Unexpected number of vertices" );
  }

  unsigned int Ntks=tks.getSize();

  //Set vertex position at T=inf as the weighted average of the tracks
  double sum_wz = 0, summ_wt = 0;
  for(unsigned int i = 0; i < Ntks; i++)
  {
    sum_wz += tks.w[i] * tks.z[i] * tks.dz2_o[i];
    sum_wt += tks.w[i] * tks.t[i] * tks.dt2_o[i];
  }
  vtx.z[0] = sum_wz / tks.sum_dz2_o;
  vtx.t[0] = sum_wt / tks.sum_dt2_o;

  // Compute the posterior distribution covariance matrix elements
  double s_zz = 0, s_tt = 0, s_zt = 0;
  for(unsigned int i = 0; i < Ntks; i++)
  {
    double dz = (tks.z[i] - vtx.z[0]) * tks.dz_o[i];
    double dt = (tks.t[i] - vtx.t[0]) * tks.dt_o[i];

    s_zz += tks.w[i] * dz * dz;
    s_tt += tks.w[i] * dt * dt;
    s_zt += tks.w[i] * dt * dz;
  }

  // Copute the max eighenvalue
  double beta_c = (s_tt - s_zz)*(s_tt - s_zz) + 4*s_zt*s_zt;
  beta_c = 1. / (s_tt + s_zz + sqrt(beta_c));

  double out;
  if (beta_c < fBetaMax)
  {
    // Cool down up to a step before the phase transition
    out = int(log(fBetaMax/beta_c) / log(fCoolingFactor)) - 1;
    out = fBetaMax / pow(fCoolingFactor, out);
  }
  else
  {
    out = fBetaMax * fCoolingFactor;
  }

  if(fVerbose)
  {
    cout << "Beta0: " << beta_c << " " << out << endl;
  }

  return out;
}

//------------------------------------------------------------------------------
// Compute the new vertexes position and mass (probability) -- mass constrained annealing without noise
// Compute and store the posterior covariance matrix elements
// Returns the squared sum of changes of vertexex position normalized by the vertex size declared in the init
double VertexFinderDAClusterizerZT::update(double beta, tracks_t &tks, vertex_t &vtx, double rho0)
{
  unsigned int nt = tks.getSize();
  unsigned int nv = vtx.getSize();

  //initialize sums
  double Z_init = rho0 * exp(-beta * fDzCutOff * fDzCutOff); // Add fDtCutOff here toghether  with this

  // Compute all the energies (aka distances) and normalization partition function
  vector<double> exp_betaE(nt * nv);
  for (unsigned int k = 0; k < nv; k++)
  {
    for (unsigned int i = 0; i < nt; i++)
    {
      if(k == 0) tks.Z[i] = Z_init;

      double aux = Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i]);
      aux = exp(-beta * aux);
      tks.Z[i] += aux;

      unsigned int idx = k*nt + i;
      exp_betaE1[idx] = aux;
    }
  }

  double delta2 = 0; // Sum of vertex displacement in this run
  for (unsigned int k = 0; k < nv; k++)
  {
    // Compute the new vertex positions and masses
    double pk_new = 0;
    double sw_z = 0, sw_t = 0;
    // Compute the posterior covariance matrix Elements
    double szz = 0, stt = 0, stz = 0;

    for (unsigned int i = 0; i < nt; i++)
    {
      unsigned int idx = k*nt + i;

      double p_ygx = vtx.pk[k] * exp_betaE[idx] / tks.Z[i]      //p(y|x), Gibbs distribution
      pk_new += tks.w[i] * p_ygx;

      sw_z += tks.w[i] * p_ygx * tks.z[i] * tks.dz2_o[i];
      sw_t += tks.w[i] * p_ygx * tks.t[i] * tks.dt2_o[i];

      // Add the track contribution to the covariance matrix
      double p_xgy = p_ygx * tks.w[i] / vtx.pk[k];
      double dz = (tks.z[i] - vtx.z[0]) * tks.dz_o[i];
      double dt = (tks.t[i] - vtx.t[0]) * tks.dt_o[i];

      szz += tks.w[i] * dz * dz;
      stt += tks.w[i] * dt * dt;
      szt += tks.w[i] * dt * dz;
    }

    double new_z = sw_z/tks.sum_dz2_o;
    double new_t = sw_t/tks.sum_dt2_o;

    double z_displ = (new_z - vtx.z[k])/fVertexZSize;
    double t_displ = (new_t - vtx.t[k])/fVertexTSize;
    delta2 += z_displ*z_displ + t_displ*t_displ;

    vtx.z[k] = new_z;
    vtx.t[k] = new_t;
    vtx.pk[k] = pk_new;
    vtx.szz[k] = szz;
    vtx.stt[k] = stt;
    vtx.szt[k] = szt;
  }

  return delta2;
}

//------------------------------------------------------------------------------
// Split critical vertices (beta_c < beta)
// Returns true if at least one cluster was split
static bool VertexFinderDAClusterizerZT::split(double beta, vertex_t &vtx, const double epsilon = 2.)
{
  bool split = false;

  auto pair_bc_k = vtx.ComputeAllBeta_c();

  // If minimum beta_c is higher than beta, no split is necessaire
  if( pair_bc_k.first > beta )
  {
    split = false;
  }
  else
  {
    for(unsigned int k = 0; k < vtx.getSize(), k++)
    {
      if(vtx.beta_c[k] <= beta)
      {
        split = true;
        // Compute splitting direction: given by the max eighenvalue eighenvector
        double z_old = vtx.z[k];
        double t_old = vtx.t[k];
        double pk_old = vtx.pk[k];

        double aux_slope = (vxt.szz[k] - vxt.stt[k])*(vxt.szz[k] - vxt.stt[k]) + 4*vtx.stz[k];
        aux_slope = vtx.szz[k] - vtx.stt[k] - sqrt(aux_slope);
        aux_slope = -2*vtx.stz[k] / aux_slope;

        vtx.t[k] = t_old + epsilon*fVertexTSize;
        vtx.z[k] = z_old + aux_slope * epsilon * fVertexZSize;
        vtx.pk[k] = pk_old/2.;

        new_t = t_old - epsilon*fVertexTSize;
        new_z = z_old - aux_slope * epsilon * fVertexZSize;
        new_pk = pk_old/2.;

        vtx.addItem(new_z, new_t, new_pk);
      }
    }
  }

  return split;
}


//------------------------------------------------------------------------------
// Merge vertexes closer than declared dimensions
static bool VertexFinderDAClusterizerZT::merge(vertex_t & vtx, double d2_merge = 2)
{
  bool merged = false;
  unsigned int k1 = 0, k2 =0;

  for(unsigned int k1 = 0; k1 < vtx.getSize(); k1++)
  {
    k2 = 0;
    for(unsigned int k2 = k1+1; k2 < vtx.getSize();)
    {
      if(vtx.DistanceSquare(k1, k2) > d2_merge)
      {
        k2++;
      }
      else
      {
        double new_pk = pk[k1] + pk[k2];
        double new_z = (z[k1]*pk[k1] + z[k2]*pk[k2])/new_pk;
        double new_t = (t[k1]*pk[k1] + t[k2]*pk[k2])/new_pk;

        vtx.removeItem(k2);
        vtx.z[k1] = new_z;
        vtx.t[k1] = new_t;
        vtx.pk[k1] = new_pk;

        merged = true;
      }
    }
  }

  return merged;
}
