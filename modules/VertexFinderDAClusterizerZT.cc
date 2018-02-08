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

#include "TAxis.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLegend.h"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>

using namespace std;

static const Double_t c_light   = 2.99792458e+8;


//------------------------------------------------------------------------------

VertexFinderDAClusterizerZT::VertexFinderDAClusterizerZT()
{
  fVerbose = 0;
  fMaxIterations = 0;
  fBetaMax = 0;
  fBetaStop = 0;
  fBetaPurge = 0;
  fVertexZSize = 0;
  fVertexTSize = 0;
  fCoolingFactor = 0;
  fDzCutOff = 0;
  fD0CutOff = 0;
  fDtCutOff = 0;
  fD2Merge = 0;
  fSplittingSize = 0;
}

//------------------------------------------------------------------------------

VertexFinderDAClusterizerZT::~VertexFinderDAClusterizerZT()
{
}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::Init()
{
  fVerbose         = GetInt("Verbose", 0);

  // !!FIX defaul values
  fMaxIterations   = GetInt("MaxIterations", 100);
  fMaxVertexNumber = GetInt("MaxVertexNumber", 500);

  fBetaMax         = GetDouble("BetaMax", 1.);
  fBetaPurge       = GetDouble("BetaPurge", 0.5);
  fBetaStop        = GetDouble("BetaStop", 0.2);

  fVertexZSize     = GetDouble("VertexZSize", 0.05); //in mm
  fVertexTSize     = 1E12*GetDouble("VertexTimeSize", 10E-12); //Convert from [s] to [ps]

  fCoolingFactor   = GetDouble("CoolingFactor", 0.8); // Multiply T so to cooldown must be <1

  // !!FIX defaul values
  // DzCutOff also used as initializer to Z_init to do outlayers, as D0CutOff must be O(1) and find anew way to the flowwing line meaning (add new var)
  fDzCutOff        = GetDouble("DzCutOff", 40);  // For the moment 3*DzCutOff is hard cut off for the considered tracks
  fD0CutOff        = GetDouble("D0CutOff", 1);  // d0/sigma_d0, used to compute the pi (weight) of the track
  fDtCutOff        = GetDouble("DtCutOff", 160);  // [ps], 3*DtCutOff is hard cut off for tracks

  fD2UpdateLim     = GetDouble("D2UpdateLim", 1.); // ((dz/ZSize)^2+(dt/TSize)^2)/nv limit for merging vertices
  fD2Merge         = GetDouble("D2Merge", 2.0); // (dz/ZSize)^2+(dt/TSize)^2 limit for merging vertices
  fSplittingSize   = GetDouble("SplittingSize", 0.5); // Size of the perturbation when splitting


  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));

  if(fVerbose)
  {
    cout << setprecision(2) << std::scientific;
  }

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

  if (fVerbose)
  {
     cout<< endl << "      Start processing vertices with VertexFinderDAClusterizerZT" << endl;
     cout<<" Found "<<fInputArray->GetEntriesFast()<<" input tracks"<<endl;
  }

  // clusterize tracks
  TObjArray *ClusterArray = new TObjArray;
  clusterize(*ClusterArray);

  if(fVerbose>10)
  {
    unsigned int N = fEnergy_rec.size();
    TGraph* gr1 = new TGraph(N, &fBeta_rec[0], &fNvtx_rec[0]);
    gr1->SetName("gr1");
    gr1->GetXaxis()->SetTitle("beta");
    gr1->GetYaxis()->SetTitle("# Vtx");
    TGraph* gr2 = new TGraph(N, &fBeta_rec[0], &fEnergy_rec[0]);
    gr2->SetName("gr2");
    gr2->GetXaxis()->SetTitle("beta");
    gr2->GetYaxis()->SetTitle("Total Energy");
    TGraph* gr3 = new TGraph(N, &fNvtx_rec[0], &fEnergy_rec[0]);
    gr3->SetName("gr3");
    gr3->GetXaxis()->SetTitle("# Vtx");
    gr3->GetYaxis()->SetTitle("Total Energy");

    auto f = new TFile("~/Desktop/debug/EnergyStat.root", "recreate");
    gr1->Write("gr1");
    gr2->Write("gr2");
    gr3->Write("gr3");

    f->Close();
  }

  if (fVerbose){std::cout <<  " clustering returned  "<< ClusterArray->GetEntriesFast() << " clusters  from " << fInputArray->GetEntriesFast() << " selected tracks" <<std::endl;}

  //loop over vertex candidates
  TIterator * ItClusterArray = ClusterArray->MakeIterator();
  ItClusterArray->Reset();
  Candidate *candidate;
  unsigned int k = 0;
  while((candidate = static_cast<Candidate*>(ItClusterArray->Next())))
  {
     if(fVerbose)
     {
       cout << Form("Cluster %d has %d tracks ", k, candidate->GetCandidates()->GetEntriesFast()) << endl;
     }

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
     k++;
   }// end of cluster loop

  delete ClusterArray;

}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::clusterize(TObjArray &clusters)
{
  tracks_t tks;
  fill(tks);
  if(fVerbose)
  {
    cout << "Tracks added: " << tks.getSize() << endl;
  }

  unsigned int nt=tks.getSize();
  double rho0=0.0;  // start with no outlier rejection

  if (nt == 0) return;

  vertex_t vtx; // the vertex prototypes
  vtx.ZSize = fVertexZSize;
  vtx.TSize = fVertexTSize;
  // initialize:single vertex at infinite temperature
  vtx.addItem(0, 0, 1);

  // Fit the vertex at T=inf and return the starting temperature
  double beta=beta0(tks, vtx);

  if( fVerbose > 1 )
  {
    cout << "Cluster position at T=inf: z = " << vtx.z[0] << " mm , t = " << vtx.t[0] << " ps" << "  pk = " << vtx.pk[0] << endl;
    cout << "Beta Start = " << setprecision(6) << beta << endl;
  }

  if( fVerbose > 10 ) plot_status(beta, vtx, tks, 0, "Ast");

  if( fVerbose > 2)
  {
    cout << "Cool down untill reaching the temperature to finish increasing the number of vertexes" << endl;
  }
  while(beta < fBetaPurge)
  {

    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0);

      if( fVerbose > 10 ) plot_status(beta, vtx, tks, niter, "Bup");
      if (fVerbose > 3)
      {
        cout << niter << ": " << delta2 << endl;
      }
      niter++;
    }
    while (delta2 > fD2UpdateLim &&  niter < fMaxIterations);


    unsigned int n_it = 0;
    while(merge(vtx, fD2Merge) && n_it < fMaxIterations)
    {
      unsigned int niter=0;
      double delta2 = 0;
      do  {
        delta2 = update(beta, tks, vtx, rho0);
        niter++;
      }
      while (delta2 > fD2UpdateLim &&  niter < fMaxIterations);
      n_it++;

      if( fVerbose > 10 ) plot_status(beta, vtx, tks, n_it, "Cme");
    }

    beta /= fCoolingFactor;

    if(beta < fBetaStop && vtx.getSize() < min(fMaxVertexNumber, tks.getSize()))
    {
      split(beta, vtx, fSplittingSize);
      if( fVerbose > 10 ) plot_status(beta, vtx, tks, 0, "Asp");
    }

    if(fVerbose > 3)
    {
      cout << endl << endl << " ----- Beta = " << beta << " --------" << endl;
      cout << "Nv: " << vtx.getSize() << endl;
    }
  }

  if(fVerbose > 2)
  {
    cout << "Cooldown untill the limit before assigning track to vertices" << endl;
  }
  while(beta < fBetaMax)
  {
    unsigned int niter=0;
    double delta2 = 0;
    do  {
      delta2 = update(beta, tks, vtx, rho0);
      niter++;
      if( fVerbose > 10 ) plot_status(beta, vtx, tks, 0, "Bup");
    }
    while (delta2 > 0.5*fD2UpdateLim &&  niter < fMaxIterations);

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
      double d2 = Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i]);

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
  return (t_z - v_z)*(t_z - v_z)* dz2_o + (t_t - v_t)*(t_t - v_t)*dt2_o;
}

//------------------------------------------------------------------------------
// Fill tks with the input candidates array
void VertexFinderDAClusterizerZT::fill(tracks_t &tks)
{
  tks.sum_w_o_dt2 = 0;
  tks.sum_w_o_dz2 = 0;
  tks.sum_w = 0;

  Candidate *candidate;

  double z, t, dz2_o, dt2_o, w;
  double p, e, bz;
  fItInputArray->Reset();

  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    z = candidate->DZ; // [mm]
    if(fabs(z) > 3*fDzCutOff) continue;

    //Temporary in v0 where the right mass is assumed
    p = candidate->Momentum.Pt() * sqrt(1 + candidate->CtgTheta*candidate->CtgTheta);
    e = sqrt(p*p + candidate->Mass*candidate->Mass);
    bz = candidate->Momentum.Pt() * candidate->CtgTheta/e;
    t = candidate->Position.T()*1.E9/c_light; // from [mm] to [ps]
    t += (z - candidate->Position.Z())*1E9/(c_light*bz);
    if(fabs(t) > 3*fDtCutOff) continue;

    dz2_o = candidate->ErrorDZ*candidate->ErrorDZ;
    dz2_o += fVertexZSize*fVertexZSize;
    // when needed add beam spot width (x-y)?? mha?
    dz2_o = 1/dz2_o; //Multipling is faster than dividing all the times

    dt2_o = candidate->ErrorT*1.E9/c_light; // [ps]
    dt2_o *= dt2_o;
    dt2_o += fVertexTSize*fVertexTSize; // [ps^2]
    // Ideally we should also add the induced uncertantiy from dz and z_out. For the moment we compensae using a high value for vertex time.
    dt2_o = 1/dt2_o;

    if(fD0CutOff > 0 && candidate->ErrorD0 > 0)
    {
      double d0_sig = candidate->D0/candidate->ErrorD0;
      w = exp(d0_sig*d0_sig - fD0CutOff*fD0CutOff);
      w = 1./(1. + w);
      if (w < 1E-10) continue;
    }
    else
    {
      w = 1;
    }

    tks.sum_w_o_dt2 += w * dt2_o;
    tks.sum_w_o_dz2 += w * dz2_o;
    tks.sum_w += w;

    tks.addItem(z, t, dz2_o, dt2_o, &(*candidate), w, candidate->PID); //PROVA: rimuovi &(*---)
  }

  if(fVerbose > 1)
  {
    cout << "----->Filled tracks" << endl;
    cout << "M\t\tz\t\tdz\t\tt\t\tdt\t\tw" << endl;
    for(unsigned int i = 0; i < tks.getSize(); i++)
    {
      cout << tks.PID[i] << "\t\t" << tks.z[i] << "\t\t" << 1/sqrt(tks.dz2_o[i]) << "\t\t" << tks.t[i] << "\t\t" << 1/sqrt(tks.dt2_o[i]) << "\t\t" << tks.w[i] << endl;
    }
  }

  return;
}

//------------------------------------------------------------------------------
// Compute higher phase transition temperature
double VertexFinderDAClusterizerZT::beta0(tracks_t & tks, vertex_t &vtx)
{
  if(vtx.getSize() != 1)
  {
    throw std::invalid_argument( "Unexpected number of vertices" );
  }

  unsigned int nt = tks.getSize();

  //Set vertex position at T=inf as the weighted average of the tracks
  double sum_wz = 0, sum_wt = 0;
  for(unsigned int i = 0; i < nt; i++)
  {
    sum_wz += tks.w[i] * tks.z[i] * tks.dz2_o[i];
    sum_wt += tks.w[i] * tks.t[i] * tks.dt2_o[i];
  }
  vtx.t[0] = sum_wt / tks.sum_w_o_dt2;
  vtx.z[0] = sum_wz / tks.sum_w_o_dz2;

  // Compute the posterior distribution covariance matrix elements
  double s_zz = 0, s_tt = 0, s_tz = 0;
  for(unsigned int i = 0; i < nt; i++)
  {
    double dz = (tks.z[i] - vtx.z[0]) * tks.dz_o[i];
    double dt = (tks.t[i] - vtx.t[0]) * tks.dt_o[i];

    s_zz += tks.w[i] * dz * dz;
    s_tt += tks.w[i] * dt * dt;
    s_tz += tks.w[i] * dt * dz;
  }
  s_tt /= tks.sum_w;
  s_zz /= tks.sum_w;
  s_tz /= tks.sum_w;

  // Copute the max eighenvalue
  double beta_c = (s_tt - s_zz)*(s_tt - s_zz) + 4*s_tz*s_tz;
  beta_c = 1. / (s_tt + s_zz + sqrt(beta_c));

  double out;
  if (beta_c < fBetaMax)
  {
    // Cool down up to a step before the phase transition
    out = beta_c * sqrt(fCoolingFactor);
  }
  else
  {
    out = fBetaMax * fCoolingFactor;
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
  vector<double> exp_betaE = Compute_exp_betaE(beta, vtx, tks, Z_init);

  double sum_pk = 0;
  double delta2 = 0; // Sum of vertex displacement in this run
  for (unsigned int k = 0; k < nv; k++)
  {
    // Compute the new vertex positions and masses
    double pk_new = 0;
    double sw_z = 0, sw_t = 0;
    // Compute the posterior covariance matrix Elements
    double szz = 0, stt = 0, stz = 0;
    // double sum_w = 0;
    double sum_wt = 0, sum_wz = 0;

    for (unsigned int i = 0; i < nt; i++)
    {
      unsigned int idx = k*nt + i;

      if(exp_betaE[idx] == 0) continue;

      double p_ygx = vtx.pk[k] * (exp_betaE[idx] / tks.Z[i]);      //p(y|x), Gibbs distribution
      if(isnan(p_ygx))
      {
        cout << Form("%1.2e    %1.2e    %1.2e", vtx.pk[k], exp_betaE[idx], tks.Z[i]);
        throw std::invalid_argument("p_ygx is nan");
      }
      // sum_w += tks.w[i];
      pk_new += tks.w[i] * p_ygx;

      double wt = tks.w[i] * p_ygx * tks.dt2_o[i];
      sw_t += wt * tks.t[i];
      sum_wt += wt;

      double wz = tks.w[i] * p_ygx * tks.dz2_o[i];
      sw_z += wz * tks.z[i];
      sum_wz += wz;

      // Add the track contribution to the covariance matrix
      double p_xgy = p_ygx * tks.w[i] / vtx.pk[k];
      double dz = (tks.z[i] - vtx.z[k]) * tks.dz_o[i];
      double dt = (tks.t[i] - vtx.t[k]) * tks.dt_o[i];

      //DEBUG
      // cout << Form("{%d,%d}: p_ygx=%1.2f, dz=%1.2f, sig_z=%1.2f", k, i, p_ygx, 1./tks.dz_o[i], dz) << endl;
      // cout << Form("Beta = %f, Energy = %f", beta, Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i])) << endl;

      szz += p_xgy * dz * dz;
      stt += p_xgy * dt * dt;
      stz += p_xgy * dt * dz;
    }
    pk_new /= tks.sum_w;
    sum_pk += pk_new;

    stt /= tks.sum_w;
    szz /= tks.sum_w;
    stz /= tks.sum_w;

    double new_t = sw_t/sum_wt;
    double new_z = sw_z/sum_wz;
    if(isnan(new_z))
    {
      cout << "pk " << k << "  " << vtx.pk[k] << endl;
      throw std::invalid_argument("new_z is nan");
    }

    double z_displ = (new_z - vtx.z[k])/fVertexZSize;
    double t_displ = (new_t - vtx.t[k])/fVertexTSize;
    delta2 += z_displ*z_displ + t_displ*t_displ;

    vtx.z[k] = new_z;
    vtx.t[k] = new_t;
    vtx.pk[k] = pk_new;
    vtx.szz[k] = szz;
    vtx.stt[k] = stt;
    vtx.stz[k] = stz;
  }

  if(fabs((sum_pk - 1.) > 1E-6))
  {
    cout << "sum_pk " << sum_pk << endl;
    for (unsigned int k = 0; k < nv; k++)
    {
      cout << Form("%d: %1.4e", k, vtx.pk[k]) << endl;
    }
    throw std::invalid_argument("Sum of masses not unitary");
  }
  if(fVerbose > 3)
  {
    cout << "===Update over" << endl;
    for (unsigned int k = 0; k < nv; k++)
    {
      cout << k << endl;
      cout << "z: " << vtx.z[k] << " , t: " << vtx.t[k] << " , p: " << vtx.pk[k] << endl;
      cout << " | " << vtx.szz[k] << "   " << vtx.stz[k] << "|" << endl;
      cout << " | " << vtx.stz[k] << "   " << vtx.stt[k] << "|" << endl << endl;
    }
    cout << "=======" << endl;
  }

  return delta2/nv;
}

//------------------------------------------------------------------------------
// Split critical vertices (beta_c < beta)
// Returns true if at least one cluster was split
bool VertexFinderDAClusterizerZT::split(double &beta, vertex_t &vtx, double epsilon)
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
    unsigned int nv = vtx.getSize();
    for(unsigned int k = 0; k < nv; k++)
    {
      if( fVerbose > 3 )
      {
        cout << "vtx " << k << "  beta_c = " << vtx.beta_c[k] << endl;
      }
      if(vtx.beta_c[k] <= beta)
      // if(k == pair_bc_k.second)
      {
        split = true;
        // Compute splitting direction: given by the max eighenvalue eighenvector
        double z_old = vtx.z[k];
        double t_old = vtx.t[k];
        double pk_old = vtx.pk[k];

        double aux_slope = (vtx.szz[k] - vtx.stt[k])*(vtx.szz[k] - vtx.stt[k]) + 4*vtx.stz[k]*vtx.stz[k];
        aux_slope = vtx.szz[k] - vtx.stt[k] - sqrt(aux_slope);
        aux_slope = -2*vtx.stz[k] / aux_slope;
        // Move the vertex of epsilon
        epsilon /= TMath::Max(1., double(floor(beta)));
        double delta = 0.5*(vtx.szz[k] + vtx.stt[k])*epsilon/sqrt(1+aux_slope*aux_slope);

        if( fVerbose > 3 )
        {
          cout << aux_slope << "   delta = " << delta << endl;
        }

        double t_displ = delta * fVertexTSize;
        double z_displ = aux_slope * delta * fVertexZSize;

        vtx.t[k] = t_old + t_displ;
        vtx.z[k] = z_old + z_displ;
        vtx.pk[k] = pk_old/2.;

        double new_t = t_old - t_displ;
        double new_z = z_old - z_displ;
        double new_pk = pk_old/2.;

        vtx.addItem(new_z, new_t, new_pk);

        if( fVerbose > 3 )
        {
          cout << "===Split happened on vtx " << k << endl;
          cout << "OLD     z: " << z_old << " , t: " << t_old << " , pk: " << pk_old << endl;
          cout << "NEW+    z: " << vtx.z[k] << " , t: " << vtx.t[k] << " , pk: " << vtx.pk[k] << endl;
          cout << "NEW-    z: " << new_z << " , t: " << new_t << " , pk: " << new_pk <<  endl;
        }
      }
    }
  }
  return split;
}


//------------------------------------------------------------------------------
// Merge vertexes closer than declared dimensions
bool VertexFinderDAClusterizerZT::merge(vertex_t & vtx, double d2_merge = 2)
{
  bool merged = false;

  if(vtx.getSize() < 2) return merged;

  for(unsigned int k1 = 0; k1 < vtx.getSize(); k1++)
  {
    for(unsigned int k2 = k1+1; k2 < vtx.getSize();)
    {
      if(vtx.DistanceSquare(k1, k2) > d2_merge)
      {
        k2++;
      }
      else
      {
        double new_pk = vtx.pk[k1] + vtx.pk[k2];
        double new_z = (vtx.z[k1]*vtx.pk[k1] + vtx.z[k2]*vtx.pk[k2])/new_pk;
        double new_t = (vtx.t[k1]*vtx.pk[k1] + vtx.t[k2]*vtx.pk[k2])/new_pk;

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

// -----------------------------------------------------------------------------
// Compute all the energies and set the partition function normalization for each track
vector<double> VertexFinderDAClusterizerZT::Compute_exp_betaE(double beta, vertex_t &vtx, tracks_t &tks, double Z_init)
{
  unsigned int nt = tks.getSize();
  unsigned int nv = vtx.getSize();

  vector<double> exp_betaE(nt * nv);
  for (unsigned int k = 0; k < nv; k++)
  {
    for (unsigned int i = 0; i < nt; i++)
    {
      if(k == 0) tks.Z[i] = Z_init;

      double aux = Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i]);
      aux = exp(-beta * aux);
      // if(aux < 1E-10) continue;
      tks.Z[i] += vtx.pk[k] * aux;

      unsigned int idx = k*nt + i;
      exp_betaE[idx] = aux;
    }
  }
  return exp_betaE;
}


// -----------------------------------------------------------------------------
// Plot status
void VertexFinderDAClusterizerZT::plot_status(double beta, vertex_t &vtx, tracks_t &tks, int n_it, const char* flag)
{
  vector<int> vtx_color = {2,4,8,1,5,6,9,14,46,3};
  while(vtx.getSize() > vtx_color.size()) vtx_color.push_back(40);

  vector<double> t_PV, dt_PV, z_PV, dz_PV;
  vector<double> t_PU, dt_PU, z_PU, dz_PU;

  double ETot = 0;
  vector<double> exp_betaE = Compute_exp_betaE(beta, vtx, tks, 0);

  for(unsigned int i = 0; i < tks.getSize(); i++)
  {
    for(unsigned int k = 0; k < vtx.getSize(); k++)
    {
      unsigned int idx = k*tks.getSize() + i;
      if(exp_betaE[idx] == 0) continue;

      double p_ygx = vtx.pk[k] * (exp_betaE[idx] / tks.Z[i]);

      ETot += tks.w[i] * p_ygx * Energy(tks.z[i], vtx.z[k], tks.dz2_o[i], tks.t[i], vtx.t[k], tks.dt2_o[i]);
    }

    if(tks.tt[i]->IsPU)
    {
      t_PU.push_back(tks.t[i]);
      dt_PU.push_back(1./tks.dt_o[i]);
      z_PU.push_back(tks.z[i]);
      dz_PU.push_back(1./tks.dz_o[i]);
    }
    else
    {
      t_PV.push_back(tks.t[i]);
      dt_PV.push_back(1./tks.dt_o[i]);
      z_PV.push_back(tks.z[i]);
      dz_PV.push_back(1./tks.dz_o[i]);
    }
  }


  ETot /= tks.sum_w;
  fEnergy_rec.push_back(ETot);
  fBeta_rec.push_back(beta);
  fNvtx_rec.push_back(vtx.getSize());

  double t_min = TMath::Min(  TMath::MinElement(t_PV.size(), &t_PV[0]), TMath::MinElement(t_PU.size(), &t_PU[0])  );
  t_min = TMath::Min(t_min, TMath::MinElement(vtx.getSize(), &(vtx.t[0]))  ) - fVertexTSize;
  double t_max = TMath::Max(  TMath::MaxElement(t_PV.size(), &t_PV[0]), TMath::MaxElement(t_PU.size(), &t_PU[0])  );
  t_max = TMath::Max(t_max, TMath::MaxElement(vtx.getSize(), &(vtx.t[0]))  ) + fVertexTSize;

  double z_min = TMath::Min(  TMath::MinElement(z_PV.size(), &z_PV[0]), TMath::MinElement(z_PU.size(), &z_PU[0])  );
  z_min = TMath::Min(z_min, TMath::MinElement(vtx.getSize(), &(vtx.z[0]))  ) - 5;
  double z_max = TMath::Max(  TMath::MaxElement(z_PV.size(), &z_PV[0]), TMath::MaxElement(z_PU.size(), &z_PU[0])  );
  z_max = TMath::Max(z_max, TMath::MaxElement(vtx.getSize(), &(vtx.z[0]))  ) + 5;

  auto c_2Dspace = new TCanvas("c_2Dspace", "c_2Dspace", 800, 600);

  TGraphErrors* gr_PVtks = new TGraphErrors(t_PV.size(), &t_PV[0], &z_PV[0], &dt_PV[0], &dz_PV[0]);
  gr_PVtks->SetTitle(Form("Clustering space - #beta = %.6f", beta));
  gr_PVtks->GetXaxis()->SetTitle("t CA [ps]");
  gr_PVtks->GetXaxis()->SetLimits(t_min, t_max);
  gr_PVtks->GetYaxis()->SetTitle("z CA [mm]");
  gr_PVtks->GetYaxis()->SetRangeUser(z_min, z_max);
  gr_PVtks->SetMarkerStyle(4);
  gr_PVtks->Draw("APE1");

  TGraphErrors* gr_PUtks = new TGraphErrors(t_PU.size(), &t_PU[0], &z_PU[0], &dt_PU[0], &dz_PU[0]);
  gr_PUtks->SetMarkerStyle(3);
  gr_PUtks->Draw("PE1");

  TGraph* gr_vtx = new TGraph(vtx.getSize(), &(vtx.t[0]), &(vtx.z[0]));
  gr_vtx->SetMarkerStyle(28);
  gr_vtx->SetMarkerColor(2);
  gr_vtx->SetMarkerSize(2.);
  gr_vtx->Draw("PE1");

  // auto leg = new TLegend(0.1, 0.1);
  // leg->AddEntry(gr_PVtks, "PV tks", "ep");
  // leg->AddEntry(gr_PUtks, "PU tks", "ep");
  // leg->AddEntry(gr_vtx, "Cluster center", "p");
  // leg->Draw();

  c_2Dspace->SetGrid();
  c_2Dspace->SaveAs(Form("~/Desktop/debug/c_2Dspace_beta%010.0f-%s%d.png", 1E7*beta, flag, n_it));

  delete c_2Dspace;
}
