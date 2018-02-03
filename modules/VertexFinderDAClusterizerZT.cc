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
  fMinTrackWeight(0),
  fUseTc(0),
  fBetaMax(0),
  fBetaStop(0),
  fBetaPurge(0),
  fVertexSpaceSize(0),
  fVertexTimeSize(0),
  fCoolingFactor(0),
  fDzCutOff(0),
  fD0CutOff(0),
  fDtCutOff(0),
  fMinPT(0),
  fUniqueTrkWeight(0),
  fDzMerge(0),
  fDtMerge(0)
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
  fMinTrackWeight  = GetDouble("MinTrackWeight", 0.5);
  fUseTc           = GetBool("UseTc", 1);

  fBetaMax         = GetDouble("BetaMax ", 0.1);
  fBetaStop        = GetDouble("BetaStop", 1.0);
  fBetaPurge       = GetDouble("BetaPurge", 1.0);

  fVertexSpaceSize = GetDouble("VertexSpaceSize", 0.5); //in mm
  fVertexTimeSize  = GetDouble("VertexTimeSize", 10E-12); //in s

  fCoolingFactor   = GetDouble("CoolingFactor", 0.8); // Multiply T so to cooldown must be <1

  // !!FIX defaul values
  // DzCutOff also used as initializer to Z_init to do outlayers, as D0CutOff must be O(1) and find anew way to the flowwing line meaning (add new var)
  fDzCutOff        = GetDouble("DzCutOff", 40);  // don't know yet. For the moment 30*DzCutOff is hard cut off for the considered tracks
  fD0CutOff        = GetDouble("D0CutOff", 3);  // d0/sigma_d0, used to compute the pi (weight) of the track
  fDtCutOff        = GetDouble("DtCutOff", 3);  // Used as first line for outlayer rejection O(1). Put only a parameter maybe and use the other for hard cutoff
  fMinPT           = GetDouble("MinPT", 0.1);

  // !!FIX defaul values
  fUniqueTrkWeight = GetDouble("UniqueTrkWeight", 1);
  fDzMerge         = GetDouble("DzMerge", 3);
  fDtMerge         = GetDouble("D0CutOff", 3);

  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));

  if (fBetaMax > fBetaPurge)
  {
    fBetaPurge = fBetaMax;
    if (fVerbose)
    {
      cout << "BetaPurge set to " << fBetaPurge << endl;
    }
  }

  if (fBetaPurge > fBetaStop)
  {
    fBetaStop = TMath::Min(1., fBetaMax);
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
  Candidate *candidate, *track;
  TObjArray *ClusterArray = new TObjArray;
  TIterator *ItClusterArray;

  Int_t ivtx = 0;

  fInputArray->Sort();

  TLorentzVector pos, mom;
  if (fVerbose)
  {
     cout<<" start processing vertices ..."<<endl;
     cout<<" Found "<<fInputArray->GetEntriesFast()<<" input tracks"<<endl;
     //loop over input tracks
     fItInputArray->Reset();
     while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
     {
        pos = candidate->InitialPosition;
        mom = candidate->Momentum;

        cout<<"pt: "<<mom.Pt()<<", eta: "<<mom.Eta()<<", phi: "<<mom.Phi()<<", z: "<<candidate->DZ/10<<endl;
     }
  }

  // clusterize tracks
  clusterize(*fInputArray, *ClusterArray);

  if (fVerbose){std::cout <<  " clustering returned  "<< ClusterArray->GetEntriesFast() << " clusters  from " << fInputArray->GetEntriesFast() << " selected tracks" <<std::endl;}
  // ----------------HERE HERE ----------------//

  //loop over vertex candidates
  ItClusterArray = ClusterArray->MakeIterator();
  ItClusterArray->Reset();
  while((candidate = static_cast<Candidate*>(ItClusterArray->Next())))
  {

     double meantime = 0.;
     double expv_x2 = 0.;
     double normw = 0.;
     double errtime = 0;

     double meanpos = 0.;
     double meanerr2 = 0.;
     double normpos = 0.;
     double errpos = 0.;

     double sumpt2 = 0.;

     int itr = 0;

     if(fVerbose)cout<<"this vertex has: "<<candidate->GetCandidates()->GetEntriesFast()<<" tracks"<<endl;

     // loop over tracks belonging to this vertex
     TIter it1(candidate->GetCandidates());
     it1.Reset();

     while((track = static_cast<Candidate*>(it1.Next())))
     {
        itr++;
        // TBC: the time is in ns for now TBC
        double t = track->InitialPosition.T()/c_light;
        double dt = track->ErrorT/c_light;
        const double time = t;
        const double inverr = 1.0/dt;
        meantime += time*inverr;
        expv_x2  += time*time*inverr;
        normw    += inverr;

        // compute error position TBC
        const double pt = track->Momentum.Pt();
        const double z = track->DZ/10.0;
        const double err_pt = track->ErrorPT;
        const double err_z = track->ErrorDZ;

        const double wi = (pt/(err_pt*err_z))*(pt/(err_pt*err_z));
        meanpos += z*wi;

        meanerr2 += err_z*err_z*wi;
        normpos += wi;
        sumpt2 += pt*pt;

        // while we are here store cluster index in tracks
        track->ClusterIndex = ivtx;
     }

     meantime = meantime/normw;
     expv_x2 = expv_x2/normw;
     errtime = TMath::Sqrt((expv_x2 - meantime*meantime)/itr);
     meanpos = meanpos/normpos;
     meanerr2 = meanerr2/normpos;
     errpos = TMath::Sqrt(meanerr2/itr);

     candidate->Position.SetXYZT(0.0, 0.0, meanpos*10.0 , meantime*c_light);
     candidate->PositionError.SetXYZT(0.0, 0.0, errpos*10.0 , errtime*c_light);
     candidate->SumPT2 = sumpt2;
     candidate->ClusterNDF = itr;
     candidate->ClusterIndex = ivtx;

     fVertexOutputArray->Add(candidate);

     ivtx++;

     if (fVerbose){
     std::cout << "x,y,z";
       std::cout << ",t";
       std::cout << "=" << candidate->Position.X()/10.0 <<" " << candidate->Position.Y()/10.0 << " " <<  candidate->Position.Z()/10.0;
       std::cout << " " << candidate->Position.T()/c_light;

       std::cout << std::endl;
       std::cout << "sumpt2 " << candidate->SumPT2<<endl;

       std::cout << "ex,ey,ez";
       std::cout << ",et";
       std::cout << "=" << candidate->PositionError.X()/10.0 <<" " << candidate->PositionError.Y()/10.0 << " " <<  candidate->PositionError.Z()/10.0;
       std::cout << " " << candidate->PositionError.T()/c_light;
       std::cout << std::endl;

      }
   }// end of cluster loop


    if(fVerbose){
      std::cout << "PrimaryVertexProducerAlgorithm::vertices candidates =" << ClusterArray->GetEntriesFast() << std::endl;
    }

    //TBC maybe this can be done later
    // sort vertices by pt**2  vertex (aka signal vertex tagging)
    /*if(pvs.size()>1){
      sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }
     */

  delete ClusterArray;

}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::clusterize(const TObjArray &tracks, TObjArray &clusters)
{
  if(fVerbose) {
    cout << "###################################################" << endl;
    cout << "# VertexFinderDAClusterizerZT::clusterize   nt="<<tracks.GetEntriesFast() << endl;
    cout << "###################################################" << endl;
  }

  vector< Candidate* > pv = vertices();

  // -----HERE HERE -----//

  if(fVerbose){ cout << "# VertexFinderDAClusterizerZT::clusterize   pv.size="<<pv.size() << endl;  }
  if (pv.size()==0){ return;  }

  // convert into vector of candidates
  //TObjArray *ClusterArray = pv.begin()->GetCandidates();
  //Candidate *aCluster = static_cast<Candidate*>(&(pv.at(0)));
  Candidate *aCluster = pv.at(0);

  // fill into clusters and merge


  if( fVerbose ) {
      std::cout << '\t' << 0;
      std::cout << ' ' << (*pv.begin())->Position.Z()/10.0 << ' ' << (*pv.begin())->Position.T()/c_light << std::endl;
    }

  for(vector<Candidate*>::iterator k=pv.begin()+1; k!=pv.end(); k++){
    if( fVerbose ) {
      std::cout << '\t' << std::distance(pv.begin(),k);
      std::cout << ' ' << (*k)->Position.Z() << ' ' << (*k)->Position.T() << std::endl;
    }


    // TBC - check units here
    if ( std::abs((*k)->Position.Z() - (*(k-1))->Position.Z())/10.0 > (2*fVertexSpaceSize) ||
         std::abs((*k)->Position.T() - (*(k-1))->Position.Z())/c_light > 2*0.010 ) {
      // close a cluster
      clusters.Add(aCluster);
      //aCluster.clear();
    }
    //for(unsigned int i=0; i<k->GetCandidates().GetEntriesFast(); i++){
      aCluster = *k;
    //}

  }
  clusters.Add(aCluster);

  if(fVerbose) { std::cout << "# VertexFinderDAClusterizerZT::clusterize clusters.size="<<clusters.GetEntriesFast() << std::endl; }

}

//------------------------------------------------------------------------------

vector< Candidate* > VertexFinderDAClusterizerZT::vertices()
{
  Candidate *candidate;
  UInt_t clusterIndex = 0;
  vector< Candidate* > clusters;

  tracks_t tks;
  fill(tks);

  unsigned int nt=tks.getSize();
  double rho0=0.0;  // start with no outlier rejection

  if (nt == 0) return clusters;

  vertex_t vtx; // the vertex prototypes
  // initialize:single vertex at infinite temperature
  vtx.addItem(0, 0, 1);

  // Fit the vertex at T=inf and return first critical temperature
  double beta=beta0(tks, vtx);

  unsigned int niter=0;
  double delta2 = 0;
  do  {
    delta2 = update(beta, tks, vtx, rho0);
    niter++;
  }
  while (delta2 > 1. &&  niter < maxIterations_)
  // --- HERE HERE ---///

  // annealing loop, stop when T<Tmin  (i.e. beta>1/Tmin)
  while(beta<fBetaMax){

    if(fUseTc){
      update1(beta, tks,y);
      while(merge(y,beta)){update1(beta, tks,y);}
      split(beta, tks,y);
      beta=beta/fCoolingFactor;
    }
    else{
      beta=beta/fCoolingFactor;
      splitAll(y);
    }

   // make sure we are not too far from equilibrium before cooling further
   niter=0; while((update1(beta, tks,y)>1.e-6)  && (niter++ < fMaxIterations)){ }

  }

  if(fUseTc){
    // last round of splitting, make sure no critical clusters are left
    update1(beta, tks,y);
    while(merge(y,beta)){update1(beta, tks,y);}
    unsigned int ntry=0;
    while( split(beta, tks,y) && (ntry++<10) ){
      niter=0;
      while((update1(beta, tks,y)>1.e-6)  && (niter++ < fMaxIterations)){}
      merge(y,beta);
      update1(beta, tks,y);
    }
  }else{
    // merge collapsed clusters
    while(merge(y,beta)){update1(beta, tks,y);}
    if(fVerbose ){ cout << "dump after 1st merging " << endl;  dump(beta,y,tks);}
  }

  // switch on outlier rejection
  rho0=1./nt; for(vector<vertex_t>::iterator k=y.begin(); k!=y.end(); k++){ k->pk =1.; }  // democratic
  niter=0; while((update2(beta, tks,y,rho0, fDzCutOff) > 1.e-8)  && (niter++ < fMaxIterations)){  }
  if(fVerbose  ){ cout << "rho0=" << rho0 <<   " niter=" << niter <<  endl; dump(beta,y,tks);}


  // merge again  (some cluster split by outliers collapse here)
  while(merge(y)){}
  if(fVerbose  ){ cout << "dump after 2nd merging " << endl;  dump(beta,y,tks);}


  // continue from freeze-out to Tstop (=1) without splitting, eliminate insignificant vertices
  while(beta<=fBetaStop){
    while(purge(y,tks,rho0, beta, fDzCutOff)){
      niter=0; while((update2(beta, tks, y, rho0, fDzCutOff) > 1.e-6)  && (niter++ < fMaxIterations)){  }
    }
    beta/=fCoolingFactor;
    niter=0; while((update2(beta, tks, y, rho0, fDzCutOff) > 1.e-6)  && (niter++ < fMaxIterations)){  }
  }


  //   // new, one last round of cleaning at T=Tstop
  //   while(purge(y,tks,rho0, beta)){
  //     niter=0; while((update2(beta, tks,y,rho0, fDzCutOff) > 1.e-6)  && (niter++ < fMaxIterations)){  }
  //   }


  if(fVerbose){
   cout << "Final result, rho0=" << rho0 << endl;
   dump(beta,y,tks);
  }


  // select significant tracks and use a TransientVertex as a container
  //GlobalError dummyError;

  // ensure correct normalization of probabilities, should make double assginment reasonably impossible
  for(unsigned int i=0; i<nt; i++){
    tks[i].Z=rho0*exp(-beta*( fDzCutOff*fDzCutOff));
    for(vector<vertex_t>::iterator k=y.begin(); k!=y.end(); k++){
      tks[i].Z += k->pk * exp(-beta*Eik(tks[i],*k));
    }
  }

  for(vector<vertex_t>::iterator k=y.begin(); k!=y.end(); k++){

    DelphesFactory *factory = GetFactory();
    candidate = factory->NewCandidate();

    //cout<<"new vertex"<<endl;
    //GlobalPoint pos(0, 0, k->z);
    double time = k->t;
    double z = k->z;
    //vector< reco::TransientTrack > vertexTracks;
    //double max_tracks_time_err2 = 0;
    double mean = 0.;
    double expv_x2 = 0.;
    double normw = 0.;
    for(unsigned int i=0; i<nt; i++){
      const double invdt = 1.0/std::sqrt(tks[i].dt2);
      if(tks[i].Z>0){
          double p = k->pk * exp(-beta*Eik(tks[i],*k)) / tks[i].Z;
          if( (tks[i].pi>0) && ( p > 0.5 ) ){
            //std::cout << "pushing back " << i << ' ' << tks[i].tt << std::endl;
            //vertexTracks.push_back(*(tks[i].tt)); tks[i].Z=0;

            candidate->AddCandidate(tks[i].tt); tks[i].Z=0;

            mean     += tks[i].t*invdt*p;
            expv_x2  += tks[i].t*tks[i].t*invdt*p;
            normw    += invdt*p;
          } // setting Z=0 excludes double assignment
      }
    }

    mean = mean/normw;
    expv_x2 = expv_x2/normw;
    const double time_var = expv_x2 - mean*mean;
    const double crappy_error_guess = std::sqrt(time_var);
    /*GlobalError dummyErrorWithTime(0,
                                   0,0,
                                   0,0,0,
                                   0,0,0,crappy_error_guess);*/
    //TransientVertex v(pos, time, dummyErrorWithTime, vertexTracks, 5);


    candidate->ClusterIndex = clusterIndex++;;
    candidate->Position.SetXYZT(0.0, 0.0, z*10.0 , time*c_light);

    // TBC - fill error later ...
    candidate->PositionError.SetXYZT(0.0, 0.0, 0.0 , crappy_error_guess*c_light);

    clusterIndex++;
    clusters.push_back(candidate);
  }


  return clusters;

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
    dz2_o += fVertexSpaceSize*fVertexSpaceSize;
    // when needed add beam spot width (x-y)?? mha??
    dz2_o = 1/dz2_o; //Multipling is faster than dividing all the times
    tks.sum_dz2_o += dz2_o;

    dt2_o = candidate->ErrorT*1.E9/c_light; // [ps]
    dt2_o *= dt2_o;
    dt2_o += fVertexTimeSize*fVertexTimeSize*1.E24; // [ps^2]
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

    double z_displ = (new_z - vtx.z[k])/fVertexSpaceSize;
    double t_displ = (new_t - vtx.t[k])/fVertexTimeSize;
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

static bool VertexFinderDAClusterizerZT::merge(vector<vertex_t> &y, double &beta)
{
  // merge clusters that collapsed or never separated,
  // only merge if the estimated critical temperature of the merged vertex is below the current temperature
  // return true if vertices were merged, false otherwise
  if(y.size()<2)  return false;

  for(vector<vertex_t>::iterator k=y.begin(); (k+1)!=y.end(); k++){
    if ( std::abs((k+1)->z - k->z) < 2.e-3 &&
         std::abs((k+1)->t - k->t) < 2.e-3    ) {
      double rho=k->pk + (k+1)->pk;
      double swE=k->swE+(k+1)->swE - k->pk * (k+1)->pk / rho * ( std::pow((k+1)->z - k->z,2.) +
                                                                 std::pow((k+1)->t - k->t,2.)   );
      double Tc=2*swE/(k->sw+(k+1)->sw);

      if(Tc*beta<1){
  if(rho>0){
    k->z = ( k->pk * k->z + (k+1)->z * (k+1)->pk)/rho;
          k->t = ( k->pk * k->t + (k+1)->t * (k+1)->pk)/rho;
  }else{
    k->z = 0.5*(k->z + (k+1)->z);
          k->t = 0.5*(k->t + (k+1)->t);
  }
  k->pk  = rho;
  k->sw += (k+1)->sw;
  k->swE = swE;
  k->Tc  = Tc;
  y.erase(k+1);
  return true;
      }
    }
  }

  return false;
}

//------------------------------------------------------------------------------

static bool VertexFinderDAClusterizerZT::purge(vector<vertex_t> &y, vector<tracks_t> &tks, double & rho0, const double beta, const double dzCutOff)
{
  // eliminate clusters with only one significant/unique track
  if(y.size()<2)  return false;

  unsigned int nt=tks.size();
  double sumpmin=nt;
  vector<vertex_t>::iterator k0=y.end();
  for(vector<vertex_t>::iterator k=y.begin(); k!=y.end(); k++){
    int nUnique=0;
    double sump=0;
    double pmax=k->pk/(k->pk+rho0*exp(-beta*dzCutOff*dzCutOff));
    for(unsigned int i=0; i<nt; i++){
      if(tks[i].Z > 0){
  double p = k->pk * std::exp(-beta*Eik(tks[i],*k)) / tks[i].Z ;
  sump+=p;
  if( (p > 0.9*pmax) && (tks[i].pi>0) ){ nUnique++; }
      }
    }

    if((nUnique<2)&&(sump<sumpmin)){
      sumpmin=sump;
      k0=k;
    }
  }

  if(k0!=y.end()){
    //cout << "eliminating prototype at " << k0->z << "," << k0->t << " with sump=" << sumpmin << endl;
    //rho0+=k0->pk;
    y.erase(k0);
    return true;
  }else{
    return false;
  }
}


//------------------------------------------------------------------------------

static bool VertexFinderDAClusterizerZT::split(double beta, vector<tracks_t> &tks, vector<vertex_t> &y)
{
  // split only critical vertices (Tc >~ T=1/beta   <==>   beta*Tc>~1)
  // an update must have been made just before doing this (same beta, no merging)
  // returns true if at least one cluster was split

  const double epsilon = 1e-3; // split all single vertices by 10 um
  bool split = false;

  // avoid left-right biases by splitting highest Tc first

  std::vector<std::pair<double, unsigned int> > critical;
  for(unsigned int ik=0; ik<y.size(); ik++){
    if (beta*y[ik].Tc > 1.){
      critical.push_back( make_pair(y[ik].Tc, ik));
    }
  }
  std::stable_sort(critical.begin(), critical.end(), std::greater<std::pair<double, unsigned int> >() );

  for(unsigned int ic=0; ic<critical.size(); ic++){
    unsigned int ik=critical[ic].second;
    // estimate subcluster positions and weight
    double p1=0, z1=0, t1=0, w1=0;
    double p2=0, z2=0, t2=0, w2=0;
    //double sumpi=0;
    for(unsigned int i=0; i<tks.size(); i++){
      if(tks[i].Z>0){
  //sumpi+=tks[i].pi;
  double p=y[ik].pk * exp(-beta*Eik(tks[i],y[ik])) / tks[i].Z*tks[i].pi;
  double w=p/(tks[i].dz2 * tks[i].dt2);
  if(tks[i].z < y[ik].z){
    p1+=p; z1+=w*tks[i].z; t1+=w*tks[i].t; w1+=w;
  }else{
    p2+=p; z2+=w*tks[i].z; t2+=w*tks[i].t; w2+=w;
  }
      }
    }
    if(w1>0){  z1=z1/w1; t1=t1/w1;} else{ z1=y[ik].z-epsilon; t1=y[ik].t-epsilon; }
    if(w2>0){  z2=z2/w2; t2=t2/w2;} else{ z2=y[ik].z+epsilon; t2=y[ik].t+epsilon;}

    // reduce split size if there is not enough room
    if( ( ik   > 0       ) && ( y[ik-1].z>=z1 ) ){ z1=0.5*(y[ik].z+y[ik-1].z); t1=0.5*(y[ik].t+y[ik-1].t); }
    if( ( ik+1 < y.size()) && ( y[ik+1].z<=z2 ) ){ z2=0.5*(y[ik].z+y[ik+1].z); t2=0.5*(y[ik].t+y[ik+1].t); }

    // split if the new subclusters are significantly separated
    if( (z2-z1)>epsilon || std::abs(t2-t1) > epsilon){
      split=true;
      vertex_t vnew;
      vnew.pk = p1*y[ik].pk/(p1+p2);
      y[ik].pk= p2*y[ik].pk/(p1+p2);
      vnew.z  = z1;
      vnew.t  = t1;
      y[ik].z = z2;
      y[ik].t = t2;
      y.insert(y.begin()+ik, vnew);

     // adjust remaining pointers
      for(unsigned int jc=ic; jc<critical.size(); jc++){
        if (critical[jc].second>ik) {critical[jc].second++;}
      }
    }
  }

  //  stable_sort(y.begin(), y.end(), clusterLessZ);
  return split;
}

//------------------------------------------------------------------------------

void VertexFinderDAClusterizerZT::splitAll(vector<vertex_t> &y)
{


  const double epsilon=1e-3;      // split all single vertices by 10 um
  const double zsep=2*epsilon;    // split vertices that are isolated by at least zsep (vertices that haven't collapsed)
  const double tsep=2*epsilon;    // check t as well

  vector<vertex_t> y1;

  for(vector<vertex_t>::iterator k=y.begin(); k!=y.end(); k++){
    if ( ( (k==y.begin())|| (k-1)->z < k->z - zsep) && (((k+1)==y.end()  )|| (k+1)->z > k->z + zsep)) {
      // isolated prototype, split
      vertex_t vnew;
      vnew.z  = k->z - epsilon;
      vnew.t  = k->t - epsilon;
      (*k).z  = k->z + epsilon;
      (*k).t  = k->t + epsilon;
      vnew.pk= 0.5* (*k).pk;
      (*k).pk= 0.5* (*k).pk;
      y1.push_back(vnew);
      y1.push_back(*k);

    }else if( y1.empty() || (y1.back().z < k->z -zsep) || (y1.back().t < k->t - tsep) ){
      y1.push_back(*k);
    }else{
      y1.back().z -= epsilon;
      y1.back().t -= epsilon;
      k->z += epsilon;
      k->t += epsilon;
      y1.push_back(*k);
    }
  }// vertex loop

  y=y1;
}
