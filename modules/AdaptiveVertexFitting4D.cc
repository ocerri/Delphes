/** \class AdaptiveVertexFitting4D
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  \authors O. Cerri

Example for delphes card:
 ##################################
 # Vertex cluster fitting
 ##################################

 module AdaptiveVertexFitting4D AdaptiveVertexFitting4D {
   set InputArray VertexFinderDAClusterizerZT/vertices

   set OutputArray tracks
   set VertexOutputArray vertices

   set Verbose 0
 }

 *
 */


#include "modules/AdaptiveVertexFitting4D.h"
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
#include "TFitter.h"
#include "TMinuit.h"
// replace ln18 of makefile
// LIBS = -lRooFitCore -lRooFit -lMinuit -lHtml -lPyROOT -lFoam -lRooStats -lTreePlayer -lTMVA
// DELPHES_LIBS = $(shell $(RC) --libs) -lEG $(SYSLIBS) $(LIBS)

#include "TAxis.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TColor.h"
#include "TLegend.h"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>

using namespace std;
using namespace AVF_4D;

//------------------------------------------------------------------------------
// Global variable to pass info to the VertexFinderDAClusterizerZTVertexFinderDAClusterizerZT
AdaptiveVertexFitting4D::tracks_t * global_tks;
double global_z_v = 0;
double global_chi2 = 0;


//------------------------------------------------------------------------------

AdaptiveVertexFitting4D::AdaptiveVertexFitting4D()
{
  fVerbose = 0;
  fMaxIterations = 0;
  fBetaMax = 0;
  fCoolingFactor = 0;
  fChi2_0 = 0;
}

//------------------------------------------------------------------------------

AdaptiveVertexFitting4D::~AdaptiveVertexFitting4D()
{
}

//------------------------------------------------------------------------------

void AdaptiveVertexFitting4D::Init()
{
  fVerbose         = GetInt("Verbose", 0);

  fMaxIterations   = GetInt("MaxIterations", 100);

  fBetaMax         = GetDouble("BetaMax", 1.);

  fCoolingFactor   = GetDouble("CoolingFactor", 0.8); // Multiply T so to cooldown must be <1
  fChi2_0   = GetDouble("Chi2_0", 4); // Multiply T so to cooldown must be <1


  fInputArray = ImportArray(GetString("InputArray", "VertexFinderDAClusterizerZT/vertices"));
  fItInputArray = fInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void AdaptiveVertexFitting4D::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void AdaptiveVertexFitting4D::Process()
{
  fItInputArray->Reset();
  Candidate * mother;
  while((mother = static_cast<Candidate*>(fItInputArray->Next())))
  {
    Candidate * cluster = static_cast<Candidate*>(mother->Clone());

    vertex_t vtx(cluster->Position.T()/c_light, cluster->Position.Z());
    if(fVerbose>2)
    {
      cout << endl << "-------------------- Vertex " << cluster->ClusterIndex << " --------------------------------"<< endl;
      cout << Form("Start: {0, 0, %.2f +/- %.2f} mm - %.0f +/- %.0f ps", vtx.z, vtx.sigma_z, vtx.t, vtx.sigma_t) << endl;
      cout << "N tracks: " << cluster->GetCandidates()->GetEntriesFast() << endl;
    }

    tracks_t tks;
    fill(cluster, tks);

    double beta = 0;

    fit(beta, tks, vtx);

    if(fVerbose>2)
    {
      cout << "chi2: " << vtx.chi2 << endl;
      cout << Form("End: {0, 0, %.2f +/- %.2f} mm - %.2f +/- %.2f ps", vtx.z, vtx.sigma_z, vtx.t, vtx.sigma_t) << endl;
    }
    cluster->Chi2 = vtx.chi2;

    cluster->Position.SetT(vtx.t * c_light);
    cluster->PositionError.SetT(vtx.sigma_t * c_light);

    cluster->Position.SetZ(vtx.z);
    cluster->PositionError.SetZ(vtx.sigma_z);

    for(unsigned int i = 0 ; i < tks.getSize(); i++)
    {
      Candidate *tr = static_cast<Candidate*>(tks.tt[i]->Clone());;
      tr->InitialPosition.SetT(cluster->Position.T());
      tr->InitialPosition.SetZ(tks.z(i, vtx.t));

      fOutputArray->Add(tr);
    }

    fVertexOutputArray->Add(cluster);
  }

}


//------------------------------------------------------------------------------
// Fill tks with the candidates associated to the vtx
void AdaptiveVertexFitting4D::fill(Candidate* vtx, tracks_t &tks)
{
  // loop over tracks belonging to this vertex
  TIter TrackIt(vtx->GetCandidates());
  TrackIt.Reset();

  Candidate *tr;
  unsigned int i = 0;
  while((tr = static_cast<Candidate*>(TrackIt.Next())))
  {
    double aux = 0;

    tks.t_out.push_back( tr->Position.T()/c_light );
    aux = tr->ErrorT / c_light;
    tks.s2t_out.push_back( aux*aux );

    tks.z_out.push_back(tr->Position.Z());
    tks.s2z_out.push_back(tr->PositionError.Z());

    tks.Pt.push_back(tr->Momentum.Pt());
    tks.s2_Pt.push_back(tr->ErrorPT * tr->ErrorPT);

    tks.ctg_theta.push_back(tr->CtgTheta);
    tks.E.push_back( sqrt(tr->Mass*tr->Mass + tr->Momentum.Pt()*tr->Momentum.Pt()*(1+ tr->CtgTheta*tr->CtgTheta)) );
    tks.beta_z.push_back( tks.Pt[i] * tks.ctg_theta[i] / tks.E[i] );
    tks.M.push_back(tr->Mass);

    // tks.wz.push_back( 1. );

    tks.tt.push_back( tr );

    if(fVerbose > 20 && vtx->GetCandidates()->GetEntriesFast()< 5) tks.dump(i);
    i++;
  }

  return;
}

// -----------------------------------------------------------------------------
// Fit the vertex and update the vertex position
void AdaptiveVertexFitting4D::fit(double beta, tracks_t &tks, vertex_t &vtx)
{
  TFitter minuit(2);
  minuit.SetDefaultFitter("Minuit");
  minuit.GetMinuit()->SetPrintLevel(-1);
  minuit.SetFCN(NegativeLogLikelihood);

  // Beta not used for the moment
  minuit.SetParameter(0, "beta", 0, 0, -beta, beta);
  minuit.FixParameter(0);

  minuit.SetParameter(1, "t_v", vtx.t, vtx.sigma_t, 5000, -5000);
  global_tks = &tks;

  minuit.ExecuteCommand("MIGRAD", new double(0), 0);

  vtx.t = minuit.GetParameter(1);
  vtx.sigma_t = minuit.GetParError(1);

  vtx.z = global_z_v;
  vtx.sigma_z = 1./sqrt(tks.sum_wz);

  vtx.chi2 = global_chi2;

  if(fVerbose > 4 && vtx.chi2 < 1.)
  {
    for(unsigned int i = 0; i < tks.getSize(); i++)
    {
      cout << "z(tv) = " << tks.z(i, vtx.t) << endl;
      cout << "dz2_(tv) = " << tks.s2z(i, vtx.t) << endl;
      double aux = c_light * tks.beta_z[i];
      aux = aux * aux * tks.s2t_out[i];
      cout << Form("Breakdown : %.2f + %.2f", tks.s2z_out[i], aux) << flush;

      aux = c_light * (vtx.t - tks.t_out[i]) * tks.ctg_theta[i] * ( 1 - tks.Pt[i]*tks.Pt[i]*(1+tks.ctg_theta[i]*tks.ctg_theta[i])/(tks.E[i]*tks.E[i]) ) / tks.E[i];
      aux = aux*aux*tks.s2_Pt[i];
      cout << Form("+ %.2f", aux) << endl;
      cout << "beta_z: " << tks.beta_z[i] << endl;
      cout << "Mass: " << tks.M[i] << endl;
      cout << "Energy: " << tks.E[i] << endl;
      cout << "CtgTheta: " << tks.ctg_theta[i] << endl;
      cout << Form("Pt: %.3f +/- %.3f ", tks.Pt[i], sqrt(tks.s2_Pt[i])) << endl;
      cout << "eta: " << tks.tt[i]->Momentum.Eta() << endl;
    }
  }
}

// -----------------------------------------------------------------------------
// Global Function for the minimization

void NegativeLogLikelihood(Int_t&, Double_t*, Double_t&f, Double_t*par, Int_t )
{
  double t_v = par[1];
  // double beta = par[0];

  unsigned int nt = global_tks->getSize();
  std::vector<double> z(nt);
  std::vector<double> wz(nt);

  global_z_v = 0;
  global_tks->sum_wz = 0;
  for(unsigned int i = 0; i< nt; i++)
  {
    z[i] = global_tks->z(i, t_v);
    wz[i] = 1./global_tks->s2z(i, t_v);

    global_z_v += wz[i]*z[i];
    global_tks->sum_wz += wz[i];
  }
  global_z_v /= global_tks->sum_wz;

  global_chi2 = 0;
  for(unsigned int i = 0; i< nt; i++)
  {
    global_chi2 += (z[i] - global_z_v)*(z[i] - global_z_v) * wz[i];
  }

  f = global_chi2;
}
