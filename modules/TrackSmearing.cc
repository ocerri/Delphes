/** \class TrackSmearing
 *
 *  Performs d0, dZ, p, Theta, Phi smearing of tracks.
 *
 *  \authors A. Hart, M. Selvaggi
 *
*/

#include "modules/TrackSmearing.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

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
#include "TFile.h"
#include "TProfile2D.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

TrackSmearing::TrackSmearing() :
  fD0Formula(0), fDZFormula(0), fPFormula(0), fCtgThetaFormula(0), fPhiFormula(0), fItInputArray(0), fZOuterResolution(0)
{
  fD0Formula = new DelphesFormula;
  fDZFormula = new DelphesFormula;
  fPFormula = new DelphesFormula;
  fCtgThetaFormula = new DelphesFormula;
  fPhiFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackSmearing::~TrackSmearing()
{
  if(fD0Formula) delete fD0Formula;
  if(fDZFormula) delete fDZFormula;
  if(fPFormula) delete fPFormula;
  if(fCtgThetaFormula) delete fCtgThetaFormula;
  if(fPhiFormula) delete fPhiFormula;
}

//------------------------------------------------------------------------------

void TrackSmearing::Init()
{

  fZOuterResolution = GetDouble("ZOuterResolution", 0.1); // [mm]
  // read resolution formula

  // !!! IF WE WANT TO KEEP ROOT INPUT !!!
  if (string (GetString("D0ResolutionFormula", "0.0")) != "0.0")
    {
      fD0Formula->Compile(GetString("D0ResolutionFormula", "0.0"));
      fUseD0Formula = true;
    }
  else
    {
      fD0ResolutionFile = GetString("D0ResolutionFile", "errors.root");
      fD0ResolutionHist = GetString("D0ResolutionHist", "d0");
      fUseD0Formula = false;
    }
  if (string (GetString("DZResolutionFormula", "0.0")) != "0.0")
    {
      fDZFormula->Compile(GetString("DZResolutionFormula", "0.0"));
      fUseDZFormula = true;
    }
  else
    {
      fDZResolutionFile = GetString("DZResolutionFile", "errors.root");
      fDZResolutionHist = GetString("DZResolutionHist", "dz");
      fUseDZFormula = false;
    }
  if (string (GetString("PResolutionFormula", "0.0")) != "0.0")
    {
      fPFormula->Compile(GetString("PResolutionFormula", "0.0"));
      fUsePFormula = true;
    }
  else
    {
      fPResolutionFile = GetString("PResolutionFile", "errors.root");
      fPResolutionHist = GetString("PResolutionHist", "p");
      fUsePFormula = false;
    }
  if (string (GetString("CtgThetaResolutionFormula", "0.0")) != "0.0")
    {
      fCtgThetaFormula->Compile(GetString("CtgThetaResolutionFormula", "0.0"));
      fUseCtgThetaFormula = true;
    }
  else
    {
      fCtgThetaResolutionFile = GetString("CtgThetaResolutionFile", "errors.root");
      fCtgThetaResolutionHist = GetString("CtgThetaResolutionHist", "ctgTheta");
      fUseCtgThetaFormula = false;
    }
  if (string (GetString("PhiResolutionFormula", "0.0")) != "0.0")
    {
      fPhiFormula->Compile(GetString("PhiResolutionFormula", "0.0"));
      fUsePhiFormula = true;
    }
  else
    {
      fPhiResolutionFile = GetString("PhiResolutionFile", "errors.root");
      fPhiResolutionHist = GetString("PhiResolutionHist", "phi");
      fUsePhiFormula = false;
    }

  fApplyToPileUp = GetBool("ApplyToPileUp", true);

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // import beamspot
  try
  {
    fBeamSpotInputArray = ImportArray(GetString("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle"));
  }
  catch(runtime_error &e)
  {
    fBeamSpotInputArray = 0;
  }

  // create output array
  fBz = GetDouble("Bz", 0.0);
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
}

//------------------------------------------------------------------------------

void TrackSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TrackSmearing::Process()
{
  Int_t iCandidate = 0;
  TLorentzVector beamSpotPosition;
  Candidate *candidate, *mother;
  Double_t pt, eta, d0, d0Error, dz, dzError, p, pError, trueP, ctgTheta, ctgThetaError, trueCtgTheta, phi, phiError, truePhi;
  Double_t x, y, z, t, px, py, pz, theta;
  TProfile2D *d0ErrorHist = NULL,
             *dzErrorHist = NULL,
             *pErrorHist = NULL,
             *ctgThetaErrorHist = NULL,
             *phiErrorHist = NULL;

  if (!fBeamSpotInputArray || fBeamSpotInputArray->GetSize () == 0)
    beamSpotPosition.SetXYZT(0.0, 0.0, 0.0, 0.0);
  else
  {
    Candidate &beamSpotCandidate = *((Candidate *) fBeamSpotInputArray->At (0));
    beamSpotPosition = beamSpotCandidate.Position;
  }

  if (!fUseD0Formula)
  {
     TFile *fin = TFile::Open (fD0ResolutionFile.c_str ());
     d0ErrorHist = (TProfile2D *) fin->Get (fD0ResolutionHist.c_str ());
     d0ErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUseDZFormula)
  {
     TFile *fin = TFile::Open (fDZResolutionFile.c_str ());
     dzErrorHist = (TProfile2D *) fin->Get (fDZResolutionHist.c_str ());
     dzErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUsePFormula)
  {
     TFile *fin = TFile::Open (fPResolutionFile.c_str ());
     pErrorHist = (TProfile2D *) fin->Get (fPResolutionHist.c_str ());
     pErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUseCtgThetaFormula)
  {
     TFile *fin = TFile::Open (fCtgThetaResolutionFile.c_str ());
     ctgThetaErrorHist = (TProfile2D *) fin->Get (fCtgThetaResolutionHist.c_str ());
     ctgThetaErrorHist->SetDirectory (0);
     fin->Close ();
  }
  if (!fUsePhiFormula)
  {
     TFile *fin = TFile::Open (fPhiResolutionFile.c_str ());
     phiErrorHist = (TProfile2D *) fin->Get (fPhiResolutionHist.c_str ());
     phiErrorHist->SetDirectory (0);
     fin->Close ();
  }

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    if (!fApplyToPileUp && candidate->IsPU)
    {
      fOutputArray->Add(candidate);
      iCandidate++;
      continue;
    }

    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->InitialPosition;

    pt = momentum.Pt();
    eta = momentum.Eta();

    d0 = candidate->D0;
    dz = candidate->DZ;

    p = trueP = candidate->P;
    ctgTheta = trueCtgTheta = candidate->CtgTheta;
    phi = truePhi = candidate->Phi;

    if (fUseD0Formula)
      d0Error = fD0Formula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < d0ErrorHist->GetXaxis ()->GetXmax () ? d0ErrorHist->GetXaxis ()->FindBin (pt)
            : d0ErrorHist->GetXaxis ()->GetBinCenter (d0ErrorHist->GetXaxis ()->GetNbins ());
       ybin = d0ErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       d0Error = d0ErrorHist->GetBinContent (xbin, ybin);
       if (!d0Error)
          d0Error = -1.0;
    }
    if (d0Error < 0.0)
      continue;

    if (fUseDZFormula)
      dzError = fDZFormula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < dzErrorHist->GetXaxis ()->GetXmax () ? dzErrorHist->GetXaxis ()->FindBin (pt)
            : dzErrorHist->GetXaxis ()->GetBinCenter (dzErrorHist->GetXaxis ()->GetNbins ());
       ybin = dzErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       dzError = dzErrorHist->GetBinContent (xbin, ybin);
       if (!dzError)
          dzError = -1.0;
    }
    if (dzError < 0.0)
      continue;

    if (fUsePFormula)
      pError = fPFormula->Eval(pt, eta) * p;
    else
    {
       Int_t xbin, ybin;

       xbin = pt < pErrorHist->GetXaxis ()->GetXmax () ? pErrorHist->GetXaxis ()->FindBin (pt)
            : pErrorHist->GetXaxis ()->GetBinCenter (pErrorHist->GetXaxis ()->GetNbins ());
       ybin = pErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       pError = pErrorHist->GetBinContent (xbin, ybin) * p;
       if (!pError)
          pError = -1.0;
    }
    if (pError < 0.0)
      continue;

    if (fUseCtgThetaFormula)
      ctgThetaError = fCtgThetaFormula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < ctgThetaErrorHist->GetXaxis ()->GetXmax () ? ctgThetaErrorHist->GetXaxis ()->FindBin (pt)
            : ctgThetaErrorHist->GetXaxis ()->GetBinCenter (ctgThetaErrorHist->GetXaxis ()->GetNbins ());
       ybin = ctgThetaErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       ctgThetaError = ctgThetaErrorHist->GetBinContent (xbin, ybin);
       if (!ctgThetaError)
         ctgThetaError = -1.0;
    }
    if (ctgThetaError < 0.0)
      continue;

    if (fUsePhiFormula)
      phiError = fPhiFormula->Eval(pt, eta);
    else
    {
       Int_t xbin, ybin;

       xbin = pt < phiErrorHist->GetXaxis ()->GetXmax () ? phiErrorHist->GetXaxis ()->FindBin (pt)
            : phiErrorHist->GetXaxis ()->GetBinCenter (phiErrorHist->GetXaxis ()->GetNbins ());
       ybin = phiErrorHist->GetYaxis ()->FindBin (TMath::Abs (eta));
       phiError = phiErrorHist->GetBinContent (xbin, ybin);
       if (!phiError)
          phiError = -1.0;
    }
    if (phiError < 0.0)
      continue;

    mother = candidate;
    candidate = static_cast<Candidate*>(mother->Clone());

    if (d0Error != 0)
    {
      d0 = gRandom->Gaus(d0, d0Error);
      candidate->D0 = d0;
    }
    candidate->ErrorD0 = d0Error;

    if (dzError != 0)
    {
      candidate->DZ = gRandom->Gaus(dz, dzError);
    }
    candidate->ErrorDZ = dzError;

    if(fZOuterResolution > 0)
    {
      candidate->Position.SetZ(candidate->Position.Z() + gRandom->Gaus(0, fZOuterResolution));
      candidate->PositionError.SetZ(fZOuterResolution);
    }

    // if(pError*ctgThetaError*phiError > 0.)
    // {
    //   cout << "To be validated" << endl;
    //   p = gRandom->Gaus(p, pError);
    //   ctgTheta = gRandom->Gaus(ctgTheta, ctgThetaError);
    //   phi = gRandom->Gaus(phi, phiError);
    //
    //   if(p < 0.0) continue;
    //   while (phi > TMath::Pi ()) phi -= TMath::TwoPi ();
    //   while (phi <= -TMath::Pi ()) phi += TMath::TwoPi ();
    //
    //
    //   candidate->P = p;
    //   candidate->CtgTheta = ctgTheta;
    //   candidate->Phi = phi;
    //
    //   theta = TMath::ACos(ctgTheta / TMath::Sqrt (1.0 + ctgTheta * ctgTheta));
    //   candidate->Momentum.SetPx (p * TMath::Cos (phi) * TMath::Sin (theta));
    //   candidate->Momentum.SetPy (p * TMath::Sin (phi) * TMath::Sin (theta));
    //   candidate->Momentum.SetPz (p * TMath::Cos (theta));
    //   candidate->Momentum.SetE (candidate->Momentum.Pt () * TMath::CosH (eta));
    //   candidate->PT = candidate->Momentum.Pt ();
    //
    //   x = position.X ();
    //   y = position.Y ();
    //   z = position.Z ();
    //   t = position.T ();
    //   px = candidate->Momentum.Px ();
    //   py = candidate->Momentum.Py ();
    //   pz = candidate->Momentum.Pz ();
    //   pt = candidate->Momentum.Pt ();
    //
    //   // -- solve for delta: d0' = ( (x+delta)*py' - (y+delta)*px' )/pt'
    //
    //   candidate->InitialPosition.SetX (x + ((px * y - py * x + d0 * pt) / (py - px)));
    //   candidate->InitialPosition.SetY (y + ((px * y - py * x + d0 * pt) / (py - px)));
    //   x = candidate->InitialPosition.X ();
    //   y = candidate->InitialPosition.Y ();
    //   candidate->InitialPosition.SetZ (z + ((pz * (px * (x - beamSpotPosition.X ()) + py * (y - beamSpotPosition.Y ())) + pt * pt * (dz - z)) / (pt * pt)));
    //
    //   candidate->InitialPosition.SetT(t);
    //
    //   candidate->ErrorP = pError;
    //   candidate->ErrorCtgTheta = ctgThetaError;
    //   candidate->ErrorPhi = phiError;
    //   candidate->ErrorPT = ptError (p, ctgTheta, pError, ctgThetaError);
    //   candidate->TrackResolution = pError/p;
    // }

    candidate->Xd = d0 * TMath::Sin(phi);
    candidate->Yd = - d0 * TMath::Cos(phi);
    candidate->Zd = candidate->DZ;

    const Double_t c_light = 2.99792458E8;
    // Recompute the new track length
    Double_t rh = candidate->Momentum.Pt() / (candidate->Charge * fBz) * 1.0E9/c_light;        // in [m]
    rh *= 1E3; // [mm]

    Double_t xc = (d0 + rh)*TMath::Sin(phi);
    Double_t yc = -(d0 + rh)*TMath::Cos(phi);


    candidate->Td = -9999*c_light*1E-3;

    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);

    iCandidate++;
  }
}

Double_t TrackSmearing::ptError (const Double_t p, const Double_t ctgTheta, const Double_t dP, const Double_t dCtgTheta)
{
  Double_t a, b;
  a = (p * p * ctgTheta * ctgTheta * dCtgTheta * dCtgTheta) / ((ctgTheta * ctgTheta + 1) * (ctgTheta * ctgTheta + 1) * (ctgTheta * ctgTheta + 1));
  b = (dP * dP) / (ctgTheta * ctgTheta + 1);
  return sqrt (a + b);
}

//------------------------------------------------------------------------------
