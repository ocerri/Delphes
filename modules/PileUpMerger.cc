/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class PileUpMerger
 *
 *  Merges particles from pile-up sample into event
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *  \author O. Cerri - Caltech, Pasadena
 *
 */

#include "modules/PileUpMerger.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesTF2.h"
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

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

PileUpMerger::PileUpMerger() :
  fFunction(0), fReader(0), fItInputArray(0)
{
  fFunction = new DelphesTF2;
}


//------------------------------------------------------------------------------

PileUpMerger::~PileUpMerger()
{
  delete fFunction;
}

//------------------------------------------------------------------------------

void PileUpMerger::Init()
{
  const char *fileName;

  fVerbose = GetInt("Verbose", 0);

  fPileUpDistribution = GetInt("PileUpDistribution", 0);

  fMeanPileUp  = GetDouble("MeanPileUp", 10);

  fZVertexSpread = GetDouble("ZVertexSpread", 0.15);
  fTVertexSpread = GetDouble("TVertexSpread", 1.5E-09);

  fInputBeamSpotX = GetDouble("InputBeamSpotX", 0.0);
  fInputBeamSpotY = GetDouble("InputBeamSpotY", 0.0);

  // read vertex smearing formula

  fFunction->Compile(GetString("VertexDistributionFormula", "0.0"));
  fFunction->SetRange(-fZVertexSpread, -fTVertexSpread, fZVertexSpread, fTVertexSpread);

  fileName = GetString("PileUpFile", "MinBias.pileup");
  fReader = new DelphesPileUpReader(fileName);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays
  fParticleOutputArray = ExportArray(GetString("ParticleOutputArray", "stableParticles"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void PileUpMerger::Finish()
{
  if(fReader) delete fReader;
}

//------------------------------------------------------------------------------

void PileUpMerger::Process()
{
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  Int_t pid, nch, nvtx = 0;
  Float_t x, y, z, t, vx, vy;
  Float_t px, py, pz, e, pt;
  Double_t dz, dphi, dt, sumpt2;
  Int_t numberOfEvents, event, numberOfParticles;
  Long64_t allEntries, entry;
  Candidate *candidate, *vertex;
  DelphesFactory *factory;

  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();

  // --- Deal with primary vertex first  ------

  fFunction->GetRandom2(dz, dt);
  dz *= 1.0E3; // necessary in order to make z in mm
  dt *= c_light*1.0E3; // necessary in order to make t in mm/c

  if(fVerbose)
  {
    cout << "-------------PU Merger---------------------" << endl;
    cout << Form("PV position: z=%.3f mm, t=%.0f ps", dz, dt*1E9/c_light) << endl;
  }

  vx = 0.0;
  vy = 0.0;
  nch = 0;
  sumpt2 = 0.0;
  Double_t sumpt2_sel = 0.0;

  factory = GetFactory();
  vertex = factory->NewCandidate();

  numberOfParticles = fInputArray->GetEntriesFast();

  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    double pt2 = candidate->Momentum.Pt() * candidate->Momentum.Pt();
    vx += candidate->Position.X() * pt2;
    vy += candidate->Position.Y() * pt2;
    z = candidate->Position.Z();
    t = candidate->Position.T();

    candidate->Position.SetZ(z + dz);
    candidate->Position.SetT(t + dt);

    candidate->IsPU = 0;
    candidate->ClusterIndex = nvtx;
    fParticleOutputArray->Add(candidate);

    if(TMath::Abs(candidate->Charge) >  1.0E-9)
    {
      sumpt2 += pt2;
      nch++;
      vertex->AddCandidate(candidate);

      if(candidate->Momentum.Pt() < 50 && candidate->Momentum.Pt() > 0.5)
      {
        if(fabs(candidate->PID) < 5000) //is SM
        {
          sumpt2_sel += pt2;
        }
      }
    }

  }

  vertex->ClusterIndex = nvtx;
  nvtx++;
  if(sumpt2 > 0)
  {
    vx /= sumpt2;
    vy /= sumpt2;
  }
  vertex->Position.SetXYZT(vx, vy, dz, dt);
  vertex->PositionError.SetXYZT(0,0,0,0);
  vertex->ClusterNDF = nch;
  vertex->SumPT2 = sumpt2_sel;
  vertex->GenSumPT2 = sumpt2;
  fVertexOutputArray->Add(vertex);

  // --- Then with pile-up vertices  ------

  switch(fPileUpDistribution)
  {
    case 0:
      numberOfEvents = gRandom->Poisson(fMeanPileUp);
      break;
    case 1:
      numberOfEvents = gRandom->Integer(2*fMeanPileUp + 1);
      break;
    case 2:
      numberOfEvents = fMeanPileUp;
      break;
    default:
      numberOfEvents = gRandom->Poisson(fMeanPileUp);
      break;
  }

  allEntries = fReader->GetEntries();


  // --- Pile-up vertex smearing--------------------
  for(event = 0; event < numberOfEvents; ++event)
  {
    do
    {
      entry = TMath::Nint(gRandom->Rndm()*allEntries);
    }
    while(entry >= allEntries);

    fReader->ReadEntry(entry);


    fFunction->GetRandom2(dz, dt);
    dz *= 1.0E3; // necessary in order to make z in mm
    dt *= c_light*1.0E3; // necessary in order to make t in mm/c
    if(fVerbose)
    {
      cout << Form("PV position: z=%.3f mm, t=%.0f ps", dz, dt*1E9/c_light) << endl;
    }

    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

    vx = 0.0;
    vy = 0.0;
    numberOfParticles = 0;
    sumpt2 = 0.0;
    Double_t sumpt2_sel = 0.0;
    nch = 0;

    vertex = factory->NewCandidate();

    while(fReader->ReadParticle(pid, x, y, z, t, px, py, pz, e))
    {
      candidate = factory->NewCandidate();

      candidate->PID = pid;

      candidate->Status = 1;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -1E6;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -1E6;

      candidate->IsPU = 1;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);
      pt = candidate->Momentum.Pt();

      x -= fInputBeamSpotX;
      y -= fInputBeamSpotY;
      candidate->Position.SetXYZT(x, y, z + dz, t + dt);
      candidate->Position.RotateZ(dphi);

      vx += x*pt*pt;
      vy += y*pt*pt;

      ++numberOfParticles;
      if(TMath::Abs(candidate->Charge) >  1.0E-9)
      {
        nch++;
        sumpt2 += pt*pt;
        vertex->AddCandidate(candidate);

        if(candidate->Momentum.Pt() < 50 && candidate->Momentum.Pt() > 0.5)
        {
          if(fabs(candidate->PID) < 5000) //is SM
          {
            sumpt2_sel += pt*pt;
          }
        }
      }
      candidate->ClusterIndex = nvtx;
      fParticleOutputArray->Add(candidate);
    }

    if(sumpt2 > 0)
    {
      vx /= sumpt2;
      vy /= sumpt2;
    }

    vertex->Position.SetXYZT(vx, vy, dz, dt);

    vertex->ClusterIndex = nvtx;
    nvtx++;
    vertex->ClusterNDF = nch;
    vertex->SumPT2 = sumpt2_sel;
    vertex->GenSumPT2 = sumpt2;

    vertex->IsPU = 1;

    fVertexOutputArray->Add(vertex);

  }
}
