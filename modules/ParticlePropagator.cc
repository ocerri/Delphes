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


/** \class ParticlePropagator
 *
 *  Propagates charged and neutral particles
 *  from a given vertex to a cylinder defined by its radius,
 *  its half-length, centered at (0,0,0) and with its axis
 *  oriented along the z-axis.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/ParticlePropagator.h"

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

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

ParticlePropagator::ParticlePropagator() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

ParticlePropagator::~ParticlePropagator()
{
}


//------------------------------------------------------------------------------

void ParticlePropagator::Init()
{
  fRadius = GetDouble("Radius", 1.0);
  fRadius2 = fRadius*fRadius;
  fHalfLength = GetDouble("HalfLength", 3.0);
  fBz = GetDouble("Bz", 0.0);
  if(fRadius < 1.0E-2)
  {
    cout << "ERROR: magnetic field radius is too low\n";
    return;
  }
  if(fHalfLength < 1.0E-2)
  {
    cout << "ERROR: magnetic field length is too low\n";
    return;
  }

  fHalfLengthMax = GetDouble("HalfLengthMax", fHalfLength);

  // import array with output from filter/classifier module

  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
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
  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
  fChargedHadronOutputArray = ExportArray(GetString("ChargedHadronOutputArray", "chargedHadrons"));
  fElectronOutputArray = ExportArray(GetString("ElectronOutputArray", "electrons"));
  fMuonOutputArray = ExportArray(GetString("MuonOutputArray", "muons"));
}

//------------------------------------------------------------------------------

void ParticlePropagator::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ParticlePropagator::Process()
{
  const Double_t c_light = 2.99792458E8;

  // // Not used for the moment
  // TLorentzVector beamSpotPosition;
  // if (!fBeamSpotInputArray || fBeamSpotInputArray->GetSize () == 0)
  //   beamSpotPosition.SetXYZT(0.0, 0.0, 0.0, 0.0);
  // else
  // {
  //   Candidate &beamSpotCandidate = *((Candidate *) fBeamSpotInputArray->At(0));
  //   beamSpotPosition = beamSpotCandidate.Position;
  // }
  // Double_t bsx = beamSpotPosition.X()*1.0E-3;
  // Double_t bsy = beamSpotPosition.Y()*1.0E-3;
  // Double_t bsz = beamSpotPosition.Z()*1.0E-3;

  fItInputArray->Reset();
  Candidate *candidate;
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    TLorentzVector candidatePosition = candidate->Position;
    TLorentzVector candidateMomentum = candidate->Momentum;
    Double_t x = candidatePosition.X()*1.0E-3;
    Double_t y = candidatePosition.Y()*1.0E-3;
    Double_t z = candidatePosition.Z()*1.0E-3;

    // check that particle position is inside the cylinder
    if(TMath::Hypot(x, y) > fRadius || TMath::Abs(z) > fHalfLengthMax)
    {
      cout << "Warning: particle produced outside the detector and will not be firther taken into account" << endl;
      continue;
    }

    Double_t q = candidate->Charge;

    Double_t px = candidateMomentum.Px();
    Double_t py = candidateMomentum.Py();
    Double_t pz = candidateMomentum.Pz();
    Double_t pt = candidateMomentum.Pt();
    Double_t pt2 = candidateMomentum.Perp2();
    Double_t e = candidateMomentum.E();

    if(pt2 < 1.0E-9)
    {
      continue;
    }

    if(TMath::Abs(q) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
    {
      // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
      Double_t tmp = px*y - py*x;
      Double_t discr2 = pt2*fRadius2 - tmp*tmp;

      if(discr2 < 0.0)
      {
        // no solutions
        continue;
      }

      tmp = px*x + py*y;
      Double_t t1 = (-tmp + TMath::Sqrt(discr2))/pt2;
      Double_t t2 = (-tmp - TMath::Sqrt(discr2))/pt2;
      Double_t t = (t1 < 0.0) ? t2 : t1;

      Double_t z_t = z + pz*t;
      if(TMath::Abs(z_t) > fHalfLength)
      {
        Double_t t3 = (+fHalfLength - z) / pz;
        Double_t t4 = (-fHalfLength - z) / pz;
        t = (t3 < 0.0) ? t4 : t3;
      }

      Double_t x_t = x + px*t;
      Double_t y_t = y + py*t;
      z_t = z + pz*t;

      Double_t l = TMath::Sqrt( (x_t - x)*(x_t - x) + (y_t - y)*(y_t - y) + (z_t - z)*(z_t - z));

      Candidate* mother = candidate;
      candidate = static_cast<Candidate*>(candidate->Clone());

      candidate->InitialPosition = candidatePosition;
      candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, candidatePosition.T() + t*e*1.0E3);
      candidate->L = l*1.0E3;

      candidate->Momentum = candidateMomentum;
      candidate->AddCandidate(mother);

      fOutputArray->Add(candidate);
      if(TMath::Abs(q) > 1.0E-9)
      {
        switch(TMath::Abs(candidate->PID))
        {
          case 11:
            fElectronOutputArray->Add(candidate);
            break;
          case 13:
            fMuonOutputArray->Add(candidate);
            break;
          default:
            fChargedHadronOutputArray->Add(candidate);
        }
      }
    }
    else
    {
      // 1.  initial transverse momentum p_{T0}: Part->pt
      //     initial transverse momentum direction phi_0 = -atan(p_X0/p_Y0)
      //     relativistic gamma: gamma = E/mc^2; gammam = gamma * m
      //     gyration frequency omega = q/(gamma m) fBz
      //     helix radius r = p_{T0} / (omega gamma m)

      Double_t gammam = e*1.0E9 / (c_light*c_light);      // gammam in [eV/c^2]
      Double_t omega = q * fBz / (gammam);                // omega is here in [89875518/s]
      Double_t r = pt / (q * fBz) * 1.0E9/c_light;        // in [m]

      Double_t phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

      // 2. helix axis coordinates
      Double_t x_c = x + r*TMath::Sin(phi_0);
      Double_t y_c = y - r*TMath::Cos(phi_0);
      Double_t r_c = TMath::Hypot(x_c, y_c);
      Double_t phi_c = TMath::ATan(y_c/x_c);
      if(x_c < 0.0) phi_c -= TMath::Sign(1., phi_c)*TMath::Pi();

      //Find the time of closest approach
      Double_t td = (phi_0 - TMath::ATan(-x_c/y_c))/omega;
      //Remove all the modulo pi that might have come from the atan
      Double_t pio = fabs(TMath::Pi()/omega);
      while(fabs(td) > 0.5*pio)
      {
        td -= TMath::Sign(1., td)*pio;
      }

      //Compute the coordinate of closed approach to z axis
      //if wants wtr beamline need to be changedto re-center with a traslation of the z axis
      Double_t phid = phi_0 - omega*td;
      Double_t xd = x_c - r*TMath::Sin(phid);
      Double_t yd = y_c + r*TMath::Cos(phid);
      Double_t zd = z + c_light*(pz/e)*td;

      //Compute momentum at closest approach (perigee??)
      px = pt*TMath::Cos(phid);
      py = pt*TMath::Sin(phid);

      candidateMomentum.SetPtEtaPhiE(pt, candidateMomentum.Eta(), phid, candidateMomentum.E());

      // calculate additional track parameters (correct for beamspot position)
      Double_t d0 = (xd*py - yd*px)/pt;
      Double_t dz = zd;
      // dz        = z - ((x - bsx) * px + (y - bsy) * py) / pt * (pz / pt);
      Double_t ctgTheta  = 1.0 / TMath::Tan (candidateMomentum.Theta());


      // 3. time evaluation t = TMath::Min(t_r, t_z)
      //    t_r : time to exit from the sides
      //    t_z : time to exit from the front or the back
      Double_t t = 0;
      Double_t t_z = 0;
      int sign_pz = (pz > 0.0) ? 1 : -1;
      if(pz == 0.0) t_z = 1.0E99;
      else t_z = gammam / (pz*1.0E9/c_light) * (-z + fHalfLength*sign_pz);

      if(r_c + TMath::Abs(r)  < fRadius)   // helix does not cross the cylinder sides
      {
        t = t_z;
      }
      // else
      // {
      //   Double_t asinrho = TMath::ASin((fRadius*fRadius - r_c*r_c - r*r) / (2*TMath::Abs(r)*r_c));
      //   Double_t delta = phi_0 - phi_c;
      //   if(delta <-TMath::Pi()) delta += 2*TMath::Pi();
      //   if(delta > TMath::Pi()) delta -= 2*TMath::Pi();
      //   Double_t t1 = (delta + asinrho) / omega;
      //   Double_t t2 = (delta + TMath::Pi() - asinrho) / omega;
      //   Double_t t3 = (delta + TMath::Pi() + asinrho) / omega;
      //   Double_t t4 = (delta - asinrho) / omega;
      //   Double_t t5 = (delta - TMath::Pi() - asinrho) / omega;
      //   Double_t t6 = (delta - TMath::Pi() + asinrho) / omega;
      //
      //   if(t1 < 0.0) t1 = 1.0E99;
      //   if(t2 < 0.0) t2 = 1.0E99;
      //   if(t3 < 0.0) t3 = 1.0E99;
      //   if(t4 < 0.0) t4 = 1.0E99;
      //   if(t5 < 0.0) t5 = 1.0E99;
      //   if(t6 < 0.0) t6 = 1.0E99;
      //
      //   Double_t t_ra = TMath::Min(t1, TMath::Min(t2, t3));
      //   Double_t t_rb = TMath::Min(t4, TMath::Min(t5, t6));
      //   Double_t t_r = TMath::Min(t_ra, t_rb);
      //   t = TMath::Min(t_r, t_z);
      //   // if(t<4E-9)
      //   // {
      //   //   cout << endl;
      //   //   cout << "------------------ here we are ---------------------" << endl;
      //   //   cout << "Initial pos: {" << x << ", " << y << ", " << z << "}" << endl;
      //   //   cout << "Initial Momentum: {" << pt << ", " << candidateMomentum.Eta() << ", " << phi_0 << "}" << endl;
      //   //   cout << "E: " << e << ", mass: " << candidate->Mass << ", PID: " << candidate->PID << endl;
      //   //   cout << "beta: " << candidateMomentum.P()/e << endl << endl;
      //   //
      //   //   cout << "Check luminal prop up to initial pos" << endl;
      //   //   cout << "T prod: " <<  candidatePosition.T()*1e-3/c_light << endl;
      //   //   cout << "T transverse minimum: " << TMath::Hypot(x,y)/c_light << endl;
      //   //   cout << endl;
      //   //
      //   //   cout << "Distances in [ps] from detector" << endl;
      //   //   double tfrom_endcap = (sign_pz*fHalfLength - z)/(c_light*pz/e);
      //   //   cout << "T from endcap: " << tfrom_endcap << endl;
      //   //   double tfrom_barrel = (fRadius - TMath::Hypot(x,y))/(c_light*pt/e);
      //   //   cout << "T from barrel: " << tfrom_barrel << endl << endl;
      //   //
      //   //   cout << "tz: " << t_z << endl;
      //   //   cout << "t1: " << t1 << endl;
      //   //   cout << "t2: " << t2 << endl;
      //   //   cout << "t3: " << t3 << endl;
      //   //   cout << "t4: " << t4 << endl;
      //   //   cout << "t5: " << t5 << endl;
      //   //   cout << "t6: " << t6 << endl;
      //   //   cout << "t:  " << t << endl << endl;
      //   // }
      // }
      // // 4. position in terms of x(t), y(t), z(t)
      // Double_t x_t = x_c + r * TMath::Sin(omega * t - phi_0);
      // Double_t y_t = y_c + r * TMath::Cos(omega * t - phi_0);
      // Double_t z_t = z + pz*1.0E9 / c_light / gammam * t;
      else
      {
        Double_t alpha = -(fRadius*fRadius - r*r - r_c*r_c)/(2*fabs(r)*r_c);
        alpha = fabs(TMath::ACos(alpha));

        t = td + alpha/fabs(omega);
      }
      Double_t x_t = x_c - r*TMath::Sin(phi_0 - omega*t);
      Double_t y_t = y_c + r*TMath::Cos(phi_0 - omega*t);
      Double_t z_t = z + c_light*t*pz/e;

      // if(t<4E-9)
      // {
      //   cout << endl;
      //   cout << "------------------ here we are ---------------------" << endl;
      //   cout << "Initial pos: {" << x << ", " << y << ", " << z << "}" << endl;
      //   cout << "Initial Momentum: {" << pt << ", " << candidateMomentum.Eta() << ", " << phi_0 << "}" << endl;
      //   cout << "E: " << e << ", mass: " << candidate->Mass << ", PID: " << candidate->PID << endl;
      //   cout << "beta: " << candidateMomentum.P()/e << endl << endl;
      //
      //   cout << "Check luminal prop up to initial pos" << endl;
      //   cout << "T prod: " <<  candidatePosition.T()*1e-3/c_light << endl;
      //   cout << "T transverse minimum: " << TMath::Hypot(x,y)/c_light << endl;
      //   cout << endl;
      //
      //   cout << "Distances in [ps] from detector" << endl;
      //   double tfrom_endcap = (sign_pz*fHalfLength - z)/(c_light*pz/e);
      //   cout << "T from endcap: " << tfrom_endcap << endl;
      //   double tfrom_barrel = (fRadius - TMath::Hypot(x,y))/(c_light*pt/e);
      //   cout << "T from barrel: " << tfrom_barrel << endl << endl;
      //
      //   cout << "tz: " << t_z << endl;
      //   cout << "t:  " << t << endl << endl;
      // }
      //
      //
      // if(fabs(TMath::Hypot(x_t,y_t) - fRadius)>0.01 && fabs(fabs(z_t) - fHalfLength)>0.01)
      // {
      //   cout << "-------------FLAG!" << endl;
      //   cout << TMath::Hypot(x_t,y_t)- fRadius << " - " << fabs(z_t) - fHalfLength << endl;
      //
      //     cout << "Distances in [ps] from detector" << endl;
      //     double tfrom_endcap = (sign_pz*fHalfLength - z)/(c_light*pz/e);
      //     cout << "T from endcap: " << tfrom_endcap << endl;
      //     double tfrom_barrel = (fRadius - TMath::Hypot(x,y))/(c_light*pt/e);
      //     cout << "T from barrel: " << tfrom_barrel << endl << endl;
      //
      //     cout << "tz: " << t_z << endl;
      //     cout << "t:  " << t << endl << endl;
      // }


      // compute path length for an helix
      Double_t vz = pz*1.0E9 / c_light / gammam;
      //lenght of the path from production to tracker
      Double_t l = t * TMath::Sqrt(vz*vz + r*r*omega*omega);
      //lenght of the path from closest approach to z axis to tracker`
      Double_t ld = (t-td) * TMath::Sqrt(vz*vz + r*r*omega*omega);


      // store these variables before cloning
      candidate->D0 = d0*1.0E3;
      candidate->DZ = dz*1.0E3;
      candidate->P  = candidateMomentum.P();
      candidate->PT = pt;
      //Momentum variables at closest approach
      candidate->CtgTheta = ctgTheta;
      candidate->Phi = phid;

      Candidate* mother = candidate;
      candidate = static_cast<Candidate*>(candidate->Clone());

      candidate->InitialPosition = candidatePosition;
      candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, candidatePosition.T() + t*c_light*1.0E3);

      //Momentum at closest approach
      candidate->Momentum = candidateMomentum;

      candidate->L = l*1.0E3;
      candidate->Ld = ld*1.0E3;

      candidate->Xd = xd*1.0E3;
      candidate->Yd = yd*1.0E3;
      candidate->Zd = zd*1.0E3;
      candidate->Td = candidatePosition.T() + td*c_light*1.0E3;

      candidate->AddCandidate(mother);

      fOutputArray->Add(candidate);
      switch(TMath::Abs(candidate->PID))
      {
        case 11:
          fElectronOutputArray->Add(candidate);
          break;
        case 13:
          fMuonOutputArray->Add(candidate);
          break;
        default:
          fChargedHadronOutputArray->Add(candidate);
      }

    }
  }
}

//------------------------------------------------------------------------------
