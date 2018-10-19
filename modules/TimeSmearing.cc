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


/** \class TimeSmearing
 *
 *  Performs transverse momentum resolution smearing.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TimeSmearing.h"

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

TimeSmearing::TimeSmearing() :
fItInputArray(0)
{
}

//------------------------------------------------------------------------------

TimeSmearing::~TimeSmearing()
{
}

//------------------------------------------------------------------------------

void TimeSmearing::Init()
{
  // read resolution formula

  fTimeResolution = GetDouble("TimeResolution", 3.0E-11);
  fEtaMax = GetDouble("EtaMax", 9999.);
  // import input array

  fInputArray = ImportArray(GetString("InputArray", "MuonMomentumSmearing/muons"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "muons"));
}

//------------------------------------------------------------------------------

void TimeSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TimeSmearing::Process()
{
  Candidate *candidate, *mother;
  Double_t ti, tf_smeared, tf;
  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    ti = candidate->InitialPosition.T()*1.0E-3/c_light;
    tf = candidate->Position.T()*1.0E-3/c_light;

    // apply smearing formula
    if(fabs(candidate->Position.Eta())<fEtaMax)
    {
      tf_smeared = tf + fTimeResolution*gRandom->Gaus(0, 1);
    }
    else continue;

    // double beta_particle = candidate->Momentum.P()/candidate->Momentum.E();
    // ti = tf_smeared - candidate->Ld*1.0E-3/(c_light*beta_particle);

    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->AddCandidate(mother);
    candidate->InitialPosition.SetT((100+ti)*1.0E3*c_light);
    candidate->Position.SetT(tf_smeared*1.0E3*c_light);
    candidate->ErrorT = fTimeResolution*1.0E3*c_light;

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
