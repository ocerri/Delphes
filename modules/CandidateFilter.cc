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

/*  \class CandidateFilter
 *
 *  Removes particles with given features
 *
 *  \author Olmo Cerri
 *
 */

#include "modules/CandidateFilter.h"

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

CandidateFilter::CandidateFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

CandidateFilter::~CandidateFilter()
{
}

//------------------------------------------------------------------------------

void CandidateFilter::Init()
{

  ExRootConfParam param;

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredCandidates"));

  fPtMin = GetDouble("PtMin", 20);
  fMassMin = GetDouble("MassMin", 0.20);
}

//------------------------------------------------------------------------------

void CandidateFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void CandidateFilter::Process()
{
  Candidate *candidate;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    Bool_t pass = candidate->Momentum.Pt() > fPtMin;
    pass *= candidate->ClusterIndex >= 0;
    pass *= candidate->Mass > fMassMin;

    if(pass) fOutputArray->Add(candidate);
  }
}
