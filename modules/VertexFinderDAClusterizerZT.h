#ifndef VertexFinderDAClusterizerZT_h
#define VertexFinderDAClusterizerZT_h

/** \class VertexFinderDAClusterizerZT
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  \authors M. Selvaggi, L. Gray
 *
 */


#include "classes/DelphesModule.h"

#include <vector>

class TObjArray;
class TIterator;
class Candidate;

class VertexFinderDAClusterizerZT: public DelphesModule
{
public:

  VertexFinderDAClusterizerZT();
  ~VertexFinderDAClusterizerZT();

  void Init();
  void Process();
  void Finish();

  void clusterize(const TObjArray &tracks, TObjArray &clusters);
  std::vector< Candidate* > vertices();

private:

  Bool_t fVerbose;

  UInt_t fMaxIterations;
  Double_t fMinTrackWeight;
  Bool_t fUseTc;

  Float_t fBetaMax;
  Float_t fBetaStop;
  Float_t fBetaPurge;


  Float_t fVertexSpaceSize;
  Float_t fVertexTimeSize;

  Double_t fCoolingFactor;

  Double_t fDzCutOff;
  Double_t fD0CutOff;
  Double_t fDtCutOff; // for when the beamspot has time
  Double_t fMinPT;

  Double_t fUniqueTrkWeight;
  Double_t fDzMerge;
  Double_t fDtMerge;

  TObjArray *fInputArray;
  TIterator *fItInputArray;

  TObjArray *fOutputArray;
  TObjArray *fVertexOutputArray;

  ClassDef(VertexFinderDAClusterizerZT, 1)
};

#endif
