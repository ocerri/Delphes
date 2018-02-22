#ifndef AdaptiveVertexFitting4D_h
#define AdaptiveVertexFitting4D_h

/** \class AdaptiveVertexFitting4D
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  Author O. Cerri
 *
 */

#include "classes/DelphesModule.h"

#include <vector>
#include <iostream>

#ifndef _c_light__
#define _c_light__
namespace AVF_4D
{
  static const Double_t c_light   = 2.99792458e-01; //[mm/ps]
}
#endif

using namespace std;
using namespace AVF_4D;

class TObjArray;
class TIterator;
class Candidate;

void NegativeLogLikelihood(Int_t&, Double_t *, Double_t&f, Double_t*par, Int_t );

class AdaptiveVertexFitting4D: public DelphesModule
{
  public:

    class tracks_t : public TObject
    {
      public:
        std::vector<double> t_out;      // t-coordinate at point of of intersection with the timing layer  [ps]
        std::vector<double> s2t_out;      // sigma square of t_out [ps]

        std::vector<double> z_out;      // z-coordinate at point of intersection with the timing layer [mm]
        std::vector<double> s2z_out;      // sigma square of z_out [mm]

        std::vector<double> Pt;    // Transverse Momentum
        std::vector<double> s2_Pt;    // sigma square of transverse Momentum

        std::vector<double> ctg_theta;    // CtgTheta = PZ/PT
        std::vector<double> E;    // Energy
        std::vector<double> beta_z;    // Beta along z
        std::vector<double> M;

        // parameters from the last iteration
        // std::vector<double> wz;     // track weight
        double sum_wz;     // track weight


        std::vector<Candidate*> tt; // a pointer to the Candidate Track

        // parameters from the last iteration

        tracks_t(){}
        ~tracks_t(){}

        void dump(unsigned int i)
        {
          cout << "Track " << i << endl;
          cout <<
          Form("t_out: %1.2e  ", t_out[i]) <<
          Form("st_out: %1.2e  ", sqrt(s2t_out[i])) <<
          Form("z_out: %1.2e  ", z_out[i]) <<
          Form("sz_out: %1.2e  ", sqrt(s2z_out[i])) <<
          Form("Pt: %1.2e  ", Pt[i]) <<
          Form("s2_Pt: %1.2e  ", sqrt(s2_Pt[i])) <<
          Form("ctg_theta: %1.2e  ", ctg_theta[i]) <<
          Form("E: %1.2e  ", E[i]) <<
          Form("beta_z: %1.2e  ", beta_z[i]) <<
          Form("M: %1.2e  ", M[i]) << endl << endl;
        }

        unsigned int getSize()
        {
          return z_out.size();
        }

        double z(unsigned int i, double t)
        {
          return z_out[i] + c_light * beta_z[i] * (t - t_out[i]);
        }

        double s2z(unsigned int i, double t)
        {
          // Solve! It overestimate the error
          double out = 0;
          out += s2z_out[i];

          double aux = c_light * beta_z[i];
          out += aux * aux * s2t_out[i];

          aux = c_light * (t - t_out[i]) * ctg_theta[i] * ( 1 - Pt[i]*Pt[i]*(1+ctg_theta[i]*ctg_theta[i])/(E[i]*E[i]) ) / E[i];
          out += aux*aux*s2_Pt[i];

          return out;
        }
    };

    class vertex_t : public TObject
    {
      public:
        double x,y,z,t;
        double sigma_t = 30;
        double sigma_z = 1;

        double chi2 = 0;

        vertex_t( double t_new, double z_new)
        {
          x = 0;
          y = 0;
          t = t_new;
          z = z_new;
        };

        ~vertex_t(){};
    };

    AdaptiveVertexFitting4D();
    ~AdaptiveVertexFitting4D();

    void Init();
    void Process();
    void Finish();


    // Fill the tracks structure from the vertex
    void fill(Candidate* vtx, tracks_t &tks);

    // Fit the vertex and update the vertex position
    void fit(double beta, tracks_t &tks, vertex_t &vtx);

  private:

    UInt_t fVerbose;

    UInt_t fMaxIterations;

    Double_t fBetaMax;

    Double_t fCoolingFactor;
    Double_t fChi2_0;


    TObjArray *fInputArray;
    TIterator *fItInputArray;

    TObjArray *fOutputArray;
    TObjArray *fVertexOutputArray;

    ClassDef(AdaptiveVertexFitting4D, 1)
};

#endif
