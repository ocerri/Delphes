#ifndef VertexFinderDAClusterizerZT_h
#define VertexFinderDAClusterizerZT_h

/** \class VertexFinderDAClusterizerZT
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  Author O. Cerri
 *
 */


#include "classes/DelphesModule.h"

#include <vector>

class TObjArray;//???? A che servono? Provare a levarli
class TIterator;
class Candidate;

class VertexFinderDAClusterizerZT: public DelphesModule
{
  public:

    struct tracks_t
    {
      std::vector<double> z;      // z-coordinate at point of closest approach to the beamline [mm]
      std::vector<double> t;      // t-coordinate at point of closest approach to the beamline  [ps]

      std::vector<double> dz2_o;    //1 over the squared error of z(pca)
      std::vector<double> dz_o;
      std::vector<double> dt2_o;    //1 over the squared error of t(pca)
      std::vector<double> dt_o;

      std::vector<Candidate*> tt; // a pointer to the Candidate Track
      std::vector<double> Z;      // Z[i] for DA clustering

      std::vector<double> w;     // track weight

      // std::vector<double> pt;
      // std::vector<double> pz;

      // std::vector<int> PID;

      double sum_w = 0;
      double sum_dz2_o = 0;
      double sum_dt2_o = 0;

      void addItem( double new_z, double new_t, double new_dz2_o, double new_dt2_o, Candidate * new_tt, double new_w)
      {
        z.push_back( new_z );
        t.push_back( new_t );
        dz2_o.push_back( new_dz2_o );
        dz_o.push_back( sqrt(new_dz2_o) );
        dt2_o.push_back( new_dt2_o );
        dt_o.push_back( sqrt(new_dt2_o) );
        tt.push_back( new_tt );

        w.push_back( new_w );
        Z.push_back( 1.0 );
      }

      unsigned int getSize()
      {
        return z.size();
      }

      void NormalizeWeights()
      {
        if(sum_w <= 0)
        {
          throw std::invalid_argument( "Non positive sum of tracks weights" );
        }
        double new_sum = 0; //In principle might be useless but prevents from numerical errors
        for(unsigned int i = 0; i < z.size(); i++)
        {
          w[i] /= sum_w;
          new_sum += w[i];
        }
        sum_w = new_sum;
      }
    };

    struct vertex_t
    {
      std::vector<double> z;  //           z coordinate
      std::vector<double> t;  //           t coordinate
      std::vector<double> pk; //          vertex weight for "constrained" clustering

      // Elements of the covariance matrix of the posterior distribution
      std::vector<double> szz, stt, stz;
      std::vector<double> beta_c;

      void addItem( double new_z, double new_t, double new_pk)
      {
        z.push_back( new_z);
        t.push_back( new_t);
        pk.push_back( new_pk);

        szz.push_back(0);
        stt.push_back(0);
        szt.push_back(0);
        beta_c.push_back(0)
      }

      void insertItem( unsigned int i, double new_z, double new_t, double new_pk   )
      {
        z.insert(z.begin() + i, new_z);
        t.insert(t.begin() + i, new_t);
        pk.insert(pk.begin() + i, new_pk);

        szz.insert(szz.begin() + i, 0 );
        stt.insert(stt.begin() + i, 0 );
        szt.insert(szt.begin() + i, 0 );
        beta_c.insert(beta_c.begin() + i, 0 );
      }

      void removeItem( unsigned int i )
      {
        z.erase( z.begin() + i );
        t.erase( t.begin() + i );
        pk.erase( pk.begin() + i );

        szz.erase(szz.begin() + i);
        stt.erase(stt.begin() + i);
        szt.erase(szt.begin() + i);
        beta_c.erase(beta_c.begin() + i);
      }

      unsigned int getSize() const
      {
        return z.size();
      }

      double DistanceSquare(unsigned int k1, unsigned int k1)
      {
        double dz = (z[k1] - z[k2])/fVertexSpaceSize;
        double dt = (t[k1] - t[k2])/fVertexTimeSize;

        return dz*dz + dt*dt;
      }

      pair<double, int> ComputeAllBeta_c()
      {
        unsigned int nv = getSize();
        unsigned int k_min = 0;
        double beta_c_min = 0;
        for(unsigned int k = 0; k < nv; k++)
        {
          if(szz[k] == 0 || stt[k] == 0)
          {
            throw std::invalid_argument( "Attempting to compute beta c for uninitialized vertex" );
          }

          double aux = (stt[k] - szz[k])*(stt[k] - szz[k]) + 4*szt[k]*szt[k];
          aux = 1. / (stt[k] + szz[k] + sqrt(aux));

          if(aux > beta_c_min)
          {
            beta_c_min = aux;
            k_min = k;
          }

          beta_c[k] = aux;
        }

        pair<double, unsigned int> out(beta_c_min, k_min);
        return out;
      }
    };

    VertexFinderDAClusterizerZT();
    ~VertexFinderDAClusterizerZT();

    void Init();
    void Process();
    void Finish();

    void clusterize(const TObjArray &tracks, TObjArray &clusters);

    // Define the distance metric between tracks and vertices
    double Energy(double t_z, double v_z, double dz2_o, double t_t, double v_t, double dt2_o);

    // Fill the tracks structure from the input array
    void fill(tracks_t &tks);

    // Compute higher phase transition temperature
    double beta0(tracks_t const & tks, vertex_t &vtx);

    // Compute the new vertexes position and mass
    double update(double beta, tracks_t &tks, vertex_t &vtx, double rho0);

    // If a vertex has beta_c lower than beta, split it
    bool split(const double beta,  vertex_t & vtx, const double epsilon = 2.);

    // Merge vertexes closer than declared dimensions
    bool merge(vertex_t & vtx, double d2_merge);



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
