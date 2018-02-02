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
      std::vector<double> z;      // z-coordinate at point of closest approach to the beamline
      std::vector<double> t;      // t-coordinate at point of closest approach to the beamline

      std::vector<double> dz_o;    //1 over square of the error of z(pca)
      std::vector<double> dt_o;    //1 over square of the error of t(pca)

      std::vector<Candidate*> tt; // a pointer to the Candidate Track
      std::vector<double> Z;      // Z[i] for DA clustering

      std::vector<double> w;     // track weight

      // std::vector<double> pt;
      // std::vector<double> pz;

      // std::vector<int> PID;

      void addItem( double new_z, double new_t, double new_dz_o, double new_dt_o, Candidate * new_tt, double new_w)
      {
        z.push_back( new_z );
        t.push_back( new_t );
        dz_o.push_back( new_dz_o );
        dt_o.push_back( new_dt_o );
        tt.push_back( new_tt );

        w.push_back( new_w );
        Z.push_back( 1.0 );
      }

      unsigned int getSize()
      {
        return z.size();
      }
    };

    struct vertex_t
    {
      std::vector<double> z;  //           z coordinate
      std::vector<double> t;  //           t coordinate
      std::vector<double> pk; //          vertex weight for "constrained" clustering

      // --- temporary numbers, used during update
      std::vector<double> ei;
      std::vector<double> swz;
      std::vector<double> swt;
      std::vector<double> se;
      std::vector<double> nuz;
      std::vector<double> nut;
      std::vector<double> szz;
      std::vector<double> stt;
      std::vector<double> szt;

      void addItem( double new_z, double new_t, double new_pk)
      {
        z.push_back( new_z);
        t.push_back( new_t);
        pk.push_back( new_pk);

        ei_cache.push_back( 0.0 );
        ei.push_back( 0.0 );
        swz.push_back( 0.0);
        swt.push_back( 0.0);
        se.push_back( 0.0);
        nuz.push_back(0.0);
        nut.push_back(0.0);
        szz.push_back(0.0);
        stt.push_back(0.0);
        szt.push_back(0.0);

        extractRaw();
      }

      void insertItem( unsigned int i, double new_z, double new_t, double new_pk   )
      {
        z.insert(z.begin() + i, new_z);
        t.insert(t.begin() + i, new_t);
        pk.insert(pk.begin() + i, new_pk);

        ei_cache.insert(ei_cache.begin() + i, 0.0 );
        ei.insert( ei.begin()  + i, 0.0 );
        swz.insert(swz.begin() + i, 0.0 );
        swt.insert(swt.begin() + i, 0.0 );
        se.insert( se.begin()  + i, 0.0 );

        nuz.insert(nuz.begin() +i, 0.0 );
        nut.insert(nut.begin() +i, 0.0 );
        szz.insert(szz.begin() + i, 0.0 );
        stt.insert(stt.begin() + i, 0.0 );
        szt.insert(szt.begin() + i, 0.0 );
      }

      void removeItem( unsigned int i )
      {
        z.erase( z.begin() + i );
        t.erase( t.begin() + i );
        pk.erase( pk.begin() + i );

        ei_cache.erase( ei_cache.begin() + i);
        ei.erase( ei.begin() + i);
        swz.erase( swz.begin() + i);
        swt.erase( swt.begin() + i);
        se.erase(se.begin() + i);

        nuz.erase(nuz.begin() + i);
        nut.erase(nut.begin() + i);
        szz.erase(szz.begin() + i);
        stt.erase(stt.begin() + i);
        szt.erase(szt.begin() + i);
      }

      unsigned int insertOrdered( double z, double t, double pk){
        // insert a new cluster according to it's z-position, return the index at which it was inserted
        unsigned int k = 0;
        for( ; k < getSize(); k++)
        {
        	if (z < z_[k]) break;
        }
        insertItem(k ,z, t, pk);
        return k;
      }

      unsigned int getSize() const
      {
        return z.size();
      }
    };

    VertexFinderDAClusterizerZT();
    ~VertexFinderDAClusterizerZT();

    void Init();
    void Process();
    void Finish();

    void clusterize(const TObjArray &tracks, TObjArray &clusters);

    std::vector< Candidate* > vertices();

    // Compute higher phase transition temperature
    double beta0(tracks_t const & tks, vertex_t &vtx);
    // --- HERE HERE ---///
    double update(double beta, tracks_t &tks, vertex_t &vtx, double rho0);

    void zorder(vertex_t & y);

    bool find_nearest(double z, double t, vertex_t & y, unsigned int & k_min, double dz, double dt);

    bool merge(vertex_t & y, double & beta);

    bool purge(vertex_t &, tracks_t &, double &, const double);

    void splitAll( vertex_t & y);

    bool split(const double beta,  tracks_t &t, vertex_t & y, double threshold = 1.);


    double get_Tc(const vertex_t & y, int k);


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

// Useful kernel functions ??? A che servono?
void kernel_calc_exp_arg(double beta, unsigned int nv, unsigned int i_tk, tracks_t &tks, vertex_t &vtx);
double kernel_add_Z(double Z_init, unsigned int nv, vertex_t &vertices);
void kernel_calc_normalization(unsigned int i_tk, tracks_t &tks, unsigned int nv, vertex_t &vtx);

#endif
