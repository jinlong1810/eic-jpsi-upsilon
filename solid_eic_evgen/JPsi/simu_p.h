#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>

class TFile;
class TTree;

class Simulator {
 
 public:
  Simulator();

  int run();  

  void set_seed( int newseed ) { seed = newseed; }

  void set_number_events( Int_t nnew ) { nevents = nnew; }

  void set_lepton_energy( Double_t enew ) { Ebeam_lab = enew; }

  void set_hadron_energy( Double_t enew ) { Etarget_lab = enew; }

  void set_output_file( TString newfile ) { output_root_file = newfile; }

  int set_meson_type( TString newtype );

  int set_process_type( TString type );

 protected:

  int init();
  int create_output_file();
  void reset_event_variables();
  int process_event();
  int end();

  Float_t t0lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);
  Float_t t1lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);
  Double_t fun_2g(Double_t x, Double_t t, Double_t M);
  Double_t fun_23g(Double_t x, Double_t t, Double_t M);

  TFile *_file_out;
  TTree *_tree_out;

  Int_t seed; //random number seed

  Int_t nevents; //number of events

  Double_t Ebeam_lab; // Electron Beam Energy in labortory frame
  Double_t Etarget_lab; // Proton Beam Energy in labortory frame

  TString output_root_file;

  TString meson_type;

  bool Is_e;
  bool Is_g;
  
  Double_t Gbeam_min;

  /* counter of actual event */
  Int_t neve;

  /* Event-wise parameters */

  /* Photon energy and flux */
  Double_t Gbeam;
  Double_t Gflux;

  /* global event paraemters */
  Double_t Q2;
  Double_t t;

  /* invariant masses */
  Double_t minv;
  Double_t minv_prest;
  Double_t minv_beam;

  /* final state particle properties in laboratory frame */
  Double_t p_e;
  Double_t pt_e;
  Double_t theta_e;
  Double_t phi_e;
  Double_t eta_e;
  Double_t p_p;
  Double_t pt_p;
  Double_t theta_p;
  Double_t phi_p;
  Double_t eta_p;
  Double_t p_jpsi;
  Double_t pt_jpsi;
  Double_t theta_jpsi;
  Double_t phi_jpsi;
  Double_t eta_jpsi;

  Double_t p_je1;
  Double_t pt_je1;
  Double_t theta_je1;
  Double_t phi_je1;
  Double_t eta_je1;
  Double_t p_je2;
  Double_t pt_je2;
  Double_t theta_je2;
  Double_t phi_je2;
  Double_t eta_je2;

  /* Weights */
  Double_t weight_decay;
  Double_t weight;
    
  /* calculation formular */
  Double_t Gamma;
  Double_t epsilon;
  Double_t Keq;
  Double_t W;
  Double_t q;
  Double_t theta_q;
  Double_t J;
  Double_t R;
  Double_t theta_cm;
  Double_t phi_cm;
  Double_t r;

  Double_t dxs;
  Double_t dxs_2g;
  Double_t dxs_23g;

  Double_t tmin;
  Double_t tmax;

  Double_t phasespace;

};
