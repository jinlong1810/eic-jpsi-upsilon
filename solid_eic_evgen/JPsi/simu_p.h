#include <TSystem.h>
#include <string>


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

  Float_t t0lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);
  Float_t t1lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);
  Double_t fun_2g(Double_t x, Double_t t, Double_t M);
  Double_t fun_23g(Double_t x, Double_t t, Double_t M);

  Int_t seed; //random number seed

  Int_t nevents; //number of events

  Double_t Ebeam_lab; // Electron Beam Energy in labortory frame
  Double_t Etarget_lab; // Proton Beam Energy in labortory frame

  TString output_root_file;

  TString meson_type;

  bool Is_e;
  bool Is_g;
  
  Double_t Gbeam_min;
};
