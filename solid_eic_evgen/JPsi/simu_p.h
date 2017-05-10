#include <TSystem.h>
#include <string>


class Simulator {
 
 public:
  Simulator();

  int run();  

  Float_t t0lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);
  Float_t t1lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);
  Double_t fun_2g(Double_t x, Double_t t, Double_t M);
  Double_t fun_23g(Double_t x, Double_t t, Double_t M);

  // protected:
  Int_t seed;

  Int_t nevents; //number of events

  Double_t Ebeam_lab; // Electron Beam Energy in labortory frame
  Double_t Etarget_lab; // Proton Beam Energy in labortory frame

  TString output_root_file;

  TString meson_type;
  TString acceptance_root_file;

  bool Is_e;
  bool Is_g;
  
  Double_t Gbeam_min;
};
