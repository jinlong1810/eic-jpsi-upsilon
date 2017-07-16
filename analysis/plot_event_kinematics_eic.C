#include "normalization.h"

int plot_event_kinematics_eic( char* suffix="new" )
{

  //float Wmax=7; // J/psi production
  float Wmax=1e8; // Upsilon production

  //TString fbase("sim_te_b5GeV_p50GeV_mjpsi");
  TString fbase("sim_te_b5GeV_p50GeV_mupsilon");
  //TString fbase("sim_te_b5GeV_p50GeV_mupsilon_Wmax20");
  TString fname("../data2/");
  fname+=fbase;
  fname+=".root";
  TFile *fin = new TFile(fname,"OPEN");
  TTree *T = (TTree*)fin->Get("T");

  TString p_base("plots/");
  p_base+=fbase;
  p_base+="_";


  gStyle->SetOptStat(0);


  TCanvas *cscratch = new TCanvas();

  /* Define acceptance */
  char acceptance_eicsphenix[200];
  sprintf(acceptance_eicsphenix,"(W<%f)",Wmax);
  cout << "Acceptance cut: " << acceptance_eicsphenix << endl;

  /* Define global weight */
  double weight_global = get_norm_eic_overall(T);

  /* Define event-by-event weight */
  char weight_event[200];
  sprintf(weight_event,"%s","dxs_2g*weight*weight_decay");


  //
  TCanvas *c1 = new TCanvas();
  T->Draw("Q2:W",Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event),"colz");
  c1->SetLogz(1);
  c1->Print(Form("plots/eic_events_Q2_W_%s.png", suffix));

  //
  TCanvas *c2 = new TCanvas();
  T->Draw("tmin:W",Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event),"colz");
  c2->SetLogz(1);
  c2->Print(Form("plots/eic_events_tmin_W_%s.png", suffix));

  //
  TCanvas *c3 = new TCanvas();
  //  TH2F* h_tdiff = new TH2F("h_tdiff","",40, 10, 18.5, 40, 0, 20);
  //T->Draw("t-tmin:W >> h_tdiff",Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event),"colz");
T->Draw("t-tmin:W ",Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event),"colz");
// h_tdiff->Draw("colz");
//  c3->SetLogz(1);
  c3->Print(Form("plots/eic_events_dt_W_%s.png", suffix));

  //
  TCanvas *c4 = new TCanvas();
  T->Draw("W",Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event),"colz");
  htemp->Scale(1./htemp->Integral());
  htemp->GetYaxis()->SetRangeUser(1e-8,1);
  //  c4->SetLogx(1);
  c4->SetLogy(1);
  c4->Print(Form("plots/eic_events_W_%s.png", suffix));


  return 0;
}
