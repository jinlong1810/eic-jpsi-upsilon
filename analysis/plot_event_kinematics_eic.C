#include "normalization.h"

int plot_event_kinematics_eic( char* dfilenmae, char* suffix="new" )
{

  TFile *fin = new TFile( dfilenmae ,"OPEN" );
  //TFile *fin = new TFile("../data/sim_te_b5GeV_p50GeV_mjpsi.root","OPEN");
  TTree *T = (TTree*)fin->Get("T");

  //  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas *cscratch = new TCanvas();

  /* Define acceptance */
  char acceptance_eicsphenix[200];
  sprintf(acceptance_eicsphenix,"%s","1");

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
  T->Draw("t-tmin:W",Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event),"colz");
  c3->SetLogz(1);
  c3->Print(Form("plots/eic_events_dt_W_%s.png", suffix));

  //
  TCanvas *c4 = new TCanvas();
  T->Draw("W",Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event),"colz");
  //  c4->SetLogx(1);
  c4->SetLogy(1);
  c4->Print(Form("plots/eic_events_W_%s.png", suffix));


  return 0;
}
