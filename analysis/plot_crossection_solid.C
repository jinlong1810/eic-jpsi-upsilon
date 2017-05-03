#include "normalization.h"

int plot_crossection_solid()
{
  TFile *fin_11GeV = new TFile("../data/sim_te_b11GeV_p0GeV_mjpsi_10k_accep.root","OPEN");
  //TFile *fin_11GeV = new TFile("../data/sim_te_b11GeV_p0GeV_mjpsi_accep.root","OPEN");
  TTree *T = (TTree*)fin_11GeV->Get("T");

  TTree *t_ref1 = new TTree();
  t_ref1->ReadFile("references/solid_pac39_xsection_data/fig16_this_proposal_mean.csv","Egamma/F:sigma/F");

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);


  /* data points from reference plot */
  t_ref1->Draw("sigma:Egamma","","same");
  TGraphErrors *g_ref1 = new TGraphErrors(t_ref1->GetEntries(), t_ref1->GetV2(), t_ref1->GetV1());
  g_ref1->SetMarkerColor(kRed);
  g_ref1->SetMarkerStyle(20);

  /* get normalization */
  Double_t overall = get_norm_solid_overall(T);

  /* define acceptance */
  char accep_normal[200];
  sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)*(accep_p_1+accep_p_2)");

  /* match factor to make cross section agree better with SOLID PAC 39 document- need to understand this and get rid of it! */
  Double_t match_factor = 1./4000.;
  Double_t overall_match = overall * match_factor;

  //  Double_t binwidth=0.2;
  TH1F* count_events = new TH1F("count_events","",45*5+1, 5, 50);
  TH1F* count_events_accep = new TH1F("count_events_accep","",45*5,5,50);
  TH1F* count_events_accep_match = new TH1F("count_events_accep_match","",45*5,5,50);

  count_events->SetMarkerStyle(20);
  count_events->GetYaxis()->SetRangeUser(0.1,2e5);

  count_events_accep->SetMarkerStyle(20);
  count_events_accep->SetMarkerColor(kGray+1);
  count_events_accep->GetYaxis()->SetTitle("#sigma (nb)");
  count_events_accep->GetXaxis()->SetTitle("E_{#gamma} (GeV)");

  count_events_accep_match->SetMarkerStyle(20);
  count_events_accep_match->SetMarkerColor(kGray+1);
  count_events_accep_match->GetYaxis()->SetRangeUser(5e-4,2e1);
  count_events_accep_match->GetYaxis()->SetTitle("#sigma (nb)");
  count_events_accep_match->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  
  cout << Form("dxs_2g*weight*weight_decay*%f",overall) << endl;

  T->Project("count_events","Keq",Form("dxs_2g*weight*weight_decay*%f",overall));
  T->Project("count_events_accep","Keq",Form("dxs_2g*weight*weight_decay*%s*%f",accep_normal,overall));
  T->Project("count_events_accep_match","Keq",Form("dxs_2g*weight*weight_decay*%s*%f",accep_normal,overall_match));

  // T->Draw("Keq >> count_events","");
  //  T->Draw("Keq >> count_events_accep","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)*(accep_p_1+accep_p_2)");

  /* create legends */
  leg1 = new TLegend(0.5,0.7,0.9,0.9);
  //leg1->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg1->AddEntry(g_ref1,"SOLID proposal","ep");
  leg1->AddEntry(count_events_accep,"this study","ep");

  leg2 = new TLegend(0.5,0.7,0.9,0.9);
  leg2->AddEntry(g_ref1,"SOLID proposal","ep");
  leg2->AddEntry(count_events_accep_match,"this study (scaled)","ep");


  TCanvas *c1 = new TCanvas();
  //  count_events->Draw("e");
  count_events_accep->Draw("e");
  g_ref1->Draw("Psame");
  leg1->Draw();
  c1->SetLogx(1);
  c1->SetLogy(1);
  c1->Print("xsection_solid.png");

  TCanvas *c2 = new TCanvas();
  count_events_accep_match->Draw("e");
  g_ref1->Draw("Psame");
  leg2->Draw();
  c2->SetLogx(1);
  c2->SetLogy(1);
  c2->Print("xsection_solid_match.png");

  return 0;
}
