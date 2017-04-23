//#include "normalization.h"

int plot_acceptance_eic()
{
  TFile *fin = new TFile("../data/sim_te_b5GeV_p50GeV_mjpsi.root","OPEN");
  //TFile *fin = new TFile("../data/sim_te_b5330GeV_p0GeV_mjpsi.root","OPEN");
  TTree *T = (TTree*)fin->Get("T");

  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  /* get normalization */
//  Double_t overall = get_norm_solid_overall(T);
//
//  /* define acceptance */
//  char accep_3fold[200];
//  char accep_4fold[200];
//  sprintf(accep_3fold,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)");
//  //sprintf(accep_4fold,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)*(accep_p_1+accep_p_2)");
//  sprintf(accep_4fold,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)*(accep_p_1)");

  /* "The scattered electron and recoil proton will be detected by the forward angle detector,
   * while the electron-positron paor from J/Psi will be mostly detected by the large angle
   * detector" (SOLID PAC39 document) */

  TCanvas *cscratch = new TCanvas();

  //  char *content[4]={"p_e:eta_e","p_p:eta_p","p_je1:eta_je1","p_je2:eta_je2"};

  float h_pmax = 55;

  TH2F* h_p_eta_e = new TH2F("h_p_eta_e","Scattered Electron",50,-10,10,50,0,h_pmax);
  h_p_eta_e->GetXaxis()->SetTitle("Pseudorapidity");
  h_p_eta_e->GetYaxis()->SetTitle("Momentum (GeV/c)");
  h_p_eta_e->GetXaxis()->SetNdivisions(505);
  h_p_eta_e->GetYaxis()->SetNdivisions(506);
  h_p_eta_e->GetXaxis()->SetTitleSize(0.05);
  h_p_eta_e->GetYaxis()->SetTitleSize(0.05);
  h_p_eta_e->GetXaxis()->SetLabelSize(0.05);
  h_p_eta_e->GetYaxis()->SetLabelSize(0.05);
  T->Project("h_p_eta_e", "p_e:eta_e");//, Form("dxs_2g*weight*weight_decay*%s*%f",accep_4fold,overall));

  TCanvas *c18_ul = new TCanvas();
  h_p_eta_e->Draw("colz");
  c18_ul->SetLogz(1);
  c18_ul->SetRightMargin(0.2);
  c18_ul->Print("new_eic_p_eta_e.png");

  cscratch->cd();



  TH2F* h_p_eta_je1 = new TH2F("h_p_eta_je1","Decay Electron",50,-10,10,50,0,h_pmax);
  h_p_eta_je1->GetXaxis()->SetTitle("Pseudorapidity");
  h_p_eta_je1->GetYaxis()->SetTitle("Momentum (GeV/c)");
  h_p_eta_je1->GetXaxis()->SetNdivisions(505);
  h_p_eta_je1->GetYaxis()->SetNdivisions(506);
  h_p_eta_je1->GetXaxis()->SetTitleSize(0.05);
  h_p_eta_je1->GetYaxis()->SetTitleSize(0.05);
  h_p_eta_je1->GetXaxis()->SetLabelSize(0.05);
  h_p_eta_je1->GetYaxis()->SetLabelSize(0.05);
  T->Project("h_p_eta_je1", "p_je1:eta_je1");//, Form("dxs_2g*weight*weight_decay*%s*%f",accep_4fold,overall));

  TCanvas *c18_ur = new TCanvas();
  h_p_eta_je1->Draw("colz");
  c18_ur->SetLogz(1);
  c18_ur->SetRightMargin(0.2);
  c18_ur->Print("new_eic_p_eta_je1.png");

  cscratch->cd();


  TH2F* h_p_eta_je2 = new TH2F("h_p_eta_je2","Decay Positron",50,-10,10,50,0,h_pmax);
  h_p_eta_je2->GetXaxis()->SetTitle("Pseudorapidity");
  h_p_eta_je2->GetYaxis()->SetTitle("Momentum (GeV/c)");
  h_p_eta_je2->GetXaxis()->SetNdivisions(505);
  h_p_eta_je2->GetYaxis()->SetNdivisions(506);
  h_p_eta_je2->GetXaxis()->SetTitleSize(0.05);
  h_p_eta_je2->GetYaxis()->SetTitleSize(0.05);
  h_p_eta_je2->GetXaxis()->SetLabelSize(0.05);
  h_p_eta_je2->GetYaxis()->SetLabelSize(0.05);
  T->Project("h_p_eta_je2", "p_je2:eta_je2");//, Form("dxs_2g*weight*weight_decay*%s*%f",accep_4fold,overall));

  TCanvas *c18_ll = new TCanvas();
  h_p_eta_je2->Draw("colz");
  c18_ll->SetLogz(1);
  c18_ll->SetRightMargin(0.2);
  c18_ll->Print("new_eic_p_eta_je2.png");

  cscratch->cd();


  TH2F* h_p_eta_p = new TH2F("h_p_eta_p","Recoil Proton",50,-10,10,50,0,h_pmax);
  h_p_eta_p->GetXaxis()->SetTitle("Pseudorapidity");
  h_p_eta_p->GetYaxis()->SetTitle("Momentum (GeV/c)");
  h_p_eta_p->GetXaxis()->SetNdivisions(505);
  h_p_eta_p->GetYaxis()->SetNdivisions(506);
  h_p_eta_p->GetXaxis()->SetTitleSize(0.05);
  h_p_eta_p->GetYaxis()->SetTitleSize(0.05);
  h_p_eta_p->GetXaxis()->SetLabelSize(0.05);
  h_p_eta_p->GetYaxis()->SetLabelSize(0.05);
  T->Project("h_p_eta_p", "p_p:eta_p");//, Form("dxs_2g*weight*weight_decay*%s*%f",accep_4fold,overall));

  TCanvas *c18_lr = new TCanvas();
  h_p_eta_p->Draw("colz");
  c18_lr->SetLogz(1);
  c18_lr->SetRightMargin(0.2);
  c18_lr->Print("new_eic_p_eta_p.png");

  cscratch->cd();

  TH2F* h_p_eta_jpsi = new TH2F("h_p_eta_jpsi","JPsi",50,-10,10,50,0,h_pmax);
  h_p_eta_jpsi->GetXaxis()->SetTitle("Pseudorapidity");
  h_p_eta_jpsi->GetYaxis()->SetTitle("Momentum (GeV/c)");
  h_p_eta_jpsi->GetXaxis()->SetNdivisions(505);
  h_p_eta_jpsi->GetYaxis()->SetNdivisions(506);
  h_p_eta_jpsi->GetXaxis()->SetTitleSize(0.05);
  h_p_eta_jpsi->GetYaxis()->SetTitleSize(0.05);
  h_p_eta_jpsi->GetXaxis()->SetLabelSize(0.05);
  h_p_eta_jpsi->GetYaxis()->SetLabelSize(0.05);
  T->Project("h_p_eta_jpsi", "p_jpsi:eta_jpsi");//, Form("dxs_2g*weight*weight_decay*%s*%f",accep_4fold,overall));

  TCanvas *cx = new TCanvas();
  h_p_eta_jpsi->Draw("colz");
  cx->SetLogz(1);
  cx->SetRightMargin(0.2);
  cx->Print("new_eic_p_eta_jpsi.png");

  cscratch->cd();

  TH2F* h_p_e_p = new TH2F("h_p_e_p","Momentum Scattered Electron vs Proton",50,0,h_pmax,50,0,h_pmax);
  h_p_e_p->GetXaxis()->SetTitle("Scattered Proton Momentum (GeV/c)");
  h_p_e_p->GetYaxis()->SetTitle("Scattered Electron Momentum (GeV/c)");
  h_p_e_p->GetXaxis()->SetNdivisions(505);
  h_p_e_p->GetYaxis()->SetNdivisions(506);
  h_p_e_p->GetXaxis()->SetTitleSize(0.05);
  h_p_e_p->GetYaxis()->SetTitleSize(0.05);
  T->Project("h_p_e_p","p_e:p_p","");

  TCanvas *ccorr1 = new TCanvas();
  h_p_e_p->Draw("colz");
  ccorr1->SetLogz(1);
  ccorr1->Print("new_eic_corr1.png");

  return 0;
}
