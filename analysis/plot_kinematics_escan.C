int plot_kinematics_escan()
{
  TFile *fin_11GeV = new TFile("../data/sim_te_b11GeV_p0GeV_mjpsi_accep.root","OPEN");
  TTree *T_11GeV = (TTree*)fin_11GeV->Get("T");

  TFile *fin_50GeV = new TFile("../data/sim_te_b50GeV_p0GeV_mjpsi.root","OPEN");
  TTree *T_50GeV = (TTree*)fin_50GeV->Get("T");

  TFile *fin_500GeV = new TFile("../data/sim_te_b500GeV_p0GeV_mjpsi.root","OPEN");
  TTree *T_500GeV = (TTree*)fin_500GeV->Get("T");

  TFile *fin_10TeV = new TFile("../data/sim_te_b10TeV_p0GeV_mjpsi.root","OPEN");
  TTree *T_10TeV = (TTree*)fin_10TeV->Get("T");

  TFile *fin_5x50GeV = new TFile("../data/sim_te_b5GeV_p50GeV_mjpsi.root","OPEN");
  TTree *T_5x50GeV = (TTree*)fin_5x50GeV->Get("T");

  TFile *fin_500GeV_ups = new TFile("../data/sim_te_b500GeV_p0GeV_mupsilon.root","OPEN");
  TTree *T_500GeV_ups = (TTree*)fin_500GeV_ups->Get("T");

  //  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  /* "The scattered electron and recoil proton will be detected by the forward angle detector,
   * while the electron-positron paor from J/Psi will be mostly detected by the large angle
   * detector" (SOLID PAC39 document) */

  TCanvas *cscratch = new TCanvas();

  /* 11 GeV */
  TH1F* h_W_11GeV = new TH1F("h_W_11GeV","11x0 GeV",100,4,4.7);
  h_W_11GeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_11GeV->GetYaxis()->SetTitle("#entries");
  h_W_11GeV->GetXaxis()->SetNdivisions(504);
  h_W_11GeV->GetYaxis()->SetNdivisions(505);

  TH1F* h_W_11GeV_3fold = (TH1F*) h_W_11GeV->Clone("h_W_11GeV_3fold");
  h_W_11GeV_3fold->SetFillColor(kGreen+1);

  T_11GeV->Draw("W >> h_W_11GeV");
  T_11GeV->Draw("W >> h_W_11GeV_3fold","accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","colz");

  TCanvas *c_11GeV = new TCanvas();
  c_11GeV->SetLogy(1);
  h_W_11GeV->Draw();
  h_W_11GeV_3fold->Draw("same");
  c_11GeV->Print("new_escan_W_11GeV.png");

  cscratch->cd();

  /* 50 GeV */
  TH1F* h_W_50GeV = new TH1F("h_W_50GeV","50x0 GeV",100,4,10);
  h_W_50GeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_50GeV->GetYaxis()->SetTitle("#entries");
  h_W_50GeV->GetXaxis()->SetNdivisions(504);
  h_W_50GeV->GetYaxis()->SetNdivisions(505);

  T_50GeV->Draw("W >> h_W_50GeV");

  TCanvas *c_50GeV = new TCanvas();
  c_50GeV->SetLogy(1);
  h_W_50GeV->Draw();
  c_50GeV->Print("new_escan_W_50GeV.png");

  cscratch->cd();


  /* 500 GeV */
  TH1F* h_W_500GeV = new TH1F("h_W_500GeV","500x0 GeV",100,4,35);
  h_W_500GeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_500GeV->GetYaxis()->SetTitle("#entries");
  h_W_500GeV->GetXaxis()->SetNdivisions(504);
  h_W_500GeV->GetYaxis()->SetNdivisions(505);

  T_500GeV->Draw("W >> h_W_500GeV");

  TCanvas *c_500GeV = new TCanvas();
  c_500GeV->SetLogy(1);
  h_W_500GeV->Draw();
  c_500GeV->Print("new_escan_W_500GeV.png");

  cscratch->cd();


  /* 10 TeV */
  TH1F* h_W_10TeV = new TH1F("h_W_10TeV","10x0 TeV",100,4,150);
  h_W_10TeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_10TeV->GetYaxis()->SetTitle("#entries");
  h_W_10TeV->GetXaxis()->SetNdivisions(504);
  h_W_10TeV->GetYaxis()->SetNdivisions(505);

  T_10TeV->Draw("W >> h_W_10TeV");

  TCanvas *c_10TeV = new TCanvas();
  c_10TeV->SetLogy(1);
  h_W_10TeV->Draw();
  c_10TeV->Print("new_escan_W_10TeV.png");

  cscratch->cd();


  /* 5x50 GeV */
  TH1F* h_W_5x50GeV = new TH1F("h_W_5x50GeV","5x50 GeV",100,0,35);
  h_W_5x50GeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_5x50GeV->GetYaxis()->SetTitle("#entries");
  h_W_5x50GeV->GetXaxis()->SetNdivisions(504);
  h_W_5x50GeV->GetYaxis()->SetNdivisions(505);
  h_W_5x50GeV->GetYaxis()->SetRangeUser(1,2e5);

  T_5x50GeV->Draw("W >> h_W_5x50GeV");
  
  TH1F* h_W_5x50GeV_ref1 = (TH1F*)h_W_500GeV->Clone("h_W_5x50GeV_ref1");
  h_W_5x50GeV_ref1->SetLineColor(kGray+1);

  TCanvas *c_5x50GeV = new TCanvas();
  c_5x50GeV->SetLogy(1);
  h_W_5x50GeV->Draw();
  h_W_5x50GeV_ref1->Draw("sames");
  c_5x50GeV->Print("new_escan_W_5x50GeV.png");

  cscratch->cd();


  /* 500 GeV Upsilon */
  TH1F* h_W_500GeV_ups = new TH1F("h_W_500GeV_ups","500x0 GeV (Upsilon)",100,4,35);
  h_W_500GeV_ups->GetXaxis()->SetTitle("W (GeV)");
  h_W_500GeV_ups->GetYaxis()->SetTitle("#entries");
  h_W_500GeV_ups->GetXaxis()->SetNdivisions(504);
  h_W_500GeV_ups->GetYaxis()->SetNdivisions(505);

  T_500GeV_ups->Draw("W >> h_W_500GeV_ups");

  TH1F* h_W_500GeV_ref1 = (TH1F*)h_W_500GeV->Clone("h_W_500GeV_ref1");
  h_W_500GeV_ref1->SetLineColor(kGray+1);

  TCanvas *c_500GeV_ups = new TCanvas();
  c_500GeV_ups->SetLogy(1);
  h_W_500GeV_ups->Draw();
  h_W_500GeV_ref1->Draw("sames");
  c_500GeV_ups->Print("new_escan_W_500GeV_ups.png");

  cscratch->cd();


  return 0;
}
