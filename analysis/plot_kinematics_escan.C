int plot_kinematics_escan()
{
  TFile *fin_11GeV = new TFile("../data/solid_e_11GeV_1M_accep.root","OPEN");
  TTree *T_11GeV = (TTree*)fin_11GeV->Get("T");

  TFile *fin_50GeV = new TFile("../data/solid_e_50GeV_1M.root","OPEN");
  TTree *T_50GeV = (TTree*)fin_50GeV->Get("T");

  TFile *fin_500GeV = new TFile("../data/solid_e_500GeV_1M.root","OPEN");
  TTree *T_500GeV = (TTree*)fin_500GeV->Get("T");

  TFile *fin_10TeV = new TFile("../data/solid_e_10TeV_1M.root","OPEN");
  TTree *T_10TeV = (TTree*)fin_10TeV->Get("T");

  //  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  /* "The scattered electron and recoil proton will be detected by the forward angle detector,
   * while the electron-positron paor from J/Psi will be mostly detected by the large angle
   * detector" (SOLID PAC39 document) */

  TCanvas *cscratch = new TCanvas();

  /* 11 GeV */
  T_11GeV->Draw("W >> h_W_11GeV");
  T_11GeV->Draw("W >> h_W_11GeV_3fold","accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","colz");

  h_W_11GeV->SetTitle("");
  h_W_11GeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_11GeV->GetYaxis()->SetTitle("#entries");
  h_W_11GeV->GetXaxis()->SetNdivisions(504);
  h_W_11GeV->GetYaxis()->SetNdivisions(505);

  h_W_11GeV_3fold->SetFillColor(kGreen+1);

  TCanvas *c_11GeV = new TCanvas();
  c_11GeV->SetLogy(1);
  h_W_11GeV->Draw();
  h_W_11GeV_3fold->Draw("same");
  c_11GeV->Print("new_escan_W_11GeV.png");
  cscratch->cd();



  /* 50 GeV */
  T_50GeV->Draw("W >> h_W_50GeV");
  h_W_50GeV->SetTitle("");
  h_W_50GeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_50GeV->GetYaxis()->SetTitle("#entries");
  h_W_50GeV->GetXaxis()->SetNdivisions(504);
  h_W_50GeV->GetYaxis()->SetNdivisions(505);

  TCanvas *c_50GeV = new TCanvas();
  c_50GeV->SetLogy(1);
  h_W_50GeV->Draw();
  c_50GeV->Print("new_escan_W_50GeV.png");
  cscratch->cd();



  /* 500 GeV */
  T_500GeV->Draw("W >> h_W_500GeV");
  h_W_500GeV->SetTitle("");
  h_W_500GeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_500GeV->GetYaxis()->SetTitle("#entries");
  h_W_500GeV->GetXaxis()->SetNdivisions(504);
  h_W_500GeV->GetYaxis()->SetNdivisions(505);

  TCanvas *c_500GeV = new TCanvas();
  c_500GeV->SetLogy(1);
  h_W_500GeV->Draw();
  c_500GeV->Print("new_escan_W_500GeV.png");
  cscratch->cd();



  /* 10 TeV */
  T_10TeV->Draw("W >> h_W_10TeV");
  h_W_10TeV->SetTitle("");
  h_W_10TeV->GetXaxis()->SetTitle("W (GeV)");
  h_W_10TeV->GetYaxis()->SetTitle("#entries");
  h_W_10TeV->GetXaxis()->SetNdivisions(504);
  h_W_10TeV->GetYaxis()->SetNdivisions(505);

  TCanvas *c_10TeV = new TCanvas();
  c_10TeV->SetLogy(1);
  h_W_10TeV->Draw();
  c_10TeV->Print("new_escan_W_10TeV.png");
  cscratch->cd();


  return 0;
}
