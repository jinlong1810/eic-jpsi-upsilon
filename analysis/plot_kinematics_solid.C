int plot_kinematics_solid()
{
  TFile *fin = new TFile("../data/solid_e_11GeV_1M_accep.root","OPEN");
  TTree *T = (TTree*)fin->Get("T");

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  /* "The scattered electron and recoil proton will be detected by the forward angle detector,
   * while the electron-positron paor from J/Psi will be mostly detected by the large angle
   * detector" (SOLID PAC39 document) */

  TCanvas *cscratch = new TCanvas();

  /* Figure 6 upper-left: 3-fold coincidence */
  TH2F* h_W_Q2 = new TH2F("h_W_Q2","",50,4.0,4.45,50,0.3,1.6);
  h_W_Q2->GetXaxis()->SetTitle("W (GeV)");
  h_W_Q2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  h_W_Q2->GetXaxis()->SetNdivisions(504);
  h_W_Q2->GetYaxis()->SetNdivisions(505);
  h_W_Q2->GetXaxis()->SetTitleSize(0.05);
  h_W_Q2->GetYaxis()->SetTitleSize(0.05);
  h_W_Q2->GetXaxis()->SetLabelSize(0.05);
  h_W_Q2->GetYaxis()->SetLabelSize(0.05);
  T->Draw("Q2:W >> h_W_Q2","accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","colz");

  TCanvas *c6_ul = new TCanvas();
  h_W_Q2->Draw("colz");
  c6_ul->Print("new_solid_W_Q2.png");

  cscratch->cd();


  /* Figure 6 lower-left: 3-fold coincidence */
  TH2F* h_W_dt_3fold = new TH2F("h_W_dt_3fold","",50,4.0,4.45,50,0.0,5.5);
  h_W_dt_3fold->GetXaxis()->SetTitle("W (GeV)");
  h_W_dt_3fold->GetYaxis()->SetTitle("|t-t_{min}| (GeV^{2})");
  h_W_dt_3fold->GetXaxis()->SetNdivisions(504);
  h_W_dt_3fold->GetYaxis()->SetNdivisions(505);
  h_W_dt_3fold->GetXaxis()->SetTitleSize(0.05);
  h_W_dt_3fold->GetYaxis()->SetTitleSize(0.05);
  h_W_dt_3fold->GetXaxis()->SetLabelSize(0.05);
  h_W_dt_3fold->GetYaxis()->SetLabelSize(0.05);
  T->Draw("abs(t-tmin):W >> h_W_dt_3fold","accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","colz");

  TCanvas *c6_ll = new TCanvas();
  h_W_dt_3fold->Draw("colz");
  c6_ll->Print("new_solid_W_dt_3fold.png");

  cscratch->cd();


  /* Figure 6 lower-right: 4-fold coincidence */
  TH2F* h_W_dt_4fold = new TH2F("h_W_dt_4fold","",50,4.0,4.45,50,0.0,5.5);
  h_W_dt_4fold->GetXaxis()->SetTitle("W (GeV)");
  h_W_dt_4fold->GetYaxis()->SetTitle("|t-t_{min}| (GeV^{2})");
  h_W_dt_4fold->GetXaxis()->SetNdivisions(504);
  h_W_dt_4fold->GetYaxis()->SetNdivisions(505);
  h_W_dt_4fold->GetXaxis()->SetTitleSize(0.05);
  h_W_dt_4fold->GetYaxis()->SetTitleSize(0.05);
  h_W_dt_4fold->GetXaxis()->SetLabelSize(0.05);
  h_W_dt_4fold->GetYaxis()->SetLabelSize(0.05);
  T->Draw("abs(t-tmin):W >> h_W_dt_4fold","accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2) && accep_p_1","colz");

  TCanvas *c6_lr = new TCanvas();
  h_W_dt_4fold->Draw("colz");
  c6_lr->Print("new_solid_W_dt_4fold.png");

  cscratch->cd();

  /* Figure 6 upper right: 3-fold coincidence */
  //  TH2F* h_W_tmin_3fold = new TH2F("h_W_tmin_3fold","",50,4.0,4.55,50,0.0,10.0);
  TH1F* h_W_t_3fold = new TH1F("h_W_t_3fold","",50,4.0,4.55);
  h_W_t_3fold->GetYaxis()->SetRangeUser(0.0,10.0);
  h_W_t_3fold->GetXaxis()->SetTitle("W (GeV)");
  h_W_t_3fold->GetYaxis()->SetTitle("-t (GeV^{2})");
  h_W_t_3fold->GetXaxis()->SetNdivisions(504);
  h_W_t_3fold->GetYaxis()->SetNdivisions(504);
  h_W_t_3fold->GetXaxis()->SetTitleSize(0.05);
  h_W_t_3fold->GetYaxis()->SetTitleSize(0.05);
  h_W_t_3fold->GetXaxis()->SetLabelSize(0.05);
  h_W_t_3fold->GetYaxis()->SetLabelSize(0.05);

  T->Draw("tmin:W","accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","");
  TGraph* g_W_tmin = new TGraph(T->GetEntries("accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)"),T->GetV2(),T->GetV1());

  T->Draw("tmax:W","accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","");
  TGraph* g_W_tmax = new TGraph(T->GetEntries("accep_e_1 && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)"),T->GetV2(),T->GetV1());
  g_W_tmax->SetMarkerColor(kRed);

  TCanvas *c6_ur_tmin = new TCanvas();
  //  h_W_tmin_3fold->Draw("colz");
  h_W_t_3fold->Draw();
  g_W_tmin->Draw("Psame");
  g_W_tmax->Draw("Psame");
  c6_ur_tmin->Print("new_solid_w_tmax_tmin_3fold.png");

  return 0;
}
