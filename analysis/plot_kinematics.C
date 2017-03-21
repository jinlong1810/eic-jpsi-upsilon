int plot_kinematics()
{

  TFile *fin = new TFile("output_e_11GeV_1M_accep.root","OPEN");
  TTree *T = (TTree*)fin->Get("T");

  // Figure 6 upper-left: 3-fold coincidence
  TCanvas *c1 = new TCanvas();
  T->Draw("Q2:W","(accep_e_1 || accep_e_2) && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","colz");

  // Figure 6 upper right: 3-fold coincidence
  TCanvas *c2 = new TCanvas();
  T->Draw("tmin:W","(accep_e_1 || accep_e_2) && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","colz");

  // Figure 6 lower-left: 3-fold coincidence
  TCanvas *c3 = new TCanvas();
  T->Draw("t-tmin:W","(accep_e_1 || accep_e_2) && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2)","colz");

  // Figure 6 lower-right: 4-fold coincidence
  TCanvas *c4 = new TCanvas();
  T->Draw("t-tmin:W","(accep_e_1 || accep_e_2) && (accep_je1_1 || accep_je1_2) && (accep_je2_1 || accep_je2_2) && accep_p_1","colz");

  return 0;
}
