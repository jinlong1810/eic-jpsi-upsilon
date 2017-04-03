int plot_kinematics_eic()
{

  TFile *fin = new TFile("../data/output_e_500GeV_1M.root","OPEN");
  TTree *T = (TTree*)fin->Get("T");

  //
  TCanvas *c1 = new TCanvas();
  T->Draw("Q2:W","","colz");

  //
  TCanvas *c2 = new TCanvas();
  T->Draw("tmin:W","","colz");

  //
  TCanvas *c3 = new TCanvas();
  T->Draw("t-tmin:W","","colz");

  return 0;
}
