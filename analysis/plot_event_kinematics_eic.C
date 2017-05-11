int plot_event_kinematics_eic()
{

  TFile *fin = new TFile("../data/sim_te_b5GeV_p50GeV_mjpsi.root","OPEN");
  TTree *T = (TTree*)fin->Get("T");

  //  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TCanvas *cscratch = new TCanvas();

  //
  TCanvas *c1 = new TCanvas();
  T->Draw("Q2:W","","colz");
  c1->SetLogz(1);
  c1->Print("new_eic_Q2_W.png");

  //
  TCanvas *c2 = new TCanvas();
  T->Draw("tmin:W","","colz");
  c2->SetLogz(1);
  c2->Print("new_eic_tmin_W.png");

  //
  TCanvas *c3 = new TCanvas();
  T->Draw("t-tmin:W","","colz");
  c3->SetLogz(1);
  c3->Print("new_eic_dt_W.png");

  return 0;
}
