int plot_compare_codes()
{

  TFile *fin1 = new TFile("../solid_evgen/JPsi/sim_te_b11GeV_p0GeV_mjpsi.root","OPEN");
  TTree *T1 = (TTree*)fin1->Get("T");

  TFile *fin2 = new TFile("../solid_eic_evgen/JPsi/sim_te_b11GeV_p0GeV_mjpsi.root","OPEN");
  TTree *T2 = (TTree*)fin2->Get("T");

  //  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  // Theta angle of J/Psi
  TCanvas *c1 = new TCanvas();
  TH1F* h1_1 = new TH1F("h1_1","J/Psi theta angle;# entries;#theta^{J/#Psi}",100,0,10);
  TH1F* h1_2 = (TH1F*)h1_1->Clone("h1_2");

  h1_1->SetLineColor(kBlue);
  h1_2->SetLineColor(kOrange+2);
  h1_2->SetLineStyle(2);

  T1->Draw("theta_jpsi >> h1_1");
  T2->Draw("theta_jpsi >> h1_2");

  h1_1->Draw();
  h1_2->Draw("same");

  c1->Print("plots/compare_codes_1.eps");


  // Total momentum of J/Psi
  TCanvas *c2 = new TCanvas();
  TH1F* h2_1 = new TH1F("h2_1","J/Psi total momentum;# entries;p^{J/#Psi}",150,0,15);
  TH1F* h2_2 = (TH1F*)h2_1->Clone("h2_2");

  h2_1->SetLineColor(kBlue);
  h2_2->SetLineColor(kOrange+2);
  h2_2->SetLineStyle(2);

  T1->Draw("p_jpsi >> h2_1");
  T2->Draw("p_jpsi >> h2_2");

  h2_1->Draw();
  h2_2->Draw("same");

  c2->SetLogy(1);

  c2->Print("plots/compare_codes_2.eps");

  return 0;
}
