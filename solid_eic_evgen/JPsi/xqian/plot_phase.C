void plot_phase(){
  gStyle->SetOptStat(1);
  gStyle->SetLabelSize(0.08,"Y");
  gStyle->SetTitleSize(0.08,"Y");
  gStyle->SetLabelSize(0.08,"X");
  gStyle->SetTitleSize(0.08,"X");
  gStyle->SetNdivisions(506,"X");
  gStyle->SetNdivisions(506,"Y");
  

  TChain *T = new TChain("T","T");
  T->AddFile("temp_p.root");
  gStyle->SetPalette(1);
  TH2F *h1 = new TH2F("h1","h1",100,7.,17,100,0.,3.);
  TH2F *h2 = new TH2F("h2","h2",100,11.,28,100,2.5,8.);
  TH2F *h3 = new TH2F("h3","h3",100,11.,28,100,2.5,8.);
  TH2F *h4 = new TH2F("h4","h4",100,7.,17,100,0.,4.5);
  T->Project("h1","p_e:theta_e","weight*weight1*dxs*accep_p*accep_p*accep_je1*accep_je2");
//   T->Project("h1","p_e:theta_e","dxs*accep_p*accep_p*accep_je1*accep_je2");  
  T->Project("h2","p_je1:theta_je1","weight*weight1*dxs*accep_p*accep_p*accep_je1*accep_je2");
  T->Project("h3","p_je2:theta_je2","weight*weight1*dxs*accep_p*accep_p*accep_je1*accep_je2");
  T->Project("h4","p_p:theta_p","weight*weight1*dxs*accep_p*accep_p*accep_je1*accep_je2");

  cout << h1->GetSum() << endl;
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  c1->cd(1);
  h1->Draw("COLZ");
  h1->SetTitle("Scattered Electron");
  h1->SetXTitle("Angle (deg)");
  h1->SetYTitle("Momentum (GeV/c)");
  
  c1->cd(2);
  h2->Draw("COLZ");
  h2->SetTitle("Decay Electron");
  h2->SetXTitle("Angle (deg)");
  h2->SetYTitle("Momentum (GeV/c)");

  c1->cd(3);
  h3->Draw("COLZ");
  h3->SetTitle("Decay Positron");
  h3->SetXTitle("Angle (deg)");
  h3->SetYTitle("Momentum (GeV/c)");
  c1->cd(4);
  h4->Draw("COLZ");
  h4->SetTitle("Recoil Proton");
  h4->SetXTitle("Angle (deg)");
  h4->SetYTitle("Momentum (GeV/c)");
}
