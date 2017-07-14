void plot_crosssection_upsilon(){

  gStyle->SetOptStat(0);

  /* reference: solid cross section projection */
  TTree *t_ref1 = new TTree();
  t_ref1->ReadFile("references/solid_pac39_xsection_data/fig16_this_proposal_mean.csv","Egamma/F:sigma/F");

  // proton mass = 0.938272; // GeV
  // W = sqrt( 2*mp*Keq + mp*mp )
  //t_ref1->Draw("sigma/1.2:sqrt(2*0.938272*Egamma + 0.938272*0.938272)","",""); // original points scaled up by factor 1.2, see SOLID PAC proposal Fig 16
  t_ref1->Draw("sigma/1.2:Egamma","",""); // original points scaled up by factor 1.2, see SOLID PAC proposal Fig 16
  TGraphErrors *g_ref1 = new TGraphErrors(t_ref1->GetEntries(), t_ref1->GetV2(), t_ref1->GetV1());
  g_ref1->SetMarkerColor(kRed);
  g_ref1->SetMarkerStyle(20);


  TH1F* hframe_jpsi = new TH1F("hframe_jpsi","",40,0,40);
  hframe_jpsi->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  hframe_jpsi->GetYaxis()->SetTitle("#sigma (nb)");
  hframe_jpsi->GetXaxis()->SetRangeUser(3e-0,3e1);
  hframe_jpsi->GetYaxis()->SetRangeUser(1e-5,2e1);

  TH1F* hframe_upsilon = new TH1F("hframe_upsilon","",100,0,1000);
  hframe_upsilon->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  hframe_upsilon->GetYaxis()->SetTitle("#sigma (nb)");
  hframe_upsilon->GetXaxis()->SetRangeUser(3e-0,3e1);
  hframe_upsilon->GetYaxis()->SetRangeUser(1e-5,5e1);

  TH1F* hframe_upsilon_W = new TH1F("hframe_upsilon_W","",100,1,150);
  hframe_upsilon_W->GetXaxis()->SetTitle("W (GeV)");
  hframe_upsilon_W->GetYaxis()->SetTitle("#sigma (nb)");
  hframe_upsilon_W->GetXaxis()->SetRangeUser(3,3e1);
  hframe_upsilon_W->GetYaxis()->SetRangeUser(1e-5,5e1);

  /* reference lines for coverage */
  TLine *W_solid = new TLine(4.05,20,4.45,20);
  TLine *W_eic_5x50 = new TLine(8,25,30,25);
  TLine *W_eic_20x250 = new TLine(40,30,140,30);

  TF1 *f1_jpsi = new TF1("f1_jpsi","fun_2g_jpsi(x)",8,42);
  TF1 *f1_jpsi_ext = new TF1("f1_jpsi_ext","fun_2g_jpsi(x)",1,1000);
  TF1 *f1_upsilon = new TF1("f1_upsilon","fun_2g_upsilon(x)",1,1000);
  // f->SetMaximum(1e2);
  // f->SetMimimum(1e-4);

  f1_jpsi->SetLineColor(kRed);
  f1_jpsi_ext->SetLineColor(kRed);
  f1_jpsi_ext->SetLineStyle(2);
  f1_jpsi_ext->SetLineWidth(2);
  f1_upsilon->SetLineColor(kBlue);

  //  TF1 *f1_jpsi_W = new TF1("f1_jpsi_W","fun_2g_jpsi_W(x)",4,22);
  TF1 *f1_jpsi_W = new TF1("f1_jpsi_W","fun_2g_jpsi_W(x)",4,7);
  //  TF1 *f1_jpsi_W_ext = new TF1("f1_jpsi_W_ext","fun_2g_jpsi_W(x)",4,1000);
  TF1 *f1_upsilon_W = new TF1("f1_upsilon_W","fun_2g_upsilon_W(x)",10.4,18);

  f1_jpsi_W->SetLineColor(kOrange);
  //  f1_jpsi_W_ext->SetLineColor(kOrange+1);
  //  f1_jpsi_W_ext->SetLineWidth(3);
  f1_upsilon_W->SetLineColor(kBlue);
  f1_upsilon_W->SetLineWidth(3);

  /* J/Psi only */
  TCanvas *c_crossection_jpsi = new TCanvas("crossection_jpsi","crossection_jpsi");//,1800,700);
  hframe_jpsi->Draw();
  f1_jpsi->Draw("same");

  g_ref1->Draw("Psame");

  c_crossection_jpsi->SetLogx(1);
  c_crossection_jpsi->SetLogy(1);

  return;
  /* J/Psi and Upsilon */
  TCanvas *c_crossection_upsilon = new TCanvas("crossection","crossection");//,1800,700);
  hframe_upsilon->Draw();
  f1_jpsi->Draw("same");
  f1_upsilon->Draw("same");

  c_crossection_upsilon->SetLogx(1);
  c_crossection_upsilon->SetLogy(1);

  /* J/Psi vs W */
  TCanvas *c_crossection_W = new TCanvas("crossection_jpsi_W","crossection_jpsi_W");//,1800,700);
  hframe_upsilon_W->Draw();
  f1_jpsi_W->Draw("same");

  W_solid->Draw();
  W_eic_5x50->Draw();
  //  W_eic_20x250->Draw();

  c_crossection_W->SetLogx(1);
  c_crossection_W->SetLogy(1);
  c_crossection_W->Print("calc_xsec_vs_W_jpsi.png");

  /* J/Psi and Upsilon vs W */
  TCanvas *c_crossection_upsilon_W = new TCanvas("crossection_W","crossection_W");//,1800,700);
  hframe_upsilon_W->Draw();
  f1_jpsi_W->Draw("same");
  f1_upsilon_W->Draw("same");

  W_solid->Draw();
  W_eic_5x50->Draw();
  //  W_eic_20x250->Draw();

  c_crossection_upsilon_W->SetLogx(1);
  c_crossection_upsilon_W->SetLogy(1);
  c_crossection_upsilon_W->Print("calc_xsec_vs_W_upsilon.png");

  //TF1 *f2 = new TF1("f2","fun_3g(x)",8,22);
  //f2->Draw("same");

  //TCanvas *c_decay = new TCanvas("decay","decay",1800,700);
  //// TF1 *f3=new TF1("f","3./4./3.1415926*(0.5 + 0.5*cos(x/180*3.1415926)*cos(x/180*3.1415926))",0,180);
  //TF1 *f3=new TF1("f","3./2.*(1 - cos(x)*cos(x))",0,3.1415926);
  //cout << "integral " << f3->Integral(0,3.1415926) << endl;
  //f3->Draw();
}

// Double_t fun_2g(Double_t x, Double_t){
Double_t fun_2g_jpsi(Double_t x){
  Double_t t=0;
  x=sqrt((x+0.938)*(x+0.938)-x*x);   //convert from E_gamma to W

  Double_t N2g = 7.5671e3;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t M = 3.097;
  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff;

  result = result / 1.13;    //total crossection is differential crossection at t=0 and divide by slope

  return result;
}

/* x is W here */
Double_t fun_2g_jpsi_W(Double_t x){
  Double_t t=0;

  Double_t N2g = 7.5671e3;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t M = 3.097;
  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff;

  result = result / 1.13;    //total crossection is differential crossection at t=0 and divide by slope

  return result;
}

// same as fun_2g except change M to upsilon mass
Double_t fun_2g_upsilon(Double_t x){
  Double_t t=0;
  x=sqrt((x+0.938)*(x+0.938)-x*x);   //convert from E_gamma to W

  Double_t ratio_upsilon_jpsi = 1;//0.0015 * 10;
  Double_t N2g = 7.5671e3 * ratio_upsilon_jpsi;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t M = 9.460;
  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff;

  result = result / 1.13;    //total crossection is differential crossection at t=0 and divide by slope

  return result;
}


/* x is W here */
Double_t fun_2g_upsilon_W(Double_t x){
  Double_t t=0;

  Double_t ratio_upsilon_jpsi = 1.; //0.0015 * 10;
  Double_t N2g = 7.5671e3 * ratio_upsilon_jpsi;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t M = 9.460;
  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff;

  result = result / 1.13;    //total crossection is differential crossection at t=0 and divide by slope

  return result;
}

// Double_t fun_3g(Double_t x, Double_t){
Double_t fun_3g(Double_t x){
  Double_t t=0;
  x=sqrt((x+0.938)*(x+0.938)-x*x);  //convert from E_gamma to W
  Double_t N3g = 2.894e3;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t M = 3.097;

  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N3g*v/R/R/R/R/M/M/M/M*pow(1-xp,0)*ff;

  result = result / 1.13; //total crossection is differential crossection at t=0 and divide by slope

  return result;
}
