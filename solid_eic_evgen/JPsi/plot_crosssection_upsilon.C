void plot_crosssection_upsilon(){

  gStyle->SetOptStat(0);

  TH1F* hframe_jpsi = new TH1F("hframe_jpsi","",40,0,40);
  hframe_jpsi->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  hframe_jpsi->GetYaxis()->SetTitle("#sigma (nb)");
  hframe_jpsi->GetXaxis()->SetRangeUser(5e-0,5e1);
  hframe_jpsi->GetYaxis()->SetRangeUser(5e-4,2e1);

  TH1F* hframe_upsilon = new TH1F("hframe_upsilon","",1000,0,1000);
  hframe_upsilon->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  hframe_upsilon->GetYaxis()->SetTitle("#sigma (nb)");
  hframe_upsilon->GetXaxis()->SetRangeUser(5e-0,5e3);
  hframe_upsilon->GetYaxis()->SetRangeUser(5e-4,5e1);

  TF1 *f1_jpsi = new TF1("f1_jpsi","fun_2g_jpsi(x)",8,42);
  TF1 *f1_jpsi_ext = new TF1("f1_jpsi_ext","fun_2g_jpsi(x)",1,1000);
  TF1 *f1_upsilon = new TF1("f1_upsilon","fun_2g_upsilon(x)",1,1000);
  // f->SetMaximum(1e2);
  // f->SetMimimum(1e-4);

  f1_jpsi->SetLineColor(kRed);
  f1_jpsi_ext->SetLineColor(kRed);
  f1_jpsi_ext->SetLineStyle(2);
  f1_upsilon->SetLineColor(kBlue);

  /* J/Psi only */
  TCanvas *c_crossection_jpsi = new TCanvas("crossection_jpsi","crossection_jpsi",1800,700);
  hframe_jpsi->Draw();
  f1_jpsi->Draw("same");

  c_crossection_jpsi->SetLogx(1);
  c_crossection_jpsi->SetLogy(1);


  /* J/Psi and Upsilon */
  TCanvas *c_crossection_upsilon = new TCanvas("crossection","crossection",1800,700);
  hframe_upsilon->Draw();
  f1_jpsi->Draw("same");
  f1_jpsi_ext->Draw("same");
  f1_upsilon->Draw("same");

  c_crossection_upsilon->SetLogx(1);
  c_crossection_upsilon->SetLogy(1);

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

// same as fun_2g except change M to upsilon mass
Double_t fun_2g_upsilon(Double_t x){
  Double_t t=0;  
  x=sqrt((x+0.938)*(x+0.938)-x*x);   //convert from E_gamma to W
  
  Double_t N2g = 7.5671e3;
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
