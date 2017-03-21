void plot_crosssection(){
TCanvas *c_crossection = new TCanvas("crossection","crossection",1800,700);
TF1 *f1 = new TF1("f1","fun_2g(x)",8,22);
// f->SetMaximum(1e2);
// f->SetMimimum(1e-4);
f1->Draw();
TF1 *f2 = new TF1("f2","fun_3g(x)",8,22);
f2->Draw("same");

TCanvas *c_decay = new TCanvas("decay","decay",1800,700);
// TF1 *f3=new TF1("f","3./4./3.1415926*(0.5 + 0.5*cos(x/180*3.1415926)*cos(x/180*3.1415926))",0,180);
TF1 *f3=new TF1("f","3./2.*(1 - cos(x)*cos(x))",0,3.1415926);
cout << "integral " << f3->Integral(0,3.1415926) << endl;
f3->Draw();
}

// Double_t fun_2g(Double_t x, Double_t){
Double_t fun_2g(Double_t x){
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
