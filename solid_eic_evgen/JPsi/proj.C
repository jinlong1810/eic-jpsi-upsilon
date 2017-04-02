void proj(char *file,string type){
  gROOT->Reset();
  gStyle->SetPalette(1);
//   gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);  
  
  //Ebeam
  Double_t cov= 1e-9 * 1e-24; //nb to cm2 coversion
  Double_t br = 5.94/100.;
  Double_t time = 50*3600*24;  //50 days in seconds
  Double_t eff = 0.85;

  Double_t lumi = 1.2e37;  // 1.2e37/cm2/s is from 3uA on 15cm long LH2 target  
  
  //Gbeam  
//   double lumi_g = 8e35; // 10e-9/1.6e-19*30*0.071*6.02e23  100nA e and 30cm LH2
  Double_t rad = 1.7e-2;    // 15cm LH2, 15/890.4=1.7e-2 radiator

  Double_t norm_experiment;  
  char accep_normal[200];
  if (type=="e")  {
    norm_experiment = cov * br * eff * time * lumi;  
//   sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1+accep_e_2)*(accep_p_1+accep_p_2)");
//   sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)*(accep_p_1)");    
  sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)*(accep_p_1+accep_p_2)");    
    
//   sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_p_1+accep_p_2)");
//     sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)"); 
  }
  else if (type=="g") {
    norm_experiment = cov * br * eff * time * lumi * rad;
    sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_p_1+accep_p_2)*Gflux"); 
//     sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*Gflux");     
  }
  else {cout << "wrong type" << endl;  return 0;}  

  TChain *T = new TChain("T","T");
  T->AddFile(file);
  cout << " open file " << file << endl;   
  Double_t phasespace;
  T->SetBranchAddress("phasespace",&phasespace);  
  Int_t neve;
  T->SetBranchAddress("neve",&neve);
  T->GetEntry(T->GetEntries()-1);
  cout << " throw events " << neve << endl;
  
  Double_t norm_simulation=phasespace/neve;
  Double_t overall=norm_experiment*norm_simulation;    
  
  const int n=13;
  TH2F *hThetaP[n][4];  //#theta(deg);P(GeV)
  TH2F *htemp=new TH2F("htemp","htemp",180,0,180,120,0,12);
  char *content[4]={"p_e:theta_e","p_p:theta_p","p_je1:theta_je1","p_je2:theta_je2"};
//     TH1F *htemp=new TH1F("htemp","htemp",80,0,8);
//     char *content[4]={"p_e","p_p","p_je1","p_je2"};
    
//   char *weight[10]={
//     "",
//     "weight*dxs","weight*dxs_2g","weight*dxs_23g",
//     "weight*dxs*weight_decay","weight*dxs_2g*weight_decay","weight*dxs_23g*weight_decay",
//     "weight*dxs*weight_decay*accep_p*accep_je1*accep_je2*accep_e*(accep_p>0.&&W<4.12)","weight*dxs_2g*weight_decay*accep_p*accep_je1*accep_je2*accep_e*(accep_p>0.&&W<4.12)","weight*dxs_23g*weight_decay*accep_p*accep_je1*accep_je2*accep_e*(accep_p>0.&&W<4.12)"
//   };
  
//     char *weight[n]={
//     "",
//     "weight*dxs","weight*dxs*weight_decay","weight*dxs*weight_decay","dxs*accep_je1*accep_je2*accep_e*accep_p",
// //     "weight*dxs","weight*dxs*weight_decay","weight*dxs*weight_decay","weight*dxs*weight_decay*accep_je1*accep_je2*accep_e*accep_p",
//     "weight*dxs_2g","weight*dxs_2g*weight_decay","weight*dxs_2g*weight_decay","weight*dxs_2g*weight_decay*accep_je1*accep_je2*accep_e*accep_p",    "weight*dxs_23g","weight*dxs_23g*weight_decay","weight*dxs_23g*weight_decay","weight*dxs_23g*weight_decay*accep_je1*accep_je2*accep_e*accep_p"
//   };

///3 fold no p
// char *weight[n]={
// "",
// "weight*dxs","weight*dxs*weight_decay","weight*dxs*weight_decay","weight*dxs*weight_decay*accep_je1*accep_je2*accep_e*%f",
// "weight*dxs_2g","weight*dxs_2g*weight_decay","weight*dxs_2g*weight_decay","weight*dxs_2g*weight_decay*accep_je1*accep_je2*accep_e*%f",    "weight*dxs_23g","weight*dxs_23g*weight_decay","weight*dxs_23g*weight_decay","weight*dxs_23g*weight_decay*accep_je1*accep_je2*accep_e*%f"
// };

///psudo 4 fold with W cut 
//     char *weight[n]={
//     "",
//     "weight*dxs","weight*dxs*weight_decay","weight*dxs*weight_decay","weight*dxs*weight_decay*(accep_p>0.&&W<4.12)*accep_je1*accep_je2*accep_e*%f",
//     "weight*dxs_2g","weight*dxs_2g*weight_decay","weight*dxs_2g*weight_decay","weight*dxs_2g*weight_decay*(accep_p>0.&&W<4.12)*accep_je1*accep_je2*accep_e*%f",    "weight*dxs_23g","weight*dxs_23g*weight_decay","weight*dxs_23g*weight_decay","weight*dxs_23g*weight_decay*(accep_p>0.&&W<4.12)*accep_je1*accep_je2*accep_e*%f"
//   };  

//general weight
char *weight[n]={
"1",
"dxs","dxs*weight","dxs*weight*weight_decay","dxs*weight*weight_decay*%s*%f",
"dxs_2g","dxs_2g*weight","dxs_2g*weight*weight_decay","dxs_2g*weight*weight_decay*%s*%f",
"dxs_23g","dxs_23g*weight","dxs_23g*weight*weight_decay","dxs_23g**weight*weight_decay*%s*%f"
};

  
  
  TCanvas *c_ThetaP = new TCanvas("ThetaP","ThetaP",800,1200);
  c_ThetaP->Divide(4,n);   
  for (Int_t i=0;i<n;i++){
      for (Int_t j=0;j<4;j++){
	T->Project("htemp",content[j],Form(weight[i],accep_normal,overall));
// 	T->Project("htemp",content[j],Form("%s%s",Form(weight[i],accep_normal,overall),"*(W<4.25)"));	
	c_ThetaP->cd(i*4+j+1);
	gPad->SetLogz(1);
	hThetaP[i][j]=(TH2F*) htemp->Clone();
	hThetaP[i][j]->Draw("colz");
      }
//       if (i==4 || i==8 || i==12) cout << hThetaP[i][0]->Integral() << endl;
//       if (i==4 || i==8 || i==12) cout << htemp->Integral() << endl;
  }      
    
    TCanvas *c_ThetaP_2g = new TCanvas("ThetaP_2g","ThetaP_2g",800,1200);
    c_ThetaP_2g->Divide(4,4);   
    for (Int_t i=5;i<9;i++){
	for (Int_t j=0;j<4;j++){
	  c_ThetaP_2g->cd((i-5)*4+j+1);
	  gPad->SetLogz(1);
	  hThetaP[i][j]->Draw("colz");
// 	  cout << hThetaP[i][j]->Integral() << endl;	  
	}
    }
//   c_ThetaP_2g->SaveAs("ThetaP_2g.png");

  char *title[4]={"scattered e^{-}","recoil p","decay e^{-}","decay e^{+}"};
    TCanvas *c_ThetaP_2g_final = new TCanvas("ThetaP_2g_final","ThetaP_2g_final",1200,1000);
    c_ThetaP_2g_final->Divide(2,2);   
    for (Int_t i=8;i<9;i++){
	for (Int_t j=0;j<4;j++){
	  c_ThetaP_2g_final->cd(j+1);
	  gPad->SetLogz(1);
	  hThetaP[i][j]->SetAxisRange(0,40,"X");
	  hThetaP[i][j]->SetTitle(Form("Counts (%s);#theta (deg);P (GeV)",title[j]));	  
// 	  gStyle->SetTitleAlign(23);	  
// 	  hThetaP[i][j]->GetXaxis()->SetLabelSize(0.1);
// 	  hThetaP[i][j]->GetYaxis()->SetLabelSize(0.1);	  
// 	  hThetaP[i][j]->GetXaxis()->SetTitleSize(0.1);
// 	  hThetaP[i][j]->GetYaxis()->SetTitleSize(0.1);	  	  
	  hThetaP[i][j]->Draw("colz");
	}
    }   
//   c_ThetaP_2g_final->SaveAs("ThetaP_2g_final.png");
  
cout << "total events with dxs_2g crossection " << hThetaP[8][3]->Integral() << endl;
  
TCanvas *c_decay = new TCanvas("decay","decay",1400,600);
c_decay->Divide(2,1);   
c_decay->cd(1);
gPad->SetLogz();	
T->Project("hdecay_gen(100,-1,1,360,0,360)","phi_cm:cos(theta_cm*3.1416/180)",Form("dxs_2g*weight*weight_decay*%f",overall));
// T->Project("hdecay_gen(100,-1,1)","cos(theta_cm*3.1416/180)",Form("dxs_2g*weight*weight_decay*%f",overall));
hdecay_gen->SetTitle("J/#psi decay (generated);cos(#theta_{GJ});#phi_{GJ} (deg)");
hdecay_gen->GetXaxis()->SetLabelSize(0.04);
hdecay_gen->GetYaxis()->SetLabelSize(0.04);
hdecay_gen->Draw("colz");
c_decay->cd(2);
gPad->SetLogz();	
T->Project("hdecay_acc(100,-1,1,360,0,360)","phi_cm:cos(theta_cm*3.1416/180)",Form("dxs_2g*weight*weight_decay*%s*%f",accep_normal,overall));
hdecay_acc->SetTitle("J/#psi decay (accepted);cos(#theta_{GJ});#phi_{GJ} (deg)");
hdecay_acc->GetXaxis()->SetLabelSize(0.04);
hdecay_acc->GetYaxis()->SetLabelSize(0.04);
hdecay_acc->Draw("colz");


TCanvas *c = new TCanvas("c","c",1800,700);
c->Divide(2,2);   
c->cd(1);
gPad->SetLogz();	
T->Project("hWt(90,0,3,100,6,16)","Keq:p_e",Form("dxs_2g*weight*weight_decay*%s*%f",accep_normal,overall));
// T->Project("hWt(60,0,6,100,4,5)","W:(t-tmin)",Form("dxs_2g*weight*weight_decay*%s*%f",accep_normal,overall));
// T->Project("hWt(60,0,6,100,4,5)","W:(t-tmin)",Form("(Q2<0.01)*(dxs_2g*weight*weight_decay*%s*%f)",accep_normal,overall));
hWt->Draw("colz");
// c->cd(2);
// gPad->SetLogz();	
// T->Project("hWKeq(50,7,12,60,3.5,5)","W:Keq","weight*dxs_2g*weight_decay*accep_p*accep_je1*accep_je2*accep_e");
// hWKeq->Draw("colz");
c->cd(2);
gPad->SetLogz();	
T->Project("htheta_je1_je2(80,0,40,80,0,40)","theta_je1:theta_je2",Form("dxs_2g*weight*weight_decay*%s*%f",accep_normal,overall));
htheta_je1_je2->Draw("colz");
c->cd(3);
gPad->SetLogz();	
T->Project("htheta_e_je1(80,0,40,80,0,40)","theta_e:theta_je1",Form("dxs_2g*weight*weight_decay*%s*%f",accep_normal,overall));
htheta_e_je1->Draw("colz");
c->cd(4);
gPad->SetLogz();
double nevent=T->GetEntries();
// T->Project("hweight_decay(100,-1,1)","cos(theta_cm*3.1416/180)",Form("weight_decay/%f",nevent));
T->Project("hweight_decay(100,0,180)","theta_cm",Form("weight_decay/%f",nevent));
hweight_decay->Draw("colz");
// cout << "decay integral " << hweight_decay->Integral() << endl;

TCanvas *c_G = new TCanvas("G","G",1800,700);
c_G->Divide(2,1);   
c_G->cd(1);
gPad->SetLogy();
// T->Project("hXGbeam(50,7,12)","Gbeam",Form("dxs_2g*weight*%f",norm_simulation));
T->Project("hXGbeam(50,7,12)","Gbeam");
hXGbeam->Draw("colz");
c_G->cd(2);
gPad->SetLogy();
T->Project("hXGbeam_Norm(50,7,12)","Gbeam",Form("dxs_2g*weight*Gflux*%f",norm_simulation));
hXGbeam_Norm->Draw("colz");


/*
  TCanvas *c = new TCanvas("c","c",800,600);
  c->Divide(3,1);
  
  TH1F *h2 = new TH1F("h2","h2",500,8.5,10.5);  
  T->Project("h2","11-p_e",Form("weight*dxs_23g*weight_decay*%f*accep_p*accep_je1*accep_je2*accep_e*(accep_p>0.&&W<4.12)",overall));
  c->cd(1);
//   h2->SetMarkerColor(kRed);
  h2->Draw();  
//   cout << h2->GetSum() << endl;
  cout << h2->Integral() << endl;  
  
  TH1F *h3 = new TH1F("h3","h3",500,8.5,10.5);  
  T->Project("h3","11-p_e",Form("weight*dxs_2g*weight_decay*%f*accep_p*accep_je1*accep_je2*accep_e*(accep_p>0.&&W<4.12)",overall));
  c->cd(2); 
  h3->Draw();    
//   cout << h2->GetSum() << endl;
  cout << h3->Integral() << endl;
  
  TH1F *h4 = new TH1F("h4","h4",500,8.5,10.5);  
  T->Project("h4","11-p_e",Form("weight*dxs*weight_decay*%f*accep_p*accep_je1*accep_je2*accep_e*(accep_p>0.&&W<4.12)",overall));
  c->cd(3); 
  h4->Draw();     
//   cout << h2->GetSum() << endl;
  cout << h4->Integral() << endl;  

  Double_t energy_limit[5];
  Double_t energy_center[4];
  Double_t count = 0,sum=0;
  Double_t total = h2->GetSum()/4.;
  Int_t abc = 0;
  energy_limit[0] = 8.5;
  for (Int_t i =0;i!=500;i++){
    count += h2->GetBinContent(i+1);
    sum += h2->GetBinCenter(i+1)*h2->GetBinContent(i+1);
    if (count > total){
      energy_center[abc] = sum/count;
      energy_limit[abc+1] = h2->GetBinCenter(i+1);
      count = 0;
      sum = 0;
      abc ++;
    }
  }
  energy_center[abc] = sum/count;
  energy_limit[abc+1] = h2->GetBinCenter(i+1);
  
  TH1F *h1 = new TH1F("h1","h1",500,0.,5.);  
  Double_t tto_ave[4][20];
  Double_t nevents[4][20];
  for (Int_t i=0;i!=4;i++){
    Double_t Emin = energy_limit[i];
    Double_t Emax = energy_limit[i+1];
    h1->Reset();
    TString cut;
    cut.Form("weight*dxs_23g*weight_decay*accep_p*accep_je1*accep_je2*accep_e*%f*(accep_p>0.&&11-p_e>%f&&11-p_e<%f)",overall,Emin,Emax);
    //out << cut << endl;
    TCut ecut = cut;
    T->Project("h1","t-tmin",ecut);
    Double_t ave = h1->GetSum()/12.;
    count = 0;
    sum = 0;
    abc = 0;
    for (Int_t j=0;j!=500;j++){
      count += h1->GetBinContent(j+1);
      sum += h1->GetBinCenter(j+1)*h1->GetBinContent(j+1);
      if ( count > ave){
	tto_ave[i][abc] = sum/count;
	nevents[i][abc] = count;
	sum = 0;
	count = 0;
	//cout << tto_ave[i][abc] << "\t" << nevents[i][abc] << endl;
	abc ++;
      }
    }
    tto_ave[i][abc] = sum/count;
    nevents[i][abc] = count;
    
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  TF1 *f1 = new TF1("f1","0.94*exp(-0.97*(x-1.5))",0.,4.);
  for (Int_t i=0;i!=4;i++){
    c1->cd(i+1);
    gPad->SetLogy(1);
    Double_t x[20],y[20],y_err[20];
    Int_t ncount = 0;
    for (Int_t j=0;j!=12;j++){
      
      if (nevents[i][j] >10.0){
	x[j] = tto_ave[i][j];
	y[j] = 0.94 * exp( -0.97*(x[j]-1.5));
	cout << nevents[i][j] << endl;
	y_err[j] = y[j] / sqrt(nevents[i][j])*sqrt(1.15);
	ncount ++;

      }
    }
    TGraphErrors *g1 = new TGraphErrors(ncount,x,y,0,y_err);
    g1->Draw("A*");
    g1->GetXaxis()->SetTitle("|t-tmin| GeV^{2}");
    g1->GetYaxis()->SetTitle("#frac{d#sigma}{dt} (nb/GeV^{2})");
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(2);
    g1->SetLineColor(2);
    TString title;
    title.Form("Virtual Photon Energy %3.2f GeV", energy_center[i]);
    g1->SetTitle(title);
    f1->Draw("same");
    f1->SetLineColor(1);
    cout << ncount << endl;
  }
  */
}
