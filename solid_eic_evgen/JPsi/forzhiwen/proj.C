void proj(char *file){
  TChain *T = new TChain("T","T");
  T->AddFile(file);
  Int_t neve;
  T->SetBranchAddress("neve",&neve);
  T->GetEntry(T->GetEntries()-1);
  Double_t lumi = 1e37 * 1e-9 * 1e-24/neve;
  Double_t br = 5.94/100.;
  Double_t accep = pow(2.*3.1415926*(cos(8./180.*3.1415926)-0.85),2)*2.5;
  Double_t overall = lumi * br * accep * 50*3600*24*0.85;

    cout << "neve" << neve << endl;
    
  // TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TH1F *h1 = new TH1F("h1","h1",500,0.,5.);
  TH1F *h2 = new TH1F("h2","h2",500,8.5,10.5);  
//   T->Project("h2","11-p_e",Form("weight*dxs_23g*weight_decay*%f*accep_p*accep_je1*accep_je2*accep_e*(accep_p>0.&&W<4.12)",overall));
    T->Project("h2","11-p_e",Form("weight*dxs_2g*weight_decay*%f*accep_p*accep_je1*accep_je2*accep_e",overall));
//     T->Project("h2","11-p_e",Form("weight*dxs_2g*weight_decay*%f*accep_je1*accep_je2*accep_e",overall));

  cout << h2->GetSum() << endl;

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
  
}
