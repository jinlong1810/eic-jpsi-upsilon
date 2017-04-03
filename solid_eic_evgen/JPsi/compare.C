#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>

using namespace std;

void compare()
{
  gROOT->Reset();
  gStyle->SetOptStat(0);
  
  Double_t cov= 1e-9 * 1e-24; //nb to cm2 coversion
  Double_t lumi = 1e37;  // 1e37/cm2/s is from 3nA on 15cm long LH2 target
  Double_t br = 5.94/100.;
//   Double_t accep = pow(2.*3.1415926*(cos(8./180.*3.1415926)-0.85),2)*2.5;
  Double_t range_angle = pow(2.*3.1415926*(cos(0./180.*3.1415926)-cos(40./180.*3.1415926)),2); // both e and proton solid angle range
    Double_t range_P= 3.-0.;  // e mom range
//   Double_t range_P= 3.-0.5;  // e mom range  
  Double_t time = 50*3600*24;  //50 days in seconds
  Double_t eff = 0.85;
  Double_t overall_NOneve = cov *lumi * br * range_angle * range_P * eff * time;  

 const int n=13;
int pos[n]={250,270,280,290,295,300,305,310,320,330,340,350,360};  


///nocut
char *weight_nocut[6][3]={
 //4fold
"dxs*weight*weight_decay*accep_je1*accep_je2*accep_e*accep_p*%f",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_e*accep_p*%f",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_e*accep_p*%f",
//3foldnp
"dxs*weight*weight_decay*accep_je1*accep_je2*accep_e*%f",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_e*%f",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_e*%f",
//3foldne
"dxs*weight*weight_decay*accep_je1*accep_je2*accep_p*%f",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_p*%f",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_p*%f",
//3foldnje
"dxs*weight*weight_decay*accep_e*accep_je2*accep_p*%f",
"dxs_2g*weight*weight_decay*accep_e*accep_je2*accep_p*%f",
"dxs_23g*weight*weight_decay*accep_e*accep_je2*accep_p*%f",
//2foldj
"dxs*weight*weight_decay*accep_je1*accep_je2*%f",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*%f",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*%f",
//2folep
"dxs*weight*weight_decay*accep_e*accep_p*%f",
"dxs_2g*weight*weight_decay*accep_e*accep_p*%f",
"dxs_23g*weight*weight_decay*accep_e*accep_p*%f",
};

///W<4.25
char *weight_Wcut[6][3]={
 //4fold
"dxs*weight*weight_decay*accep_je1*accep_je2*accep_e*accep_p*%f*(W<4.25)",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_e*accep_p*%f*(W<4.25)",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_e*accep_p*%f*(W<4.25)",
//3foldnp
"dxs*weight*weight_decay*accep_je1*accep_je2*accep_e*%f*(W<4.25)",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_e*%f*(W<4.25)",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_e*%f*(W<4.25)",
//3foldne
"dxs*weight*weight_decay*accep_je1*accep_je2*accep_p*%f*(W<4.25)",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_p*%f*(W<4.25)",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_p*%f*(W<4.25)",
//3foldnje
"dxs*weight*weight_decay*accep_e*accep_je2*accep_p*%f*(W<4.25)",
"dxs_2g*weight*weight_decay*accep_e*accep_je2*accep_p*%f*(W<4.25)",
"dxs_23g*weight*weight_decay*accep_e*accep_je2*accep_p*%f*(W<4.25)",
//2foldj
"dxs*weight*weight_decay*accep_je1*accep_je2*%f*(W<4.25)",
"dxs_2g*weight*weight_decay*accep_je1*accep_je2*%f*(W<4.25)",
"dxs_23g*weight*weight_decay*accep_je1*accep_je2*%f*(W<4.25)",
//2folep
"dxs*weight*weight_decay*accep_e*accep_p*%f*(W<4.25)",
"dxs_2g*weight*weight_decay*accep_e*accep_p*%f*(W<4.25)",
"dxs_23g*weight*weight_decay*accep_e*accep_p*%f*(W<4.25)",
};

///electron at forwad and large
// char *weight[4][3]={
//  //4fold
// "dxs*weight*weight_decay*accep_je1*accep_je2*accep_e_both*accep_p*%f",
// "dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_e_both*accep_p*%f",
// "dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_e_both*accep_p*%f",
// //3foldnp
// "dxs*weight*weight_decay*accep_je1*accep_je2*accep_e_both*%f",
// "dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_e_both*%f",
// "dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_e_both*%f",
// //3foldne
// "dxs*weight*weight_decay*accep_je1*accep_je2*accep_p*%f",
// "dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_p*%f",
// "dxs_23g*weight*weight_decay*accep_je1*accep_je2*accep_p*%f",
// //2fold
// "dxs*weight*weight_decay*accep_je1*accep_je2*%f",
// "dxs_2g*weight*weight_decay*accep_je1*accep_je2*%f",
// "dxs_23g*weight*weight_decay*accep_je1*accep_je2*%f"
// };

char *hstname[6][3]={
"dxs_t_4fold","dxs_2g_4fold","dxs_23g_4fold",
"dxs_t_3foldnp","dxs_2g_3foldnp","dxs_23g_3foldnp",
"dxs_t_3foldne","dxs_2g_3foldne","dxs_23g_3foldne",
"dxs_t_3foldnje","dxs_2g_3foldnje","dxs_23g_3foldnje",
"dxs_t_2foldj","dxs_2g_2foldj","dxs_23g_2foldj",
"dxs_t_2foldep","dxs_2g_2foldep","dxs_23g_2foldep"
};

char hWname[100];
TH1F *hW[n];
for (Int_t i=0;i<n;i++) {
  sprintf(hWname,"target%i",pos[i]);  
  hW[i] = new TH1F(hWname,hWname,100,3.8,4.8);
}


TH1F *hFOM_nocut[6][3],*hFOM_Wcut[6][3];
for (Int_t i=0;i<6;i++){
for (Int_t j=0;j<3;j++){
  hFOM_nocut[i][j] = new TH1F("","",23,-362.5,-247.5);
  hFOM_Wcut[i][j] = new TH1F("","",23,-362.5,-247.5);  
}
}

for (Int_t k=0;k<n;k++){  
      char file[100];
      sprintf(file,"simu_p_11GeV_geant4_target%i_1e6.root",pos[k]);
      TChain *T = new TChain("T","T");
      T->AddFile(file);
      cout << " open file " << file << endl;   
      Int_t neve;
      T->SetBranchAddress("neve",&neve);
      T->GetEntry(T->GetEntries()-1);
      cout << " throw events " << neve << endl;
      Double_t overall=overall_NOneve/neve;

//       TH1F *hW_temp = new TH1F("hW_temp","hW_temp",100,3.8,4.8);      
      T->Project("hW_temp","W",Form("dxs_2g*weight*weight_decay*accep_je1*accep_je2*accep_e*accep_p*%f",overall));
      hW[k]=(TH1F*) hW_temp->Clone();
      
  for (Int_t i=0;i<6;i++){      
    for (Int_t j=0;j<3;j++){      
//       TH2F *htemp=new TH2F("htemp","htemp",80,0,40,110,0,11);
      T->Project("htemp1","p_e:theta_e",Form(weight_nocut[i][j],overall));
      T->Project("htemp2","p_e:theta_e",Form(weight_Wcut[i][j],overall));      
//       cout << htemp->Integral() << endl;
      double FOM_nocut=htemp1->Integral();
      double FOM_Wcut=htemp2->Integral();      
      double targetpos=-pos[k];
      hFOM_nocut[i][j]->Fill(targetpos,FOM_nocut);
      hFOM_Wcut[i][j]->Fill(targetpos,FOM_Wcut);      
      if (k==5) cout << FOM_nocut << " " << FOM_Wcut << endl;
    }
  }
}


TCanvas *c_FOM_nocut = new TCanvas("FOM_nocut","FOM_nocut",1000,900);
c_FOM_nocut->Divide(3,6);
for (Int_t i=0;i<6;i++){
 for (Int_t j=0;j<3;j++){
  c_FOM_nocut->cd(i*3+j+1);
  hFOM_nocut[i][j]->SetMarkerStyle(8);    
  hFOM_nocut[i][j]->Draw("P");
  hFOM_nocut[i][j]->SetTitle(hstname[i][j]);
 }
}
c_FOM_nocut->SaveAs("FOM_nocut.png");
c_FOM_nocut->SaveAs("FOM_nocut.C");

TCanvas *c_FOM_Wcut = new TCanvas("FOM_Wcut","FOM_Wcut",1000,900);
c_FOM_Wcut->Divide(3,6);
for (Int_t i=0;i<6;i++){
 for (Int_t j=0;j<3;j++){
  c_FOM_Wcut->cd(i*3+j+1);
  hFOM_Wcut[i][j]->SetMarkerStyle(8);  
  hFOM_Wcut[i][j]->Draw("P");
  hFOM_Wcut[i][j]->SetTitle(hstname[i][j]);  
 }
}
c_FOM_Wcut->SaveAs("FOM_Wcut.png");
c_FOM_Wcut->SaveAs("FOM_Wcut.C");

TCanvas *c_W = new TCanvas("W","W",1800,800);
c_W->Divide(5,3);
for (Int_t k=0;k<n;k++){
c_W->cd(k+1);
gPad->SetLogz();
hW[k]->Draw();
TLine *lpos0;
lpos0=new TLine(4.25,0,4.25,1E6);
lpos0->SetLineColor(kRed);
lpos0->SetLineWidth(5);
lpos0->Draw();
}

}