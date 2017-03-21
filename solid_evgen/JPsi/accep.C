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
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void accep(string input_filename,string acceptance_filename)
{
gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

const double DEG=180./3.1415926;

TFile *file1;
TH2F *hf,*hl;
file1 = new TFile(acceptance_filename.c_str());    
if (file1->IsZombie()) {
  cout << "Error opening acceptance file" << acceptance_filename << endl;       
  return;
}
else {
  cout << "open acceptance file " << acceptance_filename << endl;    
  hf = (TH2F*)file1->Get("acceptance_ThetaP_forwardangle");
  hl = (TH2F*)file1->Get("acceptance_ThetaP_largeangle");	
//   hf = (TH2F*)file1->Get("acceptance_forwardangle");
//   hl = (TH2F*)file1->Get("acceptance_largeangle");	  
}	

char input_filename_new[200];
sprintf(input_filename_new,"%s_accep.root",input_filename.substr(0,input_filename.rfind(".")).c_str());

gSystem->Exec(Form("cp -f %s %s",input_filename.c_str(),input_filename_new));
cout << "copy " << input_filename << " to " << input_filename_new << endl;

TFile *file=new TFile(input_filename_new,"update");
if (file->IsZombie()) {
  cout << "Error opening input file" << input_filename_new << endl;
  exit(-1);
}
else cout << "open input file " << input_filename_new << endl;

TTree *T = (TTree*) file->Get("T");

Double_t accep_e_1,accep_je1_1,accep_je2_1,accep_p_1;
Double_t accep_e_2,accep_je1_2,accep_je2_2,accep_p_2;
TBranch *bnew1 = T->Branch("accep_e_1",&accep_e_1,"data/D");
TBranch *bnew2 = T->Branch("accep_je1_1",&accep_je1_1,"data/D");
TBranch *bnew3 = T->Branch("accep_je2_1",&accep_je2_1,"data/D");
TBranch *bnew4 = T->Branch("accep_p_1",&accep_p_1,"data/D");
TBranch *bnew5 = T->Branch("accep_e_2",&accep_e_2,"data/D");
TBranch *bnew6 = T->Branch("accep_je1_2",&accep_je1_2,"data/D");
TBranch *bnew7 = T->Branch("accep_je2_2",&accep_je2_2,"data/D");
TBranch *bnew8 = T->Branch("accep_p_2",&accep_p_2,"data/D");

    Double_t p_e,theta_e,phi_e;
    Double_t p_p,theta_p,phi_p;
    Double_t p_jpsi, theta_jpsi,phi_jpsi;
    
    Double_t p_je1, theta_je1, phi_je1;
    Double_t p_je2, theta_je2, phi_je2;

T->SetBranchStatus("*",1);
T->SetBranchAddress("p_e",&p_e);
T->SetBranchAddress("theta_e",&theta_e);
T->SetBranchAddress("phi_e",&phi_e);
T->SetBranchAddress("p_p",&p_p);
T->SetBranchAddress("theta_p",&theta_p);
T->SetBranchAddress("phi_p",&phi_p);
T->SetBranchAddress("p_jpsi",&p_jpsi);
T->SetBranchAddress("theta_jpsi",&theta_jpsi);
T->SetBranchAddress("phi_jpsi",&phi_jpsi);
T->SetBranchAddress("p_je1",&p_je1);
T->SetBranchAddress("theta_je1",&theta_je1);
T->SetBranchAddress("phi_je1",&phi_je1);
T->SetBranchAddress("p_je2",&p_je2);
T->SetBranchAddress("theta_je2",&theta_je2);
T->SetBranchAddress("phi_je2",&phi_je2);

//  loop on entries
Int_t nentries = T->GetEntries();
for (Int_t i=0;i<nentries;i++) {
  cout << i << "\r";
  T->GetEvent(i);

  int bin,binx,biny,binz;

  bin=hf->FindBin(theta_e,p_e);  
  hf->GetBinXYZ(bin,binx,biny,binz);
  accep_e_1 = hf->GetBinContent(binx,biny);
  accep_e_2 = hl->GetBinContent(binx,biny);

  bin=hf->FindBin(theta_je1,p_je1);
  hf->GetBinXYZ(bin,binx,biny,binz);
  accep_je1_1 = hf->GetBinContent(binx,biny);  
  accep_je1_2 = hl->GetBinContent(binx,biny);
  
  bin=hf->FindBin(theta_je2,p_je2);
  hf->GetBinXYZ(bin,binx,biny,binz);
  accep_je2_1 = hf->GetBinContent(binx,biny);  
  accep_je2_2 = hl->GetBinContent(binx,biny);

  bin=hf->FindBin(theta_p,p_p);  
  hf->GetBinXYZ(bin,binx,biny,binz);
  accep_p_1 = hf->GetBinContent(binx,biny);
  accep_p_2 = hl->GetBinContent(binx,biny);  
  
  //HallD
  // 		  if (3 < theta_je1 && theta_je1 <140)  accep_je1 = 1.;
// 		  else accep_je1 = 0;
// 		  if (3 < theta_je2 && theta_je2 <140)  accep_je2 = 1.; 		    
// 		  else accep_je2 = 0;		  
// 		  if (3 < theta_p && theta_p <150)  accep_p = 1.;
// 		  else accep_p = 0;		  
// // 		  if (8 < theta_je1 && theta_je1 <25)  accep_je1 = 1.;
// // 		  else accep_je1 = 0;
// // 		  if (8 < theta_je2 && theta_je2 <25)  accep_je2 = 1.; 		    
// // 		  else accep_je2 = 0;		  
// // 		  if (8 < theta_p && theta_p <25)  accep_p = 1.;
// // 		  else accep_p = 0;
		  
   bnew1->Fill();
   bnew2->Fill();
   bnew3->Fill();
   bnew4->Fill();
   bnew5->Fill();
   bnew6->Fill(); 
   bnew7->Fill();
   bnew8->Fill();    
}
T->Write();
        
    file->Write();
    file->Close();
}