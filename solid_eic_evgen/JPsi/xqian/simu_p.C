#include "TApplication.h"
#include "Rtypes.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TString.h"
#include "TChain.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCut.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TGenPhaseSpace.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

const Double_t mass_e = 0.511e-3;
const Double_t mass_p = 0.93827;
const Double_t mass_n = 0.939566;
const Double_t mass_jpsi = 3.097;
const Double_t width_jpsi = 93.2e-6;
TLorentzVector pSpec;

Float_t t0lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);

// void simu_p(){
  // Double_t mass_e = 0.511e-3;
//   Double_t mass_p = 0.93827;
//   Double_t mass_jpsi = 3.097;
//   Double_t width_jpsi = 93.2e-6;
int main (Int_t argc, char *argv[])
{
  Double_t mass[10];
  
  Double_t Ebeam = 11.0;

  //use acceptance from geant3 simulation
//   TFile *file1 = new TFile("./acceptance/accep.root");
//   TH1F *hf = (TH1F*)file1->Get("h2");
//   TH1F *hl = (TH1F*)file1->Get("h3");  

  //use acceptance from geant4 simulation  
  TLorentzVector *pBeam = new TLorentzVector(0.,0.,Ebeam,sqrt(Ebeam*Ebeam+mass_e*mass_e));
  TLorentzVector *pTarget = new TLorentzVector(0.,0.,0.,mass_p);
  TGenPhaseSpace *gen = new TGenPhaseSpace();
  TGenPhaseSpace *gen1 = new TGenPhaseSpace();
  mass[0] = mass_e;
  mass[1] = mass_p;
  mass[2] = mass_jpsi;

  mass[3] = mass_e;
  mass[4] = mass_e;
  TLorentzVector *ps = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *ps1 = new TLorentzVector(0.,0.,0.,0.);

  TLorentzVector *pq = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *pt = new TLorentzVector(0.,0.,0.,0.);
  *ps = *pBeam + *pTarget;
  
  TFile *file = new TFile("temp_p.root","RECREATE");
  TTree *T = new TTree("T","T");
  T->SetDirectory(file);
  
  hf->SetDirectory(file);
  hl->SetDirectory(file);
  file1->Close();
  

  Double_t weight,weight1;
  Double_t Q2,t;
  
  Double_t p_e,theta_e,phi_e;
  Double_t p_p,theta_p,phi_p;
  Double_t p_jpsi, theta_jpsi,phi_jpsi;
  
  Double_t p_je1, theta_je1, phi_je1;
  Double_t p_je2, theta_je2, phi_je2;
  
  T->Branch("weight",&weight,"data/D");
  T->Branch("weight1",&weight1,"data/D");

  T->Branch("Q2",&Q2,"data/D");
  T->Branch("t",&t,"data/D");
  
  T->Branch("p_e",&p_e,"data/D");
  T->Branch("theta_e",&theta_e,"data/D");
  T->Branch("phi_e",&phi_e,"data/D");

  T->Branch("p_p",&p_p,"data/D");
  T->Branch("theta_p",&theta_p,"data/D");
  T->Branch("phi_p",&phi_p,"data/D");

  T->Branch("p_jpsi",&p_jpsi,"data/D");
  T->Branch("theta_jpsi",&theta_jpsi,"data/D");
  T->Branch("phi_jpsi",&phi_jpsi,"data/D");

  T->Branch("p_je1",&p_je1,"data/D");
  T->Branch("theta_je1",&theta_je1,"data/D");
  T->Branch("phi_je1",&phi_je1,"data/D");

  T->Branch("p_je2",&p_je2,"data/D");
  T->Branch("theta_je2",&theta_je2,"data/D");
  T->Branch("phi_je2",&phi_je2,"data/D");

  // calculation formular
  Double_t Gamma;
  Double_t epsilon,Keq,W,q,theta_q;
  Double_t J,R,theta_cm,r;
  
  T->Branch("Gamma",&Gamma,"data/D");
  T->Branch("epsilon",&epsilon,"data/D");
  T->Branch("Keq",&Keq,"data/D");
  T->Branch("W",&W,"data/D");
  T->Branch("q",&q,"data/D");
  T->Branch("theta_q",&theta_q,"data/D");
  T->Branch("J",&J,"data/D");
  T->Branch("R",&R,"data/D");
  T->Branch("tehta_cm",&theta_cm,"data/D");
  T->Branch("r",&r,"data/D");

  Double_t dxs;
  T->Branch("dxs",&dxs,"data/D");
  Double_t tmin;
  T->Branch("tmin",&tmin,"data/D");
    
  Double_t accep_e,accep_je1,accep_je2,accep_p;
  T->Branch("accep_e",&accep_e,"data/D");
  T->Branch("accep_je1",&accep_je1,"data/D");
  T->Branch("accep_je2",&accep_je2,"data/D");
  T->Branch("accep_p",&accep_p,"data/D");
  

  Double_t A = 0.94;
  Double_t b = -0.97;
  Double_t alpha = 1./137.;
  Double_t a = 2.164;
  Double_t n = 2.131;
  
  Double_t aa,bb,cc;

  


  for (Int_t i=0;i!=20000000;i++){
    if (ps->M() > mass[0] + mass[1] + mass[2]){
      gen->SetDecay(*ps,3,mass);
      weight = gen->Generate();
      
      TLorentzVector *p4_ep = gen->GetDecay(0);
      TLorentzVector *p4_recoil = gen->GetDecay(1); 
      TLorentzVector *p4_jpsi = gen->GetDecay(2); 
      
      gen1->SetDecay(*p4_jpsi,2,&mass[3]);
      weight1 = gen1->Generate();
      
      TLorentzVector *p4_je1 = gen1->GetDecay(0);
      TLorentzVector *p4_je2 = gen1->GetDecay(1);
      *pq = *pBeam - *p4_ep;
      Q2 = -pq->M2();
      *pt = *p4_recoil - *pTarget;
      t = -pt->M2();
      
      p_e = p4_ep->P();
      theta_e = p4_ep->Theta()/3.1415926*180.;
      phi_e = p4_ep->Phi()/3.1415926*180.;
      
      p_p = p4_recoil->P();
      theta_p = p4_recoil->Theta()/3.1415926*180.;
      phi_p = p4_recoil->Phi()/3.1415926*180.;
      
      p_jpsi = p4_jpsi->P();
      theta_jpsi = p4_jpsi->Theta()/3.1415926*180.;
      phi_jpsi = p4_jpsi->Phi()/3.1415926*180.;
      
      p_je1 = p4_je1->P();
      theta_je1 = p4_je1->Theta()/3.1415926*180.;
      phi_je1 = p4_je1->Phi()/3.1415926*180.;
      
      p_je2 = p4_je2->P();
      theta_je2 = p4_je2->Theta()/3.1415926*180.;
      phi_je2 = p4_je2->Phi()/3.1415926*180.;
      
      
      R = pow((a*mass_jpsi*mass_jpsi-Q2)/(a*mass_jpsi*mass_jpsi),n) -1;
      theta_q = pq->Theta()/3.1415926*180.;
      q = pq->P();
      W = sqrt(pow(mass_p + pq->E(),2)-pow(pq->P(),2));
      Keq = (W*W-mass_p*mass_p)/2./mass_p;
      epsilon = 1./(1+2*q*q/Q2*pow(tan(theta_e/180.*3.1415926/2.),2));
      Gamma = alpha/2./3.1415926/3.1415926*p_e/Ebeam*Keq/Q2/(1.-epsilon);
      r= epsilon*R/(1.+epsilon*R);
      J = fabs((pq->E()+mass_p-q*p4_recoil->E()/p_p*(cos(theta_q/180.*3.1415926)*cos(theta_p/180.*3.1415926)+sin(theta_q/180.*3.1415926)*sin(theta_p/180.*3.1415926)*sin((phi_p-phi_e+180.)/180.*3.1415926))*tan(theta_q/180.*3.1415926))
	       /(2.*mass_p*q*p_p*(cos(theta_q*3.1415926/180.)*tan(theta_p/180.*3.1415926)-sin(theta_q/180.*3.1415926)*sin((phi_p-phi_e+180.)/180.*3.1415926))));
      TVector3 beta = p4_jpsi->Vect();
      beta *= -1./p4_jpsi->E();
      p4_je1->Boost(beta);
      theta_cm = p4_je1->Theta()/3.1415926*180.;
      //theta_cm = 0.;
      
      weight1 = 3./2.*(1-r + (3*r-1)*pow(cos(theta_cm/180.*3.1415926),2));
      
      // aa = 2*mass_p*mass_p + 2.*(Ebeam-p_e)*mass_p-mass_jpsi*mass_jpsi-Q2;
      //     bb = 2*(mass_p+Ebeam-p_e);
      //     cc = 2*q;
//     aa = (2.*aa*bb-sqrt(4*aa*aa*bb*bb-4.*(bb*bb-cc*cc)*(aa*aa+cc*cc*mass_p*mass_p)))/2./(bb*bb-cc*cc);
//     tmin = fabs(2.*mass_p*mass_p-2.*mass_p*aa);
      
      *ps1 = *pq + *pTarget;
      tmin = -1*t0lim(-sqrt(Q2),mass_p, mass_jpsi, mass_p, ps1->M2());
      //cout << tmin << " " <<  t0lim(-sqrt(Q2),mass_p, mass_jpsi, mass_p, ps1->M2()) << endl;
      
      
      dxs = J/2./3.1415926*Gamma*A*exp(b*(t-tmin));
      

      Int_t tmp_x,tmp_y;
      tmp_x = 10.*(theta_e)+1;
      tmp_y = 100.*(p_e)/8.+1;
      if (tmp_x>310) tmp_x = 310;
      if (tmp_y>100) tmp_y = 100;
      
      accep_e = hf->GetBinContent(tmp_x,tmp_y);
      
      tmp_x = 10.*(theta_je1)+1;
      tmp_y = 100.*(p_je1)/8.+1;
      if (tmp_x>310) tmp_x = 310;
      if (tmp_y>100) tmp_y = 100;
      
      accep_je1 = hf->GetBinContent(tmp_x,tmp_y) + hl->GetBinContent(tmp_x,tmp_y);

      tmp_x = 10.*(theta_je2)+1;
      tmp_y = 100.*(p_je2)/8.+1;
      if (tmp_x>310) tmp_x = 310;
      if (tmp_y>100) tmp_y = 100;
      
      accep_je2 = hf->GetBinContent(tmp_x,tmp_y) + hl->GetBinContent(tmp_x,tmp_y);

      tmp_x = 10.*(theta_p)+1;
      tmp_y = 100.*(p_p)/8.+1;
      if (tmp_x>310) tmp_x = 310;
      if (tmp_y>100) tmp_y = 100;
      
      accep_p = hf->GetBinContent(tmp_x,tmp_y);

      if (accep_e>0.&&accep_je1>0.&&accep_je2>0.&&accep_p>0.)
	T->Fill();
    }
  }
  file->Write();
  file->Close();
}


Float_t t0lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s)
{
  Float_t t1,t2,t3,t4;
  
  if (m1>=0){
     t1 = (m1*m1 - m3*m3 - m2*m2 + m4*m4)/(2*sqrt(s));
     t2 = (s + m1*m1 - m2*m2)/(2*sqrt(s));    //E1cm
    //if (t2 < 0.) {return 1.;}
    t2 = sqrt(t2*t2 - m1*m1);                        //p1cm
     t3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));    //E3cm
    //if (t3 < 0.) {return 1.;}
    t3 = sqrt(t3*t3 - m3*m3);                        //p3cm
     t4 = t2 - t3;                            //p1cm-p3cm
  }else{
     t1 = (-m1*m1 - m3*m3 - m2*m2 + m4*m4)/(2*sqrt(s));
     t2 = (s - m1*m1 - m2*m2)/(2*sqrt(s));    //E1cm
    //if (t2 < 0.) {return 1.;}
    t2 = sqrt(t2*t2 + m1*m1);                        //p1cm
     t3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));    //E3cm
    //if (t3 < 0.) {return 1.;}
    t3 = sqrt(t3*t3 - m3*m3);                        //p3cm
     t4 = t2 - t3;                            //p1cm-p3cm
  }
  return  t1*t1 - t4*t4;
}
