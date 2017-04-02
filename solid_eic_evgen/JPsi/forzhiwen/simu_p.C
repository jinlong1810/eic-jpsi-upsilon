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
#include "TRandom.h"

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
Float_t t1lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s);
Double_t fun_2g(Double_t x, Double_t t);
Double_t fun_23g(Double_t x, Double_t t);

int main (Int_t argc, char *argv[])
{
  gRandom->SetSeed(0);
  if (argc==1){
    cout << "./simu_p -a[nevents] -b[save_flag] -c[Ebeam] -d[target_flag]" << endl;
    cout << "save_flag   = 0: only save the events passing acceptance and kinematics" << endl;
    cout << "save_flag   = 1: save the events passing kinematics" << endl;
    cout << "target_flag = 1: Proton" << endl;
    cout << "target_flag = 2: Deuteron" << endl;
  }else{
    Int_t nevents=500000;
    Int_t save_flag=0; 
    Double_t Ebeam=11.0; // Beam Energy
    Int_t target_flag = 1; // proton
    
    for(Int_t i = 1; i != argc; i++){
      switch(argv[i][1]){
      case 'a':
        nevents= atoi(&argv[i][2]); 
	break;
      case 'b':
	save_flag = atoi(&argv[i][2]); 
	break;
      case 'c':
	Ebeam = atof(&argv[i][2]); 
	break;
      case 'd':
	target_flag = atoi(&argv[i][2]); 
	break;	
      default:
	cout << "Warning!!!! Unknown option: " << &argv[i][1] << endl;
	break;
      }
    }

    Double_t mass[10];
    
    TFile *file1 = new TFile("./acceptance/accep.root");
    TH1F *hf = (TH1F*)file1->Get("h2");
    TH1F *hl = (TH1F*)file1->Get("h3");

    TLorentzVector *pBeam = new TLorentzVector(0.,0.,Ebeam,sqrt(Ebeam*Ebeam+mass_e*mass_e));
    TLorentzVector *pTarget = new TLorentzVector(0.,0.,0.,mass_p);
    
    TGenPhaseSpace *gen1 = new TGenPhaseSpace();
    mass[0] = mass_e;
    mass[1] = mass_e;
    
    TLorentzVector *ps = new TLorentzVector(0.,0.,0.,0.);
    *ps = *pBeam + *pTarget;
    TLorentzVector *pq = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *pt = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *ps1 = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *ps2 = new TLorentzVector(0.,0.,0.,0.);

    TLorentzVector *p4_ep = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_recoil = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_jpsi = new TLorentzVector(0.,0.,0.,0.);


    TLorentzVector *p4_recoil_tmp = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_q_tmp = new TLorentzVector(0.,0.,0.,0.);

    TLorentzVector *p4_je1 = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_je2 = new TLorentzVector(0.,0.,0.,0.);
        
    TFile *file = new TFile("temp_p.root","RECREATE");
    TTree *T = new TTree("T","T");
    T->SetDirectory(file);
    
    hf->SetDirectory(file);
    hl->SetDirectory(file);
    file1->Close();

    Double_t weight_decay;
    Double_t Q2,t;
    
    Double_t p_e,theta_e,phi_e;
    Double_t p_p,theta_p,phi_p;
    Double_t p_jpsi, theta_jpsi,phi_jpsi;
    
    Double_t p_je1, theta_je1, phi_je1;
    Double_t p_je2, theta_je2, phi_je2;
    
    T->Branch("weight_decay",&weight_decay,"data/D");
    Double_t weight;
    T->Branch("weight",&weight,"data/D");
    
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

    Int_t neve=0,neve1=0;
    T->Branch("neve",&neve,"data/I");
    
    // calculation formular
    Double_t Gamma;
    Double_t epsilon,Keq,W,q,theta_q;
    Double_t J,R,theta_cm,r,phi_cm;
    
    T->Branch("Gamma",&Gamma,"data/D");
    T->Branch("epsilon",&epsilon,"data/D");
    T->Branch("Keq",&Keq,"data/D");
    T->Branch("W",&W,"data/D");
    T->Branch("q",&q,"data/D");
    T->Branch("theta_q",&theta_q,"data/D");
    T->Branch("J",&J,"data/D");
    T->Branch("R",&R,"data/D");
    T->Branch("theta_cm",&theta_cm,"data/D");
    T->Branch("phi_cm",&phi_cm,"data/D");
    T->Branch("r",&r,"data/D");
    
    Double_t dxs,dxs_2g,dxs_23g;
    T->Branch("dxs",&dxs,"data/D");
    T->Branch("dxs_2g",&dxs_2g,"data/D");
    T->Branch("dxs_23g",&dxs_23g,"data/D");
    
    Double_t tmin,tmax;
    T->Branch("tmin",&tmin,"data/D");
    T->Branch("tmax",&tmax,"data/D");
    
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
    
    Int_t flag_kin1,flag_kin2,flag_accep;

    
    Int_t qflag = 1;
    // start to generate particles
    while(qflag){
      //set the flag for kinematics 
      flag_kin1 = 0;
      flag_kin2 = 0;
      

      //sample electron's angle and momentum
      p_e = gRandom->Uniform(0.5,3.0);
      theta_e = acos(gRandom->Uniform(0.85,cos(8./180.*3.1415926)));
      phi_e = gRandom->Uniform(0.,2.*3.141596);
      
      p4_ep->SetPxPyPzE(p_e*sin(theta_e)*cos(phi_e),p_e*sin(theta_e)*sin(phi_e),p_e*cos(theta_e),sqrt(p_e*p_e+mass_e*mass_e));
      *pq = *pBeam - *p4_ep;
      
      Q2 = -pq->M2();
      *ps1 = *ps - *p4_ep;

      // ps1->SetPxPyPzE(0,0,11,sqrt(11*11+mass_p*mass_p)+mass_jpsi);
      
      //judge whether the event pass the kinematics
      if (ps1->M2() > pow(mass_jpsi+mass_p,2)){
	flag_kin1 = 1;
	//sample proton solid angle;
	theta_p = acos(gRandom->Uniform(0.85,cos(8./180.*3.1415926)));
	phi_p = gRandom->Uniform(0.,2.*3.141596);
	
	TVector3 p3_p(sin(theta_p)*cos(phi_p),sin(theta_p)*sin(phi_p),cos(theta_p));
	Double_t aa = (ps1->M2()+mass_p*mass_p-mass_jpsi*mass_jpsi)/2.;
	Double_t bb = ps1->E();
	TVector3 p3_s(ps1->Px(),ps1->Py(),ps1->Pz());
	Double_t cc = -p3_p.Dot(p3_s);
	Double_t sol[2];
	
	if (pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*mass_p*mass_p-aa*aa)>=0.){
	  sol[0] = (-2*aa*cc+sqrt(pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*mass_p*mass_p-aa*aa)))/2./(bb*bb-cc*cc);
	  sol[1] = (-2*aa*cc-sqrt(pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*mass_p*mass_p-aa*aa)))/2./(bb*bb-cc*cc);

	  if (sol[0]>=0||sol[1]>=0){
	    flag_kin2 = 1;
	    if (sol[0]>=0&&sol[1]<0){
	      weight = 1;
	    }else if (sol[1]>=0&&sol[0]<0){
	      weight = 1;
	    }else if (sol[0]>=0&&sol[1]>=0){
	      weight = 0.5;
	    }

	    //	    cout << sol[0] << " " << sol[1] << endl;
	    
	    Double_t theta_ps = theta_p;
	    Double_t phi_ps = phi_p;

	    for (Int_t j=0;j!=2;j++){
	      flag_accep = 0.;
	      p_p = sol[j];
	      if (p_p>0){
		p4_recoil->SetPxPyPzE(p_p*sin(theta_ps)*cos(phi_ps),p_p*sin(theta_ps)*sin(phi_ps),p_p*cos(theta_ps),sqrt(p_p*p_p+mass_p*mass_p));
		*p4_jpsi = *ps1-*p4_recoil;
		
		//	cout << p4_recoil->Theta() << " " << theta_ps << endl;

		if (p4_jpsi->M()>0.){
		  gen1->SetDecay(*p4_jpsi,2,&mass[0]);
		  weight_decay = gen1->Generate();
		  p4_je1 = gen1->GetDecay(0);
		  p4_je2 = gen1->GetDecay(1);
		  
		  theta_e = p4_ep->Theta()/3.1415926*180.;
		  phi_e = p4_ep->Phi()/3.1415926*180.;
		  
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
		  
		  *pt = *p4_recoil - *pTarget;
		  t = -pt->M2();
		  
		  R = pow((a*mass_jpsi*mass_jpsi-Q2)/(a*mass_jpsi*mass_jpsi),n) -1;
		  theta_q = pq->Theta()/3.1415926*180.;
		  q = pq->P();
		  W = sqrt(pow(mass_p + pq->E(),2)-pow(pq->P(),2));
		  Keq = (W*W-mass_p*mass_p)/2./mass_p;
		  epsilon = 1./(1+2*q*q/Q2*pow(tan(theta_e/180.*3.1415926/2.),2));
		  Gamma = alpha/2./3.1415926/3.1415926*p_e/Ebeam*Keq/Q2/(1.-epsilon);
		  r= epsilon*R/(1.+epsilon*R);
		  
		  // J = fabs((pq->E()+mass_p-q*p4_recoil->E()/p_p*(cos(theta_q/180.*3.1415926)*cos(theta_p/180.*3.1415926)+sin(theta_q/180.*3.1415926)*sin(theta_p/180.*3.1415926)*sin((phi_p-phi_e+180.)/180.*3.1415926))*tan(theta_q/180.*3.1415926))
		  //	   /(2.*mass_p*q*p_p*(cos(theta_q*3.1415926/180.)*tan(theta_p/180.*3.1415926)-sin(theta_q/180.*3.1415926)*sin((phi_p-phi_e+180.)/180.*3.1415926))));
		  // cout << aa - bb*sqrt(p_p*p_p+mass_p*mass_p)-cc*p_p << endl;

		  Double_t dEpdcp=fabs((ps1->Pz()-(ps1->Px()*cos(phi_p/180.*3.1415926)+ps1->Py()*sin(phi_p/180.*3.1415926))*cos(theta_p/180.*3.1415926)/sin(theta_p/180.*3.1415926))*p_p/(bb+cc*sqrt(p_p*p_p+mass_p*mass_p)/p_p));
		  //cout << J << " " << 1./mass_p/2./dEpdcp<< " " << dEpdcp << " " << p_p << " " << bb+cc*sqrt(p_p*p_p+mass_p*mass_p)/p_p << endl;
		  J = 2 * mass_p * dEpdcp;

		  //go to JPsi at rest frame
		  TVector3 beta = p4_jpsi->Vect();
		  beta *= -1./p4_jpsi->E();
		  p4_je1->Boost(beta);
		  //get recoil proton in the same frame
		  *p4_recoil_tmp = *p4_recoil;
		  p4_recoil_tmp->Boost(beta);

		  *p4_q_tmp = *pq;
		  p4_q_tmp->Boost(beta);
		  // calculate the theta angle between these two
		  TVector3 a1 = p4_je1->Vect();
		  TVector3 a2 = p4_recoil_tmp->Vect();
		  TVector3 a3 = p4_q_tmp->Vect();

		  a1.RotateUz(a2.Unit());
		  a3.RotateUz(a2.Unit());
		  theta_cm = a1.Theta()/3.1415926*180.;//p4_je1->Theta()/3.1415926*180.;
		  phi_cm = a1.Phi()/3.1415926*180.-a3.Phi()/3.1415926*180.;
		  if (phi_cm<0.) phi_cm+=360;
		  if (phi_cm>360) phi_cm-=360;
		  // phi_cm = p4_je1->Phi()/3.1415926*180.;
		  //theta_cm = 0.;
		  
		  weight_decay = 3./2.*(1-r + (3*r-1)*pow(cos(theta_cm/180.*3.1415926),2));
		  
		  // calculate tmin
		  *ps2 = *pq + *pTarget;
		  tmin = -1*t0lim(-sqrt(Q2),mass_p, mass_jpsi, mass_p, ps2->M2());
		  tmax = -1*t1lim(-sqrt(Q2),mass_p, mass_jpsi, mass_p, ps2->M2());
		  
		  

		  dxs     = J/2./3.1415926*Gamma*A*exp(b*(t));
		  dxs_2g  = J/2./3.1415926*Gamma*fun_2g(W,t);
		  dxs_23g = J/2./3.1415926*Gamma*fun_23g(W,t);
		  
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
		  
		  if (accep_e>0.&&accep_je1>0.&&accep_je2>0.) flag_accep=1;
		  if (flag_accep>0.&&flag_kin1>0.&&flag_kin2>0.){
		    T->Fill();
		    if (neve1%1000==0) cout << neve1 << endl;
		    neve1++;
		    if (neve1 > nevents){
		      qflag = 0;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      neve++;
    }

    file->Write();
    file->Close();
  }
  return 0;
}

Double_t fun_2g(Double_t x, Double_t t){
  Double_t N2g = 7.5671e3;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t M = 3.097;
  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff;
  return result;
}

Double_t fun_23g(Double_t x,Double_t t){
  Double_t N2g = 6.499e3;
  Double_t N3g = 2.894e3;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t M = 3.097;
 
  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff
    +N3g*v/R/R/R/R/M/M/M/M*pow(1-xp,0)*ff;
  return result;
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

Float_t t1lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s)
{
  Float_t t1,t2,t3,t4;
  if (m1>=0){
    t1 = (m1*m1 - m3*m3 - m2*m2 + m4*m4)/(2*sqrt(s));
    t2 = (s + m1*m1 - m2*m2)/(2*sqrt(s));
    //if (t2 < 0.) {return 1.;}
    t2 = sqrt(t2*t2 - m1*m1);
    t3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
    //if (t3 < 0.) {return 1.;}                      
    t3 = sqrt(t3*t3 - m3*m3);                         
    t4 = t2 + t3;                            //p1cm+p3cm 
  }else{
     t1 = (-m1*m1 - m3*m3 - m2*m2 + m4*m4)/(2*sqrt(s));
    t2 = (s - m1*m1 - m2*m2)/(2*sqrt(s));
    //if (t2 < 0.) {return 1.;}
    t2 = sqrt(t2*t2 - m1*m1);
    t3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));
    //if (t3 < 0.) {return 1.;}                      
    t3 = sqrt(t3*t3 - m3*m3);                         
    t4 = t2 + t3;                            //p1cm+p3cm 
  }
  return  t1*t1 - t4*t4;
}
