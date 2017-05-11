#include "simu_p.h"
#include "globals.h"

#include <TF1.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom.h>
#include <TVector3.h>

using namespace std;

Simulator::Simulator(){

  _file_out = NULL;
  _tree_out = NULL;

  seed = 0;

  nevents=1000000; //number of events

  Ebeam_lab=11.0; // Electron Beam Energy in labortory frame
  Etarget_lab=0.0; // Proton Beam Energy in labortory frame

  output_root_file = "output.root";

  meson_type=("jpsi");

  Is_e=false;
  Is_g=false;

  Gbeam_min = 7.5;

  neve = 0;

  /* set all event-wise variables to 0 */
  reset_event_variables();

}


int Simulator::set_process_type( TString type )
{
  if (type=="e")
    {
      cout << "Select process type: Electroproduction" << endl;
      Is_e=true;
    }
  else if (type=="g")
    {
      cout << "Select process type: Photoproduction" << endl;
      Is_g=true;
    }
  else
    {
      cout << "wrong type" << endl;
      return 1;
    }

  return 0;
}


int Simulator::set_meson_type( TString newtype )
{
  meson_type = newtype;

  /* Check selected meson type */
  if ( meson_type == "jpsi" )
    cout << "Select meson production: J/Psi" << endl;
  else if ( meson_type == "upsilon" )
    cout << "Select meson production: Upsilon" << endl;
  else
    {
      cout << "ERROR: Unknown meson type " << meson_type << endl;
      return(1);
    }

  return 0;
}


int Simulator::run ()
{
  cout << "Running Simulator..." << endl;

  gRandom->SetSeed(seed);

  init();
  process_event();

  /* run everything */

  /* Transition from 'laboratory frame' to 'target rest' frame*/
  TLorentzVector *pBeam_lab = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *pTarget_lab = new TLorentzVector(0.,0.,0.,0.);

  pBeam_lab->SetPxPyPzE(0.,0.,Ebeam_lab,sqrt(Ebeam_lab*Ebeam_lab+simglobals::mass_e*simglobals::mass_e));
  pTarget_lab->SetPxPyPzE(0.,0.,Etarget_lab,sqrt(Etarget_lab*Etarget_lab+simglobals::mass_p*simglobals::mass_p));

  TLorentzVector *pBeam_prest = (TLorentzVector*)pBeam_lab->Clone(); //new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *pTarget_prest = (TLorentzVector*)pTarget_lab->Clone(); //new TLorentzVector(0.,0.,0.,0.);

  /* determine beta to move to 'proton at rest' reference frame */
  TVector3 beta_lab_prest = pTarget_lab->Vect();
  beta_lab_prest *= -1./pTarget_lab->E();
  pBeam_prest->Boost(beta_lab_prest);
  pTarget_prest->Boost(beta_lab_prest);

  Double_t Ebeam = pBeam_prest->E();
  /* pTarget_prest->Vect().Mag() may give very small but non-0 number (rounding etc.), so force Etarget to 0 */
  Double_t Etarget = 0;
  pTarget_prest->SetPxPyPzE(0.,0.,0.,simglobals::mass_p);

  /* Monitoring output statements */
  cout << "Energies in LABORATORY  FRAME: " << pBeam_lab->E() << " GeV (e) -> " << pTarget_lab->E() << " GeV (p) " << endl;
  cout << "Energies in TARGET REST FRAME: " << pBeam_prest->E() << " GeV (e) -> " << pTarget_prest->E() << " GeV (p) " << endl;

  cout << "********************" << endl;

  cout << "Lab frame electron 4-vector before boosts: ("
       << pBeam_lab->Px() << ", "
       << pBeam_lab->Py() << ", "
       << pBeam_lab->Pz() << ", "
       << pBeam_lab->E() << ")" << endl;

    cout << "Proton Rest frame electron 4-vector: ("
	 << pBeam_prest->Px() << ", "
	 << pBeam_prest->Py() << ", "
	 << pBeam_prest->Pz() << ", "
	 << pBeam_prest->E() << ")" << endl;
 
    TLorentzVector *pBeam_lab_xcheck = new TLorentzVector(0.,0.,0.,0.);
    *pBeam_lab_xcheck = *pBeam_prest;
    pBeam_lab_xcheck->Boost(-beta_lab_prest);

    cout << "Lab frame electron 4-vector after boosts: ("
	 << pBeam_lab_xcheck->Px() << ", "
	 << pBeam_lab_xcheck->Py() << ", "
	 << pBeam_lab_xcheck->Pz() << ", "
	 << pBeam_lab_xcheck->E() << ")" << endl;

    cout << "********************" << endl;

    cout << "Lab frame proton 4-vector before boosts: ("
	 << pTarget_lab->Px() << ", "
	 << pTarget_lab->Py() << ", "
	 << pTarget_lab->Pz() << ", "
	 << pTarget_lab->E() << ")" << endl;

    cout << "Proton Rest frame proton 4-vector: ("
	 << pTarget_prest->Px() << ", "
	 << pTarget_prest->Py() << ", "
	 << pTarget_prest->Pz() << ", "
	 << pTarget_prest->E() << ")" << endl;
 
    TLorentzVector *pTarget_lab_xcheck = new TLorentzVector(0.,0.,0.,0.);
    *pTarget_lab_xcheck = *pTarget_prest;
    pTarget_lab_xcheck->Boost(-beta_lab_prest);

    cout << "Lab frame proton 4-vector after boosts: ("
	 << pTarget_lab_xcheck->Px() << ", "
	 << pTarget_lab_xcheck->Py() << ", "
	 << pTarget_lab_xcheck->Pz() << ", "
	 << pTarget_lab_xcheck->E() << ")" << endl;

    cout << "********************" << endl;
    /* END Monitoring output statements */

    /* END transition to 'target rest' frame */

    if (Is_e)  cout << "Ebeam " << Ebeam << " GeV, Etarget " << Etarget << " GeV" << endl;
    else if (Is_g) cout << "Gbeam " << Gbeam_min << " - " << Ebeam << " GeV, Etarget " << Etarget << " GeV" << endl;

    TF1 *fbr = new TF1("fbr","[1]/(x/[0])*(4./3.-4./3.*(x/[0])+(x/[0])*(x/[0]))",Gbeam_min,Ebeam);
    fbr->SetParameter(0,Ebeam);
    fbr->SetParameter(1,0.00243);
    //brem photon produced by one electron for 100% radiator
    //a=0.000146 for a 6% radiator, according to p251 of PDG 2012 and Jixie Zhang's estimation,  it only works for thin radiator
    //a=0.00243=0.000146/0.06 for a 100% radiator

    Double_t mass[10];

    TGenPhaseSpace *gen1 = new TGenPhaseSpace();
    mass[0] = simglobals::mass_e;
    mass[1] = simglobals::mass_e;

    Double_t mass_meson = 0;
    if ( meson_type == "jpsi")
      mass_meson = simglobals::mass_jpsi;
    else if ( meson_type == "upsilon" )
      mass_meson = simglobals::mass_upsilon;

    /* kinematics 4-vectors in "proton at rest" frame */
    TLorentzVector *ps_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *pq_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *pt_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *ps1_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *ps2_prest = new TLorentzVector(0.,0.,0.,0.);

    /* 4-vectors in "proton at rest" frame */
    TLorentzVector *p4_ep_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_recoil_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_jpsi_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_je1_prest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_je2_prest = new TLorentzVector(0.,0.,0.,0.);

    /* 4-vectors in "laboratory" frame */
    TLorentzVector *p4_ep_lab = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_recoil_lab = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_jpsi_lab = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_je1_lab = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_je2_lab = new TLorentzVector(0.,0.,0.,0.);

    /* 4-vectors in "J/Psi at rest" frame */
    TLorentzVector *p4_je1_jpsirest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_recoil_jpsirest = new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *p4_q_jpsirest = new TLorentzVector(0.,0.,0.,0.);

    Int_t neve1=0;

    /* calculation formular */
    Double_t A = 0.94;
    Double_t b = -0.97;
    Double_t alpha = 1./137.;
    Double_t a = 2.164;
    Double_t n = 2.131;

    /* start to generate particles */
    Int_t qflag = 1;
    Int_t counter[4]={0,0,0,0};

    if(Is_e){

      /* choose maximum energy of scattered electron based on electron beam energy */
      Double_t p_e_max = Ebeam - 6.0;

      phasespace = (p_e_max-0.)*(4*3.1415926)*(4.*3.1415926)*(4.*3.1415926);

      while(qflag){

        //sample electron's angle and momentum
        //       p_e = gRandom->Uniform(0.5,3.0);
        p_e = gRandom->Uniform(0.,p_e_max);
        //       theta_e = acos(gRandom->Uniform(0.85,cos(8./simglobals::DEG)));
        //       theta_e = acos(gRandom->Uniform(cos(40./simglobals::DEG),cos(0./simglobals::DEG))); //random selection in solid angle need to go with cos(theta)
        theta_e = acos(gRandom->Uniform(-1,1)); //random selection in solid angle need to go with cos(theta)
	eta_e = -1*log( tan( theta_e / 2.0 ) );
        phi_e = gRandom->Uniform(0.,2.*3.1415926);

        p4_ep_prest->SetPxPyPzE(p_e*sin(theta_e)*cos(phi_e),p_e*sin(theta_e)*sin(phi_e),p_e*cos(theta_e),sqrt(p_e*p_e+simglobals::mass_e*simglobals::mass_e));

        *ps_prest = *pBeam_prest + *pTarget_prest;

        *pq_prest = *pBeam_prest - *p4_ep_prest;

        Q2 = -pq_prest->M2();
        *ps1_prest = *ps_prest - *p4_ep_prest;

        // ps1_prest->SetPxPyPzE(0,0,11,sqrt(11*11+simglobals::mass_p*simglobals::mass_p)+mass_meson);

        //judge whether the event pass the kinematics
        if (ps1_prest->M2() > pow(mass_meson+simglobals::mass_p,2)){

          //sample proton solid angle;
          //    theta_p = acos(gRandom->Uniform(0.85,cos(8./simglobals::DEG)));
          //    theta_p = acos(gRandom->Uniform(cos(40./simglobals::DEG),cos(0./simglobals::DEG))); //random selection in solid angle need to go with cos(theta)
          theta_p = acos(gRandom->Uniform(-1,1)); //random selection in solid angle need to go with cos(theta)
          phi_p = gRandom->Uniform(0.,2.*3.1415926);

          //solve recoil proton mom
          TVector3 p3_p(sin(theta_p)*cos(phi_p),sin(theta_p)*sin(phi_p),cos(theta_p));
          Double_t aa = (ps1_prest->M2()+simglobals::mass_p*simglobals::mass_p-mass_meson*mass_meson)/2.;
          Double_t bb = ps1_prest->E();
          TVector3 p3_s(ps1_prest->Px(),ps1_prest->Py(),ps1_prest->Pz());
          Double_t cc = -p3_p.Dot(p3_s);
          Double_t sol[2];

          if (pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*simglobals::mass_p*simglobals::mass_p-aa*aa)>=0.){
            sol[0] = (-2*aa*cc+sqrt(pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*simglobals::mass_p*simglobals::mass_p-aa*aa)))/2./(bb*bb-cc*cc);
            sol[1] = (-2*aa*cc-sqrt(pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*simglobals::mass_p*simglobals::mass_p-aa*aa)))/2./(bb*bb-cc*cc);

            if (sol[0]>=0||sol[1]>=0){

              if (sol[0]>=0&&sol[1]<0){
                weight = 1;
              }else if (sol[1]>=0&&sol[0]<0){
                weight = 1;
              }else if (sol[0]>=0&&sol[1]>=0){
                weight = 0.5;
              }

              Double_t theta_ps_prest = theta_p;
              Double_t phi_ps_prest = phi_p;

              for (Int_t j=0;j!=2;j++){
                p_p = sol[j];
                if (p_p>0){
                  p4_recoil_prest->SetPxPyPzE(p_p*sin(theta_ps_prest)*cos(phi_ps_prest),p_p*sin(theta_ps_prest)*sin(phi_ps_prest),p_p*cos(theta_ps_prest),sqrt(p_p*p_p+simglobals::mass_p*simglobals::mass_p));
                  *p4_jpsi_prest = *ps1_prest-*p4_recoil_prest;

                  if (p4_jpsi_prest->M()>0.){
                    gen1->SetDecay(*p4_jpsi_prest,2,&mass[0]);
                    weight_decay = gen1->Generate();  //decay by phasespace, return weight as 1, replace it later with actual distribution
                    p4_je1_prest = gen1->GetDecay(0);
                    p4_je2_prest = gen1->GetDecay(1);

                    theta_e = p4_ep_prest->Theta()*simglobals::DEG;
		    eta_e = p4_ep_prest->PseudoRapidity();
                    phi_e = p4_ep_prest->Phi()*simglobals::DEG;

                    theta_p = p4_recoil_prest->Theta()*simglobals::DEG;
                    phi_p = p4_recoil_prest->Phi()*simglobals::DEG;

                    p_jpsi = p4_jpsi_prest->P();
                    theta_jpsi = p4_jpsi_prest->Theta()*simglobals::DEG;
                    phi_jpsi = p4_jpsi_prest->Phi()*simglobals::DEG;

                    p_je1 = p4_je1_prest->P();
                    theta_je1 = p4_je1_prest->Theta()*simglobals::DEG;
                    phi_je1 = p4_je1_prest->Phi()*simglobals::DEG;

                    p_je2 = p4_je2_prest->P();
                    theta_je2 = p4_je2_prest->Theta()*simglobals::DEG;
                    phi_je2 = p4_je2_prest->Phi()*simglobals::DEG;

                    *pt_prest = *p4_recoil_prest - *pTarget_prest;
                    t = -pt_prest->M2();

                    R = pow((a*mass_meson*mass_meson+Q2)/(a*mass_meson*mass_meson),n) -1;
                    //R defination and parameter a and n are from eq 18 of "R. Fiore et al. Exclusive Jpsi electroproduction in a dual model. Phys. Rev.,D80:116001, 2009"
                    theta_q = pq_prest->Theta()*simglobals::DEG;
                    q = pq_prest->P();
                    W = sqrt(pow(simglobals::mass_p + pq_prest->E(),2)-pow(pq_prest->P(),2));
                    Keq = (W*W-simglobals::mass_p*simglobals::mass_p)/2./simglobals::mass_p;
                    epsilon = 1./(1+2*q*q/Q2*pow(tan(theta_e/simglobals::DEG/2.),2));
                    Gamma = alpha/2./3.1415926/3.1415926*p_e/Ebeam*Keq/Q2/(1.-epsilon);
                    r= epsilon*R/(1.+epsilon*R);

                    // J = fabs((pq_prest->E()+simglobals::mass_p-q*p4_recoil_prest->E()/p_p*(cos(theta_q/simglobals::DEG)*cos(theta_p/simglobals::DEG)+sin(theta_q/simglobals::DEG)*sin(theta_p/simglobals::DEG)*sin((phi_p-phi_e+180.)/simglobals::DEG))*tan(theta_q/simglobals::DEG))
                    //     /(2.*simglobals::mass_p*q*p_p*(cos(theta_q/simglobals::DEG)*tan(theta_p/simglobals::DEG)-sin(theta_q/simglobals::DEG)*sin((phi_p-phi_e+180.)/simglobals::DEG))));
                    // cout << aa - bb*sqrt(p_p*p_p+simglobals::mass_p*simglobals::mass_p)-cc*p_p << endl;

                    Double_t dEpdcp=fabs((ps1_prest->Pz()-(ps1_prest->Px()*cos(phi_p/simglobals::DEG)+ps1_prest->Py()*sin(phi_p/simglobals::DEG))*cos(theta_p/simglobals::DEG)/sin(theta_p/simglobals::DEG))*p_p/(bb+cc*sqrt(p_p*p_p+simglobals::mass_p*simglobals::mass_p)/p_p));

                    J = 2 * simglobals::mass_p * dEpdcp;

                    //go to JPsi at rest frame
                    TVector3 beta = p4_jpsi_prest->Vect();
                    beta *= -1./p4_jpsi_prest->E();
		    *p4_je1_jpsirest = *p4_je1_prest;
                    p4_je1_jpsirest->Boost(beta);
                    //get recoil proton in the same frame
                    *p4_recoil_jpsirest = *p4_recoil_prest;
                    p4_recoil_jpsirest->Boost(beta);

                    *p4_q_jpsirest = *pq_prest;
                    p4_q_jpsirest->Boost(beta);
                    // calculate the theta angle between these two
                    TVector3 a1 = p4_je1_jpsirest->Vect();
                    TVector3 a2 = p4_recoil_jpsirest->Vect();
                    TVector3 a3 = p4_q_jpsirest->Vect();

                    a1.RotateUz(a2.Unit());
                    a3.RotateUz(a2.Unit());
                    theta_cm = a1.Theta()*simglobals::DEG;//p4_je1_jpsirest->Theta()*simglobals::DEG;
                    phi_cm = a1.Phi()*simglobals::DEG-a3.Phi()*simglobals::DEG;
                    if (phi_cm<0.) phi_cm+=360;
                    if (phi_cm>360) phi_cm-=360;
                    // phi_cm = p4_je1_jpsirest->Phi()*simglobals::DEG;
                    //theta_cm = 0.;

                    weight_decay = 3./2./4./3.1415926*(1-r + (3*r-1)*pow(cos(theta_cm/simglobals::DEG),2));
                    //            cout << weight_decay << endl;
                    //eq 92 of "K. Schilling and G. Wolf. How to analyze vector meson production in inelastic lepton scattering. Nucl. Phys., B61:381, 1973.",
                    // eq 32 and 31 from K. Schilling, P. Seyboth and G. Wolf, "ON THE ANALYSIS OF VECTOR-MESON PRODUCTION BY POLARIZED PHOTONS"  NucL Phys. B15 (1970) 397, B18 (1970) 332.
                    // eq 32 and 31 from K. Schilling, P. Seyboth and G. Wolf, "ON THE ANALYSIS OF VECTOR-MESON PRODUCTION BY POLARIZED PHOTONS"  NucL Phys. B15 (1970) 397, B18 (1970) 332.
                    //after removing all phi related term
                    //theta_cm and phi_cm as two degree of freedom, their phasespace size is 4pi, integral over them will give 1 with r cancls out

                    // calculate tmin
                    *ps2_prest = *pq_prest + *pTarget_prest;
                    tmin = -1*t0lim(-sqrt(Q2),simglobals::mass_p, mass_meson, simglobals::mass_p, ps2_prest->M2());
                    tmax = -1*t1lim(-sqrt(Q2),simglobals::mass_p, mass_meson, simglobals::mass_p, ps2_prest->M2());


                    //differential crossection in nb/(phasespace cell)		      
                    dxs     = J/2./3.1415926*Gamma*A*exp(b*(t));
                    dxs_2g  = J/2./3.1415926*Gamma*fun_2g(W,t,mass_meson);
                    dxs_23g = J/2./3.1415926*Gamma*fun_23g(W,t,mass_meson);


		    /* Calculate invariant mass from beam */
		    TLorentzVector p4_sum_beam(0.,0.,0.,0.);
		    p4_sum_beam += ( *pBeam_lab + *pTarget_lab );
		    minv_beam = p4_sum_beam.M();

		    /* Calculate invariant mass in 'proton at rest' frame */
		    TLorentzVector p4_sum_prest(0.,0.,0.,0.);
		    p4_sum_prest += ( *p4_ep_prest + *p4_recoil_prest + *p4_je1_prest + *p4_je2_prest );
		    minv_prest = p4_sum_prest.M();

		    /* Create vectors in laboratory frame */
		    *p4_ep_lab = *p4_ep_prest;
		    *p4_recoil_lab = *p4_recoil_prest;
		    *p4_jpsi_lab = *p4_jpsi_prest;
		    *p4_je1_lab = *p4_je1_prest;
		    *p4_je2_lab = *p4_je2_prest;

		    /* Boost final state particles to laboratory frame */
		    p4_ep_lab->Boost(-1*beta_lab_prest);
		    p4_recoil_lab->Boost(-1*beta_lab_prest);
		    p4_jpsi_lab->Boost(-1*beta_lab_prest);
		    p4_je1_lab->Boost(-1*beta_lab_prest);
		    p4_je2_lab->Boost(-1*beta_lab_prest);

		    /* Calculate invariant mass in laboratory frame */
		    TLorentzVector p4_sum(0.,0.,0.,0.);
		    p4_sum += ( *p4_ep_lab + *p4_recoil_lab + *p4_je1_lab + *p4_je2_lab );
		    minv = p4_sum.M();

		    /* re-calculate final state particle parameters */
		    p_e = p4_ep_lab->P();
		    pt_e = p4_ep_lab->Perp();
                    theta_e = p4_ep_lab->Theta()*simglobals::DEG;
                    eta_e = p4_ep_lab->PseudoRapidity();
                    phi_e = p4_ep_lab->Phi()*simglobals::DEG;

		    p_p = p4_recoil_lab->P();
		    pt_p = p4_recoil_lab->Perp();
                    theta_p = p4_recoil_lab->Theta()*simglobals::DEG;
                    eta_p = p4_recoil_lab->PseudoRapidity();
                    phi_p = p4_recoil_lab->Phi()*simglobals::DEG;

                    p_jpsi = p4_jpsi_lab->P();
                    pt_jpsi = p4_jpsi_lab->Perp();
                    theta_jpsi = p4_jpsi_lab->Theta()*simglobals::DEG;
                    eta_jpsi = p4_jpsi_lab->PseudoRapidity();
                    phi_jpsi = p4_jpsi_lab->Phi()*simglobals::DEG;

                    p_je1 = p4_je1_lab->P();
                    pt_je1 = p4_je1_lab->Perp();
                    theta_je1 = p4_je1_lab->Theta()*simglobals::DEG;
                    eta_je1 = p4_je1_lab->PseudoRapidity();
                    phi_je1 = p4_je1_lab->Phi()*simglobals::DEG;

                    p_je2 = p4_je2_lab->P();
                    pt_je2 = p4_je2_lab->Perp();
                    theta_je2 = p4_je2_lab->Theta()*simglobals::DEG;
                    eta_je2 = p4_je2_lab->PseudoRapidity();
                    phi_je2 = p4_je2_lab->Phi()*simglobals::DEG;
		    /* END boost final state particles back to laboratory frame */

		    /* Fill tree */
                    _tree_out->Fill();
                    if (neve1%(nevents/10)==0) cout << neve1 << endl;
                    neve1++;
                    if (neve1 > nevents) qflag = 0;

                  } //check if jpsi mass is positive
                  else counter[3]++;
                }//choose only positive recoil p mom
              } //end of loop through two solution
            }  //check if we can find positive solution for recoil p mom
            else counter[2]++;
          } //check if we can find solution for recoil p mom
          else counter[1]++;
        } //check if cross mass threshold
        else counter[0]++;
        neve++;
      }
    } // end of e beam

    if(Is_g){

      phasespace = (Ebeam-Gbeam_min)*(4*3.1415926)*(4*3.1415926);
      //        phasespace = (4*3.1415926);

      while(qflag){

        //sample photon energy
        Gbeam = fbr->GetRandom();
        //       Gbeam = gRandom->Uniform(7.5,Ebeam);
        //       Gbeam = 11;
        Gflux = fbr->Eval(Gbeam);
        //       cout << Gbeam << " " << Gflux << endl;

        pBeam_prest->SetPxPyPzE(0.,0.,Gbeam,Gbeam);
        pTarget_prest->SetPxPyPzE(0.,0.,0.,simglobals::mass_p);

        *ps_prest = *pBeam_prest + *pTarget_prest;

        *pq_prest = *pBeam_prest;

        Q2 = -pq_prest->M2();
        *ps1_prest = *ps_prest;

        //judge whether the event pass the kinematics
        if (ps1_prest->M2() > pow(mass_meson+simglobals::mass_p,2)){
          //sample proton solid angle;
          //    theta_p = acos(gRandom->Uniform(0.85,cos(8./simglobals::DEG)));
          theta_p = acos(gRandom->Uniform(-1,1)); //random selection in solid angle need to go with cos(theta)
          phi_p = gRandom->Uniform(0.,2.*3.1415926);

          //solve recoil proton mom
          TVector3 p3_p(sin(theta_p)*cos(phi_p),sin(theta_p)*sin(phi_p),cos(theta_p));
          Double_t aa = (ps1_prest->M2()+simglobals::mass_p*simglobals::mass_p-mass_meson*mass_meson)/2.;
          Double_t bb = ps1_prest->E();
          TVector3 p3_s(ps1_prest->Px(),ps1_prest->Py(),ps1_prest->Pz());
          Double_t cc = -p3_p.Dot(p3_s);
          Double_t sol[2];

          if (pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*simglobals::mass_p*simglobals::mass_p-aa*aa)>=0.){
            sol[0] = (-2*aa*cc+sqrt(pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*simglobals::mass_p*simglobals::mass_p-aa*aa)))/2./(bb*bb-cc*cc);
            sol[1] = (-2*aa*cc-sqrt(pow(2.*aa*cc,2)-4.*(bb*bb-cc*cc)*(bb*bb*simglobals::mass_p*simglobals::mass_p-aa*aa)))/2./(bb*bb-cc*cc);

            if (sol[0]>=0||sol[1]>=0){

              if (sol[0]>=0&&sol[1]<0){
                weight = 1;
              }else if (sol[1]>=0&&sol[0]<0){
                weight = 1;
              }else if (sol[0]>=0&&sol[1]>=0){
                weight = 0.5;
              }

              //                    cout << sol[0] << " " << sol[1] << endl;

              Double_t theta_ps_prest = theta_p;
              Double_t phi_ps_prest = phi_p;

              for (Int_t j=0;j!=2;j++){
                p_p = sol[j];
                if (p_p>0){
                  p4_recoil_prest->SetPxPyPzE(p_p*sin(theta_ps_prest)*cos(phi_ps_prest),p_p*sin(theta_ps_prest)*sin(phi_ps_prest),p_p*cos(theta_ps_prest),sqrt(p_p*p_p+simglobals::mass_p*simglobals::mass_p));
                  *p4_jpsi_prest = *ps1_prest-*p4_recoil_prest;

                  if (p4_jpsi_prest->M()>0.){
                    gen1->SetDecay(*p4_jpsi_prest,2,&mass[0]);
                    weight_decay = gen1->Generate();   //decay by phasespace, return weight as 1, replace it later with actual distribution
                    p4_je1_prest = gen1->GetDecay(0);
                    p4_je2_prest = gen1->GetDecay(1);

                    theta_e = p4_ep_prest->Theta()*simglobals::DEG;
                    phi_e = p4_ep_prest->Phi()*simglobals::DEG;

                    theta_p = p4_recoil_prest->Theta()*simglobals::DEG;
                    phi_p = p4_recoil_prest->Phi()*simglobals::DEG;

                    p_jpsi = p4_jpsi_prest->P();
                    theta_jpsi = p4_jpsi_prest->Theta()*simglobals::DEG;
                    phi_jpsi = p4_jpsi_prest->Phi()*simglobals::DEG;

                    p_je1 = p4_je1_prest->P();
                    theta_je1 = p4_je1_prest->Theta()*simglobals::DEG;
                    phi_je1 = p4_je1_prest->Phi()*simglobals::DEG;

                    p_je2 = p4_je2_prest->P();
                    theta_je2 = p4_je2_prest->Theta()*simglobals::DEG;
                    phi_je2 = p4_je2_prest->Phi()*simglobals::DEG;

                    *pt_prest = *p4_recoil_prest - *pTarget_prest;
                    t = -pt_prest->M2();

                    R = pow((a*mass_meson*mass_meson+Q2)/(a*mass_meson*mass_meson),n) -1;
                    theta_q = pq_prest->Theta()*simglobals::DEG;
                    q = pq_prest->P();
                    W = sqrt(pow(simglobals::mass_p + pq_prest->E(),2)-pow(pq_prest->P(),2));
                    Keq = (W*W-simglobals::mass_p*simglobals::mass_p)/2./simglobals::mass_p;
                    //            epsilon = 1./(1+2*q*q/Q2*pow(tan(theta_e/simglobals::DEG/2.),2));
                    //            Gamma = alpha/2./3.1415926/3.1415926*p_e/Ebeam*Keq/Q2/(1.-epsilon);
                    //            r= epsilon*R/(1.+epsilon*R);

                    // J = fabs((pq_prest->E()+simglobals::mass_p-q*p4_recoil->E()/p_p*(cos(theta_q/simglobals::DEG)*cos(theta_p/simglobals::DEG)+sin(theta_q/simglobals::DEG)*sin(theta_p/simglobals::DEG)*sin((phi_p-phi_e+180.)/simglobals::DEG))*tan(theta_q/simglobals::DEG))
                    //     /(2.*simglobals::mass_p*q*p_p*(cos(theta_q/simglobals::DEG)*tan(theta_p/simglobals::DEG)-sin(theta_q/simglobals::DEG)*sin((phi_p-phi_e+180.)/simglobals::DEG))));
                    // cout << aa - bb*sqrt(p_p*p_p+simglobals::mass_p*simglobals::mass_p)-cc*p_p << endl;

                    Double_t dEpdcp=fabs((ps1_prest->Pz()-(ps1_prest->Px()*cos(phi_p/simglobals::DEG)+ps1_prest->Py()*sin(phi_p/simglobals::DEG))*cos(theta_p/simglobals::DEG)/sin(theta_p/simglobals::DEG))*p_p/(bb+cc*sqrt(p_p*p_p+simglobals::mass_p*simglobals::mass_p)/p_p));
                    //cout << J << " " << 1./simglobals::mass_p/2./dEpdcp<< " " << dEpdcp << " " << p_p << " " << bb+cc*sqrt(p_p*p_p+simglobals::mass_p*simglobals::mass_p)/p_p << endl;
                    J = 2 * simglobals::mass_p * dEpdcp;

                    //go to JPsi at rest frame
                    TVector3 beta = p4_jpsi_prest->Vect();
                    beta *= -1./p4_jpsi_prest->E();
		    *p4_je1_jpsirest = *p4_je1_prest;
                    p4_je1_jpsirest->Boost(beta);
                    //get recoil proton in the same frame
                    *p4_recoil_jpsirest = *p4_recoil_prest;
                    p4_recoil_jpsirest->Boost(beta);

                    *p4_q_jpsirest = *pq_prest;
                    p4_q_jpsirest->Boost(beta);
                    // calculate the theta angle between these two
                    TVector3 a1 = p4_je1_jpsirest->Vect();
                    TVector3 a2 = p4_recoil_jpsirest->Vect();
                    TVector3 a3 = p4_q_jpsirest->Vect();

                    a1.RotateUz(a2.Unit());
                    a3.RotateUz(a2.Unit());
                    theta_cm = a1.Theta()*simglobals::DEG;//p4_je1_jpsirest->Theta()*simglobals::DEG;
                    phi_cm = a1.Phi()*simglobals::DEG-a3.Phi()*simglobals::DEG;
                    if (phi_cm<0.) phi_cm+=360;
                    if (phi_cm>360) phi_cm-=360;
                    // phi_cm = p4_je1_jpsirest->Phi()*simglobals::DEG;
                    //theta_cm = 0.;

                    r=0; //real photon has no transverse component
                    weight_decay = 3./2./4./3.1415916*(1-r + (3*r-1)*pow(cos(theta_cm/simglobals::DEG),2));
                    //            cout << weight_decay << endl;
                    //eq 92 of "K. Schilling and G. Wolf. How to analyze vector meson production in inelastic lepton scattering. Nucl. Phys., B61:381, 1973.",
                    // eq 32 and 31 from K. Schilling, P. Seyboth and G. Wolf, "ON THE ANALYSIS OF VECTOR-MESON PRODUCTION BY POLARIZED PHOTONS"  NucL Phys. B15 (1970) 397, B18 (1970) 332.
                    // eq 32 and 31 from K. Schilling, P. Seyboth and G. Wolf, "ON THE ANALYSIS OF VECTOR-MESON PRODUCTION BY POLARIZED PHOTONS"  NucL Phys. B15 (1970) 397, B18 (1970) 332.
                    //after removing all phi related term
                    //theta_cm and phi_cm as two degree of freedom, their phasespace size is 4pi, integral over them will give 1 with r cancls out

                    // calculate tmin
                    *ps2_prest = *pq_prest + *pTarget_prest;
                    tmin = -1*t0lim(-sqrt(Q2),simglobals::mass_p, mass_meson, simglobals::mass_p, ps2_prest->M2());
                    tmax = -1*t1lim(-sqrt(Q2),simglobals::mass_p, mass_meson, simglobals::mass_p, ps2_prest->M2());


                    //differential crossection in nb/(phasespace cell)
                    dxs     = A*exp(b*(t));
                    dxs_2g  = fun_2g(W,t,mass_meson);
                    dxs_23g = fun_23g(W,t,mass_meson);

                    _tree_out->Fill();
                    if (neve1%(nevents/10)==0) cout << neve1 << endl;
                    neve1++;
                    if (neve1 > nevents)  qflag = 0;

                  } //check if jpsi mass is positive
                  else counter[3]++;
                }//choose only positive recoil p mom
              } //end of loop through two solution
            }  //check if we can find positive solution for recoil p mom
            else counter[2]++;
          } //check if we can find solution for recoil p mom
          else counter[1]++;
        } //check if cross mass threshold
        else counter[0]++;
        neve++;
      }
    } // end of g beam

    cout << nevents << " events obtained after " <<  neve << " trials" << endl;

    cout << "counter " << counter[0] << " " << counter[1] << " " << counter[2] << " " << counter[3] << endl;


    end();
    return 0;
}



/////////////


int Simulator::init()
{

  create_output_file();

  cout << "init() done" << endl;
  return 0;
}

int Simulator::create_output_file()
{
  /* Outut file */
  _file_out = new TFile(output_root_file,"RECREATE");
  _tree_out = new TTree("T","T");
  _tree_out->SetDirectory(_file_out);

  cout << _file_out << endl;
    
  /* Set branches */
  _tree_out->Branch("Etarget",&Etarget_lab,"data/D");
  _tree_out->Branch("Ebeam",&Ebeam_lab,"data/D");

  _tree_out->Branch("Gbeam",&Gbeam,"data/D");
  _tree_out->Branch("Gflux",&Gflux,"data/D");

  _tree_out->Branch("weight_decay",&weight_decay,"data/D");
  _tree_out->Branch("weight",&weight,"data/D");

  _tree_out->Branch("Q2",&Q2,"data/D");
  _tree_out->Branch("t",&t,"data/D");

  _tree_out->Branch("m_inv",&minv,"data/D");
  _tree_out->Branch("m_inv_prest",&minv_prest,"data/D");
  _tree_out->Branch("m_inv_beam",&minv_beam,"data/D");

  _tree_out->Branch("p_e",&p_e,"data/D");
  _tree_out->Branch("pt_e",&pt_e,"data/D");
  _tree_out->Branch("theta_e",&theta_e,"data/D");
  _tree_out->Branch("phi_e",&phi_e,"data/D");
  _tree_out->Branch("eta_e",&eta_e,"data/D");

  _tree_out->Branch("p_p",&p_p,"data/D");
  _tree_out->Branch("pt_p",&pt_p,"data/D");
  _tree_out->Branch("theta_p",&theta_p,"data/D");
  _tree_out->Branch("eta_p",&eta_p,"data/D");
  _tree_out->Branch("phi_p",&phi_p,"data/D");

  _tree_out->Branch("p_jpsi",&p_jpsi,"data/D");
  _tree_out->Branch("pt_jpsi",&pt_jpsi,"data/D");
  _tree_out->Branch("theta_jpsi",&theta_jpsi,"data/D");
  _tree_out->Branch("eta_jpsi",&eta_jpsi,"data/D");
  _tree_out->Branch("phi_jpsi",&phi_jpsi,"data/D");

  _tree_out->Branch("p_je1",&p_je1,"data/D");
  _tree_out->Branch("pt_je1",&pt_je1,"data/D");
  _tree_out->Branch("theta_je1",&theta_je1,"data/D");
  _tree_out->Branch("eta_je1",&eta_je1,"data/D");
  _tree_out->Branch("phi_je1",&phi_je1,"data/D");

  _tree_out->Branch("p_je2",&p_je2,"data/D");
  _tree_out->Branch("pt_je2",&pt_je2,"data/D");
  _tree_out->Branch("theta_je2",&theta_je2,"data/D");
  _tree_out->Branch("eta_je2",&eta_je2,"data/D");
  _tree_out->Branch("phi_je2",&phi_je2,"data/D");

  _tree_out->Branch("neve",&neve,"data/I");

  _tree_out->Branch("Gamma",&Gamma,"data/D");
  _tree_out->Branch("epsilon",&epsilon,"data/D");
  _tree_out->Branch("Keq",&Keq,"data/D");
  _tree_out->Branch("W",&W,"data/D");
  _tree_out->Branch("q",&q,"data/D");
  _tree_out->Branch("theta_q",&theta_q,"data/D");
  _tree_out->Branch("J",&J,"data/D");
  _tree_out->Branch("R",&R,"data/D");
  _tree_out->Branch("theta_cm",&theta_cm,"data/D");
  _tree_out->Branch("phi_cm",&phi_cm,"data/D");
  _tree_out->Branch("r",&r,"data/D");

  _tree_out->Branch("dxs",&dxs,"data/D");
  _tree_out->Branch("dxs_2g",&dxs_2g,"data/D");
  _tree_out->Branch("dxs_23g",&dxs_23g,"data/D");

  _tree_out->Branch("tmin",&tmin,"data/D");
  _tree_out->Branch("tmax",&tmax,"data/D");

  _tree_out->Branch("phasespace",&phasespace,"data/D");

  cout << "create_output_file() done" << endl;
  return 0;
}

void Simulator::reset_event_variables()
{

  Gbeam = 0;
  Gflux = 0;

  Q2 = 0;
  t = 0;

  minv = 0;
  minv_prest = 0;
  minv_beam = 0;

  p_e = 0;
  pt_e = 0;
  theta_e = 0;
  phi_e = 0;
  eta_e = 0;
  p_p = 0;
  pt_p = 0;
  theta_p = 0;
  phi_p = 0;
  eta_p = 0;
  p_jpsi = 0;
  pt_jpsi = 0;
  theta_jpsi = 0;
  phi_jpsi = 0;
  eta_jpsi = 0;

  p_je1 = 0;
  pt_je1 = 0;
  theta_je1 = 0;
  phi_je1 = 0;
  eta_je1 = 0;
  p_je2 = 0;
  pt_je2 = 0;
  theta_je2 = 0;
  phi_je2 = 0;
  eta_je2 = 0;

  weight_decay = 0;
  weight = 0;
    
  Gamma = 0;
  epsilon = 0;
  Keq = 0;
  W = 0;
  q = 0;
  theta_q = 0;
  J = 0;
  R = 0;
  theta_cm = 0;
  phi_cm = 0;
  r = 0;

  dxs = 0;
  dxs_2g = 0;
  dxs_23g = 0;

  tmin = 0;
  tmax = 0;

  phasespace = 0;

  return;
}

int Simulator::process_event()
{
  cout << "process_event() done" << endl;
  return 0;
}

int Simulator::end()
{
  _file_out->Write();
  _file_out->Close();

  cout << "end() done" << endl;
  return 0;
}

//////////////

Double_t Simulator::fun_2g(Double_t x, Double_t t, Double_t M){
  Double_t N2g = 7.5671e3;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;
  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff;
  return result;
}

Double_t Simulator::fun_23g(Double_t x,Double_t t, Double_t M){
  Double_t N2g = 6.499e3;
  Double_t N3g = 2.894e3;
  Double_t v = 1./16/3.1415926;
  Double_t R = 1;

  Double_t xp = (2*0.938*M+M*M)/(x*x-0.938*0.938);
  Double_t ff = exp(-1.13 * t);

  Double_t result = N2g*v/R/R/M/M*pow(1-xp,2)*ff
    +N3g*v/R/R/R/R/M/M/M/M*pow(1-xp,0)*ff;
  return result;
}


Float_t Simulator::t0lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s)
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

Float_t Simulator::t1lim(Float_t m1, Float_t m2,Float_t m3, Float_t m4,Float_t s)
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



int main (Int_t argc, char *argv[])
{
  Simulator sim;

  if (argc==1)
    {
      cout << "./simu_p -n[nevents as integer like 1000000] -t[e for electroproduction,g for photoproduction] -b[Ebeam in GeV] -o[output_root_file]" << endl;
    }
  else{

    string type;

    for(Int_t i = 1; i != argc; i++){
      switch(argv[i][1]){
      case 'n':
	sim.set_number_events( int(atof(&argv[i][2])) );
	break;
      case 't':
	if ( sim.set_process_type( &argv[i][2] ) )
	  return 1;
	break;
      case 'b':
	sim.set_lepton_energy( atof(&argv[i][2]) );
	break;
      case 'p':
	sim.set_hadron_energy( atof(&argv[i][2]) );
	break;
      case 'm':
	if ( sim.set_meson_type( &argv[i][2] ) )
	  return 1;
	break;
      case 'o':
	sim.set_output_file( TString(&argv[i][2]) );
	break;
      default:
	cout << "Warning!!!! Unknown option: " << &argv[i][1] << endl;
	return 1;
	break;
      }
    }
  }

  sim.run();
  return 0;
}



