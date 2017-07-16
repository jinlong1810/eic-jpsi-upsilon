#include "normalization.h"

int plot_particle_kinematics_eic( bool detailed=false )
{
  float p_max = 55;
  //  float p_max = 110;
  //  float p_max = 270;

  //float Wmax=7; // J/psi production
  float Wmax=18; // Upsilon production

  //TString fbase("sim_te_b5GeV_p50GeV_mjpsi");
  TString fbase("sim_te_b5GeV_p50GeV_mupsilon");
  //TString fbase("sim_te_b5GeV_p50GeV_mupsilon_Wmax20");
  TString fname("../data2/");
  fname+=fbase;
  fname+=".root";
  TFile *fin = new TFile(fname,"OPEN");
  TTree *T = (TTree*)fin->Get("T");

  TString p_base("plots/");
  p_base+=fbase;
  p_base+="_";


  gStyle->SetOptStat(0);


  TCanvas *cscratch = new TCanvas();

  /* Define acceptance */
  char acceptance_eicsphenix[200];
  sprintf(acceptance_eicsphenix,"( (W<%f) )",Wmax);
  cout << "Acceptance cut: " << acceptance_eicsphenix << endl;

  /* Define global weight */
  double weight_global = get_norm_eic_overall(T);

  /* Define event-by-event weight */
  /* select only dxs_2g < 1 to avoid bias by single high-crossection points */
  char weight_event[200];
  sprintf(weight_event,"%s","dxs_2g*weight*weight_decay");


  /* 1D and 2D histograms */
  TObjArray *a_th12f = new TObjArray();


  /* Momentum vs Pseudorapidity - Electron*/
  a_th12f->Add( new TH2D( "h_e_p_vs_eta",
                          "Scattered electron;#eta^{e'};p^{e'} (GeV);",
                          50,-10,10,50,0,p_max ) );
  T->Project( ((TH2D*)a_th12f->Last())->GetName(),
              "p_e:eta_e",
              Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
              );

  /* Momentum vs Pseudorapidity - J/Psi Electron*/
  a_th12f->Add( new TH2D( "h_je1_p_vs_eta",
                          "Scattered electron;#eta^{e1};p^{e1} (GeV);",
                          50,-10,10,50,0,p_max ) );
  T->Project( ((TH2D*)a_th12f->Last())->GetName(),
              "p_je1:eta_je1",
              Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
              );

  /* Momentum vs Pseudorapidity - J/Psi Positron*/
  a_th12f->Add( new TH2D( "h_je2_p_vs_eta",
                          "Scattered positron;#eta^{e2};p^{e2} (GeV);",
                          50,-10,10,50,0,p_max ) );
  T->Project( ((TH2D*)a_th12f->Last())->GetName(),
              "p_je2:eta_je2",
              Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
              );

  /* Momentum vs Pseudorapidity - Proton*/
  a_th12f->Add( new TH2D( "h_p_p_vs_eta",
                          "Scattered proton;#eta^{p'};p^{p'} (GeV);",
                          50,-10,10,50,0,p_max ) );
  T->Project( ((TH2D*)a_th12f->Last())->GetName(),
              "p_p:eta_p",
              Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
              );

  /* Momentum vs Pseudorapidity - J/Psi*/
  a_th12f->Add( new TH2D( "h_jpsi_p_vs_eta",
                          "JPsi;#eta^{VM};p^{VM} (GeV);",
                          50,-10,10,50,0,p_max ) );
  T->Project( ((TH2D*)a_th12f->Last())->GetName(),
              "p_jpsi:eta_jpsi",
              Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
              );


  /* SUPPORTING (DETAILED) PLOTS */

  /* 1D - Invariant mass */
  if (detailed)
    {
      a_th12f->Add( new TH1F( "h_minv",
                              "Invariant mass difference (final state) - (intial state);#Delta m_{inv};# events",
                              100,-1e-10,1e-10) );
      T->Project( ((TH1F*)a_th12f->Last())->GetName(),
                  "m_inv - m_inv_beam",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* 1D - Azimuthal angle */
  if ( detailed )
    {
      a_th12f->Add( new TH1F( "h_dphi",
                              "Azimuthal angle difference (fs J/#Psi) - (fs electron + fs proton);#Delta #phi;# events",
                              100,-200,200) );
      T->Project( ((TH1F*)a_th12f->Last())->GetName(),
                  "phi_jpsi-asin((pt_e*sin(phi_e*TMath::Pi()/180)+pt_p*sin(phi_p*TMath::Pi()/180))/sqrt(pow(pt_e*sin(phi_e*TMath::Pi()/180)+pt_p*sin(phi_p*TMath::Pi()/180),2)+pow(pt_e*cos(phi_e*TMath::Pi()/180)+pt_p*cos(phi_p*TMath::Pi()/180),2)))/TMath::Pi()*180",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* 2D - Azimuthal angle */
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_phi_jpsi_vs_epcombi",
                              "Azimuthal angle;#phi_{e'+p'};#phi_{J/#Psi};",
                              50,-200,200,50,-200,200) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "phi_jpsi:asin((pt_e*sin(phi_e*TMath::Pi()/180)+pt_p*sin(phi_p*TMath::Pi()/180))/sqrt(pow(pt_e*sin(phi_e*TMath::Pi()/180)+pt_p*sin(phi_p*TMath::Pi()/180),2)+pow(pt_e*cos(phi_e*TMath::Pi()/180)+pt_p*cos(phi_p*TMath::Pi()/180),2)))/TMath::Pi()*180",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* 2D - Azimuthal angle */
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_phi_e_vs_p",
                              "Azimuthal angle;#phi_{p'};#phi_{e'};",
                              50,-200,200,50,-200,200) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "phi_e:phi_p",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* 2D - Transverse momentum */
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_pt_jpsi_vs_epcombi",
                              "Transverse momentum;p_{T}^{e'+p'};p_{T}^{J/#Psi};",
                              50,0,20,50,0,20) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "pt_jpsi:sqrt( pow( pt_e*sin(phi_e*TMath::Pi()/180) + pt_p*sin(phi_p*TMath::Pi()/180), 2) + pow( pt_e*cos(phi_e*TMath::Pi()/180) + pt_p*cos(phi_p*TMath::Pi()/180) , 2))",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }


  /* Transverse momentum vs Pseudorapidity - Electron*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_e_pt_vs_eta",
                              "Scattered electron;#eta^{e'};p_{T}^{e'};",
                              50,-10,10,50,0,55 ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "pt_e:eta_e",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* Longitudinal momentum vs Pseudorapidity - Electron*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_e_pl_vs_eta",
                              "Scattered electron;#eta^{e'};p_{z}^{e'};",
                              50,-10,10,50,0,p_max ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "sqrt(p_e*p_e-pt_e*pt_e):eta_e",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* Transverse momentum / momentum vs Pseudorapidity - Electron*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_e_pt_over_p_vs_eta",
                              "Scattered electron;#eta^{e'};p_{T}^{e'}/p^{e'};",
                              50,-10,10,50,0,1 ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "pt_e/p_e:eta_e",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }


  /* Transverse momentum vs Pseudorapidity - J/Psi*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_jpsi_pt_vs_eta",
                              "JPsi;#eta^{J/#Psi};p_{T}^{J/#Psi};",
                              50,-10,10,50,0,55 ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "pt_jpsi:eta_jpsi",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* Longitudinal momentum vs Pseudorapidity - J/Psi*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_jpsi_pl_vs_eta",
                              "JPsi;#eta^{J/#Psi};p_{z}^{J/#Psi};",
                              50,-10,10,50,0,p_max ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "sqrt(p_jpsi*p_jpsi-pt_jpsi*pt_jpsi):eta_jpsi",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* Transverse momentum / momentum vs Pseudorapidity - J/Psi*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_jpsi_pt_over_p_vs_eta",
                              "JPsi;#eta^{J/#Psi};p_{T}^{J/#Psi}/p^{J/#Psi};",
                              50,-10,10,50,0,1 ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "pt_jpsi/p_jpsi:eta_jpsi",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* Pseudorapidity Electron vs Proton*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_e_eta_vs_p_eta",
                              "Electron eta vs proton eta;#eta^{p'};eta^{e'};",
                              50,-10,10,50,-10,10 ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "eta_e:eta_p",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* Pseudorapidity J/Psi vs Proton*/
  if ( detailed )
    {
      a_th12f->Add( new TH2D( "h_jpsi_eta_vs_p_eta",
                              "JPsi eta vs proton eta;#eta^{p'};eta^{J/#Psi};",
                              50,-10,10,50,-10,10 ) );
      T->Project( ((TH2D*)a_th12f->Last())->GetName(),
                  "eta_jpsi:eta_p",
                  Form("1*%s*%f*%s",acceptance_eicsphenix, weight_global, weight_event)
                  );
    }

  /* END DETAILED HISTOGRAMS */



  /* 1D - Plot and write all 1D histograms */
  for ( unsigned i = 0; i < a_th12f->GetEntries(); i++ )
    {

      /* scale entries */
      //      ((TH1F*)a_th12f->At(i))->Scale(1. / ((TH1F*)a_th12f->At(i))->Integral() );

      /* plot histogram on canvas */
      TCanvas *c_1d = new TCanvas();
      ((TH1F*)a_th12f->At(i))->SetLineColor(kBlue);
      ((TH1F*)a_th12f->At(i))->SetMarkerColor(kBlue);
      ((TH1F*)a_th12f->At(i))->SetMarkerSize(0.5);
      a_th12f->At(i)->Draw("colz");
      c_1d->SetLogz(1);
      c_1d->SetRightMargin(0.2);

      /* print canvas to file */
      TString p_name(p_base);
      p_name+=a_th12f->At(i)->GetName();
      p_name+=".eps";
      c_1d->Print(p_name);
    }

  cscratch->cd();

  return 0;
}
