#include "normalization.h"

int plot_dsigma_dt_solid()
{
  TFile *fin_11GeV = new TFile("../data/sim_te_b11GeV_p0GeV_mjpsi_v2_accep.root","OPEN");
  TTree *T = (TTree*)fin_11GeV->Get("T");

  //  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);

  /* get normalization */
  Double_t overall = get_norm_solid_overall(T);


  /* define acceptance */
  char accep_normal[200];
  sprintf(accep_normal,"%s","(accep_je1_1+accep_je1_2)*(accep_je2_1+accep_je2_2)*(accep_e_1)*(accep_p_1+accep_p_2)");

  TH1F* count_dsigma_dt = new TH1F("count_dsigma_dt","",30,0,3);
  count_dsigma_dt->SetMarkerStyle(20);
  count_dsigma_dt->SetMarkerColor(kRed+1);
  count_dsigma_dt->GetYaxis()->SetTitle("# events");//#frac{d#sigma}{dt} (nb/GeV^{2})");
  count_dsigma_dt->GetXaxis()->SetTitle("|t-t_{min}| (GeV^{2})");

  const unsigned nbins_Keq = 5;
  float a_Keqmin[nbins_Keq] = {9.09, 9.29, 9.49, 9.69, 0};
  float a_Keqmax[nbins_Keq] = {9.29, 9.49, 9.69, 9.89, 1000};

  for ( unsigned i = 0; i < nbins_Keq; i++ )
    {
      count_dsigma_dt->Reset();

      float Keqmin = a_Keqmin[i];
      float Keqmax = a_Keqmax[i];

      T->Project("count_dsigma_dt","abs(t-tmin)",Form("dxs_2g*weight*weight_decay*%s*%f*(Keq>=%f && Keq<%f)",accep_normal,overall,Keqmin,Keqmax));

      TCanvas *c1 = new TCanvas();
      c1->SetTitle( Form("Keq range: %.2f <= Keq < %.2f", Keqmin, Keqmax) );

      cout << "Entries: " << count_dsigma_dt->GetEntries() << endl;

      count_dsigma_dt->DrawClone("");
      c1->SetLogx(0);
      c1->SetLogy(1);
      c1->Print(Form("plots/solid_dsigma_dt_KeqBin_%.2f_%.2f.png", Keqmin, Keqmax) );
    }

  return 0;
}
