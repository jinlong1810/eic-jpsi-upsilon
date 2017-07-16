#include "normalization.h"

int plot_dsigma_dt_eic()
{
  TFile *fin = new TFile("../data2/sim_te_b5GeV_p50GeV_mupsilon.root","OPEN");
  TTree *T = (TTree*)fin->Get("T");

  //  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1);

  /* get normalization */
  Double_t overall = get_norm_eic_overall(T);

  /* define acceptance */
  char accep_normal[200];
  sprintf(accep_normal,"%s","((W<18)&&(Q2<2))");

  TH1F* count_dsigma_dt = new TH1F("count_dsigma_dt","",30,0,3);
  count_dsigma_dt->SetMarkerStyle(20);
  count_dsigma_dt->SetMarkerColor(kRed+1);
  count_dsigma_dt->GetYaxis()->SetTitle("# events");//#frac{d#sigma}{dt} (nb/GeV^{2})");
  count_dsigma_dt->GetXaxis()->SetTitle("|t-t_{min}| (GeV^{2})");

  const unsigned nbins_Keq = 1;
  float a_Keqmin[nbins_Keq] = {0};
  float a_Keqmax[nbins_Keq] = {10000};

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
      c1->Print(Form("plots/eic_dsigma_dt_KeqBin_%.2f_%.2f.png", Keqmin, Keqmax) );
    }

  return 0;
}
