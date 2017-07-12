///// WORK IN PROGRESS

int
plot_tmax_tmin_analytical()
{
  gStyle->SetOptStat(0);

  /* vector meson masses */
  //float m_jpsi = 3.096916; // GeV
  //float m_upsi = 9.46030; // GeV

  /* vector meson production thresholds */
  float t_jpsi = 4.035; // GeV
  float t_upsilon = 10.4; // GeV

  /* frame for J/psi at SoLID */
  TH1F* hframe_jpsi = new TH1F("hframe_jpsi","",1000,0,100);
  hframe_jpsi->SetLineColor(kWhite);
  hframe_jpsi->GetXaxis()->SetNdivisions(505);
  hframe_jpsi->GetYaxis()->SetNdivisions(504);
  hframe_jpsi->GetXaxis()->SetTitle("W (GeV)");
  hframe_jpsi->GetYaxis()->SetTitle("-t (GeV^{2})");
  hframe_jpsi->GetXaxis()->SetRangeUser(4,4.6);
  hframe_jpsi->GetYaxis()->SetRangeUser(0.001,10);

  /* frame for Upsilon at EIC */
  TH1F* hframe_upsilon = (TH1F*)hframe_jpsi->Clone("hframe_upsilon");
  hframe_upsilon->GetXaxis()->SetRangeUser(9,20);
  hframe_upsilon->GetYaxis()->SetRangeUser(0.001,30);

  /* Functions for J/PSi at SoLID */
  TF1 *f_t0 = new TF1("ft0","fun_tlimit(x, 0)",t_jpsi,4.55);
  TF1 *f_t1 = new TF1("ft1","fun_tlimit(x, 1)",t_jpsi,4.55);
  f_t1->SetLineColor(kRed);

  /* Functions for Upsilon at EIC */
  TF1 *f_t0_upsilon = new TF1("fupsit0","fun_tlimit_upsilon(x, 0)",t_upsilon,18.9);
  TF1 *f_t1_upsilon = new TF1("fupsit1","fun_tlimit_upsilon(x, 1)",t_upsilon,18.9);
  f_t1_upsilon->SetLineColor(kRed);

  /* SOLID W covereage */
  TLine *ref_solid_Wmin = new TLine(4.05,0,4.05,10);
  TLine *ref_solid_Wmax = new TLine(4.45,0,4.45,10);

  ref_solid_Wmin->SetLineColor(kBlue);
  ref_solid_Wmax->SetLineColor(kBlue);

  ref_solid_Wmin->SetLineWidth(2);
  ref_solid_Wmax->SetLineWidth(2);

  /* EIC W covereage */
  TLine *ref_eic_Wmin = new TLine(8.0,0,8.0,10);
  TLine *ref_eic_Wmax = new TLine(30.0,0,30.0,10);

  ref_eic_Wmin->SetLineColor(kBlue);
  ref_eic_Wmax->SetLineColor(kBlue);

  ref_eic_Wmin->SetLineWidth(2);
  ref_eic_Wmax->SetLineWidth(2);

  /* J/Psi vs W at SoLID*/
  TCanvas *c_t0_t1_jpsi_W = new TCanvas("c_t0_t1_jpsi_W","c_t0_t1_jpsi_W");
  hframe_jpsi->Draw();
  f_t0->Draw("same");
  f_t1->Draw("same");
  ref_solid_Wmin->Draw("same");
  ref_solid_Wmax->Draw("same");
  c_t0_t1_jpsi_W->RedrawAxis();

  /* Upsilon vs W at EIC */
  TCanvas *c_t0_t1_upsilon_W = new TCanvas("c_t0_t1_upsilon_W","c_t0_t1_upsilon_W");
  hframe_upsilon->Draw();
  f_t0_upsilon->Draw("same");
  f_t1_upsilon->Draw("same");
  //  ref_eic_Wmin->Draw("same");
  //  ref_eic_Wmax->Draw("same");
  c_t0_t1_upsilon_W->RedrawAxis();

  /* print some reference numbers */
  cout << "J/psi, W = 15 GeV, t0 = " << fun_tlimit( 15, 0 ) << ", t1 = " << fun_tlimit( 15, 1 ) << endl;
  cout << "Upsilon, W = 15 GeV, t0 = " << fun_tlimit_upsilon( 15, 0 ) << ", t1 = " << fun_tlimit_upsilon( 15, 1 ) << endl;

  return 0;
}


float
fun_tlimit( float x, int t01 )
{
  /* Scattering:
   * particle 1 = gamma*
   * particle 2 = proton
   * particle 3 = J/psi
   * particle 4 = proton'
   */

  /* result */
  float tlimit = 0;

  float W = x;
  float s = W*W;

  float m_jpsi = 3.096916; // GeV
  float m_vmeson = m_jpsi;
  float m_proton = 0.938272; // GeV

  //  float Q2 = 0;
  float m_gammax = 0;

  float m1 = m_gammax;
  float m2 = m_proton;
  float m3 = m_vmeson;
  float m4 = m_proton;

  float E1cm = ( s + m1*m1 - m2*m2 ) / ( 2*sqrt(s) );
  float p1cm = sqrt( E1cm * E1cm - m1*m1 );

  float E3cm = ( s + m3*m3 - m4*m4 ) / ( 2*sqrt(s) );
  float p3cm = sqrt( E3cm*E3cm - m3*m3 );

  float temp1 = ( m1*m1 - m2*m2 - m3*m3 + m4*m4 )/( 2 * sqrt( s ) );

  float temp2 = 0;
  if ( t01 == 0 )
    temp2 = p1cm - p3cm;
  else if ( t01 == 1 )
    temp2 = p1cm + p3cm;
  else
    temp2 = -999;

  tlimit = temp1*temp1 - temp2*temp2;

  return -1 * tlimit;

}


float
fun_tlimit_upsilon( float x, int t01 )
{
  /* Scattering:
   * particle 1 = gamma*
   * particle 2 = proton
   * particle 3 = Upsilon
   * particle 4 = proton'
   */

  /* result */
  float tlimit = 0;

  float W = x;
  float s = W*W;

  float m_upsilon = 9.46030; // GeV
  float m_vmeson = m_upsilon;
  float m_proton = 0.938272; // GeV

  //  float Q2 = 0;
  float m_gammax = 0;

  float m1 = m_gammax;
  float m2 = m_proton;
  float m3 = m_vmeson;
  float m4 = m_proton;

  float E1cm = ( s + m1*m1 - m2*m2 ) / ( 2*sqrt(s) );
  float p1cm = sqrt( E1cm * E1cm - m1*m1 );

  float E3cm = ( s + m3*m3 - m4*m4 ) / ( 2*sqrt(s) );
  float p3cm = sqrt( E3cm*E3cm - m3*m3 );

  float temp1 = ( m1*m1 - m2*m2 - m3*m3 + m4*m4 )/( 2 * sqrt( s ) );

  float temp2 = 0;
  if ( t01 == 0 )
    temp2 = p1cm - p3cm;
  else if ( t01 == 1 )
    temp2 = p1cm + p3cm;
  else
    temp2 = -999;

  tlimit = temp1*temp1 - temp2*temp2;

  return -1 * tlimit;

}
