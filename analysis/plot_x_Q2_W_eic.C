int plot_x_Q2_W_eic()
{
  gStyle->SetOptStat(0);

  /* calc W */
  TF1* f_W1 = new TF1("f_W1", "([0]*[0] - 0.934*0.934) / (1/x - 1)",1e-6,0.999);
  f_W1->SetLineColor(kBlue);
  f_W1->SetParameter(0,4);
  TF1 *f_W2 = (TF1*)f_W1->Clone("f_W2");
  f_W2->SetParameter(0,10);
  TF1 *f_W3 = (TF1*)f_W1->Clone("f_W3");
  f_W3->SetParameter(0,30);
  TF1 *f_W4 = (TF1*)f_W1->Clone("f_W4");
  f_W4->SetParameter(0,130);

  /* calc y */
  TF1 *f_y1a = new TF1("f_y1a", "4*x*[0]*[1]*[2]", 1e-6, 0.999);
  f_y1a->SetLineColor(kGray+2);
  f_y1a->SetParameter( 0, 5);
  f_y1a->SetParameter( 1, 50);
  f_y1a->SetParameter( 2, 0.95);
  TF1 *f_y1b = (TF1*)f_y1a->Clone("f_y1b");
  f_y1b->SetParameter( 2, 0.01);

  TF1 *f_y2a = (TF1*)f_y1a->Clone("f_y2a");
  f_y2a->SetLineStyle(2);
  f_y2a->SetParameter( 0, 20);
  f_y2a->SetParameter( 1, 250);
  f_y2a->SetParameter( 2, 0.95);
  TF1 *f_y2b = (TF1*)f_y2a->Clone("f_y2a");
  f_y2b->SetParameter(2 , 0.01);

  /* frame */
  TH1F* hframe_eic = new TH1F("hframe_eic","EIC;x;Q^{2} (GeV^{2})",10,4e-5,1);
  hframe_eic->GetYaxis()->SetRangeUser(7e-1,2e3);

  /* plot canvas and lines */
  TCanvas *c1 = new TCanvas();
  hframe_eic->Draw();

  f_W1->Draw("same");
  f_W2->Draw("same");
  f_W3->Draw("same");
  f_W4->Draw("same");

  f_y1a->Draw("SAME");
  f_y1b->Draw("SAME");

  f_y2a->Draw("SAME");
  f_y2b->Draw("SAME");

  c1->SetLogx(1);
  c1->SetLogy(1);

  gPad->RedrawAxis();

  c1->Print("coverage_x_Q2_y_W_eic.png");

  return 0;
}
