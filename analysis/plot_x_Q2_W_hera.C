int plot_x_Q2_W_hera()
{
  gStyle->SetOptStat(0);

  /* calc W */
  TF1* f_W1 = new TF1("f_W1", "([0]*[0] - 0.934*0.934) / (1/x - 1)",1e-6,0.999);
  f_W1->SetParameter(0,200);
  TF1 *f_W2 = (TF1*)f_W1->Clone("f_W2");
  f_W2->SetParameter(0,100);
  TF1 *f_W3 = (TF1*)f_W1->Clone("f_W3");
  f_W3->SetParameter(0,50);
  TF1 *f_W4 = (TF1*)f_W1->Clone("f_W4");
  f_W4->SetParameter(0,10);
  TF1 *f_W5 = (TF1*)f_W1->Clone("f_W5");
  f_W5->SetParameter(0,3);

  /* calc y */
  TF1 *f_y1 = new TF1("f_y1", "4*x*[0]*[1]*[2]", 1e-6, 0.999);
  f_y1->SetParameter( 0, 26.7);
  f_y1->SetParameter( 1, 820);
  f_y1->SetParameter( 2, 1);
  TF1 *f_y01 = (TF1*)f_y1->Clone("f_y01");
  f_y01->SetParameter(2 , 0.1);
  TF1 *f_y001 = (TF1*)f_y1->Clone("f_y001");
  f_y001->SetParameter(2 , 0.01);
  TF1 *f_y0001 = (TF1*)f_y1->Clone("f_y0001");
  f_y0001->SetParameter(2 , 0.001);
  TF1 *f_y00001 = (TF1*)f_y1->Clone("f_y00001");
  f_y00001->SetParameter(2 , 0.0001);

  /* frame */
  TH1F* hframe_hera = new TH1F("hframe_hera","Hera;x;Q^{2} (GeV^{2})",10,1e-6,1);
  hframe_hera->GetYaxis()->SetRangeUser(1e-1,1e5);

  /* plot canvas and lines */
  TCanvas *c1 = new TCanvas();
  hframe_hera->Draw();

  f_W1->Draw("same");
  f_W2->Draw("same");
  f_W3->Draw("same");
  f_W4->Draw("same");
  f_W5->Draw("same");

  f_y1->Draw("SAME");
  f_y1->SetLineColor(1);
  f_y01->Draw("SAME");
  f_y01->SetLineColor(1);
  f_y001->Draw("SAME");
  f_y001->SetLineColor(1);
  f_y0001->Draw("SAME");
  f_y0001->SetLineColor(1);
  f_y00001->Draw("SAME");
  f_y00001->SetLineColor(1);


  c1->SetLogx(1);
  c1->SetLogy(1);

  return 0;
}
