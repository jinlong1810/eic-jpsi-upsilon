void plot_brem(){
// TCanvas *c_G = new TCanvas("G","G",1800,700);
// TF1 *f1 = new TF1("f1","fun_2g(x)",8,22);
TF1 *fbr = new TF1("fbr","0.001/(1-x/[0])*(4./3.-4./3.*(1-x/[0])+(1-x/[0])*(1-x/[0]))",0,11);
fbr->SetParameter(0,11);
// f->SetMaximum(1e2);
// f->SetMimimum(1e-4);
fbr->Draw();
}

