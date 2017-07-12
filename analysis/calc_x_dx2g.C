/** Calculate x (-> Brodsky's paper for 2-gluon cross section)
    Output crosschecked with paper, it matches
    x = 0.82 at E_gamma_lab = 10 GeV
    x = 0.69 at E_gamma_lab = 12 GeV */

float
calc_x_dx2g( )
{

  float x_calc = 0;

  float mp = 0.938272; // GeV

  float mJPsi = 3.097; // GeV
  float mUpsilon = 9.460; // GeV

  float M = mJPsi;
  //float M = mUpsilon;

  vector<float> v_Egamma_lab;

  for ( int i = 0; i < 100; i++ )
    v_Egamma_lab.push_back(10.0+i*2.0);

  for ( int i = 0; i < v_Egamma_lab.size(); i++ )
    {

      float Egamma_lab = v_Egamma_lab.at(i);

      float W = sqrt((Egamma_lab+0.938)*(Egamma_lab+0.938)-Egamma_lab*Egamma_lab);   //convert from E_gamma to W
      float Ecm = W;
      float s = Ecm * Ecm;

      x_calc = ( 2 * mp * M + M*M ) / ( s - mp * mp );

      cout << "E_gamma_lab = " << Egamma_lab << " , W = " << W << " , x = " << x_calc << endl;

    }
}
