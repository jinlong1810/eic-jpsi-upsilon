
float
calc_Keq_to_W( float Keq = 5 )
{

  float mp = 0.938272; // GeV

  return ( sqrt( 2*mp*Keq + mp*mp ) );
}
