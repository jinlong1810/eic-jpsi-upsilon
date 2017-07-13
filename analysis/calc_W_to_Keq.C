
float
calc_W_to_Keq( float W = 5 )
{

  float mp = 0.938272; // GeV

  return ( ( W*W - mp*mp ) / ( 2*mp ) );
}
