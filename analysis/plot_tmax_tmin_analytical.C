///// WORK IN PROGRESS

int
plot_tmax_tmin_analytical()
{

  Float_t temp1,temp2,temp3,temp4;

  float Q2 = 10;
  //  float W = ;
  float s = 141;

  float m1 = Q2;
  float m2 = 1;
  float m3 = 4;
  float m4 = m2;

  temp1 = (m1*m1 - m3*m3 - m2*m2 + m4*m4)/(2*sqrt(s));

  temp2 = (s + m1*m1 - m2*m2)/(2*sqrt(s));    //E1cm
  temp2 = sqrt(temp2*temp2 - m1*m1);                //p1cm

  temp3 = (s + m3*m3 - m4*m4)/(2*sqrt(s));    //E3cm
  temp3 = sqrt(temp3*temp3 - m3*m3);                //p3cm

  temp4 = temp2 - temp3;                            //p1cm-p3cm

  cout << temp1*temp1 - temp4*temp4 << endl;

  return 0;
}
