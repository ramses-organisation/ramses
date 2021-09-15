/* This is code designed to.... */
/* But what will the types be for each of the inputs? Furthermore, */
/* what will even be the inputs? */
/* Which kick is the "first" and which is the "second"? */
/* Declare our packages */
#include <iostream>
#include <cmath>
#include <tuple>
#include <utility>
using namespace std;

void field_interpolate() /* This is here just in case we need it */
{
  return;
}

tuple<float,float,float> BorisKick( /* "kick" is either 0 or 1, referring to the first
  or second kick. */
  int kick, float dt,
  float ctm, float ts, /* Charge to mass ratio, stopping time */
  float b1, float b2, float b3, /* Magnetic field */
  float u1, float u2, float u3, /* Fluid velocity */
  float v1, float v2, float v3 /* grain velocity */
)
{
  float v1h, v2h, v3h; /* "half" velocity to be output */
  float v1p, v2p, v3p;

  if (kick==2) /* The "second"" kick */
  {
    v1p = v1 + (2*ctm*dt*(-(b2*b2*ctm*dt*v1) + b2*(b1*ctm*dt*v2
          - 2*v3) + b3*(-(b3*ctm*dt*v1) + 2*v2
          + b1*ctm*dt*v3)))/(4 +
          (b1*b1 + b2*b2 + b3*b3)*ctm*ctm*dt*dt);
    v2p = v2 + (2*ctm*dt*(-(b3*b3*ctm*dt*v2) + b1*(b2*ctm*dt*v1
          - b1*ctm*dt*v2 + 2*v3) + b3*(-2*v1
          + b2*ctm*dt*v3)))/(4
          + (b1*b1 + b2*b2 + b3*b3)*ctm*ctm*dt*dt);
    v3p = v3 + (2*ctm*dt*(2*b2*v1 + b1*b3*ctm*dt*v1 - 2*b1*v2
          + b2*b3*ctm*dt*v2 - (b1*b1
          + b2*b2)*ctm*dt*v3))/(4 +
          (b1*b1 + b2*b2 + b3*b3)*ctm*ctm*dt*dt);
    /* The vp's advance us to the next step, but only with rotation. */
    /* They don't add in the energy changes in any way */
    v1h = (v1p - 0.5*dt*(ctm*(u2*b3-u3*b2)-u1/ts))/(1.+0.5*dt/ts);
    v2h = (v2p - 0.5*dt*(ctm*(u3*b1-u1*b3)-u2/ts))/(1.+0.5*dt/ts);
    v3h = (v3p - 0.5*dt*(ctm*(u1*b2-u2*b1)-u3/ts))/(1.+0.5*dt/ts);
  }
  else /* The "first" kick */
  {
    v1h = (v1 - 0.5*dt*(ctm*(u2*b3-u3*b2)-u1/ts))/(1.+0.5*dt/ts);
    v2h = (v2 - 0.5*dt*(ctm*(u3*b1-u1*b3)-u2/ts))/(1.+0.5*dt/ts);
    v3h = (v3 - 0.5*dt*(ctm*(u1*b2-u2*b1)-u3/ts))/(1.+0.5*dt/ts);

  }

 return make_tuple(v1h,v2h,v3h);

}
// Tomorrow, I should turn main() into a loop and test the decaying particle
// velocity. 
int main()
{
  float v1, v2, v3;
  float vgr1, vgr2, vgr3;
  float u1, u2, u3;
  float b1, b2, b3;
  float ctm, ts;
  float dt;
  int kick;

  kick = 1; dt=1.e-5;
  ctm=2.e4; ts=1.e-4;
  b1=0.;b2=0.;b3=10.;
  u1=0.;u2=0.;u3=0.;
  vgr1=100.;vgr2=0.;vgr3=0.;
  tie(v1,v2,v3) = BorisKick(1,dt,ctm,ts,b1,b2,b3,u1,u2,u3,vgr1,vgr2,vgr3);
  tie(v1,v2,v3) = BorisKick(2,dt,ctm,ts,b1,b2,b3,u1,u2,u3,v1,v2,v3);
  cout << '(' << v1 << ',' << v2 << ',' << v3 << ')' << endl;
  cout << sqrt(v1*v1+v2*v2+v3*v3) << endl;
  return 0;
}
