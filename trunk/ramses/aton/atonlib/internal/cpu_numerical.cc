// CPU version of ATON.
// Ported from the GPU version.
// Warning: Needs cleaning.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "gpu.h"
#include "constants.h"
#include "aton_cpp.h"

#define ONE 0.99999f
#define ZERO 0.00001f
#define NCELL4 ((NCELLZ+NBOUND2)*(NCELLY+NBOUND2)*(NCELLX+NBOUND2))
#define TINY 1e-26
#define FLAGSOURCE 5.
#define KBOLTZ 1.3806e-23
#define EPSCOOL 0.0001

#define CUERR() printf("\n %s \n",cudaGetErrorString(cudaGetLastError()))

#ifdef TEST7_RAYTRACING_TRICKS

// Maximum eigenvalue in the y and z directions. This would normally be 1 for
// the GLF scheme but we can reduce it to 0.1 in the case of Test 7 because the
// photon flux is mostly directed in the x direction. The effect is to reduce
// diffusion in the y and z directions so that the clump casts a sharp shadow.
#define LAMBDA_YZ 0.1

#define OTSA

#else
#define LAMBDA_YZ 1.0
#endif

// This was introduced to prevent NaNs in Test 5 with 10^5 K blackbody spectrum.
// It could possibly be decreased futher.
#define EDDINGTON_EPSILON 0.0 //1.0e-5


// Bounds on physical quantities
#define MIN_TEMP 1.0e-20
#define MIN_EGY 0
#define MIN_X 0
#define MAX_X 0.99999999999999f
#define MIN_EINT 0


static int CellIndex(int i, int j, int k) {
  return ((i + NBOUND) +
	  (j + NBOUND)*(NCELLX + NBOUND2) +
	  (k + NBOUND)*(NCELLX + NBOUND2)*(NCELLY + NBOUND2));
}

static int CellIndex4(int i, int j, int k, int component) {
  return CellIndex(i, j, k) + component * NCELL4;
}


// cuegy_new: photon number density [1/m^3]
// cu_photon_source: photon number density source [1/m^3/s]
// dt: time step [seconds]
static void AddSourceField(double *cuegy_new,
			   double *cu_photon_source,
			   double dt) {
  double max_photon_source = 0.0;
  // Loop through all cells, except the boundary cells.
  for (int i = 0; i < NCELLX; i++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int k = 0; k < NCELLZ; k++) {
	int cell = CellIndex(i, j, k);
	cuegy_new[cell] += cu_photon_source[cell] * dt;
	max_photon_source = std::max(max_photon_source, cu_photon_source[cell]);
      }
    }
  }

  printf("AddSourceField: max_photon_source=%e, dt=%e\n",
	 max_photon_source, dt);
}







static void ComputeELFSingleCell(double *cuegy, double *cuflx,
				 double *cuegy_new,
				 double c, double dx, double dt, int i, int j, int k)
{
  double C = c*dt/dx; // Courant number (~1)
  
  double um1,up1,fm1,fp1,u0;

  // REMINDER LF flux : (fl+fr-ur+ul)*0.5f;
  //  f_icjcks_p =cuCompute_FaceLF(f[2+idx*3],f[2+idxp*3],c*e[idx],c*e[idxp]);

  double res = 0.0;

  // Divergence along Z

  um1 = cuegy[CellIndex(i, j, k-1)];
  u0  = cuegy[CellIndex(i, j, k  )];
  up1 = cuegy[CellIndex(i, j, k+1)];
  
  fm1 = cuflx[CellIndex4(i, j, k-1, 2)];
  fp1 = cuflx[CellIndex4(i, j, k+1, 2)];

  res = u0 - 0.5*(fp1-fm1)/dx*dt - 0.5*LAMBDA_YZ*C*(2*u0-um1-up1);
  if (isnan(res) || isinf(res) || res <= 0.0) {
    printf("ComputeELF bad 1:\n");
    printf("  res = %e, u0 = %e, fp1 = %e, fm1 = %e, um1 = %e, up1 = %e\n",
	   res, u0, fp1, fm1, um1, up1);
  }

  // Divergence along Y

  um1 = cuegy[CellIndex(i, j-1, k)];
  u0  = cuegy[CellIndex(i, j,   k)];
  up1 = cuegy[CellIndex(i, j+1, k)];

  fm1 = cuflx[CellIndex4(i, j-1, k, 1)];
  fp1 = cuflx[CellIndex4(i, j+1, k, 1)];
  
  res = res - 0.5*(fp1-fm1)/dx*dt - 0.5*LAMBDA_YZ*C*(2*u0-um1-up1);
  if (isnan(res) || isinf(res) || res <= 0.0) {
    printf("ComputeELF bad 2:\n");
    printf("  res = %e, u0 = %e, fp1 = %e, fm1 = %e, um1 = %e, up1 = %e\n",
	   res, u0, fp1, fm1, um1, up1);
  }
  
  //Divergence along X
  
  um1 = cuegy[CellIndex(i-1, j, k)];
  u0  = cuegy[CellIndex(i,   j, k)];
  up1 = cuegy[CellIndex(i+1, j, k)];

  fm1 = cuflx[CellIndex4(i-1, j, k, 0)];
  fp1 = cuflx[CellIndex4(i+1, j, k, 0)];

  res = res - 0.5*(fp1-fm1)/dx*dt - 0.5*C*(2*u0-um1-up1);
  if (isnan(res) || isinf(res) || res <= 0.0) {
    printf("ComputeELF bad 3:\n");
    printf("  res = %e, u0 = %e, fp1 = %e, fm1 = %e, um1 = %e, up1 = %e\n",
	   res, u0, fp1, fm1, um1, up1);
  }

  cuegy_new[CellIndex(i, j, k)] = res;
}

static void ComputeELF(double *cuegy, double *cuflx, double *cuegy_new,
		       double c, double dx, double dt) {
  for (int i = 0; i < NCELLX; i++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int k = 0; k < NCELLZ; k++) {
	ComputeELFSingleCell(cuegy, cuflx, cuegy_new,
			     c, dx, dt, i, j, k);
      }
    }
  }
}

// Returns the P tensor component.
static double Eddington(double fx, double fy, double fz, double ee, double c,int i,int j)
{
  double ff=0.;
  double arg,chi,res=0.;
  double n[3];
#ifdef ISOTROP
  
  if(i==j) res=1./3.;

#else
  n[0]=0.;n[1]=0.;n[2]=0.;

  if(ee > EDDINGTON_EPSILON)
    {
      ff=sqrtf(fx*fx+fy*fy+fz*fz);
      if(ff > EDDINGTON_EPSILON)
	{
	  n[0]=fx/ff;
	  n[1]=fy/ff;
	  n[2]=fz/ff;
	}
      ff=ff/(c*ee);
    }
  
  arg=fmaxf(4.-3.*ff*ff,0.);
  chi=(3.+4.*ff*ff)/(5.+2.*sqrtf(arg));

  if(i==j) res=(1.-chi)/2.;
  arg=(3.*chi-1.)/2.;
  res+=arg*n[i]*n[j];
#endif

  // These factors are applied in three steps to avoid double overflow issues.
  res *= ee;
  //res *= c;
  //res *= c;

  if (isnan(res) || isinf(res)) {
    printf("Eddington bad:\n");
    printf(" fx = %e, fy = %e, fz = %e\n", fx, fy, fz);
    printf(" ee = %e\n", ee);
    printf(" c = %e, i = %d, j = %d\n", c, i, j);
    printf(" res = %e\n", res);
    printf(" arg = %e, chi = %e, ff = %e\n", arg, chi, ff);
    printf(" n=[%e, %e, %e]\n", n[0], n[1], n[2]);
    printf(" ee*c*c = %e * %e * %e = %e\n", ee, c, c, ee*c*c);
  }

  return res;
}

static void ComputeF_TOTAL_LF_SingleCell(const double *cuflx, double *cuflx_new, double c, double dx, double dt, const double *cuegy, int i, int j, int k)
{
  double C = c*dt/dx; // Courant number (~1)

  double fm1,fp1;

  // REMINDER LF flux : (fl+fr-ur+ul)*0.5f;
  //  f_icjcks_p =cuCompute_FaceLF(f[2+idx*3],f[2+idxp*3],c*e[idx],c*e[idxp]);

  double resfx, resfy, resfz;

  double u[3], fp[3],fm[3], ep, em;

  //================================================ Z DIRECTION =============================================
  
  u[0]=cuflx[CellIndex4(i, j, k, 0)];
  u[1]=cuflx[CellIndex4(i, j, k, 1)];
  u[2]=cuflx[CellIndex4(i, j, k, 2)];

  ep=cuegy[CellIndex(i, j, k+1)];
  em=cuegy[CellIndex(i, j, k-1)];

  fm[0]=cuflx[CellIndex4(i, j, k-1, 0)];
  fm[1]=cuflx[CellIndex4(i, j, k-1, 1)];
  fm[2]=cuflx[CellIndex4(i, j, k-1, 2)];

  fp[0]=cuflx[CellIndex4(i, j, k+1, 0)];
  fp[1]=cuflx[CellIndex4(i, j, k+1, 1)];
  fp[2]=cuflx[CellIndex4(i, j, k+1, 2)];


  // FX Divergence along Z
  
  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,0,2);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,0,2);
  resfx = u[0] - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[0]-fm[0]-fp[0]);
  if (isnan(resfx) || isinf(resfx)) {
    printf("ComputeF bad resfx 1:\n"
	   "resfx=%e, fp1=%e, fm1=%e, u[0]=%e, fm[0]=%e, fp[0]=%e, "
	   "ep = %e, em = %e\n",
	   resfx, fp1, fm1, u[0], fm[0], fp[0], ep, em);
  }
  
  // FY Divergence along Z

  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,1,2);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,1,2);
  resfy = u[1] - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[1]-fm[1]-fp[1]);

  // FZ Divergence along Z

  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,2,2);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,2,2);
  resfz = u[2] - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[2]-fm[2]-fp[2]);


  //================================================ Y DIRECTION =============================================
  
  u[0]=cuflx[CellIndex4(i, j, k, 0)];
  u[1]=cuflx[CellIndex4(i, j, k, 1)];
  u[2]=cuflx[CellIndex4(i, j, k, 2)];

  ep=cuegy[CellIndex(i, j+1, k)];
  em=cuegy[CellIndex(i, j-1, k)];

  fm[0]=cuflx[CellIndex4(i, j-1, k, 0)];
  fm[1]=cuflx[CellIndex4(i, j-1, k, 1)];
  fm[2]=cuflx[CellIndex4(i, j-1, k, 2)];

  fp[0]=cuflx[CellIndex4(i, j+1, k, 0)];
  fp[1]=cuflx[CellIndex4(i, j+1, k, 1)];
  fp[2]=cuflx[CellIndex4(i, j+1, k, 2)];

  // FX Divergence along Y
  
  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,0,1);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,0,1);
  resfx = resfx - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[0]-fm[0]-fp[0]);
  if (isnan(resfx) || isinf(resfx)) {
    printf("ComputeF bad resfx 2:\n"
	   "resfx=%e, fp1=%e, fm1=%e, u[0]=%e, fm[0]=%e, fp[0]=%e, "
	   "ep = %e, em = %e\n",
	   resfx, fp1, fm1, u[0], fm[0], fp[0], ep, em);
  }

  // FY Divergence along Y

  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,1,1);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,1,1);
  resfy = resfy - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[1]-fm[1]-fp[1]);

  // FZ Divergence along Y

  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,2,1);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,2,1);
  resfz = resfz - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[2]-fm[2]-fp[2]);


  //================================================ X DIRECTION =============================================

  u[0]=cuflx[CellIndex4(i, j, k, 0)];
  u[1]=cuflx[CellIndex4(i, j, k, 1)];
  u[2]=cuflx[CellIndex4(i, j, k, 2)];

  ep=cuegy[CellIndex(i+1, j, k)];
  em=cuegy[CellIndex(i-1, j, k)];

  fm[0]=cuflx[CellIndex4(i-1, j, k, 0)];
  fm[1]=cuflx[CellIndex4(i-1, j, k, 1)];
  fm[2]=cuflx[CellIndex4(i-1, j, k, 2)];

  fp[0]=cuflx[CellIndex4(i+1, j, k, 0)];
  fp[1]=cuflx[CellIndex4(i+1, j, k, 1)];
  fp[2]=cuflx[CellIndex4(i+1, j, k, 2)];

  // FX Divergence along X

  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,0,0);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,0,0);
  resfx = resfx - 0.5*C*(fp1-fm1) - 0.5*C*(2*u[0]-fm[0]-fp[0]);
  if (isnan(resfx) || isinf(resfx)) {
    printf("ComputeF bad resfx 3:\n"
	   "resfx=%e, fp1=%e, fm1=%e, u[0]=%e, fm[0]=%e, fp[0]=%e, "
	   "ep = %e, em = %e\n",
	   resfx, fp1, fm1, u[0], fm[0], fp[0], ep, em);
  }

  // FY Divergence along X

  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,1,0);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,1,0);
  resfy = resfy - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[1]-fm[1]-fp[1]);

  // FZ Divergence along X

  fp1 = c * Eddington(fp[0],fp[1],fp[2],ep,c,2,0);
  fm1 = c * Eddington(fm[0],fm[1],fm[2],em,c,2,0);
  resfz = resfz - 0.5*C*(fp1-fm1) - 0.5*LAMBDA_YZ*C*(2*u[2]-fm[2]-fp[2]);


  cuflx_new[CellIndex4(i, j, k, 0)] = resfx;
  cuflx_new[CellIndex4(i, j, k, 1)] = resfy;
  cuflx_new[CellIndex4(i, j, k, 2)] = resfz;
}

static void ComputeF_TOTAL_LF(double *cuflx, double *cuflx_new,
		       double c, double dx, double dt, double *cuegy) {
  for (int i = 0; i < NCELLX; i++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int k = 0; k < NCELLZ; k++) {
	ComputeF_TOTAL_LF_SingleCell(cuflx, cuflx_new,
				     c, dx, dt, cuegy, i, j, k);
      }
    }
  }
}

static double cucompute_alpha_b(double temp) {
  // CASE B recombination rate m**3 s*-1
  // temperature should be given in Kelvin
  
  // Protection from divide-by-zero errors
  temp = fmax(temp, MIN_TEMP);

  double lambda = 2e0*157807e0 / temp;
  double alpha_b = 2.753e-14*powf(lambda,1.5)/powf(1e0+powf(lambda/2.740,0.407),2.242); //cm3/s
  alpha_b=alpha_b*1e-6; //m3/s
  return alpha_b;
}

static double cucompute_alpha_a(double temp) {
  // CASE A recombination rate m**3 s*-1
  // temperature should be given in Kelvin
 
  // Protection from divide-by-zero errors
  temp = fmax(temp, MIN_TEMP);
 
  double lambda=2e0*157807e0/temp;
  double alpha_a=1.269e-13*powf(lambda,1.503)/powf(1e0+powf(lambda/0.522,0.470),1.923); //cm3/s
  alpha_a=alpha_a*1e-6; //m3/s
  return alpha_a;
}

static double cucompute_beta(double temp) {
  // Protection from divide-by-zero errors
  temp = fmax(temp, MIN_TEMP);

  // Collizional ionization rate m**3 s*-1
  // temperature in Kelvin
  double beta,T5;
  T5=temp/1e5;
  beta=5.85e-11*sqrtf(temp)/(1+sqrtf(T5))*expf(-(157809e0/temp)); //cm3/s
  beta=beta*1e-6; // !m3/s
  return beta;
}

static void CompCooling(double temp, double x, double nH, double *lambda, double *tcool, double aexp) {
  // Protection from divide-by-zero errors
  temp = fmax(temp, MIN_TEMP);
 

  double nh2 = nH*1e-6;// ! m-3 ==> cm-3


  // Collisional Ionization Cooling
  double c1 = expf(-157809.1e0/temp)*1.27e-21*sqrtf(temp)/(1e0+sqrtf(temp/1e5))*x*(1-x)*nh2*nh2;
  
  // Case A Recombination Cooling
  double c2 = 1.778e-29*temp*powf(2e0*157807e0/temp,1.965e0)/powf(1e0+powf(2e0*157807e0/temp/0.541e0,0.502e0),2.697e0)*x*x*nh2*nh2;
    
  // Case B Recombination Cooling
  //c6=3.435e-30*temp*powf(2e0*157807e0/temp,1.970e0)/powf(1e0+(powf(2e0*157807e0/temp/2.250e0,0.376e0)),3.720e0)*x*x*nh2*nh2;
  double c6 = 0.0;
  
  // Collisional excitation cooling
  double c3 = expf(-118348e0/temp)*7.5e-19/(1+sqrtf(temp/1e5))*x*(1-x)*nh2*nh2;  
  
  // Bremmsstrahlung
  double c4 = 1.42e-27*1.5e0*sqrtf(temp)*x*x*nh2*nh2;
  
  // Compton Cooling and Heating
  double c5 = 1.017e-37*powf(2.727/aexp,4)*(temp-2.727/aexp)*nh2*x;

  
  // Overall Cooling
  *lambda = c1+c2+c3+c4+c5+c6; // erg*cm-3*s-1
  
  // Unit Conversion
  *lambda=(*lambda)*1e-7*1e6; // J*m-3*s-1


  // cooling times
  double unsurtc = fmaxf(c1,c2);
  unsurtc=fmaxf(unsurtc,c3);
  unsurtc=fmaxf(unsurtc,c4);
  unsurtc=fmaxf(unsurtc,fabs(c5));
  unsurtc=fmaxf(unsurtc,c6)*1e-7;// ==> J/cm3/s
  *tcool=1.5e0*nh2*(1+x)*1.3806e-23*temp/unsurtc; //Myr
}

void ComputeTempSingleCell(double *cuxion,
			   double *cudensity,
			   double *cutemperature,
			   double *cuegy_new,
			   double *cuflx_new_x,
			   double *cuflx_new_y,
			   double *cuflx_new_z,
			   double fudgecool,
			   double c,
			   double dt,
			   double aexp,
			   double hubblet,
			   int *num_iterations) {
  double hnu,hnu0;
  double Cool;
  double tcool1;
  
  hnu=AVG_EGY*1.6022e-19;//     ! Average Photon Energy (J)
  hnu0=13.6*1.6022e-19;
  double sige=AVG_CSE;
  double sign=AVG_CSN;

  //double alpha,beta,alphai=1.6279e-22*c,alphae=1.0970e-22*c;
  double alpha,alphab,beta;
  double alphai=sign * 2.998e8;
  double alphae=sige * 2.998e8;

  double egyloc;
  double floc[3];
  double x0;
  double nH;
  double tloc;

  x0 = *cuxion;
  nH = *cudensity;
  egyloc = *cuegy_new;
  tloc = *cutemperature;
  floc[0] = *cuflx_new_x;
  floc[1] = *cuflx_new_y;
  floc[2] = *cuflx_new_z;

  double currentcool_t = 0.f;
  int nitcool = 0;
  double dtlimit = dt;
  while (currentcool_t < dt) {
    nitcool++;
    if (nitcool > 10000) {
      printf("bad bailing out nitcool=%d\n", nitcool);
      break;
    }

    double eint = 1.5*nH*KBOLTZ*(1.+x0)*tloc;

    //== Getting a timestep
    // 1. Don't evolve past the full time step.
    double dtcool = dt - currentcool_t;
    // 2. Maybe subcycle based on the cooling time.
    CompCooling(tloc,x0,nH,&Cool,&tcool1,aexp);
    const double Heat = egyloc*nH*(1-x0) *
      (alphae*hnu-alphai*hnu0);
    const double eint_rate = Heat - Cool;
    if (fabsf(eint_rate) > 1.0e-30) {
      double tcool = fabsf(eint/eint_rate);
      dtcool = fmin(dtcool, tcool);
    }
    dtcool = fmin(dtcool, dtlimit);

    //== Cross sections
#ifndef OTSA
    alpha=cucompute_alpha_a(tloc);
#else
    alpha=cucompute_alpha_b(tloc);
#endif
    alphab=cucompute_alpha_b(tloc);
    beta=cucompute_beta(tloc);


    //== Update
    double q = (alpha-alphab)*x0*x0*nH*nH;
    double p = alphai*nH*(1-x0);
    double r = 1.0;

    //== Update photon density
    double new_egy = (q*dtcool + egyloc) / (1 + dtcool*(p + 3*hubblet));

    //== Update the flux
    double new_fx = r*floc[0] / (1.f + dtcool*(p + 2*hubblet));
    double new_fy = r*floc[1] / (1.f + dtcool*(p + 2*hubblet));
    double new_fz = r*floc[2] / (1.f + dtcool*(p + 2*hubblet));
    double F = sqrtf(new_fx*new_fx + new_fy*new_fy + new_fz*new_fz);

    //== Update ionized fraction
    double new_x = 1 - (alpha*x0*x0*nH*dtcool+1.f -x0)/(1.f+dtcool*(beta*x0*nH+alphai*new_egy+3*hubblet));

    //== Update internal energy
    // The heating and cooling rates are computed using the new density values.
    CompCooling(tloc,new_x,nH,&Cool,&tcool1,aexp);
    eint = (eint + dtcool*(new_egy*nH*(1.f-new_x)*(alphae*hnu-alphai*hnu0)-Cool)) / (1 + 3*hubblet*dtcool);

    //printf("   eint=%e dtcool=%e Cool=%e\n", eint, dtcool, Cool);

    if (new_egy < MIN_EGY) {
      dtlimit = fudgecool * dtcool;
      //printf("continue new_egy\n");
      continue;
    }
    new_egy = fmaxf(new_egy, MIN_EGY);
    if ((new_x < MIN_X || new_x > MAX_X)) {
      dtlimit = fudgecool * dtcool;
      //printf("continue new_x\n");
      continue;
    }
    new_x = fminf(fmaxf(new_x, MIN_X), MAX_X);
    if (eint < MIN_EINT) {
      dtlimit = fudgecool * dtcool;
      //printf("continue eint\n");
      continue;
    }
    eint = fmaxf(eint, MIN_EINT);

    /*
    if (F > c*new_egy && fudgecool > 1.0e-5) {
      fudgecool = fudgecool / 10.0;
      continue;
    } else if (F > c*new_egy) {
      double f = c*new_egy / F;
      new_fx *= f;
      new_fy *= f;
      new_fz *= f;
    }
    */

    //== Save the new values
    egyloc = new_egy;
    floc[0] = new_fx;
    floc[1] = new_fy;
    floc[2] = new_fz;
    x0 = new_x;
    tloc = eint / (1.5f*nH*KBOLTZ*(1 + x0));

    currentcool_t += dtcool;

    //printf("nitcool=%d t=%e egyloc=%e x0=%e tloc=%e\n",
    //nitcool, currentcool_t, egyloc, x0, tloc);
  }

  if (nitcool > 10000 || isnan(tloc) || isinf(tloc) || tloc <= 0.0 ||
      egyloc <= 0.0 || x0 <= 0.0) {
    double Fold = sqrtf((*cuflx_new_x)*(*cuflx_new_x) +
			(*cuflx_new_y)*(*cuflx_new_y) +
			(*cuflx_new_z)*(*cuflx_new_z));
    double Fnew = sqrtf(floc[0]*floc[0] +
			floc[1]*floc[1] +
			floc[2]*floc[2]);

    printf("ComputeTemp bad:\n");
    printf(" nitcool = %d\n", nitcool);
    printf(" nH: %e\n", *cudensity);
    printf(" T: %e -> %e\n", *cutemperature, tloc);
    printf(" x: %e -> %e\n", *cuxion, x0);
    printf(" N: %e -> %e\n", *cuegy_new, egyloc);
    printf(" F: [%e,%e,%e] -> [%e,%e,%e]\n",
	   *cuflx_new_x, *cuflx_new_y, *cuflx_new_z,
	   floc[0], floc[1], floc[2]);
    printf(" F: %e -> %e\n", Fold, Fnew);
  }

  *num_iterations = nitcool;

  *cutemperature = tloc;
  *cuxion = x0;
  *cuegy_new = egyloc;
  *cuflx_new_x = floc[0];
  *cuflx_new_y = floc[1];
  *cuflx_new_z = floc[2];
}

static void ComputeTemp(double *cuxion,
			double *cudensity,
			double *cutemperature,
			double *cuegy_new,
			double *cuflx_new,
			double fudgecool,
			double c,
			double dt,
			double aexp,
			double hubblet) {

  double sum_iterations = 0.0;
  int max_iterations = 0;

  // Loop through all cells, except the boundary cells.
  for (int i = 0; i < NCELLX; i++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int k = 0; k < NCELLZ; k++) {
	int cell = CellIndex(i, j, k);
	int num_iterations;
	ComputeTempSingleCell(&cuxion[cell],
			      &cudensity[cell],
			      &cutemperature[cell],
			      &cuegy_new[cell],
			      &cuflx_new[cell+0*NCELL4],
			      &cuflx_new[cell+1*NCELL4],
			      &cuflx_new[cell+2*NCELL4],
			      fudgecool,
			      c,
			      dt,
			      aexp,
			      hubblet,
			      &num_iterations);
	sum_iterations += num_iterations;
	max_iterations = std::max(max_iterations, num_iterations);
      }
    }
  }

  double avg_iterations = sum_iterations / NCELL4;
  printf("avg_iterations = %f\n", avg_iterations);
  printf("max_iterations = %d\n", max_iterations);

}

static void Validate(int label, double c,
		     double *cpu_e, double* cpu_d, double* cpu_t, double* cpu_x,
		     double *cpu_f) {
  printf("Validate %d\n", label);

  for (int i = 0; i < NCELLX; i++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int k = 0; k < NCELLZ; k++) {
	int b = CellIndex(i, j, k);

	bool bad = false;

	if (cpu_d[b] == 0.0) {
	  printf("(%d) zero cpu_d\n", label);
	  bad = true;
	}

	if (isnan(cpu_e[b])) {
	  printf("(%d) nan cpu_e\n", label);
	  bad = true;
	}
	if (isinf(cpu_e[b])) {
	  printf("(%d) inf cpu_e\n", label);
	  bad = true;
	}
	if (cpu_e[b] <= 0.0) {
	  printf("(%d) cpu_e <= 0: %e\n", label, cpu_e[b]);
	  bad = true;
	}

	if (isnan(cpu_t[b])) {
	  printf("(%d) nan cpu_t\n", label);
	  bad = true;
	}
	if (isinf(cpu_t[b])) {
	  printf("(%d) inf cpu_t\n", label);
	  bad = true;
	}
	if (cpu_t[b] <= 0.0) {
	  printf("(%d) cpu_t <= 0: %e\n", label, cpu_t[b]);
	  bad = true;
	}

	if (isnan(cpu_x[b])) {
	  printf("(%d) nan cpu_x\n", label);
	  bad = true;
	}
	if (isinf(cpu_x[b])) {
	  printf("(%d) inf cpu_x\n", label);
	  bad = true;
	}
	if (cpu_x[b] <= 0.0) {
	  printf("(%d) cpu_x <= 0: %e\n", label, cpu_x[b]);
	  bad = true;
	}

	if (isnan(cpu_f[b+0*NCELL4])) {
	  printf("(%d) nan cpu_f_x\n", label);
	  bad = true;
	}
	if (isnan(cpu_f[b+1*NCELL4])) {
	  printf("(%d) nan cpu_f_y\n", label);
	  bad = true;
	}
	if (isnan(cpu_f[b+2*NCELL4])) {
	  printf("(%d) nan cpu_f_z\n", label);
	  bad = true;
	}
	if (isinf(cpu_f[b+0*NCELL4])) {
	  printf("(%d) inf cpu_f_x\n", label);
	  bad = true;
	}
	if (isinf(cpu_f[b+1*NCELL4])) {
	  printf("(%d) inf cpu_f_y\n", label);
	  bad = true;
	}
	if (isinf(cpu_f[b+2*NCELL4])) {
	  printf("(%d) inf cpu_f_z\n", label);
	  bad = true;
	}

	double F = sqrtf(cpu_f[b+0*NCELL4]*cpu_f[b+0*NCELL4] +
			 cpu_f[b+1*NCELL4]*cpu_f[b+1*NCELL4] +
			 cpu_f[b+2*NCELL4]*cpu_f[b+2*NCELL4]);
	if (F > c*cpu_e[b]) {
	  printf("(%d) F > cN. %e > %e * %e\n",
		 label, F, c, cpu_e[b]);
	  bad = true;
	}
	

	if (bad) {
	  printf("(%d) bad at %d,%d,%d, idx=%d\n", label, i, j, k, b);
	  return;
	}

      }
    }
  }
}



static void CheckCells(int label, double c,
		       const double *cpu_e, const double *cpu_f,
		       const double *cpu_new_e, const double *cpu_new_f) {

  for (int i = 0; i < NCELLX; i++) {
    for (int j = 0; j < NCELLY; j++) {
      for (int k = 0; k < NCELLZ; k++) {
	int b = CellIndex(i, j, k);

	double old_F = sqrtf(cpu_f[b+0*NCELL4]*cpu_f[b+0*NCELL4] +
			     cpu_f[b+1*NCELL4]*cpu_f[b+1*NCELL4] +
			     cpu_f[b+2*NCELL4]*cpu_f[b+2*NCELL4]);
	double new_F = sqrtf(cpu_new_f[b+0*NCELL4]*cpu_new_f[b+0*NCELL4] +
			     cpu_new_f[b+1*NCELL4]*cpu_new_f[b+1*NCELL4] +
			     cpu_new_f[b+2*NCELL4]*cpu_new_f[b+2*NCELL4]);

	bool old_bad = old_F > c*cpu_e[b];
	bool new_bad = new_F > c*cpu_new_e[b];

	if (new_bad && !old_bad) {
	  printf("%d: Cell check fail at %d,%d,%d\n", label, i, j, k);
	  printf("%d: begin\n", label);

	  for (int ii=i-1; ii<=i+1; ii++) {
	    for (int jj=j-1; jj<=j+1; jj++) {
	      for (int kk=k-1; kk<=k+1; kk++) {
		int index = CellIndex(ii, jj, kk);
		printf("%d:   %d %d %d  N: %e -> %e\n", label, ii, jj, kk,
		       cpu_e[index], cpu_new_e[index]);
		printf("%d:   %d %d %d  Fx: %e -> %e\n", label, ii, jj, kk,
		       cpu_f[index+0*NCELL4], cpu_new_f[index+0*NCELL4]);
		printf("%d:   %d %d %d  Fy: %e -> %e\n", label, ii, jj, kk,
		       cpu_f[index+1*NCELL4], cpu_new_f[index+1*NCELL4]);
		printf("%d:   %d %d %d  Fz: %e -> %e\n", label, ii, jj, kk,
		       cpu_f[index+2*NCELL4], cpu_new_f[index+2*NCELL4]);
		printf("\n");
	      }
	    }
	  }
	  printf("%d: end\n", label);
	}
      }
    }
  }
}


namespace aton {
  
void cpu_transport(State state, double c_light, double dx, double dt) {
  // FIXME: allocate these only once
  double* cpu_e_new = (double*) malloc(NCELL4*sizeof(double));
  double* cpu_f_new = (double*) malloc(NCELL4*sizeof(double)*3);

  ComputeELF(state.E, state.F, cpu_e_new, c_light, dx, dt);
  ComputeF_TOTAL_LF(state.F, cpu_f_new, c_light, dx, dt, state.E);
  
  memcpy(state.E, cpu_e_new, NCELL4*sizeof(double));
  memcpy(state.F, cpu_f_new, NCELL4*sizeof(double)*3);
}


void cpu_add_sources(State state, double dx, double dt) {
  AddSourceField(state.E, state.photon_source, dt);
}


void cpu_cooling(State state, double c_light, double dx, double dt,
                 double aexp, double hubblet, double fudgecool) {
  ComputeTemp(state.xHII, state.nH, state.T, state.E, state.F,
	      fudgecool, c_light, dt, aexp, hubblet);
}

}


extern "C" void aton_cpu_rad_(
    int *myid,
    double *c, double *dx, double *dt, int *nsource,
    double *fudgecool, double *aexp, double *hubblet,
    double *cpu_e, double* cpu_d, double* cpu_t, double* cpu_x,
    double *cpu_photon_source, double *cpu_f) {

  printf("myid=%d aton_cpu_rad_ begin\n", *myid);
  Validate(1, *c, cpu_e, cpu_d, cpu_t, cpu_x, cpu_f);

  // FIXME: allocate these only once
  double* cpu_e_new = (double*) malloc(NCELL4*sizeof(double));
  double* cpu_f_new = (double*) malloc(NCELL4*sizeof(double)*3);

  printf("myid=%d ComputeELF:\n", *myid);
  ComputeELF(cpu_e, cpu_f, cpu_e_new, *c, *dx, *dt);
  Validate(2, *c, cpu_e_new, cpu_d, cpu_t, cpu_x, cpu_f);

  printf("myid=%d ComputeF_TOTAL_LF:\n", *myid);
  ComputeF_TOTAL_LF(cpu_f, cpu_f_new, *c, *dx, *dt, cpu_e);
  Validate(3, *c, cpu_e_new, cpu_d, cpu_t, cpu_x, cpu_f_new);

  CheckCells(50000+*myid, *c, cpu_e, cpu_f, cpu_e_new, cpu_f_new);

  printf("myid=%d AddSourceField:\n", *myid);
  AddSourceField(cpu_e_new, cpu_photon_source, *dt);
  Validate(4, *c, cpu_e_new, cpu_d, cpu_t, cpu_x, cpu_f_new);

  printf("myid=%d ComputeTemp:\n", *myid);
  ComputeTemp(cpu_x, cpu_d, cpu_t, cpu_e_new, cpu_f_new,
	      *fudgecool, *c, *dt, *aexp, *hubblet);
  Validate(5, *c, cpu_e_new, cpu_d, cpu_t, cpu_x, cpu_f_new);

  printf("myid=%d memcpu:\n", *myid);
  memcpy(cpu_e, cpu_e_new, NCELL4*sizeof(double));
  memcpy(cpu_f, cpu_f_new, NCELL4*sizeof(double)*3);

  printf("myid=%d aton_cpu_rad_ done.\n", *myid);
}
