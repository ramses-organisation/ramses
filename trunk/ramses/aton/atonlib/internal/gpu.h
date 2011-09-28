// GPU definitions.

// Unfortunately the grid size is fixed at compile time.
// You can set the dimensions here.
// For example, if you want a 256^3 simulation on 16 nodes, then set NCELLX=128, NCELLY=128, NCELLZ=64.

#define NCELLX 32
#define NCELLY 32
#define NCELLZ 32

#define NBOUND 1  // This must match nbnd in coupling.f90
#define NBOUND2 (2*NBOUND)


void gpu_rad_transport(double c_light, double dx, double dt);
void gpu_rad_add_sources(double dx, double dt, int nsource);
void gpu_rad_cooling(double c_light, double dx, double dt,
                     double aexp, double hubblet, double fudgecool);


// GPU global arrays
// They are defined in gpu_memory.cc.
// See aton_cpp.h or aton.f90 for descriptions.
// FIXME: These need more descriptive names.

extern double *cuegy, *cuegy_new;
extern double *cuflx, *cuflx_new;

extern double *cudedd;
extern double *cusrc0;
extern int  *cusrc0pos;

extern double *cutemperature;
extern double *cuxion;
extern double *cudensity;

extern double *cu_photon_source;

// Boundary values of cuegy and cuflx for efficient transfer at radiation
// substeps.
// length: 6 * 4 * max(ncellx,ncelly,ncellz)^2
// contents: [surface][{e,fx,fy,fz}][i][j]
extern double *cu_boundary_values;
