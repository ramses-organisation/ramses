# /patch/phantom

## About

This is the Phantom of RAMSES (POR) QUMOND patch.

For details, see Lueghausen et al. (2014), [arxiv:1405.5963](http://arxiv.org/abs/1405.5963)

Contact:
Fabian Lueghausen <lueghausen@gmail.com>



## Compile


### Makefile

In the RAMSES Makefile, make sure to set

    NDIM=3
    PATCH=../patch/phantom

The patch requires `NDIM=3`.


### Units

Because the POR patch makes use of the _poisson_ part, you must choose the units such that `G=1` in these units. This is a requirement of the RAMSES Poisson solver.

I.e., in the _units.f90_ file, make sure that

    scale_t = 1.0/sqrt(6.674D-8 * scale_d)

Always handle unit conversions with care!

Note. _In the code, the phantom dark matter (PDM) density is stored in user units (defined by `scale_d`)_.


### Code

The code changes are minimized, and new (non-standard-RAMSES) routines are added to the _poisson_commons.f90_ file.

In particular, see the subroutine `get_nu(x, nu)` in which the nu-function is defined:

	subroutine get_nu(x, nu)
	   use amr_commons
	   implicit none
	   real(dp),intent(in) :: x
	   real(dp),intent(out) :: nu

	   nu = -0.5 + 0.5*sqrt(x*x + 4.0*x)/x      ! Simple nu function
    end subroutine get_nu

It is **IMPORTANT** to ensure that nu(x) -> `0` if x >> 0.

If you change the nu-function, do not forget to re-compile the code!


## Namelist files

### General
To enable the QUMOND extension, set

	mond=.true.

in the `RUN_PARAMS` section. Otherwise, Newtonian forces are computed.

Set

	a0=...

in user units, **_or_**

	a0_ms2=...

in m/s^2 in the `POISSON_PARAMS` section. E.g., `a0_ms2=1.0e-10`.


### Boundary conditions

RAMSES uses the following definitions of boundary conditions:

    bound_type=0   Periodic boundary
    bound_type=1   Reflexive boundary
    bound_type=2   Gas outflow
    bound_type=3   Gas inflow

The POR patch uses periodic boundary conditions in the case of `bound_type=0`,
and Dirichlet boundary conditions (as defined in the referenced paper) if `bound_type>0`

**By default, if no boundary conditions are defined, periodic boundary conditions are applied.
If isolated objects are modelled, use** `bound_type>0`**.**

E.g., set

	&BOUNDARY_PARAMS
	nboundary=6
	ibound_min=-1, 1, 0, 0, 0, 0,
	ibound_max=-1, 1, 0, 0, 0, 0,
	jbound_min= 0, 0,-1, 1, 0, 0,
	jbound_max= 0, 0,-1, 1, 0, 0,
	kbound_min= 0, 0, 0, 0,-1, 1,
	kbound_max= 0, 0, 0, 0,-1, 1,
	bound_type= 1, 1, 1, 1, 1, 1,
	/

in the `BOUNDARY_PARAMS` section in order to apply reflexive boundary conditions.
