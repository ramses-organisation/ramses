Initial conditions
--------------------

In this chapter, we cover how initial conditions are prepared and read.
There are several ways to provide the starting state of a RAMSES
simulation. For basic geometries, the user can use the parameters in the
``init_params`` namelist block. For more advanced setups that can be
described analytically, there are the ``condinit`` subroutines. Finally,
RAMSES also supports input from files with specified formats.

.. contents::


1. Analytical initial conditions on the grid with ``condinit``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RAMSES users may already have experience with implementing initial
conditions using the ``condinit`` routine. This routine is called by
``init_flow_fine`` during the setup phase of the simulation. At the time
of writing, there are two versions available in the public version: one
for hydro and one for mhd. Remark that recently, the system has been
reworked to support various default setups (rather than replying on the
patch system).

As input, the ``condinit`` routine receives the cell center positions
(in code units) of the ``nn`` cells in the vector sweep, as well as the
cell size of the current level. From this information, the primitive
variables are then calculated. Several prescriptions are available by
default and can be selected by setting ``condinit_kind`` in the
namelist. One can easily add their own analytical initial conditions as
a new subroutine following the existing examples. Finally, the primitive
variables are converted to the conservative ones. The conservative
variables are then returned to ``init_flow_fine`` through the array
``u``, where they are written to ``uold``.

.. container:: info

   **Exercise:** Implement a new ``condinit_type`` that adds a
   sinusoidal perturbation on a uniform density background in 1D:
   :math:`\rho(x) = \rho_0 [1 + A \cos(\frac{2\pi x}{\lambda})]` The
   pressure is set to the same value as the density. :::spoiler
   **Solution**

   .. code:: fortran=

      ...
        case('jeans_instability_cos')
           call jeans_instability_cos_condinit(x, q, dx, nn)
      ...
      !================================================================
      subroutine jeans_instability_cos_condinit(x,q,dx,nn)
        use amr_parameters
        use hydro_parameters
        use constants, only:pi

        implicit none
        integer ::nn                            ! Number of cells
        real(dp)::dx                            ! Cell size
        real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
        real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
        !================================================================
        ! sinusoidal perturbation of amplitude A
        !================================================================
        integer::i
        real(dp),parameter::A=1d-4, lambda=0.5

        ! Call built-in initial condition generator to init the fields
        call region_condinit(x,q,dx,nn)

        do i=1,nn
          ! density
          q(i,1) = 1+A*cos(2*pi*x(i,1)/lambda)
          ! pressure
          q(i,3)=q(i,1)
        end do

      end subroutine jeans_instability_cos_condinit

2. Input file formats
~~~~~~~~~~~~~~~~~~~~~~

Another way to provide the initial conditions is through files, by
setting the namelist parameters ``initfile`` and ``filetype``. For the
variables on the grid, supported formats are ascii and grafic (see
*init_flow_fine.f90*), while for particles ascii, grafic, and gadget are
available (see *init_part.f90*, and `Section
1.3 <https://codimd.math.cnrs.fr/z55fgvBcTjiGxrWt7VBUVw?both#13-Initialising-particles>`__
in the chapter on particles).

If you want to add your own input file, good luck.
