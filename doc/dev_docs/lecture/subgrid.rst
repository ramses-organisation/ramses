Lecture: Subgrid modelling utilities
====================================

Introduction
------------

Many processes happen on unresolved scales in any simulation: star
formation, supernova explosions, accretion onto black-holes are typical
examples of that in cosmological or galaxy-scale simulations. The
effects of these unresolved processes on the surrounding, resolved,
medium need to be modeled explicitly, as they do not result from the set
of equations (Euler+…) solved by RAMSES. This is where **subgrid
models** come in.

From a practical viewpoint, many subgrid models deal with the formation
or growth of particles, and with the injection of matter, momentum,
energy, etc. from these particles into the surrounding medium. This
means **updating particle properties based on grid properties**, and
conversely **updating grid quantities based on particle properties**.

Here, we discuss basic algorithms that are commonly used to perform
these operations in RAMSES. We will base the explanation on the examples
of *star formation* and *SN feedback* subgrid models as a way to
illustrate the different operations we need to perform on the grid and
on particles.

A **star formation subgrid model** operates on leaf cells in typically
in 3 steps: 1. Decide whether star formation may occur in the current
cell. This depends on a set of conditions (defined by the modeller).
These conditions may be local (e.g., a density threshold), or non-local,
i.e., depend on the surrounding gas properties (e.g., a converging
flow). In the latter case we need to collect information beyond the
current leaf cell. 2. Decide how much gas mass will turn into stars.
This depends on the star formation model and may again be a function of
local properties (e.g., with a SF efficiency set by the leaf cell’s
free-fall time) or non-local properties (e.g., with the multi-free-fall
models). 3. Remove gas from the mesh and spawn new star particles.

In turn, **SN feedback subgrid models** need to: 1. Decide whether a
stellar particle releases matter and/or momentum and energy. This may
happen as a one-shot event (typically 5-10 Myr after the stellar
population formed), or may be modelled as a rate of events. 2. Figure
out in which cell(s) the star particle will inject stuff. This may be a
single cell (for metal release, thermal energy) or a number of cells
(e.g. for mechanical feedback). 3. Update the star particle and the
cell(s) according to the feedback model. Again this can be done in
several ways (including through the use of *debris* particles… ).

.. container:: info

   **Exercise –** These models apply to star formation and supernovae.
   Can you think of other examples where similar operations can apply?
   :::spoiler Solution - Injection of radiation from stars: this depends
   on the age/metallicity/mass of the star particle, and adds photons to
   the cells where star particles are located. Have a look at
   ``star_RT_feedback`` in ``rt/rt_spectra.f90``. - Most of AGN physics:
   \* BH growth through gas accretion requires reading properties in a
   region around the BH particles, removing some of the gas in that
   region, and adding it to the BH \* AGN feedback requires depositing
   mass, momentum, and energy in cells surrounding the BH particles -
   Other “sink” modules (e.g., sink particles to model stars in
   high-resolution simulations)

Generic operations
------------------

For the models described above, can you define a few **generic
operations** that we need to perform: - Find the host cell of a particle
- Find the 26 neighbours of a cell (:math:`3^3` cells minus the central
cell) - Identify all particles (e.g., stars) within a cell. - Collect
particles / cells within some distance of a point - Remove stuff from
cells/particles and add it to particles/cells (without breaking
conservation laws)

We will now look at places in the code where these operations are done.

Finding the host cell of a particle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A typical example of this operation is for supernova feedback: we need
to find the cell in which energy can be dumped.

Let’s look at the ``thermal_feedback`` routine in ``feedback.f90``. As
we have seen in the lecture on particles, the way to loop over particles
*already* requires iterating over grids. This is the loop that looks
like this:

.. code:: fortran

     do icpu=1,ncpu
     ! Loop over cpus
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           if(npart1>0)then
              ipart=headp(igrid)
              ! Loop over particles
              do jpart=1,npart1
                 ! Save next particle   <--- Very important !!!
                 next_part=nextp(ipart)
                 
                 !----
                 ! Do something with particle ipart
                 !----
                 
                 ipart=next_part  ! Go to next particle
              end do
           endif
           igrid=next(igrid)   ! Go to next grid
        end do
     end do

In the ``thermal_feedback`` routine, this is performed in two steps:
first, we count how many particles are stars:

.. code:: fortran

       if ( is_star(typep(ipart)) ) then
           npart2=npart2+1
       endif

Then, in a second step, we do something with the star particles by
repeating the loop in a *vectorized* way (see the vector sweep for the
AMR structure, for example). This calls an *inner* function, ``feedbk``,
where the actual feedback is done.

It is only in this ``feedbk`` routine that the actual **cell** is
identified, using a simpler version of the CIC scheme we have seen
previously. Indeed, for this, RAMSES uses the NGP scheme (for Nearest
Grid Point).

.. container:: info

   **Exercise –** Go through the ``feedbk`` routine and explain how NGP
   differs from CIC.

Finding the neighbours of a cell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the ``star_formation`` routine of ``pm/star_formation.f90``, some
star formation models require to identify the neighbouring cells of
star-forming sites in order to compute the velocity dispersion around a
given cell.

For this, we use the ``getnbor`` routine defined in
``pm/star_formation.f90``, and that heavily relies on routines developed
in ``amr/nbors_utils.f90``. The routine returns the (global) index,
called ``ind_nbor`` of all ``2*ndim`` cells around a cell indexed by
``ind_cell``.

This index can then be used to access the properties of the neighbours.
For example, the density in the six cells around the target one are
defined through

.. code:: fortran

       d1 = uold(ind_nbor(1,1),1)
       d2 = uold(ind_nbor(1,2),1)
       d3 = uold(ind_nbor(1,3),1)
       d4 = uold(ind_nbor(1,4),1)
       d5 = uold(ind_nbor(1,5),1)
       d6 = uold(ind_nbor(1,6),1)

Note that this ``ind_nbor`` hides sometimes complex situations. For
example: - the target cell and its neighbours can be at different levels
- in some cases, a neighbour can be non-existent

The ``getnbor`` routine itself follows the logic described in the
`lecture on the AMR
structure <https://codimd.math.cnrs.fr/s/S9DoFQ2c6#16-Finding-neighbors-WIP>`__.

.. container:: info

   **Exercise –** Go through the ``pm/move_tracer.f90`` file to identify
   how tracer particles are identified and moved to neighbouring cells.

Identifying particles within a cell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As we have seen several times now, the particle linked list is defined
using the *grid* structure, through the ``headp`` and ``nextp`` arrays:
``headp(igrid)`` points to the index of the first particle in a given
grid ``igrid``, and ``nextp(ipart)`` returns the next particle in the
linked list.

To identify the particles within a given *cell*, the easiest way is then
to simply identify the grid corresponding to a cell, and loop over the
particles in that grid. Functions like ``is_tracer(typep(ipart))``,
``is_star(typep(ipart))``, etc. can then be used to select specifically
all the particles of a given type.

Collect particles and cells in a sphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to identify particles or cells in a given region is to
use the particle/grid iteration mechanism, and add a criterion based on
the position.

For example, if we want to select all the particles with a radius
``rmax`` of a position ``xtarget``, we can just add

.. code:: fortran

   dist2 = (xp(ipart, 1)-xtarget(1))**2 + (xp(ipart, 2)-xtarget(2))**2 + (xp(ipart, 3)-xtarget(3))**2

   if (dist2 .le. rmax**2) then
       ...
   endif

in the particle loop.

Let’s now look at the *kinetic feedback* implementation for a practical
example of selecting *cells* in a given region. The idea of the model is
to add momentum to cells within a given radius ``rmax`` around
supernovae sites. This is implemented in ``pm/feedback.f90`` using three
routines: - ``kinetic_feedback``, the main driver, which loops over
“debris particles” (you can think of these as “old stars with wind
particles attached to them”) and fills ``xSN``, ``vSN``, ``mSN``, and
other “supernova” arrays, - ``average_SN``, which accounts for the grid
discretization effects in the way momentum is scaled in the SN region, -
``Sedov_blast``, which actually performs the feedback, and injects
momentum and energy to the grid.

In ``Sedov_blast``, the main loop is the following:

.. code:: fortran

   do i=1,ngrid
       if(ok(i))then
           ! Get gas cell position
           x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
           do iSN=1,nSN
              ! Check if the cell lies within the SN radius
              dxx=x-xSN(iSN,1)
              dyy=y-xSN(iSN,2)
              dzz=z-xSN(iSN,3)
              dr_SN=dxx**2+dyy**2+dzz**2
              if(dr_SN.lt.rmax2)then
                 ! Do things with uold(ind_cell(i), :) and SN variables
              endif
          end do
       endif
   end do

Once again, the main idea is to select cells based on their position.

Updating cells and particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The final, and perhaps most important, operation that is required for
subgrid models is the ability to update the properties of cells and
particles.

For grid quantities (e.g., gas density or momentum), we just need to
update ``uold`` or ``unew`` with the right values.

.. warning::

   **Careful:** The choice of ``uold`` vs ``unew`` depends on where the
   subgrid model is called in the main ``amr_step`` loop.

   For example, ``thermal_feedback`` affects ``unew`` while
   ``kinetic_feedback`` affects ``uold``. This is because
   ``thermal_feedback`` is called *after* the ``set_unew`` call, while
   ``kinetic_feedback`` is called before.

Either way, the logic is similar: once the indices of the cells affected
by feedback are identified (in the ``indp`` array for thermal feedback
and in the ``indSN`` for kinetic feedback), we can directly change their
value. For example, at the end of the ``Sedov_blast`` routine, we have
the following loop:

.. code:: fortran

     do iSN=1,nSN
        if(vol_gas(iSN)==0d0)then
           u=vSN(iSN,1)
           v=vSN(iSN,2)
           w=vSN(iSN,3)
           if(indSN(iSN)>0)then
              uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas(iSN)
              uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas(iSN)*u
              uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas(iSN)*v
              uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas(iSN)*w
              uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas(iSN)*0.5d0*(u*u+v*v+w*w)+p_gas(iSN)
              if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_metal(iSN)
           endif
        endif
     end do

.. warning::

   These updates usually must be done for the **conservative**
   variables: this is because they describe extensive quantities, which
   can be added up (e.g., momentum, energy, etc). When developing a
   module, make sure that you are using the right variable set.

Implementation pitfalls
-----------------------

When writing a new subgrid model (or any extra module that deals with
cells and particles), the main task is to connect all the generic
operations discussed above in a single, consistent, piece of code.

Beyond just accessing properties, updating cells, and looping over
grids, this requires extra care to deal with the fact that the physical
model lives inside the rest of the code.

Units
~~~~~

One of the *most* bug-inducing part of writing subgrid models is dealing
with units. In RAMSES, the units can be accessed through the ``units``
subroutine in ``amr/units.f90``: this short routine defines a series of
``scale_X`` conversion factors, and most physical routines in RAMSES
start by calling it

.. code:: fortran

     ! Conversion factor from user units to cgs units
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

This defines for example the length conversion factor, ``scale_l``,
which is fixed for cosmological runs, and a namelist parameter
otherwise.

While this presents no fundamental difficulty, special care must be
taken to actually *use* these conversion factors.

For example, the ``sf_trelax`` namelist parameter is converted to **code
units** using

.. code:: fortran

     trel=sf_trelax*Myr2sec/scale_t ! relaxation timescale

Here, ``sf_trelax`` is set in the namelist in Myr, ``Myr2sec`` converts
it to seconds, and the ``scale_t`` division converts everything to code
units.

Empty cells
~~~~~~~~~~~

RAMSES reacts badly when cells are completely devoid of gas. Indeed,
when converting from conservative to primitive variable, we operate a
division by the density: schematically speaking,
``velocity = momentum / density``.

To avoid this issue, RAMSES introduces two main strategies: - using a
floor values when dividing by the density (usually called ``smallr`` in
the code), - making sure that subgrid models that remove gas from cells
**never** remove all of the gas.

The ``star_formation`` routine, for example, uses this strategy:

.. code:: fortran

   ! Security to prevent more than 90% of gas depletion
   if (mgas > 0.9d0*mcell) then
       nstar_corrected=int(0.9d0*mcell/mstar)
       mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar
       nstar(i)=nstar_corrected
   endif

TODO: others
~~~~~~~~~~~~

Supplementary exercises
-----------------------

We are now equipped to deal with *most* subgrid model requirements.
*You* can then try apply this in a series of problems. We will explore
this in more practical details in the problem session later in the week.

Star formation
~~~~~~~~~~~~~~

Complete the following metacode to describe the implementation of a star
formation model that uses only local criteria, where stars are allowed
to form in cells with a density exceeding ``density_star`` and
temperature below ``temperature_star``, which are expected to be in CGS.

.. code:: fortran

   subroutine star_formation()

       ! TO BE COMPLETED ... 
       
       ! getting physical info ... (from uold or unew ?, etc.)
       ! units, time, cell sizes, position, density etc
       ! decide how many stars you form 
       target_mstar = dt * density / tff * epsilon 
       ! Poisson sampling
       ! create particles and add to lists 
       ! remove matter from cell

In a second step, modify the code to implement a star formation model
that forms stars where the local density is above a threshold value and
the flow of gas is converging.

TODO: SN feedback?
~~~~~~~~~~~~~~~~~~
