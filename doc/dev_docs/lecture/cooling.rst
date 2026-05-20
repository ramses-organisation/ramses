.. math::


   \newcommand{\CH} {{\Lambda}}         % Cooling/heating function
   \newcommand{\CHp} {{\Lambda^\prime}} % Time derivative of coolheat function
   \newcommand{\Cool} {{\mathcal{L}}}   % Cooling function
   \newcommand{\coolU} {\rm{erg} \, \rm{cm}^{3} \, \rm{s}^{-1}}         % erg cm3 s-1
   \newcommand{\CoolZ} {{\mathcal{L}_Z}}% Metal cooling function
   \newcommand{\dt} {\Delta t}          % Discrete time interval
   \newcommand{\dtcool} {\Delta t_{\rm{cool}}}% Discrete time interval
   \newcommand{\dthyd} {\Delta t_{\rm{hydro}}}% Discrete time interval
   \newcommand{\ele}    {{\rm{e}}}
   \newcommand{\eps} {{\varepsilon}}    % Energy density symbol
   \newcommand{\ergs} {\rm{erg} \, \rm{s}^{-1}}                         %     erg s-1
   \newcommand{\eTconv} {\frac{(\gamma-1)\mh}{\rho \kb}}
   \newcommand{\Heat} {{\mathcal{H}}}   % Heating function
   \newcommand{\heisub} {\rm{He \scriptscriptstyle I}}
   \newcommand{\heiisub} {\rm{He \scriptscriptstyle II}}
   \newcommand{\heiiisub} {\rm{He \scriptscriptstyle III}}
   \newcommand{\hisub} {\rm{H \scriptscriptstyle I}}
   \newcommand{\hiisub} {\rm{H \scriptscriptstyle II}}
   \newcommand{\kb} {k_{\rm{B}}}        % Bolzmann constant
   \newcommand{\mh} {m_{\rm{H}}}        % Proton mass
   \newcommand{\nh} {n_{\rm{H}}}
   \newcommand{\nel} {n_{\rm{e}}}
   \newcommand{\nhe} {n_{\rm{He}}}
   \newcommand{\nhei} {n_{\heisub}}
   \newcommand{\nheii} {n_{\heiisub}}
   \newcommand{\nheiii} {n_{\heiiisub}}
   \newcommand{\nhi} {n_{\hisub}}
   \newcommand{\nhii} {n_{\hiisub}}
   \newcommand{\sm} {\rm{s}^{-1}}       % s-1
   \newcommand{\Tmu} {T_{\mu}}          % T_m
   \newcommand{\xhi} {x_{\rm{H \scriptscriptstyle I}}}
   \newcommand{\xhii} {x_{\rm{H \scriptscriptstyle II}}}
   \newcommand{\xhei} {x_{\rm{He \scriptscriptstyle I}}}
   \newcommand{\xheii} {x_{\rm{He \scriptscriptstyle II}}}
   \newcommand{\xheiii} {x_{\rm{He \scriptscriptstyle III}}}
   \newcommand{\Zsun} {Z_{\odot}}       % Solar metallicity



Radiative cooling and heating
=============================

.. contents::

.. admonition:: Further reading

   More information on cooling and heating in RAMSES can be found in `Rosdahl’s PhD
   thesis <https://theses.fr/2012LYO10075>`__, or in the `RAMSES-RT
   paper <https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R/abstract>`__.
   A classical paper on the topic is `Katz, Weinberg, & Hernquist
   (1996) <https://ui.adsabs.harvard.edu/abs/1996ApJS..105...19K/abstract>`__.


In most astrophysical contexts, gas dissipates thermal energy through
collisional processes, and this energy is carried away by radiation. Gas
may also gain energy from radiation, for example through
photo-ionisations. In RAMSES, these radiative cooling and heating terms
are accounted for in the energy conservation equation with the net
cooling rate :math:`\Lambda`:

.. math::
   \frac{\partial E}{\partial t} + \nabla \cdot \left( (E + P) \boldsymbol{u} \right) = -\rho \boldsymbol{u} \cdot \nabla \phi + \Lambda(\rho, e).

Here, :math:`t` is time, :math:`\boldsymbol{u}` is the
gas bulk velocity, :math:`\phi` is the gravitational potential, :math:`E = \frac{1}{2}\rho u^2 + e` is the gas
total energy, :math:`e` is the internal energy, and :math:`P = (\gamma -1)e`, with :math:`\gamma`
the ratio of specific heats.

The net cooling rate :math:`\Lambda` is a sum over many microscopic processes, involving collisions
between ions and electrons, photo-ionisations, etc. In principle, it is thus a function of the temperature and density of the gas,
of its detailed chemical composition and ionisation state, and of the local radiation field. The different cooling models implemented in RAMSES
mostly differ on how they approximate this function.

The following approximations are available in RAMSES:

1. **ISM cooling**: a simplified prescription for cold, dense interstellar medium conditions.
2. **Equilibrium cooling with a UVB** (default): assumes photo-ionisation equilibrium with a redshift-dependent UV background.
3. **Equilibrium cooling with Grackle**: same physical approximation using the external Grackle library.
4. **Non-equilibrium cooling without RHD**: tracks ion fractions explicitly with a homogeneous UVB.
5. **Non-equilibrium cooling with RHD**: full radiative transfer coupled to chemistry of H and He.

Below, we will detail the implementation of the default model (equilibrium cooling) and
of the non-equilibrium cooling (with a UVB) as illustrative examples of all cooling models implementations.

Step by step *default cooling model*
-------------------------------------

initialisations
~~~~~~~~~~~~~~~

Some models require initialisations, and these are done in ``init_time``
(in ``amr/init_time.f90``) which is called at the beginning of a
simulation (see code snippet below). For the default cooling model, this involves a call to
``set_model`` in ``hydro/cooling_module.f90``. The ``set_model`` routine
does two things: a) it initialises UV background parameters, and b), for
cosmological simulations, it computes the (homogeneous) temperature of
the Universe at the starting redshift of the simulation (typically read
from the cosmological initial conditions, but possible to override with
the aexp_ini namelist parameter).

.. code:: fortran

   subroutine adaptive_loop
      ...
      call init_time
      ...
      if(cooling) call set_table(dble(aexp))
      ...
      do ! Main time loop
         ...
         call amr_step(levelmin,1)
         ...
      end do
      ...
   end subroutine adaptive_loop


A second initialisation here is a call to ``set_table()`` (from ``hydro/cooling_module.f90``) before entering the
main time-evolution loop of RAMSES (see above). This is a first computation of the cooling
and heating rates tables.

Indeed, the default cooling module assumes photoionization equilibrium (PIE) with a redshift-dependent UVB. In such conditions,
the abundances of H and He ions are direct functions of temperature, density, and redshift only. In turn,
the cooling and heating rates are also reduced to being functions only of
temperature, density, redshift and metallicity, so that the cooling
function becomes :math:`\Lambda \simeq \Lambda(T,n_H, Z, z)`. The default
cooling method simplifies things further by assuming a linear dependence
with metallicity so that the cooling function becomes
:math:`\CH = (\Heat + \Cool + Z/\Zsun \, \CoolZ)\, \nh^2`, where

* :math:`\Heat(\Tmu,\nh)` is the heating rate from the UV background at the current redshift.
* :math:`\Cool(\Tmu,\nh)` is the cooling rate due to H and He (and electrons) assuming PIE with the UVB at the current redshift.
* :math:`\CoolZ(\Tmu,\nh)` is the cooling rate due to metals, at solar metallicity, and assuming PIE with the UVB at the current redshift.

With only two dimensions to the problem, it is much more efficient to pre-compute these rates and use linear interpolations than to calculate them on-the-fly
(exponents, powers).
The call to ``set_table`` computes and saves into tables these rates. As we'll see below, this is done at each coarse timestep to track the evolution of the UVB.



Filling in the rate tables with ``set_table``
~~~~~~~~~~~~~~~~~~~~~~~~~~

``set_table`` first computes the PIE ionisation state for H and He, with a simple iterative process that involves equating the
rates of photo-ionization, collisional ionization and recombination
(this is done in the ``cmp_chem_eq`` routine in
``hydro/cooling_module.f90``). Then the cooling and heating rates are computed and stored in tables
for a range of :math:`(\Tmu,\nh)`-bins, where :math:`T_\mu=T/\mu`
and :math:`\mu` is the mean particle mass in units of :math:`m_H`, :math:`\nh=X\rho/m_H` is
hydrogen number density, and :math:`X` is the hydrogen mass fraction in
the gas (a global constant, typically set to the value of :math:`0.76`).

**The abundances table:**

The abundances :math:`n_i`, of each of the 6
primordial species (*e*, HI, HII, HeI, HeII, HeIII) are computed and stored for a range of
(:math:`\Tmu, \nh`) values. Here, rates are used
for all the possible interactions involving these species - these are
functions of photoionization rates :math:`\Gamma_i`, :math:`T` and :math:`n_i`, and amount to
a closed set of equations that are converged iteratively to an
equilibrium solution, such that the creation rate equals the destruction
rate for each species. The species abundances also give the value of
:math:`\mu` per (:math:`\Tmu, \nh`)-bin, which can be retrieved by

.. math::

     \mu = \left[\, X(1+\xhii) + Y/4 (1+\xheii+2\xheiii) \, \right]^{-1},

where :math:`Y=1-X` is the helium mass fraction,
:math:`\xhii \equiv \nhii/\nh`, :math:`\xheii \equiv \nheii/\nhe` and
:math:`\xheiii \equiv \nheiii/\nhe`. The tabulation of the abundances
also provides a direct mapping between :math:`\Tmu` and :math:`T` which
is useful in generating the rest of the tables.

Note that these tables are written to each ramses output directory,
allowing the user to easily extract the equilibrium ionisation fractions
and :math:`\mu` in postprocessing (by performing the same kind of
interpolation for every cell as is done in RAMSES).

**The cooling rates table:**

Given the abundances, it is then a straightforward matter to calculate
and tabulate :math:`\Cool(\Tmu,\nh)` – the cooling rate is a sum of
bremsstrahlung, collisional excitation, collisional ionization,
recombination, dielectric recombination and Compton cooling rates,
all fitted analytic functions of temperature and abundances, and fetched
from various sources in the literature (e.g. Cen, 1992).

**The heating rates table:**

It is also straightforward to tabulate :math:`\Heat(\Tmu,\nh)`. Each bin
contains

.. math::


     \Heat=\sum_i \Heat_i n_i,

where the sum is over the primordial ion species, and :math:`\Heat_i`
are photoheating rates for the individual species (:math:`\nhi`,
:math:`\nhei`, :math:`\nheii`, and *e*).

**The metal-cooling-contributions table:**

RAMSES keeps a hard-coded table of a precomputed metal-cooling rate
contribution, :math:`\Cool_{Z}^{CIE}(T)`, which is the difference
between zero metallicity and solar metallicity cooling rates calculated
assuming collisional ionization equilibrium (CIE, i.e. chemical
equilibrium in zero ionizing radiation), with the CLOUDY software
package (Ferland et al., 1998). That is,

.. math::

     \Cool_{Z}^{CIE}(T) \; = \; \Cool_{Z}^{Cloudy}(T, \, Z=Z_{\odot}, \,  UV=0)
     \; - \; \Cool_{Z}^{Cloudy}(T, \, Z=0, \, UV=0).

These rates are computed for a photoionization-free environment so they
don’t depend on gas density. Using this, the photoionization equilibrium
(PIE) metal-cooling rates are approximated and tabulated as

.. math::

     \CoolZ(\Tmu,\nh) \; = \; \Cool_{Z}^{CIE}(T)
     \; \times \; f(T, \, \nh, \, z),

where :math:`f` is a dimensionless analytic function that corrects for
density :math:`\nh` and UV background photoionization at redshift
:math:`z`.

**Implementation:**
For the default cooling model, the implementation of these computations is all in ``hydro/cooling_module.f90``,
with the following structure.

.. code:: fortran

   subroutine set_table(aexp)
      ... ! define boundaries of tables in terms of nH and T
      call cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp) ! compute the tables.
   end subroutine set_table

   subroutine cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
      ! Compute radiative ionization and heating rates
      call set_rates(t_rad_spec,h_rad_spec,aexp)
      ! Create the table
      do i_n = myid,nbin_n,ncpu
         call iterate(i_n,t_rad_spec,h_rad_spec,nbin_T,aexp) ! compute cooling/heating for a range of temperatures
      end do
   end subroutine cmp_table

   subroutine iterate(i_n,t_rad_spec,h_rad_spec,nbin_T,aexp)
      nH = ... ! define value of nH from i_n
      do i_T = 1,nbin_T
         ! Compute cooling, heating and mean molecular weight
         T2=10d0**table%T2(i_T)      ! define the temperature of current bin
         call cmp_cooling(T2,nH ...) ! compute equilibrium cooling solution
         ...
         ! Compute cooling and heating derivatives
         T2_eps=10d0**(table%T2(i_T)+0.01d0)  ! define slightly larger temperature ("slightly = 0.01")
         call cmp_cooling(T2_eps,nH, ...)
         table%cool_prime(i_n,i_T)=(log10(cool_tot_eps)-log10(cool_tot))/0.01d0   ! evaluate the derivative and store in _prime tables.
         ...
         ! Compute metal contribution for solar metallicity
         call cmp_metals(T2,nH,mu,metal_tot,metal_prime,aexp)
         ...
      end do

   end subroutine iterate

   subroutine cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,cool_com,heat_com,mu_out,aexp,n_spec)

      ! Iteration to find mu (rates depend on T, but we only know T/mu)
      do while (error is too large)
         T = ... ! get T from T/mu and a guessed value of mu
         call cmp_chem_eq(T,nH,t_rad_spec,n_spec,n_TOT,mu)   ! get equilibrium solution at fixed T, including mu.
         ...
      end do

      ! Get equilibrium abundances
      n_E     = n_spec(1) ! electrons
      n_HI    = n_spec(2) ! H
      n_HII   = n_spec(3) ! H+
      ...
      ! Bremstrahlung
      cb1 = cool_bre(HI  ,T)*n_E*n_HII  /nH**2
      ...
      ! Ionization cooling
      ci1 = cool_ion(HI  ,T)*n_E*n_HI   /nH**2
      ... etc ...

   end subroutine cmp_cooling


Operator-splitting and temperature update
~~~~~~~~~~~~~~~~

For all cooling models, RAMSES uses **operator splitting** to solve the Euler equations in two
steps. In a first step, gravity is computed, and the gas is advected,
basically solving the Euler equations with :math:`\Lambda = 0`. In a
second **thermo-chemistry** step, the heating and cooling terms are
computed, and the internal energy is updated with :math:`\dot{e} = \Lambda`.
This strategy is shown in the code  snippet below taken from the main code loop ``amr_step``.
RAMSES first computes the hydrodynamics and updates ``uold`` with a call to ``godunov_fine``.
Then, RAMSES computes cooling on the updated state ``uold`` with a call to ``cooling_fine`` (more details below).

.. code:: fortran

   subroutine amr_step
      ...
      call godunov_fine ! do the hydro and update uold.
      ...
      call cooling_fine ! do the cooling and update uold.
      ...
   end subroutine amr_step


The thermo-chemistry involves the interaction of radiation and matter,
and the implementation differs if the code is used in radiation-hydrodynamics (RHD)
mode or in hydro (HD) mode illustrated above. With RHD, the
radiation is transfered accross the grid and the thermo-chemistry is computed on the fly.

In the general case, RAMSES implements a **semi-implicit scheme to evolve the temperature**
due to cooling during a timestep.
With :math:`\rho` the mass density, :math:`\gamma` the ratio of specific heats (usually given the
value of :math:`5/3` in RAMSES, corresponding to monatomic gas),
:math:`m_H` the proton mass, :math:`\kb` the Boltzmann constant,
and :math:`\mu m_H` the average particle mass, one can write the internal energy as

.. math::

     e = \frac{T}{\mu} \times \frac{\rho \kb}{(\gamma-1) m_H} ,

It is clear here that after advection, when the density (and metallicity) is fixed, evolving the
internal energy is equivalent to evolving the ratio :math:`\Tmu \equiv T/\mu`. This is what is done
in RAMSES, by solving the equation

.. math::

     \frac{\partial \Tmu}{\partial t} =
     \frac{(\gamma-1) \mh}{\rho \kb} \ \CH.


Starting at time :math:`t` with temperature :math:`\Tmu^{t}`, we compute
the evolved temperature :math:`\Tmu^{t+\Delta t}` with a semi-implicit
Euler formulation (See Press et al., 1992):

.. math::

   \Tmu^{t+\dt}= \Tmu^{t} + \frac{\CH K \dt}{1-\CHp K \dt},

where :math:`K=\eTconv` is the conversion factor between :math:`e`
and :math:`\Tmu` and
:math:`\CHp \equiv \frac{\partial \Lambda}{\partial \Tmu}` can be
estimated by finite-differencing the rate tables.




The thermochemistry is called in amr_step with a call to
``cooling_fine`` (in ``hydro/cooling_fine.f90)``, after the gravity
source term and hydro advection. This evolves the temperature of all
cells at the given level, over the current hydrodynamical timestep
length, :math:`\dthyd`.


Details of the ``cooling_fine`` subroutine.
~~~~~~~~~~~~~~~~~~~

In ``cooling_fine()``, cells at the given level are collected into
vectors of size ``NVECTOR`` and then each vector of cells is processed
in ``coolfine1()``. This is a generic operation which is useful for
vectorisation efficiency, and we illustrate it in the snippet below.
Notice that at the end of ``cooling_fine``, we again call ``set_table``
if we're at a coarse level (``ilevel==levelmin``) to prepare the next timestep.

.. code:: fortran

   subroutine cooling_fine(ilevel)

      ...

      ! Operator splitting step for cooling source term
      ! by vector sweeps
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
         call coolfine1(ind_grid,ngrid,ilevel)
      end do

      if((cooling.and..not.neq_chem.and..not.cooling_ism).and.ilevel==levelmin.and.cosmo)then
         call set_table(dble(aexp))
      endif

   end subroutine cooling_fine


The ``coolfine1`` subroutine basically collects the density,
temperature, and metallity for a batch of cells, in CGS units, and then calls
``solve_cooling()`` for those cells, which returns the (positive or
negative) change in the temperature, in Kelvin, over :math:`\dthyd`.
This is then converted back a change in internal energy densities, in
code units, and then ``uold`` is updated accordingly for each of the
cells.

``coolfine1`` follows these steps in order:

   1. **Unit conversion** — calls ``units()`` to obtain scale factors
      from code units to CGS.

   2. **loop over cells from ind_grid list**

      2.1. **Select leaf cells** — only cells with no children (``son == 0``)
      are processed. An index array ``ind_leaf`` is built; the routine
      returns immediately if ``nleaf == 0``.

      2.2. **Extract physical quantities** — reads ``nH``, metallicity
      ``Zsolar``, and total energy from ``uold``, using one loop per
      quantity.

      .. warning::

         Do **not** combine multiple field extractions into a single
         loop:

         .. code-block:: fortran

            ! BAD: hurts vectorisation
            do i = 1, nleaf
               nH(i)     = MAX(uold(ind_leaf(i), 1), smallr)
               T2(i)     = uold(ind_leaf(i), neul)
               Zsolar(i) = uold(ind_leaf(i), imetal) / nH(i) / 0.02d0
            end do

         Use a separate loop for each quantity so the compiler can
         vectorise the inner loop efficiently.

      2.3. **Compute thermal energy** — subtracts kinetic energy, radiation
      energy (if ``NENER > 0``), and magnetic energy from the total to
      isolate the thermal part. Converts to :math:`\Tmu` in Kelvin via
      ``scale_T2``.

      2.4. **Self-shielding** — if ``self_shielding`` is enabled, the UV
      background factor is attenuated exponentially at high densities with a boost factor:

      .. math::

         \mathrm{boost}(i)
         = \max\!\left(
             \exp\!\left(-\frac{\nh}{0.01\,\mathrm{cm}^{-3}}\right),
             10^{-20}
           \right).

      This is a simple "poor man's" model: dense gas that would be
      self-shielded sees a suppressed UVB. Note that since cooling/heating tables
      are already computed as a function of :math:`\Tmu` and :math:`n_H`, the boost
      is in fact applied to density. This is a good approximation, as the dependence
      of cooling and heating on the radiation intensity is really through a term :math:`\Gamma/n_H`,
      so decreasing :math:`\Gamma` or increasing :math:`n_H` by the same factor has the same effect.

      2.5. **Polytrope and temperature floor** — A minimum thermal energy is in general included. This
      may be because ``barotropic_eos`` is true, or because ``jeans_ncells>0``,
      or because ``T2_star >0``. This energy is expressed as a temperature floor (``T2min``),
      and needs to be computed and subtracted from the thermal energy before computing cooling.
      It is then inserted back in step 2.8. below.

      2.6. **Cooling integration** — calls ``solve_cooling``, which returns
      ``delta_T2``, the change in :math:`\Tmu` over the hydro timestep.

      2.7. **Delayed cooling** — if ``delayed_cooling`` is active, cooling
      is suppressed in blast-wave regions. When the delay tracer scalar
      (``idelay``) exceeds a threshold, ``delta_T2`` is clamped to zero
      (no cooling applied).

      2.8. **Energy update** — adds the polytrope floor energy back, as well as other energy terms (kinetic etc.) and
       writes the updated total energy to ``uold``:

       .. code-block:: fortran

          uold(...) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)


The call to ``solve_cooling`` (in ``hydro/cooling_module.f90``) is really the heart of
the thermochemistry. Here the temperature-change of every cell in the
vector is evolved, sub-cycling if needed with timestep lengths
:math:`\dtcool\propto\frac{T}{\Lambda}`. Within each
sub-step:

   1. Obtain :math:`\CH` and :math:`\CHp` by bilinear interpolation in
      the pre-computed 2D tables at the current :math:`(\Tmu, \nh)`.
   2. Compute the sub-step size :math:`\delta t` from the local cooling
      time :math:`\dtcool(\Tmu)`.
   3. Apply the semi-implicit update:

      .. math::

         \Tmu^{t+\delta t}
         = \Tmu^{t} + \frac{\CH K \delta t}{1 - \CHp K \delta t}.

   4. Advance the cell clock: :math:`t_{\mathrm{cell}} \leftarrow t_{\mathrm{cell}} + \delta t`.
      When :math:`t_{\mathrm{cell}} \geq \dthyd`, the cell is finished.

   At the end, return
   :math:`\Delta\Tmu = \Tmu^{\mathrm{final}} - \Tmu^{\mathrm{init}}`.












------------


Non-equilibrium cooling
-----------------------

This works quite similarly as the default equilibrium cooling, except,
here, the non-equilibrium fraction of ionised hydrogen, and, optionally,
HeII, HeIII, and neutral hydrogen (implicitly evolving molecular
hydrogen) is evolved and tracked in every cell. These ionisation
fractions are stored in passive scalars, usually right after the
metallicity scalar. Here, the ionization fractions are evolved along
with the temperature in each cooling sub-cycling step in a
quasi-implicit fashion (see the RAMSES-RT paper by Rosdahl et al. 2013).

The non-equilibrium cooling module was written specifically for
radiative transfer (see again the RAMSES-RT paper) and is found in
``rt/rt_cooling_module.f90``. It contains the routine
``rt_solve_cooling``, called instead of the default ``solve_cooling``
from ``cooling_fine`` if the non-equilibrium cooling is activated. In
this case, cooling fine updates not only the pressure variable in
``uold``, but also the passive scalars corresponding to the ionization
fractions of hydrogen and helium (and, if radiative transfer is also
activated, the momentum of the gas from radiation-gas interactions and
the photon fluxes and densities). The non-equilibrium cooling can only
be activated if RAMSES is compiled with the ``-RT`` flag. However, if
the flag is used, non-equilibrium cooling can still be activated without
radiative transfer, simply by compiling with zero radiation groups, and
using ``cooling=.true.`` and ``neq_chem=.true.`` in the
``COOLING_PARAMS`` namelist.

As for the equilibrium cooling module, cooling and heating rates are
tabulated, though in this case only against temperature (in Kelvin),
whereas the equilibrium cooling rates can be tabulated against both
temperature and gas density. The reason for this is that the cooling
rates depend on the ionisation fractions of the gas species, which are
direct functions of temperature in equilibrium, but not in
non-equilibrium (and tabulating against density is preferrable if
possible, since it reduces the computational cost compared to
multiplying the cooling rates with densities).

The non-equilibrium heating and cooling rates tables are initialised in
``rt_set_table`` (called during ``init_time``) and updated every coarse
timestep from ``amr_step``. For the cooling rates, only Compton cooling
actually needs to be updated, since the others are not
redshift-dependent (with equilibrium cooling all the cooling rates are
redshift-dependent trough the PIE ionization fractions).

With non-equilibrium cooling, the homogeneous UV background is more
flexible than the hardcoded variants in equilibrium cooling. It can be
read from files.

.. raw:: html

   <!-- ## GRACKLE cooling
   ...

   -->

ISM-cooling
-----------

Exercise: figure out how it works.

.. raw:: html

   <!--
   :::warning
   Describe 2 versions (1 with HD and 1 with RHD). Then exercise: (1) figure out how ISM cooling works, (2) implement a new version of cooling ...
   :::

   -->
