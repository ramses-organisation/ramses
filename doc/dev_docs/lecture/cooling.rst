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

Interlude: More Things
======================

-  Passive scalars (metallicity example, and refinement flag)
-  Cooling [Joki ?]

   -  General picture, interface, units, implementation choices in
      RAMSES (different options, …)
   -  operator splitting, source term, timesteps … (galaxy vs
      dense-core-collapse viewpoints…)
   -  **EXERCISE:** set time step of simulation to resolve cooling time
      (use some analytic approx).

      -  add field to hydro variable (``ivar+n``)
      -  store the cooling time
      -  output it … (in well understood units…)

   -  **EXERCISE:** introduce a source term evolving as t/tau…

-  ``NENER`` : what it is and what it can be used for.
-  Boundary conditions (review them)
-  Cosmology (:math:`a(t)`) [JS]

[TOC]

Radiative Cooling (and photo-heating)
=====================================

.. warning::

   More on cooling and heating in RAMSES can be found in `Rosdahl’s PhD
   thesis <https://theses.fr/2012LYO10075>`__, or in the `RAMSES-RT
   paper <https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R/abstract>`__.
   A classical paper on the topic is `Katz, Weinberg, & Hernquist
   (1996) <https://ui.adsabs.harvard.edu/abs/1996ApJS..105...19K/abstract>`__.

.. warning::

   -  General picture, interface, units, implementation choices in
      RAMSES (different options, …)
   -  operator splitting, source term, timesteps … (galaxy vs
      dense-core-collapse viewpoints…)
   -  **EXERCISE:** set time step of simulation to resolve cooling time
      (use some analytic approx). Another possibility is to add a
      heating term.

In most astrophysical contexts, gas dissipates thermal energy through
collisional processes, and this energy is carried away by radiation. Gas
may also gain energy from radiation, for example through
photo-ionisations. In RAMSES, these radiative cooling and heating terms
are accounted for in the energy conservation equation with the net
cooling term :math:`\Lambda`:

.. math::


   \frac{\partial E}{\partial t} + \nabla \cdot \left( (E + P) \boldsymbol{u} \right) = -\rho \boldsymbol{u} \cdot \nabla \phi + \Lambda(\rho, \varepsilon). 

Here, :math:`t` is time, :math:`E = \frac{1}{2}\rho u^2 + e` is the gas
internal energy, :math:`P` is pressure, :math:`\boldsymbol{u}` is the
gas bulk velocity, :math:`\phi` is the gravitational potential, and
:math:`\Lambda` is the net cooling term. The pressure is related to the
internal energy through :math:`P = (\gamma -1)e`, where :math:`\gamma`
is the ratio of specific heats.

RAMSES uses **operator splitting** to solve the Euler equations in
steps. In a first step, gravity is computed, and the gas is advected,
basically solving the Euler equations with :math:`\Lambda = 0`. In a
second **thermo-chemistry** step, the heating and cooling terms are
computed, and the energy is updated. The thermo-chemistry involves the
interaction of radiation and matter, and the implementation differs if
the code is used in radiation-hydrodynamics (RHD) mode or in hydro (HD)
mode. With RHD, the radiation is transfered accross the grid and the
thermo-chemistry is computed on the fly. In the HD case, the
thermo-chemistry is computed assuming the ionisation equilibrium
(possibly in the presence of a uniform UV background).

There are a number of different cooling implementations in RAMSES: -
equilibrium cooling with a UVB, using RAMSES (default). - equilibrium
cooling with a UVB, using Grackle. - ISM cooling (with or without RHD).
- non-equilibrium cooling without RHD. - non-equilibrium cooling with
RHD.

Below, we will first focus on the default equilibrium cooling and then
comment on the alternatives.

Interface
---------

Some models require initialisations, and these are done in ``init_time``
(in ``amr/init_time.f90``) which is called at the beginning of a
simulation. For the default cooling model, this involves a call to
``set_model`` in ``hydro/cooling_module.f90``. The ``set_model`` routine
does two things: a) it initialises UV background parameters, and b), for
cosmological simulations, it computes the (homogeneous) temperature of
the Universe at the starting redshift of the simulation (typically read
from the cosmological initial conditions, but possible to override with
the aexp_ini namelist parameter).

Furthermore, there is a call in ``adaptive_loop``, before entering the
main time-evolution loop of RAMSES, to ``set_table()`` in the default
``hydro/cooling_module.f90``. This is a first computation of the cooling
and heating rates tables (see below).

The thermochemistry is called in amr_step with a call to
``cooling_fine`` (in ``hydro/cooling_fine.f90)``, after the gravity
source term and hydro advection. This evolves the temperature of all
cells at the given level, over the current hydrodynamical timestep
length, :math:`\dthyd`.

In ``cooling_fine()``, cells at the given level are collected into
vectors of size ``NVECTOR`` and then each vector of cells is processed
in ``coolfine1()``. Basically, this routine collects the density,
temperature, and metallity for the cells, in CGS units, and then calls
``solve_cooling()`` for those cells, which returns the (positive or
negative) change in the temperature, in Kelvin, over :math:`\dthyd`.
This is then converted back a change in internal energy densities, in
code units, and then ``uold`` is updated accordingly for each of the
cells.

``solve_cooling`` in ``hydro/cooling_module.f90`` is really the heart of
the thermochemistry. Here the temperature-change of every cell in the
vector is evolved, sub-cycling if needed with timestep lengths
:math:`\dtcool\propto\frac{T}{\Lambda}`.

How to evolve the photo-ionization equilibrium temperature of a single cell
---------------------------------------------------------------------------

We take a cell of gas with a given internal energy density :math:`\eps`,
mass density :math:`\rho` and metallicity :math:`Z`, and we want to
update the cell energy over a time-step :math:`\Delta t`. :math:`Z` and
:math:`\rho` are updated elsewhere (during the advection step) and can
be taken as constants during the thermochemistry step. First, we convert
the internal energy to a temperature. The temperature of the cell is
given by

.. math::


     T = e \ \frac{(\gamma-1) m_H}{\rho \kb} \ \mu ,

where :math:`\gamma` is the ratio of specific heats (usually given the
value of :math:`5/3` in RAMSES, corresponding to monatomic gas),
:math:`m_H` the proton mass, :math:`\kb` the Boltzmann constant and
:math:`\mu` is the average mass per particle in the gas, in units of
:math:`m_H`. Since the ion species are not stored when using the default
cooling, :math:`\mu` is not readily available. Thus, the quantity that
is really evolved is

.. math::


   \Tmu \equiv \frac{T}{\mu},

which can be directly extracted from the variables stored in the cell.
We start at time :math:`t` with temperature :math:`\Tmu^{t}` and wish to
evolve it to :math:`\Tmu^{t+\Delta t}`, where the superscript is for
denoting *when* :math:`\Tmu` is evaluated.

Cooling and heating rates of the gas are functions of temperature,
density, redshift (via redshift-dependent UV background and CMB
radiation), metallicity, and the abundances of each primordial ion
species, :math:`\nhi`, :math:`\nhii`, :math:`\nhei`, :math:`\nheii`,
:math:`\nheii`, and :math:`\nel`. However, the default cooling module
assumes photoionization equilibrium (PIE), such that the primordial ion
abundances are direct functions of temperature, density, and redshift,
calculated with a simple iterative process that involves equating the
rates of photo-ionization, collisional ionization and recombination
(this is done in the ``cmp_chem_eq`` routine in
``hydro/cooling_module.f90``). The cooling and heating rates in the
equilibrium thermochemistry are then reduced to being functions only of
temperature, density, redshift and metallicity.

These rates are pre-computed and stored in tables every coarse time-step
for a range of :math:`(\Tmu,\nh)`-bins, where :math:`\nh=X\rho/m_H` is
hydrogen number density and :math:`X` is the hydrogen mass fraction in
the gas (a global constant, typically set to the value of :math:`0.76`).
The reason for pre-computing and storing in tables is that these cooling
and heating rates are numerically expensive to calculate on-the-fly
(exponents, powers), and table interpolation is much faster. A
redshift-dependent UV background of ionizing radiation is assumed, but
since it is homogeneous, i.e. exactly the same in every cell, it
suffices to recompute the cooling tables every coarse time-step, as the
redshift changes. Hence, the main operation in the ``solve_cooling``
routine is linear interpolation of tables that store logarithms of
cooling and heating rates.

Using :math:`\Tmu` and :math:`\nh` as look-up indexes, the following
rates, all given in [:math:`\coolU`], are fetched on-the-fly from the
precomputed tables: - Heating rate :math:`\Heat(\Tmu,\nh)`: The heating
contribution of the UV background at the current redshift. - Primordial
cooling rate :math:`\Cool(\Tmu,\nh)`, i.e. cooling rate of a mixture of
H and He (and electrons) of primordial composition. - Metal cooling rate
contribution :math:`\CoolZ(\Tmu,\nh)`, containing the
per-solar-metallicity cooling contribution of metals in the gas.

With these three rates in hand the temperature is updated by solving:

.. math::


     \frac{\partial \Tmu}{\partial t} = 
     \frac{(\gamma-1) \mh}{\rho \kb} \ \CH,

where
:math:`\CH\equiv\dot{e}=(\Heat + \Cool + Z/\Zsun \, \CoolZ)\, \nh^2`.
The update is done with a semi-implicit Euler formulation (See Press et
al., 1992):

.. math::

     
   \Tmu^{t+\dt}= \Tmu^{t} + \frac{\CH K \dt}{1-\CHp K \dt},

where :math:`K=\eTconv` is the conversion factor between :math:`\eps`
and :math:`\Tmu` and
:math:`\CHp \equiv \frac{\partial \Lambda}{\partial \Tmu}` can be
estimated by finite-differencing the rate tables.

So the integration of :math:`\Tmu` over the time-step is quick and
painless, and basically just involves interpolations of the rate-tables
to retrieve :math:`\Cool`, :math:`\Heat` and :math:`\CoolZ` for the
given temperature and density. The expensive calculations happen between
the time-steps, where the rate tables are updated, the advantage being
they are done a lot less frequently than if they were done on-the-fly
for every cell.

Filling in the rate tables
~~~~~~~~~~~~~~~~~~~~~~~~~~

First photoionization rates :math:`\Gamma_i \; [\sm]` and photoheating
rates :math:`\Heat_i \; [\ergs]` are calculated for HI, HeI, HeII and
*e* (Compton heating in the case of electrons). These rates depend
*only* on the UV background, so they are homogeneous functions of
redshift, i.e. apply to every cell.

The abundances table:
^^^^^^^^^^^^^^^^^^^^^

Then, since PIE is assumed, calculation and tabulation is performed per
(:math:`\Tmu, \nh`)-bin of the abundances :math:`n_i`, of each of the 6
primordial species (*e*, HI, HII, HeI, HeII, HeIII). Here rates are used
for all the possible interactions involving these species - these are
functions of :math:`\Gamma_i`, :math:`T` and :math:`n_i`, and amount to
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

The cooling rates table:
^^^^^^^^^^^^^^^^^^^^^^^^

Given the abundances, it is then a straightforward matter to calculate
and tabulate :math:`\Cool(\Tmu,\nh)` – the cooling rate is a sum of
bremsstrahlung-, collisional excitation-, collisional ionization-,
recombination-, dielectric recombination- and Compton- cooling rates,
all fitted analytic functions of temperature and abundances, and fetched
from various sources in the literature (e.g. Cen, 1992).

The heating rates table:
^^^^^^^^^^^^^^^^^^^^^^^^

It is also straightforward to tabulate :math:`\Heat(\Tmu,\nh)`. Each bin
contains

.. math::


     \Heat=\sum_i \Heat_i n_i,

where the sum is over the primordial ion species, and :math:`\Heat_i`
are photoheating rates for the individual species (:math:`\nhi`,
:math:`\nhei`, :math:`\nheii`, and *e*).

The metal-cooling-contributions table:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
read from files

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
