Implementation details
======================

This section of the documentation discusses the structure of the code and the implementation of the different modules. It aims to provide developers with a better understanding of the inner workings of RAMSES.

This documentation is written in the style of lectures adressing
specified topics and includes exercises. Its first version was
written in 2025 for the first RAMSES developer school organised by the
`RAMSES SNO <https://ramses.cnrs.fr>`__ (see credits below). 




.. toctree::
   :maxdepth: 1
   :caption: General overview of the code:

   general
   hydro
   gravity
   amr
   particles
   mpi
   cooling
   subgrid
   refinement


.. toctree::
   :maxdepth: 1
   :caption: Input/Output and User Interaction:

   parameters
   input
   output

.. toctree::
   :maxdepth: 1
   :caption: Methods and algoritms:

   godunov

Credits
-----

The general organization and conception of the first RAMSES Developer
School was done by: J. Blaizot (Coord.), N. Brucy, C. Cadiou, T. Colman, B. Commerçon, M. Farcy, M. Gonzalez
Rey, J. Rosdahl, J. Sorce, M. Trebitsch. The lecture notes were written by 
**J. Blaizot** (How to use Git, Radiative cooling and
heating, subgrid modelling utilities), **N. Brucy** (Refinement schemes and
implementation), **C. Cadiou** (MPI communications), **T. Colman** (overview of the code, hydrodynamics,
mesh data structures, MPI communications, namelist parameters, initial
conditions, outputting, Godunov solver, gravity), **B. Commercon** (Godunov
solver, gravity), **J. Rosdahl** (Radiative cooling and heating), **M. Trebitsch** (Particles, subgrid
modelling utilities).
