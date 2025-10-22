.. meta::
    :keywords: ramses, AMR, Adaptive Mesh  Refinement, self-gravity, magnetic field, MPI, galaxy, simulation, ISM

######################
Ramses Documentation
######################

The `Ramses` code is intended to be a versatile platform to develop
applications using  Adaptive Mesh  Refinement (AMR)  for computational
astrophysics.  The current implementation allows solving the classical
and relativistic Euler equations in presence of self-gravity, magnetic
field  and radiation  field.   The  `ramses` code  can  be used  on
massively  parallel  architectures,  if  properly linked  to  the  MPI
library.  It  can also  be used on  single processor  machines without
MPI.  Output  data are generated  using Fortran unformatted  files.  A
suite  of post-processing  routines  is delivered  within the  present
release,  allowing  the user  to  perform  a  simple analysis  of  the
generated output files. A pdf version of this documentation can be found
`here <Ramses.pdf>`_.

**New users of Ramses are invited to follow the ramses's** `tutorials <https://ramses-tutorials.readthedocs.io/en/latest/>`_,
**that cover the basics to setup a simulation and help with getting started
with some common applications.**


Table of Contents
*****************

.. toctree::
  :caption: User Documentation
  :maxdepth: 2

  wiki/Start.md
  wiki/Runtime_Parameters.md
  wiki/Advanced_Simulations.md
  wiki/Testing.md
  wiki/User_Tools.md

.. toctree::
  :caption: Developer Documentation
  :maxdepth: 2

  dev_docs/developer_guide.md
  dev_docs/implementation.md
  dev_docs/docs.md
  dev_docs/contributors.md

.. |github tag| image:: https://img.shields.io/badge/GitHub-black.svg?style=flat&logo=github
    :target: https://github.com/ramses-organisation/ramses
    :alt: ramses GitHub
