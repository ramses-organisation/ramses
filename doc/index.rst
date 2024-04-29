Ramses Documentation
=====================

The [ramses][1] code is intended to be a versatile platform to develop
applications using  Adaptive Mesh  Refinement (AMR)  for computational
astrophysics.  The current implementation allows solving the classical
and relativistic Euler equations in presence of self-gravity, magnetic
field  and radiation  field.   The  [ramses][1] code  can  be used  on
massively  parallel  architectures,  if  properly linked  to  the  MPI
library.  It  can also  be used on  single processor  machines without
MPI.  Output  data are generated  using Fortran unformatted  files.  A
suite  of post-processing  routines  is delivered  within the  present
release,  allowing  the user  to  perform  a  simple analysis  of  the
generated output files.

[1]: https://bitbucket.org/rteyssie/ramses

[2]: https://bitbucket.org/ohahn/music

[3]: ./Content

[4]: ./AutoTests

[6]: ./User%20Tools

.. toctree::
  :maxdepth: 3

  wiki/Content.md

