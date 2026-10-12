3 Physical boundaries [TODO]
============================

Physical boundaries around the computational domain
---------------------------------------------------

The strategy to set up boundary conditions in RAMSES is based on using
“ghost regions” outside the computational domain, where flow variables
are carefully specified in order to mimic the effect of the chosen type
of boundary. While we go deeper into the subject of boundary condition
in a later section, we explain here how boundary grids are positioned
with respect to the main computational domain.
