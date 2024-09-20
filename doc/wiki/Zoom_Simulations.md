

# Zoom Simulations

It might be desirable to run a simulation to analyze the properties of an individual halo within a certain mass scale (e.g. the one of a galaxy cluster, a galaxy group or an individual galaxy). For this purpose it is of course essential to maximize the resolution on this object. However, within a cosmological framework it is necessary to simulate a sufficiently large region of the universe, in other words to choose a sufficiently large Boxsize (e.g.  100 Mpc/h), to obtain an accurate result. Because of the limitations of computational resources (RAM memory and computation time) a larger box will necessarily need to be simulated with lower resolution than a smaller one, or in turn a higher resolved astrophysical object will allow for a smaller simulation volume, than a lower resolved object. So in principle there is a trade-off between resolution and size of the simulated volume. 
A way out of this offers the concept of zoom simulations. The word zoom hereby means, that within the chosen simulation volume, the specific region of interest is simulated with a higher resolution, than the rest of the box (the zoom volume might for example be 0.001 times the box size).    

To make use of this concept, two major steps are required: First, to run a simulation with uniform resolution (unigrid run) of a box of the desired volume (as described above), and then to determine the region of interest. Given this information the second step, a RAMSES run with a higher resolution in the volume of interest, can be undertaken. This is called the zoom run. Before it can be carried out it is also required to create a new set of initial conditions, which contain the desired resolution in the zoom region.

[1. Run a zoom with RAMSES](./Ramses_zoom)

2. Zoom initial conditions with grafic files
  
[3. Generate zoom initial conditions with MUSIC](./Music)

