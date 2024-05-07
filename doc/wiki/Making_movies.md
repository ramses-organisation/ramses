

There is an easy way to visually inspect your RAMSES simulations. If you specified `movie=.true.` in your namelist together with movie_params [Movie Params Wiki](./Movies)), then RAMSES produced directories named movieX, where X runs from 1 to 5. Each directory contains snapshots from the simulation. 

In order to produce a movie from these snapshots, one has to convert the Fortran binary maps into images and then combine them to a video file. This can be done with [RAMSES Animation Maker](https://bitbucket.org/biernacki/ram), which has full documentation and video examples.

If you have your own methods, please share them!