[1]: https://bitbucket.org/rteyssie/ramses/wiki/Content
[2]: https://bitbucket.org/rteyssie/ramses/wiki/AutoTests
[3]: http://www.ics.uzh.ch/~teyssier/ramses/RAMSES.html
[4]: https://bitbucket.org/rteyssie/ramses/wiki/ramses_ug.pdf
[5]: https://bitbucket.org/vperret/dice
[6]: https://bitbucket.org/ohahn/music
[7]: https://github.com/osyris-project/osyris
[8]: https://github.com/pynbody/pynbody


## This is the [ramses](https://github.com/ramses-organisation/ramses/) GitHub repository.

Ramses is an open source code to model astrophysical systems, featuring self-gravitating, magnetised, compressible, radiative fluid flows. It is based  on the Adaptive Mesh Refinement (AMR)  technique on a  fully-threaded graded octree. 
[ramses](https://github.com/ramses-organisation/ramses/) is written in  Fortran 90 and is making intensive use of the Message Passing Interface (MPI) library.

You can go to the user's guide using the [WIKI here][1].

Check regularly the [automatic test page][2].

Download the code by cloning the git repository using 
```
$ git clone https://github.com/ramses-organisation/ramses
```
You will get the latest stable version. To get the develop branch, do the following
```
$ git clone --branch develop https://github.com/ramses-organisation/ramses
```
Please register also to the [mailing list](http://groups.google.com/group/ramses_users).

To generate idealised initial conditions of galaxies, check out the [DICE][5] code.

To generate cosmological initial conditions, check out the [MUSIC][6] code.

To visualize RAMSES data, we encourage you to use the [OSYRIS][7] or the [PYNBODY][8] codes.
