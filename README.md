[1]: https://bitbucket.org/rteyssie/ramses/wiki/Content
[2]: https://bitbucket.org/rteyssie/ramses/wiki/AutoTests
[3]: http://www.ics.uzh.ch/~teyssier/ramses/RAMSES.html
[4]: https://bitbucket.org/rteyssie/ramses/wiki/ramses_ug.pdf
[5]: https://bitbucket.org/vperret/dice
[6]: https://bitbucket.org/ohahn/music

## This is the [ramses](https://bitbucket.org/rteyssie/ramses) bitbucket repository.

Ramses is an open source code to model astrophysical systems, featuring self-gravitating, magnetised, compressible,
radiative fluid flows. It is based  on the Adaptive Mesh Refinement (AMR)  technique on a  fully-threaded graded octree. 
[ramses](https://bitbucket.org/rteyssie/ramses) is written in  Fortran 90 and is making intensive use of the Message 
Passing Interface (MPI) library.

You can go to the user's guide in [PDF here][4] and in the [WIKI here][1].

Check regularly the [automatic test page][2]. Visit the code web site [here][3].

Download the code by cloning the git repository using 
```
$ git clone https://bitbucket.org/rteyssie/ramses
```
This is the development branch. To get the stable release do the following
```
$ git clone --branch stable_17_09 https://bitbucket.org/rteyssie/ramses
```
Please register also to the [mailing list](http://groups.google.com/group/ramses_users).

To generate idealised initial conditions of galaxies, check out the [DICE][5] code.

To generate cosmological initial conditions, check out the [MUSIC][6] code.