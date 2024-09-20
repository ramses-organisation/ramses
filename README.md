[1]: https://ramses-organisation.readthedocs.io/en/latest
[2]: https://bitbucket.org/rteyssie/ramses/wiki/AutoTests
[3]: http://www.ics.uzh.ch/~teyssier/ramses/RAMSES.html
[4]: https://bitbucket.org/rteyssie/ramses/wiki/ramses_ug.pdf
[5]: https://bitbucket.org/vperret/dice
[6]: https://bitbucket.org/ohahn/music
[7]: https://github.com/osyris-project/osyris
[8]: https://github.com/pynbody/pynbody
[9]: https://yt-project.org

![GitHub logo dark-mode-only](./doc/img/full_project_logo_dark.svg#gh-dark-mode-only)
![GitHub logo light-mode-only](./doc/img/full_project_logo.svg#gh-light-mode-only)

## This is the [ramses](https://github.com/ramses-organisation/ramses/) GitHub repository.

Ramses is an open source code to model astrophysical systems, featuring self-gravitating, magnetised, compressible, radiative fluid flows. It is based  on the Adaptive Mesh Refinement (AMR)  technique on a  fully-threaded graded octree.
[ramses](https://github.com/ramses-organisation/ramses/) is written in  Fortran 90 and is making intensive use of the Message Passing Interface (MPI) library.

You can go to the user's guide using [Read The Docs][1].

Download the code by cloning the git repository using
```
$ git clone https://github.com/ramses-organisation/ramses
```
You will get the latest stable version. To get the development branch, do the following
```
$ git clone --branch dev https://github.com/ramses-organisation/ramses
```
Please register also to the [mailing list](http://groups.google.com/group/ramses_users).

To generate idealised initial conditions of galaxies, check out [DICE][5].

To generate cosmological initial conditions, check out [MUSIC][6].

To visualize RAMSES data, we encourage you to use [YT][9], [OSYRIS][7] or [PYNBODY][8].
