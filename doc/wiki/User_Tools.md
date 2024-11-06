
# External tools

This page is designed to promote tools, datasets or other useful things that may be interesting to other RAMSES users but don't belong on the main RAMSES repository.

Feel free to edit this page and add your own tools with details on how to obtain them.

## Initial Conditions

### DICE

DICE is designed to set up one or more galaxies in isolation. It is included in RAMSES.

- [Code repository](https://bitbucket.org/vperret/dice/wiki/RAMSES_simulation).

### MUSIC

MUSIC generates cosmological initial conditions and can be used to set up RAMSES cosmo simulations.

- [Code repository](https://bitbucket.org/ohahn/music),
- [Doc page](./Music).

### MPgrafic

MPgrafic generates (large) cosmological initial conditions and can be used to set up RAMSES cosmo simulations.

- [Code repository](https://bitbucket.org/broukema/mpgrafic),
- [Reference](http://adsabs.harvard.edu/abs/2013ascl.soft04014P).

## Analysis and Post-Processing

### RAMSES tools
RAMSES comes with a decent number of Fortran routines and programs to extract information from RAMSES outputs.

You can find some documentation on the relevant [doc page](./RAMSES_utils).


### OSYRIS
OSIRIS is a simple python interface developed by Neil Vaytet to visualise RAMSES outputs
- [Code repository](https://osyris.readthedocs.io/en/stable/)

### YT

YT is a large community code supporting a number of simulation codes, including RAMSES.
- [Website](http://yt-project.org/)

### PYNBODY
Pynbody is a python interface developed by Andrew Pontzen to visualise particle data.
It also works for RAMSES outputs by turning cells into particles.
- [Code repository](https://github.com/pynbody/pynbody)


### Pymses
Pymses is an analysis library written in Python for RAMSES outputs.
- [Webpage](http://irfu.cea.fr/Projets/PYMSES/) (outdated)

### MERA
Mera is a Julia package developed by Manuel Behrendt to efficiently read/store/analyse RAMSES outputs.
- [Code repository](https://github.com/ManuelBehrendt/Mera.jl)

## Visualisation

### GLnemo2

GLnemo2 is an interactive visualisation 3D program using OpenGL. This software, developed by Jean-Charles Lambert, can help you visualise particles from many different codes, including RAMSES.

It works on your laptop or on larger servers.
- [Code repository](https://projets.lam.fr/projects/glnemo2)
