

This page is designed for RAMSES users to promote tools, datasets or other useful things that may be interesting to other RAMSES users but don't belong on the main RAMSES repository.

Feel free to edit this page and add your own tools with details on how to obtain them.

# Initial Conditions #

## DICE ##

DICE is designed to set up one or more galaxies in isolation. It is included in RAMSES. 

For more information, see: https://bitbucket.org/vperret/dice/wiki/RAMSES\ simulation

## MUSIC ##

MUSIC generates cosmological initial conditions and can be used to set up RAMSES cosmo simulations.

Code: https://bitbucket.org/ohahn/music

Wiki page: ./Music

## MPgrafic ##

MPgrafic generates (large) cosmological initial conditions and can be used to set up RAMSES cosmo simulations.

Code : https://bitbucket.org/broukema/mpgrafic

Reference: http://adsabs.harvard.edu/abs/2013ascl.soft04014P

# Analysis and Post-Processing #

## RAMSES tools ##
RAMSES comes with a decent number of Fortran routines and programs to extract information from RAMSES outputs. 

You can find some documentation on the [wiki page](./RAMSES\ utils)

## MERA ##
Mera is a Julia package developed by Manuel Behrendt to efficiently read/store/analyse RAMSES outputs: https://github.com/ManuelBehrendt/Mera.jl

## OSIRIS ##
OSIRIS is a simple python interface developed by Neil Vaytet to visualise RAMSES outputs: https://bitbucket.org/nvaytet/osiris

## PYNBODY ##
Pynbody is a python interface developed by Andrew Pontzen to visualise particle data. 

It also works for RAMSES outputs by turning cells into particles: https://github.com/pynbody/pynbody

## Pymses ##
Pymses is an analysis library written in Python for RAMSES outputs: http://irfu.cea.fr/Projets/PYMSES/

## YT ##

YT fills a similar role to Pymses. It is a large community code supporting a number of simulation codes, including RAMSES: http://yt-project.org/

## Hamu ##
Hamu lets you organise your simulation results and save Python analysis routine outputs automatically. 

It wraps around Pymses, although options such as YT are also possible: https://github.com/samgeen/Hamu

# Visualisation #

## GLnemo2 ##

GLnemo2 is an interactive visualisation 3D program using OpenGL. This software, developed by Jean-Charles Lambert, can help you visualise particles from many different codes, including RAMSES. 

It works on your laptop or on larger servers. Download the code here: https://projets.lam.fr/projects/glnemo2

## Hegelian ##

[Work in progress!] Hegelian is a 3D interactive visualisation with Python & OpenGL. This is a heavily work-in-progress project that was uploaded here to facilitate sharing of development code: https://github.com/samgeen/Hegelian