#! /usr/bin/env python
#import matplotlib
#matplotlib.use('WXAgg')
from sys import argv
import math
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import glob
import re
import os

print argv[1]
nbins=int(argv[2])
print nbins

filename=glob.glob(argv[1]+'/sink*.out')
print filename
try:
    f = open(filename[0], 'r')
    strcase=' sinks'
except IOError:
    f = open('clump_masses.txt','r')
    strcase=' clumps'
lstr=f.readline()
print lstr
n=int(re.findall(r'\d+(?:\.\d+E-\d+)?',lstr)[0])
lstr=f.readline()
lstr=f.readline()
lstr=f.readline()

    
tot=0.
a=[]
zeros=0
for i in range(1,n+1):
    lstr=f.readline()
    m=float(re.findall(r'\d+(?:\.\d\d+)?',lstr)[1])
    print m
    tot=tot+m
    if m>0.:
        a.append(math.log(m,10.))
    else:
        zeros=zeros+1

n=n-zeros        
print str(tot)+' M_sol in '+str(n)+strcase+'. Average sink mass: ',str(tot/n)

histn,bins,patch = plt.hist(a,nbins,log=True,alpha=0.5)
spacing=(bins[nbins]-bins[0])/nbins
#histn=np.append(histnn,histnn[nbins-1])
#chabrier2005 IMF
log_x=np.arange(-3.,0,0.01)
log_x2=np.arange(0.,2,0.01)
y=0.093*np.exp(-1.*(log_x-np.log10(0.2))**2./(2.*0.55**2.))*tot/0.0795648*spacing
y2=0.041*10**(-1.35*log_x2)*tot/0.0795648*spacing
z=0.093*np.exp(-1.*(log_x-np.log10(0.2))**2./(2.*0.55**2.))*n/0.12834*spacing
z2=0.041*10**(-1.35*log_x2)*n/0.12834*spacing

plt.plot(log_x,y,'g--')
plt.plot(log_x2,y2,'g--')

plt.plot(log_x,z,'g')
plt.plot(log_x2,z2,'g')

#Kroupa IMF
log_v1=np.arange(-3.,np.log10(0.08),0.01)
log_v2=np.arange(np.log10(0.08),np.log10(0.5),0.01)
log_v3=np.arange(np.log10(0.5),2.,0.01)
w1=1./0.08*10**(0.7*log_v1)*n/3.04105*spacing
w2=10**(-.3*log_v2)*n/3.04105*spacing
w3=0.5*10**(-1.3*log_v3)*n/3.04105*spacing

plt.plot(log_v1,w1,'r')
plt.plot(log_v2,w2,'r')
plt.plot(log_v3,w3,'r')

plt.ylim([.8,100])

plt.savefig('plot.png',format='png')
os.system('display plot.png')
os.system('rm plot.png')
#plt.show()
