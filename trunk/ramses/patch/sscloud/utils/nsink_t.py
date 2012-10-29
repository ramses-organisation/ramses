#! /usr/bin/env python
from sys import argv
import glob
import re
import matplotlib.pyplot as plt
import os
import numpy as np


#create file list by simple ls command
os.system('ls run* > files.txt')

t=[]
n=[]
tf = open('files.txt','r')
for filename in tf:
    f = open(filename.rstrip(),'r')
    print filename.rstrip()
    gettime=False
    for line in f:
        if gettime==False:
            if 'produced' in line:
                nn=int(re.findall(r'\d+(?:\.\d\d+)?',line)[2])
                gettime=True

        else:
            if ' t=' in line:
                if len(re.findall(r'\d+(?:\.\d\d+)?',line))==11:
                    tt=float(re.findall(r'\d+(?:\.\d\d+)?',line)[1])*0.1**float(re.findall(r'\d+(?:\.\d\d+)?',line)[2])
                else:
                    tt=float(re.findall(r'\d+(?:\.\d\d+)?',line)[0])*0.1**float(re.findall(r'\d+(?:\.\d\d+)?',line)[1])
                t.append(tt)
                n.append(nn)
                gettime=False

nmax=n[::-1][0]

for i in range(nmax+1):
    done=False
    while done==False:
        if n.count(i)>1:
            t.remove(t[n.index(i)])        
            n.remove(i)
            done=False
        else:
            done=True

plt.plot(t,n,'o-',linewidth=0.2)
tn=[]
tf2 = open('files.txt','r')
for filename in tf2:
    f = open(filename.rstrip(),'r')
    done=False
    for line in f:
        if' t=' in line:
            if len(re.findall(r'\d+(?:\.\d\d+)?',line))==3:
                tt=float(re.findall(r'\d+(?:\.\d\d+)?',line)[0])*0.1**float(re.findall(r'\d+(?:\.\d\d+)?',line)[1])
            else:
                tt=float(re.findall(r'\d+(?:\.\d\d+)?',line)[0])
            if done==False:
                done=True
                if tt>t[0]: plt.plot([tt,tt],[0,n[::-1][0]],'r',linewidth=0.1)
                    

plt.savefig('plot.png',format='png')
os.system('display plot.png')
os.system('rm plot.png files.txt')
