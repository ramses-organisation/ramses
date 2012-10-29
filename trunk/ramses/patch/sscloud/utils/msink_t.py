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
m=[]

tf = open('files.txt','r')
for filename in tf:
    print filename.rstrip()
    f = open(filename.rstrip(),'r')
    data_point=False
    for line in f:
        if 'acc_rate'in line:
            data_point=True
            getmass=False
            gettime=False
            mm=0.

        elif data_point==True:
            if '========' in line and getmass==False and gettime==False:
                getmass=True
            elif '========' in line and getmass==True and gettime==False:
                getmass=False
                gettime=True
            elif gettime==True and ' t=' in line:
                ok=False
                if len(re.findall(r'\d+(?:\.\d\d+)?',line))==11:
                    tt=float(re.findall(r'\d+(?:\.\d\d+)?',line)[1])*0.1**float(re.findall(r'\d+(?:\.\d\d+)?',line)[2])
                    ok=True
                elif len(re.findall(r'\d+(?:\.\d\d+)?',line))==10:
                    tt=float(re.findall(r'\d+(?:\.\d\d+)?',line)[0])*0.1**float(re.findall(r'\d+(?:\.\d\d+)?',line)[1])
                    ok=True
                if ok==True:
                    t.append(tt)
                    m.append(mm)
                data_point=False
            elif getmass==True:
                mm=mm+float(re.findall(r'\d+(?:\.\d\d+)?',line)[1])

            
    if len(t) != len(m): m.pop()

clean=False
while clean==False:
    clean=True
    i=0
    while i<(len(t)-1) and clean==True:        
        if m[i+1]<m[i]:
            t.pop(i)
            m.pop(i)
            clean=False
        i=i+1    
clean=False
while clean==False:
    clean=True
    i=0
    while i<(len(t)-1) and clean==True:        
        if t[i+1]<t[i]:
            t.pop(i)
            m.pop(i)
            clean=False
        i=i+1    


plt.plot(t,m,'o-',markersize=0.4,linewidth=0.2)
plt.savefig('plot.png',format='png')
os.system('display plot.png')
os.system('rm plot.png files.txt')
