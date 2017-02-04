#!/usr/bin/env python
from flipper import *
import speckMisc

lmin=numpy.zeros(100)
lmax=numpy.zeros(100)
lcenter=numpy.zeros(100)


lmin[0],lmax[0],lcenter[0]=0.0,110,110./2

for i in range(99):
    
    
    if lmax[i]<2000:
        lmin[i+1]=lmax[i]+1
        lmax[i+1]=lmax[i]+50
    if ((lmax[i]>=2000) & (lmax[i]<2500)):
        lmin[i+1]=lmax[i]+1
        lmax[i+1]=lmax[i]+100
    if ((lmax[i]>=2500) & (lmax[i]<3000)):
        lmin[i+1]=lmax[i]+1
        lmax[i+1]=lmax[i]+200
    if ((lmax[i]>=3000) & (lmax[i]<6000)):
        lmin[i+1]=lmax[i]+1
        lmax[i+1]=lmax[i]+400
    if lmax[i]>=6000:
        lmin[i+1]=lmax[i]+1
        lmax[i+1]=lmax[i]+800
    lcenter[i+1]=(lmax[i+1]+lmin[i+1])/2.


    print lmin[i],lmax[i],lcenter[i]

pylab.plot(lmin)
pylab.plot(lmax)
pylab.plot(lcenter)
pylab.show()

speckMisc.writeBinnedSpectrum(lmin,lmax,lcenter,'BIN_ACTPOL_50_5')
