#!/usr/bin/env python
from flipper import *
import speckMisc
import pickle
import liteMapPol


def apoCorrection(patchDir,i):
    
    print 'ok routine a bit out to date as we now have two window functions'
    print 'also be careful if you wanna pad'
    
    Ra0Array= p['Ra0Array']
    Ra1Array= p['Ra1Array']
    Dec0Array = p['Dec0Array']
    Dec1Array = p['Dec1Array']
    
    m = liteMap.liteMapFromFits("%s/T_map_%03d_0"%(patchDir,i))
    area=numpy.float(m.Nx*m.Ny)
    area/=(m.Nx*m.Ny)
    m.data[:] = 1.0
    c=liteMapPol.initializeCosineWindow(m,cA['lenApod'],cA['pad'])
    m=liteMap.liteMapFromFits('auxMaps/finalWindow%03d_T.fits'%i)
    #m.plot(show=False,saveFig='plot/window%03d.png'%i)
    pylab.imshow(m.data,origin="down")
    pylab.savefig('auxMaps/window%03d.png'%i)
    pylab.clf()
    pylab.close()
    
    if p['doPadding']==True:
        m=liteMap.liteMapFromFits('auxMaps/finalWindow%03d_T_pad.fits'%i)
        #m.plot(show=False,saveFig='plot/window%03d.png'%i)
        pylab.imshow(m.data,origin="down")
        pylab.savefig('auxMaps/window%03d_pad.png'%i)
        pylab.clf()
        pylab.close()
        m=m.selectSubMap(Ra0Array[i],Ra1Array[i],Dec0Array[i],Dec1Array[i])
    
    factor= numpy.mean(m.data**4)/(numpy.mean(m.data**2))**2/numpy.mean(c.data[:])
    facs = [factor]#*area]
    print i,  numpy.sqrt(factor)#*area)
    return numpy.sqrt(facs)





specDir = 'spectra/'
dataDir= 'DataProducts/'



l_data,Cl_TB_data,error_TB_data=numpy.loadtxt(specDir+'/spectrum_TB_000.dat',unpack=True)
l_data,Cl_EB_data,error_EB_data=numpy.loadtxt(specDir+'/spectrum_EB_000.dat',unpack=True)
l_data,Cl_BB_data,error_BB_data=numpy.loadtxt(specDir+'/spectrum_BB_000.dat',unpack=True)




cl_TB_Sim=[]
cl_EB_Sim=[]
cl_BB_Sim=[]
Cov_TBTB=[]
Cov_EBEB=[]
Cov_BBBB=[]


for k in range(720):
    l,cl_TB,error=numpy.loadtxt('spectra_all/crossSpectrum/clBinCrossMean_TB_%03d.dat'%k,unpack=True)
    l,cl_EB,error=numpy.loadtxt('spectra_all/crossSpectrum/clBinCrossMean_EB_%03d.dat'%k,unpack=True)
    l,cl_BB,error=numpy.loadtxt('spectra_all/crossSpectrum/clBinCrossMean_BB_%03d.dat'%k,unpack=True)
    id=numpy.where(l>210)
    l=l[id]
    cl_TB=cl_TB[id]
    cl_EB=cl_EB[id]
    cl_BB=cl_BB[id]

    cl_TB_Sim+=[cl_TB]
    cl_EB_Sim+=[cl_EB]
    cl_BB_Sim+=[cl_BB]

    Cov_TBTB+=[cl_TB*cl_TB]
    Cov_EBEB+=[cl_EB*cl_EB]
    Cov_BBBB+=[cl_BB*cl_BB]
    




Cov_TBTB=numpy.mean(Cov_TBTB,axis=0)-numpy.mean(cl_TB_Sim,axis=0)**2
Cov_EBEB=numpy.mean(Cov_EBEB,axis=0)-numpy.mean(cl_EB_Sim,axis=0)**2
Cov_BBBB=numpy.mean(Cov_BBBB,axis=0)-numpy.mean(cl_BB_Sim,axis=0)**2



speckMisc.writeBinnedSpectrum(l_data[id],Cl_TB_data[id],numpy.sqrt(Cov_TBTB),dataDir+'/spectrum_TB_simError.dat')
speckMisc.writeBinnedSpectrum(l_data[id],Cl_EB_data[id],numpy.sqrt(Cov_EBEB),dataDir+'/spectrum_EB_simError.dat')
speckMisc.writeBinnedSpectrum(l_data[id],Cl_BB_data[id],numpy.sqrt(Cov_BBBB),dataDir+'/spectrum_BB_simError.dat')
