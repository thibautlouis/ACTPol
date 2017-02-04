#!/usr/bin/env python
from flipper import *
import speckMisc
import pickle
import liteMapPol






p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])

theoryFile = p['theoryFile']

specDir = 'spectra/'
patchDir = 'patches/'
plotDir= 'plot/'
seasons = p['seasonTags']
arrays = p['arrayTags']


if len(sys.argv)> 2:
    specDir = 'spectra_%s/'%sys.argv[2]
    patchDir = p['mapDir']+'patches_%s/'%sys.argv[2]
    plotDir = 'plot_%s/'%sys.argv[2]


l = os.listdir(patchDir)
nDivs = 0
nPatches = 0

for il in l:
    if 'all' in il:
        continue
    if 'T_map_%s_%s'%(arrays[0],seasons[0]) in il:
        nDivs += 1
    if 'T_map_%s_%s_0'%(arrays[0],seasons[0]) in il and '_0' in il[-2:]:
        nPatches += 1

seasonPairs = []
for i in xrange(len(seasons)):
    for j in xrange(len(seasons)):
        if arrays[0] == arrays[-1] and i>j: continue
        seasonPairs += [[seasons[i],seasons[j]]]



try:
    os.makedirs(plotDir)
except:
    pass


fields=['T','E','B']


clTh={}
clth_b={}

logXRange={}
logYRange={}
LinpowerOfL={}
meanCrossSpec={}
error_cross={}
XRange={}
YRange={}


LinpowerOfL['TT']=4
LinpowerOfL['EE']=2
LinpowerOfL['EB']=2
LinpowerOfL['TE']=2
LinpowerOfL['TB']=2
LinpowerOfL['BB']=2

logXRange['TT']=[100,5000]
logYRange['TT']=[3*10**-1,10000]
logXRange['EE']=[100,5000]
logYRange['EE']=[10**-1,5*10**1]
logXRange['BB']=[100,5000]
logYRange['BB']=[10**-4,10**1]

XRange['TT']=[0,4500]
YRange['TT']=[0,2.5*10**9]
XRange['EE']=[100,3500]
YRange['EE']=[-20,50]
XRange['BB']=[100,3500]
YRange['BB']=[-20,50]

XRange['EB']=[100,3500]
YRange['EB']=[-20,50]

XRange['TE']=[100,3500]
YRange['TE']=[-150,160]
XRange['TB']=[100,3500]
YRange['TB']=[-150,160]





X = numpy.loadtxt(theoryFile)

lmax=p['trimAtL']
#lmax=3000
lTh = X[:,0][:lmax]
clTh['TT'] = X[:,1][:lmax]
clTh['EE'] = X[:,2][:lmax]
clTh['BB'] = X[:,3][:lmax]
clTh['TE'] = X[:,4][:lmax]
clTh['TB'] = X[:,4][:lmax]*0
clTh['EB'] = X[:,4][:lmax]*0


Ap=23

count1=0
for l1 in fields:
    
    count1+=1
    count2=0
    
    for l2 in fields:
        count2+=1
        if count2<count1: continue
        
        print '%s%s'%(l1,l2)
        
        for sps in seasonPairs:
            spTag = '%sx%s'%(sps[0],sps[1])

            lbin,meanCrossSpec[l1+l2],error_cross[l1+l2]=numpy.loadtxt(specDir+'spectrum_%s%s_%sx%s_%s.dat'%(l1,l2,arrays[0],arrays[-1],spTag),unpack=True)
            if l1=='T':
                if l2=='E':
                    meanCrossSpec[l1+l2]=-meanCrossSpec[l1+l2]
        
        
            id=numpy.where((lbin>100))
            lbin=lbin[id]
            meanCrossSpec[l1+l2]=meanCrossSpec[l1+l2][id]
            error_cross[l1+l2]=error_cross[l1+l2][id]
            
            
            if l1==l2:
                
                if l1=='T':
                    pylab.semilogy()
                    pylab.errorbar(lbin,meanCrossSpec[l1+l2]*lbin**2/(2*numpy.pi),error_cross[l1+l2]*lbin**2/(2*numpy.pi),fmt='.',color='red')
                    pylab.plot(lTh, Ap*(lTh/3000.)**2,color='lightblue',label='Poisson %d'%Ap)
                    pylab.plot(lTh, clTh[l1+l2]+Ap*(lTh/3000.)**2,color='blue',label='LCDM+Poisson %d'%Ap)
                    pylab.plot(lTh, clTh[l1+l2],color='purple',label='LCDM')
                    pylab.ylim(10**(-1),10**(4))
                               
                    pylab.xlabel(r'$\ell$',fontsize=22)
                    pylab.ylabel(r'$\ell(\ell+1) C_\ell/(2 \pi)$',fontsize=22)
                    pylab.legend()
                    pylab.savefig(plotDir+'Log_spectrum_%s%s_%sx%s_%s.png'%(l1,l2,arrays[0],arrays[-1],spTag))
                    pylab.clf()
                    pylab.close()
                
                    
                    pylab.errorbar(lbin,meanCrossSpec[l1+l2]*lbin**(LinpowerOfL[l1+l2])/(2*numpy.pi),error_cross[l1+l2]*lbin**(LinpowerOfL[l1+l2])/(2*numpy.pi),fmt='.')
                    pylab.plot(lTh,(clTh[l1+l2]+Ap*(lTh/3000.)**2)*lTh**(LinpowerOfL[l1+l2]-2))
                    pylab.plot(lTh, clTh[l1+l2]*lTh**(LinpowerOfL[l1+l2]-2),color='purple',label='LCDM')
                    pylab.xlabel(r'$\ell$',fontsize=22)
                    pylab.ylabel(r'$\ell^{%s} C_\ell/(2 \pi)$'%LinpowerOfL[l1+l2],fontsize=22)
                    pylab.xlim(XRange[l1+l2])
                    pylab.ylim(YRange[l1+l2])
                    pylab.legend()
                    pylab.savefig(plotDir+'spectrum_%s%s_%sx%s_%s.png'%(l1,l2,arrays[0],arrays[-1],spTag))
                    pylab.clf()
                    pylab.close()
                
                
            
                else:
                    pylab.semilogy()
                    pylab.errorbar(lbin,meanCrossSpec[l1+l2]*lbin**2/(2*numpy.pi),error_cross[l1+l2]*lbin**2/(2*numpy.pi),fmt='.')
                    pylab.plot(lTh,clTh[l1+l2])
                    pylab.xlabel(r'$\ell$',fontsize=22)
                    pylab.ylabel(r'$\ell(\ell+1) C_\ell/(2 \pi)$',fontsize=22)
                    pylab.savefig(plotDir+'Log_spectrum_%s%s_%sx%s_%s.png'%(l1,l2,arrays[0],arrays[-1],spTag))
                    pylab.clf()
                    pylab.close()
        
        
                    pylab.errorbar(lbin,meanCrossSpec[l1+l2]*lbin**(LinpowerOfL[l1+l2])/(2*numpy.pi),error_cross[l1+l2]*lbin**(LinpowerOfL[l1+l2])/(2*numpy.pi),fmt='.')
                    pylab.plot(lTh,clTh[l1+l2]*lTh**(LinpowerOfL[l1+l2]-2))
                    pylab.xlabel(r'$\ell$',fontsize=22)
                    pylab.ylabel(r'$\ell^{%s} C_\ell/(2 \pi)$'%LinpowerOfL[l1+l2],fontsize=22)
                    #pylab.xlim(XRange[l1+l2])
                    #pylab.ylim(YRange[l1+l2])
                    pylab.savefig(plotDir+'spectrum_%s%s_%sx%s_%s.png'%(l1,l2,arrays[0],arrays[-1],spTag))
                    pylab.clf()
                    pylab.close()
            else:

                pylab.errorbar(lbin,meanCrossSpec[l1+l2]*lbin**(LinpowerOfL[l1+l2])/(2*numpy.pi),error_cross[l1+l2]*lbin**(LinpowerOfL[l1+l2])/(2*numpy.pi),fmt='.')
                pylab.plot(lTh,clTh[l1+l2]*lTh**(LinpowerOfL[l1+l2]-2))

                pylab.xlabel(r'$\ell$',fontsize=22)
                pylab.xlim(XRange[l1+l2])
                pylab.ylim(YRange[l1+l2])
                pylab.ylabel(r'$\ell^{%s} C_\ell/(2 \pi)$'%LinpowerOfL[l1+l2],fontsize=22)
                pylab.savefig(plotDir+'spectrum_%s%s_%sx%s_%s.png'%(l1,l2,arrays[0],arrays[-1],spTag))
                pylab.clf()
                pylab.close()
        

