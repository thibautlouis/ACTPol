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




p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])



specDir = 'spectra/'
patchDir = 'patches/'
dataDir= 'DataProducts/'

if len(sys.argv)> 2:
    specDir = 'spectra_%s/'%sys.argv[2]
    patchDir = 'patches_%s/'%sys.argv[2]
    dataDir= 'DataProducts_%s/'%sys.argv[2]



try:
    os.makedirs(dataDir)
except:
    pass



l = os.listdir(patchDir)
nDivs = 0
nPatches = 0

for il in l:
    if 'all' in il:
        continue
    if 'T_map_000' in il:
        nDivs += 1
    if 'T_map_0' in il and '_0' in il[-2:]:
        nPatches += 1

cA = p['cosineApodization']
apoCorr =1.
fac = 2.0

apoCorr=numpy.zeros(nPatches)
for i in range(nPatches):
    
    apoCorr[i] = apoCorrection(patchDir,i)

meanCrossSpec={}
Theta={}
meanAutoSpec={}
meanNoise={}
binWeight={}
error={}

fields=['T','E']

count1=0
for l1 in fields:
    count1+=1
    count2=0
    
    for l2 in fields:
        count2+=1
        if count2<count1: continue
        
        for i in range(nPatches):
            
            lbin,meanCrossSpec[l1+l2,i],binWeight[l1+l2,i] = numpy.loadtxt('%s/clBinCrossMean_%s%s_%03d.dat'%(specDir,l1,l2,i),unpack=True)
            lbin,meanAutoSpec[l1+l2,i],binWeight[l1+l2,i] = numpy.loadtxt('%s/clBinAutoMean_%s%s_%03d.dat'%(specDir,l1,l2,i),unpack=True)
            meanNoise[l1+l2,i]=meanAutoSpec[l1+l2,i]-meanCrossSpec[l1+l2,i]
            
            id=numpy.where((lbin>210))
            lbin=lbin[id]
            meanCrossSpec[l1+l2,i]=meanCrossSpec[l1+l2,i][id]
            binWeight[l1+l2,i]=binWeight[l1+l2,i][id]
            meanAutoSpec[l1+l2,i]=meanAutoSpec[l1+l2,i][id]
            meanNoise[l1+l2,i]=meanNoise[l1+l2,i][id]
            

            meanAutoSpec[l2+l1,i]=meanAutoSpec[l1+l2,i]
            meanCrossSpec[l2+l1,i]=meanCrossSpec[l1+l2,i]
            meanNoise[l2+l1,i]=meanNoise[l1+l2,i]
            


Nbin=len(meanNoise[l2+l1,i])



theoryFile = p['theoryFile']


X= numpy.loadtxt(theoryFile)
lTh = X[:,0][:9000]
clTh_TT= X[:,1][:9000]
clTh_EE= X[:,2][:9000]
clTh_TE= X[:,3][:9000]

Bbl_T=pickle.load(open('mcm/Bbl.mcm_mask10amin_000.pkl',mode="r"))
Bbl_Pol=pickle.load(open('mcm/Bbl.mcm_mask10amin_000_cos2.pkl',mode="r"))

numpy.savetxt(dataDir+'/BblMean.dat', Bbl_T)
numpy.savetxt(dataDir+'/BblMean_Pol.dat',  Bbl_Pol)

Bbl_T=Bbl_T[:,:8999]
Bbl_Pol=Bbl_Pol[:,:8999]


print Bbl_T.shape
print clTh_TT.shape
Cb_th_TT=numpy.dot(Bbl_T,clTh_TT)
Cb_th_EE=numpy.dot(Bbl_Pol,clTh_EE)
Cb_th_TE=numpy.dot(Bbl_Pol,clTh_TE)


l,Cl_Mean_TT,error_Mean_TT=numpy.loadtxt('spectra_all/spectrum_TT.dat',unpack=True)
l,Cl_Mean_EE,error_Mean_EE=numpy.loadtxt('spectra_all/spectrum_EE.dat',unpack=True)
l,Cl_Mean_TE,error_Mean_TE=numpy.loadtxt('spectra_all/spectrum_TE.dat',unpack=True)

id=numpy.where((l>210))
l=l[id]
Cl_Mean_TT=Cl_Mean_TT[id]
error_Mean_TT=error_Mean_TT[id]
Cl_Mean_EE=Cl_Mean_EE[id]
error_Mean_EE=error_Mean_EE[id]
Cl_Mean_TE=Cl_Mean_TE[id]
error_Mean_TE=error_Mean_TE[id]
Cb_th_TT=Cb_th_TT[id]
Cb_th_EE=Cb_th_EE[id]
Cb_th_TE=Cb_th_TE[id]

print l


pylab.errorbar(l,Cl_Mean_TT*l**4/(2*numpy.pi),error_Mean_TT*l**4/(2*numpy.pi))
pylab.plot(l,Cb_th_TT*l**2)
pylab.show()

pylab.errorbar(l,Cl_Mean_EE*l**2/(2*numpy.pi),error_Mean_EE*l**2/(2*numpy.pi))
pylab.plot(l,Cb_th_EE)
pylab.show()

pylab.errorbar(l,Cl_Mean_TE*l**2/(2*numpy.pi),error_Mean_TE*l**2/(2*numpy.pi))
pylab.plot(l,Cb_th_TE)
pylab.show()


c=0
count1=0
for l1 in fields:
    count1+=1
    count2=0
    
    for l2 in fields:
        count2+=1
        if count2<count1: continue
        
        count3=0
        for l3 in fields:
            count3+=1
            count4=0
            
            for l4 in fields:
                count4+=1
                
                if count4<count3: continue
                
                for i in range(nPatches):
                    
                    print l1,l2,l2,l4
                    
                    Theta[l1+l2+l3+l4] =  meanCrossSpec[l1+l3,i]*meanCrossSpec[l2+l4,i]+meanCrossSpec[l1+l4,i]*meanCrossSpec[l2+l3,i]
                    Theta[l1+l2+l3+l4] += 1./nDivs*(meanCrossSpec[l1+l3,i]*meanNoise[l2+l4,i]+meanCrossSpec[l2+l4,i]*meanNoise[l1+l3,i]+meanCrossSpec[l1+l4,i]*meanNoise[l2+l3,i]+meanCrossSpec[l2+l3,i]*meanNoise[l1+l4,i])
                    Theta[l1+l2+l3+l4] += 1./(nDivs*(nDivs-1))*(meanNoise[l1+l3,i]*meanNoise[l2+l4,i]+meanNoise[l1+l4,i]*meanNoise[l2+l3,i])
                    
                    Theta[l1+l2+l3+l4] /= (fac*binWeight['EE',i])
                    Theta[l1+l2+l3+l4] *=apoCorr[i]**2

Ncross=3
CovMat=numpy.zeros((Ncross*Nbin,Ncross*Nbin))

for i in range(Nbin):
    CovMat[i,i]=Theta['TTTT'][i]
    CovMat[i,Nbin+i]=Theta['TTTE'][i]
    CovMat[i,2*Nbin+i]=Theta['TTEE'][i]

    CovMat[Nbin+i,i]=Theta['TTTE'][i]
    CovMat[Nbin+i,Nbin+i]=Theta['TETE'][i]
    CovMat[Nbin+i,2*Nbin+i]=Theta['TEEE'][i]

    CovMat[2*Nbin+i,i]=Theta['TTEE'][i]
    CovMat[2*Nbin+i,Nbin+i]=Theta['TEEE'][i]
    CovMat[2*Nbin+i,2*Nbin+i]=Theta['EEEE'][i]



print numpy.linalg.eigvals(CovMat)
pylab.matshow(numpy.log(numpy.abs(CovMat)))
pylab.show()

numpy.savetxt(dataDir+'/CovMat.dat', CovMat)

l_data,Cl_TT_data,error_TT_data=numpy.loadtxt(specDir+'/spectrum_TT_000.dat',unpack=True)
l_data,Cl_TE_data,error_TE_data=numpy.loadtxt(specDir+'/spectrum_TE_000.dat',unpack=True)
l_data,Cl_EE_data,error_EE_data=numpy.loadtxt(specDir+'/spectrum_EE_000.dat',unpack=True)

id=numpy.where(l_data>210)

speckMisc.writeBinnedSpectrum(l_data[id],Cl_TT_data[id],error_TT_data[id],dataDir+'/spectrum_TT.dat')
speckMisc.writeBinnedSpectrum(l_data[id],Cl_TE_data[id],error_TE_data[id],dataDir+'/spectrum_TE.dat')
speckMisc.writeBinnedSpectrum(l_data[id],Cl_EE_data[id],error_EE_data[id],dataDir+'/spectrum_EE.dat')






cl_TT_Sim=[]
cl_TE_Sim=[]
cl_EE_Sim=[]
Cov_TTTT=[]
Cov_TETE=[]
Cov_EEEE=[]
Cov_TTTE=[]
Cov_TTEE=[]
Cov_TEEE=[]

for k in range(720):
    l,cl_TT,error=numpy.loadtxt('spectra_all/crossSpectrum/clBinCrossMean_TT_%03d.dat'%k,unpack=True)
    l,cl_TE,error=numpy.loadtxt('spectra_all/crossSpectrum/clBinCrossMean_TE_%03d.dat'%k,unpack=True)
    l,cl_EE,error=numpy.loadtxt('spectra_all/crossSpectrum/clBinCrossMean_EE_%03d.dat'%k,unpack=True)
    id=numpy.where(l>210)
    l=l[id]
    cl_TT=cl_TT[id]
    cl_TE=cl_TE[id]
    cl_EE=cl_EE[id]

    cl_TT_Sim+=[cl_TT]
    cl_TE_Sim+=[cl_TE]
    cl_EE_Sim+=[cl_EE]

    Cov_TTTT+=[cl_TT*cl_TT]
    Cov_TETE+=[cl_TE*cl_TE]
    Cov_EEEE+=[cl_EE*cl_EE]
    
    Cov_TTTE+=[cl_TT*cl_TE]
    Cov_TTEE+=[cl_TT*cl_EE]
    Cov_TEEE+=[cl_TE*cl_EE]

cov = numpy.cov(cl_TT_Sim,rowvar=0,bias=1)

nBins = len(l)
i, j  = numpy.mgrid[0:nBins,0:nBins]
i = i.flatten()
j = j.flatten()
covNorm = (cov).copy()
#print numpy.diag(cov)
covNorm[i[:],j[:]] = cov[i[:],j[:]]/numpy.sqrt(cov[i[:],i[:]]*cov[j[:],j[:]])

for i in range(nBins):
    covNorm[i,i]=1-covNorm[i,i]


pylab.matshow(covNorm)
pylab.colorbar()
pylab.xlabel('Bin index')
pylab.ylabel('Bin index')
pylab.savefig('covarianceMatrixNorm_TT.png')
pylab.clf()
pylab.close()

cov = numpy.cov(cl_EE_Sim,rowvar=0,bias=1)

nBins = len(l)
i, j  = numpy.mgrid[0:nBins,0:nBins]
i = i.flatten()
j = j.flatten()
covNorm = (cov).copy()
#print numpy.diag(cov)
covNorm[i[:],j[:]] = cov[i[:],j[:]]/numpy.sqrt(cov[i[:],i[:]]*cov[j[:],j[:]])
print covNorm

for i in range(nBins):
    covNorm[i,i]=1-covNorm[i,i]


pylab.matshow(covNorm)
pylab.colorbar()
pylab.xlabel('Bin index')
pylab.ylabel('Bin index')
pylab.savefig('covarianceMatrixNorm_EE.png')
pylab.clf()
pylab.close()


cov = numpy.cov(cl_TE_Sim,rowvar=0,bias=1)

nBins = len(l)
i, j  = numpy.mgrid[0:nBins,0:nBins]
i = i.flatten()
j = j.flatten()
covNorm = (cov).copy()
#print numpy.diag(cov)
covNorm[i[:],j[:]] = cov[i[:],j[:]]/numpy.sqrt(cov[i[:],i[:]]*cov[j[:],j[:]])
print covNorm

for i in range(nBins):
    covNorm[i,i]=1-covNorm[i,i]


pylab.matshow(covNorm)
pylab.colorbar()
pylab.xlabel('Bin index')
pylab.ylabel('Bin index')
pylab.savefig('covarianceMatrixNorm_TE.png')
pylab.clf()
pylab.close()



Cov_TTTT=numpy.mean(Cov_TTTT,axis=0)-numpy.mean(cl_TT_Sim,axis=0)**2
Cov_TETE=numpy.mean(Cov_TETE,axis=0)-numpy.mean(cl_TE_Sim,axis=0)**2
Cov_EEEE=numpy.mean(Cov_EEEE,axis=0)-numpy.mean(cl_EE_Sim,axis=0)**2
Cov_TTTE=numpy.mean(Cov_TTTE,axis=0)-numpy.mean(cl_TT_Sim,axis=0)*numpy.mean(cl_TE_Sim,axis=0)
Cov_TTEE=numpy.mean(Cov_TTEE,axis=0)-numpy.mean(cl_TT_Sim,axis=0)*numpy.mean(cl_EE_Sim,axis=0)
Cov_TEEE=numpy.mean(Cov_TEEE,axis=0)-numpy.mean(cl_EE_Sim,axis=0)*numpy.mean(cl_TE_Sim,axis=0)


for i in range(Nbin):
    CovMat[i,i]=Cov_TTTT[i]
    CovMat[i,Nbin+i]=Cov_TTTE[i]
    CovMat[i,2*Nbin+i]=Cov_TTEE[i]
    
    CovMat[Nbin+i,i]=Cov_TTTE[i]
    CovMat[Nbin+i,Nbin+i]=Cov_TETE[i]
    CovMat[Nbin+i,2*Nbin+i]=Cov_TEEE[i]
    
    CovMat[2*Nbin+i,i]=Cov_TTEE[i]
    CovMat[2*Nbin+i,Nbin+i]=Cov_TEEE[i]
    CovMat[2*Nbin+i,2*Nbin+i]=Cov_EEEE[i]


speckMisc.writeBinnedSpectrum(l_data[id],Cl_TT_data[id],numpy.sqrt(Cov_TTTT),dataDir+'/spectrum_TT_simError.dat')
speckMisc.writeBinnedSpectrum(l_data[id],Cl_TE_data[id],numpy.sqrt(Cov_TETE),dataDir+'/spectrum_TE_simError.dat')
speckMisc.writeBinnedSpectrum(l_data[id],Cl_EE_data[id],numpy.sqrt(Cov_EEEE),dataDir+'/spectrum_EE_simError.dat')


numpy.savetxt(dataDir+'/CovMat_Sim.dat', CovMat)


print numpy.linalg.eigvals(CovMat)
pylab.matshow(numpy.log(numpy.abs(CovMat)))
pylab.show()





pylab.semilogy()
pylab.plot(l,Theta['TTTT'])
pylab.plot(l,Cov_TTTT)
pylab.show()

pylab.semilogy()
pylab.plot(l,Theta['TETE'])
pylab.plot(l,Cov_TETE)
pylab.show()

pylab.semilogy()
pylab.plot(l,Theta['EEEE'])
pylab.plot(l,Cov_EEEE)
pylab.show()


pylab.semilogy()
pylab.plot(l,numpy.abs(Theta['TTEE']))
pylab.plot(l,numpy.abs(Cov_TTEE))
pylab.show()


pylab.semilogy()
pylab.plot(l,numpy.abs(Theta['TTTE']))
pylab.plot(l,numpy.abs(Cov_TTTE))
pylab.show()


pylab.semilogy()
pylab.plot(l,numpy.abs(Theta['TEEE']))
pylab.plot(l,numpy.abs(Cov_TEEE))
pylab.show()



