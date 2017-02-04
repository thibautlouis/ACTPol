#!/usr/bin/env python
from flipper import *
import speckMisc
import pickle
import liteMapPol

def apoCorrection(patchDir,sps,arrays,spTag):
    
    print 'ok routine a bit out to date as we now have two window functions'
    print 'also be careful if you wanna pad'


    
    m = liteMap.liteMapFromFits('%s/T_map_%s_%s_0'%(patchDir,arrays[0],sps))
    
    area=numpy.float(m.Nx*m.Ny)
    area/=(m.Nx*m.Ny)
    
    m.data[:] = 1.0
    
    c=liteMapPol.initializeCosineWindow(m,cA['lenApod'],cA['pad'])
    
    m=liteMap.liteMapFromFits('auxMaps/finalWindow_%sx%s_%s_T.fits'%(arrays[0],arrays[-1],spTag))
    
    m.data/=numpy.max(m.data)
    
    factor= numpy.mean(m.data**4)/(numpy.mean(m.data**2))**2/numpy.mean(c.data[:])
    facs = [factor]
    
    print spTag,arrays[0],arrays[-1], numpy.sqrt(factor)
    
    return numpy.sqrt(facs)



p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])

theoryFile = p['theoryFile']
seasons = p['seasonTags']
arrays = p['arrayTags']

specDir = 'spectra/'
patchDir = 'patches/'

if len(sys.argv)> 2:
    specDir = 'spectra_%s/'%sys.argv[2]
    patchDir = p['mapDir']+'patches_%s/'%sys.argv[2]


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


cA = p['cosineApodization']
apoCorr =1.
fac = 2.0

apoCorr={}

for sps in seasonPairs:
    
    spTag = '%sx%s'%(sps[0],sps[1])
    
    apoCorr[spTag] = apoCorrection(patchDir,sps[0],arrays,spTag)


meanCrossSpec={}
binWeight={}

meanAutoSpec_A={}
meanAutoSpec_B={}
meanAutoSpec_AB={}

meanNoise_A={}
meanNoise_B={}
meanNoise_AB={}

Theta={}

fields=['T','E','B']



count1=0
for l1 in fields:
    count1+=1
    count2=0
    
    for l2 in fields:
        count2+=1
        if count2<count1: continue
        
        for sps in seasonPairs:
            spTag = '%sx%s'%(sps[0],sps[1])
        
            lbin,meanCrossSpec[l1+l2,spTag],binWeight[l1+l2,spTag] = numpy.loadtxt('%s/clBinCrossMean_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag),unpack=True)
            meanCrossSpec[l2+l1,spTag]=meanCrossSpec[l1+l2,spTag]
            
            #A seasoni AR1 x seasoni AR1
            lbin,meanAutoSpec_A[l1+l2,spTag],binWeight[l1+l2,spTag] = numpy.loadtxt('%s/clBinAutoMean_%s%s_%sx%s_%sx%s.dat'%(specDir,l1,l2,arrays[0],arrays[0],sps[0],sps[0]),unpack=True)
            
            #B seasonj AR2 x seasonj AR2
            lbin,meanAutoSpec_B[l1+l2,spTag],binWeight[l1+l2,spTag] = numpy.loadtxt('%s/clBinAutoMean_%s%s_%sx%s_%sx%s.dat'%(specDir,l1,l2,arrays[-1],arrays[-1],sps[1],sps[1]),unpack=True)
            
            #AB is only used for seasoni AR1 x seasoni AR2
            lbin,meanAutoSpec_AB[l1+l2,spTag],binWeight[l1+l2,spTag] = numpy.loadtxt('%s/clBinAutoMean_%s%s_%sx%s_%sx%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],sps[0],sps[0]),unpack=True)

            meanAutoSpec_A[l2+l1,spTag]=meanAutoSpec_A[l1+l2,spTag]
            meanAutoSpec_B[l2+l1,spTag]=meanAutoSpec_B[l1+l2,spTag]
            meanAutoSpec_AB[l2+l1,spTag]=meanAutoSpec_AB[l1+l2,spTag]


            meanNoise_A[l1+l2,spTag]=meanAutoSpec_A[l1+l2,spTag]-meanCrossSpec[l1+l2,spTag]
            meanNoise_A[l2+l1,spTag]=meanNoise_A[l1+l2,spTag]

            meanNoise_B[l1+l2,spTag]=meanAutoSpec_B[l1+l2,spTag]-meanCrossSpec[l1+l2,spTag]
            meanNoise_B[l2+l1,spTag]=meanNoise_B[l1+l2,spTag]

            meanNoise_AB[l1+l2,spTag]=meanAutoSpec_AB[l1+l2,spTag]-meanCrossSpec[l1+l2,spTag]
            meanNoise_AB[l2+l1,spTag]=meanNoise_AB[l1+l2,spTag]


            fName = '%s/noise_%s%s_%sx%s_%sx%s.dat'%(specDir,l1,l2,arrays[0],arrays[0],sps[0],sps[0])
            speckMisc.writeBinnedSpectrum(lbin,meanNoise_A[l1+l2,spTag]/nDivs,binWeight[l1+l2,spTag],fName)

            fName = '%s/noise_%s%s_%sx%s_%sx%s.dat'%(specDir,l1,l2,arrays[-1],arrays[-1],sps[1],sps[1])
            speckMisc.writeBinnedSpectrum(lbin,meanNoise_B[l1+l2,spTag]/nDivs,binWeight[l1+l2,spTag],fName)


            fName = '%s/noise_%s%s_%sx%s_%sx%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],sps[0],sps[1])
            speckMisc.writeBinnedSpectrum(lbin,meanNoise_AB[l1+l2,spTag]/nDivs,binWeight[l1+l2,spTag],fName)




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
                
                for sps in seasonPairs:
                    
                    spTag = '%sx%s'%(sps[0],sps[1])
                    
                    if sps[0]==sps[1] and arrays[0] == arrays[-1]:
                        #theta^{\alpha A \alpha A}=1/nu_b(Cl_WY Cl_XZ+ Cl_WZ Cl_XY)+1/(N_S nu_b)[Cl_WY Nl^{\alpha A}_XZ+ Cl_XZ Nl^{\alpha A}_WY+ Cl_WZ  Nl^{\alpha A}_XY+ Cl_XY  Nl^{\alpha A}_WZ]
                        #+1/(N_S(N_S-1) nu_b)[Nl^{\alpha A}_WY Nl^{\alpha A}_XZ+ Nl^{\alpha A}_WZ Nl^{\alpha A}_XY
                        
                        Theta[l1+l2+l3+l4] =meanCrossSpec[l1+l3,spTag]*meanCrossSpec[l2+l4,spTag]+meanCrossSpec[l1+l4,spTag]*meanCrossSpec[l2+l3,spTag]
                        Theta[l1+l2+l3+l4] += 1./nDivs*(meanCrossSpec[l1+l3,spTag]*meanNoise_A[l2+l4,spTag]+meanCrossSpec[l2+l4,spTag]*meanNoise_A[l1+l3,spTag]+meanCrossSpec[l1+l4,spTag]*meanNoise_A[l2+l3,spTag]+meanCrossSpec[l2+l3,spTag]*meanNoise_A[l1+l4,spTag])
                        Theta[l1+l2+l3+l4] += 1./(nDivs*(nDivs-1))*(meanNoise_A[l1+l3,spTag]*meanNoise_A[l2+l4,spTag]+meanNoise_A[l1+l4,spTag]*meanNoise_A[l2+l3,spTag])
                    
                    elif sps[0]==sps[1] and arrays[0] != arrays[-1]:
                        #theta^{\alpha A \alpha B}=1/nu_b(Cl_WY Cl_XZ+ Cl_WZ Cl_XY)+1/(N_S nu_b)[Cl_WY Nl^{\alpha B}_XZ+ Cl_XZ Nl^{\alpha A}_WY+ Cl_WZ  Nl^{\alpha\alpha AB}_XY+ Cl_XY  Nl^{\alpha\alpha AB}_WZ]
                        #+1/(N_S(N_S-1) nu_b)[Nl^{\alpha A}_WY Nl^{\alpha B}_XZ+ Nl^{\alpha AB}_WZ Nl^{\alpha AB}_XY


                        Theta[l1+l2+l3+l4] =meanCrossSpec[l1+l3,spTag]*meanCrossSpec[l2+l4,spTag]+meanCrossSpec[l1+l4,spTag]*meanCrossSpec[l2+l3,spTag]
                        Theta[l1+l2+l3+l4] += 1./nDivs*(meanCrossSpec[l1+l3,spTag]*meanNoise_B[l2+l4,spTag]+meanCrossSpec[l2+l4,spTag]*meanNoise_A[l1+l3,spTag]+meanCrossSpec[l1+l4,spTag]*meanNoise_AB[l2+l3,spTag]+meanCrossSpec[l2+l3,spTag]*meanNoise_AB[l1+l4,spTag])
                        Theta[l1+l2+l3+l4] += 1./(nDivs*(nDivs-1))*(meanNoise_A[l1+l3,spTag]*meanNoise_B[l2+l4,spTag]+meanNoise_AB[l1+l4,spTag]*meanNoise_AB[l2+l3,spTag])
                            
                    else:
                        #theta^{\alpha A \beta B}=1/nu_b(Cl_WY Cl_XZ+ Cl_WZ Cl_XY)+1/(N_S nu_b)[Cl_WY Nl^{\beta B}_XZ+Cl_XZ Nl^{\alpha A}_WY]+1/N_S**2[Nl^{\alpha A}_WY Nl^{\beta B}_XZ]
                        
                        Theta[l1+l2+l3+l4] =meanCrossSpec[l1+l3,spTag]*meanCrossSpec[l2+l4,spTag]+meanCrossSpec[l1+l4,spTag]*meanCrossSpec[l2+l3,spTag]
                        Theta[l1+l2+l3+l4] +=1./nDivs*(meanCrossSpec[l1+l3,spTag]*meanNoise_B[l2+l4,spTag]+meanCrossSpec[l2+l4,spTag]*meanNoise_A[l1+l3,spTag])
                        Theta[l1+l2+l3+l4] +=1./(nDivs**2)*(meanNoise_A[l1+l3,spTag]*meanNoise_B[l2+l4,spTag])

                    Theta[l1+l2+l3+l4] /= (fac*binWeight['BB',spTag])
                    Theta[l1+l2+l3+l4] *=apoCorr[spTag]**2
                
                
                

                    if count1==count3 and count2==count4:
                        
                        fName = '%s/spectrum_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag)
                        speckMisc.writeBinnedSpectrum(lbin,meanCrossSpec[l1+l2,spTag],numpy.sqrt(Theta[l1+l2+l1+l2]),fName)





