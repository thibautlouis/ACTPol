#!/usr/bin/env python

from flipper import *
from flipperPol import *
from scipy.interpolate import splrep,splev
import scipy
import os
import pickle
import time
import speckMisc
import fft


def readMapAndCalibrate(patchDir,ar,season,i):
    
    T0 = liteMap.liteMapFromFits('%s/T_map_%s_%s_%d'%(patchDir,ar,season,i))
    Q0 = liteMap.liteMapFromFits('%s/Q_map_%s_%s_%d'%(patchDir,ar,season,i))
    U0 = liteMap.liteMapFromFits('%s/U_map_%s_%s_%d'%(patchDir,ar,season,i))
    
    T0.data[:]*=p['calibration_T_%s_%s'%(ar,season)]
    Q0.data[:]*=p['calibration_Pol_%s_%s'%(ar,season)]
    U0.data[:]*=p['calibration_Pol_%s_%s'%(ar,season)]
    
    return T0,Q0,U0


def weightedBinInAnnuli(p2d,weightMap,binningFile,trimAtL,powerOfL):
    binLo,binHi,binCent = fftTools.readBinningFile(binningFile)
    id = numpy.where(binHi<trimAtL)
    binHi = binHi[id]
    binLo = binLo[id]
    binCent = binCent[id]

    binnedPower = binCent.copy()*0.
    binCount = binCent.copy()*0.
    weightedBincount = binCent.copy()
    modIntLMap = numpy.array(p2d.modLMap + 0.5,dtype='int64')
    for ibin in xrange(len(binHi)):
        loc = numpy.where((modIntLMap >= binLo[ibin]) & (modIntLMap <= binHi[ibin]))
        binMap = p2d.powerMap.copy()*0.
        binMap[loc] = weightMap[loc]
        binnedPower[ibin] = numpy.sum(p2d.powerMap*binMap*p2d.modLMap**powerOfL)/numpy.sum(binMap)
        binCount[ibin] = len(loc[0])
        weightedBincount[ibin] = 1./(numpy.sum(weightMap[loc]**2)/(numpy.sum(weightMap[loc]))**2)
        #print binCount[ibin]/weightedBincount
    return binLo,binHi,binCent,binnedPower, weightedBincount/2.


def makeMbbInvArray(mbb,mbb_cos2,mbb_sin2,mbb_cos1,mbb_sin1,mbb_sincos):
    
    N=mbb.shape[0]
    mbbAll= numpy.zeros((6*N,6*N))
        
    mbbAll[:N,:N]=mbb

    mbbAll[N:2*N,N:2*N]=mbb_cos1
    mbbAll[2*N:3*N,2*N:3*N]=mbb_cos1
    mbbAll[N:2*N,2*N:3*N]=-mbb_sin1
    mbbAll[2*N:3*N,N:2*N]=mbb_sin1
    
    mbbAll[3*N:4*N,3*N:4*N]=mbb_cos2
    mbbAll[3*N:4*N,4*N:5*N]=-2*mbb_sincos
    mbbAll[3*N:4*N,5*N:6*N]=mbb_sin2
    
    mbbAll[4*N:5*N,3*N:4*N]=mbb_sincos
    mbbAll[4*N:5*N,4*N:5*N]=mbb_cos2-mbb_sin2
    mbbAll[4*N:5*N,5*N:6*N]=-mbb_sincos

    mbbAll[5*N:6*N,3*N:4*N]=mbb_sin2
    mbbAll[5*N:6*N,4*N:5*N]=2*mbb_sincos
    mbbAll[5*N:6*N,5*N:6*N]=mbb_cos2

    mbbInv_All= numpy.linalg.inv(mbbAll)
    
    return(mbbInv_All)

def ProcessAndwriteSpectra(cl_vec,filterArray,name,fields,ar1,ar2,spTag,binnedBeamDict=None,iar=None):
    
    countAll=0
    count1=0
    for l1 in fields:
        count1+=1
        count2=0
        for l2 in fields:
            count2+=1
            if count2<count1: continue
            
            # remove filter
            cl=cl_vec[countAll*Nbin:(countAll+1)*Nbin]*filterArray[l1+l2]**2
            
            # There is an additional correction for the autos as MCM had a transfer
            # function B_l_AR1*B_l_AR_2
            if iar !=None:
                cl*= binnedBeamDict[sps[0]][iar-1]/binnedBeamDict[sps[0]][iar]
            
            
            gName = '%s/%s_%s%s_%sx%s_%s.dat'%(specDir,name,l1,l2,ar1,ar2,spTag)
            
            speckMisc.writeBinnedSpectrum(lBin,cl,binCount,gName)
            
            countAll+=1







if __name__ == "__main__":
    t=time.time()
    p = flipperDict.flipperDict()
    p.read_from_file(sys.argv[1])

    specDir = 'spectra/'
    spec2dDir = 'spectra2d/'
    patchDir = 'patches'
    if len(sys.argv)> 2:
        specDir = 'spectra_%s/'%sys.argv[2]
        spec2dDir = 'spectra2d_%s/'%sys.argv[2]
        patchDir = p['mapDir']+'patches_%s'%sys.argv[2]
        
    cosineApod = p['cosineApodization']
    mask = p['mask']
    Write2dSpectrum = p['Write2dSpectrum']
    binningFile=p['binningFile']
    healpixConv = p['healpixConv']
    pwDict = p['prewhitener']

    arrays = p['arrayTags']
    seasons = p['seasonTags']

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


    print "Found %d patches with %d sub-season divisions in each"%(nPatches, nDivs)
    
    trimAtL = p['trimAtL']
    
    try:
        os.makedirs(specDir)
    except:
        pass
    
    if  Write2dSpectrum['apply']==True:
        try:
            os.makedirs(spec2dDir)
        except:
            pass

    hpfDict = p['highPassCosSqFilter']
    filter_Pol = 1.0
    filter_T= 1.0

    if hpfDict['apply']:
        print "Will take off the cos^2 high pass filter"
        filter_T = speckMisc.getBinnedInvCosSqFilter(hpfDict['lMin'],hpfDict['lMax'],p['binningFile'],trimAtL)
    
    filter_Cross=numpy.sqrt(filter_T*filter_Pol)                                                                                                        

    fields=['T','E','B']

    filterArray={}
    filterArray['TT']=filter_T
    filterArray['EE'],filterArray['BB'],filterArray['EB'],filterArray['BE']=filter_Pol,filter_Pol,filter_Pol,filter_Pol
    filterArray['TE'],filterArray['TB'],filterArray['ET'],filterArray['BT']=filter_Cross,filter_Cross,filter_Cross,filter_Cross

    arraysPair = [arrays[0],arrays[-1]]
    seasonPairs = []

    for i in xrange(len(seasons)):
        for j in xrange(len(seasons)):
            if arrays[0] == arrays[-1] and i>j: continue
            seasonPairs += [[seasons[i],seasons[j]]]

    print arraysPair
    print seasonPairs

    binnedBeamDict = {}

    for seas in seasons:
        binnedBeamWindow = []
        for ar in arraysPair:
            Bb = speckMisc.getBinnedBeamTransfer(p['beam_%s_%s'%(ar,seas)],p['binningFile'],trimAtL)
            binnedBeamWindow += [Bb]
        binnedBeamDict[seas] = binnedBeamWindow

    mbbInvArray={}

    for sps in seasonPairs:
        spTag = '%sx%s'%(sps[0],sps[1])
        arrayTag = '%sx%s'%(arraysPair[0],arraysPair[1])
        
        mbb = numpy.loadtxt('mcm/'+p['mcmFileRoot']+'_%s_%s.dat'%(arrayTag,spTag))
        mbb_cos1 = numpy.loadtxt('mcm/'+p['mcmFileRoot']+'_cos1_%s_%s.dat'%(arrayTag,spTag))
        mbb_sin1 = numpy.loadtxt('mcm/'+p['mcmFileRoot']+'_sin1_%s_%s.dat'%(arrayTag,spTag))
        mbb_cos2 =numpy.loadtxt('mcm/'+p['mcmFileRoot']+'_cos2_%s_%s.dat'%(arrayTag,spTag))
        mbb_sin2 = numpy.loadtxt('mcm/'+p['mcmFileRoot']+'_sin2_%s_%s.dat'%(arrayTag,spTag))
        mbb_sincos = numpy.loadtxt('mcm/'+p['mcmFileRoot']+'_sincos_%s_%s.dat'%(arrayTag,spTag))

        mbbInvArray[spTag]=makeMbbInvArray(mbb,mbb_cos2,mbb_sin2,mbb_cos1,mbb_sin1,mbb_sincos)


    for sps in seasonPairs:
        
            spTag = '%sx%s'%(sps[0],sps[1])
            arrayTag = '%sx%s'%(arraysPair[0],arraysPair[1])
            
            p2dWeight_TT = numpy.load('noiseAndWeights/weightMap_%s_%s_TT.npy'%(arrayTag,spTag))
            p2dWeight_Cross = numpy.load('noiseAndWeights/weightMap_%s_%s_Cross.npy'%(arrayTag,spTag))
            p2dWeight_Pol = numpy.load('noiseAndWeights/weightMap_%s_%s_Pol.npy'%(arrayTag,spTag))
            
            p2dWeightArray={}

            p2dWeightArray['TT'],p2dWeightArray['EE'],p2dWeightArray['TE'],p2dWeightArray['ET']=p2dWeight_TT,p2dWeight_Pol,p2dWeight_Cross,p2dWeight_Cross
            p2dWeightArray['TB'],p2dWeightArray['EB'],p2dWeightArray['BT'],p2dWeightArray['BE']=p2dWeight_Cross,p2dWeight_Pol,p2dWeight_Cross,p2dWeight_Pol
            p2dWeightArray['BB']=p2dWeight_Pol

            window_T = liteMap.liteMapFromFits('auxMaps/finalWindow_%s_%s_T.fits'%(arrayTag,spTag))
            window_Pol = liteMap.liteMapFromFits('auxMaps/finalWindow_%s_%s_Pol.fits'%(arrayTag,spTag))
            
            fftTemp=window_T.copy()
            fft.fft(fftTemp.data,axes=[-2,-1],flags=['FFTW_MEASURE'])
            del fftTemp
    
            modLMap,angLMap=fftPol.makeEllandAngCoordinate(window_T,bufferFactor=1)
        
            p2dTemp=fftTools.powerFromLiteMap(window_T)
            p2dTemp.powerMap[:]=0
            p2dTemp=p2dTemp.trimAtL(trimAtL+500)
            
            clAutoPatch = {}
            clCrossPatch = {}
            Sump2dPatch = {}
        
  


            # Get the auto array first
            iar = 0
            for ar in arrays:
                if sps[0] == sps[1]:
                    
                    for l1 in fields:
                        for l2 in fields:
                            clAutoPatch[l1+l2] = []
                
                    for i in xrange(nDivs):
                        print 'computing spectrum %sx%s %s: %d%d '%(ar,ar,spTag, i,i)
                        
                        T0,Q0,U0=readMapAndCalibrate(patchDir,ar,sps[0],i)
                        
                        fT0,fE0,fB0=fftPol.TQUtoPureTEB(T0,Q0,U0,window_T,window_Pol, modLMap,angLMap,method='standard',fftType='fftw3')
                        
                        area = T0.Nx*T0.Ny*T0.pixScaleX*T0.pixScaleY
                        pA={}
                            
                        pA['TT'],pA['TE'],pA['ET'],pA['TB'],pA['BT'],pA['EE'],pA['EB'],pA['BE'],pA['BB'] = fftPol.fourierTEBtoPowerTEB(fT0,fE0,fB0,fT0,fE0,fB0,trimAtL+500)
                            
                        for l1 in fields:
                            for l2 in fields:
                                
                                lLower,lUpper,lBin,clBin,binCount = weightedBinInAnnuli(pA[l1+l2],p2dWeightArray[l1+l2],binningFile,trimAtL,powerOfL=0)
                                clbinDecoup = clBin*area
                                
                                fName = '%s/clBin_%s%s_%sx%s_%s_%d%d.dat'%(specDir,l1,l2,ar,ar,spTag,i,i)
                                            
  #                              speckMisc.writeBinnedSpectrum(lBin,clbinDecoup,binCount,fName)
                                    
                                clAutoPatch[l1+l2] += [clbinDecoup]
            

                    clAuto=[]
                    
                    count1=0
                    for l1 in fields:
                        count1+=1
                        count2=0
                        for l2 in fields:
                            count2+=1
    
                            if count2<count1: continue
                            print l1,l2
        
                            clAutoPatchMean = numpy.mean(clAutoPatch[l1+l2],axis=0)
                            
                            fName = '%s/window_clBinAutoMean_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,ar,ar,spTag)
                            
                            speckMisc.writeBinnedSpectrum(lBin,clAutoPatchMean,binCount,fName)
                            clAuto= numpy.append(clAuto,clAutoPatchMean)
                
                    clAuto_deconv= numpy.dot(mbbInvArray[spTag],clAuto)

                    Nbin=len(clAutoPatchMean)
                    name='clBinAutoMean'
                    
                    ProcessAndwriteSpectra(clAuto_deconv,filterArray,name,fields,ar,ar,spTag,binnedBeamDict=binnedBeamDict,iar=iar)
                    
                    iar += 1


            # Now do the cross-array spectra
            for l1 in fields:
                for l2 in fields:
                    clAutoPatch[l1+l2] = []   #cross-array auto spec
                    clCrossPatch[l1+l2] = []
                    Sump2dPatch[l1+l2] = []
   
            for i in xrange(nDivs):
                
                T0,Q0,U0=readMapAndCalibrate(patchDir,arrays[0],sps[0],i)
                
                fT0,fE0,fB0=fftPol.TQUtoPureTEB(T0,Q0,U0,window_T,window_Pol, modLMap,angLMap,method='standard',fftType='fftw3')
                
                del T0,Q0,U0

                for j in xrange(nDivs):
                    if arrays[0]==arrays[-1] and sps[0] == sps[1] and i<j: continue
                    
                    print 'computing spectrum %sx%s %s: %d%d '%(arrays[0],arrays[-1],spTag, i,j)

                    T1,Q1,U1=readMapAndCalibrate(patchDir,arrays[-1],sps[1],j)
                    fT1,fE1,fB1=fftPol.TQUtoPureTEB(T1,Q1,U1,window_T,window_Pol, modLMap,angLMap,method='standard',fftType='fftw3')

                    pA['TT'],pA['TE'],pA['ET'],pA['TB'],pA['BT'],pA['EE'],pA['EB'],pA['BE'],pA['BB'] = fftPol.fourierTEBtoPowerTEB(fT0,fE0,fB0,fT1,fE1,fB1,trimAtL+500)
                        
                    for l1 in fields:
                        for l2 in fields:


                            lLower,lUpper,lBin,clBin,binCount = weightedBinInAnnuli(pA[l1+l2],p2dWeightArray[l1+l2],binningFile,trimAtL,powerOfL=0)
                            clbinDecoup = clBin*area

                            fName = '%s/clBin_%s%s_%sx%s_%s_%d%d.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag,i,j)
             

#                            speckMisc.writeBinnedSpectrum(lBin,clbinDecoup,binCount,fName)


                            if i == j  and sps[0] == sps[1]:
                                print 'auto', l1,l2,i,j,sps[0],sps[1]
                                clAutoPatch[l1+l2] += [clbinDecoup]
                            else:
                                print 'cross', l1,l2,i,j,sps[0],sps[1]
                                clCrossPatch[l1+l2] += [clbinDecoup]
                                Sump2dPatch[l1+l2]+=[pA[l1+l2].powerMap[:]]
                                    
                                    
            clCrossPatch['TE']+=clCrossPatch['ET']
            Sump2dPatch['TE']+=Sump2dPatch['ET']


            clCrossPatch['TB']+=clCrossPatch['BT']
            Sump2dPatch['TB']+=Sump2dPatch['BT']
                                        
            clCrossPatch['EB']+=clCrossPatch['BE']
            Sump2dPatch['EB']+=Sump2dPatch['BE']
                                        
            if sps[0] == sps[1] :
                clAutoPatch['TE']+=clAutoPatch['ET']
                clAutoPatch['TB']+=clAutoPatch['BT']
                clAutoPatch['EB']+=clAutoPatch['BE']


            clCross=[]
            clAuto=[]
            count1=0
            for l1 in fields:
                count1+=1
                count2=0
                for l2 in fields:
                    count2+=1
        
                    if count2<count1: continue
                    print l1,l2,'mean done on %d cross spectra'%(len(clCrossPatch[l1+l2]))
                    clCrossPatchMean = numpy.mean(clCrossPatch[l1+l2],axis=0)
        
                    fName = '%s/window_clBinCrossMean_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag)
                    speckMisc.writeBinnedSpectrum(lBin,clCrossPatchMean,binCount,fName)
                    clCross= numpy.append(clCross,clCrossPatchMean)
        
        
                    if sps[0] == sps[1] :
                        print l1,l2,'mean done on %d auto spectra'%(len(clAutoPatch[l1+l2]))
                        clAutoPatchMean = numpy.mean(clAutoPatch[l1+l2],axis=0)
                        fName = '%s/window_clBinAutoMean_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag)
                        speckMisc.writeBinnedSpectrum(lBin,clAutoPatchMean,binCount,fName)
                        clAuto= numpy.append(clAuto,clAutoPatchMean)

                    Sump2dPatchMean=   numpy.mean(Sump2dPatch[l1+l2],axis=0)
        
                    if Write2dSpectrum['apply']:
                        p2dTemp.powerMap[:]= Sump2dPatchMean
                        pickle.dump(p2dTemp,open('%s/p2d_%s%s_%sx%s_%s.pkl'%(spec2dDir,l1,l2,arrays[0],arrays[-1],spTag),mode="w"))


            clCross_deconv= numpy.dot(mbbInvArray[spTag],clCross)
    
            if sps[0] == sps[1] :
                clAuto_deconv= numpy.dot(mbbInvArray[spTag],clAuto)
                name='clBinAutoMean'
                ProcessAndwriteSpectra(clAuto_deconv,filterArray,name,fields,arrays[0],arrays[-1],spTag)
    
            Nbin=len(clCrossPatchMean)
    
            name='clBinCrossMean'
    
            ProcessAndwriteSpectra(clCross_deconv,filterArray,name,fields,arrays[0],arrays[-1],spTag)
            print time.time()-t
