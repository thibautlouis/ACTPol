#!/usr/bin/env python
from flipper import *
from flipperPol import *
from mpi4py import MPI
import pickle


def cosineSqFilter(map,lMin,lMax,vkMaskLimits=None):
    filteredMap = map.copy()
    ft = fftTools.fftFromLiteMap(map)
    ell = ft.modLMap
    idSub = numpy.where((ell> lMin) & (ell <lMax))
    idLow = numpy.where(ell<lMin)  
    filter = (ell*0.+1.0)
    filter[idLow] = 0.
    filter[idSub] *= (numpy.cos((lMax - ell[idSub])/(lMax-lMin)*numpy.pi/2.))**2
    ft.kMap[:] *= filter[:]
    if vkMaskLimits != None:
        #Yank the k-mode
        idvk = numpy.where((ft.lx >vkMaskLimits[0]) & (ft.lx<vkMaskLimits[1]))
        ft.kMap[:,idvk] = 0.
    
    filteredMap.data[:] = ft.mapFromFFT()
    return filteredMap


print "Reading dict file"
p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])


buffer=p['buffer']
filter = p['highPassCosSqFilter']
pwDict = p['prewhitener']
pixWeight=p['pixWeight']
beam1d=p['beam1d']
mask = p['mask']

arrayTags = p['arrayTags']
seasonTags = p['seasonTags']


ra0= p['ra0']
ra1=  p['ra1']
dec0 =  p['dec0']
dec1 =  p['dec1']


patchDir='patches'
    
try:
    os.makedirs(patchDir)
except:
    pass
    
try:
    os.makedirs(patchDir+'/mapPlots')
except:
    pass


fields=['T','Q','U']

if p['applyMask']==True:
    mask=liteMap.liteMapFromFits(p['mask'])
    mask= mask.selectSubMap(ra0-buffer,ra1+buffer,dec0-buffer,dec1+buffer)
    #The double selection here is to fix 1/2 pixels nightmare.
    mask= mask.selectSubMap(ra0,ra1,dec0,dec1)
    mask.writeFits(patchDir+os.path.sep+'mask',overWrite=True)

    if p['maskPol']!=None:
        maskPol=liteMap.liteMapFromFits(p['maskPol'])
        maskPol= maskPol.selectSubMap(ra0-buffer,ra1+buffer,dec0-buffer,dec1+buffer)
        #The double selection here is to fix 1/2 pixels nightmare.
        maskPol= maskPol.selectSubMap(ra0,ra1,dec0,dec1)
        maskPol.writeFits(patchDir+os.path.sep+'maskPol',overWrite=True)


for array in arrayTags:
    
    for season in seasonTags:
        
        for l1 in fields:

            mapFiles = p['mapFiles_%s_%s_%s'%(array,season,l1)]
    
            count=0
            
            for mapFile in mapFiles:
        
                print mapFile
    
                m0 = liteMap.liteMapFromFits(mapFile)
                m= m0.selectSubMap(ra0-buffer,ra1+buffer,dec0-buffer,dec1+buffer)

                if l1=='T' and filter['apply']:
                    print "filtering subMap..."
                    m = cosineSqFilter(m,filter['lMin'],filter['lMax'])

            
                
                if l1=='T' and pwDict['apply']:
                    print "prewhitening subMap ..."
                    pw = prewhitener.prewhitener(pwDict['radius'],addBackFraction=pwDict['addBackFraction'],smoothingFWHM=pwDict['gaussFWHM'],map = m)
                    m = pw.apply(m)

                    
                
                if p['applyPixelWeights']:
                    weightMap = liteMap.liteMapFromFits(p['weightFiles_%s_%s'%(array,season)][count])
                    weightSubMap = weightMap.selectSubMap(ra0-buffer,ra1+buffer,dec0-buffer,dec1+buffer)
                    weightSubMap = weightSubMap.selectSubMap(ra0,ra1,dec0,dec1)

                    weightSubMap.writeFits(patchDir+os.path.sep+'weight_%s_%s_%d'%(array,season,count),overWrite=True)


                if buffer !=0:
                    m=m.selectSubMap(ra0,ra1,dec0,dec1)


                m.writeFits(patchDir+os.path.sep+'%s_map_%s_%s_%d'%(l1,array,season,count),overWrite=True)
            
                        
                if l1 != 'T':
                    m = cosineSqFilter(m,filter['lMin'],filter['lMax'])

                
                count+=1


				

