#!/usr/bin/env python
from flipper import *
import pickle
p = flipperDict.flipperDict()
p.readFromFile(sys.argv[1])

fields=['T','Q','U']

outDir='noiseTemplate'

try:
    os.mkdir(outDir)
except:
    pass


arrayTags = p['arrayTags']
seasonTags = p['seasonTags']



for array in arrayTags:
    for season in seasonTags:
        for l1 in fields:

            m0 = liteMap.liteMapFromFits('patches/%s_map_%s_%s_0'%(l1,array,season))
            m0.data[:] *= p['calibration_T_%s_%s'%(array,season)]
        
            m1 = liteMap.liteMapFromFits('patches/%s_map_%s_%s_1'%(l1,array,season))
            m1.data[:] *= p['calibration_T_%s_%s'%(array,season)]

            w0 = liteMap.liteMapFromFits('patches/weight_%s_%s_0'%(array,season))
            w1 = liteMap.liteMapFromFits('patches/weight_%s_%s_1'%(array,season))
        
            m0.data[:] -= m1.data[:]
        
            m2 = liteMap.liteMapFromFits('patches/%s_map_%s_%s_2'%(l1,array,season))
            m2.data[:] *= p['calibration_T_%s_%s'%(array,season)]
    
            m3 = liteMap.liteMapFromFits('patches/%s_map_%s_%s_3'%(l1,array,season))
            m3.data[:] *= p['calibration_T_%s_%s'%(array,season)]
    
            w2 = liteMap.liteMapFromFits('patches/weight_%s_%s_2'%(array,season))
            w3 = liteMap.liteMapFromFits('patches/weight_%s_%s_3'%(array,season))
    
            m2.data[:] -= m3.data[:]
        
            m0.data[:] /= numpy.sqrt(1./w0.data[:]+1./w1.data[:])
            m2.data[:] /= numpy.sqrt(1./w2.data[:]+1./w3.data[:])
        
            nn0 = fftTools.powerFromLiteMap(m0, applySlepianTaper=False)
            nn2 = fftTools.powerFromLiteMap(m2, applySlepianTaper=False)
        
            nn = nn0.copy()
            nn.powerMap[:] += nn2.powerMap[:]
            nn.powerMap[:] /= 2.
        
            pickle.dump(nn,open('%s/noiseTemplate_%s_%s_%s.pkl'%(outDir,l1,array,season),"w"))
    
            w0.writeFits('%s/weight_%s_%s_0'%(outDir,array,season),overWrite=True)
            w1.writeFits('%s/weight_%s_%s_1'%(outDir,array,season),overWrite=True)
            w2.writeFits('%s/weight_%s_%s_2'%(outDir,array,season),overWrite=True)
            w3.writeFits('%s/weight_%s_%s_3'%(outDir,array,season),overWrite=True)
        
        
