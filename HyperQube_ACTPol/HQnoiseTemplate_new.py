#!/usr/bin/env python
from flipper import *
import pickle
p = flipperDict.flipperDict()
p.readFromFile(sys.argv[1])

fields=['T','Q','U']

def getNoiseMatrix(templateDir,array,season):
    
    def svd_pow(A, exponent):
        E, V = numpy.linalg.eigh(A)
        E[E<numpy.max(E,-1)[...,None]*1e-12] = 0
        return numpy.einsum("...ab,...b,...cb->...ac",V,E**exponent,V)

    
    print 'Template dir %s'%templateDir
    
    nn_TT = pickle.load(open('%s/noiseTemplate_TT_%s_%s.pkl'%(templateDir,array,season)))
    nn_TQ = pickle.load(open('%s/noiseTemplate_TQ_%s_%s.pkl'%(templateDir,array,season)))
    nn_TU = pickle.load(open('%s/noiseTemplate_TU_%s_%s.pkl'%(templateDir,array,season)))
    nn_QQ = pickle.load(open('%s/noiseTemplate_QQ_%s_%s.pkl'%(templateDir,array,season)))
    nn_QU = pickle.load(open('%s/noiseTemplate_QU_%s_%s.pkl'%(templateDir,array,season)))
    nn_UU = pickle.load(open('%s/noiseTemplate_UU_%s_%s.pkl'%(templateDir,array,season)))
    
    w = liteMap.liteMapFromFits('%s/weight_%s_%s_0'%(templateDir,array,season))
    noise = w.copy()
    noise.data[:] = 0.0
    random=numpy.random.randn(3,w.Ny,w.Nx)+1j*numpy.random.randn(3,w.Ny,w.Nx)
    
    cov=numpy.zeros((3,3,w.Ny,w.Nx))
    cov[0,0,:,:]=nn_TT.powerMap[:]
    cov[0,1,:,:]=nn_TQ.powerMap[:]
    cov[0,2,:,:]=nn_TU.powerMap[:]
    cov[1,0,:,:]=nn_TQ.powerMap[:]
    cov[1,1,:,:]=nn_QQ.powerMap[:]
    cov[1,2,:,:]=nn_QU.powerMap[:]
    cov[2,0,:,:]=nn_TU.powerMap[:]
    cov[2,1,:,:]=nn_QU.powerMap[:]
    cov[2,2,:,:]=nn_UU.powerMap[:]
    area = w.Nx*w.Ny*w.pixScaleX*w.pixScaleY
    cov=cov/area * (w.Nx*w.Ny)**2
    
    print "cov mat created"
    #chol = numpy.linalg.cholesky(cov.T)
    
    #    sim = numpy.einsum("xyab,byx->ayx",numpy.linalg.cholesky(cov.T),random)
    
    #  print sim[:,10,10]
    M=svd_pow(cov.T, 0.5)
    numpy.save('%s/svd_cov_%s_%s.npy'%(templateDir,array,season), M)



outDir='noiseTemplateAll'

try:
    os.mkdir(outDir)
except:
    pass



arrayTags = p['arrayTags']
seasonTags = p['seasonTags']



for array in arrayTags:
    for season in seasonTags:
        count1=0
        for l1 in fields:
            count1+=1
            count2=0
            for l2 in fields:
                count2+=1
                if count2<count1: continue
                
                m0 = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_0'%(l1,array,season))
                m0.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                m0_f = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_0'%(l2,array,season))
                m0_f.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                m1 = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_1'%(l1,array,season))
                m1.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                m1_f = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_1'%(l2,array,season))
                m1_f.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                w0 = liteMap.liteMapFromFits('patchesForNoise/weight_%s_%s_0'%(array,season))
                w1 = liteMap.liteMapFromFits('patchesForNoise/weight_%s_%s_1'%(array,season))
                
                m0.data[:] -= m1.data[:]
                m0_f.data[:] -= m1_f.data[:]
                
                m2 = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_2'%(l1,array,season))
                m2.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                m2_f = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_2'%(l2,array,season))
                m2_f.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                m3 = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_3'%(l1,array,season))
                m3.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                m3_f = liteMap.liteMapFromFits('patchesForNoise/%s_map_%s_%s_3'%(l2,array,season))
                m3_f.data[:] *= p['calibration_T_%s_%s'%(array,season)]
                
                w2 = liteMap.liteMapFromFits('patchesForNoise/weight_%s_%s_2'%(array,season))
                w3 = liteMap.liteMapFromFits('patchesForNoise/weight_%s_%s_3'%(array,season))
                
                m2.data[:] -= m3.data[:]
                m2_f.data[:] -= m3_f.data[:]
                
                m0.data[:] /= numpy.sqrt(1./w0.data[:]+1./w1.data[:])
                m2.data[:] /= numpy.sqrt(1./w2.data[:]+1./w3.data[:])
                m0_f.data[:] /= numpy.sqrt(1./w0.data[:]+1./w1.data[:])
                m2_f.data[:] /= numpy.sqrt(1./w2.data[:]+1./w3.data[:])
                
                nn0 = fftTools.powerFromLiteMap(m0,m0_f, applySlepianTaper=False)
                nn2 = fftTools.powerFromLiteMap(m2,m2_f, applySlepianTaper=False)
                
                nn = nn0.copy()
                nn.powerMap[:] += nn2.powerMap[:]
                nn.powerMap[:] /= 2.
                
                pickle.dump(nn,open('%s/noiseTemplate_%s%s_%s_%s.pkl'%(outDir,l1,l2,array,season),"w"))
                
                w0.writeFits('%s/weight_%s_%s_0'%(outDir,array,season),overWrite=True)
                w1.writeFits('%s/weight_%s_%s_1'%(outDir,array,season),overWrite=True)
                w2.writeFits('%s/weight_%s_%s_2'%(outDir,array,season),overWrite=True)
                w3.writeFits('%s/weight_%s_%s_3'%(outDir,array,season),overWrite=True)
                del nn

        getNoiseMatrix(outDir,array,season)