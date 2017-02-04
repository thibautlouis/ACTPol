from flipper import *
from flipperPol import *

if __name__=="__main__":
    global p
    p = flipperDict.flipperDict()
    p.read_from_file(sys.argv[1])

    patchDir = "patches"
    Ra0Array= p['Ra0Array']                                                                                                                    
    Ra1Array=  p['Ra1Array']                                                                                                                   
    Dec0Array =  p['Dec0Array']                                                                                                                
    Dec1Array =  p['Dec1Array']

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

    for i in xrange(nPatches):
        print i
        wt = liteMap.liteMapFromFits(patchDir+os.path.sep+'T_map_%03d_0'%(i)) 
        wt.data[:] = 0.
        for j in range(nDivs):

            wtmap = liteMap.liteMapFromFits(patchDir+'/weight_%03d_%d'%(i,j))
            wt.data[:] += wtmap.data[:]

        mBig=liteMap.liteMapFromFits(patchDir+os.path.sep+'Template_%03d'%i)
        wt=liteMapPol.padLiteMap(wt,mBig)
        print "Convolve"
        wt=wt.convolveWithGaussian(fwhm=20,nSigma=2.0)
        wt=wt.selectSubMap(Ra0Array[i],Ra1Array[i],Dec0Array[i],Dec1Array[i])
        wt.writeFits(patchDir+"/totalSmoothedWeightMap_%03d"%i,overWrite=True)
