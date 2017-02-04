#!/usr/bin/env python

from flipper import *
from flipperPol import *
import HQmodeCouplingPol as mcmPol
import pickle
import scipy

def makePrewithenerTransfer(window,trimAtL,alpha=0.02):
    disk1=numpy.loadtxt('diskKern1.dat')
    disk2=numpy.loadtxt('diskKern3.dat')
    disk=window.copy()
    disk.data[:]=0
    ny1,nx1,ny2,nx2=disk1.shape[0],disk1.shape[1],disk2.shape[0],disk2.shape[1]
    Ny,Nx=disk.Ny,disk.Nx
    disk.data[Ny/2.-ny1/2.:Ny/2.+ny1/2.,Nx/2.-nx1/2.:Nx/2.+nx1/2.]=disk1
    disk.data[Ny/2.-ny2/2.:Ny/2.+ny2/2.,Nx/2.-nx2/2.:Nx/2.+nx2/2.]-=disk2
    transfer_pw=fftTools.powerFromLiteMap(disk)
    transfer_pw=transfer_pw.trimAtL(trimAtL)
    transfer_pw.powerMap[:]*=Ny*Nx/(disk.pixScaleX*disk.pixScaleY)
    transfer_pw.powerMap[:]+=alpha**2+2*alpha*numpy.sqrt(transfer_pw.powerMap[:])
    return(transfer_pw)





if __name__=="__main__":
    global p
    p = flipperDict.flipperDict()
    p.read_from_file(sys.argv[1])

    patchDir = "patches"
    if len(sys.argv)> 2:
        patchDir = p['mapDir']+'/patches_000'
    
    
    cosineApod = p['cosineApodization']
    mask = p['mask']
    pwDict = p['prewhitener']
    beam1d=p['beam1d']
    nfork=p['nfork']
    buffer=p['buffer']
    doPadding = p['doPadding']
    
    useSmoothedWeight= p['useSmoothedWeight']
    
    arrays = p['arrayTags']
    seasons = p['seasonTags']

    try:
        os.makedirs("noiseAndWeights")
    except:
        pass
    
    try:
        os.makedirs("mcm")
    except:
        pass
    
    if not(os.path.exists('auxMaps')): os.mkdir('auxMaps')
    
    
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

    binLower,binUpper,binCenter = fftTools.readBinningFile(p['binningFile'])
    id = numpy.where(binUpper <p['trimAtL'])

    # if apply taper, read in a patch and create the order zero taper taper

    transfers = {}
    pixWeightsAll = {}

    arraysPair = [arrays[0],arrays[-1]]
    seasonPairs = []

    for i in xrange(len(seasons)):
        for j in xrange(len(seasons)):
            if arrays[0]==arrays[-1] and i>j: continue
            seasonPairs += [[seasons[i],seasons[j]]]

    print arraysPair
    print seasonPairs


    for sps in seasonPairs:
        
        spTag = '%sx%s'%(sps[0],sps[1])
        arrayTag = '%sx%s'%(arraysPair[0],arraysPair[1])
        
        m = liteMap.liteMapFromFits(patchDir+os.path.sep+'T_map_%s_%s_0'%(arrays[0],seasons[0]))
        window_T = m.copy()
        window_T.data[:] = 1.0

    	l, f_l0 = numpy.transpose(numpy.loadtxt(os.environ['SPECK_DIR']+'/data/'+p['beam_%s_%s'%(arraysPair[0],sps[0])]))
        l, f_l1 = numpy.transpose(numpy.loadtxt(os.environ['SPECK_DIR']+'/data/'+p['beam_%s_%s'%(arraysPair[1],sps[1])]))

        f_l = numpy.sqrt(f_l0*f_l1)
        
        transfer= [l,f_l**2]

        if pwDict != None and pwDict['apply']:
            
            transfer_pw=makePrewithenerTransfer(window_T,p['trimAtL'],alpha=pwDict['addBackFraction'])
            transfer += [transfer_pw]
            pickle.dump(transfer,open('mcm/transfer_%s_%s.pkl'%(arrayTag,spTag),'w'))

        else:
            pickle.dump(transfer,open('mcm/transfer_%s_%s.pkl'%(arrayTag,spTag),'w'))


        if cosineApod['apply']:
            window_T=liteMapPol.initializeCosineWindow(m,cosineApod['lenApod'],cosineApod['pad'])


        if p['applyPixelWeights']:
            if useSmoothedWeight==True:
                wt=liteMap.liteMapFromFits(patchDir+'/totalSmoothedWeightMap_%s_%s'%(arraysPair[0],sps[0]))
                wt1=liteMap.liteMapFromFits(patchDir+'/totalSmoothedWeightMap_%s_%s'%(arraysPair[1],sps[1]))
                wt.data[:]+=wt1.data[:]
            
            else:
                wt = liteMap.liteMapFromFits(patchDir+os.path.sep+'T_map_%s_%s_0'%(arrays[0],seasons[0]))
                wt.data[:] = 0.
                for j in xrange(nDivs):
                    wtmap0 = liteMap.liteMapFromFits(patchDir+'/weight_%s_%s_%d'%(arraysPair[0],sps[0],j))
                    wtmap1 = liteMap.liteMapFromFits(patchDir+'/weight_%s_%s_%d'%(arraysPair[1],sps[1],j))
                    wt.data[:] += wtmap0.data[:]+wtmap1.data[:]
                    
                print 'Writing  %s'%(patchDir+'/totalWeightMap_%s_%s'%(arrayTag,spTag))
                wt.writeFits(patchDir+'/totalWeightMap_%s_%s'%(arrayTag,spTag),overWrite=True)
                window_T.data[:] *= wt.data[:]


        window_Pol=window_T.copy()
        
        if  p['applyMask'] ==True:
            m=liteMap.liteMapFromFits(patchDir+os.path.sep+'mask')
            print "applying mask in "
            window_T.data[:] *= m.data[:]
        if p['maskPol']!=None:
            maskPol=liteMap.liteMapFromFits(patchDir+os.path.sep+'maskPol')
            window_Pol.data[:] *= maskPol.data[:]

        window_T.writeFits('auxMaps/finalWindow_%s_%s_T.fits'%(arrayTag,spTag), overWrite=True)
        window_Pol.writeFits('auxMaps/finalWindow_%s_%s_Pol.fits'%(arrayTag,spTag), overWrite=True)

		
        p2dWeight_TT = fftTools.powerFromLiteMap(window_T)
        p2dWeight_TT.powerMap[:] = 1.

        if (p['verticalkMaskLimits'] or p['horizontalkMaskLimits']) != None:
            p2dWeight_TT.createKspaceMask(verticalStripe=p['verticalkMaskLimits'],horizontalStripe=p['horizontalkMaskLimits'])
            p2dWeight_TT.powerMap[:] *= p2dWeight_TT.kMask[:]


        p2dWeight_TT=p2dWeight_TT.trimAtL(2*p['trimAtL']+500)
        #No special treatment of polarization for the moment
        p2dWeight_Cross= p2dWeight_TT.copy()
        p2dWeight_Pol= p2dWeight_TT.copy()


        numpy.save('noiseAndWeights/weightMap_%s_%s_TT.npy'%(arrayTag,spTag), p2dWeight_TT.trimAtL(p['trimAtL']+500).powerMap[:])
        numpy.save('noiseAndWeights/weightMap_%s_%s_Cross.npy'%(arrayTag,spTag), p2dWeight_Cross.trimAtL(p['trimAtL']+500).powerMap[:])
        numpy.save('noiseAndWeights/weightMap_%s_%s_Pol.npy'%(arrayTag,spTag), p2dWeight_Pol.trimAtL(p['trimAtL']+500).powerMap[:])


        print " ... done"
        print "use multiple windows"


        mcmPol.generateMCM(window_T,window_Pol,binLower[id],binUpper[id],p['trimAtL'],mbbFilename="mcm/"+p['mcmFileRoot'],transfer = transfer,binningWeightMap = p2dWeight_TT.powerMap,type=0,pixWin=p['pixWin'],arrayTag=arrayTag,spTag=spTag,ellStart=p['ellStart'])

            


            



