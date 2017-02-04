#!/usr/bin/env python
from flipper import *
import liteMapPol 
import speckMisc

        
    
print "Reading dict file"
p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])

specDir='spectra'


spec=p['spec']




seasons = p['seasonTags']
arrays = p['arrayTags']

seasonPairs = []
for i in xrange(len(seasons)):
    for j in xrange(len(seasons)):
        if arrays[0] == arrays[-1] and i>j: continue
        seasonPairs += [[seasons[i],seasons[j]]]




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
        
            clBin_crossMean=[]
            clBin_autoMean=[]
            
            spec1=spec
            
            if spec=='d5xd56':
                if spTag=='s1xs1':
                    spec1='d5xd5'
                if spTag=='s1xs2_a':
                    spec1='d5xd56_p5_1'
                if spTag=='s1xs2_b':
                    spec1='d5xd56_p5_2'
                if spTag=='s2_axs2_a':
                    spec1='d56_p5_1xd56_p5_1'
                if spTag=='s2_bxs2_b':
                    spec1='d56_p5_2xd56_p5_2'
                if spTag=='s2_axs2_b':
                    spec1='d56_p5_1xd56_p5_2'
                        
            if spec=='d6xd56':
                if spTag=='s1xs1':
                    spec1='d6xd6'
                if spTag=='s1xs2_a':
                    spec1='d6xd56_p6_1'
                if spTag=='s1xs2_b':
                    spec1='d6xd56_p6_2'
                if spTag=='s2_axs2_a':
                    spec1='d56_p6_1xd56_p6_1'
                if spTag=='s2_bxs2_b':
                    spec1='d56_p6_2xd56_p6_2'
                if spTag=='s2_axs2_b':
                    spec1='d56_p6_1xd56_p6_2'

            try:
                os.makedirs('spectra_%s'%spec1)
            except:
                pass


            lBin,clBin_cross,binWeight=numpy.loadtxt('%s/clBinCrossMean_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag),unpack=True)
            fName = 'spectra_%s/clBinCrossMean_%s%s_%s.dat'%(spec1,l1,l2,spec1)
            speckMisc.writeBinnedSpectrum(lBin,clBin_cross,binWeight,fName)
                
            if sps[0] == sps[1]:
                lBin,clBin_auto,binWeight=numpy.loadtxt('%s/clBinAutoMean_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag),unpack=True)
                fName = 'spectra_%s/clBinAutoMean_%s%s_%s.dat'%(spec1,l1,l2,spec1)
                speckMisc.writeBinnedSpectrum(lBin,clBin_auto,binWeight,fName)

            lBin,nlBin,binWeight=numpy.loadtxt('%s/noise_%s%s_%sx%s_%s.dat'%(specDir,l1,l2,arrays[0],arrays[-1],spTag),unpack=True)
            fName = 'spectra_%s/nlBin_%s%s_%s.dat'%(spec1,l1,l2,spec1)
            speckMisc.writeBinnedSpectrum(lBin,nlBin,binWeight,fName)


