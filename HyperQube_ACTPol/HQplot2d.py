#!/usr/bin/env python
from flipper import *
import speckMisc
import pickle
import liteMapPol






p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])


specDir = 'spectra2d/'
patchDir= 'patches'

if len(sys.argv)> 2:
    specDir = 'spectra2d_%s/'%sys.argv[2]
    patchDir= 'patches_%s/'%sys.argv[2]


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


PlotDir='plot2d'
    
try:
    os.makedirs(PlotDir)
except:
    pass


fields=['T','E','B']



count1=0
for l1 in fields:
    
    count1+=1
    count2=0
    
    for l2 in fields:
        count2+=1
        if count2<count1: continue
        
        print '%s%s'%(l1,l2)
        
        for i in range(nPatches):

            p2d = pickle.load(open('%s/p2d_%s%s_%03d.pkl'%(specDir,l1,l2,i),mode="r"))
            
            
            if l1==l2:
                p2d.plot(log=True,show=False,pngFile='%s/p2d_%s%s_%03d.png'%(PlotDir,l1,l2,i))
                pylab.clf()
                pylab.close()
            
            else:
                p2d.plot(show=False,pngFile='%s/p2d_%s%s_%03d.png'%(PlotDir,l1,l2,i))
                pylab.clf()
                pylab.close()


