#!/usr/bin/env python
from flipper import *


mapDir = 'patches'
try:
    os.makedirs(mapDir)
    
except:
    pass



patchDir = '../deep5_groundsub_iter1000_wpoly_az20/patches/'

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


chi2Array=[]
angArray=[]
count=0
angMatt=-2



for phi in [angMatt]:
    
    ang=phi*numpy.pi/(180)

    for j in range(nDivs):
    
        Q_map=liteMap.liteMapFromFits('../deep5_groundsub_iter1000_wpoly_az20/patches/Q_map_000_%d'%j)
        U_map=liteMap.liteMapFromFits('../deep5_groundsub_iter1000_wpoly_az20/patches/U_map_000_%d'%j)
        T_map=liteMap.liteMapFromFits('../deep5_groundsub_iter1000_wpoly_az20/patches/T_map_000_%d'%j)

        Q_map_new=Q_map.copy()
        U_map_new=U_map.copy()

        Q_map_new.data=Q_map.data*numpy.cos(2*ang)+U_map.data*numpy.sin(2*ang)
        U_map_new.data=-Q_map.data*numpy.sin(2*ang)+U_map.data*numpy.cos(2*ang)

        Q_map_new.writeFits('patches/Q_map_000_%d'%j,overWrite=True)
        U_map_new.writeFits('patches/U_map_000_%d'%j,overWrite=True)
