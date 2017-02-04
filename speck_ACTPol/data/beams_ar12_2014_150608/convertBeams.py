from flipper import *

def writeSpectrum(lbin,clbin,fileName):
    
    g = open(fileName,mode="w")
    
    for k in xrange(len(lbin)):
        g.write("%f %e \n"%(lbin[k],clbin[k]))
        
    g.close()

for i in [1,2]:
    data=numpy.loadtxt('beam_tform_150608_2014_pa%d_jitter.txt'%i)
    l,bl=data[:,0],data[:,1]

    writeSpectrum(l[:4001],bl[:4001],'beam_tform_150608_2014_pa%d_jitter_Planck.txt'%i)
    writeSpectrum(l[:10000],bl[:10000],'beam_tform_150608_2014_pa%d_jitter_Thibaut.txt'%i)
