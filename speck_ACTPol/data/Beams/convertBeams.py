from flipper import *


def writeBinnedSpectrum(lbin,clbin,fileName):
    
    g = open(fileName,mode="w")
    
    for k in xrange(len(lbin)):
        g.write("%f %e \n"%(lbin[k],clbin[k]))
    
    g.close()




l1,bl1,error1=numpy.loadtxt('uranus_pa1_2014.txt',unpack=True)
l2,bl2,error2=numpy.loadtxt('uranus_pa2_2014.txt',unpack=True)
l,bl=numpy.loadtxt('beam_tform_140206_jitter_deep6_thibaut.txt',unpack=True)
beam6=numpy.loadtxt('beam_tform_140206_jitter_deep6.txt')
beam5=numpy.loadtxt('beam_tform_140206_jitter_deep5.txt')
l6,bl6=beam6[:,0],beam6[:,1]
l5,bl5=beam5[:,0],beam5[:,1]


#pylab.plot(l,bl)
pylab.plot(l1,bl1,label='AR1 marius beam (2014)')
pylab.plot(l2,bl2,label='AR2 marius beam (2014)')
pylab.plot(l5,bl5,label='AR1 matthew beam deep5 (2013)')
pylab.plot(l6,bl6,label='AR1 matthew beam deep6 (2013)')
pylab.xlabel(r'$\ell$',fontsize=22)
pylab.ylabel(r'$B_\ell$',fontsize=22)

pylab.legend()
pylab.show()


id=numpy.where(l1<=10000)
l1,bl1=l1[id],bl1[id]

id=numpy.where(l2<=10000)
l2,bl2=l2[id],bl2[id]

writeBinnedSpectrum(l1,bl1,'uranus_pa1_2014_thibaut.txt')
writeBinnedSpectrum(l2,bl2,'uranus_pa2_2014_thibaut.txt')
