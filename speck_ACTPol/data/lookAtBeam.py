from flipper import *
import healpy

l,bl=numpy.loadtxt('beam_Planck_353_Ghz_pixWin.dat',unpack=True)


bl1,bl2,bl3=numpy.loadtxt('HFI_beam_window.dat',unpack=True)
l3=numpy.arange(len(bl3))

pw=healpy.sphtfunc.pixwin(2048)
pw=pw[:len(bl3)]
pylab.plot(l,bl)
pylab.plot(l3,(bl3*pw))
pylab.show()

g = open('beamSteve.dat',mode="w")

for k in xrange(len(l3)):
    g.write("%f %e\n"%(l3[k],(bl3*pw)[k]))
        
g.close()
              
