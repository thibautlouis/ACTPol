from flipper import *

l,beam=numpy.loadtxt('planck_dust_beam.txt',unpack=True)

arcmin_to_rad= (1./60.)*(numpy.pi/180.)
beam_size = 6.1 * arcmin_to_rad

bl=numpy.exp(-l*(l+1)*(beam_size**2)/(16.*numpy.log(2)))

g = open('BeamFDS.dat',mode="w")
for k in xrange(len(l)):
    g.write("%f %e \n"%(l[k],bl[k]))
        
g.close()

pylab.plot(l,beam**2)
pylab.plot(l,bl**2)
pylab.show()
