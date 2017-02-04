from flipper import *

pl=healpy.pixwin(4096)
l=numpy.arange(len(pl))

pylab.plot(l,pl)
pylab.show()

g = open('pixWin4096.dat',mode="w")
    
for k in xrange(len(l)):
    g.write("%f %e\n"%(l[k],pl[k]))
        
g.close()
              
