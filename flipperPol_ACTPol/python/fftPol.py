from flipper import *
import liteMapPol
import numpy
from scipy.interpolate import splrep,splev


def TQUtoPureTEB(T_map,Q_map,U_map,window_T,window_Pol, modLMap,angLMap,method='standard',fftType=None):

    print 'ONLY STANDARD METHOD IMPLEMENTED WITH MULTIPLE WINDOW'
    
    T_temp=T_map.copy()
    Q_temp=Q_map.copy()
    U_temp=U_map.copy()

    T_temp.data=T_map.data*window_T.data
    fT=fftTools.fftFromLiteMap(T_temp,fftType=fftType)
    
    Q_temp.data=Q_map.data*window_Pol.data
    fQ=fftTools.fftFromLiteMap(Q_temp,fftType=fftType)
    
    U_temp.data=U_map.data*window_Pol.data
    fU=fftTools.fftFromLiteMap(U_temp,fftType=fftType)
    
    fE=fT.copy()
    fB=fT.copy()
    
    fE.kMap[:]=fQ.kMap[:]*numpy.cos(2.*angLMap)+fU.kMap[:]*numpy.sin(2.*angLMap)
    fB.kMap[:]=-fQ.kMap[:]*numpy.sin(2.*angLMap)+fU.kMap[:]*numpy.cos(2.*angLMap)
    
    if method=='standard':
        return fT, fE, fB



def fourierTQU(T_map,Q_map,U_map):
    
    fT=fftTools.fftFromLiteMap(T_map)
    fQ=fftTools.fftFromLiteMap(Q_map)
    fU=fftTools.fftFromLiteMap(U_map)
    
    return(fT, fQ, fU)
    
    

def TQUtoFourierTEB(T_map,Q_map,U_map,window,modLMap,angLMap):

    T_map.data*=window.data
    Q_map.data*=window.data
    U_map.data*=window.data

    fT=fftTools.fftFromLiteMap(T_map)
    
    fQ=fftTools.fftFromLiteMap(Q_map)
        
    fU=fftTools.fftFromLiteMap(U_map)
    
    fE=fT.copy()
    fB=fT.copy()
    fE.kMap[:]=fQ.kMap[:]*numpy.cos(2.*angLMap)+fU.kMap[:]*numpy.sin(2.*angLMap)
    fB.kMap[:]=-fQ.kMap[:]*numpy.sin(2.*angLMap)+fU.kMap[:]*numpy.cos(2.*angLMap)
    
    return(fT, fE, fB)
    
    
def generateKspacePower(liteMap,bufferFactor,l,Cl,ll,Gauss=False):
    
    Ny = liteMap.Ny*bufferFactor
    Nx = liteMap.Nx*bufferFactor
    s=splrep(l,Cl,k=3)
    kk = splev(ll,s)
    id = numpy.where(ll>l.max())
    kk[id] = 0.
    area = Nx*Ny*liteMap.pixScaleX*liteMap.pixScaleY
    p = numpy.reshape(kk,[Ny,Nx])/area * (Nx*Ny)**2
    
    if Gauss==True:
        
        p2d=fftTools.powerFromLiteMap(liteMap)
        p2d_new=p2d.copy()
        n=1000
        
        p2d=p2d.trimAtL(l.max()+n-1)
        
        Diff_x=p2d_new.Nx-p2d.Nx
        Diff_y=p2d_new.Ny-p2d.Ny
        Low_x=Diff_x/2
        Low_y=Diff_y/2
        
        if Diff_y/2.!=Diff_y/2:
            Low_y=Low_y+1
        if Diff_x/2.!=Diff_x/2:
            Low_x=Low_x+1
        
        Nx=liteMap.Nx
        Ny=liteMap.Ny
        
        area = Nx*Ny*liteMap.pixScaleX*liteMap.pixScaleY
        larray = numpy.arange(numpy.int(l.max()-1))
        deltaLx = numpy.abs(p2d.modLMap[0,1] - p2d.modLMap[0,0])
        deltaLy = numpy.abs(p2d.modLMap[1,0] - p2d.modLMap[0,0])
        delta = numpy.min([deltaLx,deltaLy])/2.0
        modMap = numpy.ravel(p2d.modLMap)
        modMap = numpy.floor(modMap)
        binMap=0

        
        big_Cl=numpy.zeros(len(Cl)+n)
        big_Cl[:len(Cl)]=Cl
        
        
        for k in xrange(len(larray)):
            print k
            #gauss = modMap*0
            #id=numpy.where(modMap==larray[k])
            #gauss[id]=1.
            #binMap+=gauss*Cl[k]
            
            gauss = numpy.exp(-(larray-larray[k])**2./(2.*delta**2.))
            sum = gauss.sum()
                    
            gauss = numpy.exp(-(modMap-larray[k])**2./(2.*delta**2.))
            gauss /= sum
            binMap+=gauss*Cl[k]


        p=numpy.reshape(binMap,[p2d.Ny,p2d.Nx])/area * (Nx*Ny)**2

        p=numpy.fft.fftshift(p)

        p_fat=numpy.zeros((p2d_new.Ny,p2d_new.Nx))

        p_fat[Low_y:p2d_new.Ny-Diff_y/2, Low_x:p2d_new.Nx-Diff_x/2]=p


        p_fat=numpy.fft.ifftshift(p_fat)

    return(p_fat)


def generateAllKspacePower(liteMap,bufferFactor,l,cl_TT,cl_EE,cl_TE,ll,Gauss=False):
    print 'use Gauss'
    
    if Gauss==False:
        print "Not implemented yet"
    else:
        p2d=fftTools.powerFromLiteMap(liteMap)
        
        p2d_new=p2d.copy()
        n=1000
        
        p2d=p2d.trimAtL(l.max()+n-1)
        
        Diff_x=p2d_new.Nx-p2d.Nx
        Diff_y=p2d_new.Ny-p2d.Ny
        Low_x=Diff_x/2
        Low_y=Diff_y/2
        
        if Diff_y/2.!=Diff_y/2:
            Low_y=Low_y+1
        if Diff_x/2.!=Diff_x/2:
            Low_x=Low_x+1
        
        Nx=liteMap.Nx
        Ny=liteMap.Ny
        
        area = Nx*Ny*liteMap.pixScaleX*liteMap.pixScaleY
        larray = numpy.arange(numpy.int(l.max()-1))
        deltaLx = numpy.abs(p2d.modLMap[0,1] - p2d.modLMap[0,0])
        deltaLy = numpy.abs(p2d.modLMap[1,0] - p2d.modLMap[0,0])
        delta = numpy.min([deltaLx,deltaLy])/2.0
        modMap = numpy.ravel(p2d.modLMap)
        modMap = numpy.floor(modMap)
        
        binMap_TT=0
        binMap_EE=0
        binMap_TE=0
        
        big_Cl_TT=numpy.zeros(len(cl_TT)+n)
        big_Cl_TT[:len(cl_TT)]=cl_TT
        
        big_Cl_EE=numpy.zeros(len(cl_EE)+n)
        big_Cl_EE[:len(cl_EE)]=cl_EE
        
        big_Cl_TE=numpy.zeros(len(cl_TE)+n)
        big_Cl_TE[:len(cl_TE)]=cl_TE
        
        
        
        
        for k in xrange(len(larray)):
            print k
            #gauss = modMap*0
            #id=numpy.where(modMap==larray[k])
            #gauss[id]=1.
            
            gauss = numpy.exp(-(larray-larray[k])**2./(2.*delta**2.))
            sum = gauss.sum()
            
            gauss = numpy.exp(-(modMap-larray[k])**2./(2.*delta**2.))
            gauss /= sum
            
            binMap_TT+=gauss*big_Cl_TT[k]
            binMap_EE+=gauss*big_Cl_EE[k]
            binMap_TE+=gauss*big_Cl_TE[k]
        
        p_TT=numpy.reshape(binMap_TT,[p2d.Ny,p2d.Nx])/area * (Nx*Ny)**2
        p_EE=numpy.reshape(binMap_EE,[p2d.Ny,p2d.Nx])/area * (Nx*Ny)**2
        p_TE=numpy.reshape(binMap_TE,[p2d.Ny,p2d.Nx])/area * (Nx*Ny)**2
        
        p_TT=numpy.fft.fftshift(p_TT)
        p_EE=numpy.fft.fftshift(p_EE)
        p_TE=numpy.fft.fftshift(p_TE)
        
        
        p_TT_fat=numpy.zeros((p2d_new.Ny,p2d_new.Nx))
        p_EE_fat=numpy.zeros((p2d_new.Ny,p2d_new.Nx))
        p_TE_fat=numpy.zeros((p2d_new.Ny,p2d_new.Nx))
        
        p_TT_fat[Low_y:p2d_new.Ny-Diff_y/2, Low_x:p2d_new.Nx-Diff_x/2]=p_TT
        p_EE_fat[Low_y:p2d_new.Ny-Diff_y/2, Low_x:p2d_new.Nx-Diff_x/2]=p_EE
        p_TE_fat[Low_y:p2d_new.Ny-Diff_y/2, Low_x:p2d_new.Nx-Diff_x/2]=p_TE
        

        p_TT_fat=numpy.fft.ifftshift(p_TT_fat)
        p_EE_fat=numpy.fft.ifftshift(p_EE_fat)
        p_TE_fat=numpy.fft.ifftshift(p_TE_fat)
    
    

    return(p_TT_fat,p_EE_fat,p_TE_fat)



def generate_EE_Power(Cl_TT,Cl_TE,Cl_EE):
    
    clCorr_EE = Cl_TE/(Cl_TT**0.5)
    tmp= Cl_EE - Cl_TE**2/Cl_TT
    loc = numpy.where(tmp<0.)
    tmp[loc] = 0.
    clUncorr_EE = (tmp)**0.5
    
    return(clCorr_EE,clUncorr_EE)

	
def fourierTEBtoFourierTQU(fT,fE,fB,modLMap,angLMap):
    fQ=fT.copy()
    fU=fT.copy()
    fQ.kMap[:]=fE.kMap[:]*numpy.cos(2.*angLMap)-fB.kMap[:]*numpy.sin(2.*angLMap)
    fU.kMap[:]=fE.kMap[:]*numpy.sin(2.*angLMap)+fB.kMap[:]*numpy.cos(2.*angLMap)
    return fT,fQ,fU
    

def fourierTQUtoPowerTEB(fT,fQ,fU,modLMap,angLMap):
    fE=fT.copy()
    fB=fT.copy()
    fE.kMap[:]=fQ.kMap[:]*numpy.cos(2.*angLMap)+fU.kMap[:]*numpy.sin(2.*angLMap)
    fB.kMap[:]=-fQ.kMap[:]*numpy.sin(2.*angLMap)+fU.kMap[:]*numpy.cos(2.*angLMap)
    
    power_TT=fftTools.powerFromFFT(fT)
    power_TE=fftTools.powerFromFFT(fT,fE)
    power_TB=fftTools.powerFromFFT(fT,fB)
    power_EE=fftTools.powerFromFFT(fE)
    power_BE=fftTools.powerFromFFT(fB,fE)
    power_BB=fftTools.powerFromFFT(fB)
    
    return(power_TT,power_TE,power_TB,power_EE,power_BE,power_BB)
    
    
    
def fourierTEBtoPowerTEB(fT0,fE0,fB0,fT1,fE1,fB1,trimAtL):
    
    TT_power=fftTools.powerFromFFT(fT0,fT1,trimAtL)
    TE_power=fftTools.powerFromFFT(fT0,fE1,trimAtL)
    ET_power=fftTools.powerFromFFT(fE0,fT1,trimAtL)
    TB_power=fftTools.powerFromFFT(fT0,fB1,trimAtL)
    BT_power=fftTools.powerFromFFT(fB0,fT1,trimAtL)
    EE_power=fftTools.powerFromFFT(fE0,fE1,trimAtL)
    BE_power=fftTools.powerFromFFT(fB0,fE1,trimAtL)
    EB_power=fftTools.powerFromFFT(fE0,fB1,trimAtL)
    BB_power=fftTools.powerFromFFT(fB0,fB1,trimAtL)
    
    return(TT_power,TE_power,ET_power,TB_power,BT_power,EE_power,EB_power,BE_power,BB_power)


def fourierTEtoPowerTE(fT0,fE0,fT1,fE1):
    
    TT_power=fftTools.powerFromFFT(fT0,fT1)
    TE_power=fftTools.powerFromFFT(fT0,fE1)
    ET_power=fftTools.powerFromFFT(fE0,fT1)
    EE_power=fftTools.powerFromFFT(fE0,fE1)
    
    return(TT_power,TE_power,ET_power,EE_power)

    
def makeTemplate(m, wl, ell, maxEll, outputFile = None):
    """
    Yanked from Toby's csFilter
    For a given map (m) return a 2D k-space template from a 1D specification wl
    ell = 2pi * i / deltaX
    (m is not overwritten)
    """

    ell = numpy.array(ell)
    wl  = numpy.array(wl)
    
    p2d = fftTools.powerFromLiteMap(m)
    p2d.powerMap[:] = 0.

    l_f = numpy.floor(p2d.modLMap)
    l_c = numpy.ceil(p2d.modLMap)
    
    for i in xrange(numpy.shape(p2d.powerMap)[0]):
        for j in xrange(numpy.shape(p2d.powerMap)[1]):
            if l_f[i,j] > maxEll or l_c[i,j] > maxEll:
                continue
            w_lo = wl[l_f[i,j]]
            w_hi = wl[l_c[i,j]]
            trueL = p2d.modLMap[i,j]
            w = (w_hi-w_lo)*(trueL - l_f[i,j]) + w_lo
            p2d.powerMap[i,j] = w



            
    if outputFile != None:
        p2d.writeFits(outputFile, overWrite = True)
    return p2d
    
    
def makeEllandAngCoordinate(liteMap,bufferFactor=1):
	
    Ny = liteMap.Ny*bufferFactor
    Nx = liteMap.Nx*bufferFactor
    ly = numpy.fft.fftfreq(Ny,d = liteMap.pixScaleY)*(2*numpy.pi)
    lx = numpy.fft.fftfreq(Nx,d = liteMap.pixScaleX)*(2*numpy.pi)
    modLMap = numpy.zeros([Ny,Nx])
    angLMap = numpy.zeros([Ny,Nx])
    iy, ix = numpy.mgrid[0:Ny,0:Nx]
    modLMap[iy,ix] = numpy.sqrt(ly[iy]**2+lx[ix]**2)
    #Trigonometric orientation
    angLMap[iy,ix]= numpy.arctan2(lx[ix],ly[iy])


    return(modLMap,angLMap)
	
	


