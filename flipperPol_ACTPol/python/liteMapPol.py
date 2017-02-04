from flipper import *
import fftPol
import numpy
import pyfits
from scipy.interpolate import splrep,splev
import systematicAndBeam


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
    return p2d.powerMap[:]

    
def padLiteMap(Small,Big):
    DeltaX=numpy.int(Big.Nx-Small.Nx)
    DeltaY=numpy.int(Big.Ny-Small.Ny)
    lowX=DeltaX/2
    lowY=DeltaY/2
    if DeltaY/2.!=DeltaY/2:
        lowY=lowY+1
    if DeltaX/2.!=DeltaX/2:
        lowX=lowX+1
    
    padMap=Big.copy()
    padMap.data[:]=0
    padMap.data[lowY:Big.Ny-DeltaY/2,lowX:Big.Nx-DeltaX/2]=Small.data[:]
    
    return(padMap)



def makeEmptyCEATemplate(raSizeDeg, decSizeDeg,meanRa = 180., meanDec = 0.,\
                      pixScaleXarcmin = 0.5, pixScaleYarcmin=0.5):
    assert(meanDec == 0.,'mean dec other than zero not implemented yet')


    cdelt1 = -pixScaleXarcmin/60.
    cdelt2 = pixScaleYarcmin/60.
    naxis1 = numpy.int(raSizeDeg/pixScaleXarcmin*60.+0.5)
    naxis2 = numpy.int(decSizeDeg/pixScaleYarcmin*60.+0.5)
    refPix1 = naxis1/2.
    refPix2 = naxis2/2.
    pv2_1 = 1.0
    cardList = pyfits.CardList()
    cardList.append(pyfits.Card('NAXIS', 2))
    cardList.append(pyfits.Card('NAXIS1', naxis1))
    cardList.append(pyfits.Card('NAXIS2', naxis2))
    cardList.append(pyfits.Card('CTYPE1', 'RA---CEA'))
    cardList.append(pyfits.Card('CTYPE2', 'DEC--CEA'))
    cardList.append(pyfits.Card('CRVAL1', meanRa))
    cardList.append(pyfits.Card('CRVAL2', meanDec))
    cardList.append(pyfits.Card('CRPIX1', refPix1+1))
    cardList.append(pyfits.Card('CRPIX2', refPix2+1))
    cardList.append(pyfits.Card('CDELT1', cdelt1))
    cardList.append(pyfits.Card('CDELT2', cdelt2))
    cardList.append(pyfits.Card('CUNIT1', 'DEG'))
    cardList.append(pyfits.Card('CUNIT2', 'DEG'))
    hh = pyfits.Header(cards=cardList)
    wcs = astLib.astWCS.WCS(hh, mode='pyfits')
    data = numpy.zeros([naxis2,naxis1])
    ltMap = liteMap.liteMapFromDataAndWCS(data,wcs)

    return ltMap
    
    

def initializeCosineWindow(liteMap,lenApod,pad):
	
	Nx=liteMap.Nx
	Ny=liteMap.Ny
	win=liteMap.copy()
	win.data[:]=1

	winX=win.copy()
	winY=win.copy()

	for j in range(pad,Ny-pad):
		for i in range(pad,Nx-pad):
			if i<=(lenApod+pad):
				r=float(i)-pad
				winX.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
			if i>=(Nx-1)-lenApod-pad:
				r=float((Nx-1)-i-pad)
				winX.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))

	for i in range(pad,Nx-pad):
		for j in range(pad,Ny-pad):
			if j<=(lenApod+pad):
				r=float(j)-pad
				winY.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
			if j>=(Ny-1)-lenApod-pad:
				r=float((Ny-1)-j-pad)
				winY.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
	
	win.data[:]*=winX.data[:,:]*winY.data[:,:]
        win.data[0:pad,:]=0
        win.data[:,0:pad]=0
        win.data[Ny-pad:Ny,:]=0
        win.data[:,Nx-pad:Nx]=0
	
	return(win)
	

def makeMask(liteMap,nHoles,holeSize,lenApodMask,show=True):
        
    pixScaleArcmin=liteMap.pixScaleX*60*360/numpy.pi
    holeSizePix=numpy.int(holeSize/pixScaleArcmin)

    mask=liteMap.copy()
    mask.data[:]=1
    holeMask=mask.copy()
    
    Nx=mask.Nx
    Ny=mask.Ny
    xList=numpy.random.rand(nHoles)*Nx
    yList=numpy.random.rand(nHoles)*Ny    
    
    for k in range(nHoles):
    	print "number of Holes",k
        holeMask.data[:]=1
        for i in range(Nx):
            for j in range(Ny):
            	rad=(i-numpy.int(xList[k]))**2+(j-numpy.int(yList[k]))**2
            	
            	if rad < holeSizePix**2:
                    holeMask.data[j,i]=0
                for pix in range(lenApodMask):
                	
                    if rad <= (holeSizePix+pix)**2 and rad > (holeSizePix+pix-1)**2:
                        	holeMask.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*float(pix)/lenApodMask))
        mask.data[:]*=holeMask.data[:]
    data=mask.data[:]
    
    if show==True:
    	pylab.matshow(data)
    	pylab.show()
    
        
    return mask


    


def initializeDerivativesWindowfuntions(liteMap):
	
    def matrixShift(l,row_shift,column_shift):	
        m1=numpy.hstack((l[:,row_shift:],l[:,:row_shift]))
        m2=numpy.vstack((m1[column_shift:],m1[:column_shift]))
        return m2
    delta=liteMap.pixScaleX
    Win=liteMap.data[:]
    
    dWin_dx=(-matrixShift(Win,-2,0)+8*matrixShift(Win,-1,0)-8*matrixShift(Win,1,0)+matrixShift(Win,2,0))/(12*delta)
    dWin_dy=(-matrixShift(Win,0,-2)+8*matrixShift(Win,0,-1)-8*matrixShift(Win,0,1)+matrixShift(Win,0,2))/(12*delta)
    d2Win_dx2=(-matrixShift(dWin_dx,-2,0)+8*matrixShift(dWin_dx,-1,0)-8*matrixShift(dWin_dx,1,0)+matrixShift(dWin_dx,2,0))/(12*delta)
    d2Win_dy2=(-matrixShift(dWin_dy,0,-2)+8*matrixShift(dWin_dy,0,-1)-8*matrixShift(dWin_dy,0,1)+matrixShift(dWin_dy,0,2))/(12*delta)
    d2Win_dxdy=(-matrixShift(dWin_dy,-2,0)+8*matrixShift(dWin_dy,-1,0)-8*matrixShift(dWin_dy,1,0)+matrixShift(dWin_dy,2,0))/(12*delta)
    
    #In return we change the sign of the simple gradient in order to agree with numpy convention
    return {'Win':Win, 'dWin_dx':-dWin_dx,'dWin_dy':-dWin_dy, 'd2Win_dx2':d2Win_dx2, 'd2Win_dy2':d2Win_dy2,'d2Win_dxdy':d2Win_dxdy}
	
	

def  simPolMapsFromEandB(Temp,l,cl_TT,cl_EE,cl_TE,cl_BB=None,fullBeamMatrix=None,beam1d=None,undoTE=False):
    
    
    buffer=1
    Ny = Temp.Ny
    Nx = Temp.Nx
    
    modLMap,angLMap=fftPol.makeEllandAngCoordinate(Temp,buffer)
	
    ll = numpy.ravel(modLMap)
    
    clCorr_EE,clUncorr_EE= fftPol.generate_EE_Power(cl_TT,cl_TE,cl_EE)
    
    
    #pCorr_EE = fftPol.generateKspacePower(Temp,buffer,l,clCorr_EE,ll,Gauss=True)
    #pUnCorr_EE = fftPol.generateKspacePower(Temp,buffer,l,clUncorr_EE,ll,Gauss=True)
    #p_TT= fftPol.generateKspacePower(Temp,buffer,l,cl_TT,ll,Gauss=True)
    #  p_EE= fftPol.generateKspacePower(Temp,buffer,l,cl_EE,ll,Gauss=True)
    
    p_TT,p_EE,p_TE=fftPol.generateAllKspacePower(Temp,buffer,l,cl_TT,cl_EE,cl_TE,ll,Gauss=True)
    
    
    randomReal_TT=numpy.random.randn(Ny,Nx)
    randomIm_TT=numpy.random.randn(Ny,Nx)
    
    randomReal_EE=numpy.random.randn(Ny,Nx)
    randomIm_EE=numpy.random.randn(Ny,Nx)
    
    realPart_T = numpy.sqrt(p_TT)*randomReal_TT
    imgPart_T = numpy.sqrt(p_TT)*randomIm_TT
    kMap_T = realPart_T+1j*imgPart_T
    
    realPart_E = numpy.sqrt(p_EE)*randomReal_EE
    imgPart_E = numpy.sqrt(p_EE)*randomIm_EE
    kMap_Eold= realPart_E+1j*imgPart_E
    
    ### add correlations - uses generate_2d_power
    
    
    loc2 = numpy.where(p_TT == 0.)
    p_TT[loc2]=1
    
    p_EE_corr = p_TE/numpy.sqrt(p_TT)
    
    
    if numpy.min(p_EE-p_TE**2/p_TT) <0:
        print "PROBLEM"
        sys.exit()
    
    
    p_EE_uncorr = numpy.sqrt(p_EE-p_TE**2/p_TT)
    
    area = Nx*Ny*Temp.pixScaleX*Temp.pixScaleY
    
    kMap_E = p_EE_corr*(randomReal_TT+1j*randomIm_TT)+ p_EE_uncorr*(randomReal_EE+1j*randomIm_EE)
    
    
    if undoTE!=None:
        if undoTE==True:
            print 'undoEE'
            kMap_E=kMap_Eold
    
    
    
    
    if cl_BB!=None:
    	p_BB= fftPol.generateKspacePower(Temp,buffer,l,cl_BB,ll,Gauss=True)
    	realPart_B = numpy.sqrt(p_BB)*numpy.random.randn(Ny,Nx)
    	imgPart_B = numpy.sqrt(p_BB)*numpy.random.randn(Ny,Nx)
    	kMap_B= realPart_B+1j*imgPart_B
    	kMap_Q=	kMap_E * numpy.cos(2*angLMap) - kMap_B*numpy.sin(2*angLMap)
    	kMap_U=	kMap_E * numpy.sin(2*angLMap) + kMap_B*numpy.cos(2*angLMap)
    
    
    else:
    	kMap_Q=	kMap_E * numpy.cos(2*angLMap)
    	kMap_U=	kMap_E * numpy.sin(2*angLMap)
    
    
    
    data_T = numpy.real(numpy.fft.ifft2(kMap_T))
    data_Q=  numpy.real(numpy.fft.ifft2(kMap_Q))
    data_U=  numpy.real(numpy.fft.ifft2(kMap_U))
    
    
    
    T_map=Temp.copy()
    Q_map=Temp.copy()
    U_map=Temp.copy()
    
    
    T_map.data=data_T
    Q_map.data=data_Q
    U_map.data=data_U
    
    
    
    return(T_map,Q_map,U_map)




def addWhiteNoise(map,rmsArcmin):
        """
        Adds white noise to a given map; returns a new map
        """
        noisyMap = map.copy()
        if rmsArcmin == 0.0:
            pass
	else:
            radToMin = 180/numpy.pi*60
            pixArea = radToMin**2 * map.pixScaleX*map.pixScaleY
            rms = rmsArcmin/numpy.sqrt(pixArea)

            noise = numpy.random.normal( scale = rms, size = map.data.shape )

            noisyMap.data[:] += noise[:]

        return noisyMap
        




def initializeDerivXCosineWindow(liteMap,lenApod):
	
	Nx=liteMap.Nx
	Ny=liteMap.Ny
	win=liteMap.copy()
	win.data[:]=1

	winX=win.copy()
	winX.data[:]=0
	winY=win.copy()

	for j in range(Ny):
		for i in range(Nx):
			if i<=lenApod:
				r=float(i)
				winX.data[j,i]=-numpy.pi/(2*lenApod*winX.pixScaleX)*(numpy.sin(-numpy.pi*r/lenApod))
			if i>=(Ny-1)-lenApod:
				r=float((Ny-1)-i)
				winX.data[j,i]=numpy.pi/(2*lenApod*winX.pixScaleX)*(numpy.sin(-numpy.pi*r/lenApod))
	
	
	for i in range(Nx):
		for j in range(Ny):
			if j<=lenApod:
				r=float(j)
				winY.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
			if j>=(Ny-1)-lenApod:
				r=float((Ny-1)-j)
				winY.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
	

	winX.data[:,:]*=winY.data[:,:]

	return(winX)
    

def initializeDerivYCosineWindow(liteMap,lenApod):
	
	Nx=liteMap.Nx
	Ny=liteMap.Ny
	win=liteMap.copy()
	win.data[:]=1

	winX=win.copy()
	winY=win.copy()
	winY.data[:]=0

	for j in range(Ny):
		for i in range(Nx):
			if i<=lenApod:
				r=float(i)
				winX.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
			if i>=(Ny-1)-lenApod:
				r=float((Ny-1)-i)
				winX.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))

	
	for i in range(Nx):
		for j in range(Ny):
			if j<=lenApod:
				r=float(j)
				winY.data[j,i]=-numpy.pi/(2*lenApod*winX.pixScaleY)*(numpy.sin(-numpy.pi*r/lenApod))
			if j>=(Ny-1)-lenApod:
				r=float((Ny-1)-j)
				winY.data[j,i]=numpy.pi/(2*lenApod*winX.pixScaleY)*(numpy.sin(-numpy.pi*r/lenApod))
	
	winX.data[:,:]*=winY.data[:,:]

	return(winX)


def initializeDerivXYCosineWindow(liteMap,lenApod):
	
	Nx=liteMap.Nx
	Ny=liteMap.Ny
	win=liteMap.copy()
	win.data[:]=1

	winX=win.copy()
	winX.data[:]=0
	winY=win.copy()
	winY.data[:]=0

	for j in range(Ny):
		for i in range(Nx):
			if i<=lenApod:
				r=float(i)
				winX.data[j,i]=-numpy.pi/(2*lenApod*winX.pixScaleX)*(numpy.sin(-numpy.pi*r/lenApod))
			if i>=(Ny-1)-lenApod:
				r=float((Ny-1)-i)
				winX.data[j,i]=numpy.pi/(2*lenApod*winX.pixScaleX)*(numpy.sin(-numpy.pi*r/lenApod))
	
	
	for i in range(Nx):
		for j in range(Ny):
			if j<=lenApod:
				r=float(j)
				winY.data[j,i]=-numpy.pi/(2*lenApod*winX.pixScaleY)*(numpy.sin(-numpy.pi*r/lenApod))
			if j>=(Ny-1)-lenApod:
				r=float((Ny-1)-j)
				winY.data[j,i]=numpy.pi/(2*lenApod*winX.pixScaleY)*(numpy.sin(-numpy.pi*r/lenApod))
	
	winX.data[:,:]*=winY.data[:,:]

	return(winX)
	


def initializeDerivXSquareCosineWindow(liteMap,lenApod):
	
	Nx=liteMap.Nx
	Ny=liteMap.Ny
	win=liteMap.copy()
	win.data[:]=1

	winX=win.copy()
	winX.data[:]=0
	winY=win.copy()

	for j in range(Ny):
		for i in range(Nx):
			if i<=lenApod:
				r=float(i)
				winX.data[j,i]=numpy.pi**2/(2*lenApod**2*winX.pixScaleX**2)*(numpy.cos(-numpy.pi*r/lenApod))
			if i>=(Ny-1)-lenApod:
				r=float((Ny-1)-i)
				winX.data[j,i]=numpy.pi**2/(2*lenApod**2*winX.pixScaleX**2)*(numpy.cos(-numpy.pi*r/lenApod))
	
	
	for i in range(Nx):
		for j in range(Ny):
			if j<=lenApod:
				r=float(j)
				winY.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
			if j>=(Ny-1)-lenApod:
				r=float((Ny-1)-j)
				winY.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
	

	winX.data[:,:]*=winY.data[:,:]

	return(winX)


def initializeDerivYSquareCosineWindow(liteMap,lenApod):
	
	Nx=liteMap.Nx
	Ny=liteMap.Ny
	win=liteMap.copy()
	win.data[:]=1

	winX=win.copy()
	winY=win.copy()
	winY.data[:]=0

	for j in range(Ny):
		for i in range(Nx):
			if i<=lenApod:
				r=float(i)
				winX.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))
			if i>=(Ny-1)-lenApod:
				r=float((Ny-1)-i)
				winX.data[j,i]=1./2*(1-numpy.cos(-numpy.pi*r/lenApod))

	
	for i in range(Nx):
		for j in range(Ny):
			if j<=lenApod:
				r=float(j)
				winY.data[j,i]=numpy.pi**2/(2*lenApod**2*winX.pixScaleX**2)*(numpy.cos(-numpy.pi*r/lenApod))
			if j>=(Ny-1)-lenApod:
				r=float((Ny-1)-j)
				winY.data[j,i]=numpy.pi**2/(2*lenApod**2*winX.pixScaleX**2)*(numpy.cos(-numpy.pi*r/lenApod))
	
	winX.data[:,:]*=winY.data[:,:]

	return(winX)





