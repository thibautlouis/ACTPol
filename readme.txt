Code for flat sky polarisation power spectrum.
It uses a modified version of MASTER, computing exact mode coupling matrix.
If you are interested in B modes at very low noise level, you should switch to a pure B mode code.

Add in your .bashrc (change 'YOURLOC')

export FLIPPER_DIR='YOURLOC/flipper_ACTPol'
export PATH=$PATH:$FLIPPER_DIR/bin
export PYTHONPATH=$PYTHONPATH:$FLIPPER_DIR/python

export FLIPPERPOL_DIR='YOURLOC/flipperPol_ACTPol'
export PATH=$PATH:$FLIPPERPOL_DIR/bin
export PYTHONPATH=$PYTHONPATH:$FLIPPERPOL_DIR/python

export SPECK_DIR='YOURLOC/speck'
export PYTHONPATH=${PYTHONPATH}:${SPECK_DIR}/python:/Users/thibaut/Desktop/PythonLibrary/lib
export PATH=${PATH}:${SPECK_DIR}/bin

The mode coupling matrix is computed in fortran, you need to compile with f2py : 
FFLAGS="-fopenmp -fPIC -Ofast -ffree-line-length-none" f2py -c -m fortran fortran.f90 -lgomp 

Then to run the code 
HQcutPatches.py global.dict
HQcomputeMCM.py global.dict
HQcompileSpectra.py global.dict
HQcomputeAnalyticCovariance.py global.dict

Enjoy.

