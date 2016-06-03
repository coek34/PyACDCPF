"""
Defines constants for named column indices to dc bus matrix.

Some examples of usage, after defining the constants using the line above,
are::

    Pd = busdc[4, PDC];        # get dc power withdrawal from dc grid at dc bus with index 4
    busdc[:, VDCMIN] = 0.95;   # set min voltage magnitude to 0.95 at all dc buses

The index, name and meaning of each column of the branch matrix is given
below:

columns 0-9 must be included in input matrix (in case file)
    0.  C{BUSDC_I}     dc bus number
    1.  C{BUSAC_I}     ac bus number (corresponding) 0 indicates no ac bus connection
    2.  C{GRIDDC}      dc grid to which the BUSDC_I is connected
    3.  C{PDC}         power withdrawn from the dc grid (MW)
    4.  C{VDC}         dc voltage (p.u.)
    5.  C{BASE_KVDC}   base dc voltage (kV)
    6.  C{VDCMAX}      max dc voltage (p.u.)
    7.  C{VDCMIN}      min dc voltage (p.u.)
    8.  C{CDC}         dc bus capacitor size (p.u.) (not used in power flow)
    
@author:Jef Beerten (KU Leuven)
@author:Roni Irnawan (Aalborg University)    
"""

# define the indices
BUSDC_I     = 0      # dc bus number
BUSAC_I     = 1      # ac bus number (corresponding) 0 indicates no ac bus connection
GRIDDC      = 2      # dc grid to which BUSDC_I is connected
PDC         = 3      # power withdrawal from the dc grid (MW)
VDC         = 4      # dc voltage (p.u.)
BASE_KVDC   = 5      # base dc voltage (kV)
VDCMAX      = 6      # max dc voltage (p.u.)
VDCMIN      = 7      # min dc voltage (p.u.)
CDC         = 8      # dc capacitor size (p.u.)
