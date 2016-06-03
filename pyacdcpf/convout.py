"""Remove converters facing outages from converter matrix.
"""

from sys import stdout, stderr

from numpy import array, where

from pyacdcpf.idx_busdc import BUSDC_I, BUSAC_I
from pyacdcpf.idx_convdc import CONVSTATUS, CONV_BUS, PCONV, QCONV

def convout(pdc):
    """
    Remove converters facing outages from converter matrix.

    The dc converter matrix is split into working (C{conv1}) and non working
    (C{conv0}) converter matrices. The corresponding ac bus connection is
    removed from the busdc matrix.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ## converter status validity check
    if ((pdc['convdc'][:,CONVSTATUS]<0).any()
            or (pdc['convdc'][:,CONVSTATUS]>1).any()):
        stderr.write('converter status flags must be either 0 or 1\n')

    ## define index
    conv0i = where(pdc['convdc'][:,CONVSTATUS] == 0)[0]
    conv1i = where(pdc['convdc'][:,CONVSTATUS] == 1)[0]

    ## define converter outage matrix
    conv0 = pdc['convdc'][conv0i,:]
    conv1 = pdc['convdc'][conv1i,:]

    ## reset converter powers and voltages
    conv0[:,PCONV] = 0.
    conv0[:,QCONV] = 0.

    ## remove ac bus of converter with outage in busdc matrix
    if conv0i.shape[0] > 0:
        idx = [(where(x == pdc['busdc'][:,BUSDC_I])[0][0]) for x in conv0[:,CONV_BUS]]
        conv0busi = array([ idx, pdc['busdc'][idx,BUSAC_I]])
        pdc['busdc'][idx,BUSAC_I] = 0
    else:
        conv0busi = array([])

    return pdc, conv0busi, conv1, conv1i, conv0, conv0i
