"""Converts dc external to internal bus numbering.
"""
import sys

from numpy import gradient, sort, where, unique, any, arange, c_, r_, zeros

from pyacdcpf.idx_busdc import GRIDDC, BUSAC_I, BUSDC_I
from pyacdcpf.idx_convdc import CONV_BUS
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC

def ext2intdc(pdc):
    """
    Converts dc external to internal bus numbering.
    
    Converts external dc bus numbers (possibly non-consecutive) to
    consecutive internal bus numbers, starting at 1.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ##-----  Check grid numbering -----
    griddc = unique(pdc['busdc'][:,GRIDDC])

    if griddc.shape[0] > 1 and any(gradient(sort(griddc))>1.):
        sys.stderr.write('Non-successive dc grid numbering detected\n')

    ##-----  Permutation of dc bus matrix  -----
    ## Part 1: Group all dc busses without ac grid connection
    noacbusi = where(pdc['busdc'][:,BUSAC_I] == 0)[0]
    acbusi = where(pdc['busdc'][:,BUSAC_I] )[0]
    i2edcpmt = r_[acbusi,noacbusi]
    pdc['busdc'] = pdc['busdc'][i2edcpmt,:]

    ## Part 2: Sort dc busses based on dc grid number
    busdcext = c_[pdc['busdc'],i2edcpmt]
    busdcext = busdcext[busdcext[:,GRIDDC].argsort()]
    pdc['busdc'] = busdcext[:,:-1]
    i2edcpmt = busdcext[:,-1].astype(int)

    ##-----  Rename dc nodes  -----
    i2edc = pdc['busdc'][:, BUSDC_I].astype(int)
    e2idc = zeros(max(i2edc) + 1)
    e2idc[i2edc] = arange(1,pdc['busdc'].shape[0]+1)
    i2edc = r_[[0],i2edc]

    pdc['busdc'][:, BUSDC_I]    = e2idc[ pdc['busdc'][:, BUSDC_I].astype(int)    ]
    pdc['convdc'][:, CONV_BUS]  = e2idc[ pdc['convdc'][:, CONV_BUS].astype(int)  ]
    pdc['branchdc'][:, F_BUSDC] = e2idc[ pdc['branchdc'][:, F_BUSDC].astype(int) ]
    pdc['branchdc'][:, T_BUSDC] = e2idc[ pdc['branchdc'][:, T_BUSDC].astype(int) ]

    return i2edcpmt, i2edc, pdc
