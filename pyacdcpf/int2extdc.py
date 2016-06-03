"""Converts dc internal to external bus numbering.
"""
import sys

from numpy import gradient, sort, where, unique, any, arange, c_, r_, zeros

from pyacdcpf.idx_busdc import GRIDDC, BUSAC_I, BUSDC_I
from pyacdcpf.idx_convdc import CONV_BUS
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC

def int2extdc(i2edcpmt, i2edc, pdc):
    """
    Converts dc internal to external bus numbering.
    
    Converts external dc bus numbers (possibly non-consecutive) to 
    consecutive internal bus numbers, starting at 1.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ##----- internal to external bus numbering -----
    ## Part 1: Bus renumbering
    pdc['busdc'][:, BUSDC_I]    = i2edc[ pdc['busdc'][:, BUSDC_I].astype(int)    ]
    pdc['convdc'][:, CONV_BUS]  = i2edc[ pdc['convdc'][:, CONV_BUS].astype(int)  ]
    pdc['branchdc'][:, F_BUSDC] = i2edc[ pdc['branchdc'][:, F_BUSDC].astype(int) ]
    pdc['branchdc'][:, T_BUSDC] = i2edc[ pdc['branchdc'][:, T_BUSDC].astype(int) ]

    ## Part 2: Change bus order of busdc matrix
    pdc['busdc'][i2edcpmt,:] = pdc['busdc']

    return pdc
