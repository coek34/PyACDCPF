"""Converts to internal bus numbering.
"""
import sys

from numpy import where, setdiff1d, unique, r_, arange, zeros

from pyacdcpf.idx_busdc import BUSAC_I, BUSDC_I
from pyacdcpf.idx_convdc import CONV_BUS
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC
from pypower.idx_bus import BUS_I
from pypower.idx_gen import GEN_BUS
from pypower.idx_brch import F_BUS, T_BUS

def int2extac(i2eac, acdum_i, pdc,ppc):
    """
    Converts to internal bus numbering.
    
    Converts consecutive internal ac bus numbers back to the original bus
    numbers using the mapping provided by I2EAC returned from EXT2INTAC and
    removes the dummy ac buses assigned to the dc buses without connection 
    to the ac grid as a result from EXT2INTAC.
    
    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ## rename bus numbers
    pdc['busdc'][:, BUSAC_I] = i2eac[pdc['busdc'][:, BUSAC_I].astype(int)]
    ppc['bus'][:, BUS_I] = i2eac[ppc['bus'][:, BUS_I].astype(int)]
    ppc['gen'][:, GEN_BUS] = i2eac[ppc['gen'][:, GEN_BUS].astype(int)]
    ppc['branch'][:, F_BUS] = i2eac[ ppc['branch'][:, F_BUS].astype(int) ]
    ppc['branch'][:, T_BUS] = i2eac[ ppc['branch'][:, T_BUS].astype(int) ]

    return pdc, ppc
