"""Converts external to internal bus numbering.
"""
import sys

from numpy import where, setdiff1d, unique, r_, arange, zeros

from pyacdcpf.idx_busdc import BUSAC_I, BUSDC_I
from pyacdcpf.idx_convdc import CONV_BUS
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC
from pypower.idx_bus import BUS_I
from pypower.idx_gen import GEN_BUS
from pypower.idx_brch import F_BUS, T_BUS

def ext2intac(pdc,ppc):
    """
    Converts external to internal bus numbering.

    Converts external ac bus numbers (possibly non-consecutive) to
    consecutive internal bus numbers. The ac grid is sorted based on the
    converter connected buses and dc buses, per dc grid. Dummy ac buses are
    assigned to dc buses without a connection to the ac grid. All other ac
    buses are grouped at the end.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ##-----  rename ac bus numbers  -----
    ## check presence of converter station on ac busses
    acdum_i = where(pdc['busdc'][:,BUSAC_I] == 0)[0]
    accnv_i = where(pdc['busdc'][:,BUSAC_I] )[0]
    accnv = pdc['busdc'][accnv_i,BUSAC_I]
    acnocnv = setdiff1d(ppc['bus'][:,BUS_I],accnv)
    acdum = acnocnv[:acdum_i.shape[0]]
    acnodum = acnocnv[acdum_i.shape[0]:]

    ## check presence of multiple converters on ac busses
    if accnv.shape[0] != unique(accnv).shape[0] :
        sys.stderr.write('More than one converter per ac node detected!\n')

    ## define index matrices
    i2eac = r_[accnv,acdum,acnodum].astype(int)
    e2iac = zeros(max(i2eac) + 1)
    e2iac[i2eac] = arange(1,ppc['bus'].shape[0]+1)
    i2eac = r_[[0],i2eac]

    ## dummy ac bus additions to busdc matrix
    pdc['busdc'][acdum_i,BUSAC_I] = acdum

    ## rename ac busses
    pdc['busdc'][:, BUSAC_I] = e2iac[pdc['busdc'][:, BUSAC_I].astype(int)]
    ppc['bus'][:, BUS_I] = e2iac[ppc['bus'][:, BUS_I].astype(int)]
    ppc['gen'][:, GEN_BUS] = e2iac[ppc['gen'][:, GEN_BUS].astype(int)]
    ppc['branch'][:, F_BUS] = e2iac[ ppc['branch'][:, F_BUS].astype(int) ]
    ppc['branch'][:, T_BUS] = e2iac[ ppc['branch'][:, T_BUS].astype(int) ]

    return acdum_i, i2eac, pdc, ppc
