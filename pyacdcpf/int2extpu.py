"""Converts internal per unit values to original values.
"""
import sys

from numpy import where

from pyacdcpf.idx_busdc import BASE_KVDC, CDC
from pyacdcpf.idx_convdc import RTF, XTF, BF, RCONV, XCONV, ICMAX
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC, BRDC_R, BRDC_L, BRDC_C

def int2extpu(baseMVA,pdc):
    """
    Converts internal per unit values to original values.
    
    Converts internal per unit values of the dc matrices (busdc, convdc
    and branchdc) to the original bases provided in the input files.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ##-----  per unit conversion  -----
    ## Only p.u. impedances and currents are changed. Voltages (p.u.) and
    ## powers (real values) are left unaltered.

    ## converter ac side conversion
    baseMVAac = pdc['baseMVAac']
    pdc['convdc'][:,RTF] *= baseMVAac/baseMVA
    pdc['convdc'][:,XTF] *= baseMVAac/baseMVA
    pdc['convdc'][:,BF] *= baseMVAac/baseMVA
    pdc['convdc'][:,RCONV] *= baseMVAac/baseMVA
    pdc['convdc'][:,XCONV] *= baseMVAac/baseMVA
    pdc['convdc'][:,ICMAX] *= baseMVA/baseMVAac

    ## ==> dc side per unit convention
    ## basekAdc = baseMVAdc/basekVdc
    ## baseRdc  = basekVdc^2/baseMVAdc

    ## converter dc side conversion
    baseMVAdc = pdc['baseMVAdc']
    baseR_busdc = pdc['busdc'][:,BASE_KVDC]**2/baseMVAdc
    baseR_busdc2ac = pdc['busdc'][:,BASE_KVDC]**2/baseMVA
    pdc['busdc'][:,CDC] *= 1/baseR_busdc2ac*baseR_busdc

    ## dc network branch conversion
    baseKVDC_brch = pdc['busdc'][F_BUSDC,BASE_KVDC]
    baseR_brchdc = baseKVDC_brch**2/baseMVAdc
    baseR_brchdc2ac = baseKVDC_brch**2/baseMVA
    pdc['branchdc'][:,BRDC_R] *= baseR_brchdc2ac/baseR_brchdc
    pdc['branchdc'][:,BRDC_L] *= baseR_brchdc2ac/baseR_brchdc
    pdc['branchdc'][:,BRDC_C] *= baseR_brchdc2ac/baseR_brchdc

    return pdc
