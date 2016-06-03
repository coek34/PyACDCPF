"""Converts external per unit inputs to internal values.
"""
import sys

from numpy import where

from pyacdcpf.idx_busdc import BASE_KVDC, CDC
from pyacdcpf.idx_convdc import RTF, XTF, BF, RCONV, XCONV, ICMAX
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC, BRDC_R, BRDC_L, BRDC_C

def ext2intpu(baseMVA,pdc):
    """
    Converts external per unit inputs to internal values.
    
    Converts external per unit quantities of the dc bus, converter and
    branch matrix into internal per unit quantities using the ac
    power network base.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ##-----  per unit conversion  -----
    ## Only p.u. impedances and currents are changed. Voltages (p.u.) and
    ## powers (real values) are left unaltered.

    ## converter ac side conversion
    baseMVAac = pdc['baseMVAac']
    pdc['convdc'][:,RTF] *= baseMVA/baseMVAac
    pdc['convdc'][:,XTF] *= baseMVA/baseMVAac
    pdc['convdc'][:,BF] *= baseMVA/baseMVAac
    pdc['convdc'][:,RCONV] *= baseMVA/baseMVAac
    pdc['convdc'][:,XCONV] *= baseMVA/baseMVAac
    pdc['convdc'][:,ICMAX] *= baseMVAac/baseMVA

    ## ==> dc side per unit convention
    ## basekAdc = baseMVAdc/basekVdc
    ## baseRdc  = basekVdc^2/baseMVAdc

    ## converter dc side conversion
    baseMVAdc = pdc['baseMVAdc']
    baseR_busdc = pdc['busdc'][:,BASE_KVDC]**2/baseMVAdc
    baseR_busdc2ac = pdc['busdc'][:,BASE_KVDC]**2/baseMVA
    pdc['busdc'][:,CDC] *= 1/baseR_busdc*baseR_busdc2ac

    ## dc network branch conversion
    ## Check for equal base voltages at two sides of a dc branch
    Vff = pdc['busdc'][where(pdc['branchdc'][:,F_BUSDC])[0],BASE_KVDC]
    Vtt = pdc['busdc'][where(pdc['branchdc'][:,T_BUSDC])[0],BASE_KVDC]

    if not all(Vff==Vtt):
        sys.stderr.write('The dc voltages at both sides of a dc branch do not match.\n')

    baseKVDC_brch = pdc['busdc'][F_BUSDC,BASE_KVDC]
    baseR_brchdc = baseKVDC_brch**2/baseMVAdc
    baseR_brchdc2ac = baseKVDC_brch**2/baseMVA
    pdc['branchdc'][:,BRDC_R] *= baseR_brchdc/baseR_brchdc2ac
    pdc['branchdc'][:,BRDC_L] *= baseR_brchdc/baseR_brchdc2ac
    pdc['branchdc'][:,BRDC_C] *= baseR_brchdc/baseR_brchdc2ac

    return pdc
