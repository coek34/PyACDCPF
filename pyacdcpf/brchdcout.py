"""Split branchdc matrix into operating and non-operating lines.
"""

import sys

from numpy import where

from pyacdcpf.idx_brchdc import BRDC_STATUS

def brchdcout(pdc):
    """
    Split branchdc matrix into operating and non-operating lines.
    
    Returns seperate branch matrices for lines in operation and those out
    of operation, as well as their indices in the original brandch matrix.
    
    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)
    """
    
    ## converter status validity check
    if ((pdc['branchdc'][:,BRDC_STATUS]<0).any() 
            or (pdc['branchdc'][:,BRDC_STATUS]>1).any()):
        sys.stderr.write('branch status flags must be either 0 or 1\n')
    
    ## define index
    brchdc0i = where(pdc['branchdc'][:,BRDC_STATUS] == 0)[0]
    brchdc1i = where(pdc['branchdc'][:,BRDC_STATUS] == 1)[0]
    
    ## define converter outage matrix
    brchdc0 = pdc['branchdc'][brchdc0i,:]
    brchdc1 = pdc['branchdc'][brchdc1i,:]
    
    return brchdc1, brchdc1i, brchdc0, brchdc0i
    