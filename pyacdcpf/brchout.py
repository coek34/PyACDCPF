"""Split branch matrix into operating and non-operating lines.
"""

import sys

from numpy import where

from pypower.idx_brch import BR_STATUS

def brchout(ppc):
    """
    Split branch matrix into operating and non-operating lines.
    
    Returns seperate branch matrices for lines in operation and those out
    of operation, as well as their indices in the original brandch matrix.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)
    """
    
    ## converter status validity check
    if ((ppc['branch'][:,BR_STATUS]<0).any() 
            or (ppc['branch'][:,BR_STATUS]>1).any()):
        sys.stderr.write('branch status flags must be either 0 or 1\n')
    
    ## define index
    brch0i = where(ppc['branch'][:,BR_STATUS] == 0)[0]
    brch1i = where(ppc['branch'][:,BR_STATUS] == 1)[0]
    
    ## define converter outage matrix
    brch0 = ppc['branch'][brch0i,:]
    brch1 = ppc['branch'][brch1i,:]
    
    return brch1, brch1i, brch0, brch0i
    