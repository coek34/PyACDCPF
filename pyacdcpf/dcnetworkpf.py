"""Runs the dc network power flow.
"""

from numpy import ones, finfo, arange, abs, zeros
from numpy import flatnonzero as find

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from pypower.idx_bus import VM, VA, PD, QD
from pypower.idx_gen import GEN_BUS, GEN_STATUS, PG, QG, QMIN, QMAX
from pypower.idx_brch import F_BUS, T_BUS, BR_STATUS, PF, PT, QF, QT

eps = finfo(float).eps


def dcnetworkpf(Ybusdc, Vdc, Pdc, slack, noslack, droop, PVdroop, Pdcset, \
        Vdcset, dVdcset, pol, tol, itmax):
    """
    Runs the dc network power flow.
    
    Runs the dc network power flow, possibly including several dc grids.
    Each dc networks can have dc slack buses or several converters in dc
    voltage control.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ## initialisation
    nb = Vdc.size ## number of dc busses
    Pdc1 = -Pdc ## convention on power flow direction
    Pdc1[droop] = -Pdcset[droop] ## droop power set-points
    drooplidx = droop-ones(droop.size)*nb+droop # linear droop indices

    ##----- dc network iteration -----
    ## initialisation
    it = 0
    converged = 0

    ## Newton-Raphson iteration
    while (not converged) and it<=itmax:
        ## update iteration counter
        it += 1

        ## calculate power injections and Jacobian matrix
        Pdccalc = pol*Vdc*(Ybusdc*Vdc)
        Pdccalc = Pdccalc.reshape(nb)
        Vdc = Vdc.reshape((nb,1)) # Reshape in order to do the matrix multiplication
        J = csr_matrix(pol*Ybusdc.multiply(Vdc@Vdc.T))
        Vdc = Vdc.reshape(nb) # Switch back to the original shape
        J[arange(J.shape[0]),arange(J.shape[0])]      = J.diagonal() + Pdccalc #replace matrix elements

        ## include droop characteristics
        Vdcsetlh = (abs(Vdc-Vdcset)<=dVdcset)*Vdc + \
                ((Vdc-Vdcset)>dVdcset)*(Vdcset + dVdcset) + \
                ((Vdc-Vdcset)<-dVdcset)*(Vdcset-dVdcset)    ## define set-point with deadband

        Pdccalc[droop] = Pdccalc[droop] + 1/PVdroop[droop]*(Vdc[droop]-Vdcsetlh[droop]) # droop addition
        J[drooplidx,drooplidx] = J[drooplidx,drooplidx] + 1/PVdroop[droop]*Vdc[droop]

        ## dc network solution
        # reduce Jacobian
        keep = zeros(J.shape[0], dtype=bool)
        keep[noslack] = True # keep the "noslack" Jacobian
        Jr = J[keep] 
        Jr = Jr[:,keep] 
        dPdcr = Pdc1[noslack] - Pdccalc[noslack] ## power mismatch vector
        dPdcr = dPdcr.reshape((dPdcr.size,1))
        dVr = spsolve(Jr,dPdcr) ##voltage corrections

        ## update dc voltages
        Vdc[noslack] = Vdc[noslack]*(ones(noslack.size)+dVr)

        ## convergence check
        if abs(dVr).max()<tol: 
            converged = 1

    ## convergence print
    if not converged:
        stdout.write('\nDC network power flow did NOT converge after %d iterations\n', it)


    ##----- Output update -----
    ## recalculate slack bus powers
    Vdc = Vdc.reshape((nb,1)) # Reshape in order to do the matrix multiplication    
    Pdc1[slack] = pol*Vdc[slack]*(Ybusdc[slack,:]@Vdc);
    Pdc[slack]  = -Pdc1[slack]

    ## recalculate voltage droop bus powers
    Pdc1[droop]        = pol*Vdc[droop]*(Ybusdc[droop,:]@Vdc);
    Pdc[droop]         = -Pdc1[droop];  
    Vdc = Vdc.reshape(nb) # Switch back to the original shape
    
    return Vdc, Pdc
