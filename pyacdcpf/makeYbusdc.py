"""Builds the dc bus admittance matrix and branch admittance matrices.
"""

from sys import stderr

from numpy import ones, conj, nonzero, any, exp, pi, r_, where
from scipy.sparse import csr_matrix

from pyacdcpf.idx_brchdc import BRDC_R, F_BUSDC, T_BUSDC


def makeYbusdc(busdc, branchdc):
    """
    Builds the dc bus admittance matrix and branch admittance matrices.

    Returns the full bus admittance matrix (i.e. for all buses) and the
    matrices C{Yfdc} and C{Ytdc} which, when multiplied by the dc voltage
    vector, yield the vector currents injected into each dc line from the
    "from" and "to" buses respectively of each line.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ## constants
    n = busdc.shape[0]          ## number of buses
    m = branchdc.shape[0]       ## number of lines

    ## for each branch, compute the elements of the branch admittance matrix where
    ##
    ##      | If |   | Yff  Yft |   | Vf |
    ##      |    | = |          | * |    |
    ##      | It |   | Ytf  Ytt |   | Vt |
    ##

    Ys = ones(m)/branchdc[:,BRDC_R]

    Ytt = Ys
    Yff = Ys
    Yft = - Ys
    Ytf = - Ys

    ## build connection matrices
    f = branchdc[:, F_BUSDC]-1 ## list of index of "from" buses
    t = branchdc[:, T_BUSDC]-1 ## list of index of "to" buses
    ## connection matrix for line & from buses
    Cf = csr_matrix((ones(m), (range(m), f)), (m, n))
    ## connection matrix for line & to buses
    Ct = csr_matrix((ones(m), (range(m), t)), (m, n))

    ## build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    ## at each branch's "from" bus, and Yt is the same for the "to" bus end
    i = r_[range(m), range(m)]                   ## double set of row indices

    Yfdc = csr_matrix((r_[Yff, Yft], (i, r_[f, t])), (m, n))
    Ytdc = csr_matrix((r_[Ytf, Ytt], (i, r_[f, t])), (m, n))

    ## build Ybus
    Ybusdc = Cf.T * Yfdc + Ct.T * Ytdc

    return Ybusdc, Yfdc, Ytdc
