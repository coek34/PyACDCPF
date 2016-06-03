"""Internal slack/droop bus power injection iteration.
"""

from sys import stdout, stderr

from numpy import ones, finfo, arange, abs, angle, zeros, real, imag, where, \
            cos, sin, r_, exp, sqrt, delete
from numpy import flatnonzero as find

from scipy.sparse import diags, hstack, vstack
from scipy.sparse.linalg import spsolve

from pypower.idx_bus import VM, VA, PD, QD
from pypower.idx_gen import GEN_BUS, GEN_STATUS, PG, QG, QMIN, QMAX
from pypower.idx_brch import F_BUS, T_BUS, BR_STATUS, PF, PT, QF, QT

## define j
## DONT USE j IN ANYWHERE ELSE!!!
j = sqrt(-1+0j)

eps = finfo(float).eps # for avoiding division by zero


def calcslackdroop(Pcspec, Qsspec, Vs, Vf, Vc, Ztf, Bf, Zc, tol, itmax):
    """
    Internal slack/droop bus power injection iteration.
    
    Internal Newton-Raphson iteration to calculate the slack buses and
    voltage droop controlled buses power injections in the ac grid using
    the converter active power injection and the ac grid state (Vs) and the
    reactive power injection as fixed values.
    
    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """
    
    ##----- initialise -----
    ## define number of slack/droop buses
    ng = Pcspec.size

    ## define matrix indices
    i1 = arange(0,ng)
    i2 = arange(ng,2*ng)
    i3 = arange(2*ng,3*ng)
    i4 = arange(3*ng,4*ng)

    ## define voltage amplitudes and angles
    Vsm = abs(Vs)      ## grid voltage amplitude
    Vsa = angle(Vs)    ## grid voltage angle
    Vfm = abs(Vf)      ## filter voltage amplitude
    Vfa = angle(Vf)    ## filter voltage angle
    Vcm = abs(Vc)      ## converter voltage amplitude
    Vca = angle(Vc)    ## converter voltage angle

    ## calculate converter admitance
    Yc = 1/Zc
    Gc = real(Yc)
    Bc = imag(Yc)

    ## determine converters without and with transformers
    tf0i = where(Ztf==0)[0]
    tf1i = where(Ztf!=0)[0]

    ## calculate transformer admitance
    Ytf = (Ztf!=0)*1/(Ztf + eps*ones(Ztf.size)) + zeros(Ztf.size)
    Gtf = real(Ytf)
    Btf = imag(Ytf)

    ##----- Vc slack bus iteration -----
    ## initialise Jacobian matrix
    # J         = spalloc(4*ng, 4*ng,14*ng);

    ## initialisation
    it = 0
    converged = 0
    cflag     = 0

    ## Newton-Raphson iteration
    while not converged and it<itmax:
        ## update iteration counter
        it += 1

        ## define cos and sin of angles
        cosfc = cos(Vfa - Vca)
        sinfc = sin(Vfa - Vca)
        cossf = cos(Vsa - Vfa)
        sinsf = sin(Vsa - Vfa)
        cossc = cos(Vsa - Vca)
        sinsc = sin(Vsa - Vca)

        ## converter side power
        Pc  =  Vcm**2*Gc - Vfm*Vcm*(Gc*cosfc - Bc*sinfc)
        Qc  = -Vcm**2*Bc + Vfm*Vcm*(Gc*sinfc + Bc*cosfc)

        ## filter side converter power
        Pcf = -Vfm**2*Gc  + Vfm*Vcm*(Gc*cosfc + Bc*sinfc)
        Qcf =  Vfm**2*Bc  + Vfm*Vcm*(Gc*sinfc - Bc*cosfc)

        ## filter reactive power
        Qf  = -Bf*Vfm**2

        ## filter side grid power
        Psf =  Vfm**2*Gtf - Vfm*Vsm*(Gtf*cossf - Btf*sinsf)
        Qsf = -Vfm**2*Btf + Vfm*Vsm*(Gtf*sinsf + Btf*cossf)

        ## grid side power
        Ps  = (-Vsm**2*Gtf + Vfm*Vsm*(Gtf*cossf + Btf*sinsf))*(Ztf!=0) + \
              (-Vsm**2*Gc + Vsm*Vcm*(Gc*cossc + Bc*sinsc))*(Ztf==0)
        Qs  = (Vsm**2*Btf + Vfm*Vsm*(Gtf*sinsf - Btf*cossf))*(Ztf!=0) + \
              (Vsm**2*(Bc+Bf) + Vsm*Vcm*(Gc*sinsc - Bc*cossc))*(Ztf==0)

        ## additional filter bus equations
        F1 = Pcf - Psf
        F2 = Qcf - Qsf - Qf

        mismatch = r_[Pcspec-Pc,
                      Qsspec-Qs,
                      -F1,
                      -F2]
        mismatch1 = mismatch - r_[zeros(3*ng),mismatch[i4]*(Ztf==0)]
        mismatch1 = mismatch1 - r_[zeros(2*ng),mismatch[i3]*(Ztf==0),zeros(ng)]
        if abs(mismatch1).max()<tol:
            cflag = 1
            break

        ## Jacobian matrix elements
        J11 = diags(-Qc - Vcm**2*Bc)                ## J(i1,i1)
        J12 = diags((Qc + Vcm**2*Bc)*(Ztf!=0))      ## J(i1,i2)
        J13 = diags(Pc + Vcm**2*Gc)                 ## J(i1,i3)
        J14 = diags((Pc - Vcm**2*Gc)*(Ztf!=0))      ## J(i1,i4)

        ## only included without transformer
        J21 = diags((-Ps - Vsm**2*Gc)*(Ztf==0))     ##J(i2,i1)
        J23 = diags((Qs - Vsm**2*(Bc+Bf))*(Ztf==0)) ##J(i2,i3)

        ## only included with transformer
        J22 = diags((-Ps - Vsm**2*Gtf)*(Ztf!=0))    ##J(i2,i2)
        J24 = diags((Qs - Vsm**2*Btf)*(Ztf!=0))     ##J(i2,i4)

        J31 = diags((Qcf - Vfm**2*Bc)*(Ztf!=0))                 ##J(i3,i1)
        J32 = diags((-Qcf + Qsf + Vfm**2*(Bc+Btf))*(Ztf!=0))    ##J(i3,i2)
        J33 = diags((Pcf + Vfm**2*Gc)*(Ztf!=0))                 ##J(i3,i3)
        J34 = diags((Pcf - Psf - Vfm**2*(Gc+Gtf))*(Ztf!=0))     ##J(i3,i4)

        J41 = diags((-Pcf - Vfm**2*Gc)*(Ztf!=0))                    ##J(i4,i1)
        J42 = diags((Pcf - Psf + Vfm**2*(Gc+Gtf))*(Ztf!=0))         ##J(i4,i2)
        J43 = diags((Qcf - Vfm**2*Bc)*(Ztf!=0))                     ##J(i4,i3)
        J44 = diags((Qcf - Qsf + Vfm**2*(Bc+Btf+2*Bf))*(Ztf!=0))    ##J(i4,i4)

        J = vstack([
                hstack([J11, J12, J13, J14]),
                hstack([J21, J22, J23, J24]),
                hstack([J31, J32, J33, J34]),
                hstack([J41, J42, J43, J44])
            ],format="csr")

        ## remove rows and colums of transformerless buses
        rem4 = ones(J.shape[0], dtype=bool)
        rem4[tf0i] = False # remove transformerless buses
        rem3 = ones(J.shape[0], dtype=bool)
        rem3[tf0i] = False # remove transformerless buses
        rem2 = ones(J.shape[0], dtype=bool)
        rem2[tf0i] = False # remove transformerless buses
        J = J[rem4,:]
        J = J[rem3,:]
        J = J[:,rem4]
        J = J[:,rem2]
        mismatch = delete(mismatch, i4[tf0i])
        mismatch = delete(mismatch, i3[tf0i])

        ## calculate correction terms
        corr1 = spsolve(J,mismatch)

        ## add rows to mismatch matrix of transformerless buses
        corr = zeros(4*ng)
        corr[r_[i1,i2[tf1i],i3,i4[tf1i]].T] = corr1

        ## update converter voltage magnitude and angle
        Vca = Vca + corr[i1]
        Vfa = Vfa + corr[i2]
        Vcm = Vcm*(ones(ng) + corr[i3])
        Vfm = Vfm*(ones(ng) + corr[i4])

    ## convergence print
    if not cflag:
        stdout.write('\nSlackbus converter power calculation did NOT converge in %d iterations\n' %it)

    ## define cos and sin of angles
    cosfc = cos(Vfa - Vca)
    sinfc = sin(Vfa - Vca)
    cossf = cos(Vsa - Vfa)
    sinsf = sin(Vsa - Vfa)


    ##----- Output update -----
    ## slack bus VSC grid injection active power
    Ps = (-Vsm**2*Gtf + Vfm*Vsm*(Gtf*cossf + Btf*sinsf))*(Ztf!=0) + \
          (-Vsm**2*Gc + Vsm*Vcm*( Gc*cossc + Bc*sinsc ))*(Ztf==0)

    ## slack bus converter side reactive power
    Qc = -Vcm**2*Bc + Vfm*Vcm*(Gc*sinfc + Bc*cosfc)

    ## slack bus converter voltage
    Vc = Vcm*exp(j*Vca)

    return Ps, Qc, Vc
