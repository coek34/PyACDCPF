"""Check for converter operation within ac voltage and current limits.
"""

from sys import stdout, stderr

from numpy import sqrt, abs, angle, real, imag, conj, pi, inf, c_, r_, ones,\
            cos, sin, arctan, array, zeros, intersect1d, where, isreal, round, finfo
from numpy.lib.scimath import arcsin, arccos

## define j
## DONT USE j IN ANYWHERE ELSE!!!
j = sqrt(-1+0j)

eps = finfo(float).eps # for avoiding division by zero

def convlim(Ss, Vs, Vc, Ztf, Bf, Zc, Icmax, Vcmax, Vcmin, convi, epslim, printopt):
    """
    Check for converter operation within ac voltage and current limits.

    Converter operation checked with the converter's PQ capability diagram
    that are calculated based on the converter limits and grid operation
    state.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ## voltage limits order check
    if Vcmax<Vcmin:
        stderr.write('Vcmin is larger than Vcmax for converter %d\n'%convi)
    elif Vcmax==Vcmin:
        stderr.write('Vcmin is equal to Vcmax for converter %d\n'%convi)

    Ssold = Ss.copy()

    ## define voltage magnitudes, angles, powers
    Vsm = abs(Vs)
    Vsa = angle(Vs)
    Vcm = abs(Vc)
    Vca = angle(Vc)
    Ps = real(Ss)
    Qs = imag(Ss)

    ##--- initialization ---
    ## load existing parameters
    Zf = 1/(j*(Bf + eps))
    Ytf = 1/(Ztf + eps) ## avoid division by zero. If Ztf=0, Ytf=inf is not used
    Yf  = j*(Bf + eps)
    Yc  = 1/(Zc + eps)

    ## pi-equivalent parameters of the converter station (voltage limits)
    # Implementation based on:
    # J. Beerten, S. Cole and R. Belmans: "Generalized Steady-State VSC MTDC
    # Model for Sequential AC/DC Power Flow Algorithms", IEEE Trans. Pow.
    # Syst., vol. 27, no. 2, 2012, pp. 821 - 829.
    #           ======
    #  --------|  Z2  |---------
    #    ||     ======     ||
    #   ====              ====
    #  |    |            |    |
    #  | Z1 |            | Z3 |
    #  |    |            |    |
    #   ====              ====
    #    ||                ||
    #
    # Expressions for Z3 and Y3 are not included, since not used in program

    if Ztf!=0 and Bf!=0:
        Z1 = (Ztf*Zc+Zc*Zf+Zf*Ztf)/Zc
        Z2 = (Ztf*Zc+Zc*Zf+Zf*Ztf)/Zf
    elif Ztf==0 and Bf!=0:
        Z1 = Zf
        Z2 = Zc
    elif Ztf!= 0 and Bf==0:
        Z1 = inf
        Z2 = Ztf+Zc
    else: ## Ztf == 0 && Bf ==0
        Z1 = inf
        Z2 = Zc

    Y1 = 1/Z1
    Y2 = 1/Z2
    G2 = real(Y2)
    B2 = imag(Y2)
    Y12 = Y1 + Y2
    G12 = real(Y12)
    B12 = imag(Y12)

    ##--- voltage and current limit parameters ----
    ## maximum current limit circle parameters
    MPL1 = -Vsm**2*(1/(conj(Zf) + conj(Ztf)))*(Bf!=0)\
            + (0+j*0)*(Bf==0)  ## center of the current limit
    rL1 = Vsm*Icmax*(\
                (Ztf!=0)*(abs(conj(Ytf)/(conj(Yf) + conj(Ytf))))\
                + (Ztf==0)*1) ## radius of current limit circle

    ## maximum and minimal active power on current limit
    PmaxL1 = real(MPL1) + rL1
    PminL1 = real(MPL1) - rL1
    QPmaxL1 = imag(MPL1)
    QPminL1 = imag(MPL1)

    ## minimum and maximum voltage limit circle parameters
    MPL2 = -Vsm**2*conj(Y1+Y2) ## center of the voltage limits
    rL2 =  Vsm*r_[Vcmin,Vcmax]*abs(Y2) ## radius of voltage limits

    ##--- voltage limit compliance of min/max active power point ---
    ## find intersection points of Vcmin/Vcmax and current limit circles
    dL12 = sqrt((real(MPL2)-real(MPL1))**2+(imag(MPL2)-imag(MPL1))**2)
    alpha = arctan((imag(MPL2)-imag(MPL1))/(real(MPL2)-real(MPL1)))
    beta = arccos(((dL12**2+rL1**2)*ones((1,2))-rL2**2)/(2*dL12*rL1))
    delta = arccos(((dL12**2-rL1**2)*ones((1,2))+rL2**2)/(2*dL12*rL2))
    gamma = c_[alpha-beta,pi-alpha-beta]
    eta = c_[alpha+delta,alpha-delta]

    ## possible intersection points (current limit) => intersection vector
    x1 = zeros((2,gamma.shape[1]),dtype=complex)
    y1 = zeros((2,gamma.shape[1]),dtype=complex)
    x1[0,:] = real(MPL1) + rL1*cos(gamma)
    x1[1,:] = real(MPL1) - rL1*cos(gamma)
    y1[0,:] = imag(MPL1) + rL1*sin(gamma)
    y1[1,:] = imag(MPL1) - rL1*sin(gamma)

    ## possible intersection points (voltage limits) => intersection vector
    x2 = zeros((2,gamma.shape[1]),dtype=complex)
    y2 = zeros((2,gamma.shape[1]),dtype=complex)
    x2[0,:] = real(MPL2) + r_[rL2,rL2]*cos(eta);
    x2[1,:] = real(MPL2) - r_[rL2,rL2]*cos(eta);
    y2[0,:] = imag(MPL2) + r_[rL2,rL2]*sin(eta);
    y2[1,:] = imag(MPL2) - r_[rL2,rL2]*sin(eta);

    ## vectorize intersection point matrices
    x1v = x1.reshape((x1.size,1),order='F')
    x2v = x2.reshape((x2.size,1),order='F')
    y1v = y1.reshape((y1.size,1),order='F')
    y2v = y2.reshape((y2.size,1),order='F')

    ## decrease accuracy to detect the intersection points
    eps2 = 1e-8;
    x1r = round(x1v/eps2)*eps2;
    x2r = round(x2v/eps2)*eps2;
    y1r = round(y1v/eps2)*eps2;
    y2r = round(y2v/eps2)*eps2;

    ## corresponding elements in intersections vectors
    x12r = intersect1d(x1r,x2r)
    x1i = array([idx for y in x12r for idx,x in enumerate(x1r) if x==y])
    x2i = array([idx for y in x12r for idx,x in enumerate(x2r) if x==y])
    y12r = intersect1d(y1r,y2r)
    y1i = array([idx for y in y12r for idx,x in enumerate(y1r) if x==y])
    y2i = array([idx for y in y12r for idx,x in enumerate(y2r) if x==y])

    x12r1 = c_[x1i,x12r]
    x12r1 = x12r1[x12r1[:,0].argsort()]
    x12r2 = c_[x2i,x12r]
    x12r2 = x12r2[x12r2[:,0].argsort()]
    y12r1 = c_[y1i,y12r]
    y12r1 = y12r1[y12r1[:,0].argsort()]
    y12r2 = c_[y2i,y12r]
    y12r2 = y12r2[y12r2[:,0].argsort()]

    ## define intersection points (full accuracy)
    VcminPQ1 = x1v[int(real(x12r1[0,0]))] + j*y1v[int(real(y12r1[0,0]))]
    VcminPQ2 = x1v[int(real(x12r1[2,0]))] + j*y1v[int(real(y12r1[2,0]))]
    VcmaxPQ1 = x1v[int(real(x12r1[1,0]))] + j*y1v[int(real(y12r1[1,0]))]
    VcmaxPQ2 = x1v[int(real(x12r1[3,0]))] + j*y1v[int(real(y12r1[3,0]))]
    
    ## Remove imaginary intersection points (no intersections, due to low/high voltage limits)
    if isreal(x12r2[0,1])==0 or isreal(y12r2[0,1])==0 or isreal(x12r2[2,1])==0 \
        or isreal(y12r2[2,1])==0 :
        VcminPQ1 = array([])     ## Imaginary intersection points are found in pairs
        VcminPQ2 = array([])
        if printopt == 1 :
            stdout.write('\n  Lower voltage limit at converter %d : No intersections with current limit were found.\n'%convi)

    if isreal(x12r2[1,1])==0 or isreal(y12r2[1,1])==0 or isreal(x12r2[3,1])==0 \
        or isreal(y12r2[3,1])==0 :
        VcmaxPQ1 = array([])
        VcmaxPQ2 = array([])
        if printopt == 1 :
            stdout.write('\n  Upper voltage limit at converter %d : No intersections with current limit were found.\n'%convi )

    ## Define maximum and minimum power points
    # Initialisation
    Pmin = PminL1
    QPmin = QPminL1
    Pmax = PmaxL1
    QPmax = QPmaxL1

    # Redefine max and min power points if min/max voltage limits are high/low
    if not VcminPQ1.size == 0 or not VcminPQ2.size == 0 :
        if printopt == 1 and (imag(VcminPQ1) > QPminL1 or imag(VcminPQ2) > QPmaxL1) :
            stdout.write('\n  High lower voltage limit detected at converter %d. \n'% convi)
        if imag(VcminPQ1) > QPminL1 :
            Pmin = real(VcminPQ1)
            QPmin = imag(VcminPQ1)
        if imag(VcminPQ2) > QPmaxL1:
            Pmax = real(VcminPQ2)
            QPmax = imag(VcminPQ2)
    if not VcmaxPQ1.size == 0 or not VcmaxPQ2.size == 0 :
        if printopt == 1 and (imag(VcmaxPQ1) < QPminL1 or imag(VcmaxPQ2) < QPmaxL1):
            stdout.write('\n  Low upper voltage limit detected at converter %d. \n '% convi)
        if imag(VcmaxPQ1) < QPminL1:
            Pmin = real(VcmaxPQ1)
            QPmin = imag(VcmaxPQ1)
        if imag(VcmaxPQ2) < QPmaxL1:
            Pmax = real(VcmaxPQ2)
            QPmax = imag(VcmaxPQ2)

    ##--- Limit check ---
    if Pmin < Ps and Ps < Pmax:
        ## maximum current limit (L1)
        if imag(MPL1)<Qs :
            Qs1 = imag(MPL1) + sqrt(rL1**2-(Ps-real(MPL1))**2)
        else: ## if Qs<imag(MPL1)
            Qs1 = imag(MPL1) - sqrt(rL1**2-(Ps-real(MPL1))**2)

        ## minimum and maximum voltage limits (L2)
        Qs2 = imag(MPL2)+sqrt(rL2**2-(Ps-real(MPL2))**2)
        a = 1+(B2/G2)**2*ones((1,2))
        b = -2*B2/G2*(Ps+Vsm**2*G12)/(Vsm*r_[Vcmin,Vcmax]*G2)
        c = ((Ps+Vsm**2*G12)/(Vsm*r_[Vcmin,Vcmax]*G2))**2-1

        ## only positive solution retained (neg. solution refers to lower part)
        sinDd = (-b + sqrt(b**2-4*a*c))/(2*a)
        Dd = arcsin(sinDd)
        cosDd = cos(Dd)
        Qs2 = Vsm**2*B12+Vsm*r_[Vcmin,Vcmax]*(G2*sinDd-B2*cosDd);

        ## adopt working point to limits
        if Qs > imag(MPL1) :
            if Qs > r_[Qs1,Qs2[:,1]].min() :
                convlimviol = 1
                Qs = r_[Qs1,Qs2[:,1]].min()
            elif Qs < Qs2[:,0] :
                convlimviol = 1
                Qs = Qs2[:,0]
            else: ## Qs < min([Qs1 Qs2[1]])
                convlimviol = 0
        else: ## Qs < imag(MPL1)
            if Qs < r_[Qs1,Qs2[:,0]].max():
                    convlimviol = 1
                    Qs = r_[Qs1,Qs2[:,0]].max()
            elif Qs > Qs2[:,1] :
                    convlimviol = 1
                    Qs = Qs2[:,1]
            else: ## Qs < min([Qs1 Qs2[1]])
                    convlimviol = 0

    ## active power outside of current limit active power range
    elif Ps <= Pmin: ## grid injected active power lower than minimum value
        ## set active power to minimum
        convlimviol = 2
        Ps = Pmin
        Qs = QPmin


    else:  ## if PmaxL1 <= Ps %% grid injected active power higher than maximum value
        ## set active power to maximum
        convlimviol = 2
        Ps = Pmax
        Qs = QPmax

    ## define output argument Ss
    Ss = Ps + j*Qs

    ## remove violation when difference is small
    if abs(Ssold-Ss)<epslim:
        convlimviol = 0

    ## define plot arguments
    if convlimviol == 1:
        SsIcmax = Ps+j*Qs1
        SsVcmin = Ps+j*Qs2.min()
        SsVcmax = Ps+j*Qs2.max()
        plotarg = r_[convlimviol, Ztf, Zf, Zc, Y1, Y2, Yf, Ytf, Vs, Icmax,\
                Vcmax, Vcmin, Ssold, Ss, SsIcmax, SsVcmax, SsVcmin]
    else:
        plotarg = r_[convlimviol, Ztf, Zf, Zc, Y1, Y2, Yf, Ytf, Vs, Icmax,\
                Vcmax, Vcmin, Ssold, Ss, 0, 0, 0];

    return convlimviol, Ss, plotarg
