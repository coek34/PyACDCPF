"""Plots converter limits and violations in PQ capability chart.
"""

from sys import stdout, stderr

from numpy import finfo, sqrt, abs, angle, real, imag, conj, pi, inf, arange, exp, \
                ones, zeros
from numpy.lib.scimath import arcsin, arccos

import matplotlib.pyplot as plt

## define j
## DONT USE j IN ANYWHERE ELSE!!!
j = sqrt(-1+0j)

eps = finfo(float).eps # for avoiding division by zero

def convlimplot(plotarg, convi):
    """
    Plots converter limits and violations in PQ capability chart.
    
    Graphical representation of converter limit violations
    PLOTARG : arguments for converter limit plot
    CONVI : converter number

        The PLOTARG elements are:
        idx - NAME      description
        ---   --------  ---------------------------------------
         0  - VIOL      converter limit violation type
                            0 - no limit hit
                            1 - reactive power violation
                            2 - active power violation
         1  - ZTF       transformer impedance
         2  - ZF        filter impedance
         3  - ZC        converter impedance
         4  - Y1        pi-equivalent admitance 1 (see convlim)
         5  - Y2        pi-equivalent admitance 2 (see convlim)
         6  - YF        filter admitance
         7  - YTF       transformer admitance
         8  - VS        grid voltage
         9  - ICMAX     converter current limit
        10  - VCMMAX    converter upper voltage limit
        11  - VCMMIN    converter lower voltage limit
        12  - SSOLD     grid side apparent power before limiting
        13  - SSNEW     grid side apparent power after limiting
        14  - SSICMAX   grid side apparent power: converter current
                        limit (Ps constant)
        15  - SSVCMAX   grid side apparent power: converter upper
                        voltage limit (Ps constant)
        16  - SSVCMIN   grid side apparent power: converter lower
                        voltage limit (Ps constant)
                        
    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ##--- initialisation ---
    ## define arguments
    viol    = plotarg[0]
    Ztf     = plotarg[1]
    Zf      = plotarg[2]
    Zc      = plotarg[3]
    Y1      = plotarg[4]
    Y2      = plotarg[5]
    Yf      = plotarg[6]
    Ytf     = plotarg[7]
    Vs      = plotarg[8]
    Icmax   = plotarg[9]
    Vcmmax  = plotarg[10]
    Vcmmin  = plotarg[11]
    Ssold   = plotarg[12]
    Ssnew   = plotarg[13]
    SsIcmax = plotarg[14]
    SsVcmax = plotarg[15]
    SsVcmin = plotarg[16]

    ## plot assumptions in code:
    ##      Vs = Vsm * exp(j*0) with Vsm constant (defined in input)
    ## converter upper and lower voltage limit:
    ##      Vc = Vcm * exp(Vca) with Vcm constant and Vca defined as variable
    ## converter current limit:
    ##      Ic = Icm * exp(Ica) with Icm constant and Ica defined as variable

    ## grid side voltage amplitude
    Vsm = abs(Vs)

    ##--- current limits ---
    ## current limit with current converter set-up
    i1 = 0.001 ## plot current angle step size
    Ica = arange(0,2*pi,i1)
    SsIlim  = Vs*Icmax*exp(-j*Ica)* \
              ((Ztf!=0 and Yf!=0)*(conj(Ytf)/(conj(Ytf) + conj(Yf))) + \
              (Ztf==0 or Yf ==0)*1) + \
              -Vsm**2*ones((1,Ica.size))*(1/(conj(Zf) + conj(Ztf)))
    SsIlim = SsIlim.reshape(SsIlim.size)
    SsIlimMP = -Vsm**2*(1/(conj(Zf) + conj(Ztf))) #middelpunt van figuur

    ## current limit without filter
    SsIlim1  = Vs*Icmax*exp(-j*Ica)   ## current limit

    ##--- voltage limits ---
    ## voltage limits with current converter set-up
    i2 = 0.001 ## plot voltage angle step size
    Vca = arange(0,2*pi,i2)
    SsVmax = -Vsm**2*conj(Y1 + Y2) + Vs*Vcmmax*exp(-j*Vca)*conj(Y2) ## upper voltage limit
    SsVmin = -Vsm**2*conj(Y1 + Y2) + Vs*Vcmmin*exp(-j*Vca)*conj(Y2) ## lower voltage limit

    ## voltage limits without filter
    SsVmax1 = -Vsm**2*conj(1/(Ztf + Zc)) + Vs*Vcmmax*exp(-j*Vca)*conj(1/(Ztf + Zc)) ## upper voltage limit (no filter bus)
    SsVmin1 = -Vsm**2*conj(1/(Ztf + Zc)) + Vs*Vcmmin*exp(-j*Vca)*conj(1/(Ztf + Zc)) ## lower voltage limit (no filter bus)


    ##--- plot converter limits ---
    ## plot options
    xmin = -1.5/1.2*Icmax.real
    xmax = 1.5/1.2*Icmax.real
    ymin = -1.5/1.2*Icmax.real
    ymax = 1.5/1.2*Icmax.real
    ix = 0.001
    iy = 0.001
    xx = arange(xmin.real,xmax.real,ix)
    yy = arange(ymin.real,ymax.real,ix)

    plt.figure()
    # hold on
    plt.axis('equal')
    plt.axis([xmin,xmax,ymin,ymax])
    plt.xlabel('$P_s$ (p.u.)')
    plt.ylabel('$Q_s$ (p.u.)')

    ## plot graph title
    if viol != 0:
        plt.title('Converter station ' + str(convi) + ' operating outside its limits.')
    elif viol == 0:
        plt.title('Normal operation of converter station ' + str(convi));

    ## plot axes
    plt.plot(xx,zeros(xx.size), 'k')
    plt.plot(zeros(xx.size),yy, 'k')

    ## center of current limit
    plt.scatter(real(SsIlimMP), imag(SsIlimMP), c='b', marker='+')

    ## old active and reactive powers
    plt.scatter(real(Ssold), imag(Ssold), c='k')

    ## new active and reactive powers
    plt.scatter(real(Ssnew), imag(Ssnew), c='k', marker='^')

    ## actual limit plots
    plt.plot(real(SsIlim), imag(SsIlim),'b')
    plt.plot(real(SsVmax), imag(SsVmax),'r')
    plt.plot(real(SsVmin), imag(SsVmin),'g')

    ## limit plots neglecting filter effect
    plt.plot(real(SsIlim1), imag(SsIlim1),'b:')
    plt.plot(real(SsVmax1), imag(SsVmax1),'r:')
    plt.plot(real(SsVmin1), imag(SsVmin1),'g:')

    ## plot limit points
    if viol == 1:
        plt.scatter(real(SsIcmax), imag(SsIcmax), c='b', marker='x');
        plt.scatter(real(SsVcmax), imag(SsVcmax), c='r', marker='x');
        plt.scatter(real(SsVcmin), imag(SsVcmin), c='g', marker='x');
    
    # show plot
    plt.show(block=False)
