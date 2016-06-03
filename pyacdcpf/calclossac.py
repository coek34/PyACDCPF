"""Converter loss calculation.
"""

from numpy import c_, r_, sign, where, sqrt, ones

def calclossac(Pc,Qc,Vc,lossa, lossb, losscr, lossci):
    """
    Converter loss calculation.

    Calculates the converter losses based on a loss model quadratically
    dependent on the converter current Ic

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)
    """

    ##----- Converter loss input data -----
    ## Per unit loss coefficients
    a = lossa                ## constant loss coefficient  [-]
    b = lossb                ## linear loss coefficient    [-]
    c = c_[losscr,lossci]      ## quadratic loss coefficient [-]

    ##----- Converter input data -----
    ## Define other coefficients
    nc = Pc.shape[0]     ## number of converters
    convmode = sign(Pc)       ## converter operation mode
    convmode[where(convmode==0)] = 1 ## in Matlab zero has positive sign
    rectifier = convmode>0
    inverter = convmode<0
    VMc = abs(Vc)
    c_mtx = rectifier*c[:,0]+inverter*c[:,1]

    ##----- Converter loss calculation -----
    ## Define other coefficients
    Ic = sqrt(Pc**2+Qc**2)/(VMc) ## reactor currents
    Ploss = a*ones(Ic.size)+b*Ic+c_mtx*Ic**2 ## reactor losses

    Ploss = Ploss.reshape(Ploss.size)

    return Ploss
