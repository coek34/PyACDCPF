"""Prints power flow results related to DC network and converters.
"""

from sys import stdout

from numpy import ones, zeros, r_, c_, sort, exp, pi, diff, arange, min, \
    max, argmin, argmax, logical_or, real, imag, any, abs

from numpy import flatnonzero as find

from pyacdcpf.idx_busdc import BUSAC_I, GRIDDC, BUSDC_I, VDC, PDC
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC, PFDC, PTDC
from pyacdcpf.idx_convdc import CONV_BUS, CONVSTATUS, CONVTYPE_DC, DCSLACK,\
        DCDROOP, DCNOSLACK, CONVTYPE_AC, PVC, PQC, PCONV, QCONV, VCONV, \
        DROOP, PDCSET, VDCSET, DVDCSET, BASEKVC, LOSSA, LOSSB, LOSSCR, \
        LOSSCI, ICMAX, VCMAX, VCMIN, RCONV, XCONV, RTF, XTF, BF, PFIL, \
        QCONVF, PCCONV, QCCONV, QCCONVF, PCLOSS, VMC, VAC
        

def printdcpf(busdc, convdc, branchdc):
    """
    Prints power flow results related to DC network and converters.

    Prints all ac/dc power flow results related to dc grid and converters.

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ##define other numbers and indices
    nconv   = convdc.shape[0]
    nbusdc  = busdc.shape[0]
        
    ## dc bus data
    stdout.write('\n================================================================================')
    stdout.write('\n|     DC bus data                                                              |')
    stdout.write('\n================================================================================')
    stdout.write('\n Bus   Bus   Voltage    Power')
    stdout.write('\n DC #  AC #  Mag(pu)    P (MW)')
    stdout.write('\n-----  ----  ---------  --------')
    for i in arange(nbusdc):
        if busdc[i,BUSAC_I]==0:
            stdout.write('\n%4d%6c%10.3f%11.3f' % (busdc[i,BUSDC_I], '-', tuple(busdc[i,[VDC,PDC]])))
        else:
            stdout.write('\n%4d%6d%10.3f%11.3f' % tuple(busdc[i,[BUSDC_I,BUSAC_I,VDC,PDC]]))
    stdout.write('\n')

    ## transformer losses
    Plosstf = abs(convdc[:,PFIL] - convdc[:,PCONV])
    Qlosstf = abs(convdc[:,QCONVF] - convdc[:,QCONV])

    ## reactor losses
    Plossc = abs(convdc[:, PCCONV] - convdc[:, PFIL])
    Qlossc = abs(convdc[:, QCCONV] - convdc[:, QCCONVF])

    ## converter data
    Plosstot = Plosstf + Plossc + convdc[:,PCLOSS]

    ## filter reactive power
    Qfilt = convdc[:,QCCONVF] - convdc[:,QCONVF]

    
    stdout.write('\n================================================================================')
    stdout.write('\n|     VSC Converter Data                                                       |')
    stdout.write('\n================================================================================')
    stdout.write('\n Bus     Bus injection           Converter Voltage                 Total loss  ' )
    stdout.write('\n DC#   P (MW)   Q (MVAr)         Mag(pu) Ang(deg)                    P (MW)    ' )
    stdout.write('\n-----  -------  --------         ------- --------                  ----------- ' )
    for i in arange(nconv):
        stdout.write('\n%4d%9.2f%9.2f%17.3f%9.3f%27.2f'%(tuple(r_[convdc[i,[CONV_BUS,PCONV,QCONV,VMC,VAC]],  Plosstot[i]])))
    stdout.write('\n                                                                      ------  ')
    stdout.write('\n                                                           Total: %9.2f '% sum(Plosstot))
    stdout.write('\n')
 
    stdout.write('\n Bus  Converter power   Filter   Transfo loss     Reactor loss    Converter loss')
    stdout.write('\n DC#  P (MW) Q (MVAr)  Q (MVAr)  P (MW) Q (MVAr)  P (MW) Q (MVAr)    P (MW) ')
    stdout.write('\n----- ------- -------  --------  ------ --------  ------ -------- --------------')
    for i in arange(nconv):
        stdout.write('\n%4d%9.2f%8.2f%8.2f%9.2f%8.2f%9.2f%8.2f%12.2f'% \
             (tuple(r_[convdc[i,[CONV_BUS,PCCONV,QCCONV]], Qfilt[i], Plosstf[i], Qlosstf[i], Plossc[i], Qlossc[i], convdc[i,PCLOSS]])))

    stdout.write('\n                                 ------ --------  ------ -------- --------------');
    stdout.write('\n                      Total: %9.2f%8.2f%9.2f%8.2f%12.2f' % \
    (sum(Plosstf), sum(Qlosstf), sum(Plossc), sum(Qlossc), sum(convdc[:,PCLOSS])))
    stdout.write('\n')

    stdout.write('\n Bus  Grid power       Traf Filt.Power  Filter    Conv Filt. Pwr   Converter Power')
    stdout.write('\n DC#  P (MW) Q (MVAr)  P (MW) Q (MVAr)  Q (MVAr)  Q (MVAr)         P (MW) Q (MVAr)')
    stdout.write('\n----- ------ --------  ------ --------  --------  --------------   ------ --------')
    for i in arange(nconv):
        stdout.write('\n%4d%9.2f%8.2f%8.2f%9.2f%8.2f%11.2f%16.2f%8.2f'% \
            (tuple(r_[convdc[i,[CONV_BUS,PCONV,QCONV,PFIL,QCONVF]], convdc[i,QCCONVF] - convdc[i,QCONVF], convdc[i,[QCCONVF,PCCONV,QCCONV]]])))

    ## dc branch data
    nbranch = branchdc.shape[0]

    P_max = max(c_[abs(branchdc[:,PFDC]),abs(branchdc[:,PTDC])],1)
    P_min = min(c_[abs(branchdc[:,PFDC]),abs(branchdc[:,PTDC])],1)
    Plossline = P_max-P_min

    stdout.write('\n================================================================================')
    stdout.write('\n|     DC branch data                                                           |')
    stdout.write('\n================================================================================')

    stdout.write('\nBrnch   From   To     From Bus    To Bus      Loss')
    stdout.write('\n  #     Bus    Bus    P (MW)      P (MW)      P (MW)  ')
    stdout.write('\n-----  -----  -----  --------    --------    -------- ')
    for i in range(nbranch):
        stdout.write('\n%4d%7d%7d%10.2f%12.2f%12.2f' % \
            (tuple(r_[i, branchdc[i,F_BUSDC], branchdc[i,T_BUSDC], branchdc[i,[PFDC,PTDC]], Plossline[i]
            ])))
    stdout.write('\n                                             --------')
    stdout.write('\n                                 Total: %12.2f' % \
            sum(Plossline))
    stdout.write('\n')
    stdout.write('\n')