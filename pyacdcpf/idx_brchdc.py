"""
Defines constants for named column indices to dc branch matrix.

Some examples of usage, after defining the constants using the line above,
are::

    Ploss = branchdc[:, PFDC] + branchdc[:, PTDC] # compute real power loss vector

The index, name and meaning of each column of the branch matrix is given
below:

columns 0-8 must be included in input matrix (in case file)
    0.  C{F_BUSDC}     f, from bus number
    1.  C{T_BUSDC}     t, to bus number
    2.  C{BRDC_R}      r, resistance (p.u.)
    3.  C{BRDC_L}      l, inductance (p.u./s) (not used in power flow)
    4.  C{BRDC_C}      c, total line charging capacity (p.u.*s) (not used in power flow)
    5.  C{RATEDC_A}    rateA, MVA rating A (long term rating)
    6.  C{RATEDC_B}    rateB, MVA rating B (short term rating)
    7.  C{RATEDC_C}    rateC, MVA rating C (emergency rating)
    8.  C{BRDC_STATUS} initial branch status, 1 - in service, 0 - out of service

columns 9-10 are added to matrix after power flow or OPF solution
they are typically not present in the input matrix
     9. C{PFDC}        real power injected at "from" bus end (MW)
    10. C{PTDC}        real power injected at "to" bus end (MW)

@author:Jef Beerten (KU Leuven)
@author:Roni Irnawan (Aalborg University)    
"""

# define the indices
F_BUSDC     = 0     # f, from bus number
T_BUSDC     = 1     # t, to bus number
BRDC_R      = 2     # r, resistance (p.u.)
BRDC_L      = 3     # inductance (p.u./s)
BRDC_C      = 4     # c, total line charging capacity (p.u.*s)
RATEDC_A    = 5     # rateA, MVA rating A (long term rating)
RATEDC_B    = 6     # rateB, MVA rating B (short term rating)
RATEDC_C    = 7     # rateC, MVA rating C (emergency rating)
BRDC_STATUS = 8     # initial branch status, 1 - in service, 0 - out of service

# included in power flow solution, not necessarily in input
PFDC        =  9    # real power injected at "from" bus end (MW)
PTDC        = 10    # real power injected at "to" bus end (MW)