"""Power flow data for 5 bus, 2 generator case.
"""

from numpy import array

def case5_stagg_MTDCslack():
    """Power flow data for 3 node DC system.
    
    3 node system (constant voltage and power controlled) can be used
    together with ac case files 'case5_stagg.py' and 'case'3_inf.py'
    
    Network data from ...
    J. Beerten, D. Van Hertem, R. Belmans, "VSC MTDC systems with a 
    distributed DC voltage control – a power flow approach", in IEEE
    Powertech2011, Trondheim, Norway, Jun 2011.
    
    @return: Power flow data for 3 node DC system.
    """
    pdc = {"version": '1'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    pdc["baseMVAac"] = 100.0
    pdc["baseMVAdc"] = 100.0
    
    ## dc grid topology
    pdc["pol"] = 2.0

    ## DC bus data
    #   busdc_i busac_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc 
    pdc["busdc"] = array([
        [1,      2,      1,     0.0,    1.0,      345.0,      1.1,    0.9,    0],
        [2,      3,      1,     0.0,    1.0,      345.0,      1.1,    0.9,    0],
        [3,      5,      1,     0.0,    1.0,      345.0,      1.1,    0.9,    0]
    ])

    ## DC converter data
    #   busdc_i type_dc type_ac P_g   Q_g   Vtar    rtf     xtf     bf     rc     xc     basekVac    Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  
    pdc["convdc"] = array([
        [1,      1,      1,    -60., -40.,   1,    0.0015, 0.1121, 0.0887,0.0001,  0.16428, 345.,     1.1,    0.9,    1.2,    1,      1.103,0.887, 2.885,   4.371],
        [2,      2,      2,      0.,   0.,   1,    0.0015, 0.1121, 0.0887,0.0001,  0.16428, 345.,     1.1,    0.9,    1.2,    1,      1.103,0.887, 2.885,   4.371],
        [3,      1,      1,     35.,   5.,   1,    0.0015, 0.1121, 0.0887,0.0001,  0.16428, 345.,     1.1,    0.9,    1.2,    1,      1.103,0.887, 2.885,   4.371]
    ])

    ## DC branch data
    #   fbusdc  tbusdc  r      l    c   rateA   rateB   rateC   status
    pdc["branchdc"] = array([
        [1,      2,      0.052,  0,  0,   100,    100,    100,    1],
        [2,      3,      0.052,  0,  0,   100,    100,    100,    1],
        [1,      3,      0.073,  0,  0,   100,    100,    100,    1]
    ])

    return pdc