"""Power flow data for 2 hvdc systems.
"""

from numpy import array

def case24_ieee_rts1996_MTDC():
    """Power flow data for 2 hvdc systems.
    3 node dc system (dc buses 1-3) and 4 node meshed dc system 
    (dc buses 4-7)
    can be used together with ac case files "case24_ieee_rts1996_....py"
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
        [1,      107,    1,      0,      1,      150,        1.1,    0.9,    0],
        [2,      204,    1,      0,      1,      150,        1.1,    0.9,    0],
        [3,      301,    1,      0,      1,      150,        1.1,    0.9,    0],
        [4,      113,    2,      0,      1,      300,        1.1,    0.9,    0],
        [5,      123,    2,      0,      1,      300,        1.1,    0.9,    0],
        [6,      215,    2,      0,      1,      300,        1.1,    0.9,    0],
        [7,      217,    2,      0,      1,      300,        1.1,    0.9,    0]
    ])

    ## DC converter data
    #   busdc_i type_dc type_ac P_g   Q_g   Vtar    rtf     xtf     bf     rc     xc     basekVac    Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  
    pdc["convdc"] = array([
        [1,      2,      1,        0  ,  50  , 1,    0.001 , 0.10, 0.09 , 0.0001,  0.16, 138,        1.2,    0.9,    1.1,    1,       1.103,0.887, 2.885,  4.371],
        [2,      1,      2,       75.3, -50  , 1,    0.001 , 0.10, 0.09 , 0.0001,  0.16, 138,        1.2,    0.9,    1.1,    1,       1.103,0.887, 2.885,  4.371],
        [3,      1,      1,     -141.9, 130  , 1,    0.001 , 0.05, 0.045, 0.0001,  0.08, 138,        1.2,    0.9,    2.2,    1,       2.206,0.887, 1.442,  2.185],
        [4,      2,      1,      131.5,  75.9, 1,    0.0005, 0.05, 0    , 0.0001,  0.08, 345,        1.2,    0.5,    2.2,    1,       2.206,1.8  , 5.94 ,  9    ],
        [5,      1,      1,      -61.7,   0  , 1,    0.001 , 0.10, 0    , 0.0001,  0.16, 345,        1.2,    0.5,    1.1,    1,       1.103,1.8  , 11.88,  18   ],
        [6,      1,      2,     -123.4, -10  , 1,    0.0005, 0.05, 0    , 0.0001,  0.08, 345,        1.2,    0.5,    2.2,    1,       2.206,1.8  , 5.94 ,  9    ],
        [7,      1,      1,       50  ,  20  , 1,    0.001 , 0.10, 0    , 0.0001,  0.16, 345,        1.2,    0.5,    1.1,    1,       1.103,1.8  , 11.88,  18   ]
    ])

    ## DC branch data
    #   fbusdc  tbusdc  r      l    c   rateA   rateB   rateC   status
    pdc["branchdc"] = array([
        [1,      3,      0.0352, 0, 0,   100,    100,    100,    1],
        [2,      3,      0.0352, 0, 0,   100,    100,    100,    1],
        [4,      5,      0.0828, 0, 0,   100,    100,    100,    1],
        [4,      7,      0.0704, 0, 0,   100,    100,    100,    1],
        [4,      6,      0.0718, 0, 0,   100,    100,    100,    1],
        [5,      7,      0.0760, 0, 0,   100,    100,    100,    1],
        [6,      7,      0.0248, 0, 0,   100,    100,    100,    1]
    ])

    return pdc