"""Power flow data for system with 3 zones consisting of 3 infinite
"""


def case3_inf():
    from numpy import array, inf, zeros
    """Power flow data for system with 3 zones consisting of 3 infinite
    buses.

    case file can be used together with dc case files "case5_stagg_....py"
    """
    ppc = {"version": '2'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## bus data
    # bus_i     type    Pd      Qd      Gs  Bs  area    Vm      Va baseKV zone      Vmax    Vmin
    ppc["bus"] = array([
        [2,      inf,      0 ,    0 ,     0,   0,  1,      1.06,   0,	345,    1,      1.1,    0.9],
        [3,      inf,      0 ,    0 ,     0,   0,  1,      1   ,   0,	345,    2,      1.1,    0.9],
        [5,      inf,      0 ,    0 ,     0,   0,  1,      1   ,   0,	345,    3,      1.1,    0.9]
    ])

    ## generator data
    #   bus Pg,Qg, Qmax, Qmin, Vg,  mBase, status, Pmax,   Pmin
    ppc["gen"] = array(
        zeros((0,10))
    )

    ## branch data
    # fbus, tbus,   r,      x,      b,    rateA, rateB, rateC, ratio, angle, status
    ppc["branch"] = array(
        zeros((0,11))
    )

    return ppc