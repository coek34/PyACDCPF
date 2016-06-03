"""Power flow data for system based on the IEEE RELIABILITY TEST SYSTEM.
"""

from numpy import array

def case24_ieee_rts1996_3zones():
    """Power flow data for system based on the IEEE RELIABILITY TEST SYSTEM.
    This system data is based on the MATPOWER case file CASE24_IEEE_RTS
    which is based on the IEEE RELIABILITY TEST SYSTEM

    The data has been adopted to corresponding with the
    IEEE Two Area RTS-96 data from...
    IEEE Reliability Test System Task Force of Applications of
    Probability Methods Subcommittee, "IEEE reliability test system-96,"
    IEEE Transactions on Power Systems, Vol. 14, No. 3, Aug. 1999,
    pp. 1010-1020.

    The IEEE Two Area RTS-96 network has been extended and now includes 3
    asynchronous zones (node numbers 1xx, 2xx and 3xx).
    Data on zone 1 and 2 taken from the IEEE Two Area RTS-96 with following
    adaptations:
    - nodes renumbered according to IEEE Two Area RTS-96 data
    - gen U100 at node 107 disabled (commented)
    - gen U76 at node 201 disabled (commented)
    - slack node zone 2: node 213
    - lines 107-203, 113-215, 123-217 removed (commented)
    Data on zone 3 added:
    - nodes 301 and 302
    - gen at node 302
    - line 301-302
    """
    ppc = {"version": '2'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## bus data
    # bus_i     type    Pd      Qd      Gs  Bs  area    Vm      Va baseKV zone      Vmax    Vmin
    ppc["bus"] = array([
        [101,     2  ,    108	,22,0	,0    , 11,     1   ,   0,  138,    1, 1.05,0.95],
        [102,     2  ,    97   , 20,0	,0    , 11,     1   ,   0,  138,    1, 1.05,0.95],
        [103,     1  ,    180	,37,0	,0    , 11,     1   ,   0,  138,    1, 1.05,0.95],
        [104,     1  ,    74   , 15,0	,0    , 11,     1   ,   0,  138,    1, 1.05,0.95],
        [105,     1  ,    71   , 14,0	,0    , 11,     1   ,   0,  138,    1, 1.05,0.95],
        [106,     1  ,    136	,28,0   ,1.00 , 12,     1   ,   0,  138,    1, 1.05,0.95],
        [107,     2  ,    125	,25,0	,0    , 12,     1   ,   0,  138,    1, 1.05,0.95],
        [108,     1  ,    171	,35,0	,0    , 12,     1   ,   0,  138,    1, 1.05,0.95],
        [109,     1  ,    175	,36,0	,0    , 11,     1   ,   0,  138,    1, 1.05,0.95],
        [110,     1  ,    195	,40,0	,0    , 12,     1   ,   0,  138,    1, 1.05,0.95],
        [111,     1  ,    0    , 0	,0	,0    , 13,     1   ,   0,  230,    1, 1.05,0.95],
        [112,     1  ,    0    , 0	,0	,0    , 13,     1   ,   0,  230,    1, 1.05,0.95],
        [113,     3  ,    265	,54,0	,0    , 13,     1   ,   0,  230,    1, 1.05,0.95],
        [114,     2  ,    194	,39,0	,0    , 13,     1   ,   0,  230,    1, 1.05,0.95],
        [115,     2  ,    317	,64,0	,0    , 14,     1   ,   0,  230,    1, 1.05,0.95],
        [116, 	  2  ,   100	,20,0	,0    , 14,     1   ,   0,  230,    1, 1.05,0.95],
        [117,     1  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    1, 1.05,0.95],
        [118,     2  ,    333	,68,0	,0    , 14,     1   ,   0,  230,    1, 1.05,0.95],
        [119,     1  ,    181	,37,0	,0    , 13,     1   ,   0,  230,    1, 1.05,0.95],
        [120,     1  ,    128	,26,0	,0    , 13,     1   ,   0,  230,    1, 1.05,0.95],
        [121,     2  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    1, 1.05,0.95],
        [122,     2  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    1, 1.05,0.95],
        [123,     2  ,    0    , 0	,0	,0    , 13,     1   ,   0,  230,    1, 1.05,0.95],
        [124,     1  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    1, 1.05,0.95],
        [201,     2  ,    108	,22,0	,0    , 11,     1   ,   0,  138,    2, 1.05,0.95],
        [202,     2  ,    97   , 20,0	,0    , 11,     1   ,   0,  138,    2, 1.05,0.95],
        [203,     1  ,    180	,37,0	,0    , 11,     1   ,   0,  138,    2, 1.05,0.95],
        [204,     1  ,    74   , 15,0	,0    , 11,     1   ,   0,  138,    2, 1.05,0.95],
        [205,     1  ,    71   , 14,0	,0    , 11,     1   ,   0,  138,    2, 1.05,0.95],
        [206,     1  ,    136	,28,0   ,1.00 , 12,     1   ,   0,  138,    2, 1.05,0.95],
        [207,     2  ,    125	,25,0	,0    , 12,     1   ,   0,  138,    2, 1.05,0.95],
        [208,     1  ,    171	,35,0	,0    , 12,     1   ,   0,  138,    2, 1.05,0.95],
        [209,     1  ,    175	,36,0	,0    , 11,     1   ,   0,  138,    2, 1.05,0.95],
        [210,     1  ,    195	,40,0	,0    , 12,     1   ,   0,  138,    2, 1.05,0.95],
        [211,     1  ,    0    , 0	,0	,0    , 13,     1   ,   0,  230,    2, 1.05,0.95],
        [212,     1  ,    0    , 0	,0	,0    , 13,     1   ,   0,  230,    2, 1.05,0.95],
        [213,     3  ,    265	,54,0	,0    , 13,     1   ,   0,  230,    2, 1.05,0.95],
        [214,     2  ,    194	,39,0	,0    , 13,     1   ,   0,  230,    2, 1.05,0.95],
        [215,     2  ,    317	,64,0	,0    , 14,     1   ,   0,  230,    2, 1.05,0.95],
        [216, 	  2  ,   100	,20,0	,0    , 14,     1   ,   0,  230,    2, 1.05,0.95],
        [217,     1  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    2, 1.05,0.95],
        [218,     2  ,    333	,68,0	,0    , 14,     1   ,   0,  230,    2, 1.05,0.95],
        [219,     1  ,    181	,37,0	,0    , 13,     1   ,   0,  230,    2, 1.05,0.95],
        [220,     1  ,    128	,26,0	,0    , 13,     1   ,   0,  230,    2, 1.05,0.95],
        [221,     2  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    2, 1.05,0.95],
        [222,     2  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    2, 1.05,0.95],
        [223,     2  ,    0    , 0	,0	,0    , 13,     1   ,   0,  230,    2, 1.05,0.95],
        [224,     1  ,    0    , 0	,0	,0    , 14,     1   ,   0,  230,    2, 1.05,0.95],
        [301,     1  ,    0    , 0	,0	,0    , 14,     1   ,	0,  230,    3, 1.05,0.95],
        [302,     3  ,    0    , 0	,0	,0    , 14,     1   ,	0,  230,    3, 1.05,0.95]
    ])

    ## generator data
    #   bus Pg,Qg, Qmax, Qmin, Vg,  mBase, status, Pmax,   Pmin
    ppc["gen"] = array([
        [101,    10  ,  0    ,  10 ,    0  ,1.035	,100,   1,  20	,16	],
        [101,    10  ,  0    ,  10 ,    0  ,1.035	,100,   1,  20	,16	],
        [101,    76  ,  14.1 ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        [101,    76  ,  14.1 ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        [102,    10  ,  0    ,  10 ,    0	,1.035	,100,   1,  20	,16	],
        [102,    10  ,  0    ,  10 ,    0	,1.035	,100,   1,  20	,16	],
        [102,    76  ,  7.0  ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        [102,    76  ,  7.0  ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        [107,    80  ,  17.2 ,  60 ,    0	,1.025	,100,   1,  100,25	],
        [107,    80  ,  17.2 ,  60 ,    0	,1.025	,100,   1,  100,25	],
        # [107,    80  ,  17.2 ,  60 ,    0	,1.025	,100,   1,  100,25	],
        [113,    95.1,  40.7 ,  80 ,    0	,1.020	,100,   1,  197,69	],
        [113,    95.1,  40.7 ,  80 ,    0	,1.020	,100,   1,  197,69	],
        [113,    95.1,  40.7 ,  80 ,    0	,1.020	,100,   1,  197,69	],
        [114,    0   ,  13.7 ,  200,    -50,0.980	,100,   1,  0	,0	],
        [115,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [115,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [115,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [115,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [115,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [115,    155 ,  0.05 ,  80 ,    -50,1.014	,100,   1,  155,54.3],
        [116,    155 ,  25.22,  80 ,    -50,1.017	,100,   1,  155,54.3],
        [118,    400 ,  137.4,  200,	-50,1.05	,100,   1,  400,100	],
        [121,    400 ,  108.2,  200,    -50,1.05	,100,   1,  400,100	],
        [122,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [122,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [122,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [122,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [122,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [122,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [123,    155 ,  31.79,  80 ,    -50,1.05	,100,   1,  155,54.3],
        [123,    155 ,  31.79,  80 ,    -50,1.05	,100,   1,  155,54.3],
        [123,    350 ,  71.78,  150,    -25,1.05	,100,   1,  350,140	],
        [201,    10  ,  0    ,  10 ,    0  ,1.035	,100,   1,  20	,16	],
        [201,    10  ,  0    ,  10 ,    0  ,1.035	,100,   1,  20	,16	],
        [201,    76  ,  14.1 ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        # [201,    76  ,  14.1 ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        [202,    10  ,  0    ,  10 ,    0	,1.035	,100,   1,  20	,16	],
        [202,    10  ,  0    ,  10 ,    0	,1.035	,100,   1,  20	,16	],
        [202,    76  ,  7.0  ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        [202,    76  ,  7.0  ,  30 ,    -25,1.035	,100,   1,  76	,15.2],
        [207,    80  ,  17.2 ,  60 ,    0	,1.025	,100,   1,  100,25	],
        [207,    80  ,  17.2 ,  60 ,    0	,1.025	,100,   1,  100,25	],
        [207,    80  ,  17.2 ,  60 ,    0	,1.025	,100,   1,  100,25	],
        [213,    95.1,  40.7 ,  80 ,    0	,1.020	,100,   1,  197,69	],
        [213,    95.1,  40.7 ,  80 ,    0	,1.020	,100,   1,  197,69	],
        [213,    95.1,  40.7 ,  80 ,    0	,1.020	,100,   1,  197,69	],
        [214,    0   ,  13.68,  200,    -50,0.980	,100,   1,  0	,0	],
        [215,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [215,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [215,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [215,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [215,    12  ,  0    ,  6  ,    0	,1.014	,100,   1,  12	,2.4	],
        [215,    155 ,  0.048,  80 ,    -50,1.014	,100,   1,  155,54.3],
        [216,    155 ,  25.22,  80 ,    -50,1.017	,100,   1,  155,54.3],
        [218,    400 ,  137.4,  200,	-50,1.05	,100,   1,  400,100	],
        [221,    400 ,  108.2,  200,    -50,1.05	,100,   1,  400,100	],
        [222,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [222,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [222,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [222,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [222,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [222,    50  ,  -4.96,  16 ,    -10,1.05	,100,   1,  50	,10	],
        [223,    155 ,  31.79,  80 ,    -50,1.05	,100,   1,  155,54.3],
        [223,    155 ,  31.79,  80 ,    -50,1.05	,100,   1,  155,54.3],
        [223,    350 ,  71.78,  150,    -25,1.05	,100,   1,  350,140	],
        [302,    150 ,     10,  150,    -25,1.05	,100,   1,  350,140	]
    ])

    ## branch data
    # fbus, tbus,   r,      x,      b,    rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["branch"] = array([
        [101,	102,    0.003,  0.014,  0.461,  175,    250,    200,    0     , 0,	1,  -360,   360],
        [101,	103,    0.055,  0.211,  0.057,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [101,	105,    0.022,  0.085,  0.023,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [102,	104,    0.033,  0.127,  0.034,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [102,	106,    0.050,  0.192,  0.052,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [103,	109,    0.031,  0.119,  0.032,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [103,	124,    0.002,  0.084,  0    ,  400,    510,    600,    1.015 , 0,	1,  -360,   360],
        [104,	109,    0.027,  0.104,  0.028,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [105,	110,    0.022,  0.088,  0.024,  175,    208,    220,	0     , 0,	1,  -360,   360],
        [106,	110,    0.014,  0.061,  2.459,  175,    193,    200,    0     , 0,	1,  -360,   360],
        [107,	108,    0.016,  0.061,  0.017,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [108,	109,    0.043,  0.165,  0.045,  175,    208,    220,	0     , 0,	1,  -360,   360],
        [108,	110,    0.043,  0.165,  0.045,  175,    208,    220,	0     , 0,	1,  -360,   360],
        [109,	111,    0.002,  0.084,  0    ,  400,    510,    600,    1.03  , 0,	1,  -360,   360],
        [109,	112,    0.002,  0.084,  0    ,  400,    510,    600,    1.03  , 0,	1,  -360,   360],
        [110,	111,    0.002,  0.084,  0    ,  400,    510,    600,    1.015 , 0,	1,  -360,   360],
        [110,	112,    0.002,  0.084,  0    ,  400,    510,    600,    1.015 , 0,	1,  -360,   360],
        [111,	113,    0.006,  0.048,  0.100,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [111,	114,    0.005,  0.042,  0.088,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [112,	113,    0.006,  0.048,  0.100,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [112,	123,    0.012,  0.097,  0.203,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [113,	123,    0.011,  0.087,  0.182,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [114,	116,    0.005,  0.059,  0.082,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [115,	116,    0.002,  0.017,  0.036,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [115,	121,	0.006,  0.049,  0.103,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [115,	121,	0.006,  0.049,  0.103,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [115,	124,    0.007,  0.052,  0.109,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [116,	117,    0.003,  0.026,  0.055,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [116,	119,    0.003,  0.023,  0.049,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [117,	118,    0.002,  0.014,  0.030,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [117,	122,    0.014,  0.105,  0.221,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [118,	121,    0.003,  0.026,  0.055,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [118,	121,    0.003,  0.026,  0.055,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [119,	120,    0.005,  0.040,  0.083,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [119,	120,    0.005,  0.040,  0.083,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [120,	123,    0.003,  0.022,  0.046,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [120,	123,    0.003,  0.022,  0.046,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [121,	122,    0.009,  0.068,  0.142,  500,    600,    625,	0     , 0,	1,  -360,   360],
        [201,	202,    0.003,  0.014,  0.461,  175,    250,    200,    0     , 0,	1,  -360,   360],
        [201,	203,    0.055,  0.211,  0.057,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [201,	205,    0.022,  0.085,  0.023,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [202,	204,    0.033,  0.127,  0.034,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [202,	206,    0.050,  0.192,  0.052,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [203,	209,    0.031,  0.119,  0.032,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [203,	224,    0.002,  0.084,  0    ,  400,    510,    600,    1.015 , 0,	1,  -360,   360],
        [204,	209,    0.027,  0.104,  0.028,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [205,	210,    0.022,  0.088,  0.024,  175,    208,    220,	0     , 0,	1,  -360,   360],
        [206,	210,    0.014,  0.061,  2.459,  175,    193,    200,    0     , 0,	1,  -360,   360],
        [207,	208,    0.016,  0.061,  0.017,  175,    208,    220,    0     , 0,	1,  -360,   360],
        [208,	209,    0.043,  0.165,  0.045,  175,    208,    220,	0     , 0,	1,  -360,   360],
        [208,	210,    0.043,  0.165,  0.045,  175,    208,    220,	0     , 0,	1,  -360,   360],
        [209,	211,    0.002,  0.084,  0    ,  400,    510,    600,    1.03  , 0,	1,  -360,   360],
        [209,	212,    0.002,  0.084,  0    ,  400,    510,    600,    1.03  , 0,	1,  -360,   360],
        [210,	211,    0.002,  0.084,  0    ,  400,    510,    600,    1.015 , 0,	1,  -360,   360],
        [210,	212,    0.002,  0.084,  0    ,  400,    510,    600,    1.015 , 0,	1,  -360,   360],
        [211,	213,    0.006,  0.048,  0.100,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [211,	214,    0.005,  0.042,  0.088,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [212,	213,    0.006,  0.048,  0.100,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [212,	223,    0.012,  0.097,  0.203,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [213,	223,    0.011,  0.087,  0.182,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [214,	216,    0.005,  0.059,  0.082,  500,    625,    625,    0     , 0,	1,  -360,   360],
        [215,	216,    0.002,  0.017,  0.036,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [215,	221,	0.006,  0.049,  0.103,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [215,	221,	0.006,  0.049,  0.103,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [215,	224,    0.007,  0.052,  0.109,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [216,	217,    0.003,  0.026,  0.055,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [216,	219,    0.003,  0.023,  0.049,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [217,	218,    0.002,  0.014,  0.030,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [217,	222,    0.014,  0.105,  0.221,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [218,	221,    0.003,  0.026,  0.055,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [218,	221,    0.003,  0.026,  0.055,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [219,	220,    0.005,  0.040,  0.083,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [219,	220,    0.005,  0.040,  0.083,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [220,	223,    0.003,  0.022,  0.046,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [220,	223,    0.003,  0.022,  0.046,  500,    600,    625,    0     , 0,	1,  -360,   360],
        [221,	222,    0.009,  0.068,  0.142,  500,    600,    625,	0     , 0,	1,  -360,   360],
        # [107,	203,    0.042,  0.161,  0.044,  175,    208,    220,	0     , 0,	1,  -360,   360],
        # [113,	215,    0.010,  0.075,  0.158,  500,    600,    625,	0     , 0,	1,  -360,   360],
        # [123,	217,    0.010,  0.074,  0.155,  500,    600,    625,	0     , 0,	1,  -360,   360]
        [301,	302,    0.000,  0.001,  0.000,  500,    600,    625,	0     , 0,	1,  -360,   360]
    ])

    return ppc