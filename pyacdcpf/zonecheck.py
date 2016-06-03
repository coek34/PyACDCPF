"""Check for non-synchronized ac zones.
"""
from sys import stdout, stderr

from numpy import sort, unique, where, inf, setdiff1d, roll, intersect1d, c_, r_

from pypower.idx_bus import ZONE, BUS_I, BUS_TYPE, REF
from pypower.idx_gen import GEN_BUS
from pypower.idx_brch import F_BUS, T_BUS

def zonecheck(bus, gen, branch, i2eac, output):
    """
    Check for non-synchronized ac zones.
    
    Check for non-synchronized ac zones. If present, the presence of one
    and only one slack bus is checked, as well as the presence of ac
    connections between the zones, in wich case an error message is
    displayed.
    
    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    aczones =  sort(unique(bus[:,ZONE]))

    if aczones.shape[0] > 1 :
        if output:
            stdout.write('Non-synchronised zones: %d AC zones detected.'%aczones.size)

    ## check interzonal connections
    branchzone = c_[bus[[where(bus[:,BUS_I]==x)[0][0] for x in branch[:,F_BUS]],ZONE],\
                    bus[[where(bus[:,BUS_I]==x)[0][0] for x in branch[:,T_BUS]],ZONE]]
    branchconfl = branchzone[:,0] != branchzone[:,1]

    if branchconfl.sum() != 0:
        brconfl = c_[branch[branchconfl,F_BUS].astype(int), branch[branchconfl,T_BUS].astype(int)]
        stdout.write('\nRemove branch between buses %d and %d.\n', i2eac[brconfl,0], i2eac[brconfl,1])
        stderr.write('Connection between different AC zones detected.\n')

    ## check for one ac slack bus (with a generator) in every zone
    ##  addition to original matpower code to avoid that MATPOWER takes a
    ##  voltage controlling converter (dummy generator) and uses it as an ac
    ##  slack bus when there is no generator present.
    acslack = where((bus[:,BUS_TYPE] == REF))[0]
    acslackz = sort(bus[acslack,ZONE])
    acinf = where((bus[:,BUS_TYPE] == inf))[0]
    acinfz = sort(bus[acinf,ZONE])
    nogenslack = setdiff1d(bus[acslack,BUS_I], gen[:,GEN_BUS])
    if acslackz.shape[0] != aczones.shape[0] :
        ## ac zones without slack bus
        if acslackz.shape[0] < aczones.shape[0] :
            zone_noslack = setdiff1d(aczones, acslackz)
            zone_noslack_noinf = setdiff1d(zone_noslack, acinfz)
            if not zone_noslack_noinf.size == 0:
                stdout.write('\nNo AC slack bus detected in AC zone %d.\n'% zone_noslack_noinf)
                stderr.write('Define an AC slack bus for every non-synchronized zone!\n')
        ## multiple slack buses in ac zone
        elif acslackz.shape[0] > aczones.shape[0] :
            ## could also be 2 in one zone
            acslackzs = sort(acslackz)
            multslack = acslackzs[(acslackzs == roll(acslackzs, 1))]
            stdout.write('\nMultiple AC slack bus detected in AC zone %d.\n'% multslack)
            stderr.write('Reduce number of AC slack buses to 1 for every non-synchronized zone!')
    else:
        ## length(acslackz) == length(aczones) %% check for multiple ac slack buses in one zone and inf bus zones
        if acinf:
            acslackzs = sort(acslackz)
            multslack = acslackzs[(acslackzs == roll(acslackzs, 1))]
            if multslack.shape[0] > 1:
                stdout.write('\nMultiple AC slack bus detected in AC zone %d.\n'% multslack)
                stderr.write('Reduce number of AC slack buses to 1 for every non-synchronized zone!')

    ## check  for a generator for every slack bus
    if not nogenslack.size == 0:
        stdout.write('\nAC slack bus without generator at bus %d.\n', i2eac[nogenslack])
        stderr.write('Add a generator for every AC slack bus!')


    ## check for connections to infinite buses
    fbusinf = intersect1d(branch[:,F_BUS], acinf)
    tbusinf = intersect1d(branch[:,T_BUS], acinf)
    brchinf = r_[fbusinf, tbusinf]
    if not brchinf.size == 0:
        stdout.write('\n Connection with an infinite bus at bus %d.\n', i2eac[brchinf])
        stderr.write('Remove connections to infinite buses!')


    ## check for one infinite bus in every zone (without other buses)
    busi = where(bus[:,BUS_I])[0]
    ## bus number of dc busses
    acninf = setdiff1d(busi, acinf).astype(int)
    acninfz = sort(bus[acninf,ZONE])
    conflinfz = intersect1d(acinfz, acninfz)

    if conflinfz:
        stdout.write('\n Infinite buses and regular buses detected in zone %d.\n', conflinfz)
        stderr.write('Remove infinite or regular buses!')


    ## check for multiple infinute buses in 1 AC zone
    multinf = acinfz[(acinfz == roll(acinfz, 1))]
    if multinf.shape[0] > 1:
        stdout.write('\nMultiple infinite buses detected in AC zone %d.\n', np.unique(multinf))
        stderr.write('Reduce number of infinite buses to 1 for every non-synchronized zone!')

    return
