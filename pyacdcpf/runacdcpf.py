"""RUNACDCPF  Runs a sequential ac/dc power flow.
"""

from sys import stdout, stderr

from os.path import dirname, join

from time import time

from numpy import r_, c_,  zeros, pi, exp, where, inf, equal, not_equal,\
                 setdiff1d, arange, intersect1d, union1d, sqrt, sort,\
                 unique, real, imag, conj, delete, complex, abs, angle, \
                 argmax

from pypower.loadcase import loadcase
from pypower.ppoption import ppoption
from pypower.runpf import runpf

from pypower.idx_bus import PD, QD, VM, VA, BUS_TYPE, PQ, REF, PV, ZONE, \
                BUS_I
from pypower.idx_brch import QT, F_BUS
from pypower.idx_gen import PG, QG, VG, QMAX, QMIN, GEN_BUS, GEN_STATUS, \
                MBASE, PMAX, PMIN

from pyacdcpf.idx_busdc import BUSAC_I, GRIDDC, BUSDC_I, VDC, PDC
from pyacdcpf.idx_brchdc import F_BUSDC, T_BUSDC
from pyacdcpf.idx_convdc import CONV_BUS, CONVSTATUS, CONVTYPE_DC, DCSLACK,\
        DCDROOP, DCNOSLACK, CONVTYPE_AC, PVC, PQC, PCONV, QCONV, VCONV, \
        DROOP, PDCSET, VDCSET, DVDCSET, BASEKVC, LOSSA, LOSSB, LOSSCR, \
        LOSSCI, ICMAX, VCMAX, VCMIN, RCONV, XCONV, RTF, XTF, BF

from pyacdcpf.pacdcoption import pacdcoption
from pyacdcpf.loadcasedc import loadcasedc
from pyacdcpf.convout import convout
from pyacdcpf.brchdcout import brchdcout
from pyacdcpf.brchout import brchout
from pyacdcpf.ext2intdc import ext2intdc
from pyacdcpf.ext2intac import ext2intac
from pyacdcpf.ext2intpu import ext2intpu
from pyacdcpf.int2extdc import int2extdc
from pyacdcpf.int2extac import int2extac
from pyacdcpf.int2extpu import int2extpu
from pyacdcpf.makeYbusdc import makeYbusdc
from pyacdcpf.zonecheck import zonecheck
from pyacdcpf.convlim import convlim
from pyacdcpf.convlimplot import convlimplot
from pyacdcpf.calclossac import calclossac
from pyacdcpf.dcnetworkpf import dcnetworkpf
from pyacdcpf.calcslackdroop import calcslackdroop
from pyacdcpf.printdcpf import printdcpf
from pyacdcpf.printpf import printpf # Small adaptation was made in order to support inf network


## define j
## DONT USE j IN ANYWHERE ELSE!!!
j = sqrt(-1+0j)


def runacdcpf(caseac=None, casedc=None, pacdcopt=None, ppopt=None):
    """
	Runs a sequential AC/DC power flow, optionally
	returning the results, a convergence flag and the time.

    Inputs (optional):
		CASEAC : ac power flow data
			either a PYPOWER case struct or a string containing the name
			of the file with the data (default ac case is 'case5_stagg',
			only used when both ac and dc power flow data are not defined)
			(see also CASEFORMAT and LOADCASE and PYPOWER)
		CASEDC : dc power flow data
			either a PYACDCPF case struct or a string containing
			the name of the file with the data (default dc case is
			'case5_stagg_MTDCslack')
			(see also LOADCASEDC)
		PACDCOPT : PYACDCPF options vector to override default ac/dc power
			flow options. Can be used to specify tolerances, inclusion of
			limits, plot options and more (see also MACDCOPTION).
		PPOPT : PYPOWER options vector to override default options
			can be used to specify the solution algorithm, output options
			termination tolerances, and more (see also MPOPTION).

	Outputs:
		RESULTSAC : results struct, with the following fields from the
		input PYPOWER case: baseMVA, bus, branch, gen  (but with solved
		voltages, power flows, etc.)
		RESULTSDC : results struct, with the following fields:
		input PYACDCPF dc case: baseMVAac, baseMVAdc, pol, busdc, convdc,
		branchdc (but with solved voltages, power flows, etc.)
		CONVERGED : converge flag, can additionally be returned
		TIMECALC : elapsed time, can additionally be returned

	Examples of usage:
		[resultsac, resultsdc, converged, te] = runacdcpf('case5_stagg', \
		'case5_stagg_MTDCdroop');

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    ## start of time calculation
    t0 = time()

    ## add subdirectories to path
    if caseac is None:
        dirac = join(dirname(__file__), 'Cases', 'PowerflowAC')
        caseac = join(dirac, 'case5_stagg')
        # caseac = join(dirac, 'case3_inf')
        # caseac = join(dirac, 'case24_ieee_rts1996_3zones')
        # caseac = join(dirac, 'case24_ieee_rts1996_3zones_inf')
    if casedc is None:
        dirdc = join(dirname(__file__), 'Cases', 'PowerflowDC')
        # casedc = join(dirdc, 'case5_stagg_HVDCptp')
        # casedc = join(dirdc, 'case5_stagg_MTDCslack')
        casedc = join(dirdc, 'case5_stagg_MTDCdroop')
        # casedc = join(dirdc, 'case24_ieee_rts1996_MTDC')

    ## default arguments
    ppopt = ppoption(ppopt)
    ppopt["VERBOSE"] = 0
    ppopt["OUT_ALL"] = 0
    pacdcopt = pacdcoption(pacdcopt)

    ## options
    tolacdc = pacdcopt["TOLACDC"]
    itmaxacdc = pacdcopt["ITMAXACDC"]
    toldc = pacdcopt["TOLDC"]
    itmaxdc = pacdcopt["ITMAXDC"]
    tolslackdroop = pacdcopt["TOLSLACKDROOP"]
    itmaxslackdroop = pacdcopt["ITMAXSLACKDROOP"]
    tolslackdroopint = pacdcopt["TOLSLACKDROOPINT"]
    itmaxslackdroopint = pacdcopt["ITMAXSLACKDROOPINT"]
    multslack = pacdcopt["MULTSLACK"]

    limac = pacdcopt["LIMAC"]
    limdc = pacdcopt["LIMDC"]
    tollim = pacdcopt["TOLLIM"]

    output = pacdcopt["OUTPUT"]
    convplotopt = pacdcopt["CONVPLOTOPT"]

    ## -----  initialise  -----
    ## read data
    pdc = loadcasedc(casedc)
    ppc = loadcase(caseac)

    ##-----  Data preparation -----
    ## converter outage are considered as stations without ac grid connection
    pdc, conv0busi, conv1, conv1i, conv0, conv0i = convout(pdc)
    pdc['convdc'] = conv1 # only use converters without outage

    ## dc branch outages (remove branches from input data)
    brchdc1, brchdc1i, brchdc0, brchdc0i = brchdcout(pdc)
    pdc['branchdc'] = brchdc1 # only include branches in operation

    ## ac branch outages (remove branches from input data)
    brch1, brch1i, brch0, brch0i = brchout(ppc)
    ppc['branch'] = brch1 # only include branches in operation

    ## generator outages (remove non-operational generators from input data)
    gon = where(ppc['gen'][:,GEN_STATUS] > 0)[0]   # which gens are on?
    goff = where(ppc['gen'][:,GEN_STATUS] == 0)[0] # which gens are off?
    gen0 = ppc['gen'][goff,:] # non-operational gens data
    ppc['gen'] = ppc['gen'][gon,:] #keep operational generators"

    ##-----  External to internal numbering  -----
    ## dc network external to internal bus numbering
    i2edcpmt, i2edc, pdc = ext2intdc(pdc)

    ## ac network external to internal bus numbering
    acdmbus, i2eac, pdc, ppc = ext2intac(pdc,ppc)

    ## sort matrices by new bus numbers
    i2ebus    = ppc['bus'][:,0].argsort()
    i2egen    = ppc['gen'][:,0].argsort()
    i2ebrch   = ppc['branch'][:,0].argsort()
    i2ebusdc  = pdc['busdc'][:,0].argsort()
    i2econvdc = pdc['convdc'][:,0].argsort()
    i2ebrchdc = pdc['branchdc'][:,0].argsort()

    ppc['bus'] = ppc['bus'][i2ebus,:]
    ppc['gen'] = ppc['gen'][i2egen,:]
    ppc['branch'] = ppc['branch'][i2ebrch,:]
    pdc['busdc'] = pdc['busdc'][i2ebusdc,:]
    pdc['convdc'] = pdc['convdc'][i2econvdc,:]
    pdc['branchdc'] = pdc['branchdc'][i2ebrchdc,:]

    ## Per unit external to internal data conversion
    pdc = ext2intpu(ppc['baseMVA'],pdc)

    ##-----  Additional data preparation & index initialisation  -----
    ## zero rows addition to convdc matrix (dc buses without converter)
    convdc1 = zeros((pdc['busdc'].shape[0]-pdc['convdc'].shape[0],
                     pdc['convdc'].shape[1]))
    pdc['convdc'] = r_[pdc['convdc'],convdc1]

    ## indices initialisation
    bdci = where(pdc['busdc'][:,BUSAC_I])[0].astype(int)
    cdci = where(pdc['convdc'][:,CONVSTATUS] == 1)[0].astype(int)
    slackdc = where(pdc['convdc'][:,CONVTYPE_DC] == DCSLACK)[0].astype(int)

    droopdc = where(pdc['convdc'][:,CONVTYPE_DC] == DCDROOP)[0].astype(int)
    ngriddc = pdc['busdc'][:,GRIDDC].max()

    ## convert to internal indexing
    baseMVA, bus, gen, branch = \
        ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"]

    baseMVAac, baseMVAdc, pol, busdc, convdc, branchdc = \
        pdc["baseMVAac"], pdc["baseMVAdc"], pdc["pol"], \
        pdc["busdc"], pdc["convdc"], pdc["branchdc"]

    ##-----  Violation check  -----
    ## dc slack bus and distributed voltage bus violation check
    gridviol = setdiff1d(arange(1,ngriddc+1),busdc[r_[slackdc, droopdc],GRIDDC])

    if not gridviol.size == 0:
        stdout.write('\nMultiple dc slack busses defined in grid %d \n' % (gridviol))
        stderr.write('No droop controlled bus or slack bus defined for every dc grid !\n')

    ## remove multiple slack buses
    if multslack == 0:
        for ii in arange(1,ngriddc+1):
            slackdcii = intersect1d(slackdc,where(busdc[:,GRIDDC] == ii)[0])
            if slackdcii.size > 1:
                convdc[slackdcii[0].astype(int),CONVTYPE_DC ] = DCSLACK
                convdc[slackdcii[1:].astype(int),CONVTYPE_DC ] = DCNOSLACK
                slackdcii = slackdcii[0]

                ##  printout changes
                stdout.write('\nMultiple dc slack busses defined in grid %d' %(ii))
                stdout.write('\n     Bus %d kept as the slack bus\n'%\
                            (i2edc[slackdcii.astype(int)+1]))

        ## redefine slack buses
        slackdc  = where(convdc[:,CONVTYPE_DC] == DCSLACK)[0].astype(int)

    ## define indices of slack, droop and power controlled buses
    slackdroopdc = union1d(slackdc, droopdc)
    noslackbdc   = setdiff1d(where(busdc[:,BUSDC_I]), slackdc)

    ## remove converter and generator V control violations
    vcontrvsc = where(convdc[:,CONVTYPE_AC] == PVC)[0]
    vcontrgen = r_[where(bus[:,BUS_TYPE] == PV)[0], \
                   where(bus[:,BUS_TYPE] == REF)[0]]
    ##  buses with V control conflicts
    vconfl = intersect1d(vcontrvsc, vcontrgen).astype(int)
    convdc[vconfl,CONVTYPE_AC] = PQC
    convdc[vconfl,QCONV] = 0
    convdc[:,QCONV] *= convdc[:,CONVTYPE_AC] == PQC
    if not vconfl.size == 0:
        stdout.write('Generator & VSC converter on the same bus')
        stdout.write('\n   Conflicting voltage control on bus %d' %(i2eac[vconfl+1]))
        stdout.write('\n=> Corresponding VSC Converter set to PQ control without Q injections.\n')

    ##-----  initialisation ac network  -----
    ## dummy generator initialisation
    Vcref = convdc[:,VCONV] # voltage setpoints
    busVSC = bus.copy()
    gendm = zeros((0,gen.shape[1]))
    genPQ = zeros((0,1)).astype(int)
    genPQi = zeros((0,1)).astype(int)
    Qcmin_dum = -99999
    Qcmax_dum   =  99999
    Pcmin_dum   =      0
    Pcmax_dum   =  99999

    ## dummy generator addition
    for ii in arange(convdc.shape[0]):
        ## change control from PQ to PV for buses with converter in PV control
        if bus[ii,BUS_TYPE] == PQ and convdc[ii,CONVTYPE_AC] == PVC:
            busVSC[ii,BUS_TYPE] = PV
            ## add dummy generator to V controlling converter bus without generator
            if not any(gen[:,GEN_BUS] == bus[ii,BUS_I]):
                gendm = r_[gendm, zeros((1,gen.shape[1]))]
                gendm[-1,[GEN_BUS,PG,QG,QMAX,QMIN,VG,MBASE,GEN_STATUS,PMAX,PMIN]] = \
                    [ii+1,0,0,Qcmax_dum,Qcmin_dum,Vcref[ii],baseMVAac,1,Pcmax_dum,Pcmin_dum]
            else:
                genPQ = r_[genPQ,bus[ii,BUS_I]]
                genPQii = where(gen[:,GEN_BUS] == bus[ii,BUS_I])[0]
                genPQi = r_[genPQi,genPQii]

    ## define buses with dummy generator
    # gdmbus = where(gendm[[where(bus[:,BUS_I]==x)[0][0] for x in \
                # gendm[:,GEN_BUS]],GEN_BUS])[0]
    if any(gendm[:,GEN_BUS]):
        gdmbus = where(gendm[:,GEN_BUS] == bus[:,BUS_I])[0]
    else:
        gdmbus = []

    ## converter stations power injections into ac network
    Pvsc = convdc[:,PCONV]/baseMVA
    Qvsc = convdc[:,QCONV]/baseMVA

    ## dc voltage droop setpoints and parameters
    PVdroop = zeros(busdc.shape[0])
    Pdcset = zeros(busdc.shape[0])
    Vdcset = zeros(busdc.shape[0])
    dVdcset = zeros(busdc.shape[0])

    PVdroop[cdci] = convdc[cdci,DROOP]*baseMVA
    Pdcset[cdci] = convdc[cdci,PDCSET]/baseMVA
    Vdcset[cdci] = convdc[cdci,VDCSET]
    dVdcset[cdci] = convdc[cdci,DVDCSET]

    ## voltage droop converter power initialisation
    Pvsc[droopdc] = Pdcset[droopdc]     ## assumption: operating in reference set-point & no converter losses

    ## dc slack converter power injection initialisation
    if slackdc.size != 0:
        for ii in arange(1,ngriddc+1):
            #slackdcii = where(convdc[busdc[:,GRIDDC]==ii,CONVTYPE_DC]==DCSLACK)[0]
            slackdcii = intersect1d(where(busdc[:,GRIDDC]==ii)[0],\
                            where(convdc[:,CONVTYPE_DC]==DCSLACK)[0])
            if slackdcii.size != 0:
                #Pvscii = Pvsc*(convdc[busdc[:,GRIDDC]==1,CONVTYPE_DC]!=DCSLACK)
                Pvscii = Pvsc*(equal(busdc[:,GRIDDC],ii)* \
                        not_equal(convdc[:,CONVTYPE_DC],DCSLACK))
                Pvsc[slackdcii] = -Pvscii.sum()/slackdcii.shape[0]

    ## Inclusion of converters as loads
    busVSC[cdci,PD] = bus[cdci,PD] - Pvsc[cdci]*baseMVA
    busVSC[cdci,QD] = bus[cdci,QD] - Qvsc[cdci]*baseMVA


    ##-----  initialisation of converter quantities -----
    ## per unit converter loss coefficients values
    basekA = baseMVA/(sqrt(3)*convdc[:,BASEKVC])
    lossa = convdc[:,LOSSA]/baseMVA
    lossb = convdc[:,LOSSB]*basekA/baseMVA
    losscr = convdc[:,LOSSCR]*basekA**2/baseMVA
    lossci = convdc[:,LOSSCI]*basekA**2/baseMVA

    ## converter reactor parameters
    Rc = convdc[:,RCONV]
    Xc = convdc[:,XCONV]
    Zc = Rc+j*Xc

    ## converter limits data
    Icmax = convdc[:,ICMAX]
    Vcmax = convdc[:,VCMAX]
    Vcmin = convdc[:,VCMIN]

    ## filter reactance
    Bf = convdc[:,BF]

    ## transformer parameters
    Rtf = convdc[:,RTF]
    Xtf = convdc[:,XTF]
    Ztf = Rtf+j*Xtf


    ##-----  initialisation of dc network quantities -----
    ## build dc bus matrix
    Ybusdc, Yfdc, Ytdc = makeYbusdc( busdc, branchdc )

    ## detect ac islands errors (non-synchronised zones => to be solved independently)
    zonecheck(bus, gen, branch, i2eac, output)
    aczones = sort(unique(bus[:,ZONE])).astype(int)


    ##-----  main iteration loop -----
    ## initialise
    Vdc = busdc[:,VDC] #dc bus voltages
    genVSC = r_[gen, gendm] #inclusion of dummy generators for ac solution

    gendmidx = where(genVSC[:,GEN_BUS] == setdiff1d(genVSC[:,GEN_BUS], gen[:, GEN_BUS]))[0] # index of dummy generators in genVSC matrix
    Ps = Pvsc #grid side converter power initialisation
    Pdc = zeros(busdc.shape[0])
    Ifdc = zeros(branchdc.shape[0])
    Pfdc = zeros(branchdc.shape[0])
    Ptdc = zeros(branchdc.shape[0])

    ## iteration options
    it = 0
    converged = 0

    ## main loop
    while (not converged) and (it <= itmaxacdc):
        ## update iteration counter
        it += 1

        ## reset grid side converter reactive power injection
        Qs = Qvsc
        Ss = Ps +j*Qs

        ##-----  ac network power flow  -----
        ## ac power flow with converters as loads (PQ mode) or load+generator (PV mode)
        busVSCext = zeros((0,bus.shape[1]))
        genVSCext = zeros((0,gen.shape[1]))
        branchext = zeros((0,QT+1))
        buszi = []
        genVSCzi = []
        brchzi = []
        for i in arange(aczones.size):
            ## select buses, generators and branches in the specified ac zone
            buszi = where(bus[:,ZONE] == aczones[i])[0]
            genVSCzi = where(bus[[where(bus[:,BUS_I]==x)[0][0] for x in \
                genVSC[:,GEN_BUS]],ZONE] == aczones[i])[0]
            brchzi = where(bus[[where(bus[:,BUS_I]==x)[0][0] for x in \
                branch[:,F_BUS]],ZONE] == aczones[i])[0]


            busVSCz = busVSC[buszi,:]
            genVSCz = genVSC[genVSCzi,:]
            branchz = branch[brchzi,:]

            ## solve ac power flow for specified ac zone (if not infinite bus)
            if busVSCz.shape[0]>1:
                accaseVSCz = {}
                accaseVSCz['baseMVA'] = baseMVA
                accaseVSCz['bus'] = busVSCz
                accaseVSCz['gen'] = genVSCz
                accaseVSCz['branch'] = branchz
                resultsz,successz = runpf(accaseVSCz, ppopt)
                busVSCz = resultsz['bus']
                genVSCz = resultsz['gen']
                branchz = resultsz['branch']

            ## store solutions for specified ac zone in extended matrices
            for k,idx in enumerate(buszi):
                if busVSCext.shape[0] <= idx:
                    busVSCext = r_[busVSCext,zeros((idx-busVSCext.shape[0]+1, \
                                    bus.shape[1]))]
                busVSCext[idx,:] = busVSCz[k,:]

            for k,idx in enumerate(genVSCzi):
                if genVSCext.shape[0] <= idx:
                    genVSCext = r_[genVSCext,zeros((idx-genVSCext.shape[0]+1, \
                                    gen.shape[1]))]
                genVSCext[idx,:] = genVSCz[k,:]

            for k,idx in enumerate(brchzi):
                if branchext.shape[0] <= idx:
                    branchext = r_[branchext,zeros((idx-branchext.shape[0]+1, \
                                    QT+1))]
                branchext[idx,:] = branchz[k,:]

        busVSC = busVSCext.copy()
        genVSC = genVSCext.copy()
        branch = branchext.copy()

        ## dummy generator update
        gendm = genVSC[gen.shape[0]:,:]

        ## dummy generator on converter V controlled bus
        Ss[gdmbus] = Ss[gdmbus] + j*gendm[:,QG]/baseMVA

        ## PQ generator on converter V controlled bus
        Ss[genPQ] = Ss[genPQ] + \
        j*(genVSC[genPQi,QG] - gen[genPQi,QG])/baseMVA

        ## update grid side converter power injections
        Ps = real(Ss)
        Qs = imag(Ss)

        ## generator reset
        genVSC[gendmidx,QG] = 0
        genVSC[genPQi,QG] = gen[genPQi,QG]

        ##----- Converter calculations -----
        ## converter reactor voltages and power
        Vs = busVSC[bdci,VM]*exp(j*busVSC[bdci,VA]*pi/180)
        Itf = conj(Ss/Vs)         ## transformer current
        Vf = Vs + Itf*Ztf       ## filter side voltage
        Ssf = Vf*conj(Itf)        ## filter side transformer complex power
        Qf = -Bf*abs(Vf)**2      ## filter reactive power
        Scf = Ssf + j*Qf             ## filter side converter complex power
        Ic = conj(Scf/Vf)        ## converter current
        Vc = Vf + Ic*Zc          ## converter side voltage
        Sc = Vc*conj(Ic)         ## converter side complex power

        ## converter active and reactive powers
        Pc = real(Sc)
        Qc = imag(Sc)
        Pcf = real(Scf)
        Qcf = imag(Scf)
        Psf = real(Ssf)
        Qsf = imag(Ssf)

        ## initialisation
        Ps_old = Ps.copy()

        if limac == 1:
            ##--- converter limit check ---
            ## initialisation
            limviol = zeros((busdc.shape[0]))
            SsL     = zeros((busdc.shape[0]),dtype=complex)
            plotarg = zeros((busdc.shape[0],17),dtype=complex)

            for ii in arange(1,ngriddc+1):
                ## remove slack converters from limit check
                cdcii = where(busdc[:,GRIDDC] == ii)[0]
                ccdcslackii = where(intersect1d(cdcii,slackdc)==cdcii)[0]
                if not ccdcslackii.size == 0:
                    cdcii = delete(cdcii,ccdcslackii) #remove slack converter
                cdci0 = where(convdc[cdcii,CONV_BUS]==0)[0]
                cdcii = delete(cdcii,cdci0) # remove zero elements (converter outages)
                ## converter limit check
                for jj in arange(cdcii.size):
                    cvjj = cdcii[jj]
                    limviol[cvjj],SsL[cvjj], plotarg[cvjj,:] = convlim(Ss[cvjj], \
                        Vs[cvjj], Vc[cvjj], Ztf[cvjj], Bf[cvjj], Zc[cvjj], \
                        Icmax[cvjj], Vcmax[cvjj], Vcmin[cvjj], i2edc[cvjj+1], \
                        tollim, convplotopt)

                ## converter limit violations (1 = Q limit, 2 = P limit)
                limviolii   = limviol*(busdc[:,GRIDDC] == ii)
                dSii  = (SsL-Ss)*(busdc[:,GRIDDC] == ii)*(convdc[:,CONVTYPE_DC] != DCSLACK)
                if (2 in limviolii) or (1 in limviolii):
                    if (2 in limviolii):
                        dSii = dSii*(limviolii==2)
                        dSiimaxi = where(abs(real(dSii)).max())[0]
                        stdout.write('\n  Active power setpoint of converter %d changed from %.2f MW to %.2f MW.'%( \
                            i2edc[dSiimaxi+1], real(Ss[dSiimaxi])*baseMVA, real(SsL[dSiimaxi])*baseMVA))
                        stdout.write('\n  Reactive power setpoint of converter %d changed from %.2f MVAr to %.2f MVAr.\n'%(\
                            i2edc[dSiimaxi+1], imag(Ss[dSiimaxi])*baseMVA, imag(SsL[dSiimaxi])*baseMVA))
                    else: ## if ismember(1, limviolii)
                        dSii = dSii*(limviolii==1)
                        dSiimaxi = argmax(abs(imag(dSii)))
                        stdout.write('\n  Reactive power setpoint of converter %d changed from %.2f MVAr to %.2f MVAr. \n'%(\
                            i2edc[dSiimaxi+1], imag(Ss[dSiimaxi])*baseMVA, imag(SsL[dSiimaxi])*baseMVA))

                    ## plot converter setpoint adaptation
                    if convplotopt != 0 :
                        convlimplot(plotarg[dSiimaxi,:], i2edc[dSiimaxi])

                    ## update converter powers
                    Ss[dSiimaxi] = SsL[dSiimaxi]
                    Pvsc[dSiimaxi] = real(Ss[dSiimaxi])
                    Qvsc[dSiimaxi] = imag(Ss[dSiimaxi])
                    busVSC[dSiimaxi,PD] = bus[dSiimaxi,PD] - \
                        Pvsc[dSiimaxi]*baseMVA  ## converter P injection from input files included as load
                    busVSC[dSiimaxi,QD] = bus[dSiimaxi,QD] - \
                        Qvsc[dSiimaxi]*baseMVA  ## only Q from input files is included, not for V control
                else:
                    dSiimaxi = []

                ## Remove voltage control on violated converter
                if convdc[dSiimaxi, CONVTYPE_AC]==PVC:
                    convdc[dSiimaxi, CONVTYPE_AC] = PQC
                    stdout.write('  Voltage control at converter bus %d removed.\n'% i2edc[dSiimaxi+1])

                    busVSC[dSiimaxi, BUS_TYPE]  = PQ
                    ## Remove dummy generator (PV bus changed to PQ bus)
                    if dSiimaxi in gdmbus:
                        dSidx = where(gdmbus == dSiimaxi)[0]
                        dSgenidx = gendmidx[dSidx]
                        gendm = delete(gendm,dSidx)
                        genVSC = delete(genVSC,dSgenidx)
                        gdmbus = delete(gdmbus,dSidx)
                        gendmidx = where(genVSC[:,GEN_BUS] == np.setdiff1d(genVSC[:, GEN_BUS],gen[:, GEN_BUS])) ## index of dummy generators in genVSC matrix

                    ## Remove VSC voltage control at genPQ bus
                    if dSiimaxi in genPQ:
                        dSidx = where(genPQ == dSiimaxi)[0]
                        genPQ = delete(genPQ,dSidx)
                        genPQi = delete(genPQi,dSidx)

                ## Remove droop control on violated converter
                if convdc[dSiimaxi, CONVTYPE_DC]==DCDROOP:
                   convdc[dSiimaxi, CONVTYPE_DC] = DCNOSLACK
                   droopdc = setdiff1d(droopdc,dSiimaxi) ## remove converter from droop converters
                   slackdroopdc = setdiff1d(slackdroopdc,dSiimaxi) ## remove converter from slack/droop converters (additional loss iteration)
                   stdout.write('  Droop control at converter bus %d disabled.\n'%i2edc[dSiimaxi+1])

            ## recalculate converter quantities after limit check
            Itf = conj(Ss/Vs)         ## transformer current
            Vf = Vs + Itf*Ztf        ## filter side voltage
            Ssf = Vf*conj(Itf)        ## filter side transformer complex power
            Qf = -Bf*abs(Vf)**2      ## filter reactive power
            Scf = Ssf + j*Qf             ## filter side converter complex power
            Ic = conj(Scf/Vf)        ## converter current
            Vc = Vf + Ic*Zc          ## converter side voltage
            Sc = Vc*conj(Ic)         ## converter side complex power

            ## converter active and reactive powers after limit check
            Ps = real(Ss)
            Qs = imag(Ss)
            Pc = real(Sc)
            Qc = imag(Sc)
            Pcf = real(Scf)
            Qcf = imag(Scf)
            Psf = real(Ssf)
            Qsf = imag(Ssf)

        ## converter losses and dc side power
        Ploss = calclossac(Pc, Qc, Vc, lossa, lossb, losscr, lossci)
        Pdc[cdci] = Pc[cdci] + Ploss[cdci]

        ##-----  dc networks power flow  -----
        ## calculate dc networks
        Vdc, Pdc = dcnetworkpf(Ybusdc, Vdc, Pdc,slackdc, noslackbdc,\
            droopdc, PVdroop, Pdcset, Vdcset, dVdcset, pol, toldc, itmaxdc)

        ## calculate dc line powers
        Ifdc = Yfdc*Vdc ## current through dc lines
        Vdcf = Vdc[[where(busdc[:,BUSDC_I]==x)[0][0] for x in \
                branchdc[:,F_BUSDC]]]
        Vdct = Vdc[[where(busdc[:,BUSDC_I]==x)[0][0] for x in \
                branchdc[:,T_BUSDC]]]
        Pfdc = pol*Vdcf*Ifdc ## power at the "from" bus
        Ptdc = pol*Vdct*(-Ifdc) ## power at the "to" bus


        ##----- slack/droop bus voltage and converter loss -----
        ## Initialisation
        Pc[slackdroopdc] = Pdc[slackdroopdc] - Ploss[slackdroopdc] ## Pc initialisation
        itslack = 0
        convergedslackdroop = 0

        ## dc slack bus loss calculation
        while not convergedslackdroop and itslack<=itmaxslackdroop:
           ## update iteration counter and convergence variable
           itslack += 1
           Pcprev = Pc.copy()

           ## update slack bus powers Ps, Qc and voltage Vc
           Ps[slackdroopdc], Qc[slackdroopdc], Vc[slackdroopdc] = calcslackdroop(
               Pc[slackdroopdc], Qs[slackdroopdc],  Vs[slackdroopdc], \
               Vf[slackdroopdc], Vc[slackdroopdc], Ztf[slackdroopdc], \
               Bf[slackdroopdc], Zc[slackdroopdc], \
               tolslackdroopint, itmaxslackdroopint)

           ## update slack bus losses
           Ploss[slackdroopdc]  = calclossac(Pc[slackdroopdc], Qc[slackdroopdc], \
                Vc[slackdroopdc], lossa[slackdroopdc], lossb[slackdroopdc], \
                losscr[slackdroopdc], lossci[slackdroopdc])

           ## update slack bus converter side power Pc
           Pc[slackdroopdc] = Pdc[slackdroopdc] - Ploss[slackdroopdc]

           ## slack bus tolerance check
           if max(abs(Pcprev[slackdroopdc] - Pc[slackdroopdc])) < tolslackdroop:
               convergedslackdroop = 1

        if not convergedslackdroop:
            stdout.write('\nSlackbus/Droop converter loss calculation of grid did NOT converge in %d iterations\n'% itslack)

        ## extended bus matrix update
        busVSC[cdci,PD] = bus[cdci,PD] - Ps[cdci]*baseMVA

        ## convergence check
        if abs(Ps_old - Ps).max() < tolacdc:
            converged = 1
    ## end of iteration
    timecalc = time() - t0

    ##-----  Post processing  -----
    ## convergence
    if converged:
        if output:
            stdout.write('\nSequential solution method converged in %d iterations\n'%it)
    else:
        stdout.write('\nSequential solution method did NOT converge after %d iterations\n'%it)

    ## converter limit check
    if limac == 1:
        for ii in arange(cdci.size):
            cvii = cdci[ii]
            limviol, _, plotarg = convlim(Ss[cvii], Vs[cvii], Vc[cvii], Ztf[cvii], \
                Bf[cvii], Zc[cvii], Icmax[cvii], Vcmax[cvii], Vcmin[cvii], i2edc[cvii+1], tollim, 1)
            if limviol != 0:     ## limits are hit
                if (convdc[cvii,CONVTYPE_DC] == DCSLACK):
                    stdout.write('\n  Slackbus converter %d is operating outside its limits.\n'%i2edc[cvii+1])
                elif (convdc[cvii,CONVTYPE_DC] == DCNOSLACK):
                    stdout.write('\n  Converter %d is operating outside its limits.\n'%i2edc[cvii+1])
            if convplotopt == 2 :
                convlimplot(plotarg, i2edc[cvii])

    ## update bus matrix
    bus[:,VM] = busVSC[:,VM]
    bus[:,VA] = busVSC[:,VA]

    ## dummy generators removal
    gen = genVSC[arange(gen.shape[0]),:]

    ## update busdc matrix
    busdc[:,PDC] = Pdc*baseMVA
    busdc[:,VDC] = Vdc

    ## update convdc matrix
    convdc[:,PCONV] = Ps*baseMVA
    convdc[:,QCONV] = Qs*baseMVA
    # new addition to convdc matrix
    convdc = c_[convdc,abs(Vc)]
    convdc = c_[convdc,angle(Vc)*180/pi]
    convdc = c_[convdc,Pc*baseMVA]
    convdc = c_[convdc,Qc*baseMVA]
    convdc = c_[convdc,Ploss*baseMVA]
    convdc = c_[convdc,abs(Vf)]
    convdc = c_[convdc,angle(Vf)*180/pi]
    convdc = c_[convdc,Psf*baseMVA]
    convdc = c_[convdc,Qsf*baseMVA]
    convdc = c_[convdc,Qcf*baseMVA]

    ## new addition to branchdc matrix
    branchdc = c_[branchdc,Pfdc*baseMVA]
    branchdc = c_[branchdc,Ptdc*baseMVA]

    #-----  internal to external bus renumbering  -----
    # remove dummy converters
    convdc = convdc[cdci,:]

    ## convert to external indexing
    ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"] = \
            baseMVA, bus, gen, branch

    pdc["baseMVAac"], pdc["baseMVAdc"], pdc["pol"], \
        pdc["busdc"], pdc["convdc"], pdc["branchdc"] = \
        baseMVAac, baseMVAdc, pol, busdc, convdc, branchdc

    ## Per unit internal to external data conversion
    pdc =int2extpu(ppc['baseMVA'],pdc);

    ## Undo the matrices sorting based on the bus numbers
    ppc['bus'] = ppc['bus'][i2ebus.argsort(),:]
    ppc['gen'] = ppc['gen'][i2egen.argsort(),:]
    ppc['branch'] = ppc['branch'][i2ebrch.argsort(),:]
    pdc['busdc'] = pdc['busdc'][i2ebusdc.argsort(),:]
    pdc['convdc'] = pdc['convdc'][i2econvdc.argsort(),:]
    pdc['branchdc'] = pdc['branchdc'][i2ebrchdc.argsort(),:]

    ## ac network internal to external bus numbering
    pdc, ppc = int2extac(i2eac, acdmbus, pdc, ppc)

    ## dc network internal to external bus numbering
    pdc = int2extdc(i2edcpmt, i2edc, pdc)

    ## generator outage inclusion
    gen1 = ppc['gen'] ## operational generators
    gen0[:,[PG, QG]] = 0 ## reset generator power injection
    ppc['gen'] = zeros((gon.shape[0]+goff.shape[0], gen1.shape[1]));
    ppc['gen'][gon,:] = gen1 ## include operational generators
    ppc['gen'][goff,:]  = gen0 ## include non-operational generators

    ## converter with outages inclusion
    conv1 = pdc['convdc']
    conv0 = c_[conv0, zeros((conv0.shape[0],conv1.shape[1] - conv0.shape[1]))]
    pdc['convdc'][conv0i, :] = conv0
    pdc['convdc'][conv1i, :] = conv1
    if conv0busi.shape[0]>0:
        pdc['busdc'][conv0busi[:,0], BUSAC_I] = conv0busi[:,1]

    ## dc branch outages inclusion
    brchdc1 = pdc['branchdc']
    brchdc0 = c_[brchdc0, zeros((brchdc0.shape[0], brchdc1.shape[1] - brchdc0.shape[1]))]
    pdc['branchdc'][brchdc0i,:] = brchdc0
    pdc['branchdc'][brchdc1i,:] = brchdc1

    ## ac branch outages inclusion
    if ppc['branch'].shape[0] == 0: ## all infinite buses
        # python start the index at 0
        brch0 = c_[brch0, zeros((brch0.shape[0], QT + 1 - brch0.shape[1]))] # not necessary anymore after rewriting the code
        ppc['branch'] = brch0;
    else:
        brch1 = ppc['branch']
        brch0 = c_[brch0, zeros((brch0.shape[0], brch1.shape[1] - brch0.shape[1]))];
        ppc['branch'][brch0i,:] = brch0
        ppc['branch'][brch1i,:] = brch1


    ##-----  output results  -----
    ## print results
    if output:
        printpf(ppc['baseMVA'], ppc['bus'], ppc['gen'], ppc['branch'],None,converged,timecalc)
        printdcpf(pdc['busdc'], pdc['convdc'], pdc['branchdc'])

    ##-----  output results  -----
    # as dict
    resultsac = {}
    resultsac['baseMVA'] = baseMVA
    resultsac['bus'] = bus
    resultsac['gen'] = gen
    resultsac['branch'] = branch

    resultsdc = {}
    resultsdc['baseMVAac'] = baseMVAac
    resultsdc['baseMVAdc'] = baseMVAdc
    resultsdc['pol'] = pol
    resultsdc['busdc'] = busdc
    resultsdc['convdc'] = convdc
    resultsdc['branchdc'] = branchdc

    # if nargout == 2 || nargout == 3 || nargout == 4
        # baseMVA = resultsac;
        # bus = resultsdc;
        # gen = converged;
        # branch =  timecalc;
    # end
    input()
    return resultsac, resultsdc, converged, timecalc

if __name__ == '__main__':
    resultsac, resultsdc, converged, timecalc =runacdcpf()
