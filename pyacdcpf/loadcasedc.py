"""Loads a PYACDCPF case dictionary.
"""

import sys

from os.path import basename, splitext, exists

from copy import deepcopy

from numpy import array, zeros, ones, c_

from scipy.io import loadmat

from pyacdcpf._compat import PY2

if not PY2:
    basestring = str


def loadcasedc(casefiledc, return_as_obj=True):
    """
    Returns the individual data matrices or an dict containing them
    as values.

    Here C{casefiledc} is either a dict containing the keys C{baseMVA}, C{baseMVAdc},
    C{pol}, C{busdc}, C{convdc}, C{branchdc}, or a string containing the name
    of the file. If C{casefiledc} contains the extension '.mat' or '.py', then
    the explicit file is searched. If C{casefiledc} containts no extension, then
    L{loadcase} looks for a '.mat' file first, then for a '.py' file.  If the
    file does not exist or doesn't define all matrices, the function returns
    an exit code as follows:

        0.  all variables successfully defined
        1.  input argument is not a string or dict
        2.  specified extension-less file name does not exist
        3.  specified .mat file does not exist
        4.  specified .py file does not exist
        5.  specified file fails to define all matrices or contains syntax
            error

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """
    
    info = 0

    # read data into case object
    if isinstance(casefiledc, basestring):
        # check for explicit extension
        if casefiledc.endswith(('.py', '.mat')):
            rootname, extension = splitext(casefiledc)
            fname = basename(rootname)
        else:
            # set extension if not specified explicitly
            rootname = casefiledc
            if exists(casefiledc + '.mat'):
                extension = '.mat'
            elif exists(casefiledc + '.py'):
                extension = '.py'
            else:
                info = 2
            fname = basename(rootname)

        lasterr = ''

        ## attempt to read file
        if info == 0:
            if extension == '.mat':       ## from MAT file
                try:
                    d = loadmat(rootname + extension, struct_as_record=True)
                    if 'pdc' in d or 'mdc' in d:    ## it's a MAT/PYACDCPF dict
                        if 'pdc' in d:
                            struct = d['pdc']
                        else:
                            struct = d['mdc']
                        val = struct[0, 0]

                        s = {}
                        for a in val.dtype.names:
                            s[a] = val[a]
                    else:                 ## individual data matrices
                        d['version'] = '1'

                        s = {}
                        for k, v in d.items():
                            s[k] = v

                    s['baseMVAac'] = s['baseMVAac'][0]  # convert array to float
                    s['baseMVAdc'] = s['baseMVAdc'][0]  # convert array to float

                except IOError as e:
                    info = 3
                    lasterr = str(e)
            elif extension == '.py':      ## from Python file
                try:
                    if PY2:
                        execfile(rootname + extension)
                    else:
                        exec(compile(open(rootname + extension).read(),
                                     rootname + extension, 'exec'))

                    try:                      ## assume it returns an object
                        s = eval(fname)()
                    except ValueError as e:
                        info = 4
                        lasterr = str(e)
                    ## if not try individual data matrices
                    if info == 0 and not isinstance(s, dict):
                        s = {}
                        s['version'] = '1'
                        if return_as_obj:
                            try:
                                s['baseMVAac'], s['baseMVAdc'], s['pol'], \
                                    s['busdc'], s['convdc'], \
                                    s['branchdc'] = eval(fname)()
                            except ValueError as e:
                                info = 4
                                lasterr = str(e)
                        else:
                            try:
                                s['baseMVAac'], s['baseMVAdc'], s['pol'], \
                                    s['busdc'], s['convdc'], \
                                    s['branchdc'] = eval(fname)()
                            except ValueError as e:
                                info = 4
                                lasterr = str(e)

                except IOError as e:
                    info = 4
                    lasterr = str(e)


                if info == 4 and exists(rootname + '.py'):
                    info = 5
                    err5 = lasterr

    elif isinstance(casefiledc, dict):
        s = deepcopy(casefiledc)
    else:
        info = 1

    # check contents of dict
    if info == 0:
        # check for required keys
        if (s['baseMVAac'] is None or s['baseMVAdc'] is None \
            or s['busdc'] is None or s['convdc'] is None \
            or s['branchdc'] is None):
            info = 5  ## missing some expected fields
            err5 = 'missing data'
        else:
            ## all fields present, copy to pdc
            pdc = deepcopy(s)
            if not hasattr(pdc, 'version'):  ## hmm, struct with no 'version' field
                if pdc['branchdc'].shape[1] < 10:    ## version 2 has 21 or 25 cols
                    pdc['version'] = '1'
                else:
                    pdc['version'] = '2'

            #if (pdc['version'] == '1'):
                # convert from version 1 to version 2
                # currently do nothing
                
    if info == 0:  # no errors
        ## add voltage droop parameters if not defined in input files
        if pdc['convdc'].shape[1] < 24:
            shift = 24-pdc['convdc'].shape[1]
            tmp = zeros((pdc['convdc'].shape[0], shift))
            pdc['convdc'] = c_[ pdc['convdc'], tmp ]
        
        if return_as_obj:
            return pdc
        else:
            result = [pdc['baseMVAac'], pdc['baseMVAdc'], pdc['pol'], pdc['busdc'], \
            pdc['convdc'], pdc['branchdc']]
            return result
    else:  # error encountered
        if info == 1:
            sys.stderr.write('Input arg should be a case or a string '
                             'containing a filename\n')
        elif info == 2:
            sys.stderr.write('Specified case not a valid file\n')
        elif info == 3:
            sys.stderr.write('Specified MAT file does not exist\n')
        elif info == 4:
            sys.stderr.write('Specified Python file does not exist\n')
        elif info == 5:
            sys.stderr.write('Syntax error or undefined data '
                             'matrix(ices) in the file\n')
        else:
            sys.stderr.write('Unknown error encountered loading case.\n')

        sys.stderr.write(lasterr + '\n')

        return info
