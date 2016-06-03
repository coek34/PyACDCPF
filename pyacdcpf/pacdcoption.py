"""Used to retrieve a PYACDCPF options vector.
"""

ACDCPF_OPTIONS = [
    ('TOLACDC', 1e-8, 'tolerance ac/dc power flow'),

    ('ITMAXACDC', 10, 'maximum iterations ac/dc power flow'),

    ('TOLDC', 1e-8, 'tolerance dc power flow (Newton\'s method)'),

    ('ITMAXDC', 10, 'maximum iterations dc power flow (Newton\'s method)'),

    ('TOLSLACKDROOP', 1e-8, 'tolerance dc slack bus iteration'),

    ('ITMAXSLACKDROOP', 10, 'maximum iterations dc slack bus iteration'),

    ('TOLSLACKDROOPINT', 1e-8, 'tolerance internal slack bus iteration (Newton\'s method)'),

    ('ITMAXSLACKDROOPINT', 10, 'maximum iterations internal slack bus iteration (Newton\'s method)'),

    ('MULTSLACK', 0, '''multiple dc slack buses (dc voltage controlling converters) per dc grid
0 - only 1 dc voltage controlling per dc grid allowed
1 - more than 1 dc voltage controlling per dc grid allowed'''),

    ('LIMAC', 0, '''enforce ac converter limits
0 - do NOT enforce limits
1 - enforce converter current and voltage limits'''),

    ('LIMDC', 0, 'enforce dc converter limits (not implemented)'),

    ('TOLLIM', 1e-2, 'maximum difference between subsequent violations'),

    ('OUTPUT', 1, 'print output'),

    ('CONVPLOTOPT', 0, '''plot converter limit violations
0 - do not plot converter limit violations
1 - plot only converter limit violations
2 - plot converter limit violations and end situation''')

]

def pacdcoption(pacdcopt=None, **kw_args):
    """
    Used to set and retrieve a PYACDCPF options vector.

    C{opt = pacdcoption()} returns the default options vector

    C{opt = pacdcoption(NAME1=VALUE1, NAME2=VALUE2, ...)} returns the default
    options vector with new values for the specified options, NAME# is the
    name of an option, and VALUE# is the new value.

    C{opt = pacdcoption(OPT, NAME1=VALUE1, NAME2=VALUE2, ...)} same as above
    except it uses the options vector OPT as a base instead of the default
    options vector.

    Examples::
        opt = pacdcoption(ITMAXACDC=2, TOLDC=1e-4);
        opt = pacdcoption(opt, MULTSLACK=0, OUTPUT=1)

    @author:Jef Beerten (KU Leuven)
    @author:Roni Irnawan (Aalborg University)    
    """

    default_pacdcopt = {}

    options = ACDCPF_OPTIONS

    for name, default, _ in options:
        default_pacdcopt[name.upper()] = default

    pacdcopt = default_pacdcopt if pacdcopt == None else pacdcopt.copy()

    pacdcopt.update(kw_args)

    return pacdcopt