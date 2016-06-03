============
Introduction
============

PYACDCPF is a a sequential ac/dc power flow solver. It is a port of
MatACDC_ to the Python_ programming language. Current
features include:

* AC and DC sequential power flow and
* HVDC converter limit plotting using matplotlib_


Installation
============

PYACDCPF depends upon:

* Python_ 3.3 or later 
* PYPOWER_ 5.0 or later
* matplotlib_ 1.0 or later 
* NumPy_ 1.0 or later and
* SciPy_ 0.9 or later.

It can be installed after unpacking the package using::

  $ python setup.py install


Using PYACDCPF
=============

PYACDCPF module can be used first by importing runacdcpf from PYACDCPF module::

  >>> from pyacdcpf import runacdcpf
  >>> resultac, resultdc, converged, te = pfacdc.runacdcpf()


Support
=======

Questions and comments regarding PYACDCPF should be directed to:

    roni.irnawan@google.com


License & Copyright
===================

Copyright (C) 2012 Jef Beerten (KU Leuven)  
Copyright (C) 2016 Roni Irnawan (Aalborg University)


Links
=====

* MatACDC_ by Jef Beerten
* PYPOWER_ from PSERC (Cornell)
* MATPOWER_ from PSERC (Cornell)


.. _Python: http://www.python.org
.. _SciPy: http://www.scipy.org
.. _NumPy: http://www.numpy.org
.. _matplotlib: http://www.matplotlib.org
.. _MATPOWER: http://www.pserc.cornell.edu/matpower/
.. _PYPOWER: http://pypi.python.org/pypi/PYPOWER
.. _MatACDC: http://www.esat.kuleuven.be/electa/teaching/matacdc
