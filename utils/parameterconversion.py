"""
parameterconversion.py
"""

import numpy as np
import sys


def ChirpMass_MassRatio_To_ComponentMasses(Mc, eta):
    """
    Parameters
    -----------------------------------
    Mc: float or numpy.ndarray
        Chirp mass(es).
    eta: float or numpy.ndarray
        Mass ratio(s) defined as eta=m1/m2 (m1>m2).

    Returns
    ------------------------------------
    m1: float or numpy.ndarray
        Larger mass
    m2: float or numpy.ndarray
        Smaller mass
    """

    if isinstance(Mc, np.ndarray) and isinstance(eta, np.ndarray):
        if Mc.shape != eta.shape:
            sys.exit("Error: Mc and eta have different sizes.")

    m1 = Mc * ((1+eta)/(eta**3.0))**(1.0/5.0)
    m2 = eta * m1
    return m1, m2


def ComponentMasses_To_ChirpMass_MassRatio(m1, m2):
    """
    Parameters
    ------------------------------------
    m1: float or numpy.ndarray
        Larger mass
    m2: float or numpy.ndarray
        Smaller mass  

    Returns
    -----------------------------------
    Mc: float or numpy.ndarray
        Chirp mass(es).
    eta: float or numpy.ndarray
        Mass ratio(s) defined as eta=m1/m2 (m1>m2).
    """


    if isinstance(m1, np.ndarray) and isinstance(m2, np.ndarray):
        if m1.shape != m2.shape:
            sys.exit("Error: m1 and m2 have different sizes.")

    Mc = ( (m1*m2)**3.0 / (m1+m2))**(1.0/5.0)
    eta = m2 / m1
    return Mc, eta