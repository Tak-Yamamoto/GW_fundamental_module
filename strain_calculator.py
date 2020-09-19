import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt
import sys

G = const.G.cgs.value
c = const.c.cgs.value
Msun = const.M_sun.cgs.value
Mpc = 3.08568e+24

def strain(D, M, f):
    """
    D: [pc]
    M: chirp mass [Msun]
    f: frequency [Hz]
    """

    dist = D * Mpc
    Mc = M * Msun
    h = (2.0 / np.pi / (c**4.0) / dist) * ((np.pi*G*Mc)**(5.0/3.0)) * f**(2.0/3.0)

    return h


def detectable_distance(M, f, SNRth, Tobs):
    """
    M: chirp mass [Msun]
    f: frequency [Hz]
    SNRth: the threshold SNR
    Tobs: observation time [sec]
    """
    Mc = M * Msun
    Sn = 1.0e-46
    distmin = (2.0 / np.pi / (c**4.0)) * ((np.pi*G*Mc)**(5.0/3.0)) * f**(2.0/3.0) / SNRth * (Tobs/Sn)**(0.5)

    return distmin / Mpc




if __name__ == '__main__':
    
    """
    D = float(sys.argv[1])
    M = float(sys.argv[2])
    f = float(sys.argv[3])
    """

    D = 1e-3
    f = 300.0

    Mlist = np.linspace(1.0e-7, 1e-4, 1000)
    strainlist = []
    for M in Mlist:
        strainlist.append(strain(D, M, f))

    plt.figure()
    plt.loglog(Mlist, strainlist)
    plt.show()


