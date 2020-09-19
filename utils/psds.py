"""
psds.py
"""


def PSD_LIGO_O1():
    """ Taken from LIGO's tutorial """
    return lambda freqs: (1.e-22*(18./(0.1+freqs))**2)**2+0.7e-23**2+((freqs/2000.)*4.e-23)**2




if __name__=='__main__':

    import numpy as np
    import matplotlib.pyplot as plt

    fmax = 2048
    df = 1.0/8.0
    N = int(fmax / df)
    f = (np.arange(N)+1)*df

    psd = PSD_LIGO_O1()

    plt.figure()
    plt.loglog(f, np.sqrt(psd(f)))
    plt.show()

