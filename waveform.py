"""
waveform.py
"""

import numpy as np
import astropy.constants as const
from scipy.signal import tukey

G = const.G.cgs.value
c = const.c.cgs.value
M_sun = const.M_sun.cgs.value
tunit = G * M_sun / c / c / c

def tsy_get_fd_25PNwaveform(m1, m2, a1, a2, fs, duration, fmin=20.0):
    """
    Generate a 2.5PN waveform in fourier domain.
    Implementation is based on T.Tanaka&H.Tagoshi, PRD 62, 082001
    """

    mtot = m1 + m2
    eta = m1*m2/mtot/mtot
    xs = 0.5 * (a1 + a2)
    xa = 0.5 * (a1 - a2)
    f = np.fft.rfftfreq(int(duration*fs), 1.0/fs)
    fmax = 80.0
    kmin = np.argmin(abs(f-fmin))
    kmax = np.argmin(abs(f-fmax))
    w = tukey(kmax - kmin, 1.0/8.0)
    x = np.pi * mtot * f * tunit

    psi1 = 3.0 * (x**(-5./3.)) / 128. / eta
    psi2 = (3715./84. + 55.*eta) / 384. / eta / x
    psi3 = ((113.-86.*eta) * xs + 113. * xa * (m2-m1)/mtot - 48.*np.pi) / 128. / eta * (x**(-2./3.))
    psi4 = 3. * (15293365./508032. + 27415.*eta/504. + 3085.*eta*eta/72. + (30.+275./4.*eta)*(xs*xs-xa*xa)) * (x**(-1./3.))
    psi5 = np.pi / 128. / eta * (38645./252. + 5.*eta) * np.log(f)
    psi5=0.0

    psi = psi1 + psi2 + psi3 + psi4 + psi5 - f*5.0

    hp = f**(-7./6.) * np.exp(1.j * psi)
    hp[:kmin] = 0.0
    hp[kmax:] = 0.0
    hp[kmin:kmax] = hp[kmin:kmax] * w
    hc = f**(-7./6.) * np.exp(1.j * psi) * 1.j
    hc[:kmin] = 0.0
    hc[kmax:] = 0.0
    hc[kmin:kmax] = hc[kmin:kmax] * w

    return f, hp, hc





if __name__=='__main__':

    import matplotlib.pyplot as plt
    from pycbc.waveform import get_fd_waveform

    fmin = 20.0
    fs = 4096
    dt = 1.0/fs
    m1 = 30.0
    m2 = 30.0
    f, hp2, hc2 = tsy_get_fd_25PNwaveform(m1, m2, 0.0, 0.0, fs, 8.0, fmin)
    hp, hc = get_fd_waveform(approximant="TaylorF2", mass1=m1, mass2=m2, delta_f=1.0/8.0, f_lower=fmin)

    plt.figure()
    plt.loglog(hp.sample_frequencies, hp)
    plt.loglog(f, np.real(hp2))
    plt.show()


    ht = np.fft.fft(hp2) / dt

    plt.figure()
    plt.plot(ht.real)
    plt.plot(ht.imag)
    #plt.loglog(f, hp2.real)
    #plt.loglog(freq, hf.real)
    plt.show()
