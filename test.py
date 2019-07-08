import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy import constants as const

G = const.G.cgs.value
c = const.c.cgs.value
Msun = const.M_sun.cgs.value
Mpc = 3.08568e+24

def get_integral(freq, psd, fmax):

    psd = psd[freq<fmax]
    freq = freq[freq<fmax]
    df = freq[1] - freq[0]

    # integration
    integral_of_chirp = np.trapz(freq**(-7.0/3.0) / (psd**2.0), dx=df)
    return integral_of_chirp

def get_FISCO(mass):

    '''
    fmax can be estimated as 2fisco.
    '''
    M = mass * Msun
    fisco = c**3.0 / (np.pi * (6**1.5)) / (G*M)
    return fisco
    

def get_SNR(mass, distance, freq, psd, fmax):
    '''
    M: chirp mass [Msun]
    D: distance [Mpc]
    '''

    integral_of_chirp = get_integral(freq, psd, fmax)
    D = distance * Mpc
    M = mass * Msun
    rhosq = 2.0 * ((G*M)**(5.0/3.0)) / c**3.0 / (D**2.0) / (15.0 * np.pi**(4.0/3.0)) * integral_of_chirp
    return np.sqrt(rhosq)


def get_SightDistance(mass, snrth, freq, pdf, fmax):

    integral_of_chirp = get_integral(freq, pdf, fmax)
    M = mass * Msun
    numfactor = (2.0/15.0)**(0.5) / np.pi**(2.0/3.0)
    dist = numfactor * (G*M)**(5.0/6.0) / c**1.5 / snrth * (integral_of_chirp**0.5)
    return dist / Mpc



# Load the psd
psddir = 'data/PSDs/'
psdfile = 'ZERO_DET_high_P.txt'
data = np.genfromtxt(psddir + psdfile)
freq = data[:,0]
psd = data[:,1]

# interpolation
psd_interp = interp1d(freq, psd)

# resampling
N = 10000
df = (freq.max() - freq.min()) / N
freq = np.arange(N) * df + freq.min()
psd = psd_interp(freq)

mass = 0.01 
distance = 400.0
snrth = 8.0
fmax = 2.0 * get_FISCO(mass)

print("# chirp mass: {}[Msun]".format(mass))
#print("# ISCO frequency: {}[Hz]".format(fmax/2.0))
print("# SNR threshold: {}".format(snrth))

dist = get_SightDistance(mass, snrth, freq, psd, fmax)
print("# sight distance: {}[Mpc]".format(dist))




