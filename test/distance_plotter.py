from strain_calculator import detectable_distance
import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy.constants as const
import astropy.units as u

# physical constants
G = const.G.cgs.value
c = const.c.cgs.value
M_sun = const.M_sun.cgs.value
sday = u.sday.cgs
syr = u.year.cgs
print(syr)


# calculate merger time
def time_to_merger(mass, freq):

    Mc = mass * M_sun
    omega = 2 * np.pi * freq
    tc = 5.0 * (c**5.0) / (G**(5.0/3.0)) / (Mc**(5.0/3.0)) / (omega**(8.0/3.0)) / 256.0
    return tc


f = float(sys.argv[1])
SNRth = float(sys.argv[2])
print("Frequency: {}[Hz]".format(f))
print("SNR threshold: {}".format(SNRth))


Mlist = 10.0 ** (np.arange(-9, 0))
Tlist = [3e+7]

plt.figure()

print("GW150914", time_to_merger(30.0, 10.0))
print("GW170817", time_to_merger(1.4, 10.0))

for Tobs in Tlist:
    distlist = []
    for M in Mlist:
        distlist.append(detectable_distance(M, f, SNRth, Tobs))
        tc = time_to_merger(M, 10.0)
        print("{:.3e} {:.3e}".format(M, tc))
    plt.loglog(Mlist, distlist, '-o', label='{:.0f}yr'.format(Tobs/(3e+7)))

plt.xlabel('chirp mass[Msun]')
plt.ylabel('detectable distance [Mpc]')
plt.legend()
plt.savefig('figure/detectable_distance.png')


