from strain_calculator import detectable_distance
import numpy as np
import matplotlib.pyplot as plt
import sys


f = float(sys.argv[1])
SNRth = float(sys.argv[2])
print("Frequency: {}[Hz]".format(f))
print("SNR threshold: {}".format(SNRth))


Mlist = 10.0 ** (np.arange(-9, 0))
Tlist = np.array([1.0, 5.0, 10.0]) * (3.0e+7)

plt.figure()


for Tobs in Tlist:
    distlist = []
    for M in Mlist:
        distlist.append(detectable_distance(M, f, SNRth, Tobs))
    plt.loglog(Mlist, distlist, '-o', label='{:.0f}yr'.format(Tobs/(3e+7)))

plt.xlabel('chirp mass[Msun]')
plt.ylabel('detectable distance [Mpc]')
plt.legend()
plt.savefig('figure/detectable_distance.png')



