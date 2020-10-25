from strain_calculator import strain
import numpy as np
import matplotlib.pyplot as plt
import sys


f = float(sys.argv[1])
print("Frequency: {}[Hz]".format(f))


Mlist = 10.0 ** (np.arange(-9, 0))
Dlist = 10.0 ** (np.arange(-5, 1, 2))

plt.figure()

for D in Dlist:

    strainlist = []
    for M in Mlist:
        strainlist.append(strain(D, M, f))

    plt.loglog(Mlist, strainlist, label='{:.1f}kpc'.format(D*1e+3))

plt.xlabel('chirp mass[Msun]')
plt.ylabel('strain')
plt.legend()
plt.show()



