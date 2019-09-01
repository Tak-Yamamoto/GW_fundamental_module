import matplotlib.pyplot as plt
from scipy.special import jv
import numpy as np


zlist = [100.0, 1000.0, 2000.0]
vrange_list = [np.arange(10.0, 150.0), np.arange(900, 1100), np.arange(1900, 2100)]
for n in range(3):
    plt.figure()
    z = zlist[n]
    vrange = vrange_list[n]
    jvlist = []
    for v in vrange:
        jvlist.append(jv(v, z))
    plt.plot(vrange, jvlist, '-o', markersize=3, label='z={:.1f}'.format(z))
    plt.vlines(z+2*np.sqrt(z), ymin=-0.15, ymax = 0.15)
    plt.grid()
plt.show()

