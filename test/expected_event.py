import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import G, c, M_sun

M_sun = M_sun.cgs.value
c = c.cgs.value
G = G.cgs.value
pc = 3.08568e+18  # [cm]
syr = 3.15576e+07

f_pbh = 0.1
M_DM = (1.7e+15) * M_sun
R = 3.0 * 1e+6 * pc  # galaxy group

def get_df(m_pbh, fgw):

    num_fac = 96.0 * (np.pi**(8.0/3.0)) / 5.0 
    mass_fac = (G * m_pbh * M_sun / (c**3.0))**(5.0/3.0)
    return num_fac * mass_fac * (fgw**(11.0/3.0))

def get_distance(m_pbh, fgw, h0):

    num_fac = 4.0 * (np.pi**(2.0/3.0))
    mass_fac = (G * m_pbh * M_sun / (c**3.0))**(5.0/3.0)
    return num_fac * mass_fac * c * (fgw**(2.0/3.0)) / h0

def get_event_num(m_pbh, fgw, h0):

    num_density = f_pbh * M_DM / (m_pbh * M_sun) / (R**3.0)
    print(m_pbh, num_density)
    dist = get_distance(m_pbh, fgw, h0)
    return num_density * (dist**3.0)

fgw = 100.0
h0 = 1.0e-25

print("df/dt: {}[Hz/s]".format(get_df(1.0e-6, fgw)))
print("get_distance: {}[pc]".format(get_distance(1.0e-6, fgw, h0) / pc))

DT1 = 30.0
DT2 = 3000.0
DF1 = 1.0 / (syr * DT1)
DF2 = 1.0 / (syr * DT2)


Nplot = 200
m_index_max = -5.0
m_index_min = -8.0
m_max = 10.0**m_index_max
m_min = 10.0**m_index_min

list_pbhmass = 10.0**(np.linspace(m_index_min, m_index_max, Nplot))
list_df = []
list_eventnum = []


for m_pbh in list_pbhmass:

    list_df.append(get_df(m_pbh, fgw))
    list_eventnum.append(get_event_num(m_pbh, fgw, h0))


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

ax1.loglog(list_pbhmass, list_eventnum, 'b', lw=2)
ax2.loglog(list_pbhmass, list_df, 'r', lw=2)

ax1.hlines(1.0, xmin=m_min, xmax=m_max, color='b', linestyle='-.')
ax2.hlines(DF2, xmin=m_min, xmax=m_max, color='r', linestyle=':')
ax2.hlines(DF1, xmin=m_min, xmax=m_max, color='r', linestyle='--')

ax1.text(1.5e-8, 1.2, ">1event expected", color="b")
ax2.text(2.6e-6, 4.0e-10, r"$\Delta T = 30$sec", color="r")
ax2.text(2.6e-6, 4.0e-12, r"$\Delta T = 3000$sec", color="r")

ax1.set_ylabel('expected number of event')
ax1.set_ylim([1.0e-6, 1.0e2])

ax2.set_xlabel(r'$m_\mathrm{PBH} [M_\odot]$')
ax2.set_ylabel(r'$df/dt$ [Hz/s]')
ax2.set_ylim([1.0e-14, 1.0e-7])


plt.xlim([m_min, m_max])
ax1.grid(which='both', axis='x')
plt.tight_layout()
plt.savefig('expected_number_of_event.pdf')
plt.show()




