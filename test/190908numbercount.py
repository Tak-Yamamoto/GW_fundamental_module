import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import G, c, M_sun

Omega_m = 0.33
Omega_k = 0.0
Omega_lmd = 0.67
H0 = 65.0  # km / s / Mpc
G = G.cgs.value * 1e-15  #[km/g/s2]
c = c.cgs.value * 1e-5  #[km/s]
M_sun = M_sun.cgs.value  #[g]
Mpc = 3.08568e+24  #[cm]


def get_chirp_mass(m1, m2):

    return (m1*m2 * ((m1+m2)**(-1.0/3.0)))**(3.0/5.0)

def get_Eofz(z):

    return (Omega_m * ((1.0+z)**3.0) + Omega_k * ((1+z)**2.0) + Omega_lmd)**0.5


def get_luminosity_distance(z):

    # luminosity distance in Mpc

    N = 500
    z_arr = np.linspace(0.0, z, N)
    E_arr = get_Eofz(z_arr)

    lum_dist = c / H0 * (1+z) * np.trapz(1.0/E_arr, z_arr, z_arr[1]-z_arr[0])

    return lum_dist

def get_dVdz(z):

    dist_m = get_luminosity_distance(z) / (1+z)
    E = get_Eofz(z)

    return 4.0 * np.pi * c / H0 * (dist_m**2.0) / E


def get_N_per_comvol(fr, Mc, merger_rate_per_comvol):

    num_fac = 5.0*np.pi/96.0
    mass_fac = c**5.0 / ((G*Mc*M_sun)**(5.0/3.0))
    return num_fac * mass_fac * merger_rate_per_comvol / ((np.pi*fr)**(11.0/3.0))


def get_N_fr_z(fr, Mc, merger_rate_per_comvol, z):

    Nper_comvol = get_N_per_comvol(fr, Mc, merger_rate_per_comvol)
    dVdz = get_dVdz(z)

    return Nper_comvol * dVdz


DT = 3000.0
DF = 1.0 / DT
DTheta = 1.0e-2
zmax = 4.0
f = 1.0e-3

z_arr = np.linspace(0.0, zmax, 1000)
MergerRate = 1000.0

N_f_z = get_N_fr_z(f*(1.0+z_arr), 1.4, MergerRate, z_arr)
integrand_arr = (Dtheta**2.0 / 4.0 / np.pi) * N_f_z * (1.0 + z_arr) * DF
N_obtained = np.trapz(integrand_arr, z_arr, z_arr[1] - z_arr[0])
print(N_obtained)



print(get_N_per_comvol(1e-3, 1.4, MergerRate))


z = 4.0
print("z {}, luminosity distance {}[Mpc]".format(z, get_luminosity_distance(z)))

    


