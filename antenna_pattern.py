import numpy as np
import sys

Omega_ES = 2.0 * np.pi / 365.25 / 24. / 60. / 60.
Omega_E = 2.0 * np.pi / 24. / 60. / 60.


# The following values are taken from Jaranowski, Krolak and Schutz PRDxx. xxxxxx
DETECTOR_CONFIGURATIONS = {

    "LIGO Hanford":{
        "lambda": 46.45 * np.pi / 180.,
        "longitude": 119.41 * np.pi / 180.,
        "gamma": 171.8 * np.pi / 180.
    },

    "LIGO Livingston":{
        "lambda": 30.56 * np.pi / 180.,
        "longitude": 90.77 * np.pi / 180.,
        "gamma": 243.0 * np.pi / 180.
    },

    "Virgo":{
        "lambda": 43.63 * np.pi / 180.,
        "longitude": -10.5 * np.pi / 180.,
        "gamma": 116.5 * np.pi / 180.
    }

}


def _antenna_a(t, alpha, delta, detparam):
    lmd = detparam["lambda"]
    gamma = detparam["gamma"]

    coslmd = np.cos(lmd)
    sinlmd = np.sin(lmd)
    cos2lmd = np.cos(2.*lmd)
    sin2lmd = np.sin(2.*lmd)
    cosalp = np.cos(alpha - Omega_E * t)
    sinalp = np.sin(alpha - Omega_E * t)
    cos2alp = np.cos(2.*(alpha - Omega_E * t))
    sin2alp = np.sin(2.*(alpha - Omega_E * t))
    cos2del = np.cos(2.*delta)
    sin2del = np.sin(2.*delta)
    cos2gam = np.cos(2.*gamma)
    sin2gam = np.sin(2.*gamma)

    return sin2gam * (3. - cos2lmd) * (3. - cos2del) * cos2alp / 16.\
        - cos2gam * sinlmd * (3. - cos2del) * sin2alp / 4.\
            + sin2gam * sin2lmd * sin2del * cosalp / 4.\
                - cos2gam * coslmd * sin2del * sinalp / 2.\
                    + 3. * sin2gam * (1. + cos2lmd) * (1. + cos2del) / 16.


def _antenna_b(t, alpha, delta, detparam):
    lmd = detparam["lambda"]
    gamma = detparam["gamma"]

    coslmd = np.cos(lmd)
    sinlmd = np.sin(lmd)
    cos2lmd = np.cos(2.*lmd)
    sin2lmd = np.sin(2.*lmd)
    cosalp = np.cos(alpha - Omega_E * t)
    sinalp = np.sin(alpha - Omega_E * t)
    cos2alp = np.cos(2.*(alpha - Omega_E * t))
    sin2alp = np.sin(2.*(alpha - Omega_E * t))
    cosdel = np.cos(delta)
    sindel = np.sin(delta)
    cos2gam = np.cos(2.*gamma)
    sin2gam = np.sin(2.*gamma)

    return cos2gam * sinlmd * sindel * cos2alp\
        + sin2gam * (3. - cos2lmd) * sindel * sin2alp / 4.\
            + cos2gam * coslmd * cosdel * cosalp\
                + sin2gam * sin2lmd * cosdel * sinalp / 2.

def antenna_Fplus(t, alpha, delta, psi, detname):

    if not (detname in DETECTOR_CONFIGURATIONS):
        sys.ext("invalid detector name.")

    detparam = DETECTOR_CONFIGURATIONS[detname]

    return _antenna_a(t, alpha, delta, detparam) * np.cos(2. * psi)\
        + _antenna_b(t, alpha, delta, detparam) * np.sin(2. * psi)


def antenna_Fcross(t, alpha, delta, psi, detname):

    if not (detname in DETECTOR_CONFIGURATIONS):
        sys.ext("invalid detector name.")

    detparam = DETECTOR_CONFIGURATIONS[detname]

    return _antenna_b(t, alpha, delta, detparam) * np.cos(2. * psi)\
        - _antenna_a(t, alpha, delta, detparam) * np.sin(2. * psi)



if __name__ == '__main__':


    print("antenna_pattern.py")
    print("----------------------------")
    print("\
    Implementation check.\n\
    The average of F+^2 over source locations and polarizations should be 1/5.\n\
    The same for Fx^2.\n\
    The average of F+*Fx vanishes.\n\
    Confirm above properties by MonteCarlo integration.\
    ")
    print("----------------------------")


    N = 100000
    alpha = np.random.uniform(0.0, 2.0* np.pi, (N,))
    delta = np.random.uniform(-np.pi/2., np.pi/2., (N,))
    psi = np.random.uniform(0., 2.*np.pi, (N,))
    Fp2list = []
    Fc2list = []
    FpFclist = []
    t = 0.0
    for n in range(N):
        a = alpha[n]
        d = delta[n]
        p = psi[n]
        Fp = antenna_Fplus(t, a, d, p, "LIGO Hanford")
        Fc = antenna_Fcross(t, a, d, p, "LIGO Hanford")
        Fp2list.append(Fp * Fp * np.cos(d))
        Fc2list.append(Fc * Fc * np.cos(d))
        FpFclist.append(Fp * Fc * np.cos(d))

    print(f"<Fp^2> - 1/5 = {np.sum(Fp2list) * np.pi / 2. / N - 0.2 :e}")
    print(f"<Fc^2> - 1/5 = {np.sum(Fc2list) * np.pi / 2. / N - 0.2 :e}")
    print(f"<FpFc> = {np.sum(FpFclist) * np.pi / 2. / N :e}")