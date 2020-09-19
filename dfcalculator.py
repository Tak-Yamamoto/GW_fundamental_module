import numpy as np
import matplotlib.pyplot as plt

c = 3.0e+10
Munit = 4.925490e-6
factor = 96 * (np.pi**(8.0/3.0)) / 5.0
f0 = 1.0e-2

for Mc in [1.0e+2]:
    mass = (Munit * Mc) ** (5.0/3.0)
    df = factor * mass * (f0**(11.0/3.0))
    print(Mc, df)

    

Mc = 1.0e-6
f0 = 100.0
massfactor = (Munit * Mc) ** (5.0/3.0)
df = factor * massfactor * (f0**(11.0/3.0))


