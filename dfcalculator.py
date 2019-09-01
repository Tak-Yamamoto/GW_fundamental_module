import numpy as np
import matplotlib.pyplot as plt



Munit = 4.925490e-6
factor = 96 * (np.pi**(8.0/3.0)) / 5.0

for Mc in (np.arange(1.0, 11.0) * 1e-7):
    mass = (Munit * Mc) ** (5.0/3.0)
    df = factor * mass * (100**(11.0/3.0))
    print(Mc, df)


