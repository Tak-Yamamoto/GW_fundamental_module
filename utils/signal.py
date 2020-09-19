import numpy as np
import sys

TSYGWModule_WindowList = ['tukey']

def get_windowedstrain(strain, windowtype, **params):

    """
    Parameters
    ---------------------------------------------
    strain: numpy.ndarray
        A time-domain (real-value) signal to be windowed.
    windowtype: string
        Specify which window we use.
    params: dict
        Parameters charachterizing a window function.
    """

    


def get_whitenedstrain(strain, dt, psd):
    """
    Parameters
    ----------------------------------------------
    strain: numpy.ndarray
        A time-domain (real-value) signal to be whitened
    dt: float
        The time step
    psd: function
        A function takes frequency array and returns psd array.

    Returns
    ------------------------------------------------
    strain_wh: numpy.ndarray
        The whitened time-domain signal.
    """

    N = len(strain)
    freq = np.fft.rfftfreq(N, dt)
    df = abs(freq[1] - freq[0])
    psdarray = psd(freq)
    hf = np.fft.rfft(strain) * dt / psdarray

    return np.fft.irfft(hf) * df