import numpy as np
from scipy.special import wofz

def voigtfwhm(x, p):
    """
    Voigt function with normalized parameters.
    
    Parameters:
    x : array-like
        Input coordinates where the function will be calculated.
    p : list or array
        Parameters [A, x0, wG, wL] where:
        A : Amplitude (area under the curve)
        x0 : Peak center position
        wG : FWHM of the Gaussian component
        wL : FWHM of the Lorentzian component
        
    Returns:
    y : array
        Voigt profile at each x position.
    """
    A, x0, wG, wL = p
    wG /= 2
    wL /= 2
    z = ((x - x0) + 1j * wL) / (wG * np.sqrt(2))
    return A * np.real(wofz(z)) / (wG * np.sqrt(2 * np.pi))