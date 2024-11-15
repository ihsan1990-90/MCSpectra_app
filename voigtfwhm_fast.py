import numpy as np
from scipy.special import wofz

def voigtfwhm_fast(x, p):
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
    A0, X, Y = p
    A = [-1.2150, -1.3509, -1.2150, -1.3509]
    B = [1.2359, 0.3786, -1.2359, -0.3786]
    C = [-0.3085, 0.5906, -0.3085, 0.5906]
    D = [0.0210, -1.1858, -0.0210, 1.1858]

    V = np.zeros(len(x))
    for i in [0,1,2,3]:
        V += np.divide((C[i]*(Y-A[i]) + D[i]*(X-B[i])),((Y-A[i])**2+np.power((X-B[i]),2)))
    
    return A0*V