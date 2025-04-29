#This is just the sig figs functions. 

import numpy as np 

def sig_figs(x: float, precision: int):
    """
    Rounds a number to number of significant figures
    Parameters:
    - x - the number to be rounded
    - precision (integer) - the number of significant figures
    Returns:
    - float
    """

    x = float(x)
    precision = int(precision)
    if x > 0: 
        return np.round(x, -int(np.floor(np.log10(abs(x)))) + (precision - 1))
    elif x < 0: 
        return -np.round(abs(x), -int(np.floor(np.log10(abs(x)))) + (precision - 1))
    else: 
        return x 
