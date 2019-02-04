def gaussian(x,A,mu,sig):
    """This produces a gaussian function on the grid x with amplitude A, mean mu
    and standard deviation sig. Will need to expand it with a version that has
    a polynomial continuum in the same way that IDL does it."""

    import numpy as np
    return A * np.exp(-0.5*(x - mu)/sig*(x - mu)/sig)

def findgen(n):
    """This is basically IDL's findgen function.
    a = findgen(5) will return an array with 5 elements from 0 to 4:
    [0,1,2,3,4]
    """
    import numpy as np
    return np.linspace(0,n-1,n)
