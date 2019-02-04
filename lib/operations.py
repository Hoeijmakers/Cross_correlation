#This package contains operations that act on spectra, such as blurring with
#various kernels and continuum normalization.
def blur_gaussian_lsf(wl,order,v):
    """This function takes a spectrum, and blurs it using a Gaussian kernel, which
    has a FWHM width of v km/s everywhere. Meaning that its width changes dynamically.
    Because the kernel needs to be recomputed on each element of the wavelength axis
    individually, this operation is (or should be) much slower than convolution with
    a constant kernel, in which a simple shifting of the array, rather than a recomputation
    of the Gaussian is sufficient.

    This program is repeated below for a rotation kernel as well (and may be depricated
    by it...)"""
  import numpy as np
