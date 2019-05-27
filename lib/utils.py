#This package contains utility functions that are mostly used for making my life
#easier, such as a wrapper for measuring elapsed time, all sorts of variable tests,
#and quickly saving fits files.

def start():
    import time
    return(time.time())

def end(start,id=''):
    import time
    end=time.time()
    print('Elapsed %s: %s' % ('on timer '+id,end-start))
    return end-start

def nantest(varname,var):
    import numpy as np
    if np.isnan(var).any()  == True:
        raise Exception("NaN error: %s contains NaNs." % varname)
    if np.isinf(var).any()  == True:
        raise Exception("Finite error: %s contains in-finite values." % varname)

def postest(a,varname=''):
    """This function tests whether a number/array is strictly positive."""
    import numpy as np
    if np.min(a) <= 0:
        raise Exception('POSTEST ERROR: Variable %s should be strictly positive' % varname)

def notnegativetest(a,varname=''):
    """This function tests whether a number/array is strictly positive."""
    import numpy as np
    if np.min(a) < 0:
        raise Exception('POSTEST ERROR: Variable %s should not be negative' % varname)



def typetest(varname,var,vartype):
    """This program tests the type of var which has the name varname against
    the type vartype, and raises an exception if either varname is not a string,
    or if type(var) is not equal to vartype.

    Example:
    a = 'ohai'
    utils.typtest('a',a,str)"""
    if isinstance(varname,str) != True:
        raise Exception("Input error in typetest: varname should be of type string.")
    if isinstance(var,vartype) != True:
        raise Exception("Type error: %s should be %s." % (varname,vartype))

def typetest_array(varname,var,vartype):
    """This program tests the type of the elements in the array or list var which has the
    name varname, against the type vartype, and raises an exception if either
    varname is not a string, type(var) is not equal to numpy.array or list, or the elements of
    var are not ALL of a type equal to vartype.

    Example:
    a = ['alef','lam','mim']
    utils.typetest_array('teststring',a,str)"""
    #NEED TO FIX: MAKE SURE THAT A RANGE OF TYPES CAN BE TESTED FOR, SUCH AS
    #float, np.float32, np.float64... should all pass as a float.
    import numpy as np
    if isinstance(varname,str) != True:
        raise Exception("Input error in typetest: varname should be of type string.")
    if (isinstance(var,list) != True) and (isinstance(var,np.ndarray) != True):
        raise Exception("Input error in typetest_array: %s should be of class list or numpy array." % varname)
    for i in range(0,len(var)):
        typetest('element %s of %s' % (i,varname),var[i],vartype)


def dimtest(var,sizes):

    """This program tests the dimensions and shape of the input array var.
    Sizes is the number of elements on each axis.
    The program uses the above type tests to make sure that the input is ok.
    If an element in sizes is set to zero, that dimension is not checked against.
    Example:
    import numpy as np
    a=[[1,2,3],[4,3,9]]
    b=np.array(a)
    dimtest(a,2,[2,3])
    dimtest(a,2,[3,10])
    """
    import numpy as np
    typetest_array('sizes',sizes,int)
    ndim=len(sizes)

    dimerror=0.0
    sizeerror=0.0
    if np.ndim(var) != ndim:
        raise Exception("Dimension error in vartest:  ndim = %s but was required to be %s." % (np.ndim(var),ndim))

    sizes_var=np.shape(var)

    for i in range(0,len(sizes)):
        if sizes[i] < 0:
            raise Exception("Sizes was not set correctly. It contains negative values. (%s)" % sizes(i))
        if sizes[i] > 0:
            if sizes[i] != sizes_var[i]:
                raise Exception("Dimension error in vartest: Axis %s contains %s elements, but %s were required." % (i,sizes_var[i],sizes[i]))

def save_stack(filename,list_of_2D_frames):
    """This code saves a stack of fits-files to a 3D cube, that you can play
    through in DS9. For diagnostic purposes."""
    import astropy.io.fits as fits
    import numpy as np
    base = np.shape(list_of_2D_frames[0])
    N = len(list_of_2D_frames)
    out = np.zeros((base[0],base[1],N))
    for i in range(N):
        out[:,:,i] = list_of_2D_frames[i]
    fits.writeto(filename,np.swapaxes(np.swapaxes(out,2,0),1,2),overwrite=True)

def writefits(filename,array):
    """This is a fast wrapper for fits.writeto, with overwrite enabled.... ! >_<"""
    import astropy.io.fits as fits
    fits.writeto(filename,array,overwrite=True)
