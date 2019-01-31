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

def typetest_list(varname,var,vartype):
    """This program tests the type of the elements in the list var which has the
    name varname, against the type vartype, and raises an exception if either
    varname is not a string, type(var) is not equal to list or the elements of
    var are not ALL equal to vartype.

    Example:
    a = 'ohai'
    utils.typtest_array('a',a,str)"""
    if isinstance(varname,str) != True:
        raise Exception("Input error in typetest: varname should be of type string.")
    if isinstance(var,list) != True:
        raise Exception("Input error in typetest_array: %s should be of class list." % varname)
    for i in range(0,len(var)):
        typetest('element %s of sizes' % i,var[i],vartype)


def vartest(var,ndim,sizes):
    import numpy as np
    typetest('ndim',ndim,int)
    typetest('sizes',sizes,list)


    dimerror=0.0
    sizeerror=0.0
    if np.dim(var) != ndim:
        raise Exception("Dimension error in vartest:  ")
