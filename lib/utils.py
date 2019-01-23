def typetest(varname,var,vartype):
    if isinstance(varname,str) != True:
        raise Exception("Unit error in unit test: varname should be of class string.")
    if isinstance(var,vartype) != True:
        raise Exception("Unit error: %s should be of class string." % varname)
