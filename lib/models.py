def get_model(name,library='models/library'):
    """This program queries a model from a library file, with predefined models
    for use in model injection, cross-correlation and plotting. These models have
    a standard format. They are 2 rows by N points, where N corresponds to the
    number of wavelength samples. The first row contains the wavelengths, the
    second the associated flux values. The library file has two columns that are
    formatted as follows:

    modelname  modelpath
    modelname  modelpath
    modelname  modelpath

    modelpath starts in the models/ subdirectory.

    Example call:
    wlm,fxm = get_model('WASP-121_TiO',library='models/testlib')
    """

    from lib.utils import typetest
    from lib.utils import dimtest
    import pdb
    from astropy.io import fits
    typetest('name',name,str)
    typetest('library',library,str)

    #First try to open the library file.
    try:
        f = open(library, 'r')
    except FileNotFoundError:
        raise Exception('Model library does not exist at %s' % library) from None
    x = f.read().splitlines()#Read everything into a big string array.
    f.close()
    n_lines=len(x)#Number of models in the library.
    models={}#This will contain the model names.

    for i in range(0,n_lines):
        line=x[i].split()
        value=(line[1])
        models[line[0]] = value

    try:
        modelpath='models/'+models[name]
    except KeyError:
        raise Exception('Model %s is not present in library at %s' % (name,library)) from None


    modelarray=fits.getdata(modelpath)#Two-dimensional array with wl on the first row and flux on the second row.
    dimtest(modelarray,[2,0])
    return(modelarray[0,:],modelarray[1,:])



def inject_model(wld,order,dp,modelname):
    """This function takes a spectral order and injects a model with library
    identifier modelname, and system parameters as defined in dp.

    It returns a copy of order with the model injected."""
    #I chose not to save
    import lib.system_parameters as sp
    import lib.models
    import lib.constants as const
    import numpy as np
    Rd=sp.paramget('resolution',dp)
    planet_radius=sp.paramget('Rp',dp)
    inclination=sp.paramget('inclination',dp)
    P=sp.paramget('P',dp)

    wlm,fxm=get_model('modelname')

    shape=np.shape(order)
    n_exp=shape[1]

    mask=sp.transit(dp)

    #Resume with blurring. First need to define a function for resolution_kernel.
    #Well actually, I want to do blurring with a dynamic kernel; such that I no
    #longer make the assumption of a small wavelength range - without the need to
    #interpolate onto a  constant velocity grid. Need to google this. Or code it.
    #It can be done in a single forloop.

    #And, either I need to do this blurring (which is slow) only once per dataset
    #(meaning, outside of this routine).
    #Or, I do it only for the chunk that is injected. I.e. selecting the model
    #just around the wavelength chunk.
