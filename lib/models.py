def get_model(name,library='models/library'):
    """This program queries a model from a library file, with predefined models
    for use in model injection, cross-correlation and plotting. These models have
    a standard format. They are 2 rows by N points, where N corresponds to the
    number of wavelength samples. The first row contains the wavelengths, the
    second the associated flux values. The library file has two columns that are
    formatted as follows:

    modelname  modelpath

    modelpath starts in the models/ subdirectory.

    Example call:
    wl,fx = get_model('WASP-121_TiO',library='models/testlib')
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
