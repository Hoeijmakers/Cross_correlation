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

    import lib.utils as ut
    import lib.system_parameters as sp
    import lib.models
    import lib.constants as const
    import numpy as np
    #import matplotlib.pyplot as plt
    import scipy
    import lib.operations as ops
    #from astropy.io import fits
    #import pdb

    ut.dimtest(order,[0,len(wld)])
    ut.typetest('dp',dp,str)
    ut.typetest('modelname',dp,str)

    Rd=sp.paramget('resolution',dp)
    planet_radius=sp.paramget('Rp',dp)
    inclination=sp.paramget('inclination',dp)
    P=sp.paramget('P',dp)

    wlm,fxm=get_model(modelname)

    if wlm[-1] <= wlm[0]:#Reverse the wl axis if its sorted the wrong way.
        wlm=np.flipud(wlm)
        fxm=np.flipud(fxm)

    #With the model and the revelant parameters in hand, now only select that
    #part of the model that covers the wavelengths of the order provided.
    #A larger wavelength range would take much extra time because the convolution
    #is a slow operation.
    if np.min(wlm) > np.min(wld)-1.0 or np.max(wlm) < np.max(wld)+1.0:
        raise Exception('ERRROR in model injection: Data grid falls (partly) outside of model range.')

    modelsel=[(wlm >= np.min(wld)-1.0) & (wlm <= np.max(wld)+1.0)]

    wlm=wlm[modelsel]
    fxm=fxm[modelsel]

    shape=np.shape(order)
    n_exp=shape[0]
    transit=sp.transit(dp)
    rv=sp.RV(dp)
    dRV=sp.dRV(dp)*5.0
    phi=sp.phase(dp)
    ut.dimtest(transit,[n_exp])
    ut.dimtest(rv,[n_exp])
    ut.dimtest(phi,[n_exp])
    ut.dimtest(dRV,[n_exp])

    mask=(transit-1.0)/(np.min(transit-1.0))

    injection=order*0.0
    injection_rot_only=order*0.0
    injection_pure=order*0.0
    #t1=ut.start()
    fxm_b=ops.blur_rotate(wlm,fxm,(const.c/1000.0)/Rd,planet_radius,P,inclination)
    #t2=ut.end(t1)
    #print('...in blur rotate.')


    wlm_cv,fxm_bcv,vstep=ops.constant_velocity_wl_grid(wlm,fxm_b,oversampling=1.5)
    print('v_step is %s km/s' % vstep)
    print('So the exp-time blurkernel has an avg width of %s px.' % (np.mean(dRV)/vstep))
    #t3=ut.start()
    #test=fxm_cv*0.0
    #test[[1000,2000,3000,4000]]=1.0
    #fxm_b2=ops.smooth(fxm_cv,np.mean(dRV)/vstep,mode='box')
    #t4=ut.end(t3)

    #plt.plot(wlm_cv,test)
    #plt.plot(wlm_cv,fxm_b2)
    #plt.show()
    #print('...in smooth.')


    #t5=ut.start()
    #fxm_b2=ops.blur_spec(wlm,fxm_b,np.mean(dRV),mode='box')
    #t6=ut.end(t5)
    #print('...in blur_spec.')
    #THIS TAKES 1.4 SECONDS FOR 29k POINTS.
    #THAT IS MUCH TOO SLOW BECAUSE IT NEEDS TO BE DONE INSIDE THE FORLOOP.
    #WOULD TAKE OVER 1M PER ORDER.
    #INSTEAD, I SHOULD AT THIS POINT PROBABLY CHOOSE TO INTERPOLATE THE MODEL
    #TO AN OVERSAMPLED CONSTANT-VELOCITY WL GRID, AND USE A BUILTIN CONVOLUTION
    #FUNCTION WITH A CONSTANT-SIZED BOX. THE DOWNSIDE OF THE BUILTIN CONVOLUTIOn
    #OPERATOR IS THAT IT SCREWS UP THE EDGES, BUT THE MODEL IS CROPPED QUITE
    #WIDELY AROUND THE DATA ANYWAY...


    for i in range(0,n_exp):
        fxm_b2=ops.smooth(fxm_bcv,dRV[i]/vstep,mode='box')
        shift=(1.0+rv[i]/(const.c/1000.0))
        fxm_i=scipy.interpolate.interp1d(wlm_cv*shift,fxm_b2) #This is a class that can be called.
        injection[i,:]=fxm_i(wld)
        injection_rot_only[i,:]=scipy.interpolate.interp1d(wlm*shift,fxm_b)(wld)
        injection_pure[i,:]=scipy.interpolate.interp1d(wlm*shift,fxm)(wld)
    #fits.writeto('test.fits',injection,overwrite=True)

    #plt.plot(wld,injection_pure[15,:])
    #plt.plot(wld,injection_rot_only[15,:])
    #plt.plot(wld,injection[15,:])
    #plt.show()
    #pdb.set_trace()
    return(injection*order)




def read_binary_mask(inpath,R=1000000.0):
    """This short function reads in a Geneva-style binary mask file and converts
    it to a highly-sampled FITS-spectrum (IS THAT WHAT YOU WANT??) (Only if the
    shifting algorithm that I will use for CCVs later conserves flux the right way...).
    (However it will be more accurate to do the CCV based on the true edges of the
    mask on the sampling of the data itself; rather than first converting to a
    spectrum and then converting to another spectrum. Would be a keyword in the
    CCV routine. However, for now I keep it like this, so that I can plot the template)
    Provide the input path of the mask file and the output file of the FITS file
    and optionally the spectral resolution of the sampling."""
    import numpy as np
    import pylab
    import pdb
    from lib.utils import typetest
    typetest('inpath',inpath,str)

    try:
        f = open(inpath, 'r')
    except FileNotFoundError:
        raise Exception('Configfile does not exist at %s' % inpath) from None
    x = f.read().splitlines()
    f.close()
    n_lines=len(x)
    wl=np.linspace(350.0,1100.0,(1100-350)*R/350.0)
    fx=wl*0.0
    for i in x:
        line=i.split()
        start=float(line[0])/10.0
        end=float(line[1])/10.0
        weight=float(line[2])
        sel=np.where((wl >= start) & (wl <= end))[0]
        fx[sel]=weight

    return wl,fx
