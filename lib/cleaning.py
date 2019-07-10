

def clean_ccf(rv,ccf,ccf_e,dp):
    """This routine normalizes the CCF fluxes and subtracts the average out of
    transit CCF, using the transit lightcurve as a mask."""

    import numpy as np
    import lib.functions as fun
    import lib.utils as ut
    from matplotlib import pyplot as plt
    import pdb
    import lib.system_parameters as sp
    import sys
    ut.typetest('rv in clean_ccf',rv,np.ndarray)
    ut.typetest('ccf in clean_ccf',ccf,np.ndarray)
    ut.typetest('ccf_e in clean_ccf',ccf_e,np.ndarray)
    ut.typetest('dp in clean_ccf',dp,str)
    ut.dimtest(ccf,[0,len(rv)])
    ut.dimtest(ccf_e,[0,len(rv)])
    ut.nantest('rv in clean_ccf',rv)
    ut.nantest('ccf in clean_ccf',ccf)
    ut.nantest('ccf_e in clean_ccf',ccf_e)


    transit=sp.transit(dp)



    meanflux=np.median(ccf,axis=1)#Normalize the baseline flux.
    meanblock=fun.rebinreform(meanflux,len(rv))
    transitblock = fun.rebinreform(transit,len(rv))
    ccf_n = ccf/meanblock.T#*transitblock.T

    ccf_ne = ccf_e/meanblock.T#*transitblock.T

    meanccf=np.mean(ccf_n[transit == 1.0,:],axis=0)
    meanblock2=fun.rebinreform(meanccf,len(meanflux))
    # ut.writefits('test.fits',meanblock)
    ccf_nn = ccf_n/meanblock2-1.0#MAY NEED TO DO SUBTRACTION INSTEAD TOGETHER W. NORMALIZATION OF LIGHTCURVE. SEE ABOVE.
    ut.save_stack('test.fits',[ccf,meanblock.T,ccf_n,meanblock2,ccf_nn,ccf_n-meanblock2,transitblock.T])

    #ONLY WORKS IF LIGHTCURVE MODEL IS ACCURATE, i.e. if Euler observations are available.
    print("---> WARNING IN CLEANING.CLEAN_CCF(): NEED TO ADD A FUNCTION THAT YOU CAN NORMALIZE BY THE LIGHTCURVE AND SUBTRACT INSTEAD OF DIVISION!")
    return(ccf_n,ccf_ne,ccf_nn)
