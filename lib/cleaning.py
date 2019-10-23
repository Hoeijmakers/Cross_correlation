

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
    transitblock = fun.rebinreform(transit,len(rv))


    meanflux=np.median(ccf,axis=1)#Normalize the baseline flux.
    meanflux_e=1.0/len(rv)*np.sqrt(np.nansum(ccf_e**2.0,axis=1))#1/N times sum of squares.
    #I validated that this is approximately equal to ccf_e/sqrt(N).
    meanblock=fun.rebinreform(meanflux,len(rv))
    meanblock_e=fun.rebinreform(meanflux_e,len(rv))

    ccf_n = ccf/meanblock.T
    ccf_ne = np.abs(ccf_n) * np.sqrt((ccf_e/ccf)**2.0 + (meanblock_e.T/meanblock.T)**2.0)#R=X/Z -> dR = R*sqrt( (dX/X)^2+(dZ/Z)^2 )
    #I validated that this is essentially equal to ccf_e/meanblock.T; as expected because the error on the mean spectrum is small compared to ccf_e.



    meanccf=np.nanmean(ccf_n[transit == 1.0,:],axis=0)
    meanccf_e=1.0/np.sum(transit==1)*np.sqrt(np.nansum(ccf_ne[transit == 1.0,:]**2.0,axis=0))#I validated that this is approximately equal
    #to sqrt(N)*ccf_ne, where N is the number of out-of-transit exposures.

    meanblock2=fun.rebinreform(meanccf,len(meanflux))
    meanblock2_e=fun.rebinreform(meanccf_e,len(meanflux))



    ccf_nn = ccf_n/meanblock2 -1.0#MAY NEED TO DO SUBTRACTION INSTEAD TOGETHER W. NORMALIZATION OF LIGHTCURVE. SEE ABOVE.
    ccf_nne = np.abs(ccf_n/meanblock2)*np.sqrt((ccf_ne/ccf_n)**2.0 + (meanblock2_e/meanblock2)**2.0)
    #I validated that this error is almost equal to ccf_ne/meanccf



    #ONLY WORKS IF LIGHTCURVE MODEL IS ACCURATE, i.e. if Euler observations are available.
    print("---> WARNING IN CLEANING.CLEAN_CCF(): NEED TO ADD A FUNCTION THAT YOU CAN NORMALIZE BY THE LIGHTCURVE AND SUBTRACT INSTEAD OF DIVISION!")
    return(ccf_n,ccf_ne,ccf_nn,ccf_nne)
