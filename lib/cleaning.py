def construct_outlier_mask(order,w,c_thresh):
    print('THIS DOESNT EXIST YET')
    return(order*0.0+1.0)


def clean_ccf(rv,ccf,ccf_e,dp):
    """This routine normalizes the CCF fluxes and subtracts the average out of
    transit CCF, using the transit lightcurve as a mask."""

    import numpy as np
    import lib.functions as fun
    import lib.utils as ut
    from matplotlib import pyplot as plt
    import pdb
    import lib.system_parameters as sp
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
    ccf_n = ccf/meanblock.T
    ccf_ne = ccf_e/meanblock.T
    meanccf=np.mean(ccf_n[transit == 1.0,:],axis=0)
    meanblock=fun.rebinreform(meanccf,len(meanflux))
    ccf_nn = ccf_n-meanblock
    return(ccf_n,ccf_ne,ccf_nn)
