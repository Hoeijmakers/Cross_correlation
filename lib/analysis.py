


def xcor(list_of_wls,list_of_orders,wlm,fxm,drv,RVrange,NaN=None,verticalmask=None,plot=False,list_of_errors=None):
    """This routine takes a combined dataset (in the form of lists of wl spaces,
    spectral orders and possible a matching list of errors on those spectal orders),
    as well as a template (wlm,fxm) to cross-correlate with, and the cross-correlation
    parameters (drv,RVrange). The code takes on the order of ~10 minutes for an entire
    HARPS dataset, which appears to be superior than my IDL pipe.

    The CCF used is the Geneva-style weighted average; not the Pearson CCF. Therefore
    it measures true 'average' planet lines, with flux on the y-axis of the CCF.
    The template must therefore be (something close to) a binary mask, with values
    inside spectral lines (the CCF is scale-invariant so their overall scaling
    doesn't matter),

    It returns the RV axis and the resulting CCF in a tuple, """
    #NOW I NEED TO SIT DOWN AND DO THE MATH OF THE CCF AGAIN.
    #SPECIFICALLY I NEED TO DETERMINE AGAIN HOW I AM GOING TO DEAL WITH DIFFERENT
    #SPECTRAL ORDERS. I wrote some beautiful math in my paper. Maybe I should try
    #To reproduce it.
    import lib.functions as fun
    import lib.constants as const
    import lib.utils as ut
    import numpy as np
    #import matplotlib.pyplot as plt
    import scipy.interpolate
    import pdb
    import astropy.io.fits as fits
    import matplotlib.pyplot as plt


#===FIRST ALL SORTS OF TESTS ON THE INPUT===
    if len(list_of_wls) != len(list_of_orders):
        raise Exception('ERROR in xcor: List of wls and list of orders have different length (%s & %s).' % (len(list_of_wls),len(list_of_orders)))
    ut.dimtest(wlm,[len(fxm)])
    ut.typetest('wlm in xcor',wlm,np.ndarray)
    ut.typetest('fxm in xcor',fxm,np.ndarray)
    ut.typetest('drv in xcor',drv,float)
    ut.typetest('RVrange in xcor',RVrange,float)
    ut.postest(RVrange,varname='RVrange in xcor')
    ut.postest(drv,varname='drv in xcor')
    ut.nantest('fxm in xcor',fxm)
    N=len(list_of_wls)#Number of orders.
    if np.ndim(list_of_orders[0]) == 1.0:
        n_exp=1
    else:
        n_exp=len(list_of_orders[0][:,0])#Number of exposures.
#===Then check that all orders indeed have n_exp exposures===
        for i in range(N):
            if len(list_of_orders[i][:,0]) != n_exp:
                raise Exception('ERROR in xcor: Not all orders have %s exposures.' % n_exp)
#===END OF TESTS. NOW DEFINE CONSTANTS===
    c=const.c/1000.0
    RV=fun.findgen(2.0*RVrange/drv+1)*drv-RVrange#..... CONTINUE TO DEFINE THE VELOCITY GRID
    shift=1.0+RV/c#Array along which the template will be shifted during xcor.



#===Define the output CCF array.
    CCF = np.zeros((n_exp,len(shift)))#*float('NaN')
    CCF_E = CCF*0.0
    Tsums = fun.findgen(len(shift))*0.0*float('NaN')
#===Then comes the big forloop.
    #The outer loop is the shifts. For each, we loop through the orders.


    counter = 0
    for i in range(len(shift)):
        T_sum = 0.0
        wlms=wlm*shift[i]
        for j in range(N):
            # if j == 0:
            #     CCF[:,i]=0.0
            wl=list_of_wls[j]
            order=list_of_orders[j]

            T_i=scipy.interpolate.interp1d(wlms[(wlms >= np.min(wl)-10.0) & (wlms <= np.max(wl)+10.0)],fxm[(wlms >= np.min(wl)-10.0) & (wlms <= np.max(wl)+10.0)])
            T=T_i(wl)#Interpolated onto wl
            T_matrix=fun.rebinreform(T_i(wl),n_exp)
            CCF[:,i]+=np.nansum(T_matrix*order,1)
            if list_of_errors != None:
                sigma=list_of_errors[j]
                CCF_E[:,i]+=np.nansum((T_matrix*sigma)**2.0,1)#SIGMA SEEMS OVERESTIMATED BY A FACTOR OF 10 DEBUG ON THE SPECTRA THEMSELVES.
            T_sum+=np.sum(np.abs(T))#ARE THESE TWO CORRECT? CHECK IT! DO THE MATH
            #Will do when the error propagation comes in.
        CCF[:,i] /= T_sum
        CCF_E[:,i] /= T_sum**2.0
        Tsums[i] = T_sum
        T_sum = 0.0
        counter += 1
        print('  '+f"{i/(len(shift)-1)*100:.1f} %", end="\r")#Statusbar.
    ut.nantest('CCF in Xcor',CCF)
    ut.nantest('CCF_E in Xcor',CCF_E)


    if plot == True:
        fig, (a0,a1,a2) = plt.subplots(3,1,gridspec_kw = {'height_ratios':[1,1,1]},figsize=(10,7))
        a02 = a0.twinx()
        for i in range(N):
            meanspec=np.nanmean(list_of_orders[i],axis=0)
            meanwl=list_of_wls[i]
            T_i=scipy.interpolate.interp1d(wlm[(wlm >= np.min(meanwl)-0.0) & (wlm <= np.max(meanwl)+0.0)],fxm[(wlm >= np.min(meanwl)-0.0) & (wlm <= np.max(meanwl)+0.0)],fill_value='extrapolate')
            T = T_i(meanwl)
        # meanspec-=min(meanspec)
        # meanspec/=np.max(meanspec)
        # T_plot-=np.min(T)
        # T_plot/=np.median(T_plot)
            a02.plot(meanwl,T,color='orange',alpha=0.3)
            a0.plot(meanwl,meanspec,alpha=0.5)
        a1.plot(RV,Tsums,'.')
        if list_of_errors != None:
            a2.errorbar(RV,np.mean(CCF,axis=0),fmt='.',yerr=np.mean(np.sqrt(CCF_E),axis=0)/np.sqrt(n_exp),zorder=3)
        else:
            a2.plot(RV,np.mean(CCF,axis=0),'.')
        a0.set_title('t-averaged data and template')
        a1.set_title('Sum(T)')
        a2.set_title('t-averaged CCF')
        a0.tick_params(axis='both',labelsize='5')
        a02.tick_params(axis='both',labelsize='5')
        a1.tick_params(axis='both',labelsize='5')
        a2.tick_params(axis='both',labelsize='5')
        fig.tight_layout()
        plt.show()
                    # a0.set_xlim((674.5,675.5))
                    # a1.set_xlim((674.5,675.5))
    if list_of_errors != None:
        return(RV,CCF,np.sqrt(CCF_E))
    return(RV,CCF)
    #plt.plot(RV,CCF[10,:])
    #plt.show()




#Kelt-9 radius ratio:
#Rp/Rs = df = 0.08228 \pm 0.00043
#Rs = 2.196 \pm 0.012 Rssun
#Rp = df*Rs

#Propagate errors:
#S_Rp^2 = (dRp/ddf sdf)^2 + (dRp/dRs sRs)^2



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



def plot_ccf(rv,ccf,dp):
    """This is a routine that does all the plotting of the two-dimensional CCF.
    I expect this to be an organic function that is adapted to my plotting needs."""

    import numpy as np
    import matplotlib.pyplot as plt
    import pdb

    plt.imshow(ccf)

    pdb.set_trace()
