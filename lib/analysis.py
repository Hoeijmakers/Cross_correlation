def xcor(list_of_wls,list_of_orders,wlm,fxm,drv,RVrange,NaN=None,verticalmask=None):

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

    t1=ut.start()
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
#===END OF TESTS. NOW DEFINE CONSTANTS===
    c=const.c/1000.0
    N=len(list_of_wls)#Number of orders.
    RV=fun.findgen(2.0*RVrange/drv+1)*drv-RVrange#..... CONTINUE TO DEFINE THE VELOCITY GRID
    shift=1.0+RV/c#Array along which the template will be shifted during xcor.
    n_exp=len(list_of_orders[0][:,0])#Number of exposures.

#===Then check that all orders indeed have n_exp exposures===
    for i in range(N):
        if len(list_of_orders[i][:,0]) != n_exp:
            raise Exception('ERROR in xcor: Not all orders have %s exposures.' % n_exp)

#===Define the output CCF array.
    CCF = np.zeros((n_exp,len(shift)))

#===Then comes the big forloop.
    #The outer loop is the shifts. For each, we loop through the orders.
    for i in range(len(shift)):
        T_sum = 0.0
        wlms=wlm*shift[i]
        for j in range(N):
            wl=list_of_wls[j]
            order=list_of_orders[j]
            T_i=scipy.interpolate.interp1d(wlms[(wlms >= np.min(wl)-1.0) & (wlms <= np.max(wl)+1.0)],fxm[(wlms >= np.min(wl)-1.0) & (wlms <= np.max(wl)+1.0)])
            T=T_i(wl)#Interpolated onto wl
            T_matrix=fun.rebinreform(T_i(wl),n_exp)
            CCF[:,i]+=np.sum(T_matrix*order,1)
            T_sum+=np.sum(T)#ARE THESE TWO CORRECT? CHECK IT! DO THE MATH
        CCF[:,i] /= T_sum
        T_sum = 0.0
    t2=ut.end(t1)
    plt.plot(RV,CCF[10,:])
    plt.show()
    pdb.set_trace()
            #ETC ETC ETC




#Kelt-9 radius ratio:
#Rp/Rs = df = 0.08228 \pm 0.00043
#Rs = 2.196 \pm 0.012 Rssun
#Rp = df*Rs

#Propagate errors:
#S_Rp^2 = (dRp/ddf sdf)^2 + (dRp/dRs sRs)^2
