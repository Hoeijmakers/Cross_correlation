#This pacakge contains routines used for analyzing CCFs and KpVsys diagrams.
#Its main routine in the cross-correlation function.
#The rest is construction of KpVsys and plotting routines.

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
    import sys

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
                print('ERROR in xcor: Not all orders have %s exposures.' % n_exp)
                sys.exit()
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

            T_i=scipy.interpolate.interp1d(wlms[(wlms >= np.min(wl)-10.0) & (wlms <= np.max(wl)+10.0)],fxm[(wlms >= np.min(wl)-10.0) & (wlms <= np.max(wl)+10.0)],bounds_error=False,fill_value=0.0)
            T = T_i(wl)
            T_matrix=fun.rebinreform(T,n_exp)
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
        ut.statusbar(i,shift)
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
        return(RV,CCF,np.sqrt(CCF_E),Tsums)
    return(RV,CCF,Tsums)




def xcor_new(list_of_wls,list_of_orders,wlm,fxm,drv,RVrange,NaN=None,verticalmask=None,plot=False,list_of_errors=None):
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

    It returns the RV axis and the resulting CCF in a tuple,

    THIS VERSION OF XCOR IS PROBLEMATIC DUE TO HOW I TRIED TO DEAL WITH NaNs IN THE INDIVIDUAL
    SPECTRA. NaNs ARE PROBLEMATIC BECAUSE THEY CHANGE THE AVERAGE. SO I TRIED TO ACCOUNT
    FOR MISSING DATAPOINTS IN THE CALCULATION OF T_SUM AS WELL (DENOMINATOR). HOWEVER,
    ACCOUNTING FOR NaNs WILL IN ALL CASES AFFECT THE MEASUREMENT OF THE AVERAGE SPECTRAL
    LINE AS A FUNCTION OF TIME; LEADING TO RESIDUALS WHEN CLEANING THE CCF.
    THIS HAPPENS ALWAYS WHEN THERE IS ONE OR MULTIPLE NaNs IN A COLUMN, UNLESS THE
    ENTIRE COLUMN IS NaNned OUT.

    THERE ARE TWO WAYS TO OVERCOME THIS.
    THE FIRST IS TO DO ALL CLEANING IN SPECTRAL SPACE RATHER THAN CCF SPACE.
    NaNs WILL THEN NOT BE PROBLEMATIC BECAUSE THE TIME-AVERAGE OF THE CCF IS NEVER USED.
    ONLY WILL NaNning OUT A CERTAIN PIXEL IN ONE SPECTRUM AFFEC THE AVERAGE OF THE
    MEASURED AVG LINE OF THE PLANET, BUT THIS EFFECT IS TINY, ESPECIALLY IN THE CASE
    OF MANY LINES. THE ONLY REASON IT POPS UP NOW IS BECAUSE ITS A TINY EFFECT ON
    VERY STRONG STELLAR LINES THAT I'M TRYING TO REMOVE. THE DOWNSIDE OF THIS IS THAT
    THE NOISE IS NOT CONSTANT; SO I WILL HAVE TO TRY TO DEAL WITH THIS (NORMALIZE BY
    STDDEV? WHAT TO DO ABOUT THE TEMPLATE THEN?).

    THE SECOND IS TO INTERPOLATE INSTEAD OF PUTTING NaNs. This is the easiest...

    UNTIL I HAVE DONE EITHER OF THESE, I WILL KEEP SIGMA-CLIPPING SWITCHED OFF.

    May need to write a small simulator to show how the SNR of the retrieved planet
    behaves as a function of choice."""
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
    import sys

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
                print('ERROR in xcor: Not all orders have %s exposures.' % n_exp)
                sys.exit()
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
        T_sum = fun.findgen(n_exp)*0.0
        wlms=wlm*shift[i]
        for j in range(N):
            # if j == 0:
            #     CCF[:,i]=0.0
            wl=list_of_wls[j]
            order=list_of_orders[j]
            mask = order*0.0+1.0
            T_i=scipy.interpolate.interp1d(wlms[(wlms >= np.min(wl)-10.0) & (wlms <= np.max(wl)+10.0)],fxm[(wlms >= np.min(wl)-10.0) & (wlms <= np.max(wl)+10.0)])
            T=T_i(wl)#Interpolated onto wl
            T_matrix=fun.rebinreform(T_i(wl),n_exp)
            CCF[:,i]+=np.nansum(T_matrix*order,1)
            if list_of_errors != None:
                sigma=list_of_errors[j]
                CCF_E[:,i]+=np.nansum((T_matrix*sigma)**2.0,1)#SIGMA SEEMS OVERESTIMATED BY A FACTOR OF 10 DEBUG ON THE SPECTRA THEMSELVES.
            T_sum+=np.nansum(np.abs(T_matrix*mask),axis=1)#ARE THESE TWO CORRECT? CHECK IT! DO THE MATH
            #Will do when the error propagation comes in.

        CCF[:,i] /= T_sum
        CCF_E[:,i] /= T_sum**2.0
        Tsums[i] = np.nanmean(T_sum)
        counter += 1
        ut.statusbar(i,shift)
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
        return(RV,CCF,np.sqrt(CCF_E),Tsums)
    return(RV,CCF,Tsums)
    #plt.plot(RV,CCF[10,:])
    #plt.show()

def plot_XCOR(list_of_wls,list_of_orders,wlm,fxm,RV,CCF,Tsums,dp,CCF_E=None,dv=12.0):
    import lib.functions as fun
    import lib.constants as const
    import lib.utils as ut
    import numpy as np
    #import matplotlib.pyplot as plt
    import scipy.interpolate
    import pdb
    import astropy.io.fits as fits
    import matplotlib.pyplot as plt
    import sys
    import lib.system_parameters as sp
    import lib.operations as ops

    if len(list_of_wls) != len(list_of_orders):
        raise Exception('ERROR in plot_XCOR: List of wls and list of orders have different length (%s & %s).' % (len(list_of_wls),len(list_of_orders)))

    ut.dimtest(wlm,[len(fxm)])
    ut.typetest('wlm in plot_xcor',wlm,np.ndarray)
    ut.typetest('fxm in plot_xcor',fxm,np.ndarray)
    ut.nantest('fxm in xcor',fxm)
    N=len(list_of_wls)#Number of orders.

    if np.ndim(list_of_orders[0]) == 1.0:
        n_exp=1
    else:
        n_exp=len(list_of_orders[0][:,0])#Number of exposures.
#===Then check that all orders indeed have n_exp exposures===
        for i in range(N):
            if len(list_of_orders[i][:,0]) != n_exp:
                print('ERROR in xcor: Not all orders have %s exposures.' % n_exp)
                sys.exit()


    # fig, (a0,a2) = plt.subplots(2,1,gridspec_kw = {'height_ratios':[1.0/3,2.0/3]},figsize=(10,7))
    fig=plt.figure(figsize=(10,7))
    a0 = plt.axes([0.05, 0.6, 0.9, 0.35])
    a2 = plt.axes([0.05, 0.05, 0.4, 0.5])
    a3 = plt.axes([0.55, 0.05, 0.4, 0.5])
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
    # a1.plot(RV,Tsums,'.')

    transit = sp.transit(dp)
    sel = (transit == 1.0)
    n_exp_out_of_transit = len(np.where(transit == 1)[0])


    for i in range(0,n_exp):
        a2.plot(RV,CCF[i,:],'.',color='gray',alpha=0.2)

    CCF1D = np.mean(CCF[sel,:],axis=0)
    if len(np.shape(CCF_E)) > 0:
        a2.errorbar(RV,CCF1D,fmt='.',yerr=np.mean(np.sqrt(CCF_E),axis=0)/np.sqrt(n_exp_out_of_transit),zorder=3)
    else:
        a2.plot(RV,CCF1D,'.')




    maxindex = (abs(CCF1D-np.nanmedian(CCF1D))).argmax()
    RV_0 = RV[maxindex]
    sel = [(RV > RV_0-dv) & (RV < RV_0+dv)]
    y = CCF1D
    fit=ops.gauss_fit(RV[tuple(sel)],y[tuple(sel)],start=[0.0,RV_0,dv,np.nanmedian(y)])
    a2.plot(RV,y,'.')
    a2.plot(RV[tuple(sel)],y[tuple(sel)],'.',color='red')#Using tuples here to avoid a deprication warning.
    a2.plot(RV[tuple(sel)],fun.gaussian(RV[tuple(sel)],fit[0],fit[1],fit[2],cont=fit[3]))
    a2.text(fit[1]+dv/2,fit[3]+fit[0],'RV = %s km/s' % round(fit[1],3),color='green')
    a0.set_title('t-averaged data and template')
    # a1.set_title('Sum(T)')
    a2.set_title('t-averaged CCF')
    a0.tick_params(axis='both',labelsize='5')
    a02.tick_params(axis='both',labelsize='5')
    # a1.tick_params(axis='both',labelsize='5')
    a2.tick_params(axis='both',labelsize='5')
    # fig.tight_layout()
    plt.show()


def plot_RV_star(dp,RV,CCF2D,RVrange=[]):
    """This program fits the RV of a ccf in a 2D ccf.  May not be useful anymore,
    #only for diagnostic?"""
    import lib.operations as ops
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import lib.system_parameters as sp
    import pdb
    rv_c = []
    n_exp = np.shape(CCF2D)[0]
    if len(RVrange) > 0:
        sel = np.squeeze([(RV > RVrange[0]) & (RV < RVrange[1])])
        RV=RV[sel]
        CCF2D = CCF2D[:,sel]
    for i in range(n_exp):
        rv_c = np.append(rv_c,ops.measure_rv(RV,CCF2D[i,:]))
    rv_K = sp.RV_star(dp)
    transit = sp.transit(dp)
    phase = sp.phase(dp)
    berv = sp.berv(dp)
    phase[phase>0.5]-=1.0
    sel = [transit == 1.0]
    plt.plot(phase[sel],rv_c[sel],'.')
    plt.plot(phase[sel],rv_c[sel]+berv[sel]-rv_K[sel],'.')
    plt.show()


def construct_KpVsys(rv,ccf,ccf_e,dp,kprange=[0,300],dkp=1.0):
    """The name says it all. Do good tests."""
    import lib.functions as fun
    import lib.operations as ops
    import numpy as np
    import lib.system_parameters as sp
    import matplotlib.pyplot as plt
    import astropy.io.fits as fits
    import lib.utils as ut
    import sys
    import pdb
    Kp = fun.findgen((kprange[1]-kprange[0])/dkp+1)*dkp+kprange[0]
    n_exp = np.shape(ccf)[0]
    KpVsys = np.zeros((len(Kp),len(rv)))
    KpVsys_e = np.zeros((len(Kp),len(rv)))
    transit = sp.transit(dp)-1.0
    transit /= np.nansum(transit)
    transitblock = fun.rebinreform(transit,len(rv)).T

    j = 0
    ccfs = []
    for i in Kp:
        dRV = sp.RV(dp,vorb=i)*(-1.0)
        ccf_shifted = ops.shift_ccf(rv,ccf,dRV)
        ccf_e_shifted = ops.shift_ccf(rv,ccf_e,dRV)
        ccfs.append(ccf_shifted)
        KpVsys[j,:] = np.nansum(transitblock * ccf_shifted,axis=0)
        KpVsys_e[j,:] = (np.nansum((transitblock*ccf_e_shifted)**2.0,axis=0))**0.5
        # plt.plot(rv,KpVsys[j,:])
        # plt.fill_between(rv, KpVsys[j,:]-KpVsys_e[j,:], KpVsys[j,:]+KpVsys_e[j,:],alpha=0.5)
        # plt.show()
        # pdb.set_trace()
        j+=1
        ut.statusbar(i,Kp)
    return(Kp,KpVsys,KpVsys_e)





def plotting_scales_2D(rv,y,ccf,xrange,yrange,Nxticks=10.0,Nyticks=10.0,nsigma=3.0):
    """This is a primer for plot_ccf, which defines the plotting ranges and colour
    scale of a 2D CCF. It can probably be used generally for any 2D image."""
    import numpy as np
    import lib.functions as fun
    import math
    import pdb
    #Define the plotting range of the images.
    nrv=len(rv)
    ny=len(y)
    drv=rv[1]-rv[0]
    dy=y[1]-y[0]
    sel = ((rv >= xrange[0]) & (rv <=xrange[1]))
    ccf_sub = ccf[:,sel]
    sely = ((y >= yrange[0]) & (y <= yrange[1]))
    ccf_sub = ccf_sub[sely,:]
    vmin,vmax = fun.sigma_clip(ccf_sub,nsigma=nsigma)#Sigma clipping.
    rvmin=rv[sel].min()
    rvmax=rv[sel].max()
    ymin=y[sely].min()
    ymax=y[sely].max()
    #This initiates the meshgrid onto which the 2D image is defined.
    xmin = rvmin; xmax = rvmax; dx = drv#Generalization away from rv.
    x2,y2 = np.meshgrid(np.arange(xmin,xmax+dx+dx,dx)-dx/2.,np.arange(ymin,ymax+dy+dy,dy)-dy/2.)

    dxt = (xmax+dx - xmin) / Nxticks
    dyt = (ymax+dy - ymin) / Nyticks

    ndigits_dxt= -1.0*(min([math.floor(np.log(dxt)),0]))#Set the rounding number of digits to either 0 (ie integers)
    ndigits_dyt= -1.0*(min([math.floor(np.log(dyt)),0]))#or a power of 10 smaller than that, otherwise.

    xticks = np.arange(xmin,xmax+dx,round(dxt,int(ndigits_dxt)))
    yticks = np.arange(ymin,ymax+dy,round(dyt,int(ndigits_dyt)))
    return(x2,y2,ccf_sub,rv[sel],y[sely],xticks,yticks,vmin,vmax)





def plot_ccf(rv,ccf,dp,xrange=[-200,200],yrange=[0,0],Nxticks=10.0,Nyticks=10.0,title='',doppler_model = [],i_legend=True,show=True):
    """This is a routine that does all the plotting of the cleaned 2D CCF.
    I expect this to be an organic function that is adapted to my plotting needs.
    THIS CURRENTLY HAS CORE USAGE IN RUN_INSTANCE. DO GOOD TESTS AND WRITE GOOD DOCUMENTATION

    If yrange is not set, it defaults to the entire y axis (see below).

    STILL NEED TO ADD FUNCTIONALITY FOR NXTICKS AND NYTICKS SEPARATELY.
    MAKE SURE THAT XRANGE AND YRANGE BEHAVE THE SAME.
    INSTEAD OF DOPPLER_MODEL=[], ADD ARBITRARY LINES TO OVERPLOT.. with labels.
    .. and tests on those."""

    import numpy as np
    import matplotlib.pyplot as plt
    import pdb
    import lib.drag_colour as dcb
    import lib.functions as fun
    import pylab as pl
    import lib.system_parameters as sp
    import lib.plotting as fancyplots
    import lib.cleaning as cleaning

    #Load necessary physics for overplotting planet velocity.
    vsys = sp.paramget('vsys',dp)
    RVp = sp.RV(dp)+vsys
    nexp = np.shape(ccf)[0]

    #Default to the entire y-axis if yrange = [0,0]
    if all(v == 0 for v in yrange):
        yrange=[0,nexp-1]

    x2,y2,z,rv_sel,y_sel,xticks,yticks,vmin,vmax = plotting_scales_2D(rv,fun.findgen(nexp),ccf,xrange,yrange,Nxticks=Nxticks,Nyticks=Nyticks,nsigma=3.0)
    #The plotting
    fig,ax = plt.subplots(figsize=(12,6))
    img=ax.pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
    ax.axis([x2.min(),x2.max(),y2.min(),y2.max()])
    line1, = ax.plot(RVp,fun.findgen(nexp),'--',color='black',label='Planet rest-frame')
    if len(doppler_model) > 0:
        line2, = ax.plot(doppler_model+vsys,fun.findgen(nexp),'--',color='black',label='Doppler shadow')
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_title(title)
    ax.set_xlabel('Radial velocity (km/s)')
    ax.set_ylabel('Exposure')
    #The colourbar
    cbar = plt.colorbar(img,format='%05.4f',aspect = 15)
    # cbar.set_norm(dcb.MyNormalize(vmin=vmin,vmax=vmax,stretch='linear'))
    cbar = dcb.DraggableColorbar_fits(cbar,img,'hot')
    cbar.connect()

    #The clickable legend.
    if len(doppler_model) > 0:
        lines = [line1, line2]
    else:
        lines = [line1]

    if i_legend == True:
        fancyplots.interactive_legend(fig,ax,lines)
    if show == True:
        plt.show()
    return(fig,ax,cbar)






class KpVsys_callback(object):
    def __init__(self,vorb):
        self.vorb=vorb


def plot_KpVsys(rv,Kp,KpVsys,dp,xrange=[-100,100],yrange=[0,0],Nxticks = 10.0,Nyticks = 10.0,title='',invert=False,injected=[0],filter=False):
    """This is a routine that does all the plotting of the KpVsys diagram.
    I expect this to be an organic function that is adapted to my plotting needs.
    THIS CURRENTLY HAS CORE USAGE IN RUN_INSTANCE, DO GOOD TESTS AND WRITE GOOD
    DOCUMENTATION.

    EDIT: THIS IS ALMOST THE SAME AS PLOT_CCF ABOVE. NEED TO MERGE THESE."""

    import numpy as np
    import matplotlib.pyplot as plt
    import pdb
    import lib.drag_colour as dcb
    import lib.functions as fun
    import lib.system_parameters as sp
    import lib.plotting as fancyplots
    import lib.operations as ops
    from matplotlib.widgets import Slider
    import math
    import sys

    if invert == True:
        KpVsys*=(-1.0)
    #Load necessary physics
    vsys = sp.paramget('vsys',dp)
    vorb = sp.v_orb(dp)
    vorb_sel = vorb*1.0
    #Define the extent of the images.
    nrv=len(rv)
    nKp=len(Kp)
    if all(v == 0 for v in yrange):
        yrange=[0,max(Kp)]

    x2,y2,z,rv_sel,y_sel,xticks,yticks,vmin,vmax = plotting_scales_2D(rv,fun.findgen(nKp),KpVsys,xrange,yrange,Nxticks=Nxticks,Nyticks=Nyticks,nsigma=3.0)

    #The plotting
    fig,ax = plt.subplots(1,2,sharex=True,figsize=(14,6))
    plt.subplots_adjust(left=0.05)#Make them more tight, we need all the space we can get.
    plt.subplots_adjust(right=0.9)
    img=ax[0].pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
    ax[0].axis([x2.min(),x2.max(),y2.min(),y2.max()])
    line1, = ax[0].plot(rv,rv*0.0+vorb,'--',color='black',label='Orbital velocity')
    line2, = ax[0].plot(Kp*0.0+vsys,Kp,'--',color='black',label='Systemic velocity')
    line3, = ax[0].plot(rv,rv*0.0+vorb_sel,'--',color='gray',label='Orbital velocity (selected)')

    ax[0].set_xticks(xticks)
    ax[0].set_yticks(yticks)
    ax[0].set_title(title)
    ax[0].set_xlabel('Radial velocity (km/s)')
    ax[0].set_ylabel('Orbital velocity (km/s)')

    dKp=Kp[1]-Kp[0]#Assumes that Kp is a linear axis.

    vorb_index = np.argmin(np.abs(Kp - vorb_sel))



    if filter == True:
        CCF_smooth = ops.smooth(KpVsys[vorb_index,:],50.0)
        CCF1D = KpVsys[vorb_index,:]-CCF_smooth
    else:
        CCF1D = KpVsys[vorb_index,:]

    img2=ax[1].plot(rv,CCF1D)
    # ax[1].plot(rv,CCF_smooth)
    sel=(np.abs(rv) > 50.0)

    threesigma = np.std(CCF1D[sel])*3.0+np.mean(CCF1D[sel])
    line4,=ax[1].plot([np.min(rv),np.max(rv)],threesigma*np.array([1.0,1.0]),color='green',alpha=0.5)

    if len(injected) != 1:
        CCF1D_i = injected[vorb_index,:]
        img2_i = ax[1].plot(rv,CCF1D_i-CCF1D,'--')
    ax[1].axvline(vsys,color='black')
    ax[1].set_xlim(x2.min(),x2.max())
    ax[1].set_ylim(z.min(),z.max())
    ax[1].set_xlabel('Radial velocity (km/s)')
    ax[1].set_ylabel('CCF')

    #The colourbar
    cbar = fig.colorbar(img,format='%05.4f',ax=[ax[0]],aspect=20)
    # cbar.set_norm(dcb.MyNormalize(vmin=vmin,vmax=vmax,stretch='linear'))
    cbar = dcb.DraggableColorbar_fits(cbar,[img],'hot')
    cbar.connect()

    #The clickable legend.
    lines = [line1, line2, line3]
    fancyplots.interactive_legend(fig,ax[0],lines)

    #Then the slider:
    rax_slider = plt.axes([0.65, 0.90, 0.15, 0.02])
    rax_slider.set_title('Orbital velocity (km/s)')
    vorb_slider = Slider(rax_slider,'', Kp.min(),Kp.max(),valinit=Kp[vorb_index],valstep=dKp)

    def update(val):
        vorb_sel = vorb_slider.val
        vorb_index = np.argmin(np.abs(Kp - vorb_sel))

        if filter == True:
            CCF_smooth = ops.smooth(KpVsys[vorb_index,:],50.0)
            CCF1D = KpVsys[vorb_index,:]-CCF_smooth
        else:
            CCF1D = KpVsys[vorb_index,:]

        img2[0].set_ydata(CCF1D)
        if len(injected) != 1:
            CCF1D_i = injected[vorb_index,:]
            img2_i[0].set_ydata(CCF1D_i-CCF1D)
        line3.set_ydata(vorb_sel)

        threesigma = np.std(CCF1D[sel])*3.0+np.mean(CCF1D[sel])
        line4.set_ydata(threesigma*np.array([1.0,1.0]))
        #ax[1].autoscale(enable=True,axis='y')
        #ylim = ax[1].get_ylim()
        #dylim = 10.0**np.max([math.floor(np.log10(ylim[0])),math.floor(np.log10(ylim[1]))])
        #if min(CCF1D) <= ylim[0] or max(CCF1D) >= ylim[1]:
    #        ax[1].set_ylim
        fig.canvas.draw_idle()

    vorb_slider.on_changed(update)



    plt.show()







def plot_KpVsys_old(rv,Kp,KpVsys,dp,xrange=[-100,100],Nticks = 10.0,title='',invert=False):
    """This is a routine that does all the plotting of the KpVsys diagram.
    I expect this to be an organic function that is adapted to my plotting needs.
    THIS CURRENTLY HAS CORE USAGE IN RUN_INSTANCE, DO GOOD TESTS AND WRITE GOOD
    DOCUMENTATION.

    EDIT: THIS IS ALMOST THE SAME AS PLOT_CCF ABOVE. NEED TO MERGE THESE."""

    import numpy as np
    import matplotlib.pyplot as plt
    import pdb
    import lib.drag_colour as dcb
    import lib.functions as fun
    import pylab as pl
    import lib.system_parameters as sp
    import lib.plotting as fancyplots


    if invert == True:
        KpVsys*=(-1.0)
    #Load necessary physics
    vsys = sp.paramget('vsys',dp)
    vorb = sp.v_orb(dp)
    #Define the extent of the images.
    nrv=len(rv)
    nKp=len(Kp)

    drv=rv[1]-rv[0]
    sel = ((rv >= xrange[0]) & (rv <=xrange[1]))
    KpVsys_sub = KpVsys[:,sel]
    m = np.nanmedian(KpVsys_sub)
    s = np.nanstd(KpVsys_sub)
    vmin = m-3.0*s
    vmax = m+3.0*s
    rvmin=rv[sel].min()
    rvmax=rv[sel].max()

    #This initiates the meshgrid.
    xmin = rvmin; xmax = rvmax; dx = drv
    ymin = 0; ymax = nKp-1; dy = 1
    x2,y2 = np.meshgrid(np.arange(xmin,xmax+dx+dx,dx)-dx/2.,np.arange(ymin,ymax+dy+dy,dy)-dy/2.)
    z = KpVsys_sub

    dxt = (xmax+dx - xmin) / Nticks
    dyt = (ymax+dy - ymin) / Nticks
    xticks = np.arange(xmin,xmax+dx,round(dxt,0))
    yticks = np.arange(ymin,ymax+dy,round(dyt,0))



    #The plotting
    fig,ax = plt.subplots(1,2,sharex=True,figsize=(14,6))
    plt.subplots_adjust(left=0.05)#Make them more tight, we need all the space we can get.
    plt.subplots_adjust(right=0.9)
    img=ax[0].pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
    ax[0].axis([x2.min(),x2.max(),y2.min(),y2.max()])
    line1, = ax[0].plot(rv,rv*0.0+vorb,'--',color='black',label='Orbital velocity')
    line2, = ax[0].plot(Kp*0.0+vsys,Kp,'--',color='black',label='Systemic velocity')
    pl.xticks(xticks)
    pl.yticks(yticks)
    ax[0].set_title(title)
    ax[0].set_xlabel('Radial velocity (km/s)')
    ax[0].set_ylabel('Exposure')

    #The colourbar
    cbar = fig.colorbar(img,format='%05.4f',ax=[ax[0]],aspect=20)
    # cbar.set_norm(dcb.MyNormalize(vmin=vmin,vmax=vmax,stretch='linear'))
    cbar = dcb.DraggableColorbar_fits(cbar,[img],'hot')
    cbar.connect()

    #The clickable legend.
    lines = [line1, line2]

    fancyplots.interactive_legend(fig,ax,lines)
    plt.show()



def combine_KpVsys(list_of_datasets,list_of_templates):
    """This small function takes the KpVsys output of multiple templates and
    datasets (likely multiple nights on the same object, because different objects
    do not have the same Kp) and averages them."""
    import lib.utils as ut
    import astropy.io.fits as fits
    ut.typetest('list_of_templates in combine_KpVsys',list_of_templates,list)
    ut.typetest('list_of_datasets in combine_KpVsys',list_of_datasets,list)

    Nd = len(list_of_datasets)
    Nt = len(list_of_templates)
    list_of_KpVsys = []
    for i in range(Nd):
        for j in range(Nt):
            inpath = 'output/'+list_of_datasets[i]+'/'+list_of_templates[j]+'/KpVsys.fits'
            list_of_KpVsys.append(fits.getdata(inpath))


    Nk=len(list_of_KpVsys)
    KpVsys_combined = list_of_KpVsys[0]*0.0
    for i in range(Nk):
        KpVsys_combined+=list_of_KpVsys[i]

    return(KpVsys_combined/Nk)
