


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

    transit = sp.transit(dp)
    sel = (transit == 1.0)
    n_exp_out_of_transit = len(np.where(transit == 1)[0])

    CCF1D = np.mean(CCF[sel,:],axis=0)
    if len(np.shape(CCF_E)) > 0:
        a2.errorbar(RV,CCF1D,fmt='.',yerr=np.mean(np.sqrt(CCF_E),axis=0)/np.sqrt(n_exp_out_of_transit),zorder=3)
    else:
        a2.plot(RV,CCF1D,'.')

    maxindex = (abs(CCF1D-np.nanmedian(CCF1D))).argmax()
    RV_0 = RV[maxindex]
    sel = [(RV > RV_0-dv) & (RV < RV_0+dv)]
    y = CCF1D
    fit=ops.gauss_fit(RV[sel],y[sel],start=[0.0,RV_0,dv,np.nanmedian(y)])
    a2.plot(RV,y,'.')
    a2.plot(RV[sel],y[sel],'.',color='red')
    a2.plot(RV[sel],fun.gaussian(RV[sel],fit[0],fit[1],fit[2],cont=fit[3]))
    a2.text(fit[1]+dv/2,fit[3]+fit[0],'RV = %s km/s' % round(fit[1],3),color='green')
    a0.set_title('t-averaged data and template')
    a1.set_title('Sum(T)')
    a2.set_title('t-averaged CCF')
    a0.tick_params(axis='both',labelsize='5')
    a02.tick_params(axis='both',labelsize='5')
    a1.tick_params(axis='both',labelsize='5')
    a2.tick_params(axis='both',labelsize='5')
    fig.tight_layout()
    plt.show()


#Kelt-9 radius ratio:
#Rp/Rs = df = 0.08228 \pm 0.00043
#Rs = 2.196 \pm 0.012 Rssun
#Rp = df*Rs

#Propagate errors:
#S_Rp^2 = (dRp/ddf sdf)^2 + (dRp/dRs sRs)^2
def plot_RV_star(dp,RV,CCF2D,RVrange=[]):
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


def construct_KpVsys(rv,ccf,dp,kprange=[0,300],dkp=1.0):
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
    transit = sp.transit(dp)-1.0
    transit /= np.nansum(transit)
    transitblock = fun.rebinreform(transit,len(rv)).T

    j = 0
    ccfs = []
    for i in Kp:
        dRV = sp.RV(dp,vorb=i)*(-1.0)
        ccf_shifted = ops.shift_ccf(rv,ccf,dRV)
        ccfs.append(ccf_shifted)
        KpVsys[j,:] = np.nansum(transitblock * ccf_shifted,axis=0)
        j+=1
    ut.save_stack('test.fits',ccfs)
    return(Kp,KpVsys)


def plot_ccf(rv,ccf,dp,xrange,Nticks = 10.0,title='',doppler_model = []):
    """This is a routine that does all the plotting of the cleaned 2D CCF.
    I expect this to be an organic function that is adapted to my plotting needs."""

    import numpy as np
    import matplotlib.pyplot as plt
    import pdb
    import lib.drag_colour as dcb
    import lib.functions as fun
    import pylab as pl
    import lib.system_parameters as sp
    import lib.plotting as fancyplots

    #Load necessary physics
    vsys = sp.paramget('vsys',dp)
    RVp = sp.RV(dp)+vsys

    #Define the extent of the images.
    nrv=len(rv)
    nexp=np.shape(ccf)[0]
    y = fun.findgen(nexp)
    drv=rv[1]-rv[0]
    sel = ((rv >= xrange[0]) & (rv <=xrange[1]))
    ccf_sub = ccf[:,sel]
    m = np.nanmedian(ccf_sub)
    s = np.nanstd(ccf_sub)
    vmin = m-3.0*s
    vmax = m+3.0*s
    rvmin=rv[sel].min()
    rvmax=rv[sel].max()

    #This initiates the meshgrid.
    xmin = rvmin; xmax = rvmax; dx = drv
    ymin = 0; ymax = nexp-1; dy = 1
    x2,y2 = np.meshgrid(np.arange(xmin,xmax+dx+dx,dx)-dx/2.,np.arange(ymin,ymax+dy+dy,dy)-dy/2.)
    z = ccf_sub

    dxt = (xmax+dx - xmin) / Nticks
    dyt = (ymax+dy - ymin) / Nticks
    xticks = np.arange(xmin,xmax+dx,round(dxt,0))
    yticks = np.arange(ymin,ymax+dy,round(dyt,0))

    #The plotting
    fig,ax = plt.subplots(figsize=(12,6))
    img=ax.pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
    ax.axis([x2.min(),x2.max(),y2.min(),y2.max()])
    line1, = ax.plot(RVp,fun.findgen(nexp),'--',color='black',label='Planet rest-frame')
    if len(doppler_model) > 0:
        line2, = ax.plot(doppler_model+vsys,fun.findgen(nexp),'--',color='black',label='Doppler shadow')
    pl.xticks(xticks)
    pl.yticks(yticks)
    ax.set_title(title)
    plt.xlabel('Radial velocity (km/s)')
    plt.ylabel('Exposure')

    #The colourbar
    cbar = plt.colorbar(img,format='%05.4f',aspect = 15)
    cbar.set_norm(dcb.MyNormalize(vmin=vmin,vmax=vmax,stretch='linear'))
    cbar = dcb.DraggableColorbar(cbar,img)
    cbar.connect()

    #The clickable legend.
    if len(doppler_model) > 0:
        lines = [line1, line2]
    else:
        lines = [line1]
    fancyplots.interactive_legend(fig,ax,lines)
    plt.show()

def plot_KpVsys(rv,Kp,KpVsys,dp,xrange=[-100,100],Nticks = 10.0,title=''):
    """This is a routine that does all the plotting of the KpVsys diagram.
    I expect this to be an organic function that is adapted to my plotting needs."""

    import numpy as np
    import matplotlib.pyplot as plt
    import pdb
    import lib.drag_colour as dcb
    import lib.functions as fun
    import pylab as pl
    import lib.system_parameters as sp
    import lib.plotting as fancyplots

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
    fig,ax = plt.subplots()
    img=ax.pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
    ax.axis([x2.min(),x2.max(),y2.min(),y2.max()])
    line1, = ax.plot(rv,rv*0.0+vorb,'--',color='black',label='Orbital velocity')
    line2, = ax.plot(Kp*0.0+vsys,Kp,'--',color='black',label='Systemic velocity')
    pl.xticks(xticks)
    pl.yticks(yticks)
    ax.set_title(title)
    plt.xlabel('Radial velocity (km/s)')
    plt.ylabel('Exposure')

    #The colourbar
    cbar = plt.colorbar(img,format='%05.4f',aspect = 15)
    cbar.set_norm(dcb.MyNormalize(vmin=vmin,vmax=vmax,stretch='linear'))
    cbar = dcb.DraggableColorbar(cbar,img)
    cbar.connect()

    #The clickable legend.
    lines = [line1, line2]

    fancyplots.interactive_legend(fig,ax,lines)
    plt.show()
