def plot_models(list_of_models,mode):
    import lib.models as models
    import sys
    import numpy as np
    import lib.constants as const
    import lib.system_parameters as sp
    import matplotlib.pyplot as plt
    import lib.operations as ops
    if mode < 0 or mode > 2:
        print('ERROR: Mode should be set to:')
        print(' 0, just plot the high res spectra.')
        print(' 1, just plot the broadened spectra as injected.')
        print(' 2, overplot the hires and broadened spectra.')
        sys.exit()

    minwl = 350.0
    maxwl = 700.0
    modellib = 'models/library_WASP_121'
    dp = 'data/Wasp-121/night1/'
    dRV=np.median(sp.dRV(dp))
    Rd=sp.paramget('resolution',dp)
    planet_radius=sp.paramget('Rp',dp)
    inclination=sp.paramget('inclination',dp)
    P=sp.paramget('P',dp)


    fig, axes = plt.subplots(1,1,figsize=(12,7))
    plt.subplots_adjust(left=0.1)#Make them more tight, we need all the space we can get.
    plt.subplots_adjust(right=0.98)
    plt.subplots_adjust(top=0.98)
    plt.subplots_adjust(bottom=0.07)
    plt.subplots_adjust(hspace=0.04)

    for i in list_of_models:
        a = models.get_model(i,library = modellib)
        wlm = a[0]
        fxm = a[1]
        if wlm[-1] <= wlm[0]:#Reverse the wl axis if its sorted the wrong way.
            wlm=np.flipud(wlm)
            fxm=np.flipud(fxm)
        #Select right part of model
        modelsel=[(wlm >= minwl) & (wlm <= maxwl)]
        wlm=wlm[tuple(modelsel)]
        fxm=fxm[tuple(modelsel)]

        if mode < 2:
            void = plt.plot(wlm,fxm,alpha=0.4)
        if mode > 0:
            #Blur th model
            print('Rotation blurring')
            fxm_b=ops.blur_rotate(wlm,fxm,(const.c/1000.0)/Rd,planet_radius,P,inclination,status=True)
            wlm_cv,fxm_bcv,vstep=ops.constant_velocity_wl_grid(wlm,fxm_b,oversampling=1.5)
            print('Exptime blurring')
            fxm_b2=ops.smooth(fxm_bcv,dRV/vstep,mode='box')
            if mode == 1:
                plt.plot(wlm_cv,fxm_b2,alpha=0.6,color=void[0].get_color())
            if mode == 2:
                plt.plot(wlm_cv,fxm_b2,alpha=0.6)
    plt.axvline(373.0)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Transit depth')
    plt.show()






def plot_species(species,skipnight=-1,skipplot=False):
    """This plots the three individual nights plus combination of a single species."""
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
    import lib.analysis as ana
    import copy

    N = 3#Number of nights
    Nxticks=10.0
    Nyticks=10.0

    fig, axes = plt.subplots(2,N+1,figsize=(12,7),sharex=True)
    plt.subplots_adjust(left=0.1)#Make them more tight, we need all the space we can get.
    plt.subplots_adjust(right=0.98)
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(bottom=0.07)
    plt.subplots_adjust(hspace=0.04)
    plots = {}#initiate dictionary.
    list_of_ccfs = []
    list_of_ccfis = []
    list_of_kpvsys = []
    list_of_ccfes = []
    vmin_all=0.0
    vmax_all=0.0
    for i in range(1,N+2):
        if i < N+1 and i != skipnight:
            inpath = ('output/Wasp-121/night%s/' % i)+species+'/'
            dp = 'data/Wasp-121/night%s/' % i

            RV = fits.getdata(inpath+'RV.fits')
            k = fits.getdata(inpath+'KpVsys.fits')*(-1.0)
            ke = fits.getdata(inpath+'KpVsys_e.fits')
            ki = fits.getdata(inpath+'KpVsys_i.fits')*(-1.0)
            Kp = fits.getdata(inpath+'Kp.fits')


            xrange=[-150,150]
            yrange=[min(Kp),max(Kp)]
            vsys = sp.paramget('vsys',dp)
            vorb = sp.v_orb(dp)
            vorb_sel = vorb*1.0
            vorb_index = np.argmin(np.abs(Kp - vorb_sel))
            ccf1d = k[vorb_index,:]
            ccf1de = ke[vorb_index,:]
            ccf1di = ki[vorb_index,:]-k[vorb_index,:]
            list_of_ccfs.append(ccf1d)
            list_of_ccfis.append(ccf1di)
            list_of_kpvsys.append(k)
            list_of_ccfes.append(ccf1de)

            x2,y2,z,rv_sel,y_sel,xticks,yticks,vmin,vmax = ana.plotting_scales_2D(RV,Kp,k,xrange,yrange,Nxticks=Nxticks,Nyticks=Nyticks,nsigma=3.0)
            vmin_all = min([vmin_all,vmin])
            vmax_all = max([vmax_all,vmax])
            plots['Kp'] = Kp
            plots['RV'] = RV
            plots['vsys'] = vsys
            plots['vorb'] = vorb
            plots['plot_KpVsys_night%s'%i]=axes[0][i-1].pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
            plots['shade_CCF_night%s'%i]=axes[1][i-1].fill_between(RV, ccf1d-ccf1de, ccf1d+ccf1de,alpha=0.5)
            plots['plot_CCF_night%s'%i]=axes[1][i-1].plot(RV,ccf1d)
            plots['plot_CCF_i_night%s'%i]=axes[1][i-1].plot(RV,ccf1di,'--',alpha=0.5)
            axes[1][i-1].axvline(vsys,linestyle='--',color='black',alpha=0.5)
            axes[0][i-1].axvline(vsys,linestyle='--',color='black',alpha=0.5)
            axes[0][i-1].axhline(vorb,linestyle='--',color='black',alpha=0.5)
            axes[0][i-1].set_title('Night%s'%i)
            axes[1][i-1].set_xlabel('Sys. radial velocity (km/s)')
            plots['KpVsys_night%s'%i] = k
            plots['CCF_night%s'%i] = ccf1d



        if i == N+1:#The last column is the combined average of the other plots.
            K = 0.0#Initialize as zeroes.
            CCF1D = 0.0
            CCF1DI = 0.0
            CCF1DE2 = 0.0#Error bars.
            for j in range(len(list_of_ccfs)):#Compute averages.
                K+= list_of_kpvsys[j]/N
                CCF1D+=list_of_ccfs[j]/N
                CCF1DE2+=list_of_ccfes[j]**2.0
                CCF1DI+=list_of_ccfis[j]/N
            CCF1DE=np.sqrt(CCF1DE2)/N
            x2,y2,z,rv_sel,y_sel,xticks,yticks,vmin,vmax = ana.plotting_scales_2D(RV,Kp,K,xrange,yrange,Nxticks=Nxticks,Nyticks=Nyticks,nsigma=4.0)
            plots['KpVsys_avg'] = z
            plots['CCF_avg'] = CCF1D
            plots['CCF_i_avg'] = CCF1DI
            plots['plot_KpVsys_avg']=axes[0][i-1].pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
            plots['shade_CCF_avg']=axes[1][i-1].fill_between(RV, CCF1D-CCF1DE, CCF1D+CCF1DE,alpha=0.5)
            plots['plot_CCF_avg']=axes[1][i-1].plot(RV,CCF1D)
            plots['plot_CCF_i_avg']=axes[1][i-1].plot(RV,CCF1DI,'--',alpha=0.5)
            axes[1][i-1].axvline(vsys,linestyle='--',color='black',alpha=0.5)
            axes[0][i-1].axvline(vsys,linestyle='--',color='black',alpha=0.5)
            axes[0][i-1].axhline(vorb,linestyle='--',color='black',alpha=0.5)
            axes[0][i-1].set_title('Combined')
            axes[1][i-1].set_xlabel('Sys. radial velocity (km/s)')
            plots['KpVsys_avg'] = K
            plots['CCF_avg'] = CCF1D


    ccf_sample = list_of_ccfs+list_of_ccfis #LIST CONCATANATION

    ymin = min(np.concatenate(tuple(ccf_sample)))#Then they become all stitched together.
    ymax = max(np.concatenate(tuple(ccf_sample)))

    for i in range(0,N+1):
        axes[1][i].set_ylim(ymin,ymax)
        if i >= 1:
            axes[1][i].set_yticklabels([])
    axes[1][0].set_xlim(xrange[0],xrange[1])
    axes[0][0].set_ylabel('Orbital velocity (km/s)')


    for i in range(1,N):
        if i != skipnight:
            plots['plot_KpVsys_night%s'%i].set_clim(vmin=vmin_all,vmax=vmax_all)
    plots['plot_KpVsys_avg'].set_clim(vmin=vmin_all,vmax=vmax_all)

    axes[1][0].set_ylabel('Excess transit depth')


    cb_ax = fig.add_axes([0.1, 0.96, 0.86, 0.03])
    cbar = fig.colorbar(plots['plot_KpVsys_avg'], cax=cb_ax,orientation='horizontal')
    cbar.set_label('Excess transit depth in '+species,fontsize=14)

    # cbar = fig.colorbar(im, cax=cb_ax,orientation='horizontal')


    if skipplot == False:
        plt.show()
    else:
        plt.close()
    return(plots)



def slide_plots_w121():
    import matplotlib.pyplot as plt
    import pdb
    import lib.analysis as ana
    species = ['FeI','CrI','TiI','VI']
    N=len(species)
    ccfs = []
    xrange=[-150,150]
    Nxticks=10.0
    Nyticks=10.0
    fig, axes = plt.subplots(2,N,figsize=(12,7),sharex=True)
    plt.subplots_adjust(left=0.1)#Make them more tight, we need all the space we can get.
    plt.subplots_adjust(right=0.98)
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(bottom=0.07)
    plt.subplots_adjust(hspace=0.04)

    for i in range(N):
        print('Doing '+species[i])
        ccfs.append(plot_species(species[i],skipplot=True))

    for i in range(N):
        dic = ccfs[i]
        Kp = dic['Kp']
        RV = dic['RV']
        vsys = dic['vsys']
        vorb = dic['vorb']
        yrange=[min(Kp),max(Kp)]
        CCF1D = dic['CCF_avg']
        KpVsys = dic['KpVsys_avg']

        x2,y2,z,rv_sel,y_sel,xticks,yticks,vmin,vmax = ana.plotting_scales_2D(RV,Kp,KpVsys,xrange,yrange,Nxticks=Nxticks,Nyticks=Nyticks,nsigma=4.0)
        axes[0][i].pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
        # axes[1][i-1].fill_between(RV, CCF1D-CCF1DE, CCF1D+CCF1DE,alpha=0.5)
        axes[1][i].plot(RV,CCF1D)
        # plots['plot_CCF_i_avg']=axes[1][i-1].plot(RV,CCF1DI,'--',alpha=0.5)
        axes[1][i].axvline(vsys,linestyle='--',color='black',alpha=0.5)
        axes[0][i].axvline(vsys,linestyle='--',color='black',alpha=0.5)
        axes[0][i].axhline(vorb,linestyle='--',color='black',alpha=0.5)
        axes[0][i].set_title(species[i])
        axes[1][i].set_xlabel('Sys. radial velocity (km/s)')
        axes[0][i].set_xlim(xrange[0],xrange[1])
        axes[1][i].set_xlim(xrange[0],xrange[1])
        axes[1][i].set_ylim(-0.0003,0.0007)
        if i > 0 :
            axes[0][i].set_yticklabels([])
            axes[1][i].set_yticklabels([])
        if i == 0:
            axes[0][i].set_ylabel('Orbital velocity (km/s)')
            axes[1][i].set_ylabel('Transit depth')

    plt.show()




def slide_plots_w189():
    import matplotlib.pyplot as plt
    import pdb
    import lib.analysis as ana
    import astropy.io.fits as fits
    import lib.system_parameters as sp
    import numpy as np
    species = ['FeI','FeII','TiI','CrI']
    N=len(species)
    ccfs = []
    xrange=[-150,150]
    Nxticks=10.0
    Nyticks=10.0
    fig, axes = plt.subplots(2,N,figsize=(12,7),sharex=True)
    plt.subplots_adjust(left=0.1)#Make them more tight, we need all the space we can get.
    plt.subplots_adjust(right=0.98)
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(bottom=0.07)
    plt.subplots_adjust(hspace=0.04)


    path = 'output/Wasp-189/night1/'
    dp = 'data/Wasp-189/night1/'
    vsys = sp.paramget('vsys',dp)
    vorb = sp.v_orb(dp)
    vorb_sel = vorb*1.0

    for i in range(N):
        inpath=path+species[i]
        Kp = fits.getdata(inpath+'/Kp.fits')
        RV = fits.getdata(inpath+'/RV.fits')
        yrange=[min(Kp),max(Kp)]
        KpVsys = fits.getdata(inpath+'/KpVsys.fits')*(-1.0)

        vorb_index = np.argmin(np.abs(Kp - vorb_sel))

        CCF1D = KpVsys[vorb_index,:]
        if species[i] == 'FeII':
            CCF1D*=0.1


        x2,y2,z,rv_sel,y_sel,xticks,yticks,vmin,vmax = ana.plotting_scales_2D(RV,Kp,KpVsys,xrange,yrange,Nxticks=Nxticks,Nyticks=Nyticks,nsigma=4.0)
        axes[0][i].pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
        # axes[1][i-1].fill_between(RV, CCF1D-CCF1DE, CCF1D+CCF1DE,alpha=0.5)
        axes[1][i].plot(RV,CCF1D)
        # plots['plot_CCF_i_avg']=axes[1][i-1].plot(RV,CCF1DI,'--',alpha=0.5)
        axes[1][i].axvline(vsys,linestyle='--',color='black',alpha=0.5)
        axes[0][i].axvline(vsys,linestyle='--',color='black',alpha=0.5)
        axes[0][i].axhline(vorb,linestyle='--',color='black',alpha=0.5)
        axes[0][i].set_title(species[i])
        axes[1][i].set_xlabel('Sys. radial velocity (km/s)')
        axes[0][i].set_xlim(xrange[0],xrange[1])
        axes[1][i].set_xlim(xrange[0],xrange[1])
        axes[1][i].set_ylim(-0.00005,0.0001)
        if i > 0 :
            axes[0][i].set_yticklabels([])
            axes[1][i].set_yticklabels([])
        if i == 0:
            axes[0][i].set_ylabel('Orbital velocity (km/s)')
            axes[1][i].set_ylabel('Transit depth')

    plt.show()


def slide_plots_m2():
    import matplotlib.pyplot as plt
    import pdb
    import lib.analysis as ana
    import astropy.io.fits as fits
    import lib.system_parameters as sp
    import numpy as np
    species = ['FeI','FeII','CrII','TiII']
    N=len(species)
    ccfs = []
    xrange=[-150,150]
    Nxticks=10.0
    Nyticks=10.0
    fig, axes = plt.subplots(2,N,figsize=(12,7),sharex=True)
    plt.subplots_adjust(left=0.1)#Make them more tight, we need all the space we can get.
    plt.subplots_adjust(right=0.98)
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(bottom=0.07)
    plt.subplots_adjust(hspace=0.04)


    path = 'output/MASCARA-2/'
    dp = 'data/MASCARA-2/'
    vsys = sp.paramget('vsys',dp)
    vorb = sp.v_orb(dp)
    vorb_sel = vorb*1.0

    for i in range(N):
        inpath=path+species[i]
        Kp = fits.getdata(inpath+'/Kp.fits')
        RV = fits.getdata(inpath+'/RV.fits')
        yrange=[min(Kp),max(Kp)]
        KpVsys = fits.getdata(inpath+'/KpVsys.fits')*(-1.0)

        vorb_index = np.argmin(np.abs(Kp - vorb_sel))

        CCF1D = KpVsys[vorb_index,:]
        if species[i] == 'FeII' or species[i] == 'CrII' or species[i]=='TiII':
            CCF1D*=0.1

        x2,y2,z,rv_sel,y_sel,xticks,yticks,vmin,vmax = ana.plotting_scales_2D(RV,Kp,KpVsys,xrange,yrange,Nxticks=Nxticks,Nyticks=Nyticks,nsigma=3.0)
        axes[0][i].pcolormesh(x2,y2,z,vmin=vmin,vmax=vmax,cmap='hot')
        # axes[1][i-1].fill_between(RV, CCF1D-CCF1DE, CCF1D+CCF1DE,alpha=0.5)
        axes[1][i].plot(RV,CCF1D)
        # plots['plot_CCF_i_avg']=axes[1][i-1].plot(RV,CCF1DI,'--',alpha=0.5)
        axes[1][i].axvline(vsys,linestyle='--',color='black',alpha=0.5)
        axes[0][i].axvline(vsys,linestyle='--',color='black',alpha=0.5)
        axes[0][i].axhline(vorb,linestyle='--',color='black',alpha=0.5)
        axes[0][i].set_title(species[i])
        axes[1][i].set_xlabel('Sys. radial velocity (km/s)')
        axes[0][i].set_xlim(xrange[0],xrange[1])
        axes[1][i].set_xlim(xrange[0],xrange[1])
        axes[1][i].set_ylim(-0.0002,0.00028)
        if i > 0 :
            axes[0][i].set_yticklabels([])
            axes[1][i].set_yticklabels([])
        if i == 0:
            axes[0][i].set_ylabel('Orbital velocity (km/s)')
            axes[1][i].set_ylabel('Transit depth')

    plt.show()
