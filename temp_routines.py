def read_noise_csv(filename,wl):
    import csv
    import pdb
    from scipy import interpolate
    import numpy as np
    import matplotlib.pyplot as plt
    with open(filename, mode='r') as csv_file:
        reader = csv.DictReader(csv_file)
        data = {}
        for row in reader:
            for header, value in row.items():
                try:
                    data[header].append(np.float(value))
                except KeyError:
                    data[header] = [np.float(value)]


# extract the variables you want
        wls = np.array(data['wavelength'])/10.0
        counts = np.array(data['total_counts'])

        i_s = interpolate.interp1d(wls,np.sqrt(counts),fill_value='extrapolate')
        SNR_i=i_s(wl)
        return(SNR_i)


def blur_HST(wlt,T,dRV):
    import lib.operations as ops
    #dRV = #For HST Medium res grating (10) plus intrinsic width (20).
    wlt_cv,T_cv,vstep=ops.constant_velocity_wl_grid(wlt,T,oversampling=1.0)
    # print('v_step is %s km/s' % vstep)
    # print('So the resolution blurkernel has an avg width of %s px.' % (dRV/vstep))
    T_b=ops.smooth(T_cv,dRV/vstep,mode='gaussian')
    return(wlt_cv,T_b)

def HST_proposal():
    import numpy as np
    import lib.models as models
    import lib.operations as ops
    import matplotlib.pyplot as plt
    import lib.functions as fun
    import pdb
    import astropy.io.fits as fits
    import os
    import lib.utils as ut
    from scipy import interpolate
    import lib.analysis as an
    mp='../Kelt-9/TEMP_new_grid_apr_1/templates/var_grav_4000K/'
    mpl=mp+'../../all_species/'
    outp='../../Proposals/HST_Cycle_27/plot_ccfs/'
    Rsun = 695510.0
    Rs = 2.196*Rsun
    atoms = ['0','H','He','Na','Sc','Cr','Co']

    atoms=['Ag','Cr','Dy','Gd','K','Lu','Na','O','Pt','Ru','Sc','Sm','V','Yb','As','Be','Ca','Cl','Ho','Mo','Pb','Rb','Tb','Al','Au','Er','Fe','Ge','Mg','Re','Sn','Tl','Zn','B','Cs','In','La','Nb','Os','Pd','S','Se','Te','Tm','W','Bi','Cd','Co','Cu','Eu','Hf','N','Ni','Rh','Si','Sr','Th','Ar','C','Ga','Ir','Li','Mn','P','Pr','U','Y','Zr','Ba','Ce','F','Hg','Nd','Sb','Ta','Ti']
    atoms=['Fe','Mn','Ca','Co','Cr','Mg','Ni','Si','Ti','V']

    ions = ['','_p','_2p']
    spec_low3 = np.flipud(np.squeeze(fun.read_binary(mpl+'transit_m_1_t_3000.dat',double=True)))
    spec_low4 = np.flipud(np.squeeze(fun.read_binary(mpl+'transit_m_1_t_4000.dat',double=True)))
    spec_low5 = np.flipud(np.squeeze(fun.read_binary(mpl+'transit_m_1_t_5000.dat',double=True)))
#The following reads the spectral_grid.dat file. I have run it before and
#written the output to a fitsfile, which reads faster.
    # import csv
    # wl=[]
    # with open(mp+'spectral_grid.dat',mode='r') as csv_file:
    #     csv_reader = csv.DictReader(csv_file,delimiter='\t')
    #     for line in csv_reader:
    #         wl.append(np.float(line['#wavelength (microns)']))
    #     wl=np.array(wl)*1000.0
    # fits.writeto(mp+'spectral_grid.fits',wl,overwrite=True)
    wl=np.flipud(fits.getdata(mp+'spectral_grid.fits'))



    # t1=ut.start()
    #Fe = np.flipud(np.squeeze(fun.read_binary(mp+'Fe_p/transit_m_1_t_4000.dat',double=True)))
    Fe = np.flipud(np.squeeze(fun.read_binary(mpl+'transit_m_1_t_4000.dat',double=True)))
    cont = np.flipud(np.squeeze(fun.read_binary(mp+'N_2p/transit_m_1_t_4000.dat',double=True)))
    dRV=(10.0**2.0 + 20.0**2.0)**0.5
    wl_b,Fe_b = blur_HST(wl,Fe,dRV)
    wlt,T = blur_HST(wl,(Fe-cont)*(-1.0),10.0)
    wlmin = 230.3
    wlmax = 311.1
        #Now I need to model the spectrograph and the noise.
    SNR = 80.7*np.sqrt(50.0*60.0/600.0*2.0*2.0)/np.sqrt(2.0) #per 1 pixel per 2 orbits per 2 transits.
    x=fun.findgen(90000.0)
    wld = np.exp(5.0/3e5 * x)*220.0#5km/s per pixel of STIS. Min wl is 230nm.
    noise = np.random.normal(0.0,1.0/SNR,len(wld))
    #Ok now interpolate and add noise....
    i_s = interpolate.interp1d(wl_b,1.0-(Fe_b/Rs)**2.0)
    Fe_STIS =i_s(wld)


    wl_b_low,spec_b_low3 = blur_HST(wl,spec_low3,600.0)#600kms is a resolution of 500.0
    wl_b_low,spec_b_low4 = blur_HST(wl,spec_low4,600.0)#600kms is a resolution of 500.0
    wl_b_low,spec_b_low5 = blur_HST(wl,spec_low5,600.0)#600kms is a resolution of 500.0
    SNR_low = 2000.0*np.sqrt(50.0*60.0/300.0*2.0/2.0)#SNR is about 2000 after a 300s exposure.
    #multiply by sqrt of number of seconds in a 50m orbit, times 2 orbits in transit, divide by sqrt 2
    #because the noise of the out-of-transit baseline is added. ADJUSTED BELOW! into an array
    #to take into the wavelength dependent output of the ETC.
    x_low=fun.findgen(1280.0)
    wld_low = np.exp(300.0/3e5 * x_low)*250.0#5km/s per pixel of STIS. Min wl is 230nm.

    i_slo3 = interpolate.interp1d(wl_b_low,1.0-(spec_b_low3/Rs)**2.0)
    i_slo4 = interpolate.interp1d(wl_b_low,1.0-(spec_b_low4/Rs)**2.0)
    i_slo5 = interpolate.interp1d(wl_b_low,1.0-(spec_b_low5/Rs)**2.0)

    spec_low_STIS3 =i_slo3(wld_low)
    spec_low_STIS4 =i_slo4(wld_low)
    spec_low_STIS5 =i_slo5(wld_low)


    wld_low_binned = []
    spec_low_binned3 = []
    spec_low_binned4 = []
    spec_low_binned5 = []
    binfactor = 10
    k = binfactor*1
    while k < np.max(wld_low):
        wld_low_bin=np.mean(wld_low[(k-binfactor):k])
        spec_low_bin3=np.mean(spec_low_STIS3[(k-binfactor):k])
        spec_low_bin4=np.mean(spec_low_STIS4[(k-binfactor):k])
        spec_low_bin5=np.mean(spec_low_STIS5[(k-binfactor):k])
        wld_low_binned.append(wld_low_bin)
        spec_low_binned3.append(spec_low_bin3)
        spec_low_binned4.append(spec_low_bin4)
        spec_low_binned5.append(spec_low_bin5)
        k+=binfactor

    wllow_min = 290.0
    wllow_max = 570.0
    wld_low_binned=np.array(wld_low_binned)
    spec_low_binned3=np.array(spec_low_binned3)
    spec_low_binned4=np.array(spec_low_binned4)
    spec_low_binned5=np.array(spec_low_binned5)

    wld_low_binned_sel = wld_low_binned[(wld_low_binned >= wllow_min) & (wld_low_binned <= wllow_max)]
    spec_low_binned_sel3 = spec_low_binned3[(wld_low_binned >= wllow_min) & (wld_low_binned <= wllow_max)]
    spec_low_binned_sel4 = spec_low_binned4[(wld_low_binned >= wllow_min) & (wld_low_binned <= wllow_max)]
    spec_low_binned_sel5 = spec_low_binned5[(wld_low_binned >= wllow_min) & (wld_low_binned <= wllow_max)]
    SNR_low=read_noise_csv(outp+'../G430L_ETC_300s.csv',wld_low_binned_sel)*np.sqrt(50.0*60.0/300.0*2.0/2.0)

    fig = plt.figure(figsize=(10,3))
    ax = fig.add_subplot(1,1,1)
    ax.plot(wld,Fe_STIS,color='lightgray',label='High-res spectrum at 4,000 K.',linewidth=0.7)
    #plt.plot(wld_low,spec_low_STIS,'.')
    ax.errorbar(wld_low_binned_sel,spec_low_binned_sel3,fmt='.',yerr=3.0/SNR_low/np.sqrt(binfactor),zorder=3,label='T = 3,000K')
    ax.errorbar(wld_low_binned_sel,spec_low_binned_sel4,fmt='.',yerr=3.0/SNR_low/np.sqrt(binfactor),zorder=3,label='T = 4,000K')
    ax.errorbar(wld_low_binned_sel,spec_low_binned_sel5,fmt='.',yerr=3.0/SNR_low/np.sqrt(binfactor),zorder=3,label='T = 5,000K')
    #plt.plot(wld_low,spec_low_STIS,'.')
    ax.set_xlim((270.0,590.0))
    ax.set_ylim((0.9892,0.9925))
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Transit depth (relative to star)')
    ax.tick_params(axis='both',labelsize='7')
    ax.legend(loc='best')

    fig.savefig('../../Proposals/HST_Cycle_27/lowres.pdf',format='pdf',bbox_inches='tight',dpi=500)
    plt.close(fig)

    pdb.set_trace()



    #Here we do the cross-correlation once for the purpose of testing/playing.
    #Below, its done in a forloop to generate the multiplot.
    nmodels=60
    # modelstack=[Fe_STIS]
    # for i in range(0,nmodels):
    #     modelstack.append(Fe_STIS+np.random.normal(0.0,1.0/SNR,len(wld)))
    # rv,ccf = an.xcor([wld],[np.stack(modelstack)],wlt,T,5.0,600.0)
    # plt.plot(rv,ccf[0,:],color='black',linewidth=1.5)
    # ccf_avg = np.mean(ccf,axis=0)
    # ccf_std = np.std(ccf,axis=0)
    # plt.errorbar(rv,ccf_avg,fmt='.',yerr=ccf_std,color='gray')
    # for i in range(1,nmodels):
    #     plt.plot(rv,ccf[i,:],'.',alpha=2.0/nmodels/2.0,color='blue')
    # plt.show()





# f, (a0, a1) = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1]})
# a0.plot(x,y)
# a1.plot(y,x)

# f.tight_layout()
# f.savefig('grid_figure.pdf')




    for i in atoms:
        Ipath=mp+i+ions[0]+'/transit_m_1_t_4000.dat'
        IIpath=mp+i+ions[1]+'/transit_m_1_t_4000.dat'
        IIIpath=mp+i+ions[2]+'/transit_m_1_t_4000.dat'
        if os.path.isfile(Ipath) == True:
            print(i)
            neutral = np.flipud(np.squeeze(fun.read_binary(Ipath,double=True)))
            wl_b,neutral_b = blur_HST(wl,neutral,dRV)
            fig, (a0,a1,a2) = plt.subplots(1,3,gridspec_kw = {'width_ratios':[3,1,1]},figsize=(10,3))
            a0.set_title = 'Spectrum '+i
            #if i != 'leendert':
                #a0.plot(wl_b,1.0-(Fe_b/Rs)**2,color='lightgray',linewidth=0.5,label='Full model')
            a0.plot(wl_b,1.0-(neutral_b/Rs)**2,alpha=0.5,linewidth=0.5,label=i+' I',color='blue')
            #if os.path.isfile(IIpath) == True:
            ion = np.flipud(np.squeeze(fun.read_binary(IIpath,double=True)))
            ion_inflated = (ion-cont)*3.0+cont

            wl_b,ion_b = blur_HST(wl,ion,dRV)
            wl_b5,ion_b5 = blur_HST(wl,ion_inflated,dRV)
            a0.plot(wl_b,1.0-(ion_b/Rs)**2,alpha=0.5,linewidth=0.5,label=i+' II',color='orange')

            wlt_neutral,T_neutral = blur_HST(wl,(neutral-cont)*(-1.0),10.0)
            wlt_ion,T_ion = blur_HST(wl,(ion-cont)*(-1.0),10.0)
            i_sn = interpolate.interp1d(wl_b,1.0-(neutral_b/Rs)**2.0)
            i_si = interpolate.interp1d(wl_b,1.0-(ion_b/Rs)**2.0)
            i_si5 = interpolate.interp1d(wl_b5,1.0-(ion_b5/Rs)**2.0)

            wldsel=wld[(wld >= wlmin) & (wld <= wlmax)]
            neutral_STIS =i_sn(wldsel)
            ion_STIS = i_si(wldsel)
            ion5_STIS = i_si5(wldsel)

            modelstack_all=[i_s(wldsel)]
            modelstack_neutral=[neutral_STIS]
            modelstack_ion=[ion_STIS]
            modelstack_ion5=[ion5_STIS]

            for j in range(0,nmodels):
                modelstack_all.append(i_s(wldsel)+np.random.normal(0.0,1.0/SNR,len(wldsel)))
                modelstack_neutral.append(neutral_STIS+np.random.normal(0.0,1.0/SNR,len(wldsel)))
                modelstack_ion.append(ion_STIS+np.random.normal(0.0,1.0/SNR,len(wldsel)))
                modelstack_ion5.append(ion5_STIS+np.random.normal(0.0,1.0/SNR,len(wldsel)))

            if i == 'Mg' or i == 'Si':
                RVrange = 3000.0
            else:
                RVrange = 300.0

            rv_n,ccf_n = an.xcor([wldsel],[np.stack(modelstack_neutral)],wlt_neutral,T_neutral,5.0,RVrange)
            rv_i,ccf_i = an.xcor([wldsel],[np.stack(modelstack_ion)],wlt_ion,T_ion,5.0,RVrange)
            rv_i5,ccf_i5 = an.xcor([wldsel],[np.stack(modelstack_ion5)],wlt_ion,T_ion,5.0,RVrange)

            ccf_avg_n = np.mean(ccf_n,axis=0)
            ccf_std_n = np.std(ccf_n,axis=0)
            ccf_avg_i = np.mean(ccf_i,axis=0)
            ccf_std_i = np.std(ccf_i,axis=0)
            a1.plot(rv_n,ccf_n[0,:],'--',color='black',linewidth=1.5)
            a2.plot(rv_i,ccf_i[0,:],'--',color='black',linewidth=1.5)
            a2.plot(rv_i5,ccf_i5[0,:]-ccf_i5[0,0]+ccf_i[0,0],color='black',linewidth=1.5)
            #a1.errorbar(rv_n,ccf_avg_n,fmt='.',yerr=ccf_std_n,color='gray')
            #a2.errorbar(rv_i,ccf_avg_i,fmt='.',yerr=ccf_std_i,color='gray')
            for k in range(1,nmodels):
                a1.plot(rv_n,ccf_n[k,:],'.',alpha=2.0/nmodels/2.0,color='blue')
                a2.plot(rv_i,ccf_i5[k,:]-ccf_i5[0,0]+ccf_i[0,0],'.',alpha=2.3/nmodels/2.0*1.5,color='orange')


            # a0.set_fontsize(16)
            # a1.set_fontsize(16)
            # a2.set_fontsize(16)
            a0.set_xlabel('Wavelength (nm)')
            a0.set_ylabel('Transit depth (relative to star)')
            a0.set_xlim((200.0,350.0))
            a1.set_ylabel('Line-averaged transit depth')
            a1.set_xlabel('RV (km/s)')
            a2.set_xlabel('RV (km/s)')
            a0.tick_params(axis='both',labelsize='7')
            a1.tick_params(axis='both',labelsize='7')
            a2.tick_params(axis='both',labelsize='7')
            #a2.set_yticklabels([])
            #a2.yaxis.set_visible(False)
            #a2.yaxis.tick_right()
            a0.legend(loc='best')
            fig.tight_layout()
            plt.subplots_adjust(wspace=0.28)
            fig.savefig(outp+i+'.pdf', bbox_inches='tight',dpi=500)
            plt.close(fig)






def plot_spec():

    import lib.functions as fun
    from astropy.io import fits
    import numpy as np
    import pdb
    from matplotlib import pyplot as plt


    dp='data/MASCARA-2/'
    n=30

    nrs=fun.findgen(1138-1073)+1073
    b=0
    ii=0
    for i in nrs:
        path=dp+'185603_180601.%s.spec.fits' % int(i)
        try:
            a=fits.getdata(path)
            b+=a
            ii+=1
        except:
            print('%s not found' % int(i))


    print(np.shape(b))
    plt.plot(b[1,n,:]/ii,b[0,n,:]/b[2,n,:])
    plt.show()

#plot_spec()
HST_proposal()
