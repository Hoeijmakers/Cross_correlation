def sector_14():
    import matplotlib.pyplot as plt
    import astropy.io.fits as fits
    import lib.utils as ut
    import numpy as np
    import lib.operations as ops
    import scipy.interpolate as interpolate
    p = 'Sector_14.fits'

    hdu=fits.open(p)
    img=hdu[1].data['FLUX']
    time=hdu[1].data['TIME']

    xc = 8-1
    yc = 6-1

    dpx = 1
    dt = 0.4

    period = 1.4811235


    print(type(img))

    # plt.plot(img[:,xc,yc])

    lc1 = []
    lc2 = []
    lc3 = []
    lc4 = []
    t1 = []
    t2 = []
    t3 = []
    t4 =[]

    for i in range(len(time)):
        if time[i] < 1692.23:
            lc1.append(np.sum(img[i,xc-1:xc+1,yc-1:yc+1]))
            t1.append(time[i])
        if time[i] > 1694.21 and time[i] < 1696.35:
            lc2.append(np.sum(img[i,xc-1:xc+1,yc-1:yc+1]))
            t2.append(time[i])
        if time[i] > 1697.66 and time[i] < 1706.41:
            lc3.append(np.sum(img[i,xc-1:xc+1,yc-1:yc+1]))
            t3.append(time[i])
        if time[i] > 1708.29:
            lc4.append(np.sum(img[i,xc-1:xc+1,yc-1:yc+1]))
            t4.append(time[i])


    plt.plot(t1,lc1,color='blue')
    plt.plot(t2,lc2,color='blue')
    plt.plot(t3,lc3,color='blue')
    plt.plot(t4,lc4,color='blue')

    tb1,eb1=ops.envelope(np.array(t1),np.array(lc1),dt,selfrac=0.8,mode='top')#Bins of the envelope
    tb2,eb2=ops.envelope(np.array(t2),np.array(lc2),dt,selfrac=0.8,mode='top')
    tb3,eb3=ops.envelope(np.array(t3),np.array(lc3),dt,selfrac=0.8,mode='top')
    tb4,eb4=ops.envelope(np.array(t4),np.array(lc4),dt,selfrac=0.8,mode='top')

    e1 = interpolate.interp1d(tb1,eb1,fill_value='extrapolate')(t1)
    e2 = interpolate.interp1d(tb2,eb2,fill_value='extrapolate')(t2)
    e3 = interpolate.interp1d(tb3,eb3,fill_value='extrapolate')(t3)
    e4 = interpolate.interp1d(tb4,eb4,fill_value='extrapolate')(t4)


    plt.plot(t1,e1,color='red')
    plt.plot(t2,e2,color='red')
    plt.plot(t3,e3,color='red')
    plt.plot(t4,e4,color='red')
    plt.xlabel('Time (days)')
    plt.ylabel('Raw flux')
    plt.show()
    fxn1 = np.array(lc1/e1)
    fxn2 = np.array(lc2/e2)
    fxn3 = np.array(lc3/e3)
    fxn4 = np.array(lc4/e4)

    # plt.plot(t1,fxn1)
    # plt.plot(t2,fxn2)
    # plt.plot(t3,fxn3)
    # plt.plot(t4,fxn4)
    # plt.show()
    tc = np.concatenate((t1,t2,t3,t4))
    fxc = np.concatenate((fxn1,fxn2,fxn3,fxn4))


    sel = (fxc > 0.95)

    tc2 = tc[(fxc >0.95)]
    fx2 = fxc[(fxc >0.95)]

    print(len(tc))
    print(len(fx2))


    phase = (tc2 % period)/period -0.6
    phase[(phase < -0.5)]+=1
    plt.plot(phase,fx2,'.',alpha=0.5)
    plt.axhline(1.0)
    plt.xlabel('Orbital phase')
    plt.ylabel('Normalized flux')
    plt.show()
    # ut.writefits('test.fits',img)
sector_14()

def plot_KpVsyses_Wasp_121():
    import astropy.io.fits as fits
    import lib.analysis as an
    Lt = ['LiI','NaI','MgI','CaI','KI','ScI','ScII','TiI','TiII','VI','VII','CrI','CrII','MnI','MnII','FeI','FeII','NiI','NiII','CuI','CoI','CoII','ZnI','SrI','SrII','YI','YII','RbI']
    Ld = ['Wasp-121/night1','Wasp-121/night2','Wasp-121/night3']
    Lm = ['KI','ScI','MnI','NiI','CuI','CoI','ZnI','YI']
    for i in range(len(Lm)):
        Lm[i]+='_2500K'
    ts = Lt[0]+'_2500K'
    Kp = fits.getdata('output/'+Ld[0]+'/'+ts+'/Kp.fits')
    RV = fits.getdata('output/'+Ld[0]+'/'+ts+'/RV.fits')
    dp = 'data/Wasp-121/night3/' #But should be the same for all three nights.
    N = len(Lt)
    KpVsys = an.combine_KpVsys(Ld,Lm)
    an.plot_KpVsys(RV,Kp,KpVsys,dp,xrange=[-100,100],Nticks = 10.0,title='Metals',invert=True)

    for i in range(N):
        KpVsys = an.combine_KpVsys(Ld,[Lt[i]+'_2500K'])
        an.plot_KpVsys(RV,Kp,KpVsys,dp,xrange=[-100,100],Nticks = 10.0,title=Lt[i],invert=True)

def plot_model(species):
    import matplotlib.pyplot as plt
    import lib.models as m
    name=species+'_2500K'
    wl,fx=m.get_model(name,library='models/library_WASP_121')
    plt.plot(wl,fx)
    plt.show()


def convert_models_daniel_wasp121():
    import lib.models as m
    p = '/Volumes/TOSHIBA/WASP-121/Grid_Daniel_Wasp-121/Xcor_templates/templates_1500K/'
    Rstar = 1.458#Rsun
    m.convert_models_daniel(p,Rstar,'models/WASP-121/templates_daniel/','1500K_1_')
    p = '/Volumes/TOSHIBA/WASP-121/Grid_Daniel_Wasp-121/Xcor_templates/templates_2500K/'
    m.convert_models_daniel(p,Rstar,'models/WASP-121/templates_daniel/','2500K_1_')
# convert_models_daniel_wasp121()

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



def plot_model_presentation():
    import lib.models as m
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.interpolate as int
    import pdb
    from scipy.signal import savgol_filter
    import lib.utils as ut

    t1 = ut.start()
    wl,fx1 = m.get_model('W121-TiO')
    wl2,fx2 = m.get_model('FeII',library='models/library_KELT_9')
    wl3,fx3 = m.get_model('OII',library='models/library_KELT_9')
    wl4,fx4 = m.get_model('NaI',library='models/library_KELT_9')

    fx1c = fx1 - np.max(fx1)
    fx2c = fx2 - np.max(fx2)
    fx3c = fx3 - np.max(fx3)
    fx4c = fx4 - np.max(fx4)

    fx2c_i = int.interp1d(wl2,fx2c)
    fx3c_i = int.interp1d(wl3,fx3c)
    fx4c_i = int.interp1d(wl4,fx4c)

    fx2c = fx2c_i(wl)
    fx3c = fx3c_i(wl)
    fx4c = fx4c_i(wl)

    fx3c = np.flip(fx3c*15.0,0)

    fx_final = []

    for i in range(len(wl)):
        fx_final.append(np.min([fx1c[i],fx2c[i],fx3c[i],fx4c[i]]))

    dx = 15000
    fx_final = np.array(fx_final)


    wlb = []
    fxb = []
    i = 0
    while (i+1)*dx < len(wl):
        wlb.append(np.mean(wl[i*dx:(i+1)*dx]))
        fxb.append(np.mean(fx_final[i*dx:(i+1)*dx]*100.0))
        i+=1

    t2 = ut.end(t1)
    plt.rcParams.update({'font.size': 60})
    fig = plt.figure(figsize=(40,30))
    ax = fig.add_subplot(111)
    ax.plot(wl,fx_final*100.0,linewidth=3.0,color='tomato')
    ax.plot(wlb,fxb,'.',color='black',ms=50)
    ax.plot(wlb,fxb,color='black',linewidth = 10)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Transit depth (%)')
    ax.set_xlim((400,700))
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(6.0)
    plt.savefig('Wasp_121_demo.png', transparent=True)


#plot_spec()
#HST_proposal()
