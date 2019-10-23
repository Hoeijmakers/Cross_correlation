#This package contains high-level wrappers for running the entire sequence.



def run_instance(dataname,modelname,templatename,shadowname,maskname,RVrange=500.0,drv=1.0,do_colour_correction=True,do_xcor=True,plot_xcor=False,do_berv_correction=True,do_keplerian_correction = True,make_doppler_model=True,inject_model = False,model_library='models/library',template_library='models/library',skip_doppler_model = False,make_mask=False,apply_mask=True,do_telluric_correction=False):
    """This is the main script that runs through the entire sequence of steps."""
    import numpy as np
    from lib import utils as ut
    from lib import operations as ops
    from lib import functions as fun
    from lib import system_parameters as sp
    from lib import models
    from lib import analysis
    from lib import constants as const
    from lib import cleaning
    from lib import read_data as rd
    from lib import system_parameters as sp
    from lib import masking as masking
    from lib import shadow as shadow
    from lib import molecfit as telcor
    from astropy.io import fits
    from matplotlib import pyplot as plt
    import scipy.interpolate as interp
    import pylab
    import pdb
    import os.path
    import os
    import sys
    import glob
    import distutils.util
    import pickle
    import copy


#Need to build in a lot of tests here on the input.


#We start by defining constants, preparing for generating output.
    c=const.c/1000.0
    dp = 'data/'+dataname+'/'
    outpath='output/'+dataname+'/'+templatename+'/'

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    list_of_wls=[]
    list_of_orders=[]
    list_of_sigmas=[]
    trigger1 = 0#These triggers are used to generate warning messages in the forloop.
    trigger2 = 0


    startorder = int(sp.paramget('startorder',dp))
    endorder = int(sp.paramget('endorder',dp))
    air = bool(distutils.util.strtobool(sp.paramget('air',dp)))#Read bool from str in config file.



#Loading the data from the datafolder.
    if do_xcor == True or plot_xcor == True or make_mask == True:
        print('---Reading orders from '+dp)

        for i in range(startorder,endorder+1):
            if air == False:
                if i == startorder:
                    trigger1=-1
                    print("------Assuming wavelengths are in vaccuum.")
                list_of_wls.append(fits.getdata(dp+'wave_%s.fits' % i))
            else:
                if i == startorder:
                    print("------Applying airtovac correction.")
                    trigger1 =-1
                list_of_wls.append(ops.airtovac(fits.getdata(dp+'wave_%s.fits' % i)))
            order_i = fits.getdata(dp+'order_%s.fits' % i)
            order_i[order_i <= 0] = np.nan
            ut.postest(order_i,varname='order %s in main' %i)
            list_of_orders.append(order_i)
            try:
                list_of_sigmas.append(fits.getdata(dp+'sigma_%s.fits' % i))
            except FileNotFoundError:
                if trigger2 == 0:
                    print('------WARNING: Sigma (flux error) files not provided. Assuming sigma = sqrt(flux). This is standard practise for HARPS data.')
                    trigger2=-1
                list_of_sigmas.append(np.sqrt(order_i))



#Apply telluric correction file or not.

    # plt.plot(list_of_wls[0],list_of_orders[56][0])
    # for i in range(57,69):
    #     plt.plot(list_of_wls[i],list_of_orders[i][0])
    # plt.show()
    if do_telluric_correction == True and len(list_of_orders) > 0:
        print('---Applying telluric correction')
        telpath = dp+'telluric_transmission_spectra.pkl'
        list_of_orders = telcor.apply_telluric_correction(telpath,list_of_wls,list_of_orders)
    # plt.plot(list_of_wls[0],list_of_orders[56][0])
    # plt.show()
    # pdb.set_trace()

#Do velocity correction of wl-solution. Explicitly after telluric correction
#but before masking! Because the cross-correlation relies on columns being masked.
#Then if you start to move the CCFs around before removing the time-average,
#each masked column becomes slanted. Bad deal.

    rv_cor = 0
    if do_berv_correction == True:
        rv_cor += sp.berv(dp)
    if do_keplerian_correction == True:
        rv_cor-=sp.RV_star(dp)*(1.0)

    if type(rv_cor) != int and len(list_of_orders) > 0:
        print('---Reinterpolating data to correct velocities')
        list_of_orders_cor = []
        list_of_sigmas_cor = []
        for i in range(len(list_of_wls)):
            order = list_of_orders[i]
            sigma = list_of_sigmas[i]
            order_cor = order*0.0
            sigma_cor = sigma*0.0
            for j in range(len(list_of_orders[0])):
                wl_i = interp.interp1d(list_of_wls[i],order[j],bounds_error=False)
                si_i = interp.interp1d(list_of_wls[i],sigma[j],bounds_error=False)
                wl_cor = list_of_wls[i]*(1.0-rv_cor[j]*1000.0/const.c)#The minus sign was tested on a slow-rotator.
                order_cor[j] = wl_i(wl_cor)
                sigma_cor[j] = si_i(wl_cor)
            list_of_orders_cor.append(order_cor)
            list_of_sigmas_cor.append(sigma_cor)
            ut.statusbar(i,fun.findgen(len(list_of_wls)))
        list_of_orders = list_of_orders_cor
        list_of_sigmas = list_of_sigmas_cor




#Do masking or not.
    if make_mask == True and len(list_of_orders) > 0:
        print('---Constructing mask ')
        masking.mask_orders(list_of_wls,list_of_orders,dp,maskname,40.0,5.0,manual=True)
        if apply_mask == False:
            print('---Warning in run_instanace: Mask was made but is not applied to data (apply_mask = False)')
    if apply_mask == True and len(list_of_orders) > 0:
        print('---Applying mask')
        list_of_orders = masking.apply_mask_from_file(dp,maskname,list_of_orders)
    print('---Healing NaNs')
    list_of_orders = masking.interpolate_over_NaNs(list_of_orders)



#Inject_model, or not.
    if inject_model == True and do_xcor == True:
        print('---Injecting model')
        list_of_orders_injected=models.inject_model(list_of_wls,list_of_orders,dp,modelname,model_library=model_library)


#Normalize the orders to their average flux in order to effectively apply
#a broad-band colour correction (colour is a function of airmass and seeing).
    if do_xcor == True and do_colour_correction == True:
        print('---Normalizing orders to common flux level')
        list_of_orders = ops.normalize_orders(list_of_orders)
        if inject_model == True:
            list_of_orders_injected = ops.normalize_orders(list_of_orders_injected)


#Construct the cross-correlation template in case we will be doing or plotting xcor.
    if do_xcor == True or plot_xcor == True:
        print('---Building template')
        wlt,T=models.build_template(templatename,binsize=0.5,maxfrac=0.01,resolution=120000.0,template_library=template_library)
        T*=(-1.0)






#CAN'T I OUTSOURCE THIS TO DANIEL AND SIMON? I mean for each spectrum they could also
#generate a binary mask by putting delta's instead of.. well... not precisely
#because of masking by line wings; that's the whole trick of opacity calculation.
#And then Daniel can subtract the continuum himself...
#I can rewrite this to include the name pointing to another model to be
#used as continuum normalization.

#The following tests XCOR on synthetic orders.
#wlm = fun.findgen(2e6)/10000.0+650.0
#fxm = wlm*0.0
#fxm[[(fun.findgen(400)*4e3+1e3).astype(int)]] = 1.0
#px_scale=wlm[1]-wlm[0]
#dldv = np.min(wlm) / c /px_scale
#T=ops.smooth(fxm,dldv * 5.0)
#fxm_b=ops.smooth(fxm,dldv * 20.0,mode='gaussian')
#plt.plot(wlm,T)
#plt.plot(wlm,fxm_b)
#plt.show()
#ii = interpolate.interp1d(wlm,fxm_b)
#dspec = ii(wl)
#order = fun.rebinreform(dspec/np.max(dspec),30)
#fits.writeto('test.fits',order,overwrite=True)
#ccf=analysis.xcor([wl,wl,wl+0.55555555555],[order,order*3.0,order*5.0],wlm,T,drv,RVrange)

# meanspec=np.mean(order1,axis=0)
# meanspec-=min(meanspec)
# meanspec/=np.max(meanspec)
# T-=np.min(T)
# T/=np.median(T)
# T-=np.max(T[(wlt >= min(wl1)) & (wlt <= max(wl1))])
# plt.plot(wl1,meanspec)
# plt.plot(wlt*(1.0-200.0/300000.0),T)
# plt.xlim((min(wl1),max(wl1)))
# plt.show()


# plt.plot(wlt,T)
# plt.show()



#Perform the cross-correlation on the entire list of orders.
    if do_xcor == True:
        print('---Cross-correlating spectra')
        rv,ccf,ccf_e,Tsums=analysis.xcor(list_of_wls,list_of_orders,np.flipud(np.flipud(wlt)),T,drv,RVrange,list_of_errors=list_of_sigmas)
        print('------Writing output to '+outpath)
        ut.writefits(outpath+'ccf.fits',ccf)
        ut.writefits(outpath+'ccf_e.fits',ccf_e)
        ut.writefits(outpath+'RV.fits',rv)
        ut.writefits(outpath+'Tsum.fits',Tsums)


        if inject_model == True:
            print('---Cross-correlating model injection')
            rv_i,ccf_i,ccf_e_i,Tsums_i=analysis.xcor(list_of_wls,list_of_orders_injected,np.flipud(np.flipud(wlt)),T,drv,RVrange,list_of_errors=list_of_sigmas)
            print('------Writing output to '+outpath)
            ut.writefits(outpath+'ccf_i.fits',ccf_i)
            ut.writefits(outpath+'ccf_e_i.fits',ccf_e_i)

    else:
        print('---Reading CCFs from '+outpath)
        try:
            f = open(outpath+'ccf.fits', 'r')
        except FileNotFoundError:
            print('------ERROR: Necessary CCF output not located at '+outpath+'. Set do_xcor to True.')
            sys.exit()
        rv=fits.getdata(outpath+'rv.fits')
        ccf = fits.getdata(outpath+'ccf.fits')
        ccf_e = fits.getdata(outpath+'ccf_e.fits')
        Tsums = fits.getdata(outpath+'Tsum.fits')

        if inject_model == True:
            print('---Reading injected CCFs from '+outpath)
            try:
                f = open(outpath+'ccf_i.fits', 'r')
            except FileNotFoundError:
                print('------ERROR: Necessary CCF output not located at '+outpath+'. Set do_xcor and inject_model to True.')
                sys.exit()
            ccf_i = fits.getdata(outpath+'ccf_i.fits')
            ccf_e_i = fits.getdata(outpath+'ccf_e_i.fits')


#Velocity corrections
    # rv_cor = 0
    # if do_berv_correction == True:
    #     rv_cor += sp.berv(dp)
    # if do_keplerian_correction == True:
    #     rv_cor-=sp.RV_star(dp)

    # if type(rv_cor) != int:
    #     print('---Performing velocity corrections')
    #     ccf_cor = ops.shift_ccf(rv,ccf,rv_cor)
    #     ccf_e_cor = ops.shift_ccf(rv,ccf_e,rv_cor)
    # else:
    ccf_cor = ccf*1.0
    ccf_e_cor = ccf_e*1.0


    if plot_xcor == True:
        print('---Plotting orders and XCOR')
        fitdv = sp.paramget('fitdv',dp)
        analysis.plot_XCOR(list_of_wls,list_of_orders,wlt,T,rv,ccf_cor,Tsums,dp,CCF_E=ccf_e_cor,dv=fitdv)

    sys.exit()
    print('---Cleaning CCFs')
    ccf_n,ccf_ne,ccf_nn,ccf_nne = cleaning.clean_ccf(rv,ccf_cor,ccf_e_cor,dp)
    ut.writefits(outpath+'ccf_normalized.fits',ccf_nn)
    ut.writefits(outpath+'ccf_ne.fits',ccf_ne)

    if inject_model == True:
        print('---Cleaning injected CCFs')
        ccf_n_i,ccf_ne_i,ccf_nn_i,ccf_nne_i = cleaning.clean_ccf(rv,ccf_i,ccf_e_i,dp)
        ut.writefits(outpath+'ccf_normalized_i.fits',ccf_nn_i)
        ut.writefits(outpath+'ccf_ne_i.fits',ccf_ne_i)

    if make_doppler_model == True and skip_doppler_model == False:
        shadow.construct_doppler_model(rv,ccf_nn,dp,shadowname,xrange=[-200,200],Nxticks=20.0,Nyticks=10.0)

    if skip_doppler_model == False:
        print('---Reading doppler shadow model from '+shadowname)
        doppler_model,maskHW = shadow.read_shadow(dp,shadowname,rv,ccf)
        ccf_clean,matched_ds_model = shadow.match_shadow(rv,ccf_nn,dp,doppler_model,maskHW)

        if inject_model == True:
                ccf_clean_i,matched_ds_model_i = shadow.match_shadow(rv,ccf_nn_i,dp,doppler_model,maskHW)
    else:
        print('---Not performing shadow correction')
        ccf_clean = ccf_nn*1.0
        matched_ds_model = ccf_clean*0.0
        if inject_model == True:
            ccf_clean_i = ccf_nn_i*1.0
            matched_ds_model_i = ccf_clean_i*0.0


    ut.save_stack(outpath+'cleaning_steps.fits',[ccf,ccf_cor,ccf_nn,ccf_clean,matched_ds_model])
    ut.writefits(outpath+'ccf_cleaned.fits',ccf_clean)

    if inject_model == True:
        ut.writefits(outpath+'ccf_cleaned_i.fits',ccf_clean_i)

    if plot_xcor == True:
        print('---Plotting 2D CCF')
        print("---THIS NEEDS TO BE REVAMPED!")
        # analysis.plot_ccf(rv,ccf_nn,dp,xrange=[-200,200],Nticks=20.0,doppler_model = doppler_rv)
        # analysis.plot_ccf(rv,ccf_ds_model,dp,xrange=[-200,200],Nticks=20.0,doppler_model = doppler_rv)
        # analysis.plot_ccf(rv,ccf_clean,dp,xrange=[-200,200],Nticks=20.0,doppler_model = doppler_rv)

    # ut.save_stack('test.fits',[ccf_n,ccf_nn,ccf_ne,ccf_nne])
    # pdb.set_trace()

    print('---Constructing KpVsys')
    Kp,KpVsys,KpVsys_e = analysis.construct_KpVsys(rv,ccf_clean,ccf_nne,dp)
    ut.writefits(outpath+'KpVsys.fits',KpVsys)
    ut.writefits(outpath+'KpVsys_e.fits',KpVsys_e)
    ut.writefits(outpath+'Kp.fits',Kp)
    if inject_model == True:
        Kp,KpVsys_i,KpVsys_e_i = analysis.construct_KpVsys(rv,ccf_clean_i,ccf_nne_i,dp)
        ut.writefits(outpath+'KpVsys_i.fits',KpVsys_i)
        # ut.writefits(outpath+'KpVsys_e_i.fits',KpVsys_e_i)


    if plot_xcor == True:
        print('---Plotting KpVsys')
        if inject_model == False:
            analysis.plot_KpVsys(rv,Kp,KpVsys,dp)
        else:
            analysis.plot_KpVsys(rv,Kp,KpVsys,dp,injected=KpVsys_i)



# analysis.plot_RV_star(dp,rv,ccf,RVrange=[-50,100]) #THis entire function is probably obsolete.

# sel = (doppler_rv > -100000.0)
# ccf_ds = ops.shift_ccf(rv,ccf_nn[sel,:],(-1.0)*doppler_rv[sel])
# vsys = sp.paramget('vsys',dp)
# sel = ((rv >= -20) & (rv <= 60))
#
# plt.plot(rv[sel],np.nanmean(ccf_ds[:,sel],axis=0))
# # plt.axvline(x=vsys)
# plt.show()
# pdb.set_trace()

#The following reads in a binary mask and plots it to search for the velocity shift
#by comparing with a Kelt-9 model.
#start = time.time()
#wlb,fxb=models.read_binary_mask('models/A0v2.mas')
#end = time.time()
#print(end-start)
#FeIImodel=fits.getdata('models/FeII_4500_c.fits')
#wlm=FeIImodel[0,:]
#fxm=FeIImodel[1,:]
#pylab.plot(wlm,300.0*(fxm-np.median(fxm)))
#pylab.plot(wlb,fxb)
#pylab.show()



#The following tests blurring.
#spec = fun.findgen(len(wl))*0.0
#spec[500] = 1
#spec[1000] = 1
#spec[3000] = 1
#spec[3012] = 1
#spec_b=ops.blur_rotate(wl,spec,3.0,1.5,1.5,90.0)
#spec_b2=ops.blur_rotate(wl,spec,3.0,2.0,1.0,90.0)
#t1=ut.start()
#spec_b=ops.blur_spec(wl,spec,20.0,mode='box')
#spec_b2=ops.blur_spec(wl,spec,20.0,mode='box')
#dt1=ut.end(t1)
#pylab.plot(wl,spec)
#pylab.plot(wl,spec_b)
#pylab.plot(wl,spec_b2)
#pylab.show()




#This times the gaussian function.
#x=fun.findgen(10000)
#start = time.time()
#for i in range(0,10000):
#    g=fun.gaussian(x,20.0,50000.0,10000.0)
#end = time.time()
#print(end - start)
#plot(x,g)
#show()
