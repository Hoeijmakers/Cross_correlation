#This package contains high-level wrappers for running the entire sequence.



def run_instance(dataname,modelname,templatename,shadowname,RVrange=500.0,drv=1.0,do_colour_correction=True,do_xcor=True,plot_xcor=False,do_berv_correction=True,do_keplerian_correction = True,make_doppler_model=True,model_library='models/library',template_library='models/library',skip_doppler_model = False):
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
    from astropy.io import fits
    from matplotlib import pyplot as plt
    from scipy import interpolate
    import pylab
    import pdb
    import os.path
    import os
    import sys
    import glob
    import distutils.util

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
    if do_xcor == True or plot_xcor == True:
        print('Reading orders from '+dp)
        t1=ut.start()
        for i in range(startorder,endorder+1):
            if air == False:
                if trigger1 ==0:
                    trigger1=-1
                    print("Assuming wavelengths are in vaccuum.")
                list_of_wls.append(fits.getdata(dp+'wave_%s.fits' % i))
            else:
                if trigger1 ==0:
                    print("Applying airtovac correction.")
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
                    print('WARNING: SIGMA FILE NOT PROVIDED. ASSUMING SIGMA = SQRT(FLUX). This is ok for HARPS data.')
                    trigger2=-1
                list_of_sigmas.append(np.sqrt(order_i))
        t2=ut.end(t1)

    t1 = ut.start()
    list_of_orders = cleaning.mask_orders(list_of_orders,40.0,5.0)
    t2 = ut.end(t1)
    print("NEED TO ADD A GUI IN HERE THAT ENABLES MASKING!")
    print("And in the end, the adopted output of all these guis needs to be collected")
    print("into a single sort of report, so that at the time of paper writing, its easy")
    print(" to see which numbers were used.")


#The following tests inject_model
# print('Injecting model')
# t1=ut.start()
# order_injected=models.inject_model(wl1,order1,dp,modelname)
# t2=ut.end(t1)


#Need to normalize the orders to their average flux in order to effectively apply
#a broad-band colour correction (colour is a function of airmass and seeing)
    if do_xcor == True and do_colour_correction == True:
        print('Normalizing orders to common flux level')
        list_of_orders = ops.normalize_orders(list_of_orders)
    #The following performs continuum normalization for template construction.
    if do_xcor == True or plot_xcor == True:
        print('Building template')
        t1=ut.start()
        wlt,T=models.build_template(templatename,binsize=0.5,maxfrac=0.01,resolution=120000.0,template_library=template_library)
        t2=ut.end(t1)
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
        print('Starting cross-correlation')
        t1=ut.start()
        rv,ccf,ccf_e,Tsums=analysis.xcor(list_of_wls,list_of_orders,np.flipud(np.flipud(wlt)),T,drv,RVrange,list_of_errors=list_of_sigmas)
        t2=ut.end(t1)
        print('Writing output to '+outpath)
        ut.writefits(outpath+'ccf.fits',ccf)
        ut.writefits(outpath+'ccf_e.fits',ccf_e)
        ut.writefits(outpath+'RV.fits',rv)
        ut.writefits(outpath+'Tsum.fits',Tsums)
    else:
        print('Reading CCFs from '+outpath)
        try:
            f = open(outpath+'ccf.fits', 'r')
        except FileNotFoundError:
            print('--- ERROR: Necessary output not located at '+outpath+'. Set do_xcor to True.')
            sys.exit()
        rv=fits.getdata(outpath+'rv.fits')
        ccf = fits.getdata(outpath+'ccf.fits')
        ccf_e = fits.getdata(outpath+'ccf_e.fits')
        Tsums = fits.getdata(outpath+'Tsum.fits')




#Velocity corrections
    rv_cor = 0
    if do_berv_correction == True:
        rv_cor += sp.berv(dp)
    if do_keplerian_correction == True:
        rv_cor-=sp.RV_star(dp)

    if type(rv_cor) != int:
        print('Performing velocity corrections')
        ccf_cor = ops.shift_ccf(rv,ccf,rv_cor)
        ccf_e_cor = ops.shift_ccf(rv,ccf_e,rv_cor)
    else:
        ccf_cor = ccf*1.0
        ccf_e_cor = ccf_e*1.0


    if plot_xcor == True:
        print('Plotting orders and XCOR')
        fitdv = sp.paramget('fitdv',dp)
        analysis.plot_XCOR(list_of_wls,list_of_orders,wlt,T,rv,ccf_cor,Tsums,dp,CCF_E=ccf_e_cor,dv=fitdv)



    print('Cleaning CCFs')
    ccf_n,ccf_ne,ccf_nn = cleaning.clean_ccf(rv,ccf_cor,ccf_e_cor,dp)
    ut.writefits(outpath+'ccf_normalized.fits',ccf_nn)
    ut.writefits(outpath+'ccf_ne.fits',ccf_ne)

    # print('Building and removing Doppler Model')
    # print('THIS IS STILL HARDCODED TO WASP-121!')
    # doppler_rv,ccf_ds_model = cleaning.construct_doppler_model_vincent(dp+'../shadow_model_vincent.dat.txt',dp,rv,ccf_nn)
    # ccf_clean = ccf_nn - ccf_ds_model


    if make_doppler_model == True and skip_doppler_model == False:
        cleaning.construct_doppler_model(rv,ccf_nn,dp,shadowname,xrange=[-200,200],Nxticks=20.0,Nyticks=10.0)

    if skip_doppler_model == False:
        print('Reading model from '+shadowname)
        doppler_model,maskHW = cleaning.read_shadow(shadowname,rv,ccf)
        ccf_clean,matched_ds_model = cleaning.match_shadow(rv,ccf_nn,dp,doppler_model,maskHW)
    else:
        ccf_clean = ccf_nn*1.0
        matched_ds_model = ccf_clean*0.0
    ut.save_stack(outpath+'cleaning_steps.fits',[ccf,ccf_cor,ccf_nn,ccf_clean,matched_ds_model])
    ut.writefits(outpath+'ccf_cleaned.fits',ccf_clean)

    if plot_xcor == True:
        print('Plotting 2D CCF')
        print("THIS NEEDS TO BE REVAMPED!")
        # analysis.plot_ccf(rv,ccf_nn,dp,xrange=[-200,200],Nticks=20.0,doppler_model = doppler_rv)
        # analysis.plot_ccf(rv,ccf_ds_model,dp,xrange=[-200,200],Nticks=20.0,doppler_model = doppler_rv)
        # analysis.plot_ccf(rv,ccf_clean,dp,xrange=[-200,200],Nticks=20.0,doppler_model = doppler_rv)


    print('Constructing KpVsys')
    t1=ut.start()
    Kp,KpVsys = analysis.construct_KpVsys(rv,ccf_clean,dp)
    ut.end(t1)
    ut.writefits(outpath+'KpVsys.fits',KpVsys)
    ut.writefits(outpath+'Kp.fits',Kp)


    if plot_xcor == True:
        print('Plotting KpVsys')
        analysis.plot_KpVsys(rv,Kp,KpVsys,dp)



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


#The following tests get_model.
#wlm,fxm=models.get_model('W121-plez')
#plot(wlm,fxm)
#show()



#This times the gaussian function.
#x=fun.findgen(10000)
#start = time.time()
#for i in range(0,10000):
#    g=fun.gaussian(x,20.0,50000.0,10000.0)
#end = time.time()
#print(end - start)
#plot(x,g)
#show()
