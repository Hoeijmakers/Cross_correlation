import copy
import sys
#===============================================================================
#===============================================================================
#            THIS DEFINES THE PARAMETERS THAT CONTROL THE RUN
#===============================================================================
#===============================================================================



#===================================
#==========DATA AND MODELS==========
#===================================
dataname='Wasp-121/night3'
# dataname='Wasp-189/night1'
# dataname='WASP-49/night3'
# dataname='MASCARA-2'
# dataname='HD209458_HARPS/night1'

model_library = 'models/library_WASP_121_new'
# model_library = 'models/library_KELT_9'
# model_library = 'models/library'
# template_library = 'models/library_KELT_9'
# template_library = 'models/library'
template_library = 'models/library_WASP_121_new'




# templatename='phoenix-6200K'
# templatename=['ScI','ScII','CaI','FeI','FeII','SiI','NaI','MgI','TiI','TiII','CrI','CrII','MnI','CoI',
# templatename=['NiI','VI','VII','ZrI','NbI','AlI','SrII','YI','YII','RuI','CeII','NdII','BaII','LaII','EuII','FeI','FeII','TiI','TiII','NaI','MgI','CrI','CrII','VI','VII']
# templatename=['LiI','NaI','MgI','CaI','KI','ScI','ScII','TiI','TiII','VI','VII','CrI','CrII','MnI','MnII','FeII','NiI','NiII','CuI','CoI','CoII','ZnI','SrI','SrII','YI','YII','RbI']
# templatename=['FeI','VI','NaI','LiI','FeII','TiI','CaI','CrI']
# templatename='MgI'
# templatename=['CaI','SiI']

# for i in range(len(modelname)): modelname[i]+='_2500K'

templatename = ['CaI','CoI','CrI','CuI','FeI','GaI','GeI','H2O','KI','LiI','MgI','MnI','NaI','NiI','RbI','SH','ScI','SrI','TiI','TiO','VI','VO','YI','ZrI']#All templates that have lines in WASP-12b

# templatename=['FeI','CrI','VI','CaI','TiI','SH']
#templatename='comb_500_0.5'

#Use the following to read output
# import astropy.io.fits as fits
# import lib.analysis as an
# l=['CaI','CoI','CrI','CuI','FeI','GaI','GeI','H2O','KI','LiI','MgI','MnI','NaI','NiI','RbI','SH','ScI','SrI','TiI','TiO','VI','VO','YI','ZrI']
# inpath='output/Wasp-121/night1/library_WASP_121_new/'
# RV=fits.getdata(inpath+'FeI/RV.fits')
# Kp=fits.getdata(inpath+'FeI/Kp.fits')
# i=0
# print(l[i]); an.plot_KpVsys(RV,Kp,fits.getdata(inpath+l[i]+'/KpVsys.fits'),'data/Wasp-121/night1/',injected=fits.getdata(inpath+l[i]+'/KpVsys_i.fits'),xrange=[-250,250],invert=False)


# modelname=copy.deepcopy(templatename)
# for i in range(len(modelname)): modelname[i]+='_2500K'

# templatename=['SH_lorenzo','TiO_lorenzo']
# templatename='SH_lorenzo'
modelname='2000K_20_no_TiO'



#========================
#========SWITCHES========
#========================
do_telluric_correction=True
do_colour_correction=True
do_xcor=True#Set this to true if you want the CCF to be recomputed.
inject_model=True
#Set to False if you have already computed the CCF in a previous run, and now you
#just want to alter some plotting parameters.
plot_xcor=True
make_mask = False
apply_mask = True
do_berv_correction=True
do_keplerian_correction = True
make_doppler_model=False #Make a new doppler model (True) / use the previously generated one (False).
skip_doppler_model = False#This is skipping the application of the doppler model altogether.
RVrange=500.0
drv=1.0
#STILL NEED TO COPY THE CONFIGFILE TO THE XCOR OUTPUT, SO THAT WHEN DO_XCOR IS SET
#TO FALSE, THE CURRENT CONFIG FILE CAN BE CHECKED AGAINST THE PREVIOUS ONE (they should
#be identical, otherwise the cross-correlation should be repeated)

#===============================================
#==========SHADOW MODEL AND MASK NAMES==========
#===============================================
#These are names of shadow and mask files located in their respective data folders (dp). Vary them only if you are experimenting with different masks/shadows, and want to keep and access both.
shadowname='shadow_FeI'
maskname='generic_mask'
#maskname='mask_night1'










#The following is typically used to prepare data for analysis:
# from lib import read_data as rd
# rd.read_EXPRES_e2ds_v5('../MASCARA-2/MASCARA-2b_fifth_reduction/','MASCARA-2/')
# rd.read_HARPS_e2ds('../Raw_data/Wasp-121/Data/2018-01-14/','Wasp-121/night3/',molecfit = True)
# rd.read_HARPS_e2ds('../Raw_data/Wasp-189/night1/','Wasp-189/night1',molecfit=True)
# rd.read_HARPS_e2ds('../Raw_data/Wasp-49/night1/data/reduced/2015-12-06/','Wasp-49/night1/',nowave=True)
# rd.read_HARPS_e2ds('../Raw_data/Wasp-49/night2/data/reduced/2015-12-31/','Wasp-49/night2/',nowave=True)
# rd.read_HARPS_e2ds('../Raw_data/Wasp-49/night3/data/reduced/2016-01-14/','Wasp-49/night3/',nowave=True)
# rd.read_HARPS_e2ds('../Raw_data/Wasp-166/data/reduced/night2/','Wasp-166/night2/',nowave=True,molecfit = True)
# rd.read_HARPS_e2ds('../Raw_data/HD209458/night1/','HD209458_HARPS/night1/',nowave=True,molecfit = True,mode='HARPSN')
# sys.exit()
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================






















import lib.run
import lib.utils as ut

if isinstance(templatename,list) == True:
    N = len(templatename)
    print('Looping though %s templates.' % N)
    if plot_xcor == True:
        plot_xcor = False
        print('Setting plot_xcor to False, otherwise the loop will hang.')
    if isinstance(modelname,list) and len(modelname) != len(templatename):
        print('ERROR: Unequal number of models and templates provided.')
        sys.exit()
    for i in range(N):
        print('Starting '+templatename[i])

        if isinstance(modelname,list):
            modelname_i = modelname[i]
        else:
            modelname_i = modelname

        t1=ut.start()
        lib.run.run_instance(dataname,modelname_i,templatename[i],shadowname,maskname,RVrange=RVrange,drv=drv,do_colour_correction=do_colour_correction,do_xcor=do_xcor,plot_xcor=plot_xcor,do_berv_correction=do_berv_correction,do_keplerian_correction = do_keplerian_correction,make_doppler_model=make_doppler_model,template_library=template_library,inject_model=inject_model,model_library=model_library,skip_doppler_model=skip_doppler_model,make_mask=make_mask,apply_mask=apply_mask,do_telluric_correction=do_telluric_correction)
        ut.end(t1)
else:
    print('Starting '+templatename)
    t1=ut.start()
    lib.run.run_instance(dataname,modelname,templatename,shadowname,maskname,RVrange=RVrange,drv=drv,do_colour_correction=do_colour_correction,do_xcor=do_xcor,plot_xcor=plot_xcor,do_berv_correction=do_berv_correction,do_keplerian_correction = do_keplerian_correction,make_doppler_model=make_doppler_model,template_library=template_library,inject_model=inject_model,model_library=model_library,skip_doppler_model=skip_doppler_model,make_mask=make_mask,apply_mask=apply_mask,do_telluric_correction=do_telluric_correction)
    ut.end(t1)
