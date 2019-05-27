import lib.run
from lib import read_data as rd
import sys
#===============================================================================
#===============================================================================
#            THIS DEFINES THE PARAMETERS THAT CONTROL THE RUN
#===============================================================================
#===============================================================================


#dataname='Wasp-121/night3'
dataname='Wasp-189'
# dataname='WASP-49/night3'
#dataname='MASCARA-2'
modelname='W121-Evans-noVO'
model_library = 'models/library'
#templatename='phoenix-6200K'
templatename=['ScI','ScII','CaI','FeI','FeII','SiI','NaI','MgI','TiI','TiII','CrI','CrII','MnI','CoI','NiI','VI','VII','ZrI','NbI','AlI','SrII','YI','YII','RuI','CeII','NdII','BaII','LaII','EuII']
# templatename=['FeI','FeII','TiI','TiII','NaI','MgI','CrI','CrII','VI','VII']#Restricted list.

#templatename='FeI'
template_library = 'models/library_KELT_9'
shadowname='models/shadow_wasp_189.pkl' #THERE MAY BE HELL TO PAY IF I SET THIS TO AN EXISTING FOLDER.
#RATHER THAN A FILE. YOU NEED TO CATCH THAT!
#templatename='W121-TiO'

do_colour_correction=True
do_xcor=True#Set this to true if you want the CCF to be recomputed.
#Set to False if you have already computed the CCF in a previous run, and now you
#just want to alter some plotting parameters.
plot_xcor=True
do_berv_correction=True
do_keplerian_correction = True
make_doppler_model=False #Make a new doppler model or use the previously generated one.
skip_doppler_model = False#This is skipping the application of the doppler model altogether.
RVrange=500.0
drv=1.0


#These two are dataset specific and need to go to CONFIG.
#air=False#Set to true if the data wl is in air.
#fitdv=12.0
#fitdv=100 #for Wasp-189
#startorder = 2#For HARPS
#endorder = 70#For HARPS
#startorder = None #For EXPRES
#endorder = None #For EXPRES

#The following is typically used to prepare data for analysis:
#Reading raw
#rd.read_HARPS_e2ds('../Wasp-121/Reduced_Data/night1/','Wasp-121/night1/')
#rd.read_HARPS_e2ds('../Wasp-189/2019-04-14/','Wasp-189')
#rd.read_EXPRES_e2ds('../MASCARA-2/MASCARA-2b_fourth_reduction/','MASCARA-2/')
# rd.read_HARPS_e2ds('../Wasp-49/night1/data/reduced/2015-12-06/','Wasp-49/night1/',nowave=True)
# rd.read_HARPS_e2ds('../Wasp-49/night2/data/reduced/2015-12-31/','Wasp-49/night2/',nowave=True)
# rd.read_HARPS_e2ds('../Wasp-49/night3/data/reduced/2016-01-14/','Wasp-49/night3/',nowave=True)
# sys.exit()
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================
#===============================================================================

if isinstance(templatename,list) == True:
    N = len(templatename)
    print('Looping though %s templates.' % N)
    if plot_xcor == True:
        plot_xcor = False
        print('Setting plot_xcor to False, otherwise the loop will hang.')
    for i in range(N):
        print('Starting '+templatename[i])
        lib.run.run_instance(dataname,modelname,templatename[i],shadowname,RVrange=RVrange,drv=drv,do_colour_correction=do_colour_correction,do_xcor=do_xcor,plot_xcor=plot_xcor,do_berv_correction=do_berv_correction,do_keplerian_correction = do_keplerian_correction,make_doppler_model=make_doppler_model,template_library=template_library,model_library=model_library,skip_doppler_model=skip_doppler_model)
else:
    print('Starting '+templatename)
    lib.run.run_instance(dataname,modelname,templatename,shadowname,RVrange=RVrange,drv=drv,do_colour_correction=do_colour_correction,do_xcor=do_xcor,plot_xcor=plot_xcor,do_berv_correction=do_berv_correction,do_keplerian_correction = do_keplerian_correction,make_doppler_model=make_doppler_model,template_library=template_library,model_library=model_library,skip_doppler_model=skip_doppler_model)
