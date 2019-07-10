import lib.run
from lib import read_data as rd
import sys
#===============================================================================
#===============================================================================
#            THIS DEFINES THE PARAMETERS THAT CONTROL THE RUN
#===============================================================================
#===============================================================================

dataname='Wasp-121/night1'
# dataname='Wasp-189'
# dataname='WASP-49/night3'
#d ataname='MASCARA-2'
modelname='W121-Evans-noVO'
model_library = 'models/library'
# templatename='phoenix-6200K'
# templatename=['ScI','ScII','CaI','FeI','FeII','SiI','NaI','MgI','TiI','TiII','CrI','CrII','MnI','CoI',
# templatename=['NiI','VI','VII','ZrI','NbI','AlI','SrII','YI','YII','RuI','CeII','NdII','BaII','LaII','EuII','FeI','FeII','TiI','TiII','NaI','MgI','CrI','CrII','VI','VII']
templatename=['LiI','NaI','MgI','CaI','KI','ScI','ScII','TiI','TiII','VI','VII','CrI','CrII','MnI','MnII','FeII','NiI','NiII','CuI','CoI','CoII','ZnI','SrI','SrII','YI','YII','RbI']
for i in range(len(templatename)): templatename[i]+='_2500K'#This line and the above is for Wasp-121.
templatename='FeI'
# templatename='FeI'
#template_library = 'models/library_WASP_121'
template_library = 'models/library_KELT_9'
#templatename='W121-TiO'

do_colour_correction=True
do_xcor=True#Set this to true if you want the CCF to be recomputed.
#Set to False if you have already computed the CCF in a previous run, and now you
#just want to alter some plotting parameters.
plot_xcor=True
make_mask = True
apply_mask = True
do_berv_correction=True
do_keplerian_correction = True
#So tomorrow:
#-Rerun FeI without mask. See if it is messed up by turning mask on or off.azizamo boos konam

make_doppler_model=True #Make a new doppler model (True) ) use the previously generated one (False).
skip_doppler_model = False#This is skipping the application of the doppler model altogether.
RVrange=500.0
drv=1.0
#STILL NEED TO COPY THE CONFIGFILE TO THE XCOR OUTPUT, SO THAT WHEN DO_XCOR IS SET
#TO FALSE, THE CURRENT CONFIG FILE CAN BE CHECKED AGAINST THE PREVIOUS ONE (they should
#be identical, otherwise the cross-correlation should be repeated)


#The following is typically used to prepare data for analysis:
#Reading raw
#rd.read_HARPS_e2ds('../Wasp-121/Reduced_Data/night1/','Wasp-121/night1/')
#rd.read_HARPS_e2ds('../Wasp-189/2019-04-14/','Wasp-189')
#rd.read_EXPRES_e2ds('../MASCARA-2/MASCARA-2b_fourth_reduction/','MASCARA-2/')
# rd.read_HARPS_e2ds('../Wasp-49/night1/data/reduced/2015-12-06/','Wasp-49/night1/',nowave=True)
# rd.read_HARPS_e2ds('../Wasp-49/night2/data/reduced/2015-12-31/','Wasp-49/night2/',nowave=True)
# rd.read_HARPS_e2ds('../Wasp-49/night3/data/reduced/2016-01-14/','Wasp-49/night3/',nowave=True)
# sys.exit()

shadowname='shadow_Kelt9_Fe'#These are generic names of shadow and mask files located in their respective data folders (dp). Vary them only if you are experimenting with different masks/shadows, and want to keep and access both.
maskname='generic_mask'
# shadowname='models/shadow_wasp121_night3.pkl'
# maskname = 'models/mask_wasp_121_night3.pkl'

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
        lib.run.run_instance(dataname,modelname,templatename[i],shadowname,maskname,RVrange=RVrange,drv=drv,do_colour_correction=do_colour_correction,do_xcor=do_xcor,plot_xcor=plot_xcor,do_berv_correction=do_berv_correction,do_keplerian_correction = do_keplerian_correction,make_doppler_model=make_doppler_model,template_library=template_library,model_library=model_library,skip_doppler_model=skip_doppler_model,make_mask=make_mask,apply_mask=apply_mask)
else:
    print('Starting '+templatename)
    lib.run.run_instance(dataname,modelname,templatename,shadowname,maskname,RVrange=RVrange,drv=drv,do_colour_correction=do_colour_correction,do_xcor=do_xcor,plot_xcor=plot_xcor,do_berv_correction=do_berv_correction,do_keplerian_correction = do_keplerian_correction,make_doppler_model=make_doppler_model,template_library=template_library,model_library=model_library,skip_doppler_model=skip_doppler_model,make_mask=make_mask,apply_mask=apply_mask)
