import lib.utils as ut
import numpy as np
import lib.operations as ops
from lib import functions as fun
from lib import system_parameters as sp
from lib import models
from lib import analysis
from lib import constants as const
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy import interpolate
import pylab
import pdb

c=const.c/1000.0
dp='data/KELT-9/night1/'
modelname='W121-plez'
# templatename='phoenix-6200K'
templatename='Fe-test'
RVrange=400.0
drv=1.0

filepath=dp+'order_1.fits' #There is an order zero...
order1=fits.getdata(filepath)#Ok that was easy.

print('Reading data')
t1=ut.start()
list_of_wls=[]
list_of_orders=[]
for i in range(15):
    list_of_wls.append(fits.getdata(dp+'wave_%s.fits' % i))
    list_of_orders.append(fits.getdata(dp+'order_%s.fits' % i))
t2=ut.end(t1)

#fits.writeto(filename,data,header=h)  #would be the opposite of getdata().

#The following tests inject_model
# print('Injecting model')
# t1=ut.start()
# order_injected=models.inject_model(wl1,order1,dp,modelname)
# t2=ut.end(t1)

#The following tests continuum normalization for template construction.
print('Building template')
t1=ut.start()
wlt,T=models.build_template(templatename,binsize=0.5,maxfrac=0.01,resolution=120000.0)
t2=ut.end(t1)
#CAN'T I OUTSOURCE THIS TO DANIEL AND SIMON? I mean for each spectrum they could also
#generate a binary mask by putting delta's instead of.. well... not precisely
#because of masking by line wings; that's the whole trick of opacity calculation.
#And then Daniel can subtract the continuum himself...

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


#THE ACTUAL CALL TO CROSS_CORRELATION
print('Starting cross-correlation')
t1=ut.start()
rv,ccf=analysis.xcor(list_of_wls,list_of_orders,wlt,T,drv,RVrange,plot=True)
t2=ut.end(t1)
fits.writeto('rv.fits',rv)
fits.writeto('test.fits',ccf,overwrite=True)


#Now comes an analysis of the CCF. Could be that this needs to be moved to
#analysis.py.
rv=fits.getdata('rv.fits')
ccf = fits.getdata('test.fits')
meanflux=np.median(ccf,axis=1)
meanblock=fun.rebinreform(meanflux,len(rv))
ccf_n = ccf/meanblock.T
fits.writeto('ccf_n.fits',ccf_n,overwrite=True)
meanccf=np.mean(ccf,axis=0)
meanblock=fun.rebinreform(meanccf,len(meanflux))
ccf_nn = ccf_n/meanblock
fits.writeto('ccf_nn.fits',ccf_nn,overwrite=True)


pdb.set_trace()

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
