import lib.utils as ut
import numpy as np
import lib.operations as ops
from lib import functions as fun
from lib import system_parameters as sp
from lib import models
from astropy.io import fits
import pylab
import pdb


dp='data/WASP-19/night1/'
#v=sp.transit(dp)

filepath=dp+'order_0.fits'

b=fits.getdata(filepath)#Ok that was easy.
wl = fits.getdata(dp+'wave_0.fits')
#fits.writeto(filename,data,header=h)  #would be the opposite of getdata().

#The following tests inject_model

#c=models.inject_model(wl,b,dp,'W121-plez')





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
spec = fun.findgen(len(wl))*0.0
spec[500] = 1
spec[1000] = 1
spec[3000] = 1
spec[3012] = 1
#spec_b=ops.blur_rotate(wl,spec,3.0,1.5,1.5,90.0)
#spec_b2=ops.blur_rotate(wl,spec,3.0,2.0,1.0,90.0)
t1=ut.start()
spec_b=ops.blur_spec(wl,spec,20.0,mode='box')
spec_b2=ops.blur_spec(wl,spec,20.0,mode='box')
dt1=ut.end(t1)
pylab.plot(wl,spec)
pylab.plot(wl,spec_b)
pylab.plot(wl,spec_b2)
pylab.show()


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
