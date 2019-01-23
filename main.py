import numpy as np
from lib import system_parameters as sp
from astropy.io import fits
from pylab import *

dp='data/WASP-19/night1/'
#v=sp.transit(dp)


filepath=dp+'order_0.fits'

b=fits.getdata(filepath)#Ok that was easy.
#fits.writeto(filename,data,header=h)  #would be the opposite of getdata().


#This is Heathers analytical formula hurray,
phase=linspace(-0.05,0.05,100)
aRs=5.0
i=radians(86.0)
l1=radians(0.0)
l2=radians(60.0)
vsini=100.0

xp=aRs*np.sin(2.0*np.pi*phase)
yp=-1.0*aRs*np.cos(2.0*np.pi*phase)*np.cos(i)

vs1=vsini*(xp*np.cos(l1)-yp*np.sin(l1))
vs2=vsini*(xp*np.cos(l2)-yp*np.sin(l2))
plot(phase,vs1)
plot(phase,vs2)
show()
