#This package contains operations that act on spectra, such as blurring with
#various kernels and continuum normalization.

def envelope(wlm,fxm,binsize,selfrac=0.05,mode='top',threshold=''):
    """This program measures the top or bottom envelope of a spectrum (wl,fx), by
    chopping it up into bins of size binsze (unit of wl), and measuring the mean
    of the top n % of values in that bin. Setting the mode to 'bottom' will do the
    oppiste: The mean of the bottom n% of values. The output is the resulting wl
    and flux points of these bins.

    Example: wle,fxe = envelope(wl,fx,1.0,selfrac=3.0,mode='top')"""
    import pdb
    import numpy as np
    import lib.utils as ut
    import lib.functions as fun
    from matplotlib import pyplot as plt
    ut.typetest('wlm in envelope',wlm,np.ndarray)
    ut.typetest('fxm in envelope',fxm,np.ndarray)
    ut.dimtest(wlm,[len(fxm)])
    ut.typetest('binsize in envelope',binsize,float)
    ut.typetest('percentage in envelope',selfrac,float)
    ut.typetest('mode in envelope',mode,str)
    ut.nantest('fxm in envelope',fxm)
    ut.nantest('wlm in envelope',wlm)
    ut.postest(wlm,varname='wlm in envelope')
    ut.postest(binsize,varname='binsize in envelope')

    if mode == 'bottom':
        fxm*=-1.0
    # wlcs=np.array([])#Center wavelengths
    # fxcs=np.array([])#Flux at center wavelengths
    # wlb = 0.5*binsize+np.min(wlm)
    # t1=ut.start()
    # while wlb+0.5*binsize < np.max(wlm):
    #     #THIS ARRAY SELECTION IS INTENSELY SLOW. NEEDS TO BE FIXED! Stackoverflow says
    #     #that a forloop could be very fast. But that sucks forloops in while-loops...
    #     #sel = np.where( np.logical_and(wlm >= (wlb-0.5*binsize),wlm < (wlb+0.5*binsize)))[0]
    #     wlsel = wlm[(wlm >= (wlb-0.5*binsize)) & (wlm < (wlb+0.5*binsize))]
    #     fxsel = fxm[(wlm >= (wlb-0.5*binsize)) & (wlm < (wlb+0.5*binsize))]
    #     maxsel = fun.selmax(fxsel,selfrac)
    #     wlcs=np.append(wlcs,np.mean(wlsel[maxsel]))
    #     fxcs=np.append(fxcs,np.mean(fxsel[maxsel]))
    #     wlb+=binsize
    # t2=ut.end(t1)
    # plt.plot(wlm,fxm)
    # plt.plot(wlcs,fxcs,'.')
    #
    #

    #The part that is commented out above is the IDL-savvy way to do binning, with
    #the where function. Unfortunately np.where is super slow on large arrays,
    #so I did it in a move naive way below: Start counting wlm at its first element,
    #and count upwards, until wlm[i]-wlm[i_start] is (larger than) the binsize.
    #Then select everything in wlm from wlm[i_start] to wlm[i], calculate the avg
    #maximum flux there and reset i_start to continue to the next bin.
    #Looping over wlm once is more efficient than looping over bins and calling
    #np.where each time, as long as the number of bins is large. In most cases all
    #of this takes (much) less than a second, down from 30-ish seonds in the
    #strategy that uses np.where(). Days like these make me miss IDL...

    wlcs=np.array([])#Center wavelengths
    fxcs=np.array([])#Flux at center wavelengths
    # t3=ut.start()
    i_start=0
    wlm_start=wlm[i_start]
    for i in range(0,len(wlm)):
        if wlm[i]-wlm_start >= binsize:
            wlsel = wlm[i_start:i]
            fxsel = fxm[i_start:i]
            maxsel = fun.selmax(fxsel,selfrac)
            wlcs=np.append(wlcs,np.mean(wlsel[maxsel]))
            fxcs=np.append(fxcs,np.mean(fxsel[maxsel]))
            i_start=i+1
            wlm_start=wlm[i+1]

    if isinstance(threshold,float) == True:
        #This means that the threshold value is set, and we set all bins less than
        #that threshold to the threshold value:
        if mode == 'bottom':
            threshold*=-1.0
        fxcs[(fxcs < threshold)] = threshold


    # t4=ut.end(t3)
    # plt.plot(wlcs,fxcs,'.')
    # plt.show()
    #
    #
    # pdb.set_trace()
    if mode == 'bottom':
        fxcs*=-1.0
        fxm*=-1.0
    return wlcs,fxcs

def normalize(wlm,fxm,binsize,emission=False,mode='linear'):
    import numpy as np
    import pdb
    import lib.functions as fun
    import lib.utils as ut


def convolve(array,kernel,edge_degree=1):
    """It's unbelievable, but I could not find the python equivalent of IDL's
    /edge_truncate keyword, which truncates the kernel at the edge of the convolution.
    Therefore, unfortunately, I need to code a convolution operation myself.
    Stand by to be slowed down by an order of magnitude #thankspython.

    Nope! Because I can just use np.convolve for most of the array, just not the edge...

    So the strategy is to extrapolate the edge of the array using a polynomial fit
    to the edge elements. I fit over a range that is twice the length of the kernel.

    Example: y_blurred = convolve(x,y,edge_degree = 2)
    For a convolution where the edge is extrapolated with a second degree polynomial.
    """
    import numpy as np
    import pdb
    import lib.functions as fun
    import lib.utils as ut
    ut.typetest('edge_degree',edge_degree,int)
    ut.typetest('array',array,np.ndarray)
    ut.typetest('kernel',kernel,np.ndarray)

    if len(kernel) >= len(array)/2:
        raise Exception("Error in convolution: Kernel length is larger than half of the array. Can't extrapolate over that length. And you probably don't want to be doing a convolution like that, anyway.")

    if len(kernel) % 2 != 1:
        raise Exception('Error in convolution: Kernel needs to have an odd number of elements.')

    #Perform polynomial fits at the edges.
    x=fun.findgen(len(array))
    fit_left=np.polyfit(x[0:len(kernel)*2],array[0:len(kernel)*2],edge_degree)
    fit_right=np.polyfit(x[-2*len(kernel)-1:-1],array[-2*len(kernel)-1:-1],edge_degree)

    #Pad both the x-grid (onto which the polynomial is defined)
    #and the data array.
    pad=fun.findgen( (len(kernel)-1)/2)
    left_pad=pad-(len(kernel)-1)/2
    right_pad=np.max(x)+pad+1
    left_array_pad=np.polyval(fit_left,left_pad)
    right_array_pad=np.polyval(fit_right,right_pad)

    #Perform the padding.
    x_padded = np.append(left_pad , x)
    x_padded = np.append(x_padded , right_pad) #Pad the array with the missing elements of the kernel at the edge.
    array_padded = np.append(left_array_pad,array)
    array_padded = np.append(array_padded,right_array_pad)

    #Reverse the kernel because np.convol does that automatically and I don't want that.
    #(Imagine doing a derivative with a kernel [-1,0,1] and it gets reversed...)
    kr = kernel[::-1]
    #The valid keyword effectively undoes the padding, leaving only those values for which the kernel was entirely in the padded array.
    #This thus again has length equal to len(array).
    return np.convolve(array_padded,kr,'valid')


def derivative(x):
    import numpy as np
    d_kernel=np.array([-1,0,1])/2.0
    return(convolve(x,d_kernel))


def smooth(fx,w,mode='box',edge_degree=1):
    """This function takes a spectrum, and blurs it using either a
    Gaussian kernel or a box kernel, which have a FWHM width of w px everywhere.
    Meaning that the width changes dynamically on a constant d-lambda grid.

    Set the mode to gaussian or box. Because in box, care is taken to correctly
    interpolate the edges, it is about twice slower than the Gaussian.
    This interpolation is done manually in the fun.box function."""

    import numpy as np
    import lib.utils as ut
    import lib.functions as fun
    import pdb
    from matplotlib import pyplot as plt
    import lib.constants as const
    import time
    ut.typetest('w',w,float)
    ut.typetest('fx',fx,np.ndarray)
    ut.typetest('mode',mode,str)
    ut.typetest('edge_degree',edge_degree,int)

    truncsize=8.0#The gaussian is truncated at 8 sigma.
    shape=np.shape(fx)

    sig_w = w / 2*np.sqrt(2.0*np.log(2)) #Transform FWHM to Gaussian sigma. In km/s.
    trunc_dist=np.round(sig_w*truncsize).astype(int)

    #First define the kernel.
    kw=int(np.round(truncsize*sig_w*2.0))
    if kw % 2.0 != 1.0:#This is to make sure that the kernel has an odd number of
    #elements, and that it is symmetric around zero.
        kw+=1

    kx=fun.findgen(kw)
    kx-=np.mean(kx)#This must be centered around zero. Doing a hardcoded check:
    if (-1.0)*kx[-1] != kx[0]:
        print(kx)
        raise Exception("ERROR in box_smooth: Kernel could not be made symmetric somehow. Attempted kernel grid is printed above. Kernel width is %s pixels." % kw)


    if mode == 'gaussian':
        k=fun.gaussian(kx,1.0,0.0,sig_w)


    if mode == 'box':
        k=fun.box(kx,1.0,0.0,w)
        kx=kx[k > 0.0]
        k=k[k > 0.0]
        if (-1.0)*kx[-1] != kx[0]:
            print(kx)
            raise Exception("ERROR in box_smooth: Kernel could not be made symmetric AFTER CROPPING OUT THE BOX, somehow. Attempted kernel grid is printed above. Kernel width is %s pixels." % kw)

    k/=np.sum(k)

    return(convolve(fx,k,edge_degree))



    if len(shape) == 1:
        #print('Do the entire thing in 1D')
        order_blurred=order*0.0
        if mode == 'gaussian':
            for i in range(0,len(wl)):
                #Im going to select wl in a bin so that I dont need to evaluate a gaussian over millions of points that are all zero
                binstart=max([0,i-trunc_dist[i]])
                binend=i+trunc_dist[i]
                k = fun.gaussian(wl[binstart:binend],1.0,wl[i],sig_wl[i])
                k_n=k/np.sum(k)
                order_blurred[i]=np.sum(k_n*order[binstart:binend])
                #To speed up, need to select wl and then append with zeroes. <= what does that mean? Jens 03 mar 18
            return(order_blurred)
        elif mode == 'box':
            for i in range(0,len(wl)):
                binstart=max([0,i-trunc_dist[i]])
                binend=i+trunc_dist[i]
                k = fun.box(wl[binstart:binend],1.0,wl[i],dwl[i])
                k_n=k/np.sum(k)
                order_blurred[i]=np.sum(k_n*order[binstart:binend])
            return(order_blurred)
        else:
            raise Exception("ERROR: Mode should be set to 'gaussian' or 'box'.")


def constant_velocity_wl_grid(wl,fx,oversampling=1.0):
    """This function will define a constant-velocity grid that is (optionally)
    sampled a number of times finer than the SMALLEST velocity difference that is
    currently in the grid.

    Example: wl_cv,fx_cv = constant_velocity_wl_grid(wl,fx,oversampling=1.5).

    WARNING: This function is hardcoded to raise an exception if wl or fx contain NaNs,
    because interp1d does not handle NaNs."""
    import lib.constants as consts
    import numpy as np
    import lib.functions as fun
    import lib.utils as ut
    from scipy import interpolate
    import pdb
    import matplotlib.pyplot as plt
    ut.typetest('oversampling',oversampling,float)
    ut.typetest('wl',wl,np.ndarray)
    ut.typetest('fx',fx,np.ndarray)
    ut.nantest('wl',wl)
    ut.nantest('fx',fx)


    if oversampling <= 0.0:
        raise Exception("ERROR in constant velocity wl grid: oversampling should be positive and finite.")

    c=consts.c/1000.0
    dl=derivative(wl)
    dv=dl/wl*c
    a=np.min(dv)/oversampling

    wl_new=0.0
    #The following while loop will define the new pixel grid.
    #It starts trying 100,000 points, and if that's not enough to cover the entire
    #range from min(wl) to max(wl), it will add 100,000 more; until it's enough.
    n=len(wl)
    while np.max(wl_new) < np.max(wl):
        x=fun.findgen(n)
        wl_new=np.exp(a/c * x)*np.min(wl)
        n+=len(wl)
    wl_new[0]=np.min(wl)#Artificially set to zero to avoid making a small round
    #off error in that exponent.

    #Then at the end we crop the part that goes too far:
    wl_new_cropped=wl_new[(wl_new <= np.max(wl))]
    x_cropped=x[(wl_new <= np.max(wl))]
    i_fx = interpolate.interp1d(wl,fx)
    fx_new_cropped =i_fx(wl_new_cropped)
    return(wl_new_cropped,fx_new_cropped,a)


def blur_rotate(wl,order,dv,Rp,P,inclination):
    """This function takes a spectrum and blurs it using a rotation x Gaussian
    kernel which has a FWHM width of dv km/s everywhere. Meaning that its width changes
    dynamically.
    Because the kernel needs to be recomputed on each element of the wavelength axis
    individually, this operation is (or should be) much slower than convolution with
    a constant kernel, in which a simple shifting of the array, rather than a recomputation
    of the rotation profile is sufficient.

    Input:
    The wavelength axis wl.
    The spectral axis order.
    The FHWM width of the resolution element in km/s.
    The Radius of the rigid body in Rj.
    The periodicity of the rigid body rotation in days.
    The inclination of the spin axis in degrees.

    Wavelength and order need to be numpy arrays and have the same number of elements.
    Rp, P and i need to be scalar floats.

    Output:
    The blurred spectral axis, with the same dimensions as wl and order.


    WARNING: THIS FUNCTION HANDLES NANS POORLY. I HAVE THEREFORE DECIDED CURRENTLY
    TO REQUIRE NON-NAN INPUT."""
    import numpy as np
    import lib.utils as ut
    import lib.functions as fun
    import pdb
    from matplotlib import pyplot as plt
    import lib.constants as const
    import time
    from scipy import interpolate
    ut.typetest('dv',dv,float)
    ut.typetest('wl',wl,np.ndarray)
    ut.typetest('order',order,np.ndarray)
    ut.typetest('P',P,float)
    ut.typetest('Rp',Rp,float)
    ut.typetest('i',inclination,float)
    ut.nantest('wl',wl)
    ut.nantest('order',order)
    ut.dimtest(wl,[0])
    ut.dimtest(order,[len(wl)])#Test that wl and order are 1D, and that
    #they have the same length.

    if np.min(np.array([dv,P,Rp])) <= 0.0:
        raise Exception("ERROR in blur_rotate: dv, P and Rp should be strictly positive.")

    #ut.typetest_array('wl',wl,np.float64)
    #ut.typetest_array('order',order,np.float64)
    #This is not possible because order may be 2D...
    #And besides, you can have floats, np.float32 and np.float64... All of these would
    #need to pass. Need to fix typetest_array some day.


    order_blurred=order*0.0#init the output.
    truncsize=5.0#The gaussian is truncated at 5 sigma from the extremest points of the RV amplitude.
    sig_dv = dv / (2*np.sqrt(2.0*np.log(2))) #Transform FWHM to Gaussian sigma. In km/s.
    deriv = derivative(wl)
    if max(deriv) < 0:
        raise Exception("ERROR in ops.blur_rotate: WL derivative is smaller than 1.0. Sort wl in ascending order.")
    sig_wl=wl*sig_dv/const.c*1000.0#in nm
    sig_px=sig_wl/deriv

    n=1000.0
    a=fun.findgen(n)/(n-1)*np.pi
    rv=np.cos(a)*2.0*np.pi*Rp*const.Rjup/P/const.day*np.sin(np.radians(inclination))/1000.0
    trunc_dist=np.round(sig_px*truncsize+np.max(rv)*wl/const.c*1000.0/deriv).astype(int)



    rvgrid_max=(np.max(trunc_dist)+1.0)*sig_dv+np.max(rv)
    rvgrid_n=rvgrid_max / dv * 100.0 #100 samples per lsf fwhm.
    rvgrid=(fun.findgen(2*rvgrid_n+1)-rvgrid_n)/rvgrid_n*rvgrid_max#Need to make sure that this is wider than the truncation bin and more finely sampled than wl - everywhere.

    lsf=rvgrid*0.0
    #We loop through velocities in the velocity grid to build up the sum of Gaussians
    #that is the LSF.
    for v in rv:
        lsf+=fun.gaussian(rvgrid,1.0,v,sig_dv)#This defines the LSF on a velocity grid wih high fidelity.

    #Now we loop through the wavelength grid to place this LSF at each wavelength position.
    for i in range(0,len(wl)):
        binstart=max([0,i-trunc_dist[i]])
        binend=i+trunc_dist[i]
        wlbin=wl[binstart:binend]

        wlgrid =   wl[i]*rvgrid/(const.c/1000.0)+wl[i]#This converts the velocity grid to a d-wavelength grid centered on wk[i]
        #print([np.min(wlbin),np.min(wlgrid),np.max(wlbin),np.max(wlgrid)])

        i_wl = interpolate.interp1d(wlgrid,lsf) #This is a class that can be called.
        lsf_wl=i_wl(wlbin)
        k_n=lsf_wl/np.sum(lsf_wl)#Normalize at each instance of the interpolation to make sure flux is conserved exactly.
        order_blurred[i]=np.sum(k_n*order[binstart:binend])

    return order_blurred


def blur_spec(wl,order,dv,mode='gaussian'):
    """This function takes a spectrum, and blurs it using either a
    Gaussian kernel or a box kernel, which have a FWHM width of dv km/s everywhere.
    Meaning that the width changes dynamically on a constant d-lambda grid.
    Because the kernel needs to be recomputed on each element of the wavelength axis
    individually, this operation is much slower than convolution with
    a constant kernel, in which a simple shifting of the array, rather than a recomputation
    of the kernel is sufficient.

    This program is repeated above for a rotation kernel as well.

    Set the mode to gaussian or box. Because in box, care is taken to correctly
    interpolate the edges, it is about twice slower than the Gaussian.
    This interpolation is done manually in the fun.box function."""

    print('I MAY NOT WANT TO USE BLUR-SPEC BECAUSE IT IS SLOW, AT LEAST IN BOX MODE.')
    print('AND I HAVE NOT THOROUGHLY BENCHMARKED IT.')
    import numpy as np
    import lib.utils as ut
    import lib.functions as fun
    import pdb
    from matplotlib import pyplot as plt
    import lib.constants as const
    import time
    ut.typetest('dv',dv,float)
    ut.typetest('wl',wl,np.ndarray)
    ut.typetest('order',order,np.ndarray)
    ut.typetest('mode',mode,str)


    truncsize=8.0#The gaussian is truncated at 5 sigma.
    shape=np.shape(order)

    sig_dv = dv / 2*np.sqrt(2.0*np.log(2)) #Transform FWHM to Gaussian sigma. In km/s.

    d_kernel=np.array([-1,0,1])/2.0
    deriv = convolve(wl,d_kernel)
    #l*dv/c=dl
    dwl=wl*dv/const.c*1000.0
    sig_wl=wl*sig_dv/const.c*1000.0#in nm
    sig_px=sig_wl/deriv
    trunc_dist=np.round(sig_px*truncsize).astype(int)

    if len(shape) == 1:
        #print('Do the entire thing in 1D')
        order_blurred=order*0.0
        if mode == 'gaussian':
            for i in range(0,len(wl)):
                #Im going to select wl in a bin so that I dont need to evaluate a gaussian over millions of points that are all zero
                binstart=max([0,i-trunc_dist[i]])
                binend=i+trunc_dist[i]
                k = fun.gaussian(wl[binstart:binend],1.0,wl[i],sig_wl[i])
                k_n=k/np.sum(k)
                order_blurred[i]=np.sum(k_n*order[binstart:binend])
                #To speed up, need to select wl and then append with zeroes. <= what does that mean? Jens 03 mar 18
            return(order_blurred)
        elif mode == 'box':
            for i in range(0,len(wl)):
                binstart=max([0,i-trunc_dist[i]])
                binend=i+trunc_dist[i]
                k = fun.box(wl[binstart:binend],1.0,wl[i],dwl[i])
                k_n=k/np.sum(k)
                order_blurred[i]=np.sum(k_n*order[binstart:binend])
            return(order_blurred)
        else:
            raise Exception("ERROR: Mode should be set to 'gaussian' or 'box'.")


    if len(shape) == 2:
        print('Do the entire thing in 2D')
        raise Exception("ERROR: THIS FUNCTIONALITY HAS NOT YET BEEN ADDED.")

    if len(shape) >=3:
        raise Exception("Error in blur_gaussian_lsf: dimension of order should be 1 or 2, but is %s" % len(shape))


def airtovac(wlnm):
    wlA=wlnm*10.0
    s = 1e4 / wlA
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    return(wlA*n/10.0)

def normalize_orders(list_of_orders):
    import numpy as np
    import lib.functions as fun
    import pdb
    N = len(list_of_orders)
    out_list_of_orders=[]
    n_px=np.shape(list_of_orders[0])[1]

    for i in range(N):
        meanflux=np.nanmean(list_of_orders[i],axis=1)#Average flux in each order.
        meanblock=fun.rebinreform(meanflux/np.nanmean(meanflux),n_px).T
        out_list_of_orders.append(list_of_orders[i]/meanblock)
    return(out_list_of_orders)




def shift_ccf(RV,CCF,drv):
        import lib.functions as fun
        import lib.utils as ut
        import numpy as np
        #import matplotlib.pyplot as plt
        import scipy.interpolate
        import pdb
        import astropy.io.fits as fits
        import matplotlib.pyplot as plt
        import sys
        import scipy.ndimage.interpolation as scint


        if np.ndim(CCF) == 1.0:
            print("ERROR in shift_ccf: CCF should be a 2D block.")
            sys.exit()
        else:
            n_exp=len(CCF[:,0])#Number of exposures.
            n_rv=len(CCF[0,:])
        if len(RV) != n_rv:
            print('ERROR in shift_ccf: RV does not have the same length as the base size of the CCF block.')
            sys.exit()
        if len(drv) != n_exp:
            print('ERROR in shift_ccf: drv does not have the same height as the CCF block.')
            sys.exit()
        dv = RV[1]-RV[0]
        CCF_new=CCF*0.0
        for i in range(n_exp):
            #C_i=scipy.interpolate.interp1d(RV,CCF[i],fill_value=(0.0,0.0))
            #CCF_new[i,:] = C_i(RV-drv[i]*2.0)
            CCF_new[i,:] = scint.shift(CCF[i],drv[i]/dv,mode='nearest',order=1)
        return(CCF_new)

def gauss_fit(x,y,start=None,plot=False):
    """This simply fits a gaussian with an offset. Can be used generally, but
    was created for measure_RV below."""
    import matplotlib.pyplot as plt
    import numpy as np
    from lmfit import Model
    import pdb

    def gaussian(x, amp, cen, wid, cont):
        """1-d gaussian: gaussian(x, amp, cen, wid)"""
        return (amp * np.exp(-(x-cen)**2 / (2*wid**2)) + cont)


    gmodel = Model(gaussian)

    if start == None:
        maxindex = (abs(y-np.nanmedian(y))).argmax()

        start = [0.0,x[maxindex],(max(x)-min(x))/10.0,np.nanmedian(y)]
        print(start)
    result = gmodel.fit(y, x=x, amp=start[0], cen=start[1], wid=start[2], cont=start[3])

    # print(result.fit_report())

    aaa=result.best_values
    fitamp = aaa['amp']
    fitcen = aaa['cen']
    fitwid = aaa['wid']
    fitcon = aaa['cont']

    if plot == True:
        plt.plot(x, y, 'bo')
        plt.plot(x, result.init_fit, 'k--')
        plt.plot(x, result.best_fit, 'r-')
        plt.show()
    return(fitamp,fitcen,fitwid,fitcon)



def measure_rv(RV,CCF1D,dv=15.0,plot=False):
    """This fits a Gaussian to a 1D CCF."""
    import numpy as np
    import sys
    import lib.functions as fun
    import lib.operations as ops
    import matplotlib.pyplot as plt
    import pdb

    maxindex = (abs(CCF1D-np.nanmedian(CCF1D))).argmax()
    RV_0 = RV[maxindex]
    sel = [(RV > RV_0-dv) & (RV < RV_0+dv)]
    y = CCF1D
    fit=ops.gauss_fit(RV[sel],y[sel],start=[0.0,RV_0,dv,np.nanmedian(y)])
    if plot == True:
        plt.plot(RV,y,'.')
        plt.plot(RV[sel],y[sel],'.')
        plt.plot(RV,fun.gaussian(RV,fit[0],fit[1],fit[2],cont=fit[3]))
        plt.show()
    return(fit[1])
