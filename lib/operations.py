#This package contains operations that act on spectra, such as blurring with
#various kernels and continuum normalization.

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
    The blurred spectral axis, with the same dimensions as wl and order."""

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
    ut.dimtest(wl,[0])
    ut.dimtest(order,[len(wl)])#Test that wl and order are 1D, and that
    #they have the same length.


    #ut.typetest_array('wl',wl,np.float64)
    #ut.typetest_array('order',order,np.float64)
    #This is not possible because order may be 2D...
    #And besides, you can have floats, np.float32 and np.float64... All of these would
    #need to pass. Need to fix typetest_array.
    order_blurred=order*0.0#init the output.
    truncsize=5.0#The gaussian is truncated at 5 sigma from the extremest points of the RV amplitude.
    sig_dv = dv / (2*np.sqrt(2.0*np.log(2))) #Transform FWHM to Gaussian sigma. In km/s.
    d_kernel=np.array([-1,0,1])/2.0
    deriv = convolve(wl,d_kernel)

    sig_wl=wl*sig_dv/const.c*1000.0#in nm
    sig_px=sig_wl/deriv

    n=1000.0
    a=fun.findgen(n)/(n-1)*np.pi
    rv=np.cos(a)*2.0*np.pi*Rp*const.Rjup/P/const.day*np.sin(np.radians(inclination))/1000.0
    trunc_dist=np.round(sig_px*truncsize+np.max(rv)*wl/const.c*1000.0/deriv).astype(int)

    rvgrid_max=(truncsize+1.0)*sig_dv+np.max(rv)
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
        print('Do the entire thing in 1D')
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
