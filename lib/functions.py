def sigma_clip(array,nsigma=3.0):
    """This returns the edge values of a sigma-clipping operation.
    Write tests and documentation please."""
    import numpy as np
    m = np.nanmedian(array)
    s = np.nanstd(array)
    vmin = m-nsigma*s
    vmax = m+nsigma*s
    return vmin,vmax

def selmax(y,p,s=0.0):
    """This program returns the p (fraction btw 0 and 1) highest points in y,
    ignoring the very top s % (default zero, i.e. no points ignored), for the
    purpose of outlier rejection."""
    import lib.utils as ut
    import numpy as np
    ut.postest(p)
    if s < 0.0:
        raise Exception("ERROR in selmax: s should be zero or positive.")
    if p >= 1.0:
        raise Exception("ERROR in selmax: p should be strictly between 0.0 and 1.0.")
    if s >= 1.0:
        raise Exception("ERROR in selmax: s should be strictly less than 1.0.")
    ut.postest(-1.0*p+1.0)
    ut.nantest('y in selmax',y)
    ut.dimtest(y,[0])#Test that it is one-dimensional.

    y_sorting = np.flipud(np.argsort(y))#These are the indices in descending order (thats why it gets a flip)
    N=len(y)
    if s == 0.0:
        max_index = np.max([int(round(N*p)),1])#Max because if the fraction is 0 elements, then at least it should contain 1.0
        return y_sorting[0:max_index]

    if s > 0.0:
        min_index = np.max([int(round(N*s)),1])#If 0, then at least it should be 1.0
        max_index = np.max([int(round(N*(p+s))),2]) #If 0, then at least it should be 1+1.
        return y_sorting[min_index:max_index]


def gaussian(x,A,mu,sig,cont=0.0):
    """This produces a gaussian function on the grid x with amplitude A, mean mu
    and standard deviation sig. Will need to expand it with a version that has
    a polynomial continuum in the same way that IDL does it."""

    import numpy as np
    return A * np.exp(-0.5*(x - mu)/sig*(x - mu)/sig)+cont


def box(x,A,c,w):
    """This function computes a box with width w, amplitude A and center c
    on grid x. It would be a simple multiplication of two heaviside functions,
    where it not that I am taking care to interpolate the edge pixels.
    And because I didn't want to do the cheap hack of oversampling x by a factor
    of many and then using np.interpolate, I computed the interpolation manually
    myself. This function guarantees* to preserve the integral of the box, i.e.
    np.sum(box(x,A,c,w))*dwl equals w, unless the box is truncated by the edge.)

    *I hope... Uncomment the two relevant lines below to make sure."""
    import numpy as np
    #import matplotlib.pyplot as plt
    #import pdb
    y=x*0.0

    #The following manually does the interpolation of the edge pixels of the box.
    #First do the right edge of the box.
    r_dist=x-c-0.5*w
    r_px=np.argmin(np.abs(r_dist))
    r_mindist=r_dist[r_px]
    y[0:r_px]=1.0
    if r_px != 0 and r_px != len(x)-1:#This only works if we are not at an edge.
        dpx=abs(x[r_px]-x[r_px-np.int(np.sign(r_mindist))])#Distance to the
        #previous or next pixel, depending on whether the edge falls to the left
        #(posive mindist) or right (negative mindist) of the pixel center.

    #If we are an edge, dpx needs to be approximated from the other neigboring px.
    elif r_px == 0:
        dpx=abs(x[r_px+1]-x[0])
    else:
        dpx=abs(x[r_px]-x[r_px-1])

    if w/dpx < 2.0:
        raise Exception("ERROR IN BOXSMOOTH: WIDTH TOO SMALL (< 2x sampling rate)")
    frac=0.5-r_mindist/dpx
    #If mindist = 0, then frac = 0.5.
    #If mindist = +0.5*dpx, then frac=0.0
    #If mindist = -0.5*dpx, then frac=1.0
    y[r_px] = np.clip(frac,0,1)#Set the pixel that overlaps with the edge to the
    #fractional distance that the edge is away from the px center. Clippig needs
    #to be done to take into account the edge pixels, in which case frac can be
    #larger than 1.0 or smaller than 0.0.

    #And do the same for the left part. Note the swapping of two signs.
    l_dist=x-c+0.5*w
    l_px=np.argmin(np.abs(l_dist))
    l_mindist=l_dist[l_px]
    y[0:l_px]=0.0
    if l_px != 0 and l_px != len(x)-1:#This only works if we are not at an edge.
        dpx=abs(x[l_px]-x[l_px-np.int(np.sign(l_mindist))])
    elif l_px == 0:
        dpx=abs(x[l_px+1]-x[0])
    else:
        dpx=abs(x[l_px]-x[l_px-1])
    frac=0.5+l_mindist/dpx
    y[l_px] = np.clip(frac,0,1)
    #If mindist = 0, then frac = 0.5.
    #If mindist = +0.5*dpx, then frac=1.0
    #If mindist = -0.5*dpx, then frac=0.0

    #print([w,np.sum(y)*dpx]) #<=This line compares the integral.
    #In perfect agreement!!

    #plt.plot(x,y,'.')
    #plt.axvline(x=c+0.5*w)
    #plt.axvline(x=c-0.5*w)
    #plt.ion()
    #plt.show()
    return y*A


def findgen(n):
    """This is basically IDL's findgen function.
    a = findgen(5) will return an array with 5 elements from 0 to 4:
    [0,1,2,3,4]
    """
    import numpy as np
    return np.linspace(0,n-1,n)


def rebinreform(a,n):
    """This program works like the rebin(reform()) trick in IDL, where you use fast
    array manipulation operations to transform a 1D array into a 2D stack of itself,
    to be able to do operations on another 2D array by multiplication/addition/division
    without having to loop through the second dimension of said array."""
    import numpy as np
    return(np.transpose(np.repeat(np.expand_dims(a,1),n,axis=1)))


def read_binary(path,double=False):
    #This code comes through Lorenzo, to read binary files.
    import numpy as np
    import struct
# Define format of the single chunk, in this case 1 float object
    struct_fmt = '=f'
    if double == True:
        struct_fmt = '=d'
# struct_fmt = '@f'  # Same result, not sure which is the difference
# Length of the single chunk, 4 bytes for us
    struct_len = struct.calcsize(struct_fmt)
    #pdb.set_trace()
# Build the interpreter of the binary format
    struct_unpack = struct.Struct(struct_fmt).unpack_from
    opacities = []
# Open with the b option to read binary
    # i=0
    with open(path, "rb") as f:
    # Until eof
        while True:
            data = f.read(struct_len)
            # if i > 1.1e7:
                # print(i)
                # print(len(data))
            if not data: break
            opacities.append(struct_unpack(data))
            # i+=1
    return(np.array(opacities))


def doppler_shift(wl_source,dv):
    """This function returns the rel. doppler shifted wavelength of a velocity
    dv applied to source wavelength wl_source."""
    import lib.constants as const
    import numpy as np
    b = (dv*1000.0)/const.c
    return(np.sqrt((1+b)/(1-b))*wl_source)

def local_v_star(phase,aRstar,inclination,vsini,l):
    """This is the rigid-body, circular-orbt approximation of the local velocity occulted
    by the planet as it goes through transit, as per Cegla et al. 2016"""
    import numpy as np
    xp = aRstar * np.sin(2.0*np.pi*phase)
    yp = (-1.0)*aRstar * np.cos(2.0*np.pi*phase) * np.cos(np.deg2rad(inclination))
    x_per = xp*np.cos(np.deg2rad(l)) - yp*np.sin(np.deg2rad(l))
    return(x_per*vsini)

def running_MAD_2D(z,w):
    """Computers a running standard deviation of a 2-dimensional array z.
    The stddev is evaluated over the vertical block with width w pixels.
    The output is a 1D array with length equal to the width of z."""
    import astropy.stats as stats
    import numpy as np
    size = np.shape(z)
    ny = size[0]
    nx = size[1]
    import pdb
    s = findgen(nx)*0.0
    for i in range(nx):
        minx = max([0,i-int(0.5*w)])
        maxx = min([nx-1,i+int(0.5*w)])
        s[i] = stats.mad_std(z[:,minx:maxx],ignore_nan=True)

    return(s)
