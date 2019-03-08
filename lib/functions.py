def gaussian(x,A,mu,sig):
    """This produces a gaussian function on the grid x with amplitude A, mean mu
    and standard deviation sig. Will need to expand it with a version that has
    a polynomial continuum in the same way that IDL does it."""

    import numpy as np
    return A * np.exp(-0.5*(x - mu)/sig*(x - mu)/sig)

    

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
