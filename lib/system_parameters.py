#This replaces calctimes and includes paramget.

def paramget(keyword,dp):
    from lib.utils import typetest
    typetest('dp',dp,str)
    typetest('keyword',keyword,str)

    try:
        f = open(dp+'config', 'r')
    except FileNotFoundError:
        raise Exception('Configfile does not exist at %s' % dp) from None
    x = f.read().splitlines()
    f.close()
    n_lines=len(x)
    keywords={}
    for i in range(0,n_lines):
        line=x[i].split()
        try:
            value=float(line[1])
        except ValueError:
            value=(line[1])
        keywords[line[0]] = value
    try:
        return(keywords[keyword])
    except KeyError:
        raise Exception('Keyword %s is not present in configfile at %s' % (keyword,dp)) from None



def v_orb(dp):
    from lib.utils import typetest
    """This program calculates the orbital velocity in km/s for the planet in the
    data sequence provided in dp, the data-path. dp starts in the root folder,
    i.e. it starts with data/projectname/, and it ends with a slash.

    Example: v=v_orb(data/Kelt-9/night1/)
    The output is a number in km/s."""
    typetest('dp',dp,str)

    import numpy as np
    import pdb
    import lib.constants as const


    P=paramget('P',dp)*const.day
    r=paramget('a',dp)*const.AU

    typetest('P',P,float)
    typetest('r',r,float)

    if P <= 0:
        raise Exception("P <= zero. Check configfile at %s." % dp)
    if r <= 0:
        raise Exception("r <= zero. Check configfile at %s." % dp)
    return 2.0*np.pi*r/P/1000.0

def berv(dp):

    """This program retrieves the BERV corrcetion tabulated in the obs_times file.
    Example: brv=berv('data/Kelt-9/night1/')
    The output is an array with length N, corresponding to N exposures. These values
    are / should be taken from the FITS header.
    """
    from lib.utils import typetest
    import numpy as np
    import pdb
    from astropy.io import ascii
    from astropy.time import Time
    from astropy import units as u, coordinates as coord
    typetest('dp',dp,str)
    d=ascii.read(dp+'obs_times',comment="#")#,names=['mjd','time','exptime','airmass'])
    #Removed the named columns because I may not know for sure how many columns
    #there are, and read-ascii breaks if only some columns are named.
    #The second column has to be a date array though.
    berv = d['col5']
    return berv.data

def phase(dp):

    """This program calculates the orbital phase of the planet in the data
    sequence provided using the parameters in dp/config and the timings in
    dp/obstimes. dp starts in the root folder, i.e. it starts with
    data/projectname/, and it ends with a slash.

    Example: p=phase('data/Kelt-9/night1/')
    The output is an array with length N, corresponding to N exposures.

    Be CAREFUL: This program provides a time difference of ~1 minute compared
    to IDL/calctimes. This likely has to do with the difference between HJD
    and BJD, and the TDB timescale. In the future you should have a thorough
    look at the time-issue, because you should be able to get this right to the
    second.

    More importantly: The transit center time needs to be provided in config
    in BJD."""
    from lib.utils import typetest
    import numpy as np
    import pdb
    from astropy.io import ascii
    from astropy.time import Time
    from astropy import units as u, coordinates as coord
    typetest('dp',dp,str)
    d=ascii.read(dp+'obs_times',comment="#")#,names=['mjd','time','exptime','airmass'])
    #Removed the named columns because I may not know for sure how many columns
    #there are, and read-ascii breaks if only some columns are named.
    #The second column has to be a date array though.
    t = Time(d['col2'],scale='utc', location=coord.EarthLocation.of_site('paranal'))
    jd = t.jd
    P=paramget('P',dp)
    RA=paramget('RA',dp)
    DEC=paramget('DEC',dp)
    Tc=paramget('Tc',dp)#Needs to be given in BJD!

    typetest('P',P,float)
    typetest('Tc',Tc,float)

    ip_peg = coord.SkyCoord(RA,DEC,unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = t.light_travel_time(ip_peg)

    n=0.0
    Tc_n=Time(Tc,format='jd',scale='tdb')
    while Tc_n.jd >= min(jd):
        Tc_n=Time(Tc-100.0*n*P,format='jd',scale='tdb')#This is to make sure that the Transit central time PRECEDES the observations (by tens or hundreds or thousands of years). Otherwise, the phase could pick up a minus sign somewhere and be flipped. I hate that.
        n+=1
    BJD = t.tdb + ltt_bary
    diff = BJD-Tc_n
    phase=((diff.jd) % P)/P
    return phase

def RV_star(dp):
    """This program calculates the radial velocity in km/s for the STAR in the
    data sequence provided in dp, the data-path. dp starts in the root folder,
    i.e. it starts with data/projectname/, and it ends with a slash.

    Example: v=RV_star('data/Kelt-9/night1/')
    The output is an array with length N, corresponding to N exposures.
    The radial velocity is provided in km/s. This is meant to be used to correct
    (align) the stellar spectra to the same reference frame. It requires K (the
    RV-semi amplitude to be provided in the config file, in km/s as well. Often
    this value is given in discovery papers.
    Like everything, this code assumes a circular orbit."""
    from lib.utils import typetest
    import numpy as np
    p=phase(dp)
    K=paramget('K',dp)
    typetest('K',K,float)
    rv=K*np.sin(2.0*np.pi*p) * (-1.0)
    return(rv)

def RV(dp):
    """This program calculates the radial velocity in km/s for the planet in the
    data sequence provided in dp, the data-path. dp starts in the root folder,
    i.e. it starts with data/projectname/, and it ends with a slash.

    Example: v=RV('data/Kelt-9/night1/')
    The output is an array with length N, corresponding to N exposures.
    The radial velocity is provided in km/s."""
    from lib.utils import typetest
    import numpy as np
    typetest('dp',dp,str)
    p=phase(dp)
    i=paramget('inclination',dp)
    typetest('i',i,float)
    vorb=v_orb(dp)
    rv=vorb*np.sin(2.0*np.pi*p)*np.sin(np.radians(i))
    return rv#In km/s.


def dRV(dp):
    """This program calculates the change in radial velocity in km/s for the
    planet in the data sequence provided in dp, the data-path. dp starts in the
    root folder,i.e. it starts with data/projectname/, and it ends with a slash.

    Example: dv=dRV('data/Kelt-9/night1/')
    The output is an array with length N, corresponding to N exposures.
    The change in radial velocity is calculated using the first derivative of the
    formula for RV, multiplied by the exposure time provided in obs_times.
    The answer is provided in units of km/s change within each exposure."""
    from lib.utils import typetest
    import numpy as np
    import lib.constants as const
    from astropy.io import ascii
    typetest('dp',dp,str)

    d=ascii.read(dp+'obs_times',names=['mjd','time','exptime','airmass'])
    Texp=d['exptime'].astype('float')
    vorb=v_orb(dp)
    p=phase(dp)
    P=paramget('P',dp)
    i=paramget('inclination',dp)
    typetest('P',P,float)
    typetest('i',i,float)

    dRV=vorb*np.cos(2.0*np.pi*p)*2.0*np.pi/(P*const.day)*np.sin(np.radians(i))
    return abs(dRV*Texp)

def transit(dp):
    """This code uses Ians astro python routines for the approximate Mandel &
    Agol transit lightcurve to produce the predicted transit lightcurve for the
    planet described by the configfile located at dp/config.
    This all assumes a circular orbit.
    ===========
    Derivation:
    ===========
    occultnonlin_small(z,p, cn) is the algorithm of the Mandel&Agol derivation.
    z = d/R_star, where d is the distance of the planet center to the LOS to the
    center of the star.
    sin(alpha) = d/a, with a the orbital distance (semi-major axis).
    so sin(alpha)*a/Rstar = d/a*a/Rstar = d/Rstar = z.
    a/Rstar happens to be a quantity that is well known from the transit light-
    curve. So z = sin(2pi phase)*a/Rstar. But this is in the limit of i = 90.

    From Cegla 2016 it follows that z = sqrt(xp^2 + yp^2). These are given
    as xp = a/Rstar sin(2pi phase) and yp = -a/Rstar * cos(2pi phase) * cos(i).

    The second quantity, p, is Rp/Rstar, also well known from the transit light-
    curve.

    cn is a four-element vector with the nonlinear limb darkening coefficients.
    If a shorter sequence is entered, the later values will be set to zero.
    By default I made it zero; i.e. the injected model does not take into
    account limb-darkening.
    """
    from lib.utils import typetest
    import lib.iansastropy as iap
    import numpy as np
    import pdb
    typetest('dp',dp,str)
    p=phase(dp)
    a_Rstar=paramget('aRstar',dp)
    Rp_Rstar=paramget('RpRstar',dp)
    i=paramget('inclination',dp)
    typetest('Rp_Rstar',a_Rstar,float)
    typetest('a_Rstar',a_Rstar,float)
    typetest('i',i,float)

    xp=np.sin(p*2.0*np.pi)*a_Rstar
    yp=np.cos(p*2.0*np.pi)*np.cos(np.radians(i))*a_Rstar
    z=np.sqrt(xp**2.0 + yp**2.0)
    transit=iap.occultnonlin_small(z,Rp_Rstar,[0.0,0.0])
    return transit
