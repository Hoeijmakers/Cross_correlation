#This replaces calctimes and includes paramget.

day=86400.0  #seconds
AU=1.496e+11 #meters
Msun=1.989e+30 #kg
Mjup=1.898e27 #kg
c=299792458.0 #m/s

def typetest(varname,var,vartype):
    if isinstance(varname,str) != True:
        raise Exception("Unit error in unit test: varname should be of class string.")
    if isinstance(var,vartype) != True:
        raise Exception("Unit error: %s should be of class string." % varname)


def paramget(keyword,dp):
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
    """This program calculates the orbital velocity in km/s for the planet in the
    data sequence provided in dp, the data-path. dp starts in the root folder,
    i.e. it starts with data/projectname/, and it ends with a slash.

    Example: v=v_orb(data/Kelt-9/night1/)
    The output is a number in km/s."""
    typetest('dp',dp,str)

    import numpy as np
    import pdb
    P=paramget('P',dp)*day
    r=paramget('a',dp)*AU

    typetest('P',P,float)
    typetest('r',r,float)

    if P <= 0:
        raise Exception("P <= zero. Check configfile at %s." % dp)
    if r <= 0:
        raise Exception("r <= zero. Check configfile at %s." % dp)
    return 2.0*np.pi*r/P/1000.0




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
    typetest('dp',dp,str)

    import numpy as np
    import pdb
    from astropy.io import ascii
    from astropy.time import Time
    from astropy import units as u, coordinates as coord

    d=ascii.read(dp+'obs_times',names=['mjd','time','exptime','airmass'])
    t = Time(d['time'],scale='utc', location=coord.EarthLocation.of_site('paranal'))
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


def RV(dp):
    """This program calculates the radial velocity in km/s for the planet in the
    data sequence provided in dp, the data-path. dp starts in the root folder,
    i.e. it starts with data/projectname/, and it ends with a slash.

    Example: v=RV('data/Kelt-9/night1/')
    The output is an array with length N, corresponding to N exposures.
    The radial velocity is provided in km/s."""
    typetest('dp',dp,str)
    import numpy as np
    p=phase(dp)
    i=paramget('inclination',dp)
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
    typetest('dp',dp,str)
    import numpy as np
    from astropy.io import ascii
    d=ascii.read(dp+'obs_times',names=['mjd','time','exptime','airmass'])
    Texp=d['exptime'].astype('float')
    vorb=v_orb(dp)
    p=phase(dp)
    P=paramget('P',dp)
    i=paramget('inclination',dp)
    dRV=vorb*np.cos(2.0*np.pi*p)*2.0*np.pi/(P*day)*np.sin(np.radians(i))
    return abs(dRV*Texp)
