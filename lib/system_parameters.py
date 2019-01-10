#This replaces calctimes and includes paramget.

def paramget(keyword,dp):
    import pdb

    try:
        f = open(dp+'config', 'r')
    except FileNotFoundError:
        raise Exception('Configfile does not exist at %s' % dp) from None
    x = f.read().splitlines()
    f.close()

    n_lines=len(x)
    keywords={}
    for i in range(0,n_lines-1):
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
    The output is a number."""
    import numpy as np
    import pdb
    from astropy.io import ascii
    from astropy.time import Time

    P=paramget('P',dp)
    Tc=paramget('Tc',dp)
    r=paramget('a',dp)
    duration=paramget('duration',dp)
    inclination=paramget('inclination',dp)



def rv(dp):
    """This program calculates the radial velocity in km/s for the planet in the
    data sequence provided in dp, the data-path. dp starts in the root folder,
    i.e. it starts with data/projectname/, and it ends with a slash.

    Example: v=rv(data/Kelt-9/night1/)
    The output is an array with length N, corresponding to N exposures."""

    import numpy as np
    import pdb
    from astropy.io import ascii
    from astropy.time import Time

    d=ascii.read(dp+'obs_times',names=['mjd','time','exptime','airmass'])
    t = Time(d['time'])
    jd = t.jd
    pdb.set_trace()
