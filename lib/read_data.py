def read_EXPRES_e2ds_v5(inpath,outname,molecfit=False):
    """This reads the 5th version of the EXPRES reduction; files not ending in .spec.fits."""
    import os
    import pdb
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import lib.utils as ut
    import scipy.interpolate
    import copy
    from scipy.ndimage.interpolation import shift
    import lib.molecfit as mol
    import lib.operations as ops
    import pickle
    #First check the input:
    ut.typetest('inpath in read_EXPRES_e2ds ',inpath,str)
    ut.typetest('outname in read_EXPRES_e2ds ',outname,str)
    # ut.typetest('air in read_EXPRES_e2ds ',air,bool)
    if os.path.exists(inpath) != True:
        print("ERROR in read_EXPRES_e2ds: Data input path (%s) does not exist." % inpath)
        sys.exit()

    filelist=os.listdir(inpath)
    N=len(filelist)

    if len(filelist) == 0:
        print("ERROR in read_EXPRES_e2ds: input folder (%s) is empty." % inpath)
        sys.exit()

    #The following variables define the lists in which all the necessary data will be stored.
    framename=[]
    types=[]
    texp=np.array([])
    date=[]
    headers=[]
    mjd=np.array([])
    ccfmjd=np.array([])
    npx=np.array([])
    nrv=np.array([])
    norders=np.array([])
    e2ds=[]
    airmass=np.array([])
    berv=np.array([])
    wave=[]
    blaze=[]
    ccfs=[]
    outpath = 'data/'+outname
    if os.path.exists(outpath) != True:
        os.makedirs(outpath)

    #ccftotal = 0 #This will hold the sum of the CCFs
    e2ds_count = 0
    sci_count = 0
    wave_count = 0
    ccf_count = 0
    for i in range(N):
        if filelist[i].endswith('.fits'):
            e2ds_count += 1
            print(filelist[i])
            hdul = fits.open(inpath+filelist[i])
            hdr = copy.deepcopy(hdul[0].header)
            # pdb.set_trace()

            hdr2 = copy.deepcopy(hdul[1].header)
            data = copy.deepcopy(hdul[1].data['spectrum'])
            cont = copy.deepcopy(hdul[1].data['blaze'])
            wl = copy.deepcopy(hdul[1].data['wavelength'])#in Angstrom
            hdul.close()
            del hdul[0].data
            # data,hdr=fits.getdata(inpath+filelist[i],header=True)
            framename.append(filelist[i])
            types.append(hdr['OBSTYPE'])
            texp=np.append(texp,float(hdr['AEXPTIME']))
            date.append(hdr['DATE-OBS'].replace(' ','T'))
            mjd=np.append(mjd,float(hdr['TELMJD']))
            npx=np.append(npx,hdr['NAXIS1'])
            norders=np.append(norders,hdr['NAXIS2'])
            headers.append(hdr)
            e2ds.append(data)
            wave.append(wl)
            blaze.append(cont)
            if hdr['OBSTYPE'] == 'Science':
                sci_count += 1
                berv=np.append(berv,float(hdr2['WS_BARY']))
                airmass=np.append(airmass,float(hdr['AIRMASS']))
            else:
                berv=np.append(berv,np.nan)
                airmass=np.append(airmass,np.nan)



    #Now we catch some errors:
    #-The above should have read a certain number of e2ds files.
    #-A certain number of these should be SCIENCE frames.
    #-There should be at least one WAVE file.
    #-All exposures should have the same number of spectral orders.
    #-All orders should have the same number of pixels (this is true for HARPS).
    #-The wave frame should have the same dimensions as the order frames.
    if e2ds_count == 0:
        print("ERROR in read_EXPRES_e2ds: The input folder (%s) does not contain files ending in .spec.fits." % inpath)
        sys.exit()
    if sci_count == 0:
        print("ERROR in read_EXPRES_e2ds: The input folder (%2) contains e2ds files, but none of them are classified as Science frames with the OBSTYPE keyword.")
        print("These are the files and their types:")
        for i in range(len(types)):
            print('   '+framename[i]+'  %s' % types[i])
        sys.exit()
    if np.max(np.abs(norders-norders[0])) == 0:
        norders=int(norders[0])
    else:
        print("ERROR in read_EXPRES_e2ds: Not all files have the same number of orders.")
        print("These are the files and their number of orders:")
        for i in range(len(types)):
            print('   '+framename[i]+'  %s' % norders[i])
        sys.exit()
    if np.max(np.abs(npx-npx[0])) == 0:
        npx=int(npx[0])
    else:
        print("ERROR IN read_EXPRES_e2ds: Not all files have the same number of pixels.")
        print("These are the files and their number of pixels:")
        for i in range(len(types)):
            print('   '+framename[i]+'  %s' % npx[i])
        sys.exit()






    #Now we need to construct the s1d files:
    # Need to deal with overlap.
    list_of_wls = []
    list_of_spectra = []

    orderfraction=0.85#This fraction of the order is parsed into the s1d, centered on the middle of the order.
    #meaning the edges are not included in the s1d.
    longpass = 5600.0#Wavelengths below which are ignored.

    N_orders = len(e2ds[0])
    N_exp = len(e2ds)
    dwl = np.inf
    minwl=np.inf
    maxwl=-1.0*np.inf

    for i in range(N_orders):
        wl_i = np.array(wave[0][i])
        minwl=np.min([np.min(wl_i),minwl])#Update minwl as we loop through the orders; update if min(wl_i) is smaller than current wl.
        maxwl=np.max([np.max(wl_i),maxwl])#Vice versa
        wl_diff = np.abs(wl_i - shift(wl_i,1,cval=np.NaN))#Search for the smallest delta-wl
        dwl  = np.min([dwl,np.nanmin(wl_diff)])


    wl_s1d=np.arange(minwl-10.0,maxwl+10.0,dwl)#Create new wlgrid for all the s1ds.

    trigger = 0
    if trigger == 1:
        print('---Building s1d spectra')
        for i in range(N_exp):
            list_of_wls.append(wl_s1d)
            fx_s1d = wl_s1d*0.0+1.0#Everywhere it is 1.0.
            for j in range(N_orders):
                nwl=len(wave[i][0])
                start_index = int((1.0-orderfraction)/2.0*nwl)
                end_index = int(start_index+0.5*nwl)
                current_wl = ops.vactoair(wave[i][j][start_index:end_index])
                current_spec = (e2ds[i][j]/blaze[i][j])[start_index:end_index]
                wle,fxe = ops.envelope(current_wl,current_spec,binsize=1.0,selfrac=0.3)
                env_i = scipy.interpolate.interp1d(wle,fxe,bounds_error=False,fill_value='extrapolate')
                current_spec_flat = current_spec/env_i(current_wl)
                fx_i=scipy.interpolate.interp1d(current_wl,current_spec_flat,bounds_error=False,fill_value=1.0)
                fx_s1d *= fx_i(wl_s1d)
            fx_s1d[(wl_s1d<longpass)] = 1.0
            # list_of_wls.append(wl_s1d)
            list_of_spectra.append(fx_s1d)
            ut.statusbar(i,N_exp)
            with open('temp.pkl', 'wb') as f: pickle.dump(list_of_spectra,f)
    else:
        pickle_in = open('temp.pkl',"rb")
        list_of_spectra=pickle.load(pickle_in)


    #Ok, so now we should have ended up with a number of lists that contain all
    #the relevant information of our science frames.
    #We determine how to sort the resulting lists in time:
    sorting = np.argsort(mjd)
    ccfsorting = np.argsort(ccfmjd)
    ccftotal = 0.0



    #First sort the s1d files for application of molecfit.
    if molecfit == True:

        #We add all necessary fits keywords to the headers.
        for i in range(len(headers)):
            headers[i]['HIERARCH ESO DRS BERV'] = 0.0 #We have already taken berv-corrected spectra, so the berv is zero.
            headers[i]['CRVAL1'] = min(wl_s1d)
            headers[i]['CDELT1'] = dwl#The above three lines are to fool mol.write_file_to_molecfit.
            alt_string = headers[i]['TELALT'].split(':')
            headers[i]['ALTDEG'] = float(alt_string[0])+float(alt_string[1])/60.0+float(alt_string[2])/3600
            headers[i]['MJD-OBS'] = float(headers[i]['TELMJD'])
            headers[i]['UTCS'] = (float(headers[i]['TELMJD'])%1.0)*86400.0
            headers[i]['ELEV'] = float(headers[i]['SITEELEV'])
            headers[i]['LONGIT'] = float(headers[i]['SITELONG'])
            headers[i]['LATITU'] = float(headers[i]['SITELAT'])
            # list_of_wls[i]=ops.vactoair(list_of_wls[i]/10.0)*10.0

        headers_sorted=[]
        s1d_sorted=[]
        for i in range(len(sorting)):
            headers_sorted.append(headers[int(sorting[i])])
            s1d_sorted.append(list_of_spectra[int(sorting[i])])

        print('Molecfit will be executed onto the files in this order:')
        for x in headers_sorted:
            print(x['DATE-OBS'])
        list_of_wls,list_of_trans = mol.do_molecfit(headers_sorted,s1d_sorted,load_previous=False,mode='EXPRES')
        mol.write_telluric_transmission_to_file(list_of_wls,list_of_trans,outpath+'telluric_transmission_spectra.pkl')




    #Now we loop over all exposures and collect the i-th order from each exposure,
    #put these into a new matrix and save them to FITS images:
    file=open(outpath+'obs_times','w',newline='\n')
    headerline = 'MJD'+'\t'+'DATE'+'\t'+'EXPTIME'+'\t'+'MEAN AIRMASS'+'\t'+'BERV (km/s)'+'\t'+'FILE NAME'
    for i in range(norders):
        order = np.zeros((sci_count,npx))
        wave_axis = wave[0][i,:]/10.0#Convert to nm.
        print('CONSTRUCTING ORDER %s' % i)
        c = 0#To count the number of science frames that have passed. The counter
        #c is not equal to j because the list of files contains not only SCIENCE
        #frames.
        for j in range(len(sorting)):#Loop over exposures
            if i ==0:
                print('---'+types[sorting[j]]+'  '+date[sorting[j]])
            if types[sorting[j]] == 'Science':
                exposure = e2ds[sorting[j]]
                wave2d = wave[sorting[j]]
                spec_i=scipy.interpolate.interp1d(wave2d[i,:]/10.0,exposure[i,:],fill_value=np.nan,bounds_error=False)
                order[c,:] = spec_i(wave_axis)
                #Now I also need to write it to file.
                if i ==0:#Only do it the first time, not for every order.
                    line = str(mjd[sorting[j]])+'\t'+date[sorting[j]]+'\t'+str(texp[sorting[j]])+'\t'+str(airmass[sorting[j]])+'\t'+str(berv[sorting[j]]/1000.0)+'\t'+framename[sorting[j]]+'\n'
                    file.write(line)
                c+=1
        fits.writeto(outpath+'/order_'+str(i)+'.fits',order,overwrite=True)
        fits.writeto(outpath+'/wave_'+str(i)+'.fits',wave_axis,overwrite=True)
    file.close()
    print('obs_times written to '+outpath+'/')
    print('WARNING: FORMATTING IS STILL SCREWED UP!')
    print('FIGURE OUT HOW TO FORMAT THOSE LINES IN A MORE HUMAN READABLE WAY')
    print('WHEN YOU HAVE INTERNET AGAIN.')







def read_EXPRES_e2ds_v4(inpath,outname,air=True):
    """THIS IS A PYTHON TRANSLATION OF READ_DATA (BELOW). IT SHOULD NOT WORK
    WITH A PATHS FILE, JUST A FOLDER THAT CONTAINS ONLY FITS FILES AND THEN
    IT WORKS FROM THE KEYWORDS TO DO EVERYTHING AUTOMATICALLY.

    WRITE GOOD TESTS AND DOCUMENTATION.

    The air keyword is depricated."""
    import os
    import pdb
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import lib.utils as ut
    import scipy.interpolate

    #First check the input:
    ut.typetest('inpath in read_EXPRES_e2ds ',inpath,str)
    ut.typetest('outname in read_EXPRES_e2ds ',outname,str)
    ut.typetest('air in read_EXPRES_e2ds ',air,bool)
    if os.path.exists(inpath) != True:
        print("ERROR in read_EXPRES_e2ds: Data input path (%s) does not exist." % inpath)
        sys.exit()

    filelist=os.listdir(inpath)
    N=len(filelist)

    if len(filelist) == 0:
        print("ERROR in read_EXPRES_e2ds: input folder (%s) is empty." % inpath)
        sys.exit()

    #The following variables define the lists in which all the necessary data will be stored.
    framename=[]
    type=[]
    texp=np.array([])
    date=[]
    mjd=np.array([])
    ccfmjd=np.array([])
    npx=np.array([])
    nrv=np.array([])
    norders=np.array([])
    e2ds=[]
    airmass=np.array([])
    berv=np.array([])
    wave=[]
    ccfs=[]
    outpath = 'data/'+outname
    if os.path.exists(outpath) != True:
        os.makedirs(outpath)

    #ccftotal = 0 #This will hold the sum of the CCFs
    e2ds_count = 0
    sci_count = 0
    wave_count = 0
    ccf_count = 0
    for i in range(N):
        if filelist[i].endswith('.spec.fits'):
            e2ds_count += 1
            print(filelist[i])
            data,hdr=fits.getdata(inpath+filelist[i],header=True)
            framename.append(filelist[i])
            type.append(hdr['OBSTYPE'])
            texp=np.append(texp,hdr['AEXPTIME'])
            date.append(hdr['DATE-OBS'].replace(' ','T'))
            mjd=np.append(mjd,hdr['TELMJD'])
            npx=np.append(npx,hdr['NAXIS1'])
            norders=np.append(norders,hdr['NAXIS2'])
            e2ds.append(data[0,:,:])
            wave.append(data[1,:,:])
            if hdr['OBSTYPE'] == 'Science':
                sci_count += 1
                berv=np.append(berv,hdr['WS_BARY'])
                airmass=np.append(airmass,hdr['AIRMASS'])
            else:
                berv=np.append(berv,np.nan)
                airmass=np.append(airmass,np.nan)


    #Now we catch some errors:
    #-The above should have read a certain number of e2ds files.
    #-A certain number of these should be SCIENCE frames.
    #-There should be at least one WAVE file.
    #-All exposures should have the same number of spectral orders.
    #-All orders should have the same number of pixels (this is true for HARPS).
    #-The wave frame should have the same dimensions as the order frames.
    if e2ds_count == 0:
        print("ERROR in read_EXPRES_e2ds: The input folder (%s) does not contain files ending in .spec.fits." % inpath)
        sys.exit()
    if sci_count == 0:
        print("ERROR in read_EXPRES_e2ds: The input folder (%2) contains e2ds files, but none of them are classified as Science frames with the OBSTYPE keyword.")
        print("These are the files and their types:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % type[i])
        sys.exit()
    if np.max(np.abs(norders-norders[0])) == 0:
        norders=int(norders[0])
    else:
        print("ERROR in read_EXPRES_e2ds: Not all files have the same number of orders.")
        print("These are the files and their number of orders:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % norders[i])
        sys.exit()
    if np.max(np.abs(npx-npx[0])) == 0:
        npx=int(npx[0])
    else:
        print("ERROR IN read_EXPRES_e2ds: Not all files have the same number of pixels.")
        print("These are the files and their number of pixels:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % npx[i])
        sys.exit()


    #Ok, so now we should have ended up with a number of lists that contain all
    #the relevant information of our science frames.
    #We determine how to sort the resulting lists in time:
    sorting = np.argsort(mjd)
    ccfsorting = np.argsort(ccfmjd)
    ccftotal = 0.0
    #Now we loop over all exposures and collect the i-th order from each exposure,
    #put these into a new matrix and save them to FITS images:
    f=open(outpath+'obs_times','w',newline='\n')
    headerline = 'MJD'+'\t'+'DATE'+'\t'+'EXPTIME'+'\t'+'MEAN AIRMASS'+'\t'+'BERV (km/s)'+'\t'+'FILE NAME'
    for i in range(norders):
        order = np.zeros((sci_count,npx))
        wave_axis = wave[0][i,:]/10.0#Convert to nm.
        print('CONSTRUCTING ORDER %s' % i)
        c = 0#To count the number of science frames that have passed. The counter
        #c is not equal to j because the list of files contains not only SCIENCE
        #frames.
        for j in range(len(sorting)):#Loop over exposures
            if i ==0:
                print('---'+type[sorting[j]]+'  '+date[sorting[j]])
            if type[sorting[j]] == 'Science':
                exposure = e2ds[sorting[j]]
                wave2d = wave[sorting[j]]
                spec_i=scipy.interpolate.interp1d(wave2d[i,:]/10.0,exposure[i,:],fill_value=np.nan)
                order[c,:] = spec_i(wave_axis)
                #Now I also need to write it to file.
                if i ==0:#Only do it the first time, not for every order.
                    line = str(mjd[sorting[j]])+'\t'+date[sorting[j]]+'\t'+str(texp[sorting[j]])+'\t'+str(airmass[sorting[j]])+'\t'+str(berv[sorting[j]]/1000.0)+'\t'+framename[sorting[j]]+'\n'
                    f.write(line)
                c+=1
        fits.writeto(outpath+'/order_'+str(i)+'.fits',order,overwrite=True)
        fits.writeto(outpath+'/wave_'+str(i)+'.fits',wave_axis,overwrite=True)
    f.close()
    print('obs_times written to '+outpath+'/')
    print('WARNING: FORMATTING IS STILL SCREWED UP!')
    print('FIGURE OUT HOW TO FORMAT THOSE LINES IN A MORE HUMAN READABLE WAY')
    print('WHEN YOU HAVE INTERNET AGAIN.')









def read_HARPS_e2ds(inpath,outname,air=True,nowave=False,molecfit=True,mode='HARPS'):
    """THIS IS A PYTHON TRANSLATION OF READ_DATA (BELOW). IT SHOULD NOT WORK
    WITH A PATHS FILE, JUST A FOLDER THAT CONTAINS ONLY FITS FILES AND THEN
    IT WORKS FROM THE KEYWORDS TO DO EVERYTHING AUTOMATICALLY.

    WRITE GOOD TESTS AND DOCUMENTATION.

    ALSO, ULTIMATELY THIS WILL NEED A WRAPPER THAT CAN SWITCH BETWEEN DIFFERENT STANDARD DATASETS.
    IN THE CASE OF UVES (AND MAYBE MOST OTHER DATASETS) IT WILL NEED TO DEAL WITH BERV CORRECTIONS.
    GREAT WAY TO DO THIS IS HERE: https://docs.astropy.org/en/stable/coordinates/velocities.html
    DID THAT WITH JEHAN FOR 55 CNC E.

    Set the nowave keyword to True if the dataset has no wave files associated with it.
    This may happen if you downloaded ESO Advanced Data Products, which include
    reduced science e2ds's but not reduced wave e2ds's. The wavelength solution
    is still encoded in the fits header however, so we take it from there, instead.


    IF IN THE FUTURE A BERV KEYWORD WOULD BE MISSING, I HAVE INCLUDED AN ASTROPY
    IMPLEMENTATION THAT ACCURATELY CALCULATES THE BERV FROM THE MJD. SEE SYSTEM_PARAMETERS.PY
    """
    import os
    import pdb
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import lib.utils as ut
    import lib.molecfit as mol
    import pyfits
    import copy
    import scipy.interpolate as interp
    import pickle

    #First check the input:
    ut.typetest('inpath in read_HARPS_e2ds ',inpath,str)
    ut.typetest('outname in read_HARPS_e2ds ',outname,str)
    ut.typetest('air in read_HARPS_e2ds ',air,bool)
    if os.path.exists(inpath) != True:
        print("ERROR in read_HARPS_e2ds: Data input path (%s) does not exist." % inpath)
        sys.exit()

    filelist=os.listdir(inpath)
    N=len(filelist)

    if len(filelist) == 0:
        print("ERROR in read_HARPS_e2ds: input folder (%s) is empty." % inpath)
        sys.exit()

    #The following variables define the lists in which all the necessary data will be stored.
    framename=[]
    header=[]
    s1dhdr=[]
    type=[]
    texp=np.array([])
    date=[]
    mjd=np.array([])
    ccfmjd=np.array([])
    s1dmjd=np.array([])
    npx=np.array([])
    nrv=np.array([])
    norders=np.array([])
    e2ds=[]
    s1d=[]
    wave1d=[]
    airmass=np.array([])
    berv=np.array([])
    wave=[]
    blaze=[]
    ccfs=[]
    wavefile_used = []
    outpath = ut.path('data/'+outname)
    if os.path.exists(outpath) != True:
        os.makedirs(outpath)

    #ccftotal = 0 #This will hold the sum of the CCFs
    e2ds_count = 0
    sci_count = 0
    wave_count = 0
    ccf_count = 0
    blaze_count = 0
    s1d_count = 0


    #MODE SWITCHING HERE:
    if mode == 'HARPS':
        catkeyword = 'HIERARCH ESO DPR CATG'
        bervkeyword = 'HIERARCH ESO DRS BERV'
        thfilekeyword = 'HIERARCH ESO DRS CAL TH FILE'
        Zstartkeyword = 'HIERARCH ESO TEL AIRM START'
        Zendkeyword = 'HIERARCH ESO TEL AIRM END'
    if mode == 'HARPSN':
        catkeyword = 'OBS-TYPE'
        bervkeyword = 'HIERARCH TNG DRS BERV'
        thfilekeyword = 'HIERARCH TNG DRS CAL TH FILE'
        Zstartkeyword = 'AIRMASS'
        Zendkeyword = 'AIRMASS'#These are the same because HARPSN doesnt have start and end keywords.
        #Down there, the airmass is averaged, so there is no problem in taking the average of the same number.



    for i in range(N):
        if filelist[i].endswith('e2ds_A.fits'):
            e2ds_count += 1
            print(filelist[i])
            #data,hdr=fits.getdata(inpath+filelist[i],header=True)

            hdul = fits.open(inpath+filelist[i])
            data = copy.deepcopy(hdul[0].data)
            hdr = hdul[0].header
            hdul.close()

            del hdul[0].data
            if hdr[catkeyword] == 'SCIENCE':
                framename.append(filelist[i])
                header.append(hdr)
                type.append(hdr[catkeyword])
                texp=np.append(texp,hdr['EXPTIME'])
                date.append(hdr['DATE-OBS'])
                mjd=np.append(mjd,hdr['MJD-OBS'])
                npx=np.append(npx,hdr['NAXIS1'])
                norders=np.append(norders,hdr['NAXIS2'])
                e2ds.append(data)

                sci_count += 1
                berv=np.append(berv,hdr[bervkeyword])
                airmass=np.append(airmass,0.5*(hdr[Zstartkeyword]+hdr[Zendkeyword]))

                if nowave == True:
                    #Record which wavefile was used by the pipeline to
                        #create the wavelength solution.
                    wavefile_used.append(hdr[thfilekeyword])

                    wavedata=read_wave_from_HARPS_header(hdr,mode=mode)
                    wave.append(wavedata)
            # else:
                # berv=np.append(berv,np.nan)
                # airmass=np.append(airmass,np.nan)
        if filelist[i].endswith('wave_A.fits'):
            print(filelist[i]+' (wave)')
            if nowave == True:
                print("WARNING in read_HARPS_e2ds: nowave was set to True but a wave_A file")
                print("was detected. This wave file is now ignored in favor of the header.")
            else:
                wavedata=fits.getdata(inpath+filelist[i])
                wave.append(wavedata)
                wave_count += 1
        if filelist[i].endswith('ccf_G2_A.fits'):
            #ccf,hdr=fits.getdata(inpath+filelist[i],header=True)
            hdul = fits.open(inpath+filelist[i])
            ccf = copy.deepcopy(hdul[0].data)
            hdr = hdul[0].header
            hdul.close()
            del hdul[0].data

            if hdr[catkeyword] == 'SCIENCE':
                #ccftotal+=ccf
                ccfs.append(ccf)
                ccfmjd=np.append(ccfmjd,hdr['MJD-OBS'])
                nrv=np.append(nrv,hdr['NAXIS1'])
                ccf_count += 1

        if filelist[i].endswith('blaze_A.fits'):
                print(filelist[i]+' (blaze)')
                blazedata=fits.getdata(inpath+filelist[i])
                blaze.append(blazedata)
                blaze_count += 1
        if filelist[i].endswith('s1d_A.fits'):
            hdul = fits.open(inpath+filelist[i])
            data_1d = copy.deepcopy(hdul[0].data)
            hdr = hdul[0].header
            hdul.close()
            del hdul[0].data
            if hdr[catkeyword] == 'SCIENCE':
                s1d.append(data_1d)
                s1dhdr.append(hdr)
                s1dmjd=np.append(s1dmjd,hdr['MJD-OBS'])
                s1d_count += 1
    #Now we catch some errors:
    #-The above should have read a certain number of e2ds files.
    #-A certain number of these should be SCIENCE frames.
    #-There should be at least one WAVE file.
    #-All exposures should have the same number of spectral orders.
    #-All orders should have the same number of pixels (this is true for HARPS).
    #-The wave frame should have the same dimensions as the order frames.
    #-If nowave is set, test that all frames used the same wave_A calibrator.
    #-The blaze file needs to have the same shape as the e2ds files.
    #-The number of s1d files should be the same as the number of e2ds files.


    if ccf_count != sci_count:
        print("ERROR in read_HARPS_e2ds: There is a different number of science CCFs as there is science frames.")
        sys.exit()
    # if e2ds_count != s1d_count:
    #     print('ERROR in read_HARPS_e2ds: The numbers of 1ds and e2ds files are different.')
    #     print("These are the files and their types:")
    #     for i in range(len(type)):
    #         print('   '+framename[i]+'  %s' % type[i])
    #     sys.exit()
    if e2ds_count == 0:
        print("ERROR in read_HARPS_e2ds: The input folder (%s) does not contain files ending in e2ds.fits." % inpath)
        sys.exit()
    if sci_count == 0:
        print("ERROR in read_HARPS_e2ds: The input folder (%2) contains e2ds files, but none of them are classified as SCIENCE frames with the HIERARCH ESO DPR CATG/OBS-TYPE keyword.")
        print("These are the files and their types:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % type[i])
        sys.exit()
    if np.max(np.abs(norders-norders[0])) == 0:
        norders=int(norders[0])
    else:
        print("ERROR in read_HARPS_e2ds: Not all files have the same number of orders.")
        print("These are the files and their number of orders:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % norders[i])
        sys.exit()
    if np.max(np.abs(npx-npx[0])) == 0:
        npx=int(npx[0])
    else:
        print("ERROR IN read_HARPS_e2ds: Not all files have the same number of pixels.")
        print("These are the files and their number of pixels:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % npx[i])
        sys.exit()
    if np.max(np.abs(nrv-nrv[0])) == 0:
        nrv=int(nrv[0])
    else:
        print("ERROR IN read_HARPS_e2ds: Not all files have the same number of pixels.")
        print("These are the files and their number of pixels:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % npx[i])
        sys.exit()
    if wave_count >= 1:
        wave=wave[0]#SELECT ONLY THE FIRST WAVE FRAME. The rest is ignored.
    else:
        if nowave == False:
            print("ERROR in read_HARPS_e2ds: No wave_A.fits file was detected.")
            print("These are the files in the folder:")
            for i in range(N):
                print(filelist[i])
            print("This may have happened if you downloaded the HARPS data from the")
            print("ADP query form, which doesn't include wave_A files (as far as I")
            print("have seen). Set the /nowave keyword in your call to read_HARPS_e2ds")
            print("if you indeed do not expect a wave_A file to be present.")
    if nowave == True:
        if all(x == wavefile_used[0] for x in wavefile_used):
            print("Nowave is set, and simple wavelength calibration extraction")
            print("works, as all files in the dataset used the same wave_A file.")
            wave=wave[0]
        else:
            print("ERROR IN read_HARPS_e2ds: Nowave is set, but not all files")
            print("in the dataset used the same wave_A file when the pipeline was")
            print("run. Catching this requres an interpolation step that is currently")
            print("not yet implemented. Exiting. These are the filenames and their")
            print("wave_A file used:")
            for i in range(N-1):
                print('   '+framename[i]+'  %s' % wavefile_used[0])
            wave=wave[0]
            print("I ALLOW YOU TO CONTINUE BUT USING ONLY THE FIRST WAVELENGTH")
            print("SOLUTION. A PART OF THE DATA MAY BE AFFECTED BY HAVING ASSUMED")
            print("THE WRONG SOLUTION. If you are doing transits, you don't need")
            print("this kind of precision.")
    if blaze_count >= 1:
        blaze=blaze[0]#SELECT ONLY THE FIRST WAVE FRAME. The rest is ignored.
    if np.shape(wave) != np.shape(e2ds[0]):
        print("ERROR in read_HARPS_e2ds: A wave file was detected but its shape (%s,%s) does not match that of the orders (%s,%s)" % (np.shape(wave)[0],np.shape(wave)[1],np.shape(e2ds[0])[0],np.shape(e2ds[0])[1]))
    if np.shape(blaze) != np.shape(e2ds[0]) and blaze_count > 0:
        print("ERROR in read_HARPS_e2ds: A blaze file was detected but its shape (%s,%s) does not match that of the orders (%s,%s)" % (np.shape(blaze)[0],np.shape(wave)[1],np.shape(e2ds[0])[0],np.shape(e2ds[0])[1]))
    if len(s1dhdr) != len(e2ds) and molecfit == True:
        print('ERROR in read_HARPS_e2ds: The number of s1d SCIENCE files and e2ds SCIENCE files is not the same. (%s vs %s)' % (len(s1dhdr),len(e2ds)))
        print('Switching off the molecfit option will suppress this error.')



    #Ok, so now we should have ended up with a number of lists that contain all
    #the relevant information of our science frames.
    #We determine how to sort the resulting lists in time:
    sorting = np.argsort(mjd)
    ccfsorting = np.argsort(ccfmjd)
    s1dsorting = np.argsort(s1dmjd)

    if mode == 'HARPSN':
        for i in range(len(header)):
            s1dhdr[i]['TELALT'] = np.degrees(float(s1dhdr[i]['EL']))
            s1dhdr[i]['UTC'] = (float(s1dhdr[i]['MJD-OBS'])%1.0)*86400.0




    # sys.exit()




    #First sort the s1d files for application of molecfit.
    if molecfit == True:
        s1dhdr_sorted=[]
        s1d_sorted=[]
        for i in range(len(s1dsorting)):
            s1dhdr_sorted.append(s1dhdr[s1dsorting[i]])
            s1d_sorted.append(s1d[s1dsorting[i]])

        # print('Molecfit will be executed onto the files in this order:')
        # for x in s1dhdr_sorted:
        #     print(x['DATE-OBS'])
        list_of_wls,list_of_trans = mol.do_molecfit(s1dhdr_sorted,s1d_sorted,load_previous=False,mode=mode)
        mol.write_telluric_transmission_to_file(list_of_wls,list_of_trans,outpath+'telluric_transmission_spectra.pkl')











    ccftotal = 0.0
    #Now we loop over all exposures and collect the i-th order from each exposure,
    #put these into a new matrix and save them to FITS images:
    f=open(outpath+'obs_times','w',newline='\n')
    headerline = 'MJD'+'\t'+'DATE'+'\t'+'EXPTIME'+'\t'+'MEAN AIRMASS'+'\t'+'BERV (km/s)'+'\t'+'FILE NAME'
    for i in range(norders):
        order = np.zeros((sci_count,npx))
        ccforder = np.zeros((ccf_count,nrv))
        wave_axis = wave[i,:]/10.0#Convert to nm.
        print('CONSTRUCTING ORDER %s' % i)
        c = 0#To count the number of science frames that have passed. The counter
        # c is not equal to j because the list of files contains not only SCIENCE
        # frames.
        cc = 0#Same for ccfs
        for j in range(len(ccfsorting)):
            ccf=ccfs[ccfsorting[j]]
            ccforder[cc,:] = ccf[i,:]
            cc+=1
        for j in range(len(sorting)):#Loop over exposures
            if i ==0:
                print('---'+type[sorting[j]]+'  '+date[sorting[j]])
            if type[sorting[j]] == 'SCIENCE':#This check may be redundant.
                exposure = e2ds[sorting[j]]
                order[c,:] = exposure[i,:]
                #T_i = interp.interp1d(list_of_wls[j],list_of_trans[j])#This should be time-sorted, just as the e2ds files.
                #Do a manual check here that the MJDs are identical.
                #Also, determiine what to do with airtovac.
                #tel_order[c,:] = T_i[wave_axis]
                #Now I also need to write it to file.
                if i ==0:#Only do it the first time, not for every order.
                    line = str(mjd[sorting[j]])+'\t'+date[sorting[j]]+'\t'+str(texp[sorting[j]])+'\t'+str(airmass[sorting[j]])+'\t'+str(berv[sorting[j]])+'\t'+framename[sorting[j]]+'\n'
                    f.write(line)
                c+=1
        ccftotal+=ccforder



        fits.writeto(outpath+'ccf_'+str(i)+'.fits',ccforder,overwrite=True)
        fits.writeto(outpath+'order_'+str(i)+'.fits',order,overwrite=True)
        fits.writeto(outpath+'wave_'+str(i)+'.fits',wave_axis,overwrite=True)
    fits.writeto(outpath+'ccftotal.fits',ccftotal,overwrite=True)
    f.close()
    print('Time-table written to '+outpath+'obs_times')
    print('WARNING: FORMATTING IS STILL SCREWED UP!')
    print('FIGURE OUT HOW TO FORMAT THOSE LINES IN A MORE HUMAN READABLE WAY')
    print('WHEN YOU HAVE INTERNET AGAIN.')


    #ADD OUTPUT OF THE TELLURIC TRANSMISSION SPECTRUM HERE
    #Also output the list of trans spectra to the datafolder while checking that the
    #dates are the same as in the just exported config file.



def read_wave_from_HARPS_header(h,mode='HARPS'):
    """This reads the wavelength solution from the HARPS header keywords that
    encode the coefficients as a 4-th order polynomial."""
    import numpy as np
    import sys
    import pdb
    import lib.functions as fun
    import matplotlib.pyplot as plt
    import lib.utils as ut
    npx = h['NAXIS1']
    no = h['NAXIS2']
    x = fun.findgen(npx)
    wave=np.zeros((npx,no))

    if mode == 'HARPS':
        coeffkeyword = 'ESO'
    if mode == 'HARPSN':
        coeffkeyword = 'TNG'
    key_counter = 0
    for i in range(no):
        l = x*0.0
        for j in range(4):
            l += h[coeffkeyword+' DRS CAL TH COEFF LL%s' %key_counter]*x**j
            key_counter +=1
        wave[:,i] = l
    wave = wave.T
    return(wave)
