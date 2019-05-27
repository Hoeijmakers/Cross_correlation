def read_EXPRES_e2ds(inpath,outname,air=True):
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









def read_HARPS_e2ds(inpath,outname,air=True,nowave=False):
    """THIS IS A PYTHON TRANSLATION OF READ_DATA (BELOW). IT SHOULD NOT WORK
    WITH A PATHS FILE, JUST A FOLDER THAT CONTAINS ONLY FITS FILES AND THEN
    IT WORKS FROM THE KEYWORDS TO DO EVERYTHING AUTOMATICALLY.

    WRITE GOOD TESTS AND DOCUMENTATION.

    Set the nowave keyword to True if the dataset has no wave files associated with it.
    This may happen if you downloaded ESO Advanced Data Products, which include
    reduced science e2ds's but not reduced wave e2ds's. The wavelength solution
    is still encoded in the fits header however, so we take it from there, instead."""
    import os
    import pdb
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import lib.utils as ut

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
    wavefile_used = []
    outpath = 'data/'+outname
    if os.path.exists(outpath) != True:
        os.makedirs(outpath)

    #ccftotal = 0 #This will hold the sum of the CCFs
    e2ds_count = 0
    sci_count = 0
    wave_count = 0
    ccf_count = 0
    for i in range(N):
        if filelist[i].endswith('e2ds_A.fits'):
            e2ds_count += 1
            print(filelist[i])
            data,hdr=fits.getdata(inpath+filelist[i],header=True)
            framename.append(filelist[i])
            header.append(hdr)
            type.append(hdr['HIERARCH ESO DPR CATG'])
            texp=np.append(texp,hdr['EXPTIME'])
            date.append(hdr['DATE-OBS'])
            mjd=np.append(mjd,hdr['MJD-OBS'])
            npx=np.append(npx,hdr['NAXIS1'])
            norders=np.append(norders,hdr['NAXIS2'])
            e2ds.append(data)
            if hdr['HIERARCH ESO DPR CATG'] == 'SCIENCE':
                sci_count += 1
                berv=np.append(berv,hdr['HIERARCH ESO DRS BERV'])
                airmass=np.append(airmass,0.5*(hdr['HIERARCH ESO TEL AIRM START']+hdr['HIERARCH ESO TEL AIRM END']))

                if nowave == True:
                    #Record which wavefile was used by the pipeline to
                        #create the wavelength solution.
                    wavefile_used.append(hdr['HIERARCH ESO DRS CAL TH FILE'])

                    wavedata=read_wave_from_HARPS_header(hdr)
                    wave.append(wavedata)
            else:
                berv=np.append(berv,np.nan)
                airmass=np.append(airmass,np.nan)
        if filelist[i].endswith('wave_A.fits'):
            if nowave == True:
                print("WARNING in read_HARPS_e2ds: nowave was set to True but a wave_A file")
                print("was detected. This wave file is now ignored in favor of the header.")
            else:
                wavedata=fits.getdata(inpath+filelist[i])
                wave.append(wavedata)
                wave_count += 1
        if filelist[i].endswith('ccf_G2_A.fits'):
            ccf,hdr=fits.getdata(inpath+filelist[i],header=True)
            if hdr['HIERARCH ESO DPR CATG'] == 'SCIENCE':
                #ccftotal+=ccf
                ccfs.append(ccf)
                ccfmjd=np.append(ccfmjd,hdr['MJD-OBS'])
                nrv=np.append(nrv,hdr['NAXIS1'])
                ccf_count += 1
    #Now we catch some errors:
    #-The above should have read a certain number of e2ds files.
    #-A certain number of these should be SCIENCE frames.
    #-There should be at least one WAVE file.
    #-All exposures should have the same number of spectral orders.
    #-All orders should have the same number of pixels (this is true for HARPS).
    #-The wave frame should have the same dimensions as the order frames.
    #-If nowave is set, test that all frames used the same wave_A calibrator.
    if ccf_count != sci_count:
        print("ERROR in read_HARPS_e2ds: There is a different number of science CCFs as there is science frames.")
        sys.exit()
    if e2ds_count == 0:
        print("ERROR in read_HARPS_e2ds: The input folder (%s) does not contain files ending in e2ds.fits." % inpath)
        sys.exit()
    if sci_count == 0:
        print("ERROR in read_HARPS_e2ds: The input folder (%2) contains e2ds files, but none of them are classified as SCIENCE frames with the HIERARCH ESO DPR CATG keyword.")
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
            for i in range(N):
                print('   '+framename[i]+'  %s' % wavefile_used[i])
            sys.exit()
    if np.shape(wave) != np.shape(e2ds[0]):
        print("ERROR in read_HARPS_e2ds: A wave file was detected but its shape (%s,%s) does not match that of the orders (%s,%s)" % (np.shape(wave)[0],np.shape(wave)[1],np.shape(e2ds[0])[0],np.shape(e2ds[0])[1]))

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
        ccforder = np.zeros((ccf_count,nrv))
        wave_axis = wave[i,:]/10.0#Convert to nm.
        print('CONSTRUCTING ORDER %s' % i)
        c = 0#To count the number of science frames that have passed. The counter
        #c is not equal to j because the list of files contains not only SCIENCE
        #frames.
        cc = 0#Same for ccfs
        for j in range(len(ccfsorting)):
            ccf=ccfs[ccfsorting[j]]
            ccforder[cc,:] = ccf[i,:]
            cc+=1
        for j in range(len(sorting)):#Loop over exposures
            if i ==0:
                print('---'+type[sorting[j]]+'  '+date[sorting[j]])
            if type[sorting[j]] == 'SCIENCE':
                exposure = e2ds[sorting[j]]
                order[c,:] = exposure[i,:]
                #Now I also need to write it to file.
                if i ==0:#Only do it the first time, not for every order.
                    line = str(mjd[sorting[j]])+'\t'+date[sorting[j]]+'\t'+str(texp[sorting[j]])+'\t'+str(airmass[sorting[j]])+'\t'+str(berv[sorting[j]])+'\t'+framename[sorting[j]]+'\n'
                    f.write(line)
                c+=1
        ccftotal+=ccforder
        fits.writeto(outpath+'/ccf_'+str(i)+'.fits',ccforder,overwrite=True)
        fits.writeto(outpath+'/order_'+str(i)+'.fits',order,overwrite=True)
        fits.writeto(outpath+'/wave_'+str(i)+'.fits',wave_axis,overwrite=True)
    fits.writeto(outpath+'/ccftotal.fits',ccftotal,overwrite=True)
    f.close()
    print('obs_times written to '+outpath)
    print('WARNING: FORMATTING IS STILL SCREWED UP!')
    print('FIGURE OUT HOW TO FORMAT THOSE LINES IN A MORE HUMAN READABLE WAY')
    print('WHEN YOU HAVE INTERNET AGAIN.')
    #pdb.set_trace()


def read_wave_from_HARPS_header(h):
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

    key_counter = 0
    for i in range(no):
        l = x*0.0
        for j in range(4):
            l += h['ESO DRS CAL TH COEFF LL%s' %key_counter]*x**j
            key_counter +=1
        wave[:,i] = l
    wave = wave.T
    return(wave)
