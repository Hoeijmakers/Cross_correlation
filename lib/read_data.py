def read_HARPS_e2ds(inpath,outname,air=True):
    """THIS IS A PYTHON TRANSLATION OF READ_DATA (BELOW). IT SHOULD NOT WORK
    WITH A PATHS FILE, JUST A FOLDER THAT CONTAINS ONLY FITS FILES AND THEN
    IT WORKS FROM THE KEYWORDS TO DO EVERYTHING AUTOMATICALLY.

    UNDER CONSTRUCTION"""
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
        if filelist[i].endswith('e2ds_A.fits'):
            e2ds_count += 1
            print(filelist[i])
            data,hdr=fits.getdata(inpath+filelist[i],header=True)
            framename.append(filelist[i])
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
            else:
                berv=np.append(berv,np.nan)
                airmass=np.append(airmass,np.nan)
        if filelist[i].endswith('wave_A.fits'):
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
        print("ERROR in read_HARPS_e2ds: No wave_A.fits file was detected.")
        print("These are the files in the folder:")
        for i in range(N):
            print(filelist[i])
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
    print('obs_times written to '+outpath+'/')
    print('WARNING: FORMATTING IS STILL SCREWED UP!')
    print('FIGURE OUT HOW TO FORMAT THOSE LINES IN A MORE HUMAN READABLE WAY')
    print('WHEN YOU HAVE INTERNET AGAIN.')
    #pdb.set_trace()





# pro read_data,path,filename,outname,oldharps=oldharps,blaze=blaze,ccf=startype,air=air
#                                 ;This program reads a HARPS timeseries
#                                 ;and reforms it order by order.
#                                 ;Specify the list of paths in a
#                                 ;textfile and use the path to that
#                                 ;file as input (with slash at the end). This program will then
#                                 ;sort out the orders from the provided
#                                 ;exposures, give the wavelength
#                                 ;solution per order and generate a
#                                 ;file with obs-time / exp-time
#                                 ;information.
#                                 ;Outname is the name of the folder in
#                                 ;./data in which the extracted orders
#                                 ;will be located. No slash at the end.
#
#                                 ;The pathfile needs 3 columns: The
#                                 ;first one is the filepath to the 2D
#                                 ;HARPS exposure. The second is the
#                                 ;type (either SCIENCE, WAVE or
#                                 ;BLAZE). The flag particular frame
#                                 ;needs to be discarded.
#
#                                 ;Set the ccf keyword to the name of
#                                 ;the template type that was used to
#                                 ;compute CCFs (e.g. K5), and then the program
#                                 ;will output the total order-averaged
#                                 ;CCF as well.
#
#                                 ;Set the blaze keyword to de-blaze the
#                                 ;orders using _blaze_A.fits files. The
#                                 ;/blaze keyword should have had to
#                                 ;have been set for this when running
#                                 ;prep_pathsfile.
#
#
#
#   readcol,path+filename,paths,type,flag,format='(A,A,A)'
#   i_w=where(type eq 'WAVE' and flag eq 1)
#   i_s=where(type eq 'SCIENCE' and flag eq 1)
#   if keyword_set(blaze) then begin
#      i_b=where(type eq 'BLAZE' and flag eq 1)
#      for ib=0,n_elements(i_b)-1,1 do begin
#         bf=readfits(path+paths[i_b[ib]]+'_blaze_A.fits',/silent)
#         if ib eq 0 then blazeframe=bf*1.0 else blazeframe+=bf
#      endfor
#      blazeframe=blazeframe/n_elements(i_b)
#   endif
#
#
#   n_s=n_elements(i_s)
#   outpath='data/'+outname
#   file_mkdir,outpath
#   openw,1,outpath+'/obs_times'
#   for i=0,n_s-1,1 do begin
#      exp=double(readfits(path+paths[i_s[i]]+'_e2ds_A.fits',h))
#      if keyword_set(startype) then ccf=double(readfits(path+paths[i_s[i]]+'_ccf_'+startype+'_A.fits',ccfheader))
#      if keyword_set(blaze) then exp=exp/blazeframe
#      if i eq 0 then begin
#         size=size(exp)          ;Size of the first exposure
#         n_orders=size[2]
#         n_pixels=size[1]
#                                 ;Wow, these are always the same for HARPS??
#         out_matrix=make_array(n_pixels,n_s,n_orders,/double)
#         if keyword_set(startype) then ccf_matrix=make_array(n_elements(ccf[*,0]),n_s,/double)
#      endif
#      median=median(exp,dimension=1)
#      m_array=rebin(reform(median,1,n_orders),n_pixels,n_orders)
#      out_matrix[*,i,*]=reform(exp,n_pixels,1,n_orders)
#      if keyword_set(startype) then ccf_matrix[*,i]=ccf[*,n_orders];The ccf has one extra line at the top which is the average.
#      date_obs=sxpar(h,'DATE-OBS')
#      mjd_obs= sxpar(h,'MJD-OBS')
#      exptime= sxpar(h,'EXPTIME')
#      if keyword_set(oldharps) then airmass=float(sxpar_eso(h,'HIERARCH ESO TEL AIRM START')) else airmass= sxpar(h,'AIRMASS')
#
#      printf,1,mjd_obs,'    '+date_obs,exptime,airmass,format='(D,A,I,F)'
#   endfor
#   close,1
#   wave_array=double(readfits(path+paths[i_w[0]]+'_wave_A.fits'))
#   if keyword_set(air) then begin
#      airtovac,wave_array,wave_array_vac
#      wave_array=wave_array_vac*1.0
#   endif
#
#   for i=0,n_orders-1 do begin
#      writefits,outpath+'/order_'+trim(i)+'.fits',out_matrix[*,*,i]
#      writefits,outpath+'/wave_'+trim(i)+'.fits',wave_array[*,i]/10.0;NANOMETERS
#   endfor
#   if keyword_set(startype) then begin
#      writefits,outpath+'/CCF.fits',ccf_matrix
#      void=readfits(outpath+'/CCF.fits',newccfheader)
#      sxaddpar,newccfheader,'CRVAL1',sxpar(ccfheader,'CRVAL1')
#      sxaddpar,newccfheader,'CTYPE1',sxpar(ccfheader,'CTYPE1')
#      sxaddpar,newccfheader,'CDELT1',sxpar(ccfheader,'CDELT1')
#      modfits,outpath+'/CCF.fits',0,newccfheader
#   endif
#
#   print,'Processed '+trim(n_orders)+' spectral orders in '+trim(n_s)+' science exposures.'
#
# end
