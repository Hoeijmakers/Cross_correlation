pro read_data,path,filename,outname,oldharps=oldharps,blaze=blaze,ccf=startype,air=air
                                ;This program reads a HARPS timeseries
                                ;and reforms it order by order.
                                ;Specify the list of paths in a
                                ;textfile and use the path to that
                                ;file as input (with slash at the end). This program will then
                                ;sort out the orders from the provided
                                ;exposures, give the wavelength
                                ;solution per order and generate a
                                ;file with obs-time / exp-time
                                ;information.
                                ;Outname is the name of the folder in
                                ;./data in which the extracted orders
                                ;will be located. No slash at the end.
  
                                ;The pathfile needs 3 columns: The
                                ;first one is the filepath to the 2D
                                ;HARPS exposure. The second is the
                                ;type (either SCIENCE, WAVE or
                                ;BLAZE). The flag particular frame
                                ;needs to be discarded.
  
                                ;Set the ccf keyword to the name of
                                ;the template type that was used to
                                ;compute CCFs (e.g. K5), and then the program
                                ;will output the total order-averaged
                                ;CCF as well.

                                ;Set the blaze keyword to de-blaze the
                                ;orders using _blaze_A.fits files. The
                                ;/blaze keyword should have had to
                                ;have been set for this when running
                                ;prep_pathsfile.
  
  

  readcol,path+filename,paths,type,flag,format='(A,A,A)'
  i_w=where(type eq 'WAVE' and flag eq 1)
  i_s=where(type eq 'SCIENCE' and flag eq 1)
  if keyword_set(blaze) then begin
     i_b=where(type eq 'BLAZE' and flag eq 1)
     for ib=0,n_elements(i_b)-1,1 do begin
        bf=readfits(path+paths[i_b[ib]]+'_blaze_A.fits',/silent)
        if ib eq 0 then blazeframe=bf*1.0 else blazeframe+=bf
     endfor
     blazeframe=blazeframe/n_elements(i_b)
  endif
  
  
  n_s=n_elements(i_s)
  outpath='data/'+outname
  file_mkdir,outpath
  openw,1,outpath+'/obs_times'
  for i=0,n_s-1,1 do begin
     exp=double(readfits(path+paths[i_s[i]]+'_e2ds_A.fits',h))
     if keyword_set(startype) then ccf=double(readfits(path+paths[i_s[i]]+'_ccf_'+startype+'_A.fits',ccfheader))
     if keyword_set(blaze) then exp=exp/blazeframe
     if i eq 0 then begin
        size=size(exp)          ;Size of the first exposure
        n_orders=size[2]
        n_pixels=size[1]
                                ;Wow, these are always the same for HARPS??
        out_matrix=make_array(n_pixels,n_s,n_orders,/double)
        if keyword_set(startype) then ccf_matrix=make_array(n_elements(ccf[*,0]),n_s,/double)
     endif
     median=median(exp,dimension=1)
     m_array=rebin(reform(median,1,n_orders),n_pixels,n_orders)
     out_matrix[*,i,*]=reform(exp,n_pixels,1,n_orders)
     if keyword_set(startype) then ccf_matrix[*,i]=ccf[*,n_orders];The ccf has one extra line at the top which is the average.
     date_obs=sxpar(h,'DATE-OBS')
     mjd_obs= sxpar(h,'MJD-OBS')
     exptime= sxpar(h,'EXPTIME')
     if keyword_set(oldharps) then airmass=float(sxpar_eso(h,'HIERARCH ESO TEL AIRM START')) else airmass= sxpar(h,'AIRMASS')

     printf,1,mjd_obs,'    '+date_obs,exptime,airmass,format='(D,A,I,F)'
  endfor
  close,1
  wave_array=double(readfits(path+paths[i_w[0]]+'_wave_A.fits'))
  if keyword_set(air) then begin
     airtovac,wave_array,wave_array_vac
     wave_array=wave_array_vac*1.0
  endif
  
  for i=0,n_orders-1 do begin
     writefits,outpath+'/order_'+trim(i)+'.fits',out_matrix[*,*,i]
     writefits,outpath+'/wave_'+trim(i)+'.fits',wave_array[*,i]/10.0;NANOMETERS
  endfor
  if keyword_set(startype) then begin
     writefits,outpath+'/CCF.fits',ccf_matrix
     void=readfits(outpath+'/CCF.fits',newccfheader)
     sxaddpar,newccfheader,'CRVAL1',sxpar(ccfheader,'CRVAL1')
     sxaddpar,newccfheader,'CTYPE1',sxpar(ccfheader,'CTYPE1')
     sxaddpar,newccfheader,'CDELT1',sxpar(ccfheader,'CDELT1')
     modfits,outpath+'/CCF.fits',0,newccfheader
  endif
  
  print,'Processed '+trim(n_orders)+' spectral orders in '+trim(n_s)+' science exposures.'
  
end

pro read_express,path,filename,outname,air=air
  ;This is for the EXPRESS MASCARA-2 data.
  readcol,path+filename,paths,type,target,format='(A,A,A)'
  n=n_elements(paths)
  outpath='data/'+outname
  file_mkdir,outpath
  openw,1,outpath+'/obs_times'
  for i=0,n-1,1 do begin
     name=path+paths[i]
     exp=double(readfits(name,h,/silent))   
     date_obs=sxpar(h,'DATE-OBS')
     mjd_obs=sxpar(h,'TELMJD')
     exptime=sxpar(h,'REXPTIME')
     airmass=sxpar(h,'AIRMASS')
     bzero=sxpar(h,'O_BZERO')
     exp-=bzero
     
    printf,1,mjd_obs,'     '+date_obs,exptime,airmass,format='(D,A,I,F)'
     if i eq 0 then begin
        n_orders=n_elements(exp[0,0,*])
        n_px=n_elements(exp[0,*,0])
        order_stack=make_array(n_px,n,n_orders)
        wave=reform(exp[0,*,*])
        order_stack[*,i,*]=exp[1,*,*]
     endif else begin
        wave_new=reform(exp[0,*,*])
        diff = abs(wave_new-wave)
        if max(diff) ne 0 then begin
           print,max(diff)
           print,'interpolating exp '+trim(i)
           for j=0,n_orders-1,1 do order_stack[*,i,j]=interpol(exp[1,*,j],wave_new[*,j],wave[*,j])
        endif else order_stack[*,i,*]=exp[1,*,*]
     endelse
     

  endfor
  close,1

  if keyword_set(air) then begin
     airtovac,wave,wave_vac
     wave=wave_vac*1.0
  endif
  
  for i=0,n_orders-1 do begin
     writefits,outpath+'/order_'+trim(i)+'.fits',order_stack[*,*,i]
     writefits,outpath+'/wave_'+trim(i)+'.fits',wave[*,i]/10.0;NANOMETERS
  endfor
     
     

  stop    

end

pro read_express_new,path,filename,outname,air=air
                                ;This is for the EXPRESS MASCARA-2
                                ;data, with the new datareduction
                                ;that contains the wl solution in each
                                ;frame as a cube.
  readcol,path+filename,paths,type,target,format='(A,A,A)'
  n=n_elements(paths)
  outpath='data/'+outname
  file_mkdir,outpath
  openw,1,outpath+'/obs_times'
  for i=0,n-1,1 do begin
     name=path+paths[i]
     exp=double(readfits(name,h,/silent))   
     date_obs=sxpar(h,'DATE-OBS')
     mjd_obs=sxpar(h,'TELMJD')
     exptime=sxpar(h,'REXPTIME')
     airmass=sxpar(h,'AIRMASS')
     
    printf,1,mjd_obs,'     '+date_obs,exptime,airmass,format='(D,A,I,F)'
     if i eq 0 then begin
        n_orders=n_elements(exp[0,*,0])
        n_px=n_elements(exp[*,0,0])
        order_stack=make_array(n_px,n,n_orders)
        wave=reform(exp[*,*,1])
        order_stack[*,i,*]=exp[*,*,0]
     endif else begin
        wave_new=reform(exp[*,*,1])
        diff = abs(wave_new-wave)
        if max(diff) ne 0 then begin
           print,max(diff)
           print,'interpolating exp '+trim(i)
           for j=0,n_orders-1,1 do order_stack[*,i,j]=interpol(exp[*,j,0],wave_new[*,j],wave[*,j])
        endif else order_stack[*,i,*]=exp[*,*,0]
     endelse
     

  endfor
  close,1

  if keyword_set(air) then begin
     airtovac,wave,wave_vac
     wave=wave_vac*1.0
  endif
  
  for i=0,n_orders-1 do begin
     writefits,outpath+'/order_'+trim(i)+'.fits',order_stack[*,*,i]
     writefits,outpath+'/wave_'+trim(i)+'.fits',wave[*,i]/10.0;NANOMETERS
  endfor
     
    
  stop    

end







pro read_1ds,path,filename,outname,chop=chop,air=air,custom=custom
                                ;This program reads the 1d spectra
                                ;that Lorenzo can use and outputs them
                                ;as if it is a single echelle order.
  readcol,path+filename,paths,type,target,flag,format='(A,A,A,A)'
  i_s=where(type eq 'SCIENCE' and flag eq 1)

  n_s=n_elements(i_s)
  outpath='data/'+outname
  file_mkdir,outpath
  openw,1,outpath+'/obs_times'
  for i=0,n_s-1,1 do begin
     name=path+paths[i_s[i]]
     if keyword_set(suffix) then name+=suffix else name+='_s1d_A.fits'
     exp=double(readfits(name,h))
     date_obs=sxpar(h,'DATE-OBS')
     mjd_obs= sxpar(h,'MJD-OBS')
     exptime= sxpar(h,'EXPTIME')
     airmass= sxpar(h,'AIRMASS')

     
     size=size(exp)       
     n_pixels=size[1]
     wld=findgen(n_pixels)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')


                                ;Uncomment the following 5 lines to
                                ;read Romains telluric corrected files instead.
     exp=mrdfits(path+paths[i_s[i]]+'_s1d_A_correct_telluric.fits',1,h)
     wld=exp.wavelength
     fx=exp.calc_flux
     exp=fx*1.0
     n_pixels=n_elements(exp)
                                ;until here.

     
     
     
     if i eq 0 then begin
        out_matrix=make_array(n_pixels,n_s,/double)
        ;n_px_master=n_pixels*1.0
        wld_master=wld*1.0
     endif
     
     ;if n_pixels gt n_px_master then n_pixels=n_px_master*1.0


     
     out_matrix[*,i]=interpol(exp,wld,wld_master)
     printf,1,mjd_obs,'    '+date_obs,exptime,airmass,format='(D,A,I,F)'
  endfor

  close,1
  
                                ;wave_array=double(readfits(path+paths[i_w[0]]+'_wave_A.fits'))

  if keyword_set(air) then begin
     wl_air=wld_master*1.0
     airtovac,wl_air,wld_master
  endif
  
     writefits,outpath+'/order_0.fits',out_matrix
     writefits,outpath+'/wave_0.fits',wld_master/10.0;NANOMETERS
     print,'Processed the 1d spectrum in '+trim(n_s)+' science exposures.'

     
     if keyword_set(chop) then begin
        wl=wld_master/10.0
        bin_min=min(wl)
        i=0
        while bin_min lt max(wl) do begin
           chunk=where(wl ge bin_min and wl lt bin_min+chop)
           writefits,outpath+'/order_'+trim(i)+'.fits',out_matrix[chunk,*]
           writefits,outpath+'/wave_'+trim(i)+'.fits',wl[chunk]
           bin_min+=chop
           i+=1
        endwhile
        print,trim(i)+' chops of '+trim(chop)+'nm produced.'
     endif
     
  
end



pro prep_pathsfile,path_to_datafolder,night_folder,night_number,harpsn=harpsn,blaze=blaze,nowave=nowave
;This code generates a pathfile in the path_to_datafolder, which
;contains one or more folders corresponding to one or more nights of
;observations. The night number is an integer/float that names the
;number of the night. This may only work for HARPS data, and may only
;work for older HARPS data (because of the way the obstype and target
;are encoded in the FITS header.

  ;Set the harpsn keyword if you are dealing with TNG keywords.

                                ;At the end of the program, you
                                ;manually need to edit the file such
                                ;that:
                                ;1)  '_e2ds_A.fits' is removed from each path.
                                ;2)  'the leading ../folder/' is
                                ;removed from each path. The pathsfile
                                ;should be in the same folder as the
                                ;datafolder afterall.
                                ;3)  Put a flag number (0 or 1) after
                                ;each science frame, depending on
                                ;whether you want the exposure to be
                                ;included (the OBJECT column comes in
                                ;handy here because there my be
                                ;multiple targets in each night).
                                ;4) Remove the OBJECT column.
                                ;5) Manually remove all WAVE frames
                                ;(rename the type to WAVE if
                                ;necessary) and set its flag to 1.


  spawn,'ls '+path_to_datafolder+'/'+night_folder+'/*_e2ds_A.fits',file_list
n=n_elements(file_list)
  if keyword_set(nowave) ne 1 then begin
     spawn,'ls '+path_to_datafolder+'/'+night_folder+'/*_wave_A.fits',wave_list  
     nw=n_elements(wave_list)
  endif
  

  if keyword_set(blaze) then begin
     spawn,'ls '+path_to_datafolder+'/'+night_folder+'/*_blaze_A.fits',blaze_list
     nb=n_elements(blaze_list)
  endif
  

  openw,1,path_to_datafolder+'/paths_night'+trim(night_number)
  for i=0,n-1,1 do begin
     file=readfits(file_list[i],h,/silent)
     if keyword_set(harpsn) then begin
        obstype = sxpar_eso(h,'HIERARCH TNG DPR CATG')
        target=sxpar_eso(h,'HIERARCH TNG OBS TARG NAME')        
     endif else begin
        obstype = sxpar_eso(h,'HIERARCH ESO DPR CATG')
        target=sxpar(h,'OBJECT')         
     endelse
     printf,1,file_list[i],' ',obstype,' ',target,format='(A,A,A,A,A)'       
  endfor

  if keyword_set(nowave) ne 1 then begin 
     for i=0,nw-1,1 do begin
        file=readfits(wave_list[i],h,/silent)
        if keyword_set(harpsn) then begin
           obstype = sxpar_eso(h,'HIERARCH TNG DPR CATG')
           target=sxpar_eso(h,'HIERARCH TNG OBS TARG NAME')        
        endif else begin
           obstype = sxpar_eso(h,'HIERARCH ESO DPR CATG')
           target=sxpar(h,'OBJECT')         
        endelse
        printf,1,wave_list[i],' ',obstype,' ',target,format='(A,A,A,A,A)'       
     endfor
  endif
  

  if keyword_set(blaze) then begin
     for i=0,nb-1,1 do begin
        file=readfits(blaze_list[i],h,/silent)
        if keyword_set(harpsn) then begin
           obstype = sxpar_eso(h,'HIERARCH TNG DPR CATG')
           target=sxpar_eso(h,'HIERARCH TNG OBS TARG NAME')        
        endif else begin
           obstype = sxpar_eso(h,'HIERARCH ESO DPR CATG')
           target=sxpar(h,'OBJECT')         
        endelse
        printf,1,blaze_list[i],' ',obstype,' ',target,format='(A,A,A,A,A)'       
     endfor
  endif
  





  
  close,1
end

  




pro check_wavelength_solution,dp
                                ;This program checks whether the
                                ;wavelength solution of a certain
                                ;dataset/night is in air or vaccuum by
                                ;overplotting a telluric and stellar
                                ;model to see by eye whether the lines
                                ;are shifted.

;pro view_orders,dataset


  skymodel=get_model('sky_norm')
  starmodel=get_model('phoenix-6600')

  wlsky=skymodel[*,0]
  fxsky=skymodel[*,1]

  wlstar=starmodel[*,0]
  fxstar=starmodel[*,1]

  
  a=1
  i=0
  rvs=[0,0,0]
  while a eq 1 do begin
     ordername=dp+'/order_'+trim(i)+'.fits'
     wavename=dp+'/wave_'+trim(i)+'.fits'
     void=file_test(ordername)
     print,ordername
     if void ne 1 then a=0 else begin
        wave=readfits(wavename,/silent)
        spec_raw=mean(readfits(ordername,/silent),dimension=2,/nan)
        cont=ladfit_continuum(wave,spec_raw,1.0,0.9)

        spec=spec_raw/cont
        cgplot,wave,spec


        skymodel_sel=where(wlsky gt min(wave)-1.0 and wlsky lt max(wave)-1.0)
        starmodel_sel=where(wlstar gt min(wave)-1.0 and wlstar lt max(wave)-1.0)
        cgplot,wlsky,fxsky,/overplot,color='blu5'

        wlstarn=wlstar[starmodel_sel]
        fxstarn=gauss_smooth(fxstar[starmodel_sel]/ladfit_continuum(wlstar[starmodel_sel],fxstar[starmodel_sel],1.0,0.9),50.0,/edge_truncate)
        
        cgplot,wlstarn,fxstarn,/overplot,color='grn5'
        cgtext,plotpos(0.95,/x),plotpos(0.92,/y),'PHOENIX model (vaccuum)',color='grn5',/align
        cgtext,plotpos(0.95,/x),plotpos(0.88,/y),'PHOENIX model (air)',/align,color='red4'       
        cgtext,plotpos(0.95,/x),plotpos(0.84,/y),'Data',/align
        cgtext,plotpos(0.95,/x),plotpos(0.80,/y),'Sky model (vaccuum)',color='blu4',/align
        vactoair,wlstarn*10,wlstarn_air
        wlstarn_air/=10.0
        cgplot,wlstarn_air,fxstarn,/overplot,color='red3'

        ;c_star=xcor_rv(wave,spec,wlstarn,smooth(fxstarn,10.0),1.0,150.0)
        ;void=max(c_star[*,1],loc)
        ;;peakrv=c_star[loc,0]
        ;cgplot,c_star[*,0],c_star[*,1]

        ;if skymodel_sel ne [-1] and min(fxsky[skymodel_sel]) lt 0.95 then begin
        ;   c_sky=xcor_rv(wave,spec,wlsky[skymodel_sel],fxsky[skymodel_sel],1.0,150.0)
        ;   void=max(c_sky[*,1],loc)
        ;   peakrv=c_sky[loc,0]
        ;   cgplot,c_sky[*,0],c_sky[*,1]

           ;rvs=[[rvs],[i,0,peakrv]]           
        ;endif
        stop
       
     endelse
     i+=1
  endwhile

end

        
        
