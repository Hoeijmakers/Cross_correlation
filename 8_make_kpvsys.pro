pro make_kpvsys,dataset,modelname,clean=clean

  inpath='output/'+dataset+'/stacking/'+modelname+'/'

  if keyword_set(clean) then begin
     cc_c=readfits(inpath+'CCF_stacked_2D_c_clean.fits',/silent)  
     cc_i=readfits(inpath+'CCF_stacked_2D_i_clean.fits',/silent)
  endif else begin
     cc_c=readfits(inpath+'CCF_stacked_2D_c.fits',/silent)  
     cc_i=readfits(inpath+'CCF_stacked_2D_i.fits',/silent)
  endelse
  
  n_exp=n_elements(cc_c[0,*])
  n_rv=n_elements(cc_c[*,0])
  ccf_v=readfits(inpath+'/RV.fits',/silent)
  exp_weights=readfits(inpath+'/weights_exp.fits',/silent)

  if n_rv ne n_elements(ccf_v) or n_exp ne n_elements(exp_weights) then begin
     print,'ERROR in make_kpvsys: Size of RV.fits or weights_exp.fits does not match the dimensions of the CCF.'
     print,'Stopping and returning.'
     stop
     return
  endif

  
  weights_matrix=rebin(reform(exp_weights,1,n_exp),n_rv,n_exp) 

                                ;=======Then unskew the
                                ;cross_correlation for different RV
                                ;amplitudes, while weighing=======

  kp=findgen(400)           ;km/s
  phase=calctimes('data/'+dataset,/alpha)
  kpvsys_c=make_array(n_elements(ccf_v),n_elements(kp))
  kpvsys_i=kpvsys_c*0.0


  for i=0,n_elements(kp)-1,1 do begin
     v_test=kp[i]*sin(phase*2d*!pi)
     cc_c_unskewed=ccf_unskew(cc_c,ccf_v,v_test);*weights_matrix
     cc_i_unskewed=ccf_unskew(cc_i,ccf_v,v_test);*weights_matrix   
     kpvsys_c[*,i]=mean(cc_c_unskewed[*,20:37],dim=2)
     kpvsys_i[*,i]=mean(cc_i_unskewed[*,20:37],dim=2)
  endfor

  

  S=max(kpvsys_i-kpvsys_c,loc)
  maxcoords=array_indices(kpvsys_c,loc)

  noisesel=where(abs(ccf_v) gt 150.0)
  N=stddev(kpvsys_c[noisesel,maxcoords[1]],/nan)


  
  if keyword_set(clean) then begin  
     writefits,'output/'+dataset+'/stacking/'+modelname+'/KpVsys_c_clean.fits',kpvsys_c
     writefits,'output/'+dataset+'/stacking/'+modelname+'/KpVsys_i_clean.fits',kpvsys_i     
     writefits,'output/'+dataset+'/stacking/'+modelname+'/Kp.fits',kp
     writefits,'output/'+dataset+'/stacking/'+modelname+'/KpVsys_diff_clean.fits',kpvsys_i-kpvsys_c
     print,'Injected peak SNR after doppler-model cleaning: '+r_trim(S/N,2)
  endif else begin  
     writefits,'output/'+dataset+'/stacking/'+modelname+'/KpVsys_c.fits',kpvsys_c
     writefits,'output/'+dataset+'/stacking/'+modelname+'/KpVsys_i.fits',kpvsys_i     
     writefits,'output/'+dataset+'/stacking/'+modelname+'/Kp.fits',kp
     writefits,'output/'+dataset+'/stacking/'+modelname+'/KpVsys_diff.fits',kpvsys_i-kpvsys_c
     print,'Injected peak SNR without doppler-model cleaning: '+r_trim(S/N,2)     
  endelse
  

  ;All this is for plotting the Kelt 9 result:
  ;xsel=where(ccf_v gt -500 and ccf_v lt 500)
  ;ysel=where(kp gt 150)
  ;xtitle='Systemic velocity (km/s)'
  ;ytitle='Orbital velocity (km/s)'
  ;cs=1.0
  ;loadct,3
  
  ;cgps_open,'KpVsys_Fe_I_3000K.ps'
  ;cgimage,kpvsys[min(xsel):max(xsel),min(ysel):max(ysel)],/axes,xrange=[min(ccf_v[xsel]),max(ccf_v[xsel])],yrange=[min(kp[ysel]),max(kp[ysel])],charsize=cs,xtitle=xtitle,ytitle=ytitle
  ;cgtext,plotpos(0.98,/x),plotpos(0.92,/y),'K$\downp$V$\downsys$ diagram of Fe I at 3000 K',/align,color='white',charsize=cs*1.4,charthick=4.0  
  ;cgps_close
  
end



