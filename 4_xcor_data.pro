pro xcor_data,dataset,n,modelname,ccv=ccv,dRV,RVrange,fast=fast,normalize_model=normalize_model,edgew=edgew
                                ;Set the edgew keyword to a number of
                                ;pixels at the edges of the order to
                                ;ignore. This is because the cleaning
                                ;program will screw up the edges a bit
                                ;over a range of width w, due to the
                                ;hipass filtering.
  Rm=0

  model=double(get_model(modelname,resolution=Rm))
  order=double(readfits('output/'+dataset+'/cleaning/order_'+trim(n)+'.fits',/silent))
  order_i=double(readfits('output/'+dataset+'/cleaning/order_'+trim(n)+'_injected.fits',/silent))
  wl=double(readfits('data/'+dataset+'/wave_'+trim(n)+'.fits',/silent))
  wlm=model[*,0]
  fxm=model[*,1]
  ;Blur the model to match the sampling of the data wl grid:
  Rs=mean(wl)*n_elements(wl)/(max(wl)-min(wl))/2.0 ;Sampling resolution
  sigma=resolution_kernel(Rs,Rm,mean(wl))
  sel=where(wlm ge min(wl) AND wlm le max(wl)) ;To get a local pixel scale
  pxscale=(max(wlm[sel])-min(wlm[sel]))/n_elements(wlm[sel]) ;model scale nm/px
  outpath='output/'+dataset+'/xcor/'+modelname
  
                                ;if sigma gt 0.0 then
                                ;fxm=gauss_smooth(fxm,sigma/pxscale,/edge_truncate)

  ;cgplot,wlm,(fxm-1.0)/10.0,/yno,xrange=[382.19,382.22]
  ;cgplot,wl,rebinw((fxm-1.0)/10.0,wlm,wl),/overplot,psym=1,color='blu4'
  ;cgplot,wl,interpol((fxm-1.0)/10.0,wlm,wl),/overplot,psym=1,color='hotpink'
  fxm=convol(fxm,psf_gaussian(npixel=sigma/pxscale*40.0,st_dev=sigma/pxscale,/normalize,ndimen=1),/edge_truncate)
  ;cgplot,wlm,(fxm-1.0),/overplot,color='red'
                                ;cgplot,wlm,interpolate_template(fxm-1.0,wlm,wl)

                                ;THIS CONVOLUTION BROADENS THE SIGNAL
                                ;SIGNIFICANTLY, REDUCING THE MEASURED
                                ;LINE DEPTHS BY ~10% if the planet
                                ;lines are not significantly
                                ;broadened. Flux should be conserved, though.

  ;cgplot,wl,interpol(fxm-1.0,wlm,wl),/overplot,color='red',psym=1
  ;cgplot,wl,interpol(fxm-1.0)
 


  print,">>> Cross-correlating order "+trim(n)
  
  if keyword_set(normalize_model) then begin

     continuum=polyfit_continuum(wlm,fxm,5.0,0.0,3.0)
     if max(continuum) eq 0 then begin
        print,'POLYFIT CONTINUUM DIDNT WORK. CRASHING.'
        stop
     endif
     
     fxm=fxm-continuum
     sel=where(wlm ge min(wl) and wlm le max(wl))
     dsel=where(abs(fxm-median(fxm)) lt 0.005*max(abs(fxm[sel])))
     
     fxm[where(fxm gt 0.0)] = 0.0
     fxm[dsel]=0.0
     ;cgplot,wlm[sel],fxm[sel]
     print,'>>> Normalize model set'
                                ;fxm=(-1.0)*fxm/total(fxm)

    
  endif
  
  ;wlm=400.0+findgen(20000)*6.0/10000.0
  ;wl=wlm[3000:7000]
  ;order=make_array(n_elements(wl),5)
  ;mask=findgen(n_elements(wlm))
  ;newdata=mask*0.0+1.0
  ;newdata[where(mask MOD 150 eq 0)] = 0.99
  ;for i=0,n_elements(order[0,*])-1,1 do order[*,i]=newdata[3000:7000]
  ;fxm=newdata*1.0-1.0


  
                                ;fxm=fxm/smooth(fxm,10000.0,/edge_truncate)
  profile=calctimes('data/'+dataset,/transit)
  
  if keyword_set(fast) then profile[where(profile gt 0.0)]=1.0 else profile=profile*0.0+1.0

  ;diff=order_i-order
  ;ccf_t=xcor2D(wl,diff,wlm,fxm,drv,RVrange,/nan,ccv=ccv,verticalmask=profile)
                                ;stop
  
  mask_trigger=paramget('mask','data/'+dataset)
  if mask_trigger eq 1 then begin
     mask=readfits('output/'+dataset+'/cleaning/order_'+trim(n)+'_mask.fits',/silent)
     masksel=where(mask eq 1)
     if masksel ne [-1] then begin
        order[masksel]=!values.f_nan
        order_i[masksel]=!values.f_nan
     endif
  endif
  
  if keyword_set(edgew) and edgew gt 0 then begin
     wl=wl[edgew:-edgew]
     order=order[edgew:-edgew,*]
     order_i=order_i[edgew:-edgew,*]
  endif
  
  file_mkdir,outpath
  ccf=xcor2D(wl,order,wlm,fxm,drv,RVrange,/nan,ccv=ccv,verticalmask=profile)
  ccf_i=xcor2D(wl,order_i,wlm,fxm,drv,RVrange,/nan,ccv=ccv,verticalmask=profile)
  writefits,outpath+'/xcor_template_order_'+trim(n)+'.fits',ccf.xcor_template

  if n_elements(where(finite(ccf.ccf_matrix))) ne n_elements(ccf.ccf_matrix) then begin
     print,'>>> Too many NaNs in CCF. Do not use order '+trim(n)
  endif
  


  

  if n_elements(ccf.template_sum) gt 0 and total(ccf.template_sum) gt 0 then begin
     print,'Template > 0'
  endif
  
  write_ccf_matrix,ccf.ccf_matrix,ccf.RV,outpath+'/ccf_order_'+trim(n),ccf.template_sum
  write_ccf_matrix,ccf_i.ccf_matrix,ccf_i.RV,outpath+'/ccf_order_'+trim(n)+'_injected',0.0
  
  ;print,ccf.template_sum

end


