pro stack_ccfs,dataset,modelname,orders,injected=injected,ccv=ccv ;,sigma=sigma,signal=signal
                                ;Set the normal_model keyword to
                                ;normalize the CCV by the total sum of
                                ;all chuncks of the template. Thks
                                ;needs to be done if the CCV is
                                ;computed by a mask that does not
                                ;integrate to 1.0 over the bandpass of
                                ;the data. I think this is the case
                                ;most of the time.
  
                                ;Set the equal norm keyword to
                                ;normalize the CCV by the sum of the
                                ;template in each chunk. This always
                                ;needs to be done for the CCV to make sense.
  no=n_elements(orders)
  dp='data/'+dataset
  vi=calctimes(dp)
  profile=calctimes(dp,/transit)
  inpath='output/'+dataset+'/xcor/'+modelname+'/'
  original_xcor_template=get_model(modelname)
  outpath='output/'+dataset+'/stacking/'+modelname+'/'
  file_mkdir,outpath

  orderid=[]
  wlt=[]                        ;This is to plot the template in each order.
  fxt=[]
  for i=0,no-1,1 do begin
     xcortemplate=readfits(inpath+'xcor_template_order_'+trim(orders[i])+'.fits',/silent)
     wlt=[wlt,xcortemplate[*,0]]
     fxt=[fxt,xcortemplate[*,1]]
     orderid=[orderid,xcortemplate[*,0]*0.0+orders[i]]
     
  endfor

  minwlt=min(wlt)
  maxwlt=max(wlt)
  minfxt=min(fxt)
  maxfxt=max(fxt)
  


  
  for i=0,no-1,1 do begin
     if keyword_set(injected) then ccf=readfits(inpath+'ccf_order_'+trim(orders[i])+'_injected.fits',h,/silent) else ccf=readfits(inpath+'ccf_order_'+trim(orders[i])+'.fits',h,/silent)

     ccf_unsk=ccf*0.0
     if i eq 0 then begin
        nRV=n_elements(ccf[*,0])
        nexp=n_elements(ccf[0,*])
        ccf_out=ccf*0.0
        ccf_stack=make_array(nRV,nexp,no)
        ccf_stack_unsk=ccf_stack*0.0
        RV=sxpar(h,'CRVAL2')+findgen(nRV)*sxpar(h,'CDELT2')
        minRV=min(RV)
        maxRV=max(RV)
        tsum_matrix=make_array(nRV,no)
        cgps_open,'output/'+dataset+'/stacking/'+modelname+'/xcor_template.ps'
        cgdisplay,1200,500
        cgplot,original_xcor_template[*,0],original_xcor_template[*,1],xrange=[minwlt,maxwlt],yrange=[minfxt,maxfxt],charsize=1.0,xtitle='Wavelength',ytitle='Weight',thick=0.5
       
     endif
     xcortemplate=readfits(inpath+'xcor_template_order_'+trim(orders[i])+'.fits',/silent)
     wlt=xcortemplate[*,0]
     fxt=xcortemplate[*,1]
     cgplot,wlt,fxt,/overplot,color='red',thick=0.5
     cgvline,max(wlt),color='grn4',thick=3
     cgtext,mean(wlt),plotpos(0.94,/y),trim(orders[i]),align=0.5,color='grn4',charsize=1
     
     if keyword_set(ccv) then tsum_matrix[*,i]=readfits('output/'+dataset+'/xcor/'+modelname+'/ccf_order_'+trim(orders[i])+'_tsum.fits',/silent)

     ccf_stack[*,*,i]=ccf*1.0
     for j=0,nexp-1 do begin
        ccf_unsk[*,j]=interpol(ccf[*,j],RV-vi[j],RV)
        sel=where(RV+vi[j] lt min(RV) or RV+vi[j] gt max(RV))
        ccf_unsk[sel,j]=!values.f_nan
     endfor
     ccf_stack_unsk[*,*,i]=ccf_unsk*1.0
  endfor
  cgplot,0,0,/noerase,xrange=[minwlt,maxwlt],yrange=[minfxt,maxfxt],charsize=1.0
  cgps_close



  if n_elements(abs(tsum_matrix)) eq nRV then tsum_RV=abs(tsum_matrix)*1.0 else tsum_RV=total(abs(tsum_matrix),2)

  
;  if keyword_set(ccv) then begin
;     for i=0,no-1,1 do begin
;        for j=0,nexp-1,1 do begin
           
;           ccf_stack[*,j,i]=ccf_stack[*,j,i]*abs(tsum_matrix[*,i])/tsum_RV
;           ccf_stack_unsk[*,j,i]=ccf_stack_unsk[*,j,i]*abs(tsum_matrix[*,i])/tsum_RV
;        endfor
;     endfor
;  endif
;This doesn't seem to work because below you multiply again
;with a weight, which will cause the areas of the template to no
;longer add up to 1.0. I think that's the problem here... Jens Sep 12 2018

  
                                ;So after this the value of each CCV
                                ;is normalized to the total area under
                                ;the template. This means that to
                                ;complete the CCV of the entire
                                ;dataset, the individual CCV values
                                ;will need to be summed. This means
                                ;that the order-weights should not be
                                ;normalized to their sum, but to their
                                ;average.
                                ;Make sure that this happens in
                                ;co_add_ccfs.
                                ;Jens, before Sep 12 2018...
  

 
  cgps_open,outpath+'Tsum.ps'
  cgplot,orders,tsum_matrix[0.5*nRV,*],xtitle='Order number',ytitle='Sum(Template)',charsize=0.9
  cgps_close
  if keyword_set(injected) then suff='_injected' else suff=''
  writefits,outpath+'/ccf_stack'+suff+'.fits',ccf_stack
  writefits,outpath+'/ccf_stack_unsk'+suff+'.fits',ccf_stack_unsk
  writefits,outpath+'/RV.fits',RV

  ;writefits,'test.fits',total(ccf_stack_unsk,3)
  
  
end
