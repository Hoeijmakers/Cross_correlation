function xcor2D,wld,order,wlm,fxm,dRV,RVrange,nan=nan,ccv=ccv,verticalmask=verticalmask
  if keyword_set(ccv) then covariance=1 else covariance=0
  c=299792d0                    ;km/s
  v=dindgen(2*RVrange/dRV+1)*dRV-RVrange
  shift=1.00+v/c
  
  size=size(order)
  nexp=size[2]
  
  CCF_out=make_array(n_elements(v),nexp)*0.0

  finite_mask=order*0.0
  finite_mask[where(finite(order))]=1.0
  if keyword_set(verticalmask) ne 1 then verticalmask=findgen(nexp)*0.0+1.0
  tsum=[]

  ;order_n=order*1.0
  ;for i=0,nexp-1,1 do order_n[*,i]=order[*,i]-mean(order[*,i],/nan)
  
  if keyword_set(nan) then begin
     for vi=0,n_elements(v)-1,1 do begin
        
        WLMS=wlm*shift[vi]
        fxm_i=interpol(fxm,WLMS,wld)
        tsum=[tsum,total(fxm_i)]
        if covariance eq 1 then fxm_i=fxm_i*n_elements(wld)
        
        for i=0,nexp-1,1 do begin
           if verticalmask[i] eq 1 then begin
              sel=where(finite(order[*,i]))
              CCF_out[vi,i]=c_correlate(fxm_i[sel],order[sel,i],0.0,covariance=covariance)
           endif
           
        endfor
        
     endfor
  endif else begin
       for vi=0,n_elements(v)-1,1 do begin
          WLMS=wlm*shift[vi]
          fxm_i=interpol(fxm,WLMS,wld)
          tsum=[tsum,total(fxm_i)]
          if covariance eq 1 then fxm_i=fxm_i*n_elements(wld) ;HERE IT IS!! See log page after april 16 on derivation. This is correct dige.

          
          for i=0,nexp-1,1 do if verticalmask[i] eq 1 then CCF_out[vi,i]=c_correlate(fxm_i,order[*,i],0.0,covariance=covariance);if covariance eq 1 then CCF_out[vi,i]=total(
       endfor
    endelse
  

  if covariance eq 1 then template_sum=mean(tsum) else template_sum=0.0
  st_out=create_struct('CCF_matrix',CCF_out,'RV',v,'template_sum',template_sum)
  return,st_out  
end
