function xcor_rv,wld,fxd,wlm,fxm,drv,rvrange,ccv=ccv
  if keyword_Set(ccv) then yes=1 else yes=0
  c=299792d0                    ;km/s
  v=dindgen(2*RVrange/dRV+1)*dRV-RVrange
  shift=1.00+v/c
  ccf=v*0.0
  for vi=0,n_elements(v)-1,1 do begin
     WLMS=wlm*shift[vi]
     fxm_i=interpol(fxm,WLMS,wld)
     ccf[vi]=c_correlate(fxd,fxm_i,0.0,covariance=yes)     
  endfor

  return,[[v],[ccf]]
end
