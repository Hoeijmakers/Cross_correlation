function local_outlier_remove,order,t,dpx,talk=talk,sigma=sigma
  npx=n_elements(order[*,0])
  order_clean=order*1.0
  x_axis=findgen(npx)
  sig=findgen(npx)
  m=findgen(npx)
  mask=order*0.0
  for i=0,npx-1,1 do begin
     if i-dpx ge 0 and i+dpx le npx-1 then sel=order[i-dpx:i+dpx,*] else begin
        if i-dpx lt 0 then sel=order[0:2*dpx,*]
        if i+dpx gt npx-1 then sel=order[npx-1-2*dpx:npx-1,*]
     endelse

     sig[i]=robust_sigma(sel)
     m[i]=median(sel)
     column=order[i,*]
     outliers=where(column gt m[i]+t*sig[i] or column lt m[i]-t*sig[i])
     ;print,outliers
     if outliers ne [-1] then begin
        order_clean[i,outliers]=!values.f_nan
        mask[i,outliers]=1.0
        ;print,i,n_elements(outliers)
     endif
  endfor
  if keyword_set(talk) then print,'n_outliers removed: '+trim(n_elements(where(mask eq 1.0)))
  sigma=sig*1.0
  return,order_clean
end

     
