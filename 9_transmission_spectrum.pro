pro transmission_spectrum,dp,n,wlrange=wlrange,model=model
  order=readfits(dp+'/order_'+trim(n)+'.fits')
  wave=readfits(dp+'/wave_'+trim(n)+'.fits')
  profile=calctimes(dp,/transit)
  rv=calctimes(dp)

  if keyword_set(model) then modelmatrix=get_model(model)

  skymodelmatrix=get_model('sky_norm')
  wls=skymodelmatrix[*,0]
  fxs=skymodelmatrix[*,1]

  
  binsize=15.0
  w=200.0
  thresh=5.0
  
  outsel=where(profile eq 0)
  insel=where(profile eq 1)

  npx=n_elements(order[*,0])
  nexp=n_elements(order[0,*])
  


  mean_out=mean(order[*,outsel],dimension=2)
  mean_in=mean(order[*,insel],dimension=2)

  order_norm=order/rebin(reform(mean_out,npx,1),npx,nexp)
  order_norm2=order_norm/rebin(reform(mean(order_norm,dim=1),1,nexp),npx,nexp)
                                ;for i=0,nexp-1,1 do
                                ;order_flat[*,i]=order_norm[*,i]/medsmooth(order_norm[*,i],w,/edge_truncate
                                ;I dont think flattening is required if the blaze is stable.
  
  outlier_mask=construct_outlier_mask(order_norm2,w,thresh,sigma=sigma_array)


  if where(outlier_mask eq 1) ne [-1] then print,n_elements(where(outlier_mask eq 1))
  order[outlier_mask]=!values.f_nan

  mean_out=mean(order[*,outsel],dimension=2,/nan)
  mean_in=mean(order[*,insel],dimension=2,/nan)

  tspec=mean_in/mean_out
  tspec/=mean(tspec)
  tspecs=smooth(tspec,binsize)
  cgplot,wave,tspec,/yno,xrange=[651.5,652],/xs
  cgplot,wave,tspecs,/overplot,color='red'
  cgplot,wls,fxs,/overplot,color='skyblue'

  telfit=ladfit(interpol(fxs,wls,wave),tspec)
  cgplot,wls,fxs*telfit[1]+telfit[0],/overplot,color='grn5'
  stop
end
