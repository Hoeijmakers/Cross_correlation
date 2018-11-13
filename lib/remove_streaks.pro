function remove_streaks,im
  im=local_outlier_remove(im,6.0,10)
  psf=psf_gaussian(npixel=[21,21],fwhm=6.0,/normalize)
  imG=convol(im,psf,/edge_truncate)
  s=stddev(imG)
  sel=where(abs(imG-mean(imG)) gt 4.0*s)
  imM=imG
  imM[sel]=!values.f_nan

  
  ;x_kernel=reform([-4,-3,-2,-1,0,1,2,3,4],9,1)
  ;y_kernel=transpose(x_kernel)
  ;nx=n_elements(im[*,0])
  ;ny=n_elements(im[0,*])

  ;xderiv=convol(imG,x_kernel,/center,/edge_truncate)
  ;yderiv=convol(imG,y_kernel,/center,/edge_truncate)

  ;angle=atan(yderiv/xderiv)

  writefits,'ohai.fits',[[[im]],[[imG]],[[imM]]]



  stop



end

