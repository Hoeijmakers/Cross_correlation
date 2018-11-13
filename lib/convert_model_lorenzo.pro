pro convert_model_lorenzo,path,outpath

  readcol,path,wl,R,format='(D,D)'

  fx=1.0-R

  writefits,outpath,[[wl/10.0],[fx]]
  cgplot,wl,fx,/yno
end
