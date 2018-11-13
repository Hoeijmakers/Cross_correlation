pro clean_ccfs,dataset,modelname,doppler_model_path
                                ;This program takes the 2D CCFs and
                                ;detrends them with some form of
                                ;Doppler shadow model that is defined
                                ;in a textfile as a gaussian fit for
                                ;each exposure.




  tic
  vsys=paramget('vsys','data/'+dataset)
  prv=calctimes('data/'+dataset)+vsys
  inpath='output/'+dataset+'/stacking/'+modelname+'/'
  ;ccf_2D_c=readfits(inpath+'CCF_stacked_2D_cw.fits')  
  ;ccf_2D_i=readfits(inpath+'CCF_stacked_2D_iw.fits')
  CCF_2D_c=readfits(inpath+'CCF_stacked_2D_c.fits',/silent)
  CCF_2D_i=readfits(inpath+'CCF_stacked_2D_i.fits',/silent)

  
  ;weights=total(readfits(inpath+'weights_matrix.fits'),3)
  
  RV=readfits(inpath+'RV.fits',/silent)
  n_exp=n_elements(ccf_2D_c[0,*])

  ;readcol,doppler_model_path,trig,amp,pos,sig,format='(I,F,F,F)'
  readcol,doppler_model_path,trig,amp,pos,sig,amp2,pos2,sig2,amp3,pos3,sig3,amp4,pos4,sig4,amp5,pos5,sig5,format='(I,F,F,F,F,F,F,F,F,F,F,F,F)'
  doppler_model=ccf_2D_c*0.0

  ccf_masked=ccf_2D_c*1.0
  m=12.0
  tsel=where(trig eq 1)
  doppler_model_masked=doppler_model*1.0 
  for i=0,n_exp-1,1 do begin
     params=[amp[i],pos[i],sig[i]]
     params2=[amp2[i],pos2[i],sig2[i]]
     params3=[amp3[i],pos3[i],sig3[i]]
     params4=[amp4[i],pos4[i],sig4[i]]
     params5=[amp5[i],pos5[i],sig5[i]]     
     if trig[i] eq 1 then begin
        doppler_model[*,i]=gaussian(RV,params)  +gaussian(RV,params2)+gaussian(RV,params3)+gaussian(RV,params4)+gaussian(RV,params5)
        outsel=where(abs(doppler_model[*,i]) lt max(abs(doppler_model[*,i]))/1d7) ;Select all parts of the doppler_model that are infinitesimally small. I don't want them to contribute.
        doppler_model_masked[*,i]=doppler_model[*,i]*1.0
        doppler_model_masked[outsel,i]=!values.f_nan
        

        planetsel=where(RV gt prv[i]-m and RV lt prv[i]+m)
        ccf_masked[planetsel,i]=!values.f_nan
     endif
     
  endfor

  doppler_model_masked[*,where(trig eq 0.0)]=!values.f_nan
  
                                ;doppler_model=doppler_model*weights


                                ;The following 6 lines make no sense?
                                ;In the case of a vertical doppler
                                ;model I would mask it all out??
                                ;Chichi?
                                ;I moved this masking into the forloop instead.
  ;m=24.0
  ;mean_star_rv=mean(pos[where(trig eq 1)])
  ;selx=where(RV gt mean_star_rv-m and RV lt mean_star_rv+m) ;Assumes that RV is monotonic.
  ;sely=where(prv gt mean_star_rv-m and prv lt mean_star_rv+m);Assumes that planet RV is monotonic.

  ;ccf_masked=ccf_2D_c*1.0
  ;ccf_masked[min(selx):max(selx),min(sely):max(sely)]=!values.f_nan
  

  
  a=dindgen(200+1)/20.0-5.0 ;Number of orders of magnitude. So we scale between -4 and +4 orders of magnitude times the current version of the doppler shadow model.

  S=a*0.0
  for i=0,n_elements(a)-1,1 do S[i]=total((ccf_masked-(10.0^a[i])*doppler_model_masked)^2.0,/nan)
  
  min=min(S,loc)
  print,a[loc]
  
  b=2.0^(dindgen(500+1)/100.0-2.0)*10^(a[loc])
  S=b*0.0
  for i=0,n_elements(b)-1,1 do S[i]=total((ccf_masked-b[i]*doppler_model_masked)^2.0,/nan)
  min=min(S,loc)

  cgplot,b,S,/yno,xs=1,title='Minimalization of Doppler model factor',xtitle='Factor',ytitle='Sum of squares difference',/xlog

  ;writefits,'test.fits',[[[ccf_2D_c]],[[doppler_model*b[loc]]],[[ccf_2D_c-doppler_model*b[loc]]]]

  print,'Doppler model removed with a factor of '+trim(b[loc])
  writefits,inpath+'CCF_stacked_2D_c_clean.fits',ccf_2D_c-doppler_model*b[loc]
  writefits,inpath+'CCF_stacked_2D_i_clean.fits',ccf_2D_i-doppler_model*b[loc] 
end






;function ccf_shadow,v,p
                                ;p[0]=position
                                ;p[1]=amplitude
  ;p=
