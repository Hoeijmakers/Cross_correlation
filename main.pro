pro run_instance,multi,paramfile,only_xcor=only_xcor,only_stacking=only_stacking
                                ;This is the main executor of the
                                ;entire cross-correlation routine. It
                                ;loops through all spectral orders,
                                ;injecting the model into, cleaning
                                ;and cross-correlating each of
                                ;them. It can be run in parallel using
                                ;the multi variable, which is a
                                ;2-element array [a,b] specifying a)
                                ;the CPU this instance is run on (1 or
                                ;above), and b) the total
                                ;number of CPU's that you are
                                ;parallelizing on. so CPU 1 out of 4
                                ;will be called  as
                                ;run_instance,[1,4],'paramfile', 2 out
                                ;of 4 as [2,4], etc. Instance [1,N]
                                ;is the primary instance, which after
                                ;cross-correlation of all orders will
                                ;wait until the completion of the
                                ;parallel instances, and stack the
                                ;resulting CCFs, producing the primary
                                ;output of the cross-correlation
                                ;routine.

                                ;Check the output/dataset/modelname
                                ;folder for the output of the various
                                ;routine to check the quality
                                ;of the cleaned data, the CCF's
                                ;and their stacked products.




  ;==========PRIMARY CODE BELOW=========;

 
  ccv=0
  fast=0.0
  normalize_model=0.0
  doppler_model='/'
  OPENR, lun, paramfile, /GET_LUN
  commands = ''
  line = ''
  WHILE NOT EOF(lun) DO BEGIN
     READF, lun, line
     commands = [commands, line]
  ENDWHILE
; Close the file and free the file unit
  FREE_LUN, lun

  for i=0,n_elements(commands)-1 do void=execute(commands[i])

  if keyword_set(only_stacking) then multi=[1,1]
  
  for i=multi[0]-1,n_elements(orders)-1,multi[1] do begin
                                ;TIC
     if keyword_set(only_xcor) ne 1 and keyword_set(only_stacking) ne 1 then begin
        inject_model,dataset,orders[i],modelname
        clean_data,dataset,orders[i],npass=npass,w=w,/outtransit
     endif
     
     if keyword_set(only_stacking) ne 1 then xcor_data,dataset,orders[i],templatename,drv,rvrange,ccv=ccv,fast=fast,normalize_model=ccv,edgew=w
     ;TOC
  endfor

  ;To end the parallel part:
  openw,multi[0]+1,paramfile+'_p_'+trim(multi[0])
  close,multi[0]+1
  if multi[0] eq 1 then begin
     tracker=findgen(multi[1])+1
     trigger=tracker*0.0
     while n_elements(where(trigger eq 1)) ne n_elements(trigger) do begin
        for i=0,n_elements(tracker)-1,1 do trigger[i]=file_test(paramfile+'_p_'+trim(tracker[i]))        
        statusline,'Waiting for '+trim(n_elements(trigger)-n_elements(where(trigger eq 1)))+' processes'
        wait,2.0
        statusline,/clear
     endwhile
   
     for i=0,n_elements(tracker)-1 do file_delete,paramfile+'_p_'+trim(tracker[i])
     stack_ccfs,dataset,templatename,orders,ccv=ccv
     stack_ccfs,dataset,templatename,orders,ccv=ccv,/injected
 
     co_add_ccfs,orders,dataset,templatename,ccv=ccv;,constantweights=normalize_model,noiseweights=normalize_model
     if doppler_model ne '/' then clean_ccfs,dataset,templatename,doppler_model
     make_kpvsys,dataset,templatename
     if doppler_model ne '/' then make_kpvsys,dataset,templatename,/clean

     
     if multi[1] gt 1 then begin
        for fff=2,multi[1] do begin
           openw,1,paramfile+'_f_'+trim(fff)
           close,1
        endfor
     endif
     
     endif else begin
        print,'Multi complete.'
        trigger=0.0
        while file_test(paramfile+'_f_'+trim(multi[0])) ne 1 do begin
           statusline,'Waiting for stacking to complete.'
           wait,1.0
           statusline,/clear
        endwhile
        file_delete,paramfile+'_f_'+trim(multi[0])
     endelse
     
end






pro run_elements_kelt9,core,night
                                ;This script creates instances files
                                ;and executes run_instance for each
                                ;element in each night. It should be
                                ;run in 4 sessions in parallel. Core 1
                                ;and 2 for night 1 and 2 respectively. 

;elements=['LalI','LoI','LoII','Leendert']
elements=['NaI','MgI','AlI','SiI','CaI','CaII','ScI','ScII','TiI','TiII','VI','VII','CrI','CrII','MnI','MnII','FeI','FeII','CoI','NiI','NiII','CuI','ZnI']


  
instancelines=['modelname = templatename','orders    =    [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]','npass     =    0.0','w         =    75.0','drv	  =    2.0','rvrange	  =    1000.0','ccv=1']

if night eq 1 then datasetline="dataset   =    'Kelt-9/night1'" else datasetline="dataset   =    'Kelt-9/night2'"


for i=0,n_elements(elements)-1 do begin
   print,'>>>>>>> RUNNING '+elements[i]+' <<<<<<<'
   instancefile='instances/kelt_9_elements/'+elements[i]+'_4k_n'+trim(night)
   if core eq 1 then begin
      lunnie=core+5*night
      openw,lunnie,instancefile
      printf,lunnie,datasetline
      printf,lunnie,"templatename =  'Dkelt9_"+elements[i]+"_4000_0'"
      for j=0,n_elements(instancelines)-1 do printf,lunnie,instancelines[j],format='(A)'
      printf,lunnie,"doppler_model= 'output/'+dataset+'/doppler_model_kelt9_n"+trim(night)+".txt'"
      close,lunnie
   endif else while file_test(instancefile) ne 1 do begin
      wait,2
      print,'Waiting to start'
   endwhile
   
   run_instance,[core,2],instancefile

   
endfor

   



END




pro compare_CCFs
  orders=findgen(69)
  dataset='WASP-33/night1'
  modelnames=['WASP-33-norlander','W33-R2-norlander','W33-R3-norlander','W33-R4-norlander','W33-R5-norlander']
  SNs=findgen(n_elements(modelnames))*0.0
  for i=0,n_elements(modelnames)-1,1 do begin
     modelname=modelnames[i]
     stack_ccfs,dataset,modelname,orders,sigma=sigma
     stack_ccfs,dataset,modelname,orders,/injected,signal=signal  
     SNs[i]=signal/sigma
  endfor

  c1=readfits(modelnames[0]+'.fits')
  c2=readfits(modelnames[1]+'.fits')
  c3=readfits(modelnames[2]+'.fits')
  c4=readfits(modelnames[3]+'.fits')
  c5=readfits(modelnames[4]+'.fits')

  stack_ccfs,dataset,modelnames[0],orders;Now run for norlander again to get the uninjected CCF.
  c6=readfits(modelnames[0]+'.fits')

  dy=0.002
  cgps_open,'CCF_compare.ps'
  cgplot,c5[*,0],c5[*,1],xtitle='Radial velocity (km/s)',ytitle='Cross-correlation',charsize=1.0,yrange=[min(c5[*,1])-dy*6,max(c5[*,1]+dy)],ys=1
  cgplot,c4[*,0],c4[*,1]-dy,/overplot
  cgplot,c3[*,0],c3[*,1]-dy*2,/overplot  
  cgplot,c2[*,0],c2[*,1]-dy*3,/overplot
  cgplot,c1[*,0],c1[*,1]-dy*4,/overplot
  cgplot,c6[*,0],c6[*,1]-dy*5,/overplot
  cs=0.9
  cgtext,plotpos(0.98,/x),plotpos(0.94,/y),'R = 5R$\downJ$ - S/N = '+r_trim(SNs[4],1),/align,charsize=cs
  cgtext,plotpos(0.98,/x),plotpos(0.88,/y),'R = 4R$\downJ$ - S/N = '+r_trim(SNs[3],1),/align,charsize=cs
  cgtext,plotpos(0.98,/x),plotpos(0.82,/y),'R = 3R$\downJ$ - S/N = '+r_trim(SNs[2],1),/align,charsize=cs
  cgtext,plotpos(0.98,/x),plotpos(0.76,/y),'R = 2R$\downJ$ - S/N = '+r_trim(SNs[1],1),/align,charsize=cs
  cgtext,plotpos(0.98,/x),plotpos(0.70,/y),'R = 1.603R$\downJ$ - S/N = '+r_trim(SNs[0],1),/align,charsize=cs
  cgtext,plotpos(0.98,/x),plotpos(0.64,/y),'Uncontaminated',/align,charsize=cs
  cgps_close
  stop
end





pro model_wasp_31

                                ;This program models the SNR of the
                                ;wasp 31 potassium lines for the
                                ;ESPRESSO proposal in cycle 102.


  a=get_model('W31-lorenzo')
  wlm=a[*,0]
  fxm=a[*,1]
  Rd=135000.0
  sigma=resolution_kernel(Rd,1e7,760.0)
  sel=where(wlm gt 700 and wlm lt 800)
  mpxs=(max(wlm[sel])-min(wlm[sel]))/n_elements(sel)

  fxms=gauss_smooth(fxm,sigma/mpxs,/edge_truncate)
  ;cgps_open,'test.ps' 
  ;cgplot,wlm,fxms,xrange=[766,767],/yno
  ;cgplot,wlm,fxms,/overplot,color='red'

  minwl=380.0
  maxwl=790.0
  l=760.0;Target wavelength
  s=2.9;Sampling rate of ESPRESSO
  dl=l/Rd/3.3
  wld=findgen((maxwl-minwl)/dl+1)*dl+minwl

  fxmi=interpol(fxms,wlm,wld)

  SNR=25.0;SNR per exposure as from ETC
  nexp=35 ;exposures in transit
  ntransit=3
  nbin=4
  ;wlb=findgen((maxwl-minwl)/dl/nbin)*dl*nbin+minwl
  noise=randomn(321,n_elements(wld))/SNR/sqrt(nexp)/sqrt(ntransit)
  data=fxmi+noise
  i=0.0
  wlb=[]
  fxb=[]
  while i+nbin-1 lt n_elements(wld) do begin
     wlb=[wlb,mean(wld[i:i+nbin-1])]
     fxb=[fxb,mean(data[i:i+nbin-1])]
     i=i+nbin
  endwhile
  
  err=wlb*0.0+1.0/SNR/sqrt(nexp)/sqrt(ntransit)/sqrt(nbin)
  xmin=769.2
  xmax=771.0
  ;xmin=588.5
  ;xmax=590.5
  sel=where(wlb gt xmin and wlb lt xmax)


  cgps_open,'Potassium_w31.ps'
    cgdisplay,2000,450
  cgplot,wlm,fxms,xrange=[xmin,xmax],xs=1,/yno,xtitle='Wavelength (nm)',ytitle='Relative flux ($\Delta$F / F)',charsize=0.8,thick=1.5,yrange=[0.96,0.9865],Title='Red-most potassium D-line',linestyle=1
  cgplot,wlb[sel],fxb[sel],/overplot,color='red',err_yhigh=err[sel],err_ylow=err[sel],psym=16,err_width=(xmax-xmin)/n_elements(sel)/7.0,thick=0.2,symsize=0.2,err_thick=0.5  
  cgps_close


  b=get_model('sky')
  wls=b[*,0]
  fxs=b[*,1]
  cgps_open,'Potassium_W31-sky.ps'
  cgplot,wlm,fxms,yrange=[0.6,1.05],xrange=[766.0,770.5],xtitle='Wavelength (nm)',ytitle='Transmission',charsize=0.8
  cgplot,wls,fxs,/overplot,color='skyblue'
  cgps_close
  stop

end


pro model_wasp_1731

                                ;This program models the SNR of the
                                ;wasp 17 and 31 potassium lines for the
                                ;ESPRESSO proposal in cycle 103.



  T31=1600.0
  T17=1700.0
  mu=2.3
  G=6.674e-11
  Rj=69911000.0;m
  Rp31=1.537*Rj
  Rp17=1.991*Rj
  Rs=1.38
  mp=0.48
  Mj=1.898e27;kg
  mH=1.6737e-27;kg
  Rsun=695508000.0              ;m
  k=1.38d-23                    ;SI

  g31=G*mp*MJ/Rp31^2.0
  H31=k*T31/(mu*mH*g31)

  g17=G*mp*MJ/Rp17^2.0
  H17=k*T17/(mu*mH*g17)

  print,2.0*[H31,H17]*[Rp31,Rp17]/(Rs*Rsun)^2.0


  
 
 
  
  
  a=get_model('W31-lorenzo')
  wlm=a[*,0]
  fxm=a[*,1]
  Rd=135000.0
  sigma=resolution_kernel(Rd,1e7,760.0)
  sel=where(wlm gt 700 and wlm lt 800)
  mpxs=(max(wlm[sel])-min(wlm[sel]))/n_elements(sel)

  fxms=gauss_smooth(fxm,sigma/mpxs,/edge_truncate)
  ;cgps_open,'test.ps' 
  ;cgplot,wlm,fxms,xrange=[766,767],/yno
  ;cgplot,wlm,fxms,/overplot,color='red'

  minwl=380.0
  maxwl=790.0
  l=760.0;Target wavelength
  s=2.9;Sampling rate of ESPRESSO
  dl=l/Rd/3.3
  wld=findgen((maxwl-minwl)/dl+1)*dl+minwl

  fxmi=interpol(fxms,wlm,wld)

  SNR=35.0;SNR per exposure as from ETC
  nexp=35 ;exposures in transit
  ntransit=2
  nbin=2
  ;wlb=findgen((maxwl-minwl)/dl/nbin)*dl*nbin+minwl
  noise=randomn(322,n_elements(wld))/SNR/sqrt(nexp)/sqrt(ntransit)
  data=fxmi+noise
  i=0.0
  wlb=[]
  fxb=[]
  while i+nbin-1 lt n_elements(wld) do begin
     wlb=[wlb,mean(wld[i:i+nbin-1])]
     fxb=[fxb,mean(data[i:i+nbin-1])]
     i=i+nbin
  endwhile
  
  err=wlb*0.0+1.0/SNR/sqrt(nexp)/sqrt(ntransit)/sqrt(nbin)
  xmin=769.9
  xmax=770.3
  ;xmin=588.5
  ;xmax=590.5
  sel=where(wlb gt xmin and wlb lt xmax)


  cgps_open,'Potassium_w31_new.ps'
    cgdisplay,700,450
  cgplot,wlm,fxms,xrange=[xmin,xmax],xs=1,/yno,xtitle='Wavelength (nm)',ytitle='Relative flux ($\Delta$F / F)',charsize=0.8,thick=1.5,yrange=[0.96,0.987],Title='WASP-31 b: Red-most potassium line',linestyle=1
  cgplot,wlb[sel],fxb[sel],/overplot,color='red',err_yhigh=err[sel],err_ylow=err[sel],psym=16,err_width=(xmax-xmin)/n_elements(sel)/7.0,thick=0.2,symsize=0.2,err_thick=0.5  
  cgps_close

  ps_to_pdf,'Potassium_w31_new'

  b=get_model('sky')
  wls=b[*,0]
  fxs=b[*,1]
  cgps_open,'Potassium_W31-sky.ps'
  cgplot,wlm,fxms,yrange=[0.6,1.05],xrange=[766.0,770.5],xtitle='Wavelength (nm)',ytitle='Transmission',charsize=0.8
  cgplot,wls,fxs,/overplot,color='skyblue'
  cgps_close
  stop

end





pro model_wasp_43_W49,sodium=sodium

                                ;This program models the SNR of the
                                ;wasp 43 and W49 sodium and potassium
                                ;lines for the ESPRESSO proposals of
                                ;Lorenzo and Aurelien in cycle 102.
  names=['W43-lorenzo','W43-10-lorenzo','W49-lorenzo']
  SNRs=[84,84,60]
  if keyword_set(sodium) then SNRs=[73,73,60]
  nexps=[1,1,12.7]
  ntransits=[3,3,2]
  nbins=[6,6,1]
  for zz=0,2,1 do begin
  
     a=get_model(names[zz])
     wlm=a[*,0]
     fxm=a[*,1]
     Rd=120000.0
     sigma=resolution_kernel(Rd,1e7,760.0)
     sel=where(wlm gt 700 and wlm lt 800)
     mpxs=(max(wlm[sel])-min(wlm[sel]))/n_elements(sel)

     fxms=gauss_smooth(fxm,sigma/mpxs,/edge_truncate)
  ;cgps_open,'test.ps' 
  ;cgplot,wlm,fxms,xrange=[766,767],/yno
  ;cgplot,wlm,fxms,/overplot,color='red'

     minwl=380.0
     maxwl=790.0
     if keyword_set(sodium) then l=589 else l=760.0;Target wavelength
     s=3.3                      ;Sampling rate of ESPRESSO
     dl=l/Rd/3.3
     wld=findgen((maxwl-minwl)/dl+1)*dl+minwl

     fxmi=interpol(fxms,wlm,wld)

  SNR=SNRs[zz];SNR per exposure as from ETC
  nexp=nexps[zz] ;exposures in transit
  ntransit=ntransits[zz]
  nbin=nbins[zz]
  ;wlb=findgen((maxwl-minwl)/dl/nbin)*dl*nbin+minwl
  noise=randomn(123,n_elements(wld))/SNR/sqrt(nexp)/sqrt(ntransit)
  data=fxmi+noise
  i=0.0
  wlb=[]
  fxb=[]
  while i+nbin-1 lt n_elements(wld) do begin
     wlb=[wlb,mean(wld[i:i+nbin-1])]
     fxb=[fxb,mean(data[i:i+nbin-1])]
     i=i+nbin
  endwhile
  
  err=wlb*0.0+1.0/SNR/sqrt(nexp)/sqrt(ntransit)/sqrt(nbin)
  xmin=763.2
  xmax=773.5
  if keyword_set(sodium) then begin
     xmin=587.0
     xmax=592.0
  endif
  
  sel=where(wlb gt xmin and wlb lt xmax)

  if keyword_set(sodium) then cgps_open,'Sodium_'+names[zz]+'.ps' else cgps_open,'Potassium_'+names[zz]+'.ps'
  cgplot,wlm,fxms,xrange=[xmin,xmax],xs=1,/yno,xtitle='Wavelength (nm)',ytitle='Relative flux ($\Delta$F / F)',charsize=0.8,thick=0.5,yrange=[min(fxb[sel]-err[sel]),max(fxb[sel]+err[sel])]
  cgplot,wlb[sel],fxb[sel],/overplot,color='red',err_yhigh=err[sel],err_ylow=err[sel],psym=16,err_width=(xmax-xmin)/n_elements(sel)/7.0,thick=0.2,symsize=0.2,err_thick=0.5  
  cgps_close

  if keyword_set(sodium) then outpath='Sodium_'+names[zz]+'.dat' else outpath='Potassium_'+names[zz]+'.dat'
  openw,1,outpath
  for ii=0,n_elements(sel)-1,1 do printf,1,wlb[sel[ii]],fxb[sel[ii]],err[sel[ii]],format='(A,A,A)'
  close,1

  
  endfor
  

;  b=get_model('sky')
;  wls=b[*,0]
;  fxs=b[*,1]
;  cgps_open,'Potassium_W31-sky.ps'
;  cgplot,wlm,fxms,yrange=[0.6,1.05],xrange=[766.0,770.5],xtitle='Wavelength (nm)',ytitle='Transmission',charsize=0.8
;  cgplot,wls,fxs,/overplot,color='skyblue'
;  cgps_close

end
