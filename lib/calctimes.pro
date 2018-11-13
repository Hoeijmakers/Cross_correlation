FUNCTION calctimes,dp,diag=diag,alpha=alpha,transit=transit,dRV=dRV,scattering=scattering,out_transit_scattering=out_transit_scattering,vsys=vsys,vsysxcor=vsysxcor,Julianday=Julianday,HJulday=HJulday,phase_angle=phase_angle,Texp=Exposure_time,altTP=altTP,airmass=airmass,custom_v_orb=custom_v_orb
;This program calculates the time differences between the exposures and
;transit center, the associated orbital phase and associated radial
;velocity of the planet.
;The end result is the radial velocities of each exposure.
;It can also output the time derivitive of the RV (i.e. what is the
;RV-change per second?); the orbital phase of the planet at each
;exposure, the transit profile (0 during transit, 1 out of transit),
;or the forward-scattering profile as expected from Henyey-Greenstein
;scattering. The latter two outputs are useful as weights for
;co-addition of the cross-correlation result. Finally, you can also
;set scattering_out_transit. This does the same as the scattering
;keyword, but it masks out the transit times.

;I grepped the header information of the exposures using readheaders.pro then
;reformatted the data manually so I can calculate the difference in
;number of hours, minutes and seconds from the transit center time.
;It reads three datafiles containing:
;The exposure timings in UT.
;The date at which each exposure is taken (UT).
;Generic information about the transit [RA, DEC, P, Tc, a] 

;See PRO readheaders for more info on these files.

;This program is used in the analysis after xcor.pro to determine vi,
;and to unskew the data.

;================================================================
                      ;FUTURE IMPROVEMENTS
;================================================================
                                ;Use a real transit model. You might
                                ;throw away 5 to 10% of your signal
                                ;this way

                                ;Improve the forward-scattering model
                                ;to include the size of the disk of
                                ;the star.

  
;if keyword_set(dp) eq 0 then dp=dataget()
  
obs_times=dp+'/obs_times'

;readcol,obs_times,year,month,day,hr,min,sec,Texp,format='D,D,D,D,D,D,D',/silent
readcol,obs_times,reduced_JD,datestring,Texp,format='(D,A,D)',/silent
if keyword_set(Airmass) then begin
   readcol,obs_times,reduced_JD,datestring,Texp,airm,format='(D,A,D,D)',/silent   
   ;readcol,obs_times,year,month,day,hr,min,sec,Texp,airm,format='D,D,D,D,D,D,D,D',/silent  
   return,airm
endif

RA=tenv(paramget('RA',dp))*15.0
DEC=tenv(paramget('DEC',dp))
P=paramget('P',dp)
Tc=paramget('Tc',dp)
r=paramget('a',dp)
duration=paramget('duration',dp)
inclination=paramget('inclination',dp)


if keyword_set(altTP) then begin
   P=altTP[1]
   Tc=altTP[0]
endif



Td=duration/60.0d/24.0d ;in days

JD=double(make_array(n_elements(datestring)))
for i=0,n_elements(datestring)-1 do JD[i]=date_conv(datestring[i],'julian')
;JULDAY(month,day,year,hr,min,sec) ;JD of each exposure.
if keyword_set(julianday) then return,JD
HJD=JD*0.0


file_mkdir,dp+'/0_calctimes'

for i=0,n_elements(JD)-1,1 do begin
   dt=helio_jd(JD[i]-2400000d,ra,dec,/time_diff);In seconds
   HJD[i]=JD[i]+dt/3600.0d/24.0
                                ;I need to add dt (4.5 minutes), because
                                ;helio_JD returns to me the amount of
                                ;lightseconds which the Earth is
                                ;CLOSER to the Sun than the
                                ;target. The HJD is thus dt seconds
                                ;later than the local JD.
endfor

if keyword_set(Hjulday) then return,HJD

Tc_n=Tc-1000000d*P;This is to make sure that the Transit central time PRECEDES the observations (by tens or hundreds or thousands of years). Otherwise, the phase could pick up a minus sign somewhere and be flipped. I hate that.
phase=((HJD-Tc_n) MOD P)/P
v_orb=2d*!pi*(r*149598000d)/(P*86400d)

if keyword_set(custom_v_orb) then v_orb=custom_v_orb
;r comes in AU, P in days. Needs to be converted to km/s.
RV=v_orb*sin(phase*2d*!pi)

returnpar=RV*sin(inclination/!radeg)
if keyword_set(alpha) then begin
   return,phase
endif

if keyword_set(phase_angle) then begin
                                ;This is the phase angle of the
                                ;planet. Yani the fraction of the
                                ;surface that you see. Phase angle 0
                                ;is full-dayside.
   PhA=acos(-1.0*sin(inclination/!radeg)*sin(2.0*!pi*phase+0.5*!pi))
   return,PhA
endif


if keyword_set(scattering) OR keyword_set(out_transit_scattering) then begin
   loadct,4,/silent
   g1=paramget('HG_scattering_g1',dp)
   g2=paramget('HG_scattering_g2',dp)
   k=paramget('HG_scattering_k',dp)
   P_s=henyey_greenstein(phase*2*!pi,g1,g2,k)
                                ;The following few lines were copied
                                ;from the transit-ifloop below.
                                ;Because I want to visualize the in
                                ;and egress times.

   if keyword_set(out_transit_scattering) then begin
      transitprofile=calctimes(/transit) ;Can I call a program from within itself? Seems so! :D
      P_n=P_s*(1.0-transitprofile)
   endif else begin
      P_n=P_s*1.0
   endelse
   
   n_orb=(mean(HJD)-Tc) /P 
   center=Tc+round(n_orb)*P
   ingresstime=(center-0.5*Td)
   egresstime=(center+0.5*Td)
   
   cgps_open,dp+'/0_calctimes/forward_scattering_profile.ps',/quiet
   cgplot,HJD,P_n/max(P_n),thick=3.0,xtitle='HJD',ytitle='Normalized scattering intensity',title='Forward scattering during transit',xcharsize=0.6,ycharsize=0.6
   vline,ingresstime,color=100
   vline,egresstime,color=100
   cgtext,0.85*(max(HJD)-min(HJD))+min(HJD),0.9*(!y.crange[1]-!y.crange[0])+!y.crange[0],'g!If!N='+trim(g1),charsize=0.9
   cgtext,0.85*(max(HJD)-min(HJD))+min(HJD),0.85*(!y.crange[1]-!y.crange[0])+!y.crange[0],'g!Ib!N='+trim(g2),charsize=0.9
   cgtext,0.85*(max(HJD)-min(HJD))+min(HJD),0.8*(!y.crange[1]-!y.crange[0])+!y.crange[0],'k='+trim(k),charsize=0.9
   ;print,0.6*max(HJD),0.8*(!y.crange[1]-!y.crange[0])+!y.crange[0]
   cgps_close   

   return,P_n/max(P_n)
endif

if keyword_set(transit) then begin
   loadct,4,/silent
   n_orb=(mean(HJD)-Tc) /P
   center=Tc+round(n_orb)*P
   profile=HJD*0.0+1.0;Transit profile
   profile[where(HJD lt center-0.4*Td)]=0.0
   profile[where(HJD gt center+0.4*Td)]=0.0
   blur=n_elements(where(profile gt 0.5))/5.0
   profile_s=smooth(profile,blur)
   timehours=((HJD-0.5) MOD 1.0) * 24.0
   xmin=min(timehours)-3.0
   xmax=max(timehours)+3.0
   ingresstime=((center-0.5*Td-0.5) MOD 1.0)*24.0
   egresstime=((center+0.5*Td-0.5) MOD 1.0)*24.0
   midtime=((center-0.5) MOD 1.0)*24
   if ingresstime gt 12 then ingresstime=ingresstime-24
   if midtime gt 12 then midtime=midtime-24
   if egresstime gt 12 then egresstime=egresstime-24
   ;print,'Ingress, mid, egress:'+trim(ingresstime)+' '+trim(midtime)+' '+trim(egresstime)
   set_plot,'ps'
   device,filename=dp+'/0_calctimes/transit_timing.ps',/color
   ;plot,HJD,RV,psym=1,xtitle='HJD',ytitle='RV(km/s)'
   plot,timehours,RV,psym=1,xtitle='UTC',ytitle='RV(km/s)',xrange=[xmin,xmax],xstyle=1
   ;vline,center,color=200
   ;vline,center-0.5*Td,color=100
   ;vline,center+0.5*Td,color=100
   vline,midtime,color=200
   vline,ingresstime,color=100
   vline,egresstime,color=100
   device,/close
   ;print,timehours
   device,filename=dp+'/0_calctimes/transit_profile.ps'
   plot,timehours,profile_s,psym=1,yrange=[-0.1,1.1],ys=1,xtitle='UTC',ytitle='Weight',xrange=[xmin,xmax],xstyle=1
   ;vline,center-0.5*Td,color=100
   vline,ingresstime,color=100
   vline,egresstime,color=100
      vline,midtime,color=200
   device,/close
   set_plot,'x'
   return,profile_s
endif
if keyword_set(dRV) then begin
   dRV=v_orb*cos(phase*2d*!pi)*2d*!pi/(P*24d*3600d)*sin(inclination/!radeg)
   return,abs(dRV*Texp)
endif

if keyword_set(Exposure_time) then return,Texp


if keyword_set(vsys) then begin
                                ;Taken from the baryvel documentation:
   ;http://www.exelisvis.com/docs/baryvel.html
   ;     Compute the radial velocity of the Earth toward Altair on 15-Feb-1994 
   ;       using both the original Stumpf algorithm and the JPL ephemeris 
   ;   IDL> jdcnv, 1994, 2, 15, 0, jd ;==> JD = 2449398.5 
   ;   IDL> baryvel, jd, 2000, vh, vb ;Original algorithm 
   ;           ==> vh = [-17.07243, -22.81121, -9.889315] ;Heliocentric km/s 
   ;           ==> vb = [-17.08083, -22.80471, -9.886582] ;Barycentric km/s 
   ;   IDL> baryvel, jd, 2000, vh, vb, /jpl ;JPL ephemeris 
   ;           ==> vh = [-17.07236, -22.81126, -9.889419] ;Heliocentric km/s 
   ;           ==> vb = [-17.08083, -22.80484, -9.886409] ;Barycentric km/s 
   ;   IDL> ra = ten(19,50,46.77)*15/!RADEG ;RA in radians 
   ;   IDL> dec = ten(08,52,3.5)/!RADEG ;Dec in radians 
   ;   IDL> v = vb[0]*cos(dec)*cos(ra) + $ ;Project velocity toward star 
                                ;           vb[1]*cos(dec)*sin(ra) +
                                ;           vb[2]*sin(dec)
   vsys=jd*0.0
   systemic_rv=paramget('rvsys')
   for i=0,n_elements(jd)-1,1 do begin
      baryvel,JD[i],2000,vh,vb
      vsys[i] = vb[0]*cos(dec/!radeg)*cos(ra/!radeg) +vb[1]*cos(dec/!radeg)*sin(ra/!radeg) + vb[2]*sin(dec/!radeg)
   endfor
                                ;Returns the speed of the EARTH TOWARD
                                ;the target star. Radial velocity is
                                ;positive if the object is moving away
                                ;from Earth. So if baryvel returns a
                                ;negative number, it means that we are
                                ;moving away from the star. I convert
                                ;that number to a positive number, to
                                ;denote radial velocity. Add this
                                ;number to the radial velocity given
                                ;by SIMBAD to get the radial velocity
                                ;that should be retrieved by the
                                ;wavelength solution and model
                                ;injection codes. If Vsys(here) +
                                ;Vsys(simbad) dont add up, there is a problem.
   return,-1.0*vsys+systemic_rv  

endif



return,returnpar
END


function henyey_greenstein,th,g1,g2,k
                                ;This short routine calculates the
                                ;bimodal Henyey-Greenstein scattering
                                ;function, as a function of scattering
                                ;angle, scattering parameters g and
                                ;assymetry factor k for assymetry
                                ;between the forward and
                                ;backscattering component. Set g to 1
                                ;for forward scattering, 0 for isotropic
                                ;scattering and -1 for backward
                                ;scattering. Set k =  1 for total
                                ;forward scattering. k=0 for total
                                ;backscattering.
  
  if k gt 1 or k lt 0 or g1 gt 1 or g1 lt 0 or g2 gt 0 or g2 lt -1 then begin
     print,'Henyey Greenstein parameter error. 0 < g1 < 1,    -1 < g2 < 0,    0 < k < 1'
     return,0
  endif
  P=k*(1.0/(4*!pi)) * (1-g1^2) / ((1.0+g1^2 - 2*g1 * cos(th))^1.5)  +  (1.0-k) * (1.0/(4*!pi)) * (1-g2^2) / ((1.0+g2^2 - 2*g2 * cos(th))^1.5)
  return,P
end
