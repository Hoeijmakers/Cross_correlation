function interpolate_template,fx,wl,wln

                                ;This function interpolates a hi-res
                                ;(sharp-line) template onto a courses
                                ;(data) wavelength grid, by
                                ;distributing the flux in each hi-res
                                ;point over the two nearest neighbors
                                ;in the lo-res grid. This needs to be
                                ;done if the cross-correlation is done
                                ;on the lo-res grid. An alternative
                                ;approach might be to interpolate the
                                ;data onto the grid of the tempate;
                                ;but that might increase the duration
                                ;of the ccf too much (though this
                                ;thing also takes long).

  wln=wln[sort(wln)]
  fx=fx[sort(wl)]
  wl=wl[sort(wl)]
  
  dwln_start=wln[1]-wln[0]
  dwln_end=wln[-1]-wln[-2]
  wlnp=[wln[0]-dwln_start,wln,wln[-1]+dwln_end] ;Extrapolated output array to be able to deal with edges.
  fxnp=wlnp*0.0

  minsweep=min(wlnp)
  maxsweep=max(wlnp)
  tic

  wlsel=wl[where(wl gt minsweep and wl lt maxsweep)]
  fxsel=fx[where(wl gt minsweep and wl lt maxsweep)]
  for i=0,n_elements(wlsel)-1,1 do begin
     ;if wl[i] gt minsweep and wl[i] lt maxsweep then begin
        diff=wlnp-wlsel[i]
        possel=where(diff gt 0.0)
        negsel=where(diff le 0.0)
        
        nearest_pos=wlnp[possel[0]]
        nearest_neg=wlnp[negsel[-1]]

        posdiff=diff[possel[0]]
        negdiff=diff[negsel[-1]]

        dwl=abs(negdiff)+posdiff
        ;print,fx[i],nearest_neg,nearest_pos,negdiff/dwl,posdiff/dwl,dwlhi

        fxnp[possel[0]]+=posdiff/dwl*fxsel[i]
        fxnp[negsel[-1]]+=abs(negdiff/dwl)*fxsel[i]
        ;posloc=where(diff eq nearest_positive_neighbor)
        ;minloc=where(diff eq nearest_negative_neighbor)

        ;print,wlnp[posloc]
        ;print,wlnp[minloc]
        ;stop
     

     
  endfor
  toc
  stop
end
