function ccf_unskew,ccf_matrix,ccf_v,rv_planet
  ccf_unsk=ccf_matrix*0.0
  nexp=n_elements(ccf_matrix[0,*])
  for j=0,nexp-1 do begin
     ccf_unsk[*,j]=interpol(ccf_matrix[*,j],ccf_v-rv_planet[j],ccf_v)
     sel=where(ccf_v+rv_planet[j] lt min(ccf_v) or ccf_v+rv_planet[j] gt max(ccf_v))
     ccf_unsk[sel,j]=!values.f_nan
  endfor
  return,ccf_unsk
end
