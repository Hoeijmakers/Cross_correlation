FUNCTION get_model,keyword,resolution=resolution,custom_lib=custom_lib
;This program retrieves a model that is indexed in the models/library
;file.


  if keyword_set(custom_lib) ne 1 then library='models/library' else library=custom_lib
  readcol,library,names,paths,res,format='(A,A,F)',/silent


paramline=where(names eq keyword)
if n_elements(paramline) gt 1.0 then begin
   print,'ERROR. Confusion in model name (multiple hits)'
   return,'ERROR'
endif
if paramline eq -1.0 then begin
   print,'ERROR. Model named "'+keyword+'" not found.'
   return,'ERROR'
endif
resolution=res[paramline]
model=readfits('models/'+paths[paramline],/silent)

return,model
END
