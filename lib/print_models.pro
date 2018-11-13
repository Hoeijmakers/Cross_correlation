PRO print_models,custom_lib=custom_lib
;This program prints the contents of the current model library.


  if keyword_set(custom_lib) ne 1 then library='models/library' else library=custom_lib
  readcol,library,names,paths,res,format='(A,A,F)',/silent

  for i=0,n_elements(names)-1,1 do print,names[i]+'   '+paths[i]

END
