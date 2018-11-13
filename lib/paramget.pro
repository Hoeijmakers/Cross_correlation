FUNCTION paramget,keyword,dp
;This program retrieves a keyword from the configuration file of a dataset.


readcol,dp+'/config',names_g,params_g,format='(A,A)',/silent
names_all=names_g
params_all=params_g

;names_all=[names,names_g];These contain the parameters from the config file, as 
;params_all=[params,params_g];well as the generics file. All still in STRING format.

paramline=where(names_all eq keyword)
if n_elements(paramline) gt 1.0 then begin
   print,'ERROR. Confusion in keywords (multiple hits)'
   return,'ERROR'
endif
if paramline eq -1.0 then begin
   print,'ERROR. Keyword "'+keyword+'" not found.'
   return,'ERROR'
endif
value_string=params_all[paramline]
isitanumber=strnumber(value_string,value_number)
if isitanumber eq 1 and valid_num(value_string) eq 1 then return,value_number[0] else return,value_string[0]
END
