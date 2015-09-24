;*************************************************
; for calculating anomaly correlation 
;
; Concepts illustrated:
;   - Calculating a two-dimensional correlation in time
;   - Reordering an array
;   - Copying attributes from one variable to another
;   - Copying coordinate arrays from one variable to another
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;************************************************
; file handling
;************************************************

 varname  = systemfunc ( " echo ${varname} " )
; varname  = "sst"
 modfn  =  systemfunc ( " echo ${Anomaly_Prefn} " )
 reffn  =  systemfunc ( " echo ${Anomaly_Perfn} " )
 t1fn   =  systemfunc ( " echo ${acct1fn} " )

 t1in   = addfile(t1fn,"r")
 modin  = addfile(modfn,"r")
 refin  = addfile(reffn,"r")
 
 t1_var  = t1in->$varname$ 
 mod_var = modin->$varname$
 ref_var = refin->$varname$
; printVarSummary(t1_var)
; printVarSummary(mod_var)
; printVarSummary(ref_var)
; exit
;************************************************
; reorder to get time as right most dimension
;***********************************************
; mod_new = rm_single_dims(mod_var)
 modv = mod_var(depth|:,latitude|:,longitude|:,time|:)
 printVarSummary(modv)
; ref_new = rm_single_dims(ref_var)
 refv = ref_var(depth|:,latitude|:,longitude|:,time|:)
 printVarSummary(refv)

;************************************************
; calculate cross correlations
;************************************************
 maxlag = 0                                    ; set lag
; note, the max lag should not be more than N/4
  
 ccr = esccr(modv,refv,maxlag)                   ; calc cross correlations

; copy meta data and coordinate variables using contributed functions
 copy_VarAtts(modv, ccr)                    
 copy_VarCoords_1(refv,ccr)
 printVarSummary(ccr)
 print(ccr(0,:,72,0))
 ccr!3 = "time" ; rename dimension
 ccrv = ccr(time|:,depth|:,latitude|:,longitude|:)
; printVarSummary(ccr)
; print(ccr(0,:,72,0))

; *********************************************************
; create target file
;
 outn = systemfunc ( " echo ${ACCoutfn} " )
 system (" rm -f " + outn )
 ncdf   = addfile(outn,"c")      ; create new netCDF file

 att_names = getvaratts(t1in)           ; get CCM file's global attributes
 if(.not.all(ismissing(att_names))) then
   do i = 0,dimsizes(att_names)-1
     print("copy_fileatts: global attributes->" + att_names(i))
     ncdf@$att_names(i)$ = t1in@$att_names(i)$ ;copy input file's global attributes
   end do
 end if

 names = getfilevarnames(t1in)
 do i=0, dimsizes(names)-1
   if ((names(i) .eq. varname ) ) then
     print("writing ncdf: i,name="+ i+ "   "+ names(i) + " modfied ")
     filevardef(ncdf, varname ,typeof(t1_var)  ,getvardims(ccrv))
     filevarattdef(ncdf,varname ,t1_var)
; printVarSummary(ncdf->T)
; printVarSummary(ccrv)
     ncdf->$names(i)$ = (/doubletofloat(ccrv)/)
   else
     print("writing ncdf: i,name="+ i+ "   "+ names(i))
     ncdf->$names(i)$ = t1in->$names(i)$
   end if
 end do

 print ("Finish calculating correlation coeff ")
exit


end






	
