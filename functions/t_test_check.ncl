;*************************************************
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;************************************************
; file handling
;************************************************

 N = 52
 dr = 0.002
 do i = 0,300
   r = dr * i 
   t_test = rtest(r, N, 0)
   print("corr: " + r + "  t_test:" + t_test(0) + "  ")
 end do

end






	
