;*************************************************
; wind_3.ncl
;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
;*************************************************
begin
;*************************************************
; open file and read in data: data are on a gaussian grid
;*************************************************
; var = "UV200_04" ;
; ufile = "/work/earnest/UPData/Bias_FF02/SelLev_NorESM/Bias_U200_04.nc"
; vfile = "/work/earnest/UPData/Bias_FF02/SelLev_NorESM/Bias_V200_04.nc"
; fu =addfile("/work/earnest/UPData/Bias_FF02/SelLev_NorESM/Bias_U200_04.nc","r")
; fv =addfile("/work/earnest/UPData/Bias_FF02/SelLev_NorESM/Bias_V200_04.nc","r")
 fu =addfile(ufile,"r")
 fv =addfile(vfile,"r")
  u    = fu->U850             ; (time,lat,lon)
  v    = fv->V850
;*************************************************
; calculate divergence on a gaussian grid
;*************************************************
  div = uv2dvG_Wrap(u,v)                ; u,v ==> divergence

;*************************************************
; calculate velocity potential
;*************************************************
  chi = ilapsG_Wrap ( div , 0)
  chi = (/chi/1e6/)                     ; arbitrary scale
  chi@long_name = "velocity potential"
  chi@units     = "m/s" 
;*************************************************
; calculate divergent wind component 
; --
; note: the calculation uses a procedure, so memory 
; must be preallocated.
;*************************************************    
  ud    = new ( dimsizes(u), typeof(u), u@_FillValue )
  vd    = new ( dimsizes(v), typeof(v), v@_FillValue )

  dv2uvg(div,ud,vd) ; div  ==> divergent  wind components

  copy_VarCoords(u, ud ) 
  copy_VarCoords(u, vd ) 
  ud@long_name  = "Zonal Divergent Wind"
  ud@units      = u@units
  vd@long_name  = "Meridional Divergent Wind"
  vd@units      = v@units
;*************************************************
; plot results
;*************************************************    
  wks  = gsn_open_wks("pdf","wind_NA_" + var)              ; open a workstation
;  cmap = read_colormap_file("BlAqGrYeOrReVi200")
  cmap = read_colormap_file("NCV_blue_red")
  nc   = dimsizes(cmap(:,0))

  res                      = True

  res@cnFillOn             = True               ; color on
  res@cnLinesOn            = False              ; turn off contour lines

  res@gsnScalarContour     = True               ; vectors over contours
  res@cnFillPalette        = cmap(:nc-4,:)
  
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = -8               ; set min contour level
  res@cnMaxLevelValF       = 8                ; set max contour level
  res@cnLevelSpacingF      = 1               ; set contour spacing
  res@lbOrientation = "Vertical"
  res@pmLabelBarOrthogonalPosF = -0.01 ; move label bar closer

;  res@mpFillDrawOrder = "PostDraw"
;  res@mpOutlineDrawOrder = "PostDraw"
;  res@mpOutlineBoundarySets = "AllBoundaries" 

  res@vcRefMagnitudeF      = 2.                 ; make vectors larger
  res@vcRefLengthF         = 0.050              ; ref vector length
  res@vcGlyphStyle         = "CurlyVector"      ; turn on curly vectors
  res@vcMinDistanceF       = 0.012              ; thin the vectors
  res@vcRefAnnoOrthogonalPosF  = -0.13         ; Move ref anno into plot

;  res@mpLandFillColor       = "Transparent"        ; change continent color 
  res@mpLandFillColor       = "gray"        
  res@mpMinLonF            =  -100.               ; select a subregion
  res@mpMaxLonF            =  40.
  res@mpMinLatF            =   10.
  res@mpMaxLatF            =  90.


  res@tiMainString         = "Velocity Potential via Spherical Harmonics"
  res@gsnCenterString      = "Chi scaled by 1e6"
  res@gsnLeftString        = "Divergent Wind"

  plot=gsn_csm_vector_scalar_map_ce(wks,ud(0,:,:),vd(0,:,:),chi(0,:,:),res)
end

