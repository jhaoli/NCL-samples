load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  file_spectral = "~/Documents/work/data/spectral/rh/data/spectral_rh4_180x91_dt100.nc"
  file_gamil    = "~/Documents/work/data/gamil-c-grid-pole-u/rh/data/gamil_state_rh_180x90_dt20_pole_u.nc"
  file_gmcore_u = "~/Documents/work/data/gmcore-pole-u/rh/data/state_rh_180x91_midpoint_classic.nc"
  file_tmcore   = "~/Documents/work/data/tmcore/data/tmcore_rh_180x91.nc"

  spectral    = addfile(file_spectral, "r")
  gd_spectral = spectral->gd(88,:,:)/9.8

  gamil    = addfile(file_gamil, "r")
  gd_gamil = gamil->gh(88,:,:)/9.8

  gmcore    = addfile(file_gmcore_u, "r")
  gd_gmcore = gmcore->gh(88,:,:)/9.8

  tmcore = addfile(file_tmcore, "r")
  gd_tmcore = tmcore->gh(39,:,:)/9.8
; create plots 
;**************************
  out_name = "rh4_case"
  wks = gsn_open_wks("ps", out_name)
  gsn_define_colormap(wks, "GMT_panoply")

  plot = new(4, graphic)
  
  res                   = True 
  res@gsnDraw           = False 
  res@gsnFrame          = False 
  ; res@gsnMaximize       = True 
  ; res@gsnPolar          = "NH"
  res@cnFillOn          = False
  res@cnLinesOn         = True
  res@cnLineLabelsOn    = True
  res@cnLineLabelBackgroundColor = "white"
  res@cnLineLabelFontHeightF = 0.012
  res@cnLineLabelFontThicknessF = 0.02
  res@cnInfoLabelOn     = False
  ; res@cnLevelSelectionMode = "ManualLevels"
  ; res@cnMinLevelValF    = -20
  ; res@cnMaxLevelValF    = 20
; res@cnLevels          = (/-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120/)
  res@cnLevelSpacingF   = 200
  res@cnLineThicknessF  = 1.2
  res@gsnMajorLonSpacing=60
  res@gsnMajorLatSpacing=30
  res@gsnMinorLonSpacing=10
  res@gsnMinorLatSpacing=10
  res@gsnContourNegLineDashPattern = 1

  res@lbLabelBarOn      = False 

  res@tmYLLabelFontHeightF = 0.012
  res@tmYLLabelFontThicknessF = 1.0
  res@tmXBLabelFontThicknessF = 0.012
  res@tmYROn            = False
  res@tmXTOn            = False
; res@mpDataBaseVersion = "MediumRes"
  res@mpGridAndLimbOn   = False
  res@mpGridLineColor   = 1
  res@mpPerimOn         = False
  res@mpFillOn          = False
  res@mpOutlineOn       = False
; res@mpFillDrawOrder   = "PreDraw"
  ; mpres@mpMinLatF       = 30
  ; mpres@mpProjection   = "Orthographic"
  ; mpres@mpCenterLatF   = 90
  ; mpres@mpCenterLonF   = 180

  ; mpres@gsnPolarLabelFontHeightF=0.016
  ; res@gsnRightString    = "day=89"
  res@gsnLeftString     = ""
; res@tiMainFontHeightF = 0.02

; 
  res = res
  res@vpXF             = 0.05
  res@vpYF             = 0.95
  res@vpWidthF         = 0.42
  res@vpHeightF        = 0.45
  plot(0) = gsn_csm_contour_map(wks, gd_spectral, res)

  res@vpXF             = 0.53
  plot(1) = gsn_csm_contour_map(wks, gd_gamil   , res)

  res@vpXF             = 0.05
  res@vpYF             = 0.67
  plot(2) = gsn_csm_contour_map(wks, gd_gmcore  , res)
  
  res@vpXF             = 0.53
  plot(3) = gsn_csm_contour_map(wks, gd_tmcore  , res)

  draw(plot)

  ;************************************************
; Attach text to plot using plot coordinates.
;************************************************
  txres               = True                     ; text mods desired
  txres@txPerimOn     = False
  ; txres@txBackgroundFillColor = "white"
  txres@txFontHeightF = 0.015                     ; font smaller. default big
  txres@txFontThicknessF = 0.02

  gsn_text_ndc(wks, "(a) day=89", 0.26, 0.846,txres)
  gsn_text_ndc(wks, "(b) day=89", 0.74, 0.846,txres)
  gsn_text_ndc(wks, "(c) day=89", 0.26, 0.565,txres)
  gsn_text_ndc(wks, "(d) day=40", 0.74, 0.565,txres)


  frame(wks)
  system("convert -trim -density 300 -colorspace RGB  -define profile:skip='* '" + out_name + ".ps" + " " + out_name + ".jpg")
end 

