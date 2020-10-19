load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

  file_gmcore_v = "../gmcore-pole-v/mz/data/state_mz_180x90_dt15_midpoint_classic_pole_v.nc"
  file_gmcore_u_apvm = "../gmcore-pole-u/mz/data/mz_c_u_01.180x91.dt30-APVM.nc"
  file_gmcore_u_stokes = "../gmcore-pole-u/mz/data/mz_c_u_01.180x91.dt30-stokes.nc"

  f_v = addfile(file_gmcore_v, "r")
  u_v = f_v->u(20,:,:)

  f_apvm = addfile(file_gmcore_u_apvm, "r")
  u_apvm = f_apvm->u(20,:,:)

  f_stokes = addfile(file_gmcore_u_stokes, "r")
  u_stokes = f_stokes->u(20,:,:) 


; create plots 
;**************************
  out_name = "mountain_wave_case_gmcore"
  wks = gsn_open_wks("eps", out_name)
  gsn_define_colormap(wks, "GMT_panoply")

  plot = new(6, graphic)
  
  res                   = True 
  res@gsnDraw           = False 
  res@gsnFrame          = False 
  ; res@gsnMaximize       = True 
  res@gsnMajorLonSpacing=60
  res@gsnMajorLatSpacing=30
  res@gsnMinorLonSpacing=10
  res@gsnMinorLatSpacing=10
  res@cnFillOn          = True 
  res@cnLinesOn         = True
  res@cnLineLabelsOn    = False
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF    = -16
  res@cnMaxLevelValF    = 38
; res@cnLevels          = (/-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120/)
  res@cnLineThicknessF  = 1.2
  res@gsnContourNegLineDashPattern = 1
  res@lbLabelBarOn      = False 
  res@tmYLLabelFontHeightF = 0.012
  res@tmYLLabelFontThicknessF = 1.0
  res@tmXBLabelFontThicknessF = 0.012
  res@tmYROn            = False
  res@tmXTOn            = False

; res@mpDataBaseVersion = "MediumRes"
  ; res@mpGridAndLimbOn   = True 
  ; res@mpGridLineColor   = 1
  res@mpPerimOn         = False
; res@mpFillDrawOrder   = "PreDraw"

  res@gsnRightString    = ""
  res@gsnLeftString     = ""
; res@tiMainFontHeightF = 0.02

; 
  res0 = res 
  res0@vpXF             = 0.05
  res0@vpYF             = 0.95
  res0@vpWidthF         = 0.5
  res0@vpHeightF        = 0.3
  res0@mpOutlineOn      = False
  plot(0) = gsn_csm_contour_map(wks, u_v, res0)

  res1 = res 
  res1@gsnPolar        = "NH" 
  res1@vpXF            = 0.6
  res1@vpYF            = 0.93
  res1@vpWidthF        = 0.26
  res1@vpHeightF       = 0.26
  res1@mpOutlineOn      = False
  res1@mpMinLatF      = 60
  res1@mpCenterLonF   = 0
  res1@gsnPolarLabelFontHeightF=0.01

  plot(1) = gsn_csm_contour_map_polar(wks, u_v, res1)

  res0@vpYF             = 0.65
  plot(2) = gsn_csm_contour_map(wks, u_apvm, res0)

  res1@vpYF             = 0.63
  res1@lbLabelBarOn     = True                ; add common colorbar
  res1@lbOrientation    = "Vertical"
  res1@lbLabelFontHeightF  = 0.012               ; make labels smaller
  res1@pmLabelBarHeightF  = 0.5
  res1@pmLabelBarWidthF   = 0.1

  plot(3) = gsn_csm_contour_map_polar(wks, u_apvm, res1)

  res0@vpYF             = 0.35
  plot(4) = gsn_csm_contour_map(wks, u_stokes, res0)

  res1@vpYF             = 0.33
  res1@lbLabelBarOn     = False
  plot(5) = gsn_csm_contour_map_polar(wks, u_stokes, res1)
  draw(plot)
  ;************************************************
; Attach text to plot using plot coordinates.
;************************************************
  txres               = True                     ; text mods desired
  txres@txPerimOn     = False
  txres@txBackgroundFillColor = "white"
  txres@txFontHeightF = 0.018                     ; font smaller. default big
  txres@txFontThicknessF = 0.017

  gsn_text_ndc(wks,"(a)",0.074,0.903,txres)
  gsn_text_ndc(wks,"(c)",0.074,0.603,txres)
  gsn_text_ndc(wks,"(e)",0.074,0.303,txres)

  gsn_text_ndc(wks,"(b)",0.605,0.903,txres)
  gsn_text_ndc(wks,"(d)",0.605,0.603,txres)
  gsn_text_ndc(wks,"(f)",0.605,0.303,txres)

  frame(wks)
  exit()
;***
; create panel
;***
  resP                     = True                ; modify the panel plot
  resP@gsnFrame            = False
  ; resP@gsnMaximize         = False
  ; resP@gsnPanelMainString = "A plot with a common label bar"
  resP@gsnPanelFigureStrings = (/"(a)", "(b)", "(c)", "(d)", "(e)", "(f)"/)
  resP@amJust              = "TopLeft"
  resP@gsnPanelFigureStringsPerimOn = False
  resP@gsnPanelFigureStringsFontHeightF = 0.015
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  resP@lbOrientation       = "Vertical"
  ; resP@lbLabelFontHeightF  = 0.01               ; make labels smaller
  ; resP@lbRightMarginF      = 0.1
  ; resP@lbTitleOffsetF    = -0.2
  ; resP@lbTitlePosition     = "Left"
  ; resP@lbTitleDirection    = "Top"
  ; resP@lbTitleString     = "(m s~S~-1~N~)"
  ; resP@pmLabelBarOrthogonalPosF = 0.5
  ; res@pmLabelBarHeightF = 0.05 ;control the height of labelbar
  ; res@lbTitleFontHeightF=0.015
; res@lbLabelFontHeightF= 0.01 ; control the size of the labels
; resP@gsnPanelYWhiteSpacePercent = 0
 ; resP@gsnPanelXWhiteSpacePercent = 0
  
  ; gsn_panel(wks, plot, (/3,2/), resP)
  frame(wks)
  system("convert -trim -density 300 " + out_name + ".eps" + " " + out_name + ".png")
end 

