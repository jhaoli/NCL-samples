load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  file_gamil = "../gamil-c-grid-pole-u/cp/cp_c_u_01.180x91.d20.nc"
  file_gmcore_u = "../gmcore-pole-u/cp/cp_c_u_01.180x91.dt20-midpoint.nc"
  file_gmcore_u_stokes = "../gmcore-pole-u/cp/cp_c_u_01.180x91.dt20-stokes.nc"
  file_gmcore_v = "../gmcore-pole-v/cp/cp_c_v_01.180x90.dt10-midpoint.nc"

  f_gamil = addfile(file_gamil, "r")
  u_gamil = f_gamil->u(10,:,:)

  f_gmcore_u = addfile(file_gmcore_u, "r")
  u_gmcore_u = f_gmcore_u->u(10,:,:)

  f_gmcore_u_stokes = addfile(file_gmcore_u_stokes, "r")
  u_gmcore_u_stokes = f_gmcore_u_stokes->u(10,:,:) 

  f_gmcore_v = addfile(file_gmcore_v, "r")
  u_gmcore_v = f_gmcore_v->u(10,:,:)
; create plots 
;**************************
  out_name = "cross_pole_case"
  wks = gsn_open_wks("eps", out_name)
  gsn_define_colormap(wks, "GMT_panoply")

  plot = new(4, graphic)
  
  res                   = True 
  res@gsnDraw           = False 
  res@gsnFrame          = False 
  ; res@gsnMaximize       = True 
  res@gsnPolar          = "NH"
  res@cnFillOn          = True 
  res@cnLinesOn         = True
  res@cnLineLabelsOn    = False
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF    = -20
  res@cnMaxLevelValF    = 20
; res@cnLevels          = (/-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120/)
  res@cnLineThicknessF  = 1.2
  res@gsnContourNegLineDashPattern = 1
  res@lbLabelBarOn      = False 

; res@mpDataBaseVersion = "MediumRes"
  res@mpGridAndLimbOn   = True 
  res@mpGridLineColor   = 1
  res@mpPerimOn         = False

; res@mpFillDrawOrder   = "PreDraw"

  res@gsnRightString    = ""
  res@gsnLeftString     = ""
; res@tiMainFontHeightF = 0.02

; 
  mpres = res 
  mpres@mpMinLatF       = 30
  ; mpres@mpProjection   = "Orthographic"
  ; mpres@mpCenterLatF   = 90
  ; mpres@mpCenterLonF   = 180
  mpres@mpOutlineOn      = False
  mpres@gsnPolarLabelFontHeightF=0.016

  plot(0) = gsn_csm_contour_map_polar(wks, u_gamil, mpres)
  plot(1) = gsn_csm_contour_map_polar(wks, u_gmcore_u, mpres)
  plot(2) = gsn_csm_contour_map_polar(wks, u_gmcore_u_stokes, mpres)
  plot(3) = gsn_csm_contour_map_polar(wks, u_gmcore_v, mpres)

;***
; create panel
;***
  resP                     = True                ; modify the panel plot
  resP@gsnFrame            = False
  ; resP@gsnPanelMainString = "A plot with a common label bar"
  resP@gsnPanelFigureStrings = (/"(a)", "(b)", "(c)", "(d)"/)
  resP@amJust              = "TopLeft"
  resP@gsnPanelFigureStringsPerimOn = False
  resP@gsnPanelFigureStringsFontHeightF = 0.02
  resP@gsnPanelLabelBar    = True                ; add common colorbar
  resP@lbLabelFontHeightF  = 0.014               ; make labels smaller
  resP@lbTopMarginF      = 0.1
  resP@lbTitleOffsetF    = -0.2
  ; resP@lbTitleString     = "u wind component (m s~S~-1~N~)"
  ; res@pmLabelBarHeightF = 0.05 ;control the height of labelbar
  ; res@lbTitleFontHeightF=0.015
; res@lbLabelFontHeightF= 0.01 ; control the size of the labels
  
  gsn_panel(wks, plot, (/2,2/), resP)
  frame(wks)
  system("convert -trim -density 300 " + out_name + ".eps" + " " + out_name + ".png")
end 

