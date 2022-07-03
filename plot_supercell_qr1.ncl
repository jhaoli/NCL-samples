begin

  datapath = "./g4-hadv33-hnrk3-vadv3-vnrk3_L80_d5_t10_m10_E/"
  glevel = "G4" 

  f3d = addfile(datapath+"GRIST.ATM."+glevel+".dtp.00000-07200.3d.h1.nc", "r")
  f2d = addfile(datapath+"GRIST.ATM."+glevel+".dtp.00000-07200.2d.h1.nc", "r")
  f1d = addfile(datapath+"GRIST.ATM."+glevel+".dtp.00000-07200.1d.h1.nc", "r")
  
  qr      = f3d->tracerMxrt(:,:,2)*1.e3 ; location_nv, nlev, ntracer(qv,qc,qr)
  gh_face = f2d->geopotentialFace / 9.80616
  lon_nv  = f1d->lon_nv
  lat_nv  = f1d->lat_nv
  nlevp   = f1d->nlevp
  nface   = dimsizes(nlevp)
  gh_full = (gh_face(:,0:nface-2) + gh_face(:,1:nface-1)) * 0.5d0

; 1st interpolation function
  linlog = 1 ; linear interpolation
  qr_5km = int2p_Wrap(gh_full, qr, (/5000/), linlog)
;  printVarSummary(qr_5km)

; 2nd interpolation function
;  qr_5 = new((/ncell/), "double")
;  iv = 0
;  do while(iv .le. ncell-1)
;    if (lon_nv(iv)>-10 .and. lon_nv(iv) < 20 .and. \
;        lat_nv(iv)>-40 .and. lat_nv(iv) < 40) then
;        qr_5(iv) = wrf_interp_1d(qr(iv,:), gh_full(iv,:), (/5000/))
;    end if 
;    iv = iv + 1
;  end do
;printVarSummary(qr_5)

;--- Start the graphic
  wks = gsn_open_wks("ps", "supercell_G4L80d5t10m10_E")
;  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  
  res                      = True
  res@gsnFrame             = False
  res@gsnDraw              = False
  res@cnFillOn             = True
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnInfoLabelOn        = False
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = 1
  res@cnMaxLevelValF       = 12
  res@cnLevelSpacingF      = 1

  res@mpMinLonF  = -15
  res@mpMaxLonF  = 30
  res@mpMinLatF  = -40
  res@mpMaxLatF  = 40
  res@mpGeophysicalLineThicknessF = 1.2

  res@trGridType   = "TriangularMesh"

  res@tmBorderThicknessF   = 4.0
  res@tmYLMinorOn          = True
  res@tmYLLabelFontHeightF = 0.02
  res@tmYLMinorLengthF     = 0.01
  res@tmYLMinorOutwardLengthF = 0.01
  res@tmYLMinorThicknessF  = 3.0
  res@tmYLMajorThicknessF  = 5.0

  res@tmXBMinorOn          = False
  res@tmXBLabelFontHeightF = 0.02
  res@tmXBMinorLengthF     = 0.01
  res@tmXBMinorOutwardLengthF = 0.01
  res@tmXBMinorThicknessF  = 3.0
  res@tmXBMajorThicknessF  = 5.0
  res@gsnAddCyclic         = False

  res@sfXArray = lon_nv
  res@sfYArray = lat_nv
  
  res@lbLabelBarOn      = True
  res@lbBottomMarginF   = 0.05
  res@lbRightMarginF    = 0.02
  res@lbOrientation  = "Vertical"

;  plot = gsn_csm_contour_map(wks,qr(:,29),res)
  plot = gsn_csm_contour_map(wks,qr_5km(:,0),res)
;  plot = gsn_csm_contour_map(wks,qr_5(:),res)

draw(plot)
frame(wks)
end
