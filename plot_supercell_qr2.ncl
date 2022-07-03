begin
  
  tags = (/"t10_m10_E", "t10_m10_IE", "t10_m10_I", \
           "t20_m20_E", "t20_m20_IE", "t20_m20_I", \
           "t30_m30_E", "t40_m40_IE", "t40_m40_I"/)
  qr_5 = new((/9, 2562/), "double")
  datapath = "./g4-hadv33-hnrk3-vadv3-vnrk3_L80_d5_"
  do itag = 0 , 8
    f3d = addfile(datapath+tags(itag)+"/GRIST.ATM.G4.dtp.00000-07200.3d.h1.nc", "r")
    f2d = addfile(datapath+tags(itag)+"/GRIST.ATM.G4.dtp.00000-07200.2d.h1.nc", "r")
    f1d = addfile(datapath+tags(itag)+"/GRIST.ATM.G4.dtp.00000-07200.1d.h1.nc", "r")
  
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
;    printVarSummary(qr_5km)
    qr_5(itag,:) = qr_5km(:,0)
  end do

;  printVarSummary(qr_5)
;--- Start the graphic
  wks = gsn_open_wks("ps", "supercell_G4L80_2")
  gsn_define_colormap(wks,"WhBlGrYeRe")
  
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
  plot_xsize = 0.25
  plot_ysize = 0.3
  plot_xspace = 0.19
  plot_yspace = 0.32
  plot_xstart = 0.1
  plot_ystart = 0.95
  
  res@vpWidthF   = plot_xsize
  res@vpHeightF  = plot_ysize
  res@mpMinLonF  = -15
  res@mpMaxLonF  = 30
  res@mpMinLatF  = -40
  res@mpMaxLatF  = 40
  res@mpOutlineOn = False
;  res@mpGeophysicalLineThicknessF = 1.0
;  res@mpGeophysicalLineColor = "Background"

  res@trGridType   = "TriangularMesh"

  res@tmBorderThicknessF   = 3.0
  res@tmYLMinorOn          = True
  res@tmYLLabelFontHeightF = 0.015
;  res@tmYLMinorLengthF     = 0.01
;  res@tmYLMinorOutwardLengthF = 0.01
;  res@tmYLMinorThicknessF  = 2.0
;  res@tmYLMajorThicknessF  = 3.0

  res@tmXBMinorOn          = True
  res@tmXBLabelFontHeightF = 0.015
 ; res@tmXBMinorLengthF     = 0.01
 ; res@tmXBMinorOutwardLengthF = 0.01
 ; res@tmXBMinorThicknessF  = 2.0
 ; res@tmXBMajorThicknessF  = 3.0
  res@gsnAddCyclic         = False
  res@gsnStringFontHeightF = 0.015
  res@lbLabelBarOn      = False
  
  res@sfXArray = lon_nv
  res@sfYArray = lat_nv
  
  res@tmYLLabelsOn      = True
  res@tmXBLabelsOn      = True

  res@gsnCenterString = "Explicit"
  res@vpXF = plot_xstart 
  res@vpYF = plot_ystart 
  plot0 = gsn_csm_contour_map(wks,qr_5(0,:),res)
  draw(plot0)

  res@tmYLLabelsOn      = False
  res@tmXBLabelsOn      = False
  res@gsnCenterString = "Adaptively-implicit"
  res@vpXF = plot_xstart + plot_xspace
  res@vpYF = plot_ystart 
  plot1 = gsn_csm_contour_map(wks, qr_5(1,:), res)
  draw(plot1)
 
  res@gsnCenterString = "Fully implicit"
  res@vpXF = plot_xstart + 2*plot_xspace
  res@vpYF = plot_ystart 
  plot2 = gsn_csm_contour_map(wks, qr_5(2,:), res)
  draw(plot2)
  
  res@gsnCenterString = ""
  res@vpXF = plot_xstart 
  res@vpYF = plot_ystart - plot_yspace 
  plot3 = gsn_csm_contour_map(wks, qr_5(3,:), res)
  draw(plot3)

  res@gsnCenterString = ""
  res@vpXF = plot_xstart + plot_xspace
  res@vpYF = plot_ystart - plot_yspace 
  plot4 = gsn_csm_contour_map(wks, qr_5(4,:), res)
  draw(plot4)

  res@vpXF = plot_xstart + 2*plot_xspace
  res@vpYF = plot_ystart - plot_yspace 
  plot5 = gsn_csm_contour_map(wks, qr_5(5,:), res)
  draw(plot5)

  res@vpXF = plot_xstart 
  res@vpYF = plot_ystart - 2*plot_yspace 
  plot6 = gsn_csm_contour_map(wks, qr_5(6,:), res)
  draw(plot6)

  res@vpXF = plot_xstart + plot_xspace 
  res@vpYF = plot_ystart - 2*plot_yspace 
  plot7 = gsn_csm_contour_map(wks, qr_5(7,:), res)
  draw(plot7)


  res@tmYLLabelsOn             = False
  res@lbLabelBarOn             = True
  res@pmLabelBarOrthogonalPosF = 0.1
  res@pmLabelBarParallelPosF   = 1.6
  res@pmLabelBarHeightF        = 0.6
  res@pmLabelBarWidthF         = 0.05
  res@lbBottomMarginF          = 0.05
  res@lbRightMarginF           = 0.02
  res@lbOrientation            = "Vertical"

  res@vpXF = plot_xstart + 2*plot_xspace
  res@vpYF = plot_ystart - 2*plot_yspace  
  plot8 = gsn_csm_contour_map(wks, qr_5(8,:), res)
  
;  do i =1,8
;    print("Min:"+min(qr_5(i,:))+ "  Max:" + max(qr_5(i,:)))
;  end do
  draw(plot8)
  txres     = True
  txres@txFontHeightF      =  0.015
  txres@txJust             = "Centerleft"
  gsn_text_ndc(wks, "~F33~D~F~t=10s", 0.035, 0.8, txres)
  gsn_text_ndc(wks, "~F33~D~F~t=20s", 0.035, 0.48, txres)
  gsn_text_ndc(wks, "~F33~D~F~t=30s/40s", 0.035, 0.18, txres)
  
  txres@txFontHeightF     = 0.014
  txres@txFontColor       = "black"
;  txres@txBackgroundFillColor = "White"
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(0,:))), 0.21, 0.66, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(1,:))), 0.40, 0.66, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(2,:))), 0.59, 0.66, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(3,:))), 0.21, 0.34, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(4,:))), 0.40, 0.34, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(5,:))), 0.59, 0.34, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(6,:))), 0.21, 0.022, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(7,:))), 0.40, 0.022, txres)
  gsn_text_ndc(wks, "max:" + sprintf("%.2f", max(qr_5(8,:))), 0.59, 0.022, txres)

frame(wks)
end
