begin
hour=24

fn = new((/3,4/), "string")
nlev=(/"60", "120", "240", "480"/)

do ilev = 0, dimsizes(nlev)-1
  fn(0,ilev) = "./G6L"+nlev(ilev)+"_dt400"+"_E/3d_latlon3_"+hour+"h.nc" 
  fn(1,ilev) = "./G6L"+nlev(ilev)+"_dt400"+"_IE/3d_latlon3_"+hour+"h.nc" 
  fn(2,ilev) = "./G6L"+nlev(ilev)+"_dt400"+"_I/3d_latlon3_"+hour+"h.nc" 
end do

tracer60  = new((/3, 60, 181/), "double")
tracer120 = new((/3,120, 181/), "double")
tracer240 = new((/3,240, 181/), "double")
tracer480 = new((/3,480, 181/), "double")

lambda = 180
do in = 0, 2
    fi = addfile(fn(in,0), "r")
    mxrt = fi->tracerMxrt(:,lambda,:)
    tracer60(in,:,:) = mxrt(nlev|:,lat|:)
    tracer60(in,:,:) = where(tracer60(in,:,:) .lt. 10.E-10, 0, tracer60(in,:,:))
end do
delete(mxrt)
do in = 0, 2
    fi = addfile(fn(in,1), "r")
    mxrt = fi->tracerMxrt(:,lambda,:)
    tracer120(in,:,:) = mxrt(nlev|:,lat|:)
    tracer120(in,:,:) = where(tracer120(in,:,:) .lt. 10.E-10, 0, tracer120(in,:,:))
end do
delete(mxrt)
do in = 0, 2
    fi = addfile(fn(in,2), "r")
    mxrt = fi->tracerMxrt(:,lambda,:)
    tracer240(in,:,:) = mxrt(nlev|:,lat|:)
    tracer240(in,:,:) = where(tracer240(in,:,:) .lt. 10.E-10, 0, tracer240(in,:,:))
end do
delete(mxrt)
do in = 1, 2
    fi = addfile(fn(in,3), "r")
    mxrt = fi->tracerMxrt(:,lambda,:)
    tracer480(in,:,:) = mxrt(nlev|:,lat|:)
    tracer480(in,:,:) = where(tracer480(in,:,:) .lt. 10.E-10, 0, tracer480(in,:,:))
end do

wks = gsn_open_wks("ps", "DCMIP1-2_G6_dt400_"+hour+"h")

gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
;************

res            = True
res@gsnDraw    = False
res@gsnFrame   = False
res@gsnCenterStringFontHeightF = 0.0128
res@gsnLeftStringFontHeightF   = 0.0128
res@gsnRightStringFontHeightF  = 0.0128
plot_xsize     = 0.25
plot_ysize     = 0.15
plot_xspace    = 0.28
plot_yspace    = 0.17
plot_xstart    = 0.1
plot_ystart    = 0.97

res@vpWidthF   = plot_xsize
res@vpHeightF  = plot_ysize
res@lbLabelBarOn = False
res@trYReverse = True
res@cnFillOn   = True
res@cnLinesOn  = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn  = False
res@cnMinLevelValF = 0.1
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels   = fspan(0, 1, 11)
res@tmYLMode       = "Explicit"
res@tmYLValues     = (/19,39,59/) ;ispan(1,60,9)
res@tmYLLabels     = (/20,40,60/) ;ispan(1,60,9)
res@tmYLLabelsOn   = True
res@tmYLMinorOn    = True
res@tmYLMinorValues= (/09,29,49/) ;ispan(4,119,5) 
res@tmXBLabelFontHeightF = 0.0125
res@tmYLLabelFontHeightF = 0.0128
res@gsnLeftString = ""
res@tiYAxisFont    = 0.0128
res@tiYAxisString  = "Model Level Index"

res@gsnLeftString  = ""
res@gsnCenterString = "Fully explicit"
res@gsnRightString  = ""
res@tmXBLabelsOn   = False
res@vpXF = plot_xstart
res@vpYF = plot_ystart 

plot0 = gsn_csm_contour(wks, tracer60(0,:,:), res)
draw(plot0)

res@gsnCenterString= "Adaptively-implicit"
res@vpXF = plot_xstart + plot_xspace
res@vpYF = plot_ystart 
res@tmYLLabelsOn   = False
res@tmYLMinorOn    = True
res@tiYAxisString  = ""
plot1 = gsn_csm_contour(wks, tracer60(1,:,:), res)
draw(plot1)


res@gsnCenterString= "Fully implicit"
res@vpXF = plot_xstart + plot_xspace + plot_xspace
res@vpYF = plot_ystart
plot2 = gsn_csm_contour(wks, tracer60(2,:,:), res)
draw(plot2)

res@gsnCenterString = ""
res@vpXF = plot_xstart
res@vpYF = plot_ystart - plot_yspace
res@tmYLLabelsOn = True
res@tmYLValues = (/39, 79, 119/)
res@tmYLLabels = (/40, 80, 120/)
res@tmYLMinorValues = (/19, 59, 99/)
res@tiYAxisString = "Model Level Index"
plot3 = gsn_csm_contour(wks, tracer120(0,:,:), res)
draw(plot3)
res@vpXF = plot_xstart + plot_xspace
res@vpYF = plot_ystart - plot_yspace
res@tmYLLabelsOn = False
res@tiYAxisString = ""
plot4 = gsn_csm_contour(wks, tracer120(1,:,:), res)
draw(plot4)
res@vpXF = plot_xstart + 2 * plot_xspace
res@vpYF = plot_ystart - plot_yspace
plot5 = gsn_csm_contour(wks, tracer120(2,:,:), res)
draw(plot5)

res@tmXBLabelsOn    = True
res@gsnCenterString = ""
res@vpXF = plot_xstart
res@vpYF = plot_ystart - 2 * plot_yspace
res@tmYLLabelsOn  = True
res@tmYLValues = (/79, 159, 239/)
res@tmYLLabels = (/80, 160, 240/)
res@tmYLMinorValues = (/39, 119, 199/)
res@tiYAxisString = "Model Level Index"
plot6 = gsn_csm_contour(wks, tracer240(0,:,:), res)
draw(plot6)

res@tmXBLabelsOn    = False
res@vpXF = plot_xstart + plot_xspace
res@vpYF = plot_ystart - 2 * plot_yspace
res@tmYLLabelsOn = False
res@tiYAxisString = ""
plot7 = gsn_csm_contour(wks, tracer240(1,:,:), res)
draw(plot7)

res@vpXF = plot_xstart + 2 * plot_xspace
res@vpYF = plot_ystart - 2 * plot_yspace
plot8 = gsn_csm_contour(wks, tracer240(2,:,:), res)
draw(plot8)

res@gsnCenterString = ""
res@tmXBLabelsOn    = True
res@vpXF = plot_xstart
res@vpYF = plot_ystart - 3 * plot_yspace
; plot9 no solution

res@vpXF = plot_xstart + plot_xspace
res@vpYF = plot_ystart - 3 * plot_yspace
res@tmYLLabelsOn  = True
delete([/res@tmYLValues, res@tmYLLabels, res@tmYLMinorValues, res@tiYAxisString/])
res@tmYLValues = (/119, 239, 359, 479/)
res@tmYLLabels = (/120, 240, 360, 480/)
res@tmYLMinorValues = (/59, 179, 299, 419/)
;res@tiYAxisString = "Model Level Index"
;res@tiYAxisString = ""
plot10 = gsn_csm_contour(wks, tracer480(1,:,:), res)
draw(plot10)

res@lbLabelBarOn     = True
res@pmLabelBarOrthogonalPosF  = 0.003
res@pmLabelBarParallelPosF    = -0.6
res@pmLabelBarHeightF         = 0.05
res@pmLabelBarWidthF          = 0.6
res@lbLabelFontHeightF        = 0.015
res@vpXF = plot_xstart + 2 * plot_xspace
res@vpYF = plot_ystart - 3 * plot_yspace
res@tmYLLabelsOn = False
plot11 = gsn_csm_contour(wks, tracer480(2,:,:), res)
draw(plot11)

delete([/res@pmLabelBarOrthogonalPosF, res@pmLabelBarParallelPosF, res@pmLabelBarHeightF, res@pmLabelBarWidthF/])

txres = True
txres@txFontHeightF = 0.015
txres@txFont = 7
txres@txJust        = "CenterCenter"
;gsn_text_ndc(wks, "nlev=60" , 0.05, 0.9, txres)
;gsn_text_ndc(wks, "nlev=120", 0.05, 0.73, txres)
;gsn_text_ndc(wks, "nlev=240", 0.05, 0.56, txres)
;gsn_text_ndc(wks, "nlev=480", 0.05, 0.39, txres)
gsn_text_ndc(wks, "No solution", 0.2, 0.35, txres)
txres@txAngleF = 90
txres@txFontHeightF = 0.008
gsn_text_ndc(wks, "Model Level Index", 0.32, 0.38, txres)

frame(wks)
end
