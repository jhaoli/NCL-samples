begin
nlev = 60
fn0 = "./G6L60_dt900_E/3d_latlon_q3_6.nc"
fn1 = "./G6L60_dt900_E/3d_latlon_q3_12.nc"

fi0 = addfile(fn0, "r")
fi1 = addfile(fn1, "r")

tracer0 = fi0->tracerMxrt
tracer1 = fi1->tracerMxrt

tracer0_reorder = tracer0(nlev|:, lat|:, lon|:)
tracer1_reorder = tracer1(nlev|:, lat|:, lon|:)

wks = gsn_open_wks("ps", "DCMIP1-1G6L60")

gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")


res            = True
res@gsnDraw    = False
res@gsnFrame   = False
res@gsnCenterStringFontHeightF = 0.0128
res@gsnLeftStringFontHeightF   = 0.0128
res@gsnRightStringFontHeightF  = 0.0128
res@mpCenterLonF = 180
plot_xsize     = 0.4
plot_ysize     = 0.2
plot_xspace    = 0.44
plot_yspace    = 0.24
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
;res@tmYLMode       = "Explicit"
;res@tmYLValues     = ispan(19,119,20)
;res@tmYLLabels     = ispan(20,120,20)
res@tmYLLabelsOn   = True
res@tmYLMinorOn    = True
;res@tmYLMinorValues= ispan(4,119,5)
res@tmXBLabelFontHeightF = 0.0128
res@tmYLLabelFontHeightF = 0.0128
res@gsnLeftString = ""
res@tiYAxisFont    = 0.0128
res@tiYAxisString  = "Model Level Index"

res@gsnLeftString   = "(a)"
res@gsnCenterString = "G5L"+nlev
res@gsnRightString  = "12h"
res@tmXBLabelsOn    = True
res@vpXF = plot_xstart
res@vpYF = plot_ystart
tracer0_reorder = where(tracer0_reorder .lt. 10.E-10, 0,  tracer0_reorder)
plot0 = gsn_csm_contour_map(wks, tracer0_reorder(35,:,:), res)

res@gsnLeftString  = "(b)"
res@gsnCenterString= "G5L"+nlev
res@gsnRightString = "24h"
res@tmYLLabelsOn   = False
res@tiYAxisString  = ""
res@vpXF = plot_xstart + plot_xspace
res@vpYF = plot_ystart
res@lbLabelBarOn     = True
res@pmLabelBarOrthogonalPosF  = 0.0003
res@pmLabelBarParallelPosF    = -0.05
res@pmLabelBarHeightF         = 0.05
res@pmLabelBarWidthF          = 0.6
res@lbLabelFontHeightF        = 0.01
tracer1_reorder = where(tracer1_reorder .lt. 10.E-10, 0,  tracer1_reorder)
plot1 = gsn_csm_contour_map(wks, tracer1_reorder(35,:,:), res)

draw(plot0)
draw(plot1)
frame(wks)
end
