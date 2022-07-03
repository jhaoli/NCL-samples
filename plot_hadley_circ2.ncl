begin
nlev=480
method="E"
fn0 = "./G6L"+nlev+"_dt600_"+method+"/3d_latlon3_12h.nc" 
fn1 = "./G6L"+nlev+"_dt600_"+method+"/3d_latlon3_24h.nc" 

fi0 = addfile(fn0, "r")
fi1 = addfile(fn1, "r")

tracer0 = fi0->tracerMxrt
tracer1 = fi1->tracerMxrt

tracer_reorder0 = tracer0(nlev|:,lat|:,lon|:)
tracer_reorder1 = tracer1(nlev|:,lat|:,lon|:)

wks = gsn_open_wks("ps", "DCMIP1-2.G5L"+nlev+"_dt600_"+method)

gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
;************

res            = True
res@gsnDraw    = False
res@gsnFrame   = False
res@gsnCenterStringFontHeightF = 0.0128
res@gsnLeftStringFontHeightF   = 0.0128
res@gsnRightStringFontHeightF  = 0.0128
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
res@tmXBLabelsOn    = False
res@vpXF = plot_xstart
res@vpYF = plot_ystart 
tracer_reorder0 = where(tracer_reorder0 .lt. 10.E-10, 0,  tracer_reorder0)
plot0 = gsn_csm_contour(wks, tracer_reorder0(:,:,90), res)


res@gsnLeftString  = "(b)"
res@gsnCenterString= "G5L"+nlev
res@gsnRightString = "24h"
res@tmYLLabelsOn   = True
res@tiYAxisString  = "Model Level Index"
res@tmXBLabelsOn   = True
res@vpXF = plot_xstart ;+ plot_xspace
res@vpYF = plot_ystart - plot_yspace


res@lbLabelBarOn     = True
;res@pmLabelBarOrthogonalPosF  = 0.0003
res@pmLabelBarParallelPosF    = 1.1
res@pmLabelBarHeightF         = 0.4
res@pmLabelBarWidthF          = 0.06
res@lbLabelFontHeightF        = 0.01
res@lbOrientation             = "Vertical"
tracer_reorder1 = where(tracer_reorder1 .lt. 10.E-10, 0,  tracer_reorder1)
plot1 = gsn_csm_contour(wks, tracer_reorder1(:,:,90), res)

draw(plot0)
draw(plot1)
frame(wks)
end

