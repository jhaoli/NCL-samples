; plotting procedure
; ploting the cross section of 2D array(lev, lat)
;
undef("contour_zm_linearp")
procedure contour_zm_linearp(wks, data[*][*]:numeric, ci[1]:numeric, \
         cmin[1]:numeric, cmax[1]:numeric, plots[*]:graphic, \
         iplot[1]:numeric, ylabelon[1]:logical, xlabelon[1]:logical) ;, colors[*][*]:numeric)

begin
  ; set up contour levels
  ncontours = toint((cmax-cmin) / ci) +1
  levels = ispan(0, ncontours-1, 1)
  flevels = tofloat(levels)
  flevels = cmin + flevels * ci

  ; place line at every second contour, exclude zero
;  lconts = new(ncontours, string)
;  lconts(:) = "Noline"
;  izero = ncontours / 2
;  lconts(izero::2) = "LineAndLabel"
;  lconts(izero:0:-2) = "LineAndLabel"

  ; set up plotting resources
  res               = True
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@trYReverse    = True
  res@sfYArray      = data&lev
  res@gsnYAxisIrregular2Linear = True
  res@cnFillOn      = True
;  res@cnFillPalette = "GMT_panoply"
  res@gsnSpreadColors = True
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels      = flevels
  res@cnLineLabelFontHeightF = 0.015
  res@cnInfoLabelOn = False
  ; contour every second line
  res@cnMonoLevelFlag = False
;  res@cnLevelFlags  = lconts
  res@gsnContourZeroLineThicknessF = 0.0
  res@gsnContourNegLineDashPattern  = 2
  res@cnLineThicknesses = 2.
  res@gsnLeftStringFontHeightF  = 0.022
  res@gsnRightStringFontHeightF = 0.022
  ; contour labelling
  res@cnLineLabelsOn  = True
  res@cnExplicitLineLabelsOn = False
  res@cnLineLabelDensityF  = 4
  ; controlling label bar
  res@lbLabelBarOn          = False
  res@lbLabelAutoStride     = False
  res@lbBoxLinesOn          = True
  res@lbTitlePosition       = "Bottom"
  res@lbTitleFontThicknessF = 0.02
  res@lbLabelFontHeightF    = 0.015
  res@pmLabelBarHeightF     = 0.1
  res@lbBoxMinorExtentF     = 0.2
  res@lbBoxSeparatorLinesOn = False
  ; ticklabels
  res@tmYLMode         = "Explicit"
  res@tmYLValues       = fspan(0.1, 0.9, 5)
  res@tmYLLabels       = sprintf("%3.1f", fspan(0.1, 0.9, 5))
  res@tmYLMinorOn      = True
;  res@tmYLTickStartF   = 0.1
;  res@tmYLTickSpacingF = 0.2
;  res@tmYLPrecision    = 2
  res@tiYAxisString    = ""
  res@tmYLLabelsOn     = ylabelon 
  res@tmYLLabelFontHeightF = 0.03
  res@tmXBLabelsOn     = xlabelon
  res@tmXBLabelFontHeightF = 0.03
  res@tmXBValues       = (/-90, -60, -30, 0, 30, 60, 90/)
  res@tmXBLabels       = (/"90S", "60S", "30S", "0", "30N", "60N", "90N"/)
  ;
  plots(iplot) = gsn_csm_contour(wks, data, res)
  return(plots)
end
