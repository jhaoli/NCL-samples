
function wrap_plot_xy(wks, res, tiXString, tiYString, LString, CString, RString, error_norm)
begin
  res@tiXAxisFontHeightF  = 0.0128
  res@tiYAxisFontHeightF  = 0.0128
  res@tiXAxisString       = tiXString
  res@tiYAxisString       = tiYString
  res@gsnLeftString       = LString
  res@gsnCenterString     = CString
  res@gsnRightString      = RString
  plot                    = gsn_csm_xy(wks, fspan(1,3,3), log10(error_norm), res)
  return(plot)
end 

function l1_norm(comp_solution, true_solution, delz, area)
; comp_solution(nv, nlev), true_solution(nv,nlev)
; delz(nv,nlev), area(nv)
begin
  area2d = conform(delz,area,0) 
  tmp1 = area2d * abs(comp_solution - true_solution) * delz
  tmp2 = area2d * abs(true_solution) * delz
  l1   = sum(tmp1) / sum(tmp2)
  return(l1)
end

function l2_norm(comp_solution, true_solution, delz, area)
; comp_solution(nv, nlev), true_solution(nv,nlev)
; delz(nv,nlev), area(nv)
begin
  area2d = conform(delz,area,0)
  tmp1 = area2d * (comp_solution - true_solution)^2 * delz
  tmp2 = area2d * true_solution^2 * delz
  l2  = sqrt(sum(tmp1) / sum(tmp2))
  return(l2)
end

function linf_norm(comp_solution, true_solution, delz, area)
; comp_solution(nv, nlev), true_solution(nv,nlev)
; delz(nv,nlev), area(nv)
begin
  diff = abs(comp_solution - true_solution)
  linf = max(diff) / max(abs(true_solution))
  return(linf)
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin
  resolution = (/"5"  , "6"  , "7"  /)
  nlev       = (/"60" , "120", "240"/)
  timestep   = (/"800", "400", "200"/)
  scheme     = (/"E" , "IE"/)  
  error_norm = new((/6,dimsizes(resolution)/), "double")

  do ir = 0, 2
  do is = 0, 1
    file_path = "./G"+resolution(ir)+"L"+nlev(ir)+"_dt"+timestep(ir)+"_"+scheme(is)+"/"
    f1d_24h = file_path+"GRIST.ATM.G"+resolution(ir)+".tracer.00000-00000.1d.h1.nc"
    f2d_24h = file_path+"GRIST.ATM.G"+resolution(ir)+".tracer.00001-00000.2d.h1.nc"
    f3d_24h = file_path+"GRIST.ATM.G"+resolution(ir)+".tracer.00001-00000.3d.h1.nc"
    f3d_00h = file_path+"GRIST.ATM.G"+resolution(ir)+".tracer.00000-00000.3d.h1.nc"
    
    filein1d = addfile(f1d_24h, "r")
    areaPC = filein1d->areaPC
    
    filein2d = addfile(f2d_24h, "r")
    hp_face = filein2d->hPressureFace
    hp_full = filein2d->hPressureFull
    density = hp_full / (287.d * 300.d)
    delz = (hp_face(:,1:toint(nlev(ir))) - hp_face(:,0:toint(nlev(ir))-1)) / (density * 9.8d)
    copy_VarCoords(hp_full, delz)
    
    filein3d = addfile(f3d_24h, "r")
    tracer = filein3d->tracerMxrt(:,:,0)
    
    filein3d0 = addfile(f3d_00h, "r")
    tracer0 = filein3d0->tracerMxrt(:,:,0)
    if (is .eq. 0) then
      error_norm(0,ir)  = l1_norm(tracer, tracer0, delz, areaPC)
      error_norm(1,ir)  = l2_norm(tracer, tracer0, delz, areaPC)
      error_norm(2,ir)  = linf_norm(tracer, tracer0, delz, areaPC)
    else 
      error_norm(3,ir)  = l1_norm(tracer, tracer0, delz, areaPC)
      error_norm(4,ir)  = l2_norm(tracer, tracer0, delz, areaPC)
      error_norm(5,ir)  = linf_norm(tracer, tracer0, delz, areaPC)
    end if 
    delete([/areaPC, hp_face, hp_full, density, delz, tracer, tracer0/])
  end do
  end do
; plot convergence order
  wks    = gsn_open_wks("ps", "dcmip12.convergence_rate")

  res            = True
  res@gsnFrame   = False
  res@gsnDraw    = False
  res@vpHeightF  = 0.5
  res@vpWidthF   = 0.5
  res@trXMinF    = 0.9
  res@trXMaxF    = 3.1
  res@trYMinF    = -2.0
  res@trYMaxF    = -0.2
  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/1,2,3/)
  res@tmXBLabels = (/"G5L60(800s)", "G6L120(400s)", "G7L240(200s)"/)
  res@tmXBLabelFontHeightF = 0.0128

  res@tmYLMode   = "Explicit"
  res@tmYLLabelFontHeightF = 0.0128
  res@gsnCenterStringFontHeightF = 0.0128
  res@gsnLeftStringFontHeightF   = 0.0128
  res@gsnRightStringFontHeightF  = 0.0128

  plot_xsize = 0.5
  plot_ysize = 0.5
  plot_xspace = 0.41
  plot_yspace = 0.41
;  plot_xstart = 0.15
;  plot_ystart = 0.95

;  res@vpHeightF = plot_ysize
;  res@vpWidthF  = plot_xsize
;  res@vpXF      = plot_xstart
;  res@vpYF      = plot_ystart 
; legend
  res@pmLegendDisplayMode    = "Always"
  res@lgJustification        = "TopRight"
  res@pmLegendOrthogonalPosF = -0.56 ; negative to upper
  res@pmLegendParallelPosF   = 0.52 ; large to right
  res@pmLegendWidthF         = 0.3
  res@pmLegendHeightF        = 0.2
  res@lgBoxMinorExtentF      = 0.16
  res@lgLabelFontHeightF     = 0.0128
  res@lgPerimOn              = False
; xy line
  res@xyMonoDashPattern  = False
  res@xyDashPatterns     = (/0,1,2,0,1,2/)
  res@xyMarkLineMode     = "Marklines"
  res@xyLineColors       = (/"royalblue1", "red", "blue", "green", "black", "Cyan"/)
  res@xyLineThicknessF   = 3.8
  res@xyMarkerSizes      = (/0.01,0.01,0.01, 0.01,0.01, 0.01/)
  res@xyMarkers          = (/16, 6, 2, 16, 6, 2/)
  res@xyExplicitLabels   = (/"Fully-EXP:L1"     , "Fully-EXP:L2"   , "Fully-EXP:Linf", \\
                            "Adaptively-Imp:L1", "Adaptive-Imp:L2", "Adaptively-Imp:Linf"/)
  plot0 = wrap_plot_xy(wks, res, "", "Log10", "", "", "", error_norm)

; convergence line: 1st order
  res_lines = True
  res_lines@gsLineDashPattern = 0.
  res_lines@gsLineThicknessF  = 3.9
  res_lines@gsLineColor       = "gray"
  
  k  = 1
  y1 = -0.2
  y2 = y1 - k * log10(4)
  xx = (/1,3/)
  yy = (/y1,y2/)
  dum1 = gsn_add_polyline(wks, plot0, xx, yy, res_lines)
  
  k  = 2
  y1 = -0.2
  y2 = y1 - k * log10(4)
  xx = (/1,3/)
  yy = (/y1,y2/)
  dum2 = gsn_add_polyline(wks, plot0, xx, yy, res_lines)
  draw(plot0)

end
