undef ("wrap_plot_map")
function wrap_plot_map(wks, res, var1, var2, var1_ini, var2_ini, lr_value, lu_value, lo_value)
begin

  res@tiXAxisString     = "q1"
  res@tiYAxisString     = "q2"
  res@xyMarkerColor     = "red"
  res@xyMarkerSizeF     = 0.0003
  plot0 = gsn_csm_xy(wks, var1, var2, res)

  
  res@xyMarkerColor     = "black"
  res@xyMarkerSizeF     = 0.0005
  plot1 = gsn_csm_xy(wks, var1_ini, var2_ini, res)
  overlay(plot0, plot1)

  res_line                   = True
  res_line@gsLineDashPattern = 0. ; solid line
  res_line@gsLineThicknessF  = 1
  res_line@gsLineColor       = "black"
  min_x                      = min(var1_ini)
  min_y                      = min(var2_ini)
  max_x                      = max(var1_ini)
  max_y                      = max(var2_ini)

  dum = new(3, graphic)
  dum(0) = gsn_add_polyline(wks, plot0, (/min_x,max_x/), (/max_y,min_y/), res_line)
  dum(1) = gsn_add_polyline(wks, plot0, (/min_x,max_x/), (/max_y,max_y/), res_line)
  dum(2) = gsn_add_polyline(wks, plot0, (/max_x,max_x/), (/max_y,min_y/), res_line)
  plot0@dum = dum

  txres                      = True
  txres@txFontHeightF        = 0.02
  text0 = gsn_add_text(wks, plot0, "lr: "+lr_value, 0.18, 0.15, txres)
  text1 = gsn_add_text(wks, plot0, "lu: "+lu_value, 0.18, 0.10, txres)
  text2 = gsn_add_text(wks, plot0, "lo: "+lo_value, 0.18, 0.05, txres)
  return(plot0)
end 

undef ("lruo_norm")
function lruo_norm(dk, phi1, phi2, area, nsize)

;-------------------
; |  \'             |
; |   \  '          |
; |    \   '  B     |
; |     \   '       |
; |      \ A '      |
; |       \   '     |
; | B      \   '    |
; |         \  '    |
; |          \ '    |
;-------------
begin
  phi1_max_ini = 1.d
  phi1_min_ini = 0.0d
  phi2_max_ini = 0.9d
  phi2_min_ini = 0.1d

  lr = 0.d0
  lu = 0.d0
  lo = 0.d0
  la = 0.d0

  do iv = 0, nsize - 1
    if(phi1(iv).ge.phi1_min_ini.and. \
       phi1(iv).le.phi1_max_ini.and. \
       phi2(iv).ge.phi2_min_ini.and. \
       phi2(iv).le.phi2_max_ini.and. \
       phi2(iv).le.(0.9d - 0.8d * phi1(iv)^2) .and. \
       phi2(iv).ge.(0.9d - 0.8d * phi1(iv))) then
       lr = lr + dk(iv) * area(iv)
    end if

    if(phi1(iv) .ge. phi1_min_ini .and. \
       phi1(iv) .le. phi1_max_ini .and. \
       phi2(iv) .ge. phi2_min_ini .and. \
       phi2(iv) .le. phi2_max_ini .and. \
       (phi2(iv) .gt.(0.9d - 0.8d * phi1(iv)^2) .or. \
        phi2(iv) .lt. (0.9d - 0.8d * phi1(iv)))) then
       lu = lu + dk(iv) * area(iv)
    end if

    if(phi1(iv) .lt. phi1_min_ini .or. \
       phi1(iv) .gt. phi1_max_ini .or. \
       phi2(iv) .lt. phi2_min_ini .or. \
       phi2(iv) .gt. phi2_max_ini ) then
       lo = lo + dk(iv) * area(iv)
    end if
    la = la + dk(iv) * area(iv)
  end do
  lr = lr / sum(area)
  lu = lu / sum(area)
  lo = lo / sum(area)
  la = la / sum(area)
  tmp = lr+lu+lo
  print("LA1="+tmp+"LA2="+la)
  return([/lr,lu,lo/])
end

undef ("compute_dk")
function compute_dk(phi1, phi2, nsize)
begin
  phi1_max_ini = 1.0
  phi1_min_ini = 0.d
  phi2_max_ini = 0.9
  phi2_min_ini = 0.1
  delta_phi1_ini = phi1_max_ini - phi1_min_ini
  delta_phi2_ini = phi2_max_ini - phi2_min_ini

  dk = phi1

  do iv = 0, nsize - 1
    part1 = 432.d0 * phi1(iv)
    part2 = 6.d0 * sqrt(750.d0 * (2.d0 * phi2(iv) - 1.d0)^3.d0 + 5184.d0 * phi1(iv)^2)
    ck = (part1 + part2)^(1.d0 / 3.d0) / 12.d0
    root = ck + 5.d0 / ck * (1.d0 / 24.d0 - 1.d0 / 12.d0 * phi2(iv))
    input = min((/max((/phi1_min_ini, root/)), phi1_max_ini/))
    finput = 0.9d0 - 0.8d0 * input^2

    part3 = ((phi1(iv) - input) / delta_phi1_ini)^2
    part4 = ((phi2(iv) - finput) / delta_phi2_ini)^2
    dk(iv) = sqrt(part3 + part4)
  end do
  return(dk)
end 


begin
  do_calc = True
  ;
  ; calculation of diagnostics needs all layers, plot needs small samping data to avoid NCL to be killd
  ;
  nv = 40962
  dt = 900
  nlev=2000
  glevel = "G6"
    filepath="./"+glevel+"L"+nlev+"_dt"+dt+"_E/"
    filein = addfile(filepath+"GRIST.ATM."+glevel+".tracer.00000-00000.1d.h1.nc", "r")
    area = new(nv*101, "double")
    do ilev = 0, 100
      area(nv*ilev: nv*(ilev+1)-1) = filein->areaPC 
    end do
    delete(filein)

    start_nlev = 1099 ; 33-L60, 999-L2000
    end_nlev = 1199 ; 37-L60, 1333-L2000
    filein = addfile(filepath+"GRIST.ATM."+glevel+".tracer.00000-00000.3d.h1.nc", "r")
    tmp1 = filein->tracerMxrt(:,start_nlev:end_nlev,0) ; 70
    tmp2 = filein->tracerMxrt(:,start_nlev:end_nlev,1) ; 70
    tmp1_r = tmp1(nlev|:,location_nv|:)
    tmp2_r = tmp2(nlev|:,location_nv|:)
    var1_ini = ndtooned(tmp1_r)
    var2_ini = ndtooned(tmp2_r)
    delete([/tmp1,tmp2,tmp1_r,tmp2_r,filein/])
    
    filein = addfile(filepath+"GRIST.ATM."+glevel+".tracer.00006-00000.3d.h1.nc", "r")
    tmp1 = filein->tracerMxrt(:,start_nlev:end_nlev,0) ; 20
    tmp2 = filein->tracerMxrt(:,start_nlev:end_nlev,1) ; 20
    tmp1_r = tmp1(nlev|:,location_nv|:)
    tmp2_r = tmp2(nlev|:,location_nv|:)
    var1 = ndtooned(tmp1_r)
    var2 = ndtooned(tmp2_r)
    delete(filein)
  if (do_calc) then
    dk = compute_dk(var1, var2, nv)
    tmp = lruo_norm(dk, var1, var2, area, nv)
    lr = tmp[0]
    lu = tmp[1]
    lo = tmp[2]
    print("lr="+lr)
    print("lu="+lu)
    print("lo="+lo)
    fbindirwrite("lr.dat", lr)
    fbindirwrite("lu.dat", lu)
    fbindirwrite("lo.dat", lo)
    exit
  else
    lr = fbindirread("lr.dat", 0, -1, "double")
    lu = fbindirread("lu.dat", 0, -1, "double")
    lo = fbindirread("lo.dat", 0, -1, "double")
  end if

  lr_value = sprintf("%8.3e", lr)
  lu_value = sprintf("%8.3e", lu)
  lo_value = sprintf("%8.3e", lo)
  
;;;;;;;;;;;;;;;;;;;;;; plot

  wks                   = gsn_open_wks("ps", "dcmip11.corr_tracer_"+glevel+"L"+nlev+"dt"+dt+"-E")
  gsn_define_colormap(wks, "WhiteBlueGreenYellowRed")
  res                   = True
  res@gsnFrame          = False
  res@gsnDraw           = False
  res@trYMaxF           = 1.0
  res@trYMinF           = 0.0
  res@trXMaxF           = 1.1
  res@trXMinF           = -0.1
  
  res@tmXBMode          = "Explicit"
  res@tmXBValues        = (/0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.1/)
  res@tmXBLabels        = (/0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.1/)
  res@tmXBLabelFontHeightF = 0.02

  res@tmYLMode          = "Explicit"
  res@tmYLValues        = (/0., 0.2, 0.4, 0.6, 0.8, 1.0/)
  res@tmYLLabels        = (/0., 0.2, 0.4, 0.6, 0.8, 1.0/)
  res@tmYLLabelFontHeightF = 0.02

  res@tiXAxisFontHeightF   = 0.02
  res@tiYAxisFontHeightF   = 0.02

  res@xyMonoDashPattern = False
  res@xyMarkLineMode    = "markers"

  plot0 = wrap_plot_map(wks, res, var1, var2, var1_ini, var2_ini, lr_value, lu_value, lo_value)

  draw(plot0)

  frame(wks)
end
