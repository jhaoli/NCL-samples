  procedure calc_energy_spectrum(rad, ugrd, vgrd, divspec, vrtspec, verbose)
  local num_lon, num_lat, lon, lat, num_glat, glat, \
      ugau, vgau, div, vrt, im, in
  begin

    dims    = dimsizes(ugrd)
    num_lat = dims(0)
    num_lon = dims(1)
    trunc_num = (num_lon - 1) / 3
    lon = ugrd&lon
    lat = vgrd&lat
    if (lat(0) .gt. lat(1)) then
      system("echo [Error]: input values must be in ascending latitude order")
      exit()
    end if 
    ;------------------------------------------------------------
    if (verbose) then 
      system("echo [Notice]: Interpolate data onto Gaussian grid.")
    end if
    num_glat   = num_lon / 2
    gau_info   = tofloat(gaus(num_glat / 2)) ; divide by 2 to get "per hemisphere"
    glat       = gau_info(:, 0) ; gaussain latitudes (1st dimension of gau_info) 
    ; gwgt      = gau_info(:, 1) ; gaussain weights (2nd dimension of gau_info)
    ugau       = new((/num_glat, num_lon/), typeof(ugrd))
    vgau       = new((/num_glat, num_lon/), typeof(vgrd))
    ugau       = linint2(lon, lat, ugrd, True, lon, glat, 0)
    vgau       = linint2(lon, lat, vgrd, True, lon, glat, 0)
    ugau!0     = "lat"
    ugau!1     = "lon"
    ugau&lat   = glat
    ugau&lon   = lon 
    copy_VarCoords(ugau, vgau)
    ;------------------------------------------------------------
    if (verbose) then
      system("echo [Notice]: Calculate divergence.")
    end if
    div        = uv2dvG_Wrap(ugau, vgau)
    ;------------------------------------------------------------
    if (verbose) then
      system("echo [Notice]: Calculate vorticity.")
    end if
    vrt        = uv2vrG_Wrap(ugau, vgau)
    ;------------------------------------------------------------
    if (verbose) then
      system("echo [Notice]: Calculate spherical harmonic of divergent components.")
    end if
    coefd      = shagC(div) ; the dimensions coefd will be 2 x nlat x N, where coefd(0,nlat,N) contains the "real" coefficients \
                       ; and coefd(1,nlat,N) contains the "imaginary" coefficients. 
                       ; N = minimum(nlat, (nlon+2)/2) if nlon is even
                       ; N = minimum(nlat, (nlon+1)/2) if nlon is odd
    ; coefd = tri_trunC(coefd, trunc_num)
    cr = coefd(0, :, :) ; real coef (nlat, trunc_num)
    ci = coefd(1, :, :) ; imaginary coef (nlat, trunc_num)
    pwr = cr^2 + ci^2   ; (nlat, nlat)
    do n = 1, trunc_num
      divspec(n) = (pwr(n,0) + 2. * sum(pwr(n,1:n))) * rad^2 / 4 / n / (n + 1)
    end do
    ;------------------------------------------------------------
    if (verbose) then
      system("echo [Notice]: Calculate spherical harmonic of rotational components.")
    end if                  
    coefr      = shagC(vrt)
    ; coefr     = tri_trunC(coefr, trunc_num)
    cr = coefr(0, :, :) ; real coef (nlat, trunc_num)
    ci = coefr(1, :, :) ; imaginary coef (nlat, trunc_num)
    pwr = cr^2 + ci^2   ; (nlat, nlat)
    
    do n = 1, trunc_num
      vrtspec(n) = (pwr(n,0) + 2. * sum(pwr(n, 1:n))) * rad^2 / 4. / n / (n + 1)
    end do
    ;------------------------------------------------------------
    ; if (verbose) then
    ;   system("echo [Notice]: Calculate kinetic energy spectrum.")
    ; end if
    ; ; 0-divergenct part; 1-rotational part;
    ; spec(0,0) = 0.0
    ; spec(1,0) = 0.0
    ; do in = 1, trunc_num
    ;   spec(0,in) = coefd(0,in,0)^2 + coefd(1,in,0)^2
    ;   spec(1,in) = coefr(0,in,0)^2 + coefr(1,in,0)^2
    ;   do im = 1, in
    ;     spec(0,in) = spec(0,in) + 2 * (coefd(0,in,im)^2 + coefd(1,in,im)^2)
    ;     spec(1,in) = spec(1,in) + 2 * (coefr(0,in,im)^2 + coefr(1,in,im)^2)
    ;   end do
    ;   spec(0,in) = spec(0,in) * rad^2 / 4 / in / (in+1)
    ;   spec(1,in) = spec(1,in) * rad^2 / 4 / in / (in+1)
    ; end do
    
    return (/divspec, vrtspec/)
  end 

  procedure save_energy_spectrum_wavelength(divspec, vrtspec, wavelength, outfile)
  begin
    str = str_join((/"/bin/rm -f ", outfile/), "")
    system(str)
    ncdf = addfile(outfile, "c")
    fAtt = True
    fAtt@title = "KE spectrum and wavelength"
    fAtt@creation_date = systemfunc("date")
    fileattdef(ncdf, fAtt)
    trunc_num  = dimsizes(divspec) - 1
    dimNames = (/"trunc_num"/)
    dimSizes = (/trunc_num/)
    dimUnlim = (/False/)
    filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

    filevardef(ncdf, "divspec"    , typeof(divspec)      , (/"trunc_num"/))
    filevardef(ncdf, "vrtspec"    , typeof(vrtspec)      , (/"trunc_num"/))
    filevardef(ncdf, "wavelength" , typeof(wavelength)   , (/"trunc_num" /))
    setfileoption(ncdf,"DefineMode",False)

    ncdf->divspec    = divspec(1:trunc_num)
    ncdf->vrtspec    = vrtspec(1:trunc_num)
    ncdf->wavelength = wavelength

  end

  procedure plot_energy_spectrum(fig_name, fin, resolution)
  begin

    wks        = gsn_open_wks("ps", fig_name)
    fi         = addfile(fin, "r")
    divspec    = fi->divspec
    vrtspec    = fi->vrtspec
    wavelength = fi->wavelength
    trunc_num  = dimsizes(divspec)
    res1               = True
    res1@gsnFrame      = False
    res1@trXLog        = True
    res1@trYLog        = True
    res1@trXMinF       = 1.0
    if (resolution .eq. 0.05) then
      res1@trXMaxF     = trunc_num + 1.e3
    else 
      res1@trXMaxF     = trunc_num + 1.e2  ; 1degree 1.e2, 0.05degree 1.e3
    end if 
    if (resolution .eq. 0.05) then
      res1@trYMinF       = 1.e-6
    else
      res1@trYMinF       = 1.e-8
    end if
    res1@trYMaxF        = 1.e3
    res1@xyLineThicknessF = 2.0
    res1@tiXAxisString  = "Wavenumber"
    res1@tiYAxisString  = "Kinetic energy (m~S~2~N~ s~S~-2~N~)"
   
    res2                = True
    res2@gsnFrame       = False
    res2@trXLog         = True
    res2@trYLog         = True
    res2@trXMinF        =  2 * 3.14159265 * 6371e3 * 1.0e-3 / res1@trXMaxF
    res2@trXMaxF        = max(wavelength)
    res2@trYMinF        = res1@trYMinF
    res2@trYMaxF        = res1@trYMaxF 
    res2@trXReverse     = True 
    res2@tmXTMode       = "Explicit"
    res2@tmXUseBottom   = False
    res2@tmXTLabelsOn   = True
    ;------------------------------------------------------------
    ; reference lines (power -3 and power -5/3)
    x1 = 10
    x2 = 4000
    y1 = 1.0e2
    y2 = exp(-3.0*log(x2/x1))*y1
    res1@tmXTOn       = False
    res1@xyLineColors = (/"grey"/)
    res1@xyLineThicknessF = 4.0
    res1@tiMainString   = ""
    plot = gsn_csm_xy(wks, (/x1,x2/), (/y1,y2/), res1)
    res1@tiMainString   = ""
    res1@xyDashPatterns = 1
    y2 = exp(-5.0/3.0*log(x2/x1))*y1
    plot = gsn_csm_xy(wks, (/x1,x2/), (/y1,y2/), res1)
    ;------------------------------------------------------------
    ; kinetic energy spetrum
    res1@xyLineThicknessF = 3.0
    delete([/res1@xyLineColors, res1@xyDashPatterns/])
    res1@xyLineColor    = "black"
    res2@xyLineColor    = "black"
    totspec = divspec + vrtspec
    plot = gsn_csm_x2y(wks, ispan(1, trunc_num, 1), wavelength, totspec, res1, res2) ; total
   
    res1@xyLineColor = "blue"
    res2@xyLineColor = "blue"
    plot = gsn_csm_x2y(wks, ispan(1, trunc_num, 1), wavelength, vrtspec, res1, res2) ; rotational

    res1@xyLineColor = "red"
    res2@xyLineColor = "red"
    res2@tiXAxisString = "Wavelength (km)"
    plot = gsn_csm_x2y(wks, ispan(1, trunc_num, 1), wavelength, divspec, res1, res2) ; divergent
    delete(res2@tiXAxisString)
    ;------------------------------------------------------------
    ; effective resolution
    res2@tmXBOn       = False
    res2@tmXTOn       = False
    res2@xyLineColor  = "cyan"
    res2@xyLineThicknessF = 3.0
    data = new((/3,10/), typeof(wavelength))
    data(0,:) = 4 * resolution * 111
    data(1,:) = 6 * resolution * 111
    data(2,:) = fspan(res1@trYMinF, res1@trYMaxF, 10)
    res2@xyDashPatterns = 0
    plot1 = gsn_csm_xy(wks, data(0,:), data(2,:), res2)
    res2@xyDashPatterns = 1
    plot1 = gsn_csm_xy(wks, data(1,:), data(2,:), res2)
    ;----------------------------------------------------------------------
    ; add legend
    lgres                     = True
    lgres@lgLabelFontHeightF  = 0.1
    lgres@vpWidthF            = 0.25
    lgres@vpHeightF           = 0.2
    lgres@lgPerimOn           = True
    lgres@lgItemOrder         = (/6,5,4,3,2,1,0/) 
    lgres@lgDashIndexes       = (/0,0,0,0,1,0,1/)
    lgres@lgLineColors        = (/"black", "blue", "red", "grey", "grey", "cyan", "cyan"/)
    lgres@lgLineThicknessF    = 4.0
    labels                    = (/" Total", " Rotational", " Divergent", " -3 power law", " -5/3 power law", " 4~F33~D~F~", " 6~F33~D~F~"/)
    legend                    = gsn_create_legend(wks, 7, labels, lgres)
    amres                     = True
    amres@amZone              = 0
    amres@amParallelPosF      = -0.25 ; left-right
    amres@amOrthogonalPosF    = 0.3 ; up-down
    annota = gsn_add_annotation(plot, legend, amres)
    draw(plot)
    frame(wks)
  end

begin
  pi = 3.14159265
  rad = 6371e3

  model = "MCV"
  if (model .eq. "MCV") then
    fn = "../MCV/20230626/uv200-0116-divampingon.nc"
    f = addfile(fn, "r")
    resolution = 0.5
    u = f->u(0,0,:,:)
    v = f->v(0,0,:,:)
  else if (model .eq. "ERA5") then
    fn = "../ERA5/era5.pl.2022101512.0p25_uv200hPa.nc"
    f = addfile(fn, "r")
    resolution = 0.25
    u = short2flt(f->u(0,:,:))
    v = short2flt(f->v(0,:,:))
  else if (model .eq. "FNL") then
    fn = "../FNL/fnl_uv_200hPa.nc"
    f = addfile(fn, "r")
    resolution = 1
    u = f->u(0,:,:)
    v = f->v(0,:,:)
  else if (model .eq. "MERRA2") then
    fn = "../MERRA2/MERRA2_300.instU_3d_ana_Np.201012.nc4"
    f = addfile(fn, "r")
    resolution = 0.5
    u = f->U(2,22,:,:)
    v = f->V(2,22,:,:)
  else if (model .eq. "GFS") then
    fn = "../GFS/gfs.0p25.2015011500.f024.grib2_uv.nc"
    f = addfile(fn, "r")
    resolution = 0.25
    u = f->u(0,:,:)
    v = f->v(0,:,:)
  else if (model .eq. "IFS") then
    fn = "../ECMWF/uv200hPa.nc"
    f = addfile(fn, "r")
    resolution = 0.2
    u = short2flt(f->u(0,:,:))
    v = short2flt(f->v(0,:,:))
  end if 
  end if 
  end if
  end if
  end if
  end if

  num_lon = getfilevardimsizes(f, "lon")
  trunc_num = (num_lon - 1) / 3

  wavelength  = new((/trunc_num/), typeof(u))
  wavelength = 2 * pi * rad * 1.0e-3 / ispan(1, trunc_num, 1)

  vrtspec = new((/trunc_num+1/), typeof(u))

  calc_energy_spectrum(rad, u, v, divspec, vrtspec, True)

  f_save = "ke_mcv_20220720_0p5_divdamp4.nc"
  save_energy_spectrum_wavelength(divspec, vrtspec, wavelength, f_save)

  system("echo [Notice]: Plot figure.")
  fo_name = str_join((/"energy_spectrum", tostring(model), tostring(resolution), "uv200-0116-tangj"/), "_")
  plot_energy_spectrum(fo_name, f_save, resolution)

end
