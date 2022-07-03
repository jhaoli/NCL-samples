tmp = new(/5,4/), float)

;assign dimension names
tmp!0 = "lon"
tmp!1 = "lat"

;initialize the values in the coordinate arrays
lon_coor = (/-90, 85, -80, -75, -70/)
lat_coor = (/0, 20, 40, 60/)

;assign the coordinates arrays to the variable
tmp&lon = lon_coor
tmp&lat = lat_coor

