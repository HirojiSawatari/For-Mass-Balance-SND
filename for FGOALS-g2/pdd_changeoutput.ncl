load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

file1 = addfile("re_pdd2.nc","r")

;print(file1)

lat=file1->lat
lon=file1->lon
pdd=file1->data



;; The following is to regridding to 12*128*256
LON  = fspan(0,360, 256)
LON@units = "degrees_east"
LON!0="lon"
NLATG= 128           ; appropriate global value
;LATG = latGau (NLATG, "LATG", "latitude", "degrees_north")      ; nominally 90S -> 90N
;LAT  = LATG({-90:90})          ; extract gaussian latitudes for region
LAT = fspan(-90,90,128)
LAT@units = "degrees_north"
LAT!0="lat"

opt = True
opt@NLATo = NLATG              
re_pddd  = area_conserve_remap_Wrap (lon,lat, pdd , LON, LAT, opt)  ; (ntim,128,256)



;; Write ASCII Files for regridded data
system("rm -f re_pdd.txt")
asciiwrite("re_pdd.txt",re_pddd)


re_pdd= new((/dimsizes(LAT),dimsizes(LON)/),float,"No_FillValue")
re_pdd!0="lat"
re_pdd!1="lon"

re_pdd&lat  = new((/dimsizes(LAT)/),double,"No_FillValue")
re_pdd&lon  = new((/dimsizes(LON)/),double,"No_FillValue")

re_pdd = re_pddd

re_pdd&lat= LAT  
re_pdd&lon= LON


re_pdd@long_name="Surface type glacier"
re_pdd@units = "Glacier"

;create a new .nc to save new_monthly file
system("rm -f re_pdd.nc")
out=addfile("re_pdd.nc","c")

;; Assign the value to out
out->sftgif = re_pdd
;out->lat  = LAT
;out->lon  = LON


end
