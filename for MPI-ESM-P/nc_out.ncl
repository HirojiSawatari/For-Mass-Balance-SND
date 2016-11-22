load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"

begin

z0 = asciiread("re_output.txt",(/128,256/),"float")
;z0 = 30 *z1
;print(z0)
z0!0 = "lat"
z0!1 = "lon"

z = z0(lat|:,lon|:)

;;create new variables for each dimensions of new_monthly
PDD = new((/128,256/),float,"No_FillValue")
PDD = z
PDD&lat  = new((/128/),double,"No_FillValue")
PDD&lon  = new((/256/),double,"No_FillValue")

PDD!0 = "lat"
PDD!1 = "lon"
PDD&lat = fspan(-90,90,128)
PDD&lat@units = "degrees_north"
PDD&lat@long_name = "lattitude"

PDD&lon = fspan(0,360,256)
PDD&lon@units = "degrees_east"
PDD&lon@long_name = "lontitude"
print(PDD&lon)

PDD@long_name = "Positive Degree Day"
PDD@units = "units"


system("rm -f re_pdd.nc")
out=addfile("re_pdd.nc","c")

out->data = PDD

open = addfile("re_pdd.nc","r")

print(open)
end

