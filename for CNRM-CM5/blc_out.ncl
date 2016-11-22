load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"

begin

z0 = asciiread("balance.txt",(/128,256/),"float")
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

PDD@long_name = "Mass Balance"
PDD@units = "mm"


system("rm -f balance.nc")
out=addfile("balance.nc","c")

out->data = PDD

open = addfile("balance.nc","r")

print(open)
end
