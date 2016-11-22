load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;; Read data from raw files
file1 = addfile("pr_Amon_MPI-ESM-P_lgm_r1i1p1_185001-194912.nc","r")
;print(file1)


lat=file1->lat
lon=file1->lon
time=file1->time
pr=file1->pr

;; Check for the dimensions
print(dimsizes(time))
print(dimsizes(lat))
print(dimsizes(lon))
print(dimsizes(pr))

;; Create new varible to save monthly average data
;; new_monthly(month,lat,lon)
new_monthly= new((/12,dimsizes(lat),dimsizes(lon)/),float,"No_FillValue")
i=0
do while (i .lt. 12)
monthly=pr(i:dimsizes(time)-1:12,:,:)
avg_monthly=dim_avg_n(monthly,0)*30*24*60*60/0.9
new_monthly(i,:,:)=avg_monthly(:,:)
i=i+1
end do


new_monthly!0="time"
new_monthly!1="lat"
new_monthly!2="lon"

;; Create new variables for each dimensions of new_monthly
new_monthly&time = new((/12/),float,"No_FillValue")
new_monthly&lat  = new((/dimsizes(lat)/),double,"No_FillValue")
new_monthly&lon  = new((/dimsizes(lon)/),double,"No_FillValue")


;; Test the assignment
;; To simplify the question, we just need to use the date with 
;; lattitude above 30N
;test_data	 = new((/12/),float,"No_FillValue")
;print(file1->lat(85:127))


;; Assign the first dimension time
;; 1-12
new_monthly&time=ispan(1,12,1)
;print(new_monthly&time)

;; Assign the second and third dimensions
new_monthly&lat=file1->lat  
new_monthly&lon=file1->lon

;; Check the assignment
;print(new_monthly&lat)

;; Name the variable and unit
new_monthly@long_name="monthly precipitation"
new_monthly@units="mm"


;; Create a new .nc to save new_monthly file
system("rm -f monthlydata_pr.nc")
out=addfile("monthlydata_pr.nc","c")


;assign the value to out
out->pr = new_monthly
out->time = ispan(1,12,1)
out->lat  = lat
out->lon  = lon



tpr = dim_sum_n(new_monthly,0)  

;print(dimsizes(tpr))

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
re_tpr  = area_conserve_remap_Wrap (lon,lat, tpr , LON, LAT, opt)  ; (ntim,128,256)

;; Write ASCII Files for regridded data
system("rm -f tpr.txt")
asciiwrite("tpr.txt",re_tpr)

tpr= new((/dimsizes(LAT),dimsizes(LON)/),float,"No_FillValue")
tpr!0="lat"
tpr!1="lon"

tpr&lat  = new((/dimsizes(LAT)/),double,"No_FillValue")
tpr&lon  = new((/dimsizes(LON)/),double,"No_FillValue")

tpr = re_tpr

tpr&lat= LAT  
tpr&lon= LON




;create a new .nc to save new_monthly file
system("rm -f tpr.nc")
out=addfile("tpr.nc","c")

;; Assign the value to out
out->pr = tpr
;out->lat  = LAT
;out->lon  = LON











end
