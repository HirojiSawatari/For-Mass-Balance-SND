load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;; Read data from raw files
file1=addfile("tas_Amon_GISS-E2-R_lgm_r1i1p150_300001-302412.nc","r")

;; Test for Reading
;print(file1)

lat=file1->lat
lon=file1->lon
time=file1->time
tas=file1->tas

;; Check for the dimensions
print(dimsizes(time))
print(dimsizes(lat))
print(dimsizes(lon))
print(dimsizes(tas))


;; No need
;; Tranfer format from double to float
;lat = doubletofloat(lat)
;lon = doubletofloat(lon)

;; Check for Assignment
print(lat)
;print(lon)
;print(time)
;print (tas)

;; Create new varible to save monthly average data
;; new_monthly(month,lat,lon)
new_monthly= new((/12,dimsizes(lat),dimsizes(lon)/),float,"No_FillValue")
i=0
do while (i .lt. 12)
	monthly=tas(i:dimsizes(time)-1:12,:,:)
	avg_monthly=dim_avg_n(monthly,0)
	new_monthly(i,:,:)=avg_monthly(:,:)-273.15
i=i+1
end do

;; Check for the assignment of monthly average data
;print(new_monthly)

;; Assign the dimensions
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
new_monthly@long_name="monthly average temperature"
new_monthly@units="C"


;; Create a new .nc to save new_monthly file
system("rm -f monthlydata_CNRM.nc")
out=addfile("monthlydata_CNRM.nc","c")

;; Assign the value to out
out->data_monthly = new_monthly
out->time = ispan(1,12,1)
out->lat  = lat
out->lon  = lon

;; Check if the output is the same as expected
test_data1 = new_monthly(:,2,2)
print(test_data1)
;test_data2 = new_monthly(:,1,1)
;print(test_data2)

;; Write AscII files
system("rm -f monthlydata_CNRM.txt")
asciiwrite("monthlydata_CNRM.txt",out->data_monthly)



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
po  = area_conserve_remap_Wrap (lon,lat, new_monthly , LON, LAT, opt)  ; (ntim,128,256)

;; Check for the output
;print(po)
;print(LAT)
;print(LON)

test_data3 = po(:,2,2)
;print(test_data3)
test_data4 = po(:,1,1)
;print(test_data4)

;; Write ASCII Files for regridded data
system("rm -f re_monthlydata.txt")
asciiwrite("re_monthlydata.txt",po)

;.........................................................
;; The following is dispense for the author has use C++
;; to accelerate the programming


;TMJ = dim_avg_n(new_tas,0)
;TMA = dim_avg_n(tas,0)

;A = 365
;pi = 3.1415926

;TCA = TMA +(TMJ-TMA)*cos(2*pi*5/A)
;print(TCA)

;asciiwrite("/lustre/delver/WEI/NCL/TMJ_CNRM-CM5.txt",TMJ)
;asciiwrite("/lustre/delver/WEI/NCL/TMA_CNRM-CM5.txt",TMA)
;asciiwrite("/lustre/delver/WEI/NCL/LON_CNRM-CM5.txt",lon)
;asciiwrite("/lustre/delver/WEI/NCL/LAT_CNRM-CM5.txt",lat)


;asciiwrite("D:/Data_Output/TMJ_CNRM-CM5.txt",TMJ)
;asciiwrite("D:/Data_Output/TMA_CNRM-CM5.txt",TMA)

;print(TMJ)
;print(TMA)
;.........................................................

end

