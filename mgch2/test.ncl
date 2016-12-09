;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;; Read data from raw files
file1=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-01_ave.nc","r")
file2=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-02_ave.nc","r")
file3=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-03_ave.nc","r")
file4=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-04_ave.nc","r")
file5=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-05_ave.nc","r")
file6=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-06_ave.nc","r")
file7=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-07_ave.nc","r")
file8=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-08_ave.nc","r")
file9=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-09_ave.nc","r")
file10=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-10_ave.nc","r")
file11=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-11_ave.nc","r")
file12=addfile("/home/sharing/Bering_strait/mgch2/b30.232mgch2.cam2.h0.2600-2999-12_ave.nc","r")

;; Test for Reading
;print(file1)

lat=file1->lat
lon=file1->lon
time=file1->time
tas=file1->TS

lat2=file2->lat
lon2=file2->lon
time2=file2->time
tas2=file2->TS

lat3=file3->lat
lon3=file3->lon
time3=file3->time
tas3=file3->TS

lat4=file4->lat
lon4=file4->lon
time4=file4->time
tas4=file4->TS

lat5=file5->lat
lon5=file5->lon
time5=file5->time
tas5=file5->TS

lat6=file6->lat
lon6=file6->lon
time6=file6->time
tas6=file6->TS

lat7=file7->lat
lon7=file7->lon
time7=file7->time
tas7=file7->TS

lat8=file8->lat
lon8=file8->lon
time8=file8->time
tas8=file8->TS

lat9=file9->lat
lon9=file9->lon
time9=file9->time
tas9=file9->TS

lat10=file10->lat
lon10=file10->lon
time10=file10->time
tas10=file10->TS

lat11=file11->lat
lon11=file11->lon
time11=file11->time
tas11=file11->TS

lat12=file12->lat
lon12=file12->lon
time12=file12->time
tas12=file12->TS

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
;print(lat)
;print(lon)
;print(time)
;print (tas)

;; Create new varible to save monthly average data
;; new_monthly(month,lat,lon)
new_monthly= new((/12,dimsizes(lat),dimsizes(lon)/),float,"No_FillValue")

new_monthly(0,:,:)=tas(0,:,:)-273.15
new_monthly(1,:,:)=tas2(0,:,:)-273.15
new_monthly(2,:,:)=tas3(0,:,:)-273.15
new_monthly(3,:,:)=tas4(0,:,:)-273.15
new_monthly(4,:,:)=tas5(0,:,:)-273.15
new_monthly(5,:,:)=tas6(0,:,:)-273.15
new_monthly(6,:,:)=tas7(0,:,:)-273.15
new_monthly(7,:,:)=tas8(0,:,:)-273.15
new_monthly(8,:,:)=tas9(0,:,:)-273.15
new_monthly(9,:,:)=tas10(0,:,:)-273.15
new_monthly(10,:,:)=tas11(0,:,:)-273.15
new_monthly(11,:,:)=tas12(0,:,:)-273.15

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

