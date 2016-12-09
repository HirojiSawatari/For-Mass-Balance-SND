;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;; Read data from raw files
file1=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-01_ave.nc","r")
file2=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-02_ave.nc","r")
file3=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-03_ave.nc","r")
file4=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-04_ave.nc","r")
file5=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-05_ave.nc","r")
file6=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-06_ave.nc","r")
file7=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-07_ave.nc","r")
file8=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-08_ave.nc","r")
file9=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-09_ave.nc","r")
file10=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-10_ave.nc","r")
file11=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-11_ave.nc","r")
file12=addfile("/home/sharing/Bering_strait/mgoh/b30.232mgoh.cam2.h0.2100-2899-12_ave.nc","r")

lat=file1->lat
lon=file1->lon
time=file1->time
precc=file1->PRECC
precl=file1->PRECL

lat2=file2->lat
lon2=file2->lon
time2=file2->time
precc2=file2->PRECC
precl2=file2->PRECL

lat3=file3->lat
lon3=file3->lon
time3=file3->time
precc3=file3->PRECC
precl3=file3->PRECL

lat4=file4->lat
lon4=file4->lon
time4=file4->time
precc4=file4->PRECC
precl4=file4->PRECL

lat5=file5->lat
lon5=file5->lon
time5=file5->time
precc5=file5->PRECC
precl5=file5->PRECL

lat6=file6->lat
lon6=file6->lon
time6=file6->time
precc6=file6->PRECC
precl6=file6->PRECL

lat7=file7->lat
lon7=file7->lon
time7=file7->time
precc7=file7->PRECC
precl7=file7->PRECL

lat8=file8->lat
lon8=file8->lon
time8=file8->time
precc8=file8->PRECC
precl8=file8->PRECL

lat9=file9->lat
lon9=file9->lon
time9=file9->time
precc9=file9->PRECC
precl9=file9->PRECL

lat10=file10->lat
lon10=file10->lon
time10=file10->time
precc10=file10->PRECC
precl10=file10->PRECL

lat11=file11->lat
lon11=file11->lon
time11=file11->time
precc11=file11->PRECC
precl11=file11->PRECL

lat12=file12->lat
lon12=file12->lon
time12=file12->time
precc12=file12->PRECC
precl12=file12->PRECL
;print(file1)

;; Check for the dimensions
print(dimsizes(time))
print(dimsizes(lat))
print(dimsizes(lon))
print(dimsizes(precl))

;; Create new varible to save monthly average data
;; new_monthly(month,lat,lon)
new_monthly= new((/12,dimsizes(lat),dimsizes(lon)/),float,"No_FillValue")

new_monthly(0,:,:)=(precc(0,:,:)+precl(0,:,:))*1000*31*24*60*60
new_monthly(1,:,:)=(precc2(0,:,:)+precl2(0,:,:))*1000*28*24*60*60
new_monthly(2,:,:)=(precc3(0,:,:)+precl3(0,:,:))*1000*31*24*60*60
new_monthly(3,:,:)=(precc4(0,:,:)+precl4(0,:,:))*1000*30*24*60*60
new_monthly(4,:,:)=(precc5(0,:,:)+precl5(0,:,:))*1000*31*24*60*60
new_monthly(5,:,:)=(precc6(0,:,:)+precl6(0,:,:))*1000*30*24*60*60
new_monthly(6,:,:)=(precc7(0,:,:)+precl7(0,:,:))*1000*31*24*60*60
new_monthly(7,:,:)=(precc8(0,:,:)+precl8(0,:,:))*1000*31*24*60*60
new_monthly(8,:,:)=(precc9(0,:,:)+precl9(0,:,:))*1000*30*24*60*60
new_monthly(9,:,:)=(precc10(0,:,:)+precl10(0,:,:))*1000*31*24*60*60
new_monthly(10,:,:)=(precc11(0,:,:)+precl11(0,:,:))*1000*30*24*60*60
new_monthly(11,:,:)=(precc12(0,:,:)+precl12(0,:,:))*1000*31*24*60*60

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
