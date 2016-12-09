;load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;; Read data from raw files
file1=addfile("/home/sharing/Bering_strait/b30.232mgch2.clm2.h0.3998.nc","r")

;; Test for Reading
;print(file1)

lat=file1->lat
lon=file1->lon
time=file1->time
soil=file1->SOILICE

;; Check for the dimensions
print(dimsizes(time))
print(dimsizes(lat))
print(dimsizes(lon))
print(dimsizes(soil))

new_monthly= new((/12,dimsizes(lat),dimsizes(lon)/),float,"No_FillValue")

new_monthly(0,:,:)=soil(0,0,:,:)
new_monthly(1,:,:)=soil(1,0,:,:)
new_monthly(2,:,:)=soil(2,0,:,:)
new_monthly(3,:,:)=soil(3,0,:,:)
new_monthly(4,:,:)=soil(4,0,:,:)
new_monthly(5,:,:)=soil(5,0,:,:)
new_monthly(6,:,:)=soil(6,0,:,:)
new_monthly(7,:,:)=soil(7,0,:,:)
new_monthly(8,:,:)=soil(8,0,:,:)
new_monthly(9,:,:)=soil(9,0,:,:)
new_monthly(10,:,:)=soil(10,0,:,:)
new_monthly(11,:,:)=soil(11,0,:,:)

sftgif=new((/dimsizes(lat),dimsizes(lon)/),float,"No_FillValue")

do i=0,63
	do j=0,127
		min_sft = min(new_monthly(0:11,i,j))
		if(ismissing(min_sft)) then
			
		else 
			sftgif(i,j) = min(new_monthly(0:11,i,j))
		end if
	end do
end do

new_monthly!0="time"
new_monthly!1="lat"
new_monthly!2="lon"

;; Create new variables for each dimensions of new_monthly
new_monthly&time = new((/12/),float,"No_FillValue")
new_monthly&lat  = new((/dimsizes(lat)/),float,"No_FillValue")
new_monthly&lon  = new((/dimsizes(lon)/),float,"No_FillValue")

;; Assign the first dimension time
;; 1-12
new_monthly&time=ispan(1,12,1)
;print(new_monthly&time)

;; Assign the second and third dimensions
new_monthly&lat=file1->lat  
new_monthly&lon=file1->lon

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
re_sftgif  = area_conserve_remap_Wrap (lon,lat, sftgif , LON, LAT, opt)  ; (ntim,128,256)

do i=0,127
	do j=0,255
		if(re_sftgif(i,j) .gt. 10.)
			re_sftgif(i,j)=1
		else
        		re_sftgif(i,j)=0
		end if
	end do
end do

;; Write ASCII Files for regridded data
system("rm -f re_sft.txt")
asciiwrite("re_sft.txt",re_sftgif)

;sftgif= new((/dimsizes(LAT),dimsizes(LON)/),integer,"No_FillValue")
;sftgif!0="lat"
;sftgif!1="lon"

;sftgif&lat  = new((/dimsizes(LAT)/),double,"No_FillValue")
;sftgif&lon  = new((/dimsizes(LON)/),double,"No_FillValue")

;sftgif= re_sftgif

;sftgif&lat= LAT  
;sftgif&lon= LON

;create a new .nc to save new_monthly file
system("rm -f re_sft.nc")
out=addfile("re_sft.nc","c")

;; Assign the value to out
out->sftgif = re_sftgif

end

