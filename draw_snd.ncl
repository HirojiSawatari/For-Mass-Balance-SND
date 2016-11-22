load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;; Read data from raw files
file1 = addfile("snd_LImon_MIROC-ESM_lgm_r1i1p1_460001-469912.nc","w")
;print(file1)


lat=file1->lat
lon=file1->lon
time=file1->time
pr=file1->snd


do i=0,63
do j=0,127
pr(0,i,j) = dim_avg_n(pr(:,i,j),0)
end do
end do


wks = gsn_open_wks("png","snowdepth")
plot = new(1,graphic)
;gsn_define_colormap(wks,"BlWhRe")

tsres = True
;tsres@gsnZonalMean               = True
tsres@cnFillOn = True
tsres@cnMinLevelValF        =  0            ; set min contour level
tsres@cnMaxLevelValF        =  3.25            ; set max contour level
tsres@cnLevelSpacingF       =  0.25
tsres@gsnSpreadColors = True
tsres@gsnSpreadColorStart = 24
tsres@gsnSpreadColorEnd = -26
tsres@tiMainString         = "MIROC-ESM Snow Depth"    ; plot title

tsres@gsnAddCyclic = False

tsres@mpCenterLonF         = 180.

tsres@cnLinesOn = False
tsres@cnLevelSelectionMode = "ManualLevels"
plot = gsn_csm_contour_map_ce(wks,pr(0,:,:),tsres)

end
