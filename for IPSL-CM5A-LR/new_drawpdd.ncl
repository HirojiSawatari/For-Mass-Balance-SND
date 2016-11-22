load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

  d0	= addfile("re_pdd.nc","r")
  lat=d0->lat
  lon=d0->lon
  data0	= d0->data

  file2 = addfile("landsea_mask.nc","r")
  land_mask = file2->mask

  LON  = fspan(0,360, 144)
  LON@units = "degrees_east"
  LON!0="lon"
  NLATG= 96           ; appropriate global value
  LAT = fspan(-90,90,96)
  LAT@units = "degrees_north"
  LAT!0="lat"

  opt = True
  opt@NLATo = NLATG
  re_data0 = area_conserve_remap_Wrap (lon,lat,data0 , LON, LAT, opt)

  re_data0 = mask(re_data0, land_mask .eq. 0, False)

  wks = gsn_open_wks("png","newpdd")
  plot = new(1,graphic)
  ;gsn_define_colormap(wks,"BlWhRe")

  tsres = True
  ;tsres@gsnZonalMean               = True
  tsres@cnFillOn = True
  tsres@cnMinLevelValF        =  0            ; set min contour level
  tsres@cnMaxLevelValF        =  10000            ; set max contour level
  tsres@cnLevelSpacingF       =  500
  tsres@gsnSpreadColors = True
  tsres@gsnSpreadColorStart = 24
  tsres@gsnSpreadColorEnd = -26
  tsres@tiMainString         = "IPSL-CM5A-LR LGM PDD"    ; plot title

  tsres@gsnAddCyclic = False

  tsres@mpCenterLonF         = 180.

  tsres@cnLinesOn = False
  tsres@cnLevelSelectionMode = "ManualLevels"
  plot = gsn_csm_contour_map_ce(wks,re_data0,tsres)

end
