load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

	d0	= addfile("./balance.nc","r")
	d1	= addfile("./re_sft.nc","r")
	d2  = addfile("./resft.nc","r")
	data0	= d0->data	
	;lon0    = d0->lon
	;lat0    = d0->lat
	;data0&lon = lon0
	;data0&lat = lat0
	data1	= d1->sftgif	
	;lon1    = d1->lon
	;lat1    = d1->lat
	;data1&lon = lon1
	;data1&lat = lat1
	lf       = d2->data(0,:,:)
	data0a = where( lf .gt. 10., data0, data0@_FillValue)
	copy_VarMeta(data0, data0a)

	;printVarSummary(data0a)
	;printVarSummary(data1)

	cnlvls = fspan(-2400,2400,25)
;-------plot the MPI contour---------------------------------
   	wks_type    = "png"
	wks_type@wkPaperWidthF  = 8  ; in inches
	wks_type@wkPaperHeightF = 6  ; in inches
	xwks = gsn_open_wks(wks_type,"overlay")      ; Open an workstation.

	res	=True
        res@tiMainString         = "IPSL-CM5A-LR Mass Balance"   
	res@gsnDraw 	= False           	; don't draw
        res@gsnFrame 	= False            	; don't advance frame
        res@cnFillOn	= True
        res@cnLevelSelectionMode  = "ExplicitLevels"
        res@cnLevels	= cnlvls
	res@gsnSpreadColors = True
	res@gsnSpreadColorStart = 24
	res@gsnSpreadColorEnd = -26
	;res@cnFillPalette	= "CBR_coldhot"
        res@cnLinesOn 	= False
        res@lbLabelBarOn	= False       	; turn off individual cb's  	
        res@mpCenterLonF         = 180.

	plot0	= gsn_csm_contour_map(xwks,data0a,res)

	resov	= True
       	resov@gsnDraw 	= False           	; don't draw
        resov@gsnFrame 	= False            	; don't advance frame
        resov@cnFillOn	= False
        resov@mpCenterLonF         = 180.
	resov@cnLevelSelectionMode  = "ExplicitLevels"
        resov@cnLevels	= (/0.1/)
	resov@cnLineThicknessF	= 2.5
	resov@cnLineColor	=	"Red"
	resov@cnLineLabelsOn	= False 
  	resov@cnInfoLabelOn	= False	
        resov@lbLabelBarOn	= False       	; turn off individual cb's
	plot0_ov = gsn_csm_contour(xwks,data1,resov)

	overlay(plot0,plot0_ov)
	resP	= True
	resP@gsnMaximize	= True        	; use full page
        resP@gsnPaperOrientation	= "portrait"
	resP@gsnPanelLabelBar	= True   	; add common color bar
	gsn_panel(xwks,plot0,(/1,1/),resP)


	print("Total CPU time: " + get_cpu_time())

end

