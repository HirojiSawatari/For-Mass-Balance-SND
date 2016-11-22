load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

	d0	= addfile("./balance.nc","r")
	d1	= addfile("./re_sft.nc","r")
	d2  = addfile("./resft.nc","r")
	data0	= d0->sftgif	
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
        res@tiMainString         = "FGOALS-g2 Mass Balance"    ; plot title
	;res@mpMinLatF   = -0.
       	;res@mpMaxLatF   = 50.
        ;res@mpMinLonF   = 0.
        ;res@mpMaxLonF 	= 360.
       ; res@mpCenterLonF	= 180.
		res@gsnDraw 	= False           	; don't draw
        res@gsnFrame 	= False            	; don't advance frame
        res@cnFillOn	= True
        res@cnLevelSelectionMode  = "ExplicitLevels"
        ;res@cnLevels	= (/-3000.,-2500.,-2000.,-1500.,-1000.,-500,-300.,-100.,100.,300.,500.,1000.,1500.,2000.,2500.,3000./)
        res@cnLevels	= cnlvls
		;colormap	= read_colormap_file("CBR_coldhot")
		;tmpc	= colormap(::-1,:)
		res@cnFillPalette	= "CBR_coldhot"
        res@cnLinesOn 	= False
        res@lbLabelBarOn	= False       	; turn off individual cb's
  	
        res@mpCenterLonF         = 180.
	
	;	res@gsnLeftString	= "Precipitation "\
;							  + " GPCP 1981-2010 LTM "
	plot0	= gsn_csm_contour_map(xwks,data0a,res)

		;res@gsnLeftString	= "Precipitation "\
;							  + " CMAP 1981-2010 LTM "
	;plot1	= gsn_csm_contour_map(xwks,p1,res)
	;-------mask(protect) the area that annual precipitation larger than 600mm
	;p0m	= p0	
	;p0m	= where( y0.gt.1.66, p0 ,0)
    ;p0m(36:71,:)	= -p0m(36:71,:)

		resov	= True
        ;resov@mpMinLatF   = -50.
       	;resov@mpMaxLatF   = 50.
        ;resov@mpMinLonF   = 0.
        ;resov@mpMaxLonF 	= 360.
        ;resov@mpCenterLonF	= 180.
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
	;gsnDraw(plot0)
	;gsnFrame(plopt0)
		resP	= True
		resP@gsnMaximize	= True        	; use full page
        resP@gsnPaperOrientation	= "portrait"
		resP@gsnPanelLabelBar	= True   	; add common color bar
	gsn_panel(xwks,plot0,(/1,1/),resP)


	print("Total CPU time: " + get_cpu_time())

end

