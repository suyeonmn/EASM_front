load "/home/suyeon/my_function.ncl"
dir = "/home/suyeon/front/data/"
syear = 2000
lyear = 2010

;; JRA55  & Weather Chart
fin1 = addfile(dir+"JRA55.front_con1n2n3.1958-2021.JJ30dy.artf.final2.nc","r")
var = fin1->front_lat({syear:lyear},:,{110:145})
nyr = dimsizes(var(:,0,0))	;[year | 11] x [JJ | 120] x [lon | 29]
nx1 = dimsizes(var(0,0,:))
dmiss = var@_FillValue 
;printVarSummary(var)	;[year | 4] x [JJ | 120] x [lon | 29]
jra = reshape(var,(/nyr*120,nx1/))
jra_med = dim_median_n(jra,1)
jra_avg = dim_avg_n_Wrap(jra,1)
;printVarSummary(jra_med)

fin2 = addfile(dir+"WC.Front_LAT_EA.2000-2010.JJ30dy.6hr.1-1.nc","r")
var2 = fin2->frontLAT({syear:lyear},:,:,{110:145})
nx2 = dimsizes(var2(0,0,0,:))
;printVarSummary(var2)	;[year | 11] x [day | 30] x [utc | 4] x [lon | 35]
wc = reshape(var2,(/nyr*30*4,nx2/))
wc_med = dim_median_n(wc,1)
wc_avg = dim_avg_n_Wrap(wc,1)
ntime = dimsizes(wc_med)
;printVarSummary(wc_med)

;;===============================================
output1 = wc_med
output2 = jra_med
nbin =25

;; PDF          ;Generates a joint probability density distribution.
opt = True
opt@biny_min = 20.
opt@binx_min = 20.
opt@biny_max = 45.
opt@binx_max = 45.
pdf2 = pdfxy_conform(output1,output2,nbin,nbin,opt)      ;0 = 25 bins being used
printVarSummary(pdf2)   ;[y | 25] x [x | 25]

xx = fspan(1,100,100)
yy = xx

pdf_reg = dble2flt(pdf2)
pdf_x = dble2flt(pdf2&x)
pdf_max = dble2flt(dim_max_n(pdf2, 0))
pdf_max@_FillValue = dmiss
printVarSummary(pdf_max)

reg_y = new((/25/),"float")
reg_y!0 = "y"
reg_y&y = pdf_x
reg_y@_FillValue = dmiss
printVarSummary(reg_y)
do i = 7,16
	do j = 0,24
	if (pdf_reg(j,i).ne. 0 .and. (pdf_reg(j,i) .eq. pdf_max(i))) then	
;	print(pdf_x(j))
	reg_y(i) =dble2flt( pdf_x(j))
	end if 
	end do
end do 


;************************************************
  plot     = new (1, "graphic")
  plot2     = new (1, "graphic")
  wks      = gsn_open_wks("x11","../../pic/Verif.Frontal_LAT.WC_JRA.2PDF")      
  lft = syear+"-"+lyear 
  res      = True
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@cnInfoLabelOn = False
  res@cnFillOn             = True     ; turn on color fill
  res@cnFillPalette        = "WhiteYellowOrangeRed"
  res@cnFillMode           = "RasterFill"       ; Raster Mode
  res@cnLinesOn            = False    ; no contour lines
  res@cnLineLabelsOn       = False    ; no contour line labels
  res@cnInfoLabelOn        = False

  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02
  res@tiYAxisString            = "LAT~B~JRA55~N~ (~S~o~N~N)"
  res@tiXAxisString            = "LAT~B~Weather chart~N~ (~S~o~N~N)" ;Latitude (~S~o~N~N)"

  res@tmXBFormat  = "f"
  res@tmYLFormat  = "f"
  res@tiYAxisFontHeightF = 0.023
  res@tiXAxisFontHeightF = 0.023
  xy_range(res,25,40,25,40)

  res@gsnLeftStringFontHeightF        =  0.025
  res@lbLabelBarOn = True ;False
  res@lbOrientation       = "vertical"
  res@lbLabelFontHeightF = 0.017
  res@lbTitleString = "Probability density"

  res@lbTitlePosition = "Right"
  res@lbTitleDirection     = "Across"
  res@lbTitleAngleF = 270.
  res@lbTitleFontHeightF   = 0.02

  res@cnLineLabelsOn      = False
  cn_range(res,0.3,3,0.3)

  cres                     = True                   ; plot mods desired
  nodraw(cres)
  cres@xyMarkLineModes     = (/"Lines"/)  ; choose which have markers
  cres@xyDashPatterns      = 1                      ; solid line
  cres@xyLineThicknesses   = (/2./)                ; set second line to 2

do t = 0,0
  res@gsnLeftString        = lft(t)
  plot(t) = gsn_csm_contour (wks,pdf2, res)
  plot2(t) = gsn_csm_xy(wks,xx,yy,cres)        ; create plot
  cres@xyLineColor       = "red"
  overlay(plot(t),plot2(t))
end do

;************************************************
; create panel
;************************************************
 resP                    = True                   ; modify the panel plot
 resP                     = True         ; modify the panel plot
 resP@gsnMaximize         = True         ; large format
 resP@gsnPanelLabelBar    = False ; True         ; add common colorbar
 resP@lbLabelAutoStride   = True         ; auto stride on labels
 resP@lbLabelFontHeightF  = 0.014
 resP@gsnPanelYWhiteSpacePercent = 8
 resP@gsnPanelXWhiteSpacePercent = 5
 gsn_panel(wks,plot,(/1,1/),resP)

