load "/home/suyeon/my_function.ncl"

nx = 288
ny = 192
nx1 = 33
dmiss = -900.
nyr = 2020-1958+1
year = fspan(1958,2020,nyr)
year0 = fspan(1958,2021,2021-1958+1)
year1 = fspan(2007,2010,4)
dir = "/home/suyeon/front/data/CESM1LE/"

;;== Setting ======================
;1le historical 
ensn = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","101","102","103","104","105"/)
Nens = dimsizes(ensn)

;1le "xghg"
ensn2 = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020"/)
Nens2 = dimsizes(ensn2)

nday = 30
var = new((/Nens,nyr,nday,nx1/),"float")
var2 = new((/Nens2,nyr,nday,nx1/),"float")

;;========================================
;; CEMS-LE | historical 
do ie = 0,Nens-1
;print(ie+" "+ensn(ie))
fin = addfile(dir+"CESM1LE.hist."+ensn(ie)+".front_con1n2n3_lat.1958-2020.JJ30dy.nc","r")
var(ie,:,:,:) = fin->front_lat(:,:,{105:145})		;[year | 63] x [JJ | 30] x [lon | 33]
delete(fin)
end do 	;ie 
;printVarSummary(var)				;[42] x [year | 48] x [JJ | 30] x [lon | 33]

;; CEMS-LE | xghg
do ie = 0,Nens2-1
;print(ie+" "+ensn2(ie))
fin2 = addfile(dir+"CESM1LE.xghg."+ensn(ie)+".front_con1n2n3_lat.1958-2020.JJ30dy.nc","r")
var2(ie,:,:,:) = fin2->front_lat(:,:,{105:145})
delete(fin2)
end do
;printVarSummary(var2)				;[20] x [year | 48] x [JJ | 30] x [lon | 33]

;; JRA55  & Weather Chart
fin3 = addfile("../../data/JRA55.front_con1n2n3.1958-2021.JJ30dy.artf.final2.nc","r")
var31 = fin3->front_lat(0:62,:,{105:145})	;[year | 63] x [JJ | 120] x [lon | 33]
var3 = reshape(var31,(/nyr,30,4,nx1/))
jra = dim_median_n(var3,2)
copy_VarCoords(var31(:,0,:),jra(:,0,:))
;printVarSummary(jra)		;[year | 63] x [30] x [lon | 33]

fin5 = addfile("../../data/JRA55.front_con1n2n3_precipitation.1958-2020.JJ30dy.nc","r")
pr = fin5->pr(1,:,:,:,:)
pr = mask(pr,pr.lt.1.,False)
printVarSummary(pr)	;[year | 63] x [JJ | 120] x [lat | 19] x [lon | 33]
PR = dim_avg_n_Wrap(pr,1)
printVarSummary(PR)	;[year | 63] x [lat | 19] x [lon | 33]

fin4 = addfile("../../data/WC.Front_LAT_EA.2000-2010.JJ30dy.6hr.1-1.nc","r")
var4 = fin4->frontLAT(:,:,:,{105:145})
nx2 = dimsizes(var4(0,0,0,:))
;printVarSummary(var2)  ;[year | 11] x [day | 30] x [utc | 4] x [lon | 35]
LON2 = var4&lon
wc = dim_avg_n_Wrap(var4,2) 
printVarSummary(wc)

;;========================================
syr = 2000
lyr = 2010

le1 = var(:,{syr:lyr},:,:)
le2 = var2(:,{syr:lyr},:,:)
jra1 = jra({syr:lyr},:,:)
wc1 = wc({syr:lyr},:,:)
PR1 = PR({syr:lyr},:,:)
selyr = dimsizes(le1(0,:,0,0))

PRmean = dim_avg_n_Wrap(PR1,0)
printVarSummary(PRmean)
var_hist = reshape(le1,(/Nens,selyr*nday,nx1/))
var_xghg = reshape(le2,(/Nens2,selyr*nday,nx1/))
var_jra = reshape(jra1,(/selyr*nday,nx1/)) 	;[time | 7680] x [lon | 33]
var_wc = reshape(wc1,(/selyr*nday,nx2/)) 	;[time | 7680] x [lon | 33]

;;========================================
; Median 
med_hist = dim_median_n(var_hist,1)
med_xghg = dim_median_n(var_xghg,1)
med_jra = dim_median_n(var_jra,0)
med_wc = dim_median_n(var_wc,0)
;printVarSummary(med_hist)		;[42] x [33]

MIN_hist = dim_min_n(med_hist,0)
MIN_xghg = dim_min_n(med_xghg,0)
MAX_hist = dim_max_n(med_hist,0)
MAX_xghg = dim_max_n(med_xghg,0)
mnmx = new((/4,33/),"float")
mnmx(0,:) = MIN_hist
mnmx(1,:) = MAX_hist
mnmx(2,:) = MIN_xghg
mnmx(3,:) = MAX_xghg
;printVarSummary(mnmx)

MED_hist = dim_median_n(med_hist,0)
MED_xghg = dim_median_n(med_xghg,0)
LON = var&lon 
med_jra!0 = "lon"
med_jra&lon = LON

;;corr.
CORR = new((/3/),"float") 
CORR(0) = escorc_n(med_jra,MED_hist,0,0)
CORR(1) = escorc_n(med_jra,MED_xghg,0,0)
CORR(2) = escorc_n(MED_hist,MED_xghg,0,0)
MED = new((/3/),"float")
MED(0) = dim_median_n(med_jra,0) 
MED(1) = dim_median_n(MED_hist,0) 
MED(2) = dim_median_n(MED_xghg,0) 
MED1 = decimalPlaces(MED,1,True)
CORR1 = decimalPlaces(CORR,1,True)
print(MED+" "+CORR)
print(MED1+" "+CORR1)

fin43 = addfile("../../data/JRA55.EPT850.1958-2021.MAY.6hr.nc","r")
ref_cor = fin43->ept(0,:,:)
OUTPUT = new((/145,288/),"float")
copy_VarCoords(ref_cor,OUTPUT)
OUTPUT!0 = "lat"
OUTPUT!1 = "lon"
printVarSummary(OUTPUT)
OUTPUT({22.5:45},{105:145}) = PRmean({22.5:45},{105:145})

;===============================================
;;; Making Plot ;;;
;===============================================
     wks_type = "x11"
;    wks_type = "ps"
    foutname = "../../pic/FigS1b_Frontline.16JUN_15JUL."+syr+"-"+lyr
    npltfrm = 1
    plot = new(npltfrm,"graphic")
    plot1 = new(npltfrm,"graphic")
    plot2 = new(npltfrm,"graphic")
    plot3 = new(npltfrm,"graphic")
    plot4 = new(npltfrm,"graphic")
    plot5 = new(npltfrm,"graphic")
    plot21 = new(npltfrm,"graphic")
    wks = gsn_open_wks(wks_type,foutname)        ; Open an ps workstation.
    nplt = 0
;===============================================
    lft=syr+"-"+lyr
    res = True
    nodraw(res)
    res@gsnMaximize         = True  
    res@cnFillOn = True
    res@cnLinesOn= False
    res@lbLabelBarOn = True

    res@lbOrientation       = "vertical"
    res@lbLabelFontHeightF  = 0.016
    res@lbBoxMinorExtentF = 0.2
    res@lbTitleString = "Precip. (mm/day)"
    res@lbTitlePosition = "Right"
    res@lbTitleDirection     = "Across"   
    res@lbTitleAngleF = 270.
    res@lbTitleFontHeightF   = 0.02    

    ;;; Font
    res@tmXBLabelFontHeightF = 0.02; 0.02;3
    res@tmYLLabelFontHeightF  = 0.02; 3
    res@tmXBTickSpacingF = 10
    res@tmYLTickSpacingF = 5
    res@cnFillPalette = "WhiteBlue" ;"CBR_wet"
    cn_range(res,12,32,4)
    res@cnLineLabelsOn      = False
;===============================================
  xyRes = True
  nodraw(xyRes)
  xyRes@gsnDraw  = False
  xyRes@gsnFrame = False
  xyRes@xyMarkLineMode     = "Lines" ;"MarkLines"
  xyRes@xyMonoDashPattern  = True
  xyRes@xyLineThicknessF   = 4
  xyRes@xyMarker           = 1

  xyRes@xyMarkerColor      = "black"
  xyRes@xyLineColor        = "black"

  do t = 0,npltfrm-1
    map_range(res,23,40,110,145)
    res@gsnLeftString        = lft(t)
    res@gsnLeftStringFontHeightF        =  0.022 ; 3 ;0.02
    plot(t)  = gsn_csm_contour_map(wks,OUTPUT,res)
  xyRes@xyLineColor        = "black"

;;mnmx - shading
  delete(res@xyLineColors)
  res@gsnXYFillColors = "LightPink"
  res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
  plot1(t)  = gsn_csm_xy (wks,LON,mnmx(0:1,:),res) 

  xyRes@xyMarkerColor      = "blue"
  xyRes@xyLineColor        = "blue"
  plot11(t) = gsn_csm_xy(wks,LON,med_xghg,xyRes)

  xyRes@xyLineOpacityF = 1.
  xyRes@xyDashPattern = 3;(/2,3,4,5/)
  xyRes@xyMarkerColor      = "green"
  xyRes@xyLineColor        = "green"
  plot2(t) = gsn_csm_xy(wks,LON,med_jra,xyRes)

  xyRes@xyDashPattern = 0; (/2,3,4,5/)
  xyRes@xyLineColor        = "black"
  plot5(t) = gsn_csm_xy(wks,LON2,med_wc,xyRes)

  xyRes@xyDashPattern = 0
  xyRes@xyMarkerColor      = "blue"
  xyRes@xyLineColor        = "blue"
  plot4(t) = gsn_csm_xy(wks,LON,MED_xghg,xyRes)

  xyRes@xyMarkerColor      = "red"
  xyRes@xyLineColor        = "red"
  plot3(t) = gsn_csm_xy(wks,LON,MED_hist,xyRes)

  xyRes@xyLineColor        = "green"
  plot21(t) = gsn_csm_xy(wks,LON,med_jra,xyRes)


overlay(plot(t),plot2(t))
overlay(plot(t),plot5(t))
overlay(plot(t),plot4(t))
overlay(plot(t),plot3(t))
overlay(plot(t),plot21(t))
end do ;t 

res_text                    = True                  ; text mods desired
res_text@txFontHeightF      = 0.016                 ; change text size
res_text@txJust             = "CenterLeft"          ; text justification
res_text@txFontColor = "black"
dum1 = gsn_add_text(wks,plot,"Weather chart",137.,28.5,res_text)
res_text@txFontColor = "green"
dum1 = gsn_add_text(wks,plot,"JRA55",137.,27,res_text)
res_text@txFontColor = "red"
dum1 = gsn_add_text(wks,plot,"HIST",137.,25.5,res_text)
res_text@txFontColor = "blue"
dum1 = gsn_add_text(wks,plot,"XGHG",137.,24,res_text)

; panel plot only resources
 resP                     = True         ; modify the panel plot
 resP@gsnMaximize         = True         ; large format
 resP@gsnPanelLabelBar    = False ;True         ; add common colorbar
 resP@lbLabelAutoStride   = True         ; auto stride on labels
 resP@lbLabelFontHeightF  = 0.013 ; 0.012
 resP@gsnPanelYWhiteSpacePercent = 8
 resP@gsnPanelXWhiteSpacePercent = 5.6
 resP@pmLabelBarHeightF    = 0.08
 resP@pmLabelBarWidthF     = 0.7 ; 0.4
 gsn_panel(wks,plot,(/2,1/),resP)

frame(wks)




