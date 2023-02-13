load "/home/suyeon/my_function.ncl"

dmiss = -900.
nyr = 2015-1958+1
year = fspan(1958,2015,nyr)
ny = 60 ; 19
nx = 100 ; 33
nt = 30

widthN = 1      ;500km
;;===============================================
fin1 = addfile("../../data/APHRO.FrontalPR.1958-2015.0616-0715.6hr.nc","r")
jra = fin1->pr(widthN,:,:,:)
jra_r = reshape(jra,(/nyr,30,4,ny,nx/))
jra_r1 = dim_avg_n_Wrap(jra_r,2)
copy_VarCoords(jra(0,:,:),jra_r1(0,0,:,:))
jra_r1!0 = "year"
jra_r1&year = year      ;[year | 58] x [30] x [lat | 60] x [lon | 100]

jra_r1  = mask(jra_r1, jra_r1 .lt. 1., False)

jra_p1 = jra_r1({1958:1982},:,:,:)
jra_p2 = jra_r1({1991:2015},:,:,:)
nyr1 = dimsizes(jra_p1(:,0,0,0))
jra_p1d = reshape(jra_p1,(/nyr1*30,ny,nx/))
jra_p2d = reshape(jra_p2,(/nyr1*30,ny,nx/))
copy_VarCoords(jra_r1(0,0,:,:),jra_p1d(0,:,:))
copy_VarCoords(jra_r1(0,0,:,:),jra_p2d(0,:,:))
printVarSummary(jra_p2d)	;[900] x [lat | 19] x [lon | 33]


;Frontline
fin2 = addfile("../../data/JRA55.front_con1n2n3.1958-2021.JJ30dy.artf.final2.nc","r")
front = fin2->front_lat(:,:,{110:145})	;[year | 64] x [JJ | 120] x [lon | 29
front2dy = reshape(front,(/64,30,4,29/))
frontDY = dim_avg_n_Wrap(front2dy,2)
copy_VarCoords(front(:,0,:),frontDY(:,0,:))
front1 = frontDY({1958:1982},:,:)
front2 = frontDY({1991:2015},:,:)
front_re1 = reshape(front1,(/nyr1*30,29/))
front_re2 = reshape(front2,(/nyr1*30,29/))
copy_VarCoords(front(0,0,:),front_re1(0,:))
copy_VarCoords(front(0,0,:),front_re2(0,:))

front_re = reshape(front2,(/nyr1,30,29/))
front_dy = dim_median_n(front_re,2)
front_re2 = reshape(front_re,(/nyr1*30,29/))
front_med = dim_median_n(front_re2,0)
printVarSummary(front_med)      ;:34:
lon = front&lon
;printVarSummary(lon)		;FRONTLINE = lon,front_med (<- lat)

;==============================================
;; Decomp. 
nn1 = dimsizes(jra_p1d(:,0,0))
nn2 = dimsizes(jra_p2d(:,0,0))
pr_p1 = dim_sum_n_Wrap(jra_p1d,0)
pr_p2 = dim_sum_n_Wrap(jra_p2d,0)
pr_p1_std = dim_variance_n_Wrap(jra_p1d,0)
pr_p2_std = dim_variance_n_Wrap(jra_p2d,0)
Term1 = pr_p2 - pr_p1
Term1p = ttest(pr_p2,pr_p2_std,nn2,pr_p1,pr_p1_std,nn1,False,False)
Term1p = 1.-Term1p
copy_VarCoords(pr_p1,Term1p)
printVarSummary(Term1p)

; Intensity effect
I_p1 = dim_avg_n_Wrap(jra_p1d,0)
I_p2 = dim_avg_n_Wrap(jra_p2d,0)
delI = I_p2 - I_p1
Ip = ttest(I_p2,pr_p2_std,nn2,I_p1,pr_p1_std,nn1,False,False)
Ip = 1.-Ip
copy_VarCoords(pr_p1,Ip)
printVarSummary(Ip)

; Frequency effect
F_p1 = int2flt(dim_num_n(.not.ismissing(jra_p1d),0))
F_p2 = int2flt(dim_num_n(.not.ismissing(jra_p2d),0))
delF = F_p2 - F_p1
copy_VarCoords(pr_p1,F_p2)
copy_VarCoords(pr_p1,F_p1)
copy_VarCoords(pr_p1,delF)
printVarSummary(F_p2)
Fp = ttest(F_p2,pr_p2_std,nn2,F_p1,pr_p1_std,nn1,False,False)
Fp = 1.-Fp

output = new((/4,ny,nx/),"float")
outputp = new((/4,ny,nx/),"float")
output(0,:,:) = Term1
output(1,:,:) = I_p1*delF
output(2,:,:) = F_p1*delI
output(3,:,:) = delI*delF
copy_VarCoords(pr_p1,output(0,:,:))

outputp(0,:,:) = Term1p
outputp(1,:,:) = Fp
outputp(2,:,:) = Ip
outputp(3,:,:) = 0.
copy_VarCoords(pr_p1,outputp(0,:,:))
printVarSummary(output)

output_am = wgt_areaave_Wrap(output(:,{28:36},{117:135}),1.0,1.0, 0)

; Contribution ratio
output_mean_abs = abs(output_am)
down = sum(output_mean_abs(1:3))
percen = output_am / down*100.
print(percen)
exit







;debug_plot(output(2,:,:),200)
;exit
;;===============================================
	mode = "Intensity"
;;===============================================
lon_aph = fspan(0.25,359.75,720)
lat_aph = fspan(-89.75,89.75,360)
lon_aph@units = "degrees_east"
lat_aph@units = "degrees_north"
dmisss = jra@_FillValue
OUTPUT = new((/4,360,720/),"float")
OUTPUT!1 = "lat"
OUTPUT!2 = "lon"
OUTPUT&lat = lat_aph
OUTPUT&lon = lon_aph
OUTPUT@_FillValue = dmisss
OUTPUTp = OUTPUT
printVarSummary(OUTPUT)
OUTPUT(:,{22:45},{105:145}) = output(:,{22:45},{105:145})
OUTPUTp(:,{22:45},{105:145}) = outputp(:,{22:45},{105:145})
OUTPUT@info = "P_IdF_FdI_dIdF"

;; WRITE for PYTHON FIGURES
odir = "/home/suyeon/front/data/Figure_data/"
fout = odir+"Fig1_APHRO_Decomp_Diff.P2-P1.nc"
print(fout)
fo = addfile(fout,"c")
fo->diff = output
fo->prob = outputp
print("GREAT!!!!!!!!!!!!")


exit



;===============================================
;;      Plot
;;===============================================
      wks_type = "x11"
;     wks_type = "ps"
    foutname = "../../pic/APHRO.decomp_delI.PR."+nyr1+"yr"
    npltfrm = 1
    plot = new(npltfrm,"graphic")
    lfplot = new(npltfrm,"graphic")
    plot3 = new(npltfrm,"graphic")
    plot4 = new(npltfrm,"graphic")
    wks = gsn_open_wks(wks_type,foutname)        ; Open an ps workstation.
    nplt = 0
;===============================================
lft = (/"APHRODITE","Mean > 1mm/day","R99p","R997p"/)
;lft= (/"~F33~D~F~PR","I~F33~D~F~F","F~F33~D~F~I","~F33~D~F~F~F33~D~F~I"/)
    res = True
    nodraw(res)
    res@cnFillOn = True
    res@cnLinesOn= False
   res@cnFillMode = "RasterFill"
   res@lbLabelBarOn = True 
    res@lbOrientation       = "vertical"
    res@cnLineLabelsOn      = False
;    res@cnFillPalette = "precip_11lev" ; "precip_diff_12lev"
    res@cnFillPalette = "precip_diff_12lev"
;    cn_range(res,-2000,2000,400)
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-2000,-1600,-1200,-800,-400,-100,0,100,400,800,1200,1600,2000/)

    ;;; Font
    res@tmXBLabelFontHeightF = 0.02; 0.02;3
    res@tmYLLabelFontHeightF  = 0.02; 3
    res@tmXBTickSpacingF = 10
    res@tmYLTickSpacingF = 5 ;10

;===============================================
  cres                      = True
  nodraw(cres)
  cres@cnFillOn = False
  cres@cnLinesOn= True
  cres@cnLineLabelsOn      = False
  cres@cnInfoLabelOn       = False
  cres@gsnContourNegLineDashPattern =1
  cres@gsnContourZeroLineThicknessF = 0; 2.1
;  cres@cnLevelSelectionMode   = "ExplicitLevels"
  cres@cnLineColor        = "black" 
;  cres@cnLevel = (/50./)
  cres@cnLineThicknessF   = 2.2 ; 2.5
;===============================================
res3                    = True
 nodraw(res3)
 res3@gsnSpreadColors      = False        ; spread out color table
 res3@gsnAddCyclic         = False        ; plotted dataa are not cyclic
 res3@gsnMaximize = True

; ... contour
 res3@cnFillOn             = True         ; turn on color fill
 res3@cnLinesOn            = False        ; True is default
 res3@cnLineLabelsOn       = False
 res3@cnInfoLabelOn        = False

 res3@lbLabelBarOn         = False
 res3@cnLevelSelectionMode = "ExplicitLevels"   ; manually specify contour levels
; res3@cnLevels = (/0.9/) ; min level
 res3@cnLevels = (/0.95/) ; min level
 res3@cnFillColors         = (/"transparent","black"/)
 res3@cnMonoFillPattern    = False            ; want multiple patterns

 res3@cnFillPatterns       = (/-1,3/) ; the patterns
; res3@cnFillPatterns       = (/17,-1/) ; the patterns
 res3@cnMonoFillScale      = False            ; want different densities
 res3@cnFillScales         = (/.1,0.8/) ; the densities
 res3@cnFillDotSizeF       = 0.004            ; increase dot size of contours
;===============================================
  xyRes = True
  nodraw(xyRes)
  xyRes@gsnDraw  = False
  xyRes@gsnFrame = False
  xyRes@xyMarkLineMode     = "Lines" ;"MarkLines"
  xyRes@xyMonoDashPattern  = True
  xyRes@xyLineThicknessF   =  7
;  xyRes@xyMarker           = 1
;  xyRes@xyMarkerColor      = "black"
  xyRes@xyLineColor        = "black"


 hres  = True
nodraw(hres)
  hres@xyLineThicknesses        = (/4,4/) ; 2.0,2.0,2.0/)
  hres@xyLineColors             = (/"black","blue"/)
  hres@xyDashPatterns = (/2,0/)
;  hres@xyMonoDashPattern        = True              ; all solid
;  hres@tiYAxisString            = "PDF (%)"
  hres@tiXAxisString            = "PDF (%)" ;Times per 10-yr" ;Latitude (~S~o~N~N)"
hres@tiXAxisLabelFontHeightF = 0.02
; hres@trYMinF   = 0
; hres@trYMaxF   = ymax(t)

;hres@trXMinF   = 0.
;hres@trXMaxF   = 1.
  hres@tmXBFormat = "f"
  hres@vpWidthF  = 0.1;7
  hres@vpHeightF = 0.2

;   hres@gsnXRefLine = (/nEA_histMEAN(t),nEA_xghgMEAN(t)/)
;   hres@gsnYRefLine = (/ref_fr1,ref_fr2/)
   hres@gsnYRefLineColors = (/"black","blue"/)
   hres@gsnXRefLineDashPattern = 2
   hres@gsnXRefLineThicknessF = 3
;hres@gsnXYAboveFillColors = "red"
;hres@gsnXYBelowFillColors = "green"
; hres@gsnXYFillColors = (/"red","green"/)

;  hres@gsnLeftString        = varn(t)+"~B~"+CASE
  hres@gsnLeftStringFontHeightF        =  0.02
;  hres@gsnRightString        = "HIST-XGHG: "+rft(t)+" "
  hres@gsnRightStringFontHeightF        =  0.016


;  lfplot = gsn_csm_xy (wks, yy,xx, hres)
;  lfplot = gsn_csm_xy (wks, xx,yy, hres)


;************************************************
; create points for box
;************************************************
 ypts = (/ 32.5,  32.5,  25.0,   25.0, 32.5/)
 xpts = (/ 117.0,  137.0,  137.0,  117.0, 117.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "deeppink"                     ; color of lines
  resp@gsLineThicknessF = 3.3; 3.30                       ; thickness of lines
;  resp@gsLineDashPattern = 2.
  resp@gsLineLabelString= "";"test"                    ; adds a line label string
;************************************************
 dum = new(npltfrm,graphic)


;************************************************
if(mode .eq. "Frequency")then
maxv = (/60,9,5/)
minv = (/10.,1,1/) ; maxv*(-1.)
intv = (/10,2,1/)
else    ;intensity
maxv = (/50, 120, 140/)
minv = (/10., 10., 10./) ; maxv*(-1.)
intv = (/5, 10, 10/)
end if

cminv = (/20, 80, 80/)
cmaxv = (/60, 160, 160/)
cintv = (/2, 60, 60/)

do t=0,npltfrm-1
;    map_range(res,22,45,105,145)
    map_range(res,25,40,110,145)

    res@mpGeophysicalLineThicknessF  = 1.5 ; 2.6
;   res@mpGeophysicalLineColor = "black"
    res@gsnLeftString        = lft(t)
    res@gsnLeftStringFontHeightF        =  0.025 ;0.02
    plot(t)  = gsn_csm_contour_map(wks,OUTPUT(2,:,:),res)
delete(res@gsnLeftString)
    plot4(t)  = gsn_csm_contour(wks,OUTPUTp(2,:,:),res3)
    overlay(plot(t),plot4(t))
    plot3(t) = gsn_csm_xy(wks,lon,front_med,xyRes)

   ;dum(t)=gsn_add_polyline(wks,plot(t),xpts,ypts,resp)
;    overlay(plot(t),plot4(t))
    overlay(plot(t),plot3(t))
end do


;newplot = gsn_attach_plots(plot,(/lfplot/),hres,res)
;draw(plot)
;frame(wks)
 

; panel plot only resources
 resP                     = True         ; modify the panel plot
; resP@txString            = "APHRODITE "+nyr1+"yr" 	;mode+"~B~JRA55"

 resP@gsnMaximize         = True         ; large format
 resP@gsnPanelLabelBar    = False ;True         ; add common colorbar
 resP@lbLabelAutoStride   = True         ; auto stride on labels
; resP@lbOrientation       = "vertical"
 resP@lbLabelFontHeightF  = 0.013 ; 0.012
 resP@gsnPanelYWhiteSpacePercent = 8
 resP@gsnPanelXWhiteSpacePercent = 5.6
 resP@pmLabelBarHeightF    = 0.07
 resP@pmLabelBarWidthF     = 0.9 ; 0.4
 gsn_panel(wks,plot,(/2,1/),resP)
; gsn_panel(wks,plot,(/1,3/),resP)

frame(wks)


