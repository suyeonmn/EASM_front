load "/home/suyeon/my_function.ncl"

nx = 288
ny = 192
nx1 = 29
dmiss = -900.
syear = 2000
lyear = 2010
nyr = lyear-syear+1
year = fspan(syear,lyear,nyr)
dir = "/home/suyeon/front/data/CESM1LE/"
odir = "/home/suyeon/front/data/"

;;== Setting ======================
regionN = (/"EA","CHI","KOR","JAP"/)
regionNAME = (/"EA","China","Korea","Japan"/)
lon1 = (/110., 115., 125., 130./)
lon2 = (/145., 125., 130., 140./)

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
;; CESM-LE | historical
do ie = 0,Nens-1
;print(ie+" "+ensn(ie))
fin = addfile(dir+"CESM1LE.hist."+ensn(ie)+".front_con1n2n3_lat.1958-2020.JJ30dy.nc","r")
var(ie,:,:,:) = fin->front_lat({syear:lyear},:,{110:145})               ;[year | 63] x [JJ | 30] x [lon | 33]
delete(fin)
end do  ;ie
;printVarSummary(var)                            ;[40] x [year | 63] x [JJ | 30] x [lon | 33]

;; CEMS-LE | xghg
do ie = 0,10 ;Nens2-1
;print(ie+" "+ensn2(ie))
fin2 = addfile(dir+"CESM1LE.xghg."+ensn(ie)+".front_con1n2n3_lat.1958-2020.JJ30dy.nc","r")
var2(ie,:,:,:) = fin2->front_lat({syear:lyear},:,{110:145})
delete(fin2)
end do
;printVarSummary(var2)                           ;[20] x [year | 63] x [JJ | 30] x [lon | 33]

;; JRA55  & Weather Chart
fin3 = addfile("../../data/JRA55.front_con1n2n3.1958-2021.JJ30dy.artf.final2.nc","r")
var31 = fin3->front_lat({syear:lyear},:,{110:145})	;[year | 63] x [JJ | 120] x [lon | 33]
var3 = reshape(var31,(/nyr,30,4,nx1/))
;jra = dim_median_n(var3,2)
jra = dim_avg_n_Wrap(var3,2)
copy_VarCoords(var31(:,0,:),jra(:,0,:))
;printVarSummary(jra)		;[year | 63] x [30] x [lon | 33]

fin4 = addfile("../../data/WC.Front_LAT_EA.2000-2010.JJ30dy.6hr.1-1.nc","r")
var4 = fin4->frontLAT({syear:lyear},:,:,{110:145})
nx2 = dimsizes(var4(0,0,0,:))
;printVarSummary(var4)  ;[year | 11] x [day | 30] x [utc | 4] x [lon | 35]

;;==========================================
HIST = ndtooned(var)
XGHG = ndtooned(var2)
JRA = ndtooned(var31)
WC = ndtooned(var4)
printVarSummary(JRA)

;;====================
bin = 8
opt = True
opt@bin_min = 20.
opt@bin_max = 44.

xx = new((/4,bin/),"double")
yy = xx
meanlat = new((/4/),"float")

;PDF 
ap = pdfx(XGHG, bin, opt)
bp = pdfx(HIST, bin, opt)
cp = pdfx(JRA, bin, opt)
dp = pdfx(WC, bin, opt)
printVarSummary(ap)     ;[x | 40] ; double

xx(0,:) = ap@bin_center
xx(1,:) = bp@bin_center
xx(2,:) = cp@bin_center
xx(3,:) = dp@bin_center
yy(0,:) = (/ap/)
yy(1,:) = (/bp/)
yy(2,:) = (/cp/)
yy(3,:) = (/dp/)

xx@_FillValue= dmiss
yy@_FillValue= dmiss
printVarSummary(xx)
printVarSummary(yy)	

meanlat = new((/4/),"float")
meanlat(0) = dim_avg_n_Wrap(XGHG,0)
meanlat(1) = dim_avg_n_Wrap(HIST,0)
meanlat(2) = dim_avg_n_Wrap(JRA,0)
meanlat(3) = dim_avg_n_Wrap(WC,0)
meanlat2 = decimalPlaces(meanlat,1,True) 

;===============================================
;;; Making Plot ;;;
;===============================================
      wks_type = "x11"
;     wks_type = "ps"
    foutname = "../../pic/FigS1c.Verif_Front_lat.2000-2010.JJ30dy.PDF1d.EA.1LE"
    npltfrm = 1
    plot = new(npltfrm,"graphic")
    wks = gsn_open_wks(wks_type,foutname)        ; Open an ps workstation.
    nplt = 0
;===============================================
lft = (/"EA [110-145E]"/)
res  = True
nodraw(res)
  res@xyLineThicknessF        = 2. ;2.8;(/3,3,3/) ; 2.0,2.0,2.0/)
  res@xyLineColors             = (/"blue","red","green","black"/)
  res@xyMonoDashPattern        = True              ; all solid
  res@tiYAxisString            = "PDF (%)"
  res@tiXAxisString            = "LAT (~S~o~N~N)" ;"Latitude (~S~o~N~N)"
res@tiXAxisFontHeightF = 0.017
res@tiYAxisFontHeightF = 0.017

res@trYMinF   = 0
res@trYMaxF   = 45.;40.

res@trXMinF   = 20.
res@trXMaxF   = 44.
res@tmXBFormat = "f"
res@vpWidthF  = 0.7
res@vpHeightF = 0.25

   res@gsnXRefLineColors =(/"blue","red","green","black"/)
   res@gsnXRefLineDashPatterns = (/13,13,13,13/) ;2,3,4/)
   res@gsnXRefLineThicknessF  = 1.5; (/1.,1,1,1/)

txres1 = True
nodraw(txres1)
txres1@txFontHeightF = 0.017
txres1@txJust        = "CenterLeft"
txres1@txFontColor = "black"

 do t= 0,npltfrm-1
 res@gsnLeftString        = lft(t) 
 res@gsnLeftStringFontHeightF        =  0.018 ;0.022
 res@gsnXRefLine = meanlat ;  WC3(t,:)
 plot(t) = gsn_csm_xy (wks, xx, yy, res)
 end do 


 txres1@txFontHeightF = 0.014
 txres1@txFontColor = "black"
 gsn_text_ndc (wks,"Weather chart ("+meanlat2(3)+")",0.73,0.88,txres1)    
 txres1@txFontColor = "green"
 gsn_text_ndc (wks,"JRA55 ("+meanlat2(2)+")",0.73,0.855,txres1)   
 txres1@txFontColor = "red"
 gsn_text_ndc (wks,"HIST ("+meanlat2(1)+")",0.73,0.83,txres1)     
 txres1@txFontColor = "blue"
 gsn_text_ndc (wks,"XGHG ("+meanlat2(0)+")",0.73,0.805,txres1)     

; panel plot only resources
 resP                     = True         ; modify the panel plot
 resP@gsnMaximize         = True         ; large format
 resP@gsnPanelLabelBar    = True         ; add common colorbar
 resP@lbLabelAutoStride   = True         ; auto stride on labels
 resP@lbLabelFontHeightF  = 0.013 ; 0.012
 resP@gsnPanelYWhiteSpacePercent = 8
 resP@gsnPanelXWhiteSpacePercent = 5.6
 resP@pmLabelBarHeightF    = 0.08
 resP@pmLabelBarWidthF     = 0.8 ; 0.4
 gsn_panel(wks,plot,(/2,1/),resP)

frame(wks)


