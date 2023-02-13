load "/home/suyeon/my_function.ncl"
dmiss = -900.
nyr = 2020-1958+1
year = fspan(1958,2020,nyr)
nx = 33
ny = 32
nt = 30 
widthN = 1	;500km 

dir = "/home/suyeon/front/data/CESM1LE/"

;1le historical
ensn = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","101","102","103","104","105"/)
Nens = dimsizes(ensn)

;1le "xghg"
ensn2 = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020"/)
Nens2 = dimsizes(ensn2)

hist = new((/Nens,nyr,nt,ny,nx/),"float")
xghg = new((/Nens2,nyr,nt,ny,nx/),"float")
xghg_fr = new((/Nens2,nyr,30,nx/),"float")

;;======================= INPUT DATA ===========================================
do ie = 0,Nens-1
fin1 = addfile(dir+"CESM1LE.hist."+ensn(ie)+".front_con1n2n3_precipitation.1958-2020.JJ30dy.nc","r")
var1 = fin1->pr		;[cri | 3] x [year | 63] x [JJ | 30] x [lat | 32] x [lon | 33]
hist(ie,:,:,:,:) = var1(widthN,:,:,:,:)
delete(fin1)
delete(var1)
end do 
hist = hist*1000.*86400.                ;Unit = m/s => mm/day
hist!0 = "ens"
;printVarSummary(hist)	;[ens | 40] x [year | 63] x [JJ | 30] x [lat | 19] x [lon | 33]

do iee = 0,Nens2-1
fin2 = addfile(dir+"CESM1LE.xghg."+ensn2(iee)+".front_con1n2n3_precipitation.1958-2020.JJ30dy.nc","r")
var2 = fin2->pr
xghg(iee,:,:,:,:) = var2(widthN,:,:,:,:)
delete(fin2)
delete(var2)

fin21 = addfile("/home/suyeon/front/data/CESM1LE/frontlineV2/CESM1LE.xghg."+ensn2(iee)+".front_con1n2n3_lat.1958-2020.JJ30dy.nc","r")
xghg_fr(iee,:,:,:) = fin21->front_lat(:,:,{105:145})
delete(fin21)
end do
xghg = xghg*1000.*86400.
xghg!0 = "ens" 		;[ens | 40] x [year | 63] x [JJ | 30] x [lat | 32] x [lon | 33]



;;======================= (1) Set the time period ==============================
; SET THE PERIOD
syear = 1991
lyear = 2015

hist_p1 = hist(:,{syear:lyear},:,:,:)
xghg_p1 = xghg(:,{syear:lyear},:,:,:)
xghg_fr_out = xghg_fr(:,{syear:lyear},:,:)
nyr1 = dimsizes(hist_p1(0,:,0,0,0))
nx1 = dimsizes(xghg_fr(0,0,0,:))
xghg_fr1 = reshape(xghg_fr_out,(/20,nyr1*30,nx1/))
xghg_fr2 = dim_median_n(xghg_fr1,0)
xghg_front = dim_avg_n_Wrap(xghg_fr2,0)
LON = xghg_fr&lon 
LAT = xghg_front

cri = "mean_1mmd"
hist_p1 = mask(hist_p1,hist_p1.lt. 1. ,False)
xghg_p1 = mask(xghg_p1,xghg_p1.lt. 1. ,False)
;printVarSummary(hist_p1)	;[ens | 40] x [year | 63] x [JJ | 30] x [lat | 19] x [lon | 33]

hist_p1r = reshape(hist_p1,(/Nens,nyr1*nt,ny,nx/))
xghg_p1r = reshape(xghg_p1,(/Nens2,nyr1*nt,ny,nx/))
nt2 = dimsizes(hist_p1r(0,:,0,0))
copy_VarCoords(hist_p1(:,0,0,:,:),hist_p1r(:,0,:,:))
copy_VarCoords(xghg_p1(:,0,0,:,:),xghg_p1r(:,0,:,:))
hist_p1r@_FillValue = dmiss
xghg_p1r@_FillValue = dmiss
;printVarSummary(hist_p1r)	;[ens | 40] x [1890] x [lat | 32] x [lon | 33]

;;======================= (2) Start decomposing ==============================
;; Del = HIST - XGHG 	| Output is an ensemble avg. 
;; P = total precipitation 
p_h = dim_sum_n_Wrap(hist_p1r,1)	;[ens | 40] x [lat | 32] x [lon | 33]
P_h = dim_avg_n_Wrap(p_h,0)		;[lat | 32] x [lon | 33]
p_x = dim_sum_n_Wrap(xghg_p1r,1)
P_x = dim_avg_n_Wrap(p_x,0)		;[lat | 32] x [lon | 33]
P_hs = dim_variance_n_Wrap(p_h,0)
P_xs = dim_variance_n_Wrap(p_x,0)
probP = ttest(P_h,P_hs,40,P_x,P_xs,20,False,False)
probP = 1. - probP

;; I = Intensity (mm/day) | I -> II	| I & DelI
I_x = dim_avg_n_Wrap(xghg_p1r,1)
I = dim_avg_n_Wrap(I_x,0)		;XGHG_Intensity 
I_h = dim_avg_n_Wrap(hist_p1r,1)
II_h = dim_avg_n_Wrap(I_h,0)		;HIST_Intensity 
DelI = II_h - I
copy_VarCoords(I,DelI)

I_hs = dim_variance_n_Wrap(I_h,0)
I_xs = dim_variance_n_Wrap(I_x,0)
probI = ttest(II_h,I_hs,40,I,I_xs,20,False,False)
probI = 1. - probI

;; N = Numb. (Freq)	| N & DelN
n_x = int2flt(dim_num_n(.not.ismissing(xghg_p1r),1))
copy_VarCoords(I,n_x(0,:,:))
N = dim_avg_n_Wrap(n_x,0)

n_h = int2flt(dim_num_n(.not.ismissing(hist_p1r),1))
copy_VarCoords(I,n_h(0,:,:))
N_h =  dim_avg_n_Wrap(n_h,0)
DelN = N_h - N
N_hs = dim_variance_n_Wrap(n_h,0)
N_xs = dim_variance_n_Wrap(n_x,0)
probN = ttest(N_h,N_hs,40,N,N_xs,20,False,False)
probN = 1. - probN
copy_VarCoords(N,DelN)

;;======================= (3) Calculate decomposing ==============================
output = new((/3,ny,nx/),"float")
output(0,:,:) = I*DelN
output(1,:,:) = N*DelI
output(2,:,:) = DelN*DelI
copy_VarCoords(DelN,output(0,:,:))
printVarSummary(output)


BOXam = wgt_areaave_Wrap(output(:,{28:35},{117:135}), 1.0, 1.0, 0)
abs_am = abs(BOXam)
total_abs = sum(abs_am)
per_am = BOXam/total_abs*(100.)
print(BOXam+" "+abs_am)
print(total_abs)
print(per_am)
exit


