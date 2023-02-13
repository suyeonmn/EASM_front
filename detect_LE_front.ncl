;CESM1-LE
load "/home/suyeon/my_function.ncl"

nx = 288
ny = 192
dmiss = -900.
nyr = 2020-1958+1
year = fspan(1958,2020,nyr)
nyr1 = 2005-1958+1
tt = nyr*365
idir = "/data2/suyeon/cesm1le/day"
odir = "/home/suyeon/front/data/CESM1LE/"

expn = "hist"
;expn = "xghg"

if (expn .eq. "hist")then
ensn = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","101","102","103","104","105"/)
else 	;"xghg"
ensn = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020"/)
end if 
Nens = dimsizes(ensn)

nday = 30
do ie = 0,Nens-1
print(ie+" "+ensn(ie))
if (expn.eq. "hist")then
fin1 = addfile(idir+"/T850/b.e11.B20TRC5CNBDRD.f09_g16."+ensn(ie)+".cam.h1.T850.19200101-20051231.nc","r")
fin2 = addfile(idir+"/Q850/b.e11.B20TRC5CNBDRD.f09_g16."+ensn(ie)+".cam.h1.Q850.19200101-20051231.nc","r")
else
fin1 = addfile(idir+"/T850/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+ensn(ie)+".cam.h1.T850.19200101-20051231.nc","r")
fin2 = addfile(idir+"/Q850/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+ensn(ie)+".cam.h1.Q850.19200101-20051231.nc","r")
end if 
T = fin1->T850(13870:31389,:,:)		;[time | 31390] x [lat | 192] x [lon | 288]	;K
Q = fin2->Q850(13870:31389,:,:)		;kg/kg
T@_FillValue = dmiss
Q@_FillValue = dmiss
Tre = reshape(T,(/nyr1,365,ny,nx/))
Qre = reshape(Q,(/nyr1,365,ny,nx/))	;[48] x [365] x [192] x [288]

Tjj = Tre(:,166:195,:,:)	;Only June 16 -July 15
Qjj = Qre(:,166:195,:,:)	
copy_VarCoords(T(0,:,:),Tjj(0,0,:,:))
copy_VarCoords(Q(0,:,:),Qjj(0,0,:,:))
delete(fin1)
delete(fin2)

if (expn .eq. "hist")then
if (ie .lt. 33)then
fin3 = addfile(idir+"/T850/b.e11.BRCP85C5CNBDRD.f09_g16."+ensn(ie)+".cam.h1.T850.20060101-20801231.nc","r")
fin4 = addfile(idir+"/Q850/b.e11.BRCP85C5CNBDRD.f09_g16."+ensn(ie)+".cam.h1.Q850.20060101-20801231.nc","r")
else 	;33~
fin3 = addfile(idir+"/T850/b.e11.BRCP85C5CNBDRD.f09_g16."+ensn(ie)+".cam.h1.T850.20060101-21001231.nc","r")
fin4 = addfile(idir+"/Q850/b.e11.BRCP85C5CNBDRD.f09_g16."+ensn(ie)+".cam.h1.Q850.20060101-21001231.nc","r")
end if 	;historical_ens

else
fin3 = addfile(idir+"/T850/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+ensn(ie)+".cam.h1.T850.20060101-20801231.nc","r")
fin4 = addfile(idir+"/Q850/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+ensn(ie)+".cam.h1.Q850.20060101-20801231.nc","r")
end if 
T2= fin3->T850(0:5474,:,:)
Q2= fin4->Q850(0:5474,:,:)
T2@_FillValue = dmiss
Q2@_FillValue = dmiss
T2re = reshape(T2,(/15,365,ny,nx/))
Q2re = reshape(Q2,(/15,365,ny,nx/))
T2jj = T2re(:,166:195,:,:)
Q2jj = Q2re(:,166:195,:,:)
copy_VarCoords(T2(0,:,:),T2jj(0,0,:,:))
copy_VarCoords(T2(0,:,:),Tjj(0,0,:,:))
delete(fin3)
delete(fin4)

TT = new((/nyr,30,ny,nx/),"float")
QQ = TT
TT(0:nyr1-1,:,:,:) = Tjj
TT(nyr1:nyr-1,:,:,:) = T2jj
copy_VarCoords(T2(0,:,:),TT(0,0,:,:))
;printVarSummary(TT)	;[63] x [30] x [lat | 192] x [lon | 288]
QQ(0:nyr1-1,:,:,:) = Qjj
QQ(nyr1:nyr-1,:,:,:) = Q2jj
copy_VarCoords(Q2(0,:,:),QQ(0,0,:,:))


; Cal_EPT =====================================
press = 85000.
p = new((/nyr,nday,ny,nx/),"float")
p = press
copy_VarCoords(TT,p)
p@units = "Pa"
p@_FillValue = dmiss
ept = pot_temp_equiv(p, TT, QQ, -1, "q") ; dim=-1 or 1
;printVarSummary(ept)	;[86] x [62] x [lat | 192] x [lon | 288]

; dEPT/dy =====================================
lat = dble2flt(T&lat)
lon = dble2flt(T&lon)	;double->float 
;printVarSummary(lat)
ept_grad = grad_latlon_cfd(ept,lat,lon,True,False)
dEPTdy = ept_grad[0]
dEPTdy= dEPTdy*1000.
dEPTdy@units = "K/km"
;printVarSummary(dEPTdy)		;[63] x [30] x [lat | 192] x [lon | 288]

delete(p)
delete(T)
delete(Q)
delete(Tre)
delete(Qre)
delete(T2re)
delete(Q2re)
delete(Tjj)
delete(Qjj)
delete(TT)
delete(T2jj)
delete(Q2jj)
delete(QQ)
delete(ept_grad)

;;================================================================
;; find the front
;;================================================================
abs_dtdy = abs(dEPTdy)
copy_VarCoords(ept,abs_dtdy)	;[86] x [62] x [lat | 192] x [lon | 288]

; Condition 1. 
abs_dtdy = where( abs_dtdy .ge. 0.04 , 1., 0.)  ;front =1 / non-front = 0
;printVarSummary(abs_dtdy)	;[63] x [30] x [lat | 192] x [lon | 288]

; Condition 2.
dtdy_ea = abs_dtdy(:,:,{22:40},{105:145})
ny1 = dimsizes(dtdy_ea(0,0,:,0))
nx1 = dimsizes(dtdy_ea(0,0,0,:))

dtdy_ea2d = reshape(dtdy_ea,(/nyr,nday,ny1*nx1/))
nEA = dim_num_n( dtdy_ea2d .eq. 1., 2)	;nyr x nt

crii = decimalPlaces(nx1*ny1*0.1,0,True)	
cri = floattointeger(crii)

do iy = 0,nyr-1
do t = 0,nday-1
if (nEA(iy,t) .lt. cri) then
	abs_dtdy(iy,t,:,:) = 0. 
end if
end do 
end do 
abs_dtdy@info1 = "dtdy = 0.4"
;printVarSummary(abs_dtdy)	;[86] x [62] x [lat | 192] x [lon | 288]

delete(dEPTdy)
delete(dtdy_ea2d)
delete(nEA)
;;================================================================
; Condition 3. 
print("Start Condition 3")
abs_dtdy(:,:,{35:40},{100:115}) = 0.    
var_ea = abs_dtdy(:,:,{22:40},{105:145})
lat2 = dble2flt(var_ea&lat)
lon2 = dble2flt(var_ea&lon)
var_ea_lat = new((/nyr,nday,ny1,nx1/),"float")
var_ea_lat@_FillValue = dmiss

; Find the median lat. 
do iy = 0,nyr-1
do t = 0,nday-1
do j = 0,ny1-1
do i = 0,nx1-1
if (var_ea(iy,t,j,i) .eq. 1.) then
        var_ea_lat(iy,t,j,i) = lat2(j)
end if
end do
end do
end do
end do

ea_lat_med = dim_median_n(var_ea_lat,2)
copy_VarCoords(var_ea(:,:,0,:),ea_lat_med)
ea_lat_med@_FillValue = dmiss
;printVarSummary(ea_lat_med) 	;[86] x [62] x [lon | 33]

; Check the std. 
ea_lat_diff = new((/nyr,nday,nx1/),"float")
copy_VarCoords(var_ea(:,:,0,:),ea_lat_diff)
ea_lat_diff@_FillValue = dmiss
do iy = 0,nyr-1
do t = 0,nday-1
do i = 0,nx1-2
 if(.not. ismissing(ea_lat_med(iy,t,i)) .and. .not. ismissing(ea_lat_med(iy,t,i+1)))then
 ea_lat_diff(iy,t,i) = abs( ea_lat_med(iy,t,i+1) - ea_lat_med(iy,t,i) )
 else
 ea_lat_diff(iy,t,i) = dmiss
 end if
end do
end do
end do
avg_lat_diff = dim_avg_n_Wrap(ea_lat_diff,2)

CRI = 1.25
do iy = 0,nyr-1
do t = 0,nday-1
if(.not. ismissing(avg_lat_diff(iy,t)) .and.  avg_lat_diff(iy,t) .gt. CRI )then
abs_dtdy(iy,t,:,:) = 0.
end if
end do
end do
abs_dtdy@info3 = "lat_vari_cri = 1.25"
abs_dtdy@info4 = "delete the island-grid"
abs_dtdy!0 = "year"
abs_dtdy!1 = "JJ"
abs_dtdy&year = year
;printVarSummary(abs_dtdy) 	;[year | 63] x [JJ | 30] x [lat | 192] x [lon | 288]

delete(var_ea_lat)
delete(ea_lat_med)
delete(ea_lat_diff)

;;================================================================
; Condition 4.
VAR = abs_dtdy(:,:,{20:50},{100:150})
VAR!0 = "year"
VAR&year = year
NX = dimsizes(VAR(0,0,0,:))
NY = dimsizes(VAR(0,0,:,0))
LAT = dble2flt(VAR&lat)
LON = dble2flt(VAR&lon)

;; Find the LAT
VAR_LAT = new((/nyr,nday,NY,NX/),"float")
do iy = 0,nyr-1
do t = 0,nday-1
do j = 0,NY-1
do i = 0,NX-1
        if (VAR(iy,t,j,i).eq. 1.)then
        VAR_LAT(iy,t,j,i) = LAT(j)
	end if
end do
end do
end do
end do
copy_VarCoords(VAR,VAR_LAT)
VAR_LAT@_FillValue  = dmiss

VAR_LAT = mask(VAR_LAT,VAR_LAT.ge.40.,False)
VAR_LAT = mask(VAR_LAT,VAR_LAT.le.22.,False)


VAR_LATm = dim_median_n(VAR_LAT,2)
VAR_LATmn = new((/nyr,nday,NX/),"float")
VAR_LATmn_num = new((/nyr,nday,NX/),"float")
do i = 3,NX-4
VAR_LATmn(:,:,i) = dim_avg_n_Wrap(VAR_LATm(:,:,i-3:i+3),2)
VAR_LATmn_num(:,:,i) = dim_num_n(.not. ismissing(VAR_LATm(:,:,i-3:i+3)),2)
end do
copy_VarCoords(VAR_LAT(:,:,0,:),VAR_LATmn)
copy_VarCoords(VAR_LAT(:,:,0,:),VAR_LATmn_num)
copy_VarCoords(VAR_LAT(:,:,0,:),VAR_LATm)

;; FINAL VERSION of frontal line
VAR_LATmn = mask(VAR_LATmn,VAR_LATmn_num.le.2,False)
result = VAR_LATmn

;; FINAL VERSION of front grid
diff_lat = new((/nyr,nday,NY,NX/),"float")
diff_lat@_FillValue = dmiss
do iy = 0,nyr-1
do t= 0,nday-1
do i = 0,NX-1
        if(.not. ismissing(result(iy,t,i)))then
        do j = 0,NY-1
        diff_lat(iy,t,j,i) = LAT(j)-result(iy,t,i)
        end do
        end if
end do
end do
end do
diff_lat2 = diff_lat
diff_lat2 = where(diff_lat .ge. -1.25 .and. diff_lat .lt. 1.25, 1.,0.)
diff_lat2 = where(ismissing(diff_lat2),0.,diff_lat2)
copy_VarCoords(VAR_LAT,diff_lat2)
;printVarSummary(diff_lat2)      ;[year | 64] x [JJ | 120] x [lat | 25] x [lon | 41]



;===============================================
;WRITE
fout = odir+"CESM1LE."+expn+"."+ensn(ie)+".front_con1n2n3_lat.1958-2020.JJ30dy.nc"
print(fout)
fo = addfile(fout,"c")
fo->front_lat = result
fo->lon=LON

fout2 = odir+"CESM1LE."+expn+"."+ensn(ie)+".front_con1n2n3_gridline.1958-2020.JJ30dy.nc"
print(fout2)
fo2 = addfile(fout2,"c")
fo2->front = diff_lat2

delete(diff_lat2)
delete(LON)
delete(LAT)
delete(VAR_LAT)
delete(VAR)
delete(VAR_LATm)
delete(VAR_LATmn)
delete(VAR_LATmn_num)

delete(var_ea)
delete(diff_lat)
delete(fout)
delete(fout2)
delete(fo)
delete(fo2)
end do ;ie

