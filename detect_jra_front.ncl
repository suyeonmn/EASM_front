load "/home/suyeon/my_function.ncl"

nx = 288
ny = 145
nyr = 2021-1958+1
year = fspan(1958,2021,nyr)
nday = 30*4 
dmiss = -999.
dir = "../data/"
 
fin1 = addfile(dir+"JRA55.ABS_dEPT850dy.1958-2021.JJ30dy.6hr.nc","r")
var = fin1->dtdy 	;[64] x [30] x [4] x [g0_lat_1 | 145] x [g0_lon_2 | 288]
abs_dtdy = reshape(var,(/nyr,30*4,ny,nx/))
copy_VarCoords(var(0,0,0,:,:),abs_dtdy(0,0,:,:))
abs_dtdy@_FillValue = dmiss
delete(fin1)

;;================================================================
;; find the front
;;================================================================
; Condition 1.
abs_dtdy = where( abs_dtdy .ge. 0.04 , 1., 0.)  ;front =1 / non-front = 0
printVarSummary(abs_dtdy)

;;================================================================
; Condition 2.
dtdy_ea = abs_dtdy(:,:,{22:40},{105:145})
ny1 = dimsizes(dtdy_ea(0,0,:,0))
nx1 = dimsizes(dtdy_ea(0,0,0,:))

dtdy_ea2d = reshape(dtdy_ea,(/nyr,nday,ny1*nx1/))
nEA = dim_num_n( dtdy_ea2d .eq. 1., 2)  ;nyr x nt

crii = decimalPlaces(nx1*ny1*0.1,0,True)       ;float
cri = floattointeger(crii)

do iy = 0,nyr-1
do t = 0,nday-1
if (nEA(iy,t) .lt. cri) then
        abs_dtdy(iy,t,:,:) = 0.
end if
end do
end do
;printVarSummary(abs_dtdy)      ;[86] x [62] x [lat | 192] x [lon | 288]
delete(dtdy_ea2d)
delete(nEA)

;;================================================================
; Condition 3.
abs_dtdy(:,:,{35:40},{105:115}) = 0.     
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
abs_dtdy!0 = "year"
abs_dtdy!1 = "JJ"
abs_dtdy&year = year
;printVarSummary(abs_dtdy)      	;[64] x [120] x [g0_lat_1 | 145] x [g0_lon_2 | 288]
delete(var_ea_lat)
delete(ea_lat_med)
delete(ea_lat_diff)

;================================================================
; Condition 4. (Find the LAT)
var_lat = new((/nyr,120,ny,nx/),"float")
do iy = 0,nyr-1
do t = 0,nday-1
do j = 0,ny-1
do i = 0,nx-1
        if (abs_dtdy(iy,t,j,i).eq. 1.)then
        var_lat(iy,t,j,i) = lat(j)
end if
end do
end do
end do
end do
copy_VarCoords(var,var_lat)
var_lat@_FillValue  = dmiss
printVarSummary(var_lat)        ;[year | 64] x [JJ | 120] x [lat | 25] x [lon | 41]

var_lat = mask(var_lat,var_lat.ge.40.,False)
var_lat = mask(var_lat,var_lat.le.22.,False)


var_latm = dim_median_n(var_lat,2)
var_latmn = new((/nyr,120,41/),"float")
var_latmn_num = new((/nyr,120,41/),"float")

do i = 3,nx-4
var_latmn(:,:,i) = dim_avg_n_Wrap(var_latm(:,:,i-3:i+3),2)	;Zonal running avg. 
var_latmn_num(:,:,i) = dim_num_n(.not. ismissing(var_latm(:,:,i-3:i+3)),2)
end do
copy_VarCoords(var_lat(:,:,0,:),var_latmn)
copy_VarCoords(var_lat(:,:,0,:),var_latmn_num)
copy_VarCoords(var_lat(:,:,0,:),var_latm)
;printVarSummary(var_latmn)      ;[year | 64] x [JJ | 120] x [lon | 41]
;printVarSummary(var_latmn_num)  ;[year | 64] x [JJ | 120] x [lon | 41]

var_latmn = mask(var_latmn,var_latmn_num.le.2,False)
result = var_latmn
;printVarSummary(lon)   ;[lon | 41]
;printVarSummary(result)        ;[year | 64] x [JJ | 120] x [lon | 41]

;; Gridded front 
diff_lat = new((/64,120,ny,nx/),"float")
diff_lat@_FillValue = dmiss
do iy = 0,nyr-1
do t= 0,120-1
do i = 0,nx-1
        if(.not. ismissing(result(iy,t,i)))then
        do j = 0,ny-1
        diff_lat(iy,t,j,i) = lat(j)-result(iy,t,i)
        end do
        end if
end do
end do
end do
diff_lat2 = diff_lat
diff_lat2 = where(diff_lat .ge. -1.25 .and. diff_lat .lt. 1.25, 1.,0.)
diff_lat2 = where(ismissing(diff_lat2),0.,diff_lat2)
copy_VarCoords(var_lat,diff_lat2)
printVarSummary(diff_lat2)      ;[year | 64] x [JJ | 120] x [lat | 25] x [lon | 41]

;===============================================
;WRITE
fout = dir+"JRA55.front.1958-2021.JJ30dy.nc"
print(fout)
fo = addfile(fout,"c")
fo->front_lat = result
fo->lon = lon

fout2 = dir+"JRA55.front_grid.1958-2021.JJ30dy.nc"
print(fout2)
fo2 = addfile(fout2,"c")
fo2->front = diff_lat2
exit


