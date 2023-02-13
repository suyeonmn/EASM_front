load "/home/suyeon/my_function.ncl"
nx = 100
ny = 60

dir = "/home/suyeon/front/data/"
fin = addfile(dir+"Frontline.APHRODITE_coordinate.1958-2015.0616-0715.6hr.nc","r")
front = fin->front
dmiss = front@_FillValue
ifront = where(ismissing(front),0.,front)	;[time | 6960] x [lat | 60] x [lon | 100]

fin2 = addfile("/data2/suyeon/APHRO/APHRO_MA_050deg_V1101_1958-2015.daily.0616-0715.nc","r")
pr = fin2->precip(:,{20:50},{100:150})
;printVarSummary(pr)	;[time | 1740] x [latitude | 60] x [longitude | 100]
nt1 = dimsizes(pr(:,0,0))
pr4d = new((/nt1,4,ny,nx/),"float")
pr4d(:,0,:,:) = pr
pr4d(:,1,:,:) = pr
pr4d(:,2,:,:) = pr
pr4d(:,3,:,:) = pr
pr3d = reshape(pr4d,(/nt1*4,ny,nx/))
copy_VarCoords(front,pr3d)
pr3d@_FillValue = dmiss
printVarSummary(pr3d)	;[time | 6960] x [lat | 60] x [lon | 100]
nt = dimsizes(pr3d(:,0,0))

; 	250   500   750
cri = (/  5,  10,   15/)

;; Detection of frontal precipitation
fpr = new((/3,nt,ny,nx/),"float")
do ic = 0,2
cc = cri(ic)
;print("cri = "+cc)
do t = 0,nt-1
do i = 0,nx-1
do j = 0,ny-1
if (front(t,j,i) .eq. 1.)then
        if (j.le. cc) then
        fpr(ic,t,0:j+cc,i) = pr3d(t,0:j+cc,i)
        else
        fpr(ic,t,j-cc:j+cc,i) = pr3d(t,j-cc:j+cc,i)
        end if
end if
end do
end do
end do
end do
copy_VarCoords(pr3d,fpr(0,:,:,:))
fpr!0 = "cri"
fpr@info = "200_500_750km"
fpr@_FillValue = dmiss  
fpr@units = "mm/day"
;printVarSummary(fpr)	;[cri | 3] x [time | 6960] x [lat | 60] x [lon | 100]

fpr5d = reshape(fpr,(/3,nt1,4,ny,nx/))
copy_VarCoords(fpr(:,0,:,:),fpr5d(:,0,0,:,:))
fpr5d@_FillValue = dmiss

P =  dim_avg_n_Wrap(fpr5d,2)	;daily mean 
P@info = "200_500_750km"
printVarSummary(P)	;[cri | 3] x [1740] x [lat | 60] x [lon | 100]

;; WRITE
fout = dir+"APHRO.FrontalPR.1958-2015.0616-0715.nc"
print(fout)
fo = addfile(fout,"c")
fo->pr = P

fout2 = dir+"APHRO.FrontalPR.1958-2015.0616-0715.6hr.nc"
print(fout2)
fo2 = addfile(fout2,"c")
fo2->pr = fpr

