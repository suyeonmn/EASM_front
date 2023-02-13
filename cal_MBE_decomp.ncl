;;calculate the Moisture Budget equation (MBE) 

load "/home/suyeon/my_function.ncl"

dmiss = -900.

dir = "/home/suyeon/front/data/CESM1LE/"
fin1 = addfile(dir+"CESM1LE.hist.40ENS.U_V_VAR_IAV_Based_FPR_Box1mmd.1958-2020.nc","r")
fin2 = addfile(dir+"CESM1LE.xghg.20ENS.U_V_VAR_IAV_Based_FPR_Box1mmd.1958-2020.nc","r")
hist1 = fin1->var(:,0:2,0:57,:,:)
xghg1 = fin2->var(:,0:2,{1958:2015},:,:)
copy_VarCoords(xghg1(0,:,:,:,:),hist1(0,:,:,:,:))
printVarSummary(hist1)	;[ens | 40] x [var | 3] x [year | 58] x [lat | 96] x [lon | 288]
printVarSummary(xghg1)	;	       u_v_q

xghg = dim_avg_n_Wrap(xghg1(:,:,{1991:2015},:,:),2)	;Climatology
hist = dim_avg_n_Wrap(hist1(:,:,{1991:2015},:,:),2)	;[ens | 40] x [var | 3] x [lat | 96] x [lon | 288]
xghg_p2 = xghg1(:,:,{1991:2015},:,:)			;IAV 1991-2015 
hist_p2 = hist1(:,:,{1991:2015},:,:)

u_xghg = xghg(:,0,:,:)	;[ens | 20] x [lat | 96] x [lon | 288]
v_xghg = xghg(:,1,:,:)
q_xghg = xghg(:,2,:,:)	;g/kg
u_hist = hist(:,0,:,:)
v_hist = hist(:,1,:,:)
q_hist = hist(:,2,:,:)	
copy_VarCoords(u_xghg,q_xghg)
copy_VarCoords(u_hist,q_hist)
U_xghg = dim_avg_n_Wrap(u_xghg,0)
V_xghg = dim_avg_n_Wrap(v_xghg,0)
Q_xghg = dim_avg_n_Wrap(q_xghg,0)	
;printVarSummary(U_xghg)	;[lat | 96] x [lon | 288]

u_xghg_ano = new((/20,25,96,288/),"float")
v_xghg_ano = u_xghg_ano
q_xghg_ano = u_xghg_ano
u_hist_ano = new((/40,25,96,288/),"float")
v_hist_ano = u_hist_ano
q_hist_ano = u_hist_ano
do iy = 0,25-1
u_xghg_ano(:,iy,:,:) = xghg_p2(:,0,iy,:,:) - u_xghg
v_xghg_ano(:,iy,:,:) = xghg_p2(:,1,iy,:,:) - v_xghg
q_xghg_ano(:,iy,:,:) = xghg_p2(:,2,iy,:,:) - q_xghg
u_hist_ano(:,iy,:,:) = hist_p2(:,0,iy,:,:) - u_hist
v_hist_ano(:,iy,:,:) = hist_p2(:,1,iy,:,:) - v_hist
q_hist_ano(:,iy,:,:) = hist_p2(:,2,iy,:,:) - q_hist
end do ;iy 		;coordinate 	[40] x [25] x [96] x [288]

del_u = dim_avg_n_Wrap(u_hist,0) - dim_avg_n_Wrap(u_xghg,0)
del_v = dim_avg_n_Wrap(v_hist,0) - dim_avg_n_Wrap(v_xghg,0)
del_q = dim_avg_n_Wrap(q_hist,0) - dim_avg_n_Wrap(q_xghg,0)
copy_VarCoords(u_xghg(0,:,:),del_u)
copy_VarCoords(u_xghg(0,:,:),del_v)
copy_VarCoords(u_xghg(0,:,:),del_q)
;printVarSummary(del_u)		;[lat | 96] x [lon | 288]


;Calculations
;; (1) TH	------------------------------------
del_q@_FillValue = dmiss
QgradLatLon  = grad_latlon_cfd (del_q, dble2flt(del_q&lat), dble2flt(del_q&lon), True, False)
dqdy = QgradLatLon[0]
dqdx = QgradLatLon[1]
TH1 = U_xghg*dqdx + V_xghg*dqdy
copy_VarCoords(u_xghg(0,:,:),TH1)
TH1 = TH1*(-1.)

xUgradLatLon  = grad_latlon_cfd (U_xghg, dble2flt(del_q&lat), dble2flt(del_q&lon), True, False)
dudx = xUgradLatLon[1]
xVgradLatLon  = grad_latlon_cfd (V_xghg, dble2flt(del_q&lat), dble2flt(del_q&lon), True, False)
dvdy = xVgradLatLon[0]
TH2 = del_q*dudx + del_q*dvdy
copy_VarCoords(u_xghg(0,:,:),TH2)
TH2 = TH2*(-1.)

TH = TH1+TH2
copy_VarCoords(u_xghg(0,:,:),TH)

TH1_am =  wgt_areaave_Wrap(TH1({28:36},{117:135}),1.,1.,0)
TH2_am =  wgt_areaave_Wrap(TH2({28:36},{117:135}),1.,1.,0)
TH_am =  wgt_areaave_Wrap(TH({28:36},{117:135}),1.,1.,0)


;; (2) DY       ------------------------------------
Q_xghg@_FillValue = dmiss
QgradLatLon_xghg  = grad_latlon_cfd (Q_xghg, dble2flt(Q_xghg&lat), dble2flt(Q_xghg&lon), True, False)
dQdy = QgradLatLon_xghg[0]
dQdx = QgradLatLon_xghg[1]
DY1 = del_u * dQdx + del_v * dQdy
copy_VarCoords(u_xghg(0,:,:),DY1)
DY1 = DY1*(-1.)

dUgradLatLon  = grad_latlon_cfd (del_u, dble2flt(del_q&lat), dble2flt(del_q&lon), True, False)
deludx = dUgradLatLon[1]
dVgradLatLon  = grad_latlon_cfd (del_v, dble2flt(del_q&lat), dble2flt(del_q&lon), True, False)
delvdy = dVgradLatLon[0]
DY2 = Q_xghg*deludx + Q_xghg * delvdy
DY2 = DY2*(-1.)

DY = DY1+DY2
copy_VarCoords(u_xghg(0,:,:),DY)


;; (3) TE       ------------------------------------
uq_p_xghg = u_xghg_ano * q_xghg_ano
vq_p_xghg = v_xghg_ano * q_xghg_ano
uq_p_hist = u_hist_ano * q_hist_ano
vq_p_hist = v_hist_ano * q_hist_ano
uq_p_xghg_clm = dim_avg_n_Wrap(uq_p_xghg,1)
uq_p_xghg_clm_mme = dim_avg_n_Wrap(uq_p_xghg_clm,0)
vq_p_xghg_clm = dim_avg_n_Wrap(vq_p_xghg,1)
vq_p_xghg_clm_mme = dim_avg_n_Wrap(vq_p_xghg_clm,0)
copy_VarCoords(u_xghg(0,:,:),uq_p_xghg_clm_mme)
copy_VarCoords(u_xghg(0,:,:),vq_p_xghg_clm_mme)
uq_p_hist_clm = dim_avg_n_Wrap(uq_p_hist,1)
uq_p_hist_clm_mme = dim_avg_n_Wrap(uq_p_hist_clm,0)
vq_p_hist_clm = dim_avg_n_Wrap(vq_p_hist,1)
vq_p_hist_clm_mme = dim_avg_n_Wrap(vq_p_hist_clm,0)
copy_VarCoords(u_hist(0,:,:),uq_p_hist_clm_mme)
copy_VarCoords(u_hist(0,:,:),vq_p_hist_clm_mme)
del_uq = uq_p_hist_clm_mme - vq_p_xghg_clm_mme 
del_vq = vq_p_hist_clm_mme - vq_p_xghg_clm_mme 
copy_VarCoords(u_hist(0,:,:),del_uq)
copy_VarCoords(u_hist(0,:,:),del_vq)

TE = uv2dv_cfd(del_uq, del_vq, dble2flt(del_uq&lat), dble2flt(del_uq&lon), 3)     ;dv/dy+du/dx-(v/a)*tan(phi)
copy_VarCoords(u_hist(0,:,:),TE)
TE = TE*(-1.)
printVarSummary(TE)

TH = TH*100000.		;10^5
DY = DY*100000.
TE = TE*100000.
MFC = TH+DY+TE
copy_VarCoords(u_hist(0,:,:),MFC)

output = new((/4,96,288/),"float")
output(0,:,:) = MFC
output(1,:,:) = TH
output(2,:,:) = DY
output(3,:,:) = TE
copy_VarCoords(u_hist(0,:,:),output(0,:,:))
output@_FillValue = dmiss


;; Percentage
output_mean = wgt_areaave_Wrap(output(:,{28:36},{117:135}),1.,1.,0)
print(output_mean)

;; Contribution ratio 
output_mean_abs = abs(output_mean)
down = sum(output_mean_abs(1:3))
percen = output_mean / down*100.
print(percen)


;; WRITE
odir = "/home/suyeon/front/data/Figure_data/"
fout = odir+"CESM1LE.MBE_850hPa.Based_FPR_Box1mmd.1991-2015.nc"
print(fout)
fo = addfile(fout,"c")
fo->MFC = output

