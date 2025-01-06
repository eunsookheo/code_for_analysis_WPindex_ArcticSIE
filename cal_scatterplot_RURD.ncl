begin

;Calculate scatterplots b/w WP & PNA/SIE/anticyclone...

;remove AMOC effect
opt	=	"(c)"
;opt	=	""

;AMOC
f_amoc_ru	=	addfile("AMOC_trend_P1.nc","r")
f_amoc_rd	=	addfile("AMOC_trend_P2_90y.nc","r")

amoc_ru	=	f_amoc_ru->var
amoc_rd	=	f_amoc_rd->var

;SIE index
f_sie	=	addfile("SIE_40to90_JJAS.nc","r")
sie	=	f_sie->var

sie_ru  =   regCoef_n(ispan(1,100,1),sie(:,0:99),0,1)
sie_rd  =   regCoef_n(ispan(1,90,1),sie(:,190:279),0,1)



;WP index

f_EN_ru	=	addfile("WP_index_RU_ONDJFM.nc","r")
f_EN_rd	=	addfile("WP_index_RD_ONDJFM.nc","r")

EN_ru	=	f_EN_ru->var
EN_rd	=	f_EN_rd->var

pc_ru	=	EN_ru
;pc_ru	=	sie_ru

pc_rd	=	EN_rd
;pc_rd	=	sie_rd

;remove amoc
pc_amoc_ru =   regCoef_n(amoc_ru,pc_ru,0,0)
pc_amoc_rd =   regCoef_n(amoc_rd,pc_rd,0,0)


do i=0,28-1
pc_y_ru    =   pc_amoc_ru*amoc_ru(i)+pc_amoc_ru@yintercept
pc_y_rd    =   pc_amoc_rd*amoc_rd(i)+pc_amoc_rd@yintercept

pc_ru(i)  =   pc_ru(i)  - pc_y_ru
pc_rd(i)  =   pc_rd(i)  - pc_y_rd
end do


;Z200

f_z200_ru	=	addfile("Z200_trend_RU_ONDJFM.nc","r")
f_z200_rd	=	addfile("Z200_trend_RD_ONDJFM_90y.nc","r")

z200_ru	=	f_z200_ru->var
z200_rd =   f_z200_rd->var


lat	=	z200_ru&lat	
lon	=	z200_ru&lon
nlat	=	dimsizes(lat)
nlon	=	dimsizes(lon)
rad =   4.0 * atan(1.0) / 180.0


;RWS,VP200
f_rws_ru	=	addfile("RWS_trend_RU_ONDJFM.nc","r")
f_rws_rd	=	addfile("RWS_trend_RD_ONDJFM_90y.nc","r")
rwsi_ru	=	dble2flt(f_rws_ru->var(:,{35:45},{120:150}))
rwsi_rd	=	dble2flt(f_rws_rd->var(:,{35:45},{120:150}))
lat_rws	=	rwsi_ru&lat
wgt_rws	=	cos(lat_rws*rad)
rws_ru	=	wgt_areaave_Wrap(rwsi_ru,wgt_rws,1.0,0)
rws_rd  =   wgt_areaave_Wrap(rwsi_rd,wgt_rws,1.0,0)


;CP
f_cp_ru	=	addfile("CP_index_RU_ONDJFM.nc","r")
f_cp_rd	=	addfile("CP_index_RD_ONDJFM.nc","r")
cp_ru	=	f_cp_ru->var
cp_rd	=	f_cp_rd->var





;Local HC
f_hc_ru	=	addfile("LocalHC_trend_RU_ONDJFM.nc","r")
f_hc_rd	=	addfile("LocalHC_trend_RD_ONDJFM.nc","r")
hci_ru	=	f_hc_ru->var(:,{15000:92500},{0:30})
hci_rd	=	f_hc_rd->var(:,{15000:92500},{0:30})

hc_ru	=	dim_avg_n(dim_avg_n(hci_ru,2),1)
hc_rd	=	dim_avg_n(dim_avg_n(hci_rd,2),1)


;FLDS
flux	=	"FLDS"
f_fx_ru	=	addfile(flux+"_trend_RU_ONDJFM.nc","r")
f_fx_rd	=	 addfile(flux+"_trend_RD_ONDJFM_90y.nc","r")
;fxi_ru =   f_fx_ru->var(:,{-10:10},{120:150})
;fxi_rd  =   f_fx_rd->var(:,{-10:10},{120:150})
fxi_ru	=	f_fx_ru->var(:,{75:90},:) ;[28,lev,lat]
fxi_rd  =   f_fx_rd->var(:,{75:90},:) 
lat_fx	=	fxi_ru&lat
wgt_fx	=	cos(lat_fx*rad)

fx_ru	=	wgt_areaave_Wrap(fxi_ru,wgt_fx,1.0,0)
fx_rd   =   wgt_areaave_Wrap(fxi_rd,wgt_fx,1.0,0)

;remove amoc

z200_amoc_ru	=	regCoef_n(amoc_ru,z200_ru,0,0)
fx_amoc_ru	=	regCoef_n(amoc_ru,fx_ru,0,0)
rws_amoc_ru	=	regCoef_n(amoc_ru,rws_ru,0,0)
hc_amoc_ru	=	regCoef_n(amoc_ru,hc_ru,0,0)
cp_amoc_ru  =   regCoef_n(amoc_ru,cp_ru,0,0)

z200_amoc_rd    =   regCoef_n(amoc_rd,z200_rd,0,0)
fx_amoc_rd  =   regCoef_n(amoc_rd,fx_rd,0,0)
rws_amoc_rd =   regCoef_n(amoc_rd,rws_rd,0,0)
hc_amoc_rd =   regCoef_n(amoc_rd,hc_rd,0,0)
cp_amoc_rd =   regCoef_n(amoc_rd,cp_rd,0,0)

do i=0,28-1
z200_y_ru    =   z200_amoc_ru*amoc_ru(i)+onedtond(z200_amoc_ru@yintercept,(/nlat,nlon/))
z200_ru(i,:,:)  =   z200_ru(i,:,:)  - z200_y_ru

z200_y_rd    =   z200_amoc_rd*amoc_rd(i)+onedtond(z200_amoc_rd@yintercept,(/nlat,nlon/))
z200_rd(i,:,:)  =   z200_rd(i,:,:)  - z200_y_rd


fx_y_ru	=	fx_amoc_ru*amoc_ru(i)+fx_amoc_ru@yintercept
fx_ru(i)	=	fx_ru(i) - fx_y_ru

fx_y_rd =   fx_amoc_rd*amoc_rd(i)+fx_amoc_rd@yintercept
fx_rd(i)    =   fx_rd(i) - fx_y_rd


rws_y_ru	=	rws_amoc_ru*amoc_ru(i)+rws_amoc_ru@yintercept
rws_ru(i)    =   rws_ru(i) - rws_y_ru

rws_y_rd    =   rws_amoc_rd*amoc_rd(i)+rws_amoc_rd@yintercept
rws_rd(i)    =   rws_rd(i) - rws_y_rd


hc_y_ru    =   hc_amoc_ru*amoc_ru(i)+hc_amoc_ru@yintercept
hc_ru(i)    =   hc_ru(i) - hc_y_ru

hc_y_rd    =   hc_amoc_rd*amoc_rd(i)+hc_amoc_rd@yintercept
hc_rd(i)    =   hc_rd(i) - hc_y_rd

cp_y_ru    =   cp_amoc_ru*amoc_ru(i)+cp_amoc_ru@yintercept
cp_ru(i)    =   cp_ru(i) - cp_y_ru

cp_y_rd    =   cp_amoc_rd*amoc_rd(i)+cp_amoc_rd@yintercept
cp_rd(i)    =   cp_rd(i) - cp_y_rd


end do


;Make PNA index

region1_ru    =   z200_ru(:,{12.5:22.5},{185:195})
region2_ru    =   z200_ru(:,{40:50},{175:185})
region3_ru    =   z200_ru(:,{65:75},{220:230})
;region1_ru    =   z200_ru(:,{12.5:22.5},{170:200})
;region2_ru    =   z200_ru(:,{40:50},{160:190})
;region3_ru    =   z200_ru(:,{65:75},{180:230})

region1_rd    =   z200_rd(:,{12.5:22.5},{185:195})
region2_rd    =   z200_rd(:,{40:50},{175:185})
region3_rd    =   z200_rd(:,{65:75},{220:230})
;region1_rd    =   z200_rd(:,{12.5:22.5},{170:200})
;region2_rd    =   z200_rd(:,{40:50},{160:190})
;region3_rd    =   z200_rd(:,{65:75},{180:230})


lat1	=	region1_ru&lat
lat2	=	region2_ru&lat
lat3	=	region3_ru&lat

rad	=	4.0 * atan(1.0) / 180.0

wgt1	=	cos(lat1*rad)
wgt2	=	cos(lat2*rad)
wgt3	=	cos(lat3*rad)

r1_ru	=	wgt_areaave_Wrap(region1_ru,wgt1,1.0,0)
r2_ru	=	wgt_areaave_Wrap(region2_ru,wgt2,1.0,0)
r3_ru  =   wgt_areaave_Wrap(region3_ru,wgt3,1.0,0)

r1_rd   =   wgt_areaave_Wrap(region1_rd,wgt1,1.0,0)
r2_rd   =   wgt_areaave_Wrap(region2_rd,wgt2,1.0,0)
r3_rd  =   wgt_areaave_Wrap(region3_rd,wgt3,1.0,0)

PNA_ru	=	(r1_ru-r2_ru+r3_ru)/3.0
PNA_rd =   (r1_rd-r2_rd+r3_rd)/3.0


;Draw scatter

var1_ru	=	pc_ru
var1_rd	=	pc_rd
var2_ru	=	rws_ru
var2_rd	=	rws_rd

varname1	=		"WP_index"
varname2	=	"RWS";flux+"_index"

unit1	=	"[K/year]"
unit2	=	"[s~S~-2~N~/year]"


cor_ru	=	 escorc_n(var1_ru,var2_ru,0,0)
cor_rd  =    escorc_n(var1_rd,var2_rd,0,0)

pr_ru  =   rtest(cor_ru,28,0)
pr_rd	=	rtest(cor_rd,28,0)

ad_ru  =   " "
ad_rd	=	" "
if(pr_ru.le.0.05) then
if(pr_ru.le.0.01) then
ad_ru  =   "**"
else
ad_ru  =   "*" 
end if
end if


if(pr_rd.le.0.05) then
if(pr_rd.le.0.01) then
ad_rd  =   "**"
else
ad_rd  =   "*"
end if
end if


rc_ru  =   regline(var1_ru,var2_ru)
rc_rd	=	regline(var1_rd,var2_rd)
data_ru      = new ( (/2,28/),"float")
data_rd      = new ( (/2,28/),"float")

rstd_ru	=	rc_ru@rstd
rstd_rd	=	rc_rd@rstd
diff	=	abs(rc_ru - rc_rd)
prob	=	diff/sqrt((rstd_ru)^2+(rstd_rd)^2)
print(prob)
add_sig	=	""
if(prob.ge.1.64) then
add_sig	=	"*"
end if

data_ru(0,:) = var2_ru(:)
data_ru(1,:) = rc_ru*(var1_ru(:)-rc_ru@xave) + rc_ru@yave

data_rd(0,:) = var2_rd(:)
data_rd(1,:) = rc_rd*(var1_rd(:)-rc_rd@xave) + rc_rd@yave


wks    =   gsn_open_wks("x11","scatter_"+varname1+"_"+varname2+"_RURD")

plotname   =  opt

xname   =   "ONDJFM "+varname1+" "+unit1
yname   =   "ONDJFM "+varname2+" "+unit2

res                     = True 
res@gsnDraw		=	False                  ; plot mods desired
res@gsnFrame	=	False
res@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
res@xyMarkers           = 16                     ; choose type of marker
res@xyMarkerColor       = "red"                  ; Marker color
res@xyMarkerSizeF       = 0.007                  ; Marker size (default 0.01)
res@xyDashPattern   =   0
res@xyLineColor =   "pink"
res@xyLineThicknesses   =   (/1,2/)
res@trYMaxF =   max(var2_ru)+stddev(var2_ru)
res@trYMinF =   min(var2_ru)-stddev(var2_ru)

res@trXMaxF =   max(var1_ru)+stddev(var1_ru)
res@trXMinF =   min(var1_ru)-stddev(var1_ru)
res@gsnLeftStringFontHeightF    =  18.0
res@gsnLeftString   =   plotname
res@tiXAxisFontHeightF  = 0.015
res@tiYAxisFontHeightF  =   0.015
res@vpWidthF    =   0.35
res@vpHeightF   =   0.35

res@tiXAxisString   =   xname
res@tiYAxisString   =   yname


plot    =   gsn_csm_xy(wks,var1_ru,data_ru,res)

res@xyLineColor =   "skyblue"

plot2	=	gsn_csm_xy(wks,var1_rd,data_rd,res)

pmres   =   True
pmres@gsMarkerIndex =   16
pmres@gsMarkerSizeF =   0.007
pmres@gsMarkerColor =   "blue"

dum =   new((/28/),graphic)

do i=0,28-1

xmarker =   var1_rd(i)
ymarker =   var2_rd(i)


dum(i)  =   gsn_add_polymarker(wks,plot,xmarker,ymarker,pmres)

end do


txres   =   True
xpos    =   0.487
ypos    =    (/0.5,0.475/);(/0.785,0.76/)
colors  =   (/"red","blue"/)
labels  =   (/"Cor="+decimalPlaces(cor_ru,2,True)+ad_ru,"Cor="+decimalPlaces(cor_rd,2,True)+ad_rd/)

do i=0,1

txres@txFontHeightF     =       0.015
txres@txFontColor   =   colors(i)

gsn_text_ndc(wks,labels(i),xpos,ypos(i),txres)
end do
overlay(plot,plot2)
draw(plot)
frame(wks)


end



