begin

expname =   "RU_ONDJFM"
;expname	=	"RD_ONDJFM_90y"
drmod   =   "x11"

varname	=	"Z200"


f_sie	=	addfile("SIE_40to90_JJAS.nc","r")
sie	=	f_sie->var


sie_ru	=	regCoef_n(ispan(1,100,1),sie(:,0:99),0,1)
sie_rd	=	regCoef_n(ispan(1,90,1),sie(:,190:279),0,1)

f_amoc =   addfile("AMOC_trend_P1.nc","r")
;f_amoc	=	addfile("AMOC_trend_P2_90y.nc","r")
amoc	=	f_amoc->var


;pc	=	sie_ru;*(-1.0)
;pc	=	sie_rd

;pc	=	gmst
;pc	=	amoc

;f_EN	=	addfile("WP_index_RU_ONDJFM.nc","r")
;f_EN   =   addfile("WP_index_RD_ONDJFM.nc","r")
f_EN	=	addfile("WP_index_RU_JJA.nc","r")
;f_EN	=	addfile("WP_index_RD_JJA_90y.nc","r")


;f_EN	=	addfile("OLR_index2_RD.nc","r")
EN	=	f_EN->var

;pc	=	EN*(-1.0)
pc	=	EN


rad  = 4.0 * atan(1.0) / 180.0

;f_RWS  =   addfile("RWS_trend_"+expname+".nc","r")
;RWSi   =   f_RWS->var(:,{35:45},{120:140})
;RWSi	=	f_RWS->var(:,{30:60},{180:230})
;wlat   =   RWSi&lat
;wgt    =   cos(wlat*rad)
;RWS    =   wgt_areaave(RWSi,wgt,1.0,1)

;pc	=	RWS


;remove amoc

pc_amoc =   regCoef_n(amoc,pc,0,0)

do i=0,28-1
pc_y    =   pc_amoc*amoc(i)+pc_amoc@yintercept
pc(i)  =   pc(i)  - pc_y
end do


;f_ts	=	addfile(varname+"_trend_"+expname+".nc","r")
f_ts	=	addfile(varname+"_trend_RU_JJA.nc","r")

tsi =   f_ts->var(:,:,:)

lat    =   tsi&lat
lon    =   tsi&lon
nlat   =   dimsizes(lat)
nlon   =   dimsizes(lon)

f_fx	=	addfile("Fx_monthly_"+expname+".nc","r")
f_fy	=	addfile("Fy_monthly_"+expname+".nc","r")

fx_i	=	f_fx->Fx
fy_i	=	f_fy->Fy

fx	=	regCoef_n(ispan(1,100,1),fx_i,0,1)
fy	=	 regCoef_n(ispan(1,100,1),fy_i,0,1) ;fx,fy trend
;fx =   regCoef_n(ispan(1,90,1),fx_i,0,1)
;fy =    regCoef_n(ispan(1,90,1),fy_i,0,1) ;fx,fy trend


;Regression

;reg_ts	=	regCoef_n(pc,tsi,0,0)
;reg_fx	=	regCoef_n(pc,fx,0,0)
;reg_fy	=	regCoef_n(pc,fy,0,0)

ts_amoc	=	regCoef_n(amoc,tsi,0,0)
fx_amoc	=	regCoef_n(amoc,fx,0,0)
fy_amoc	=	regCoef_n(amoc,fy,0,0)

do i=0,28-1

ts_y	=	ts_amoc*amoc(i)+onedtond(ts_amoc@yintercept,(/nlat,nlon/))
fx_y	=	fx_amoc*amoc(i)+onedtond(fx_amoc@yintercept,(/nlat,nlon/))
fy_y    =   fy_amoc*amoc(i)+onedtond(fy_amoc@yintercept,(/nlat,nlon/))

tsi(i,:,:)	=	tsi(i,:,:)	- ts_y
fx(i,:,:)  =   fx(i,:,:)  - fx_y
fy(i,:,:)  =   fy(i,:,:)  - fy_y

end do

reg_ts =   regCoef_n(pc,tsi,0,0)
reg_fx =   regCoef_n(pc,fx,0,0)
reg_fy =   regCoef_n(pc,fy,0,0)
;reg_ts	=	escorc_n(pc,tsi,0,0)
;reg_fx	=	escorc_n(pc,fx,0,0)
;reg_fy  =   escorc_n(pc,fy,0,0)

;====remove zonal mean========

;z_m	=	dim_avg_n(reg_ts,1) ;[lat]

;do i=0,nlat-1
;reg_ts(i,:)	=	reg_ts(i,:) - z_m(i)
;end do


reg_ts!0	=	"lat"
reg_ts!1	=	"lon"
reg_ts&lat	=	lat
reg_ts&lon	=	lon

copy_VarCoords(reg_ts,reg_fx)
copy_VarCoords(reg_ts,reg_fy)

tval    =   onedtond(reg_ts@tval,dimsizes(reg_ts))
df  =   onedtond(reg_ts@nptxy,dimsizes(reg_ts))
b   =   tval
b   =   0.5 
prob    =   (1-betainc(df/(df+tval^2),df/2.0,b))

;prob	=	rtest(reg_ts,28,0)

copy_VarCoords(reg_ts,prob)




wks        =   gsn_open_wks(drmod,"REG_"+varname+"trend_WP_"+expname+"_rAMOC")
;wks        =   gsn_open_wks(drmod,"Cor_"+varname+"trend_WP_"+expname+"_rAMOC")

res1     =       True
res1@cnFillOn    =       True
res1@cnLinesOn   =       False
res1@gsnDraw     =       False
res1@gsnFrame    =       False
res1@cnInfoLabelOn       =       False
res1@lbLabelBarOn   =   False
res1@cnLineLabelsOn      =       False
res1@cnMonoFillColor = True
res1@cnMonoFillPattern  =   False
res1@cnLevelSelectionMode= "ExplicitLevels"
res1@cnLevels   =   (/0.95/)
;res1@cnLevels	=	(/0.05/)
res1@cnFillPatterns =   (/-1,3/)
;res1@cnFillPatterns =   (/3,-1/)

plot1   =       gsn_csm_contour(wks,prob(:,:),res1)


maxval	=50
minval	=-50
space	=10

 res2=True
 res2@gsnDraw	=	False
 res2@gsnFrame	=	False
 res2@cnFillOn=True
 res2@cnLinesOn=False
 res2@cnLineLabelsOn=False
 res2@gsnAddCyclic=True
 res2@mpMinLonF=0
 res2@mpMaxLonF=360
 res2@mpMinLatF=-90
 res2@mpMaxLatF=90
 res2@gsnMaximize       =       True
 res2@mpCenterLonF=180.
 res2@cnLevelSelectionMode =    "ManualLevels"
; res2@cnMaxLevelValF    =maxval;500;3.*10^(-8);500;40;2.6; 
; res2@cnMinLevelValF    =minval;-500;-3.*10^(-8);-500;-40;-2.6;
; res2@cnLevelSpacingF   =space;50;0.5*10^(-8);50;5;0.2


 res2@tmXBLabelFontHeightF   = .018         ; change font height
 res2@tmYLLabelFontHeightF   = .018
 res2@lbAutoManage      =       False
 res2@lbOrientation      =       "Vertical"
 res2@pmLabelBarOrthogonalPosF   =       -0.01
 res2@lbLabelFontHeightF        =       0.012
 res2@cnLineLabelsOn    =       False
; res2@gsnLeftString =   "REG["+varname+",(-)SIE trend index], P1"
;res2@gsnLeftString =   "(a)REG["+varname+"/WAF trend,WP index], P1"
res2@gsnLeftString =   "(b)REG["+varname+"/WAF trend,WP index], RD"
;
 res2@gsnLeftStringFontHeightF  =   20.0

 plot2=gsn_csm_contour_map(wks,reg_ts(:,:),res2)



res    =   True
res@gsnDraw    =   False
res@gsnFrame   =   False
res@vcMaxMagnitudeF	=	0.04
res@vcMinMagnitudeF	=	0.0
res@vcRefMagnitudeF          = 0.03            ; define vector ref mag
res@vcRefLengthF             = 0.025            ; define length of vec ref
res@vfXCStride =   10   
res@vfYCStride =   8    
res@vcMinDistanceF =   0.02
res@vcLineArrowThicknessF  =   1.5
res@vcGlyphStyle             = "LineArrow"    ; turn on curly vectors
res@vcRefAnnoOrthogonalPosF  = 0.2;.2              ; move ref vector down


plotv	=	gsn_csm_vector(wks,reg_fx,reg_fy,res)


txres	=	True
txres@txFontHeightF =   0.015

unit	=	"[mK~S~-1~N~/m~S~2~N~s~S~-2~N~K~S~-1~N~]"
gsn_text_ndc(wks,unit, 0.87, 0.72, txres)


lonL1	=	185
lonR1	=	195
latS1	=	12.5
latN1	=	22.5

lonL2	=	175
lonR2	=	185
latS2	=	40
latN2	=	50

lonL3	=	220
lonR3	=	230
latS3	=	65
latN3	=	75

lon1    =       (/lonL1,lonR1,lonR1,lonL1,lonL1/)
lat1    =       (/latS1,latS1,latN1,latN1,latS1/)

lon2    =       (/lonL2,lonR2,lonR2,lonL2,lonL2/)
lat2    =       (/latS2,latS2,latN2,latN2,latS2/)

lon3    =       (/lonL3,lonR3,lonR3,lonL3,lonL3/)
lat3    =       (/latS3,latS3,latN3,latN3,latS3/)



pres1   =       True
pres1@gsLineThicknessF  =       4.0
pres1@gsLineColor       =       "green" ;R1
pres1@gsLineDashPattern =       0

dum1    =       new(4,graphic)
dum2    =       new(4,graphic)
dum3    =       new(4,graphic)

npts    =       2

do i=0,4-1

gcdist_R1       =       gc_latlon(lat1(i),lon1(i),lat1(i+1),lon1(i+1),npts,2)
gcdist_R2       =       gc_latlon(lat2(i),lon2(i),lat2(i+1),lon2(i+1),npts,2)
gcdist_R3       =       gc_latlon(lat3(i),lon3(i),lat3(i+1),lon3(i+1),npts,2)


dum1(i)         =       gsn_add_polyline(wks,plot2,gcdist_R1@gclon,gcdist_R1@gclat,pres1)
dum2(i)         =       gsn_add_polyline(wks,plot2,gcdist_R2@gclon,gcdist_R2@gclat,pres1)
dum3(i)         =       gsn_add_polyline(wks,plot2,gcdist_R3@gclon,gcdist_R3@gclat,pres1)


delete(gcdist_R1)
delete(gcdist_R2)
delete(gcdist_R3)
end do




overlay(plot2,plot1)
;overlay(plot2,plotv)
draw(plot2)
frame(wks)

end
