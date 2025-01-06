begin

;expname	=	"RU_ONDJFM"
expname =   "RD_ONDJFM_90y"
drmod   =   "png"

varname	=	"RWS"

;re   = 6.37122e06
rad  = 4.0 * atan(1.0) / 180.0 
;wgt = cos(lat2d*rad) 

f_sie	=	addfile("SIE_40to90_JJAS.nc","r")
sie	=	f_sie->var


sie_ru	=	regCoef_n(ispan(1,100,1),sie(:,0:99),0,1)
sie_rd	=	regCoef_n(ispan(1,90,1),sie(:,190:279),0,1)


f_amoc	=	addfile("AMOC_trend_P2_90y.nc","r")
;f_amoc	=	addfile("AMOC_trend_P1.nc","r")
amoc	=	f_amoc->var


;f_EN    =   addfile("WP_index_RU_ONDJFM.nc","r")
f_EN    =   addfile("WP_index_RD_ONDJFM.nc","r")

EN  =   f_EN->var


;f_RWS	=	addfile("RWS_trend_"+expname+".nc","r")
;RWSi	=	f_RWS->var(:,{35:45},{120:140})
;wlat	=	RWSi&lat
;wgt	=	cos(wlat*rad)
;RWS	=	wgt_areaave(RWSi,wgt,1.0,1)

;pc	=	RWS
pc =   EN


;pc	=	sie_ru;*(-1.0)
;pc	=	sie_rd

;pc	=	gmst
;pc	=	amoc


;remove amoc


pc_amoc =   regCoef_n(amoc,pc,0,0)

do i=0,28-1
pc_y    =   pc_amoc*amoc(i)+pc_amoc@yintercept
pc(i)  =   pc(i)  - pc_y
end do


f_ts	=	addfile(varname+"_trend_"+expname+".nc","r")
;tsi =   f_ts->var(:,:,:)
tsi	=	f_ts->term1


lat    =   tsi&lat
lon    =   tsi&lon
nlat   =   dimsizes(lat)
nlon   =   dimsizes(lon)



;Regression

ts_amoc	=	regCoef_n(amoc,tsi,0,0)

do i=0,28-1
ts_y    =   ts_amoc*amoc(i)+onedtond(ts_amoc@yintercept,(/nlat,nlon/))
tsi(i,:,:)  =   tsi(i,:,:)  - ts_y
end do


reg_ts	=	regCoef_n(pc,tsi,0,0)
;reg_ts	=	escorc_n(pc,tsi,0,0)

;rstd	=	onedtond(reg_ts@rstd,dimsizes(reg_ts))

;f_regg	=	addfile("REG_RWS_RD.nc","r")
;reg_rd	=	f_regg->var

;f_r1	=	addfile("RSTD_RWS_RD.nc","r")
;rstdi	=	f_r1->var


reg_ts!0	=	"lat"
reg_ts!1	=	"lon"
reg_ts&lat	=	lat
reg_ts&lon	=	lon


;pr_ts	=	rtest(reg_ts,28,0)
;copy_VarCoords(reg_ts,pr_ts)

tval    =   onedtond(reg_ts@tval,dimsizes(reg_ts))
df  =   onedtond(reg_ts@nptxy,dimsizes(reg_ts))
b   =   tval
b   =   0.5 
prob    =   (1-betainc(df/(df+tval^2),df/2.0,b))

;draw_reg	=	reg_ts - reg_rd
;prob	=	draw_reg/sqrt((rstd)^2+(rstdi)^2)
draw_reg	=	reg_ts


;prob!0  =   "lat"
;prob!1  =   "lon"
;prob&lat    =   lat 
;prob&lon    =   lon 
copy_VarCoords(reg_ts,prob)
copy_VarCoords(reg_ts,draw_reg)


;wks		=	gsn_open_wks(drmod,"Cor_"+varname+"trend_SIE_"+expname+"_rAMOC")
;wks        =   gsn_open_wks(drmod,"REG_"+varname+"trend_WP_"+expname+"_rAMOC")
wks        =   gsn_open_wks(drmod,"REG_"+varname+"term1_WP_"+expname+"_rAMOC")

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
res1@cnLevels   = (/0.95/) ;(/1.64/)
;res1@cnLevels	=	(/0.05/)
res1@cnFillPatterns =   (/-1,3/)
;res1@cnFillPatterns =   (/3,-1/)

plot1   =       gsn_csm_contour(wks,prob(:,:),res1)
;plot1   =       gsn_csm_contour(wks,pr_ts(:,:),res1)



maxval	=2e-10;2.8e+6;5e-08;
minval	=-2e-10;-2.8e+6;-5e-08;
space	=0.4e-10;0.4e+6;1e-08

 res2=True
 res2@gsnDraw	=	False
 res2@gsnFrame	=	False
 res2@cnFillOn=True
 res2@cnLinesOn=False
 res2@cnLineLabelsOn=False
 res2@gsnAddCyclic=True
 res2@mpMinLonF=120
 res2@mpMaxLonF=270
 res2@mpMinLatF=0
 res2@mpMaxLatF=70
; res2@gsnMaximize       =       True
 res2@mpCenterLonF=195
 res2@cnLevelSelectionMode =    "ManualLevels"
;res2@cnMaxLevelValF    =maxval;0.7;500;3.*10^(-8);500;40;2.6; 
;res2@cnMinLevelValF    =minval;-0.7;-500;-3.*10^(-8);-500;-40;-2.6;
;res2@cnLevelSpacingF   =space;0.1;50;0.5*10^(-8);50;5;0.2

res2@cnMaxLevelValF    = maxval;
res2@cnMinLevelValF    = minval
res2@cnLevelSpacingF   =space

 ;res2@tmXBLabelFontHeightF   = .018         ; change font height
 ;res2@tmYLLabelFontHeightF   = .018
 ;res2@lbAutoManage      =       False
 ;res2@lbOrientation      =       "Vertical"
 ;res2@pmLabelBarOrthogonalPosF   =       -0.01
 ;res2@lbLabelFontHeightF        =       0.012
 ;res2@cnLineLabelsOn    =       False


;res2@gsnLeftString =   "(b) REG["+varname+" trend,WP index], RD"
res2@gsnLeftString =   "(b) REG["+varname+"_term1 trend,WP index], RD"

;res2@gsnLeftString =   "REG["+varname+" trend,WP index], RD"
;
; res2@gsnLeftString =   "(a) Cor["+varname+" trend,SIA trend], RU"
;res2@gsnLeftString =   "(b)Cor["+varname+" trend,SIA trend index], RD"
 res2@gsnLeftStringFontHeightF  =   20.0


 plot2=gsn_csm_contour_map(wks,reg_ts(:,:),res2)
;plot2=gsn_csm_contour_map(wks,draw_reg(:,:),res2)
;plot2=gsn_csm_contour_map(wks,cor_ts(:,:),res2)

;Polyline

lonL1   =   120
lonR1   =   150
latS1   =   35
latN1   =   45

lon1    =       (/lonL1,lonR1,lonR1,lonL1,lonL1/)
lat1    =       (/latS1,latS1,latN1,latN1,latS1/)

pres1   =       True
pres1@gsLineThicknessF  =       4.0
pres1@gsLineColor       =       "red" ;R1
pres1@gsLineDashPattern =       0

dum1    =       new(4,graphic)

npts    =       2

do i=0,4-1

gcdist_R1       =       gc_latlon(lat1(i),lon1(i),lat1(i+1),lon1(i+1),npts,2)

;dum1(i)         =       gsn_add_polyline(wks,plot2,gcdist_R1@gclon,gcdist_R1@gclat,pres1)

delete(gcdist_R1)
end do





txres   =   True
txres@txFontHeightF =   0.015


unit    =   "[s~S~-2~N~K~S~-1~N~]"
;gsn_text_ndc(wks,unit, 0.87, 0.72, txres)




overlay(plot2,plot1)
draw(plot2)
frame(wks)

end
