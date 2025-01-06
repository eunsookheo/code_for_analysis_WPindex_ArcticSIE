begin

expname1 =   "RU_ONDJFM"
expname2	=	"RD_ONDJFM_90y"
drmod   =   "x11"
varname	=	"FLNT" 

f_sie	=	addfile("SIE_40to90_JJAS.nc","r")
sie	=	f_sie->var

f_amoc1	=	addfile("AMOC_trend_P1.nc","r")
f_amoc2  =   addfile("AMOC_trend_P2_90y.nc","r")
amoc1    =   f_amoc1->var
amoc2	=	f_amoc2->var

sie_ru	=	regCoef_n(ispan(1,100,1),sie(:,0:99),0,1)
sie_rd	=	regCoef_n(ispan(1,100,1),sie(:,180:279),0,1)
;sie_ru	=	 regCoef_n(ispan(1,40,1),sie(:,100:139),0,1)
;sie_rd  =   regCoef_n(ispan(1,40,1),sie(:,140:179),0,1)


;pc	=	sie_ru*(-1.0)
;pc	=	sie_rd




f_EN1	=	addfile("WP_index_RU_ONDJFM.nc","r")
f_EN2	=	addfile("WP_index_RD_ONDJFM.nc","r")
EN1	=	f_EN1->var
EN2	=	f_EN2->var

pc1		=	EN1
pc2		=	EN2
;pc1	=	sie_ru*(-1.0)
;pc2	=	sie_rd*(-1.0)




;remove amoc

pc1_amoc =   regCoef_n(amoc1,pc1,0,0)
pc2_amoc	=	regCoef_n(amoc2,pc2,0,0)

do i=0,28-1
pc1_y    =   pc1_amoc*amoc1(i)+pc1_amoc@yintercept
pc2_y    =   pc2_amoc*amoc2(i)+pc2_amoc@yintercept

pc1(i)  =   pc1(i)  - pc1_y
pc2(i)  =   pc2(i)  - pc2_y

end do



;============================

f_ts1	=	addfile(varname+"_trend_"+expname1+".nc","r")
tsi1 =   f_ts1->var(:,:,:)

f_ts2   =   addfile(varname+"_trend_"+expname2+".nc","r")
tsi2 =   f_ts2->var(:,:,:)

;RD -> OLR index

olr_i	=	wgt_areaave_Wrap(tsi2(:,{10:20},{120:140}),1.0,1.0,0)
f_out	=	addfile("OLR_index2_RD.nc","c")
f_out->var	=	olr_i
exit


lat    =   tsi1&lat
lon    =   tsi1&lon
nlat   =   dimsizes(lat)
nlon   =   dimsizes(lon)

ts1_amoc	=	regCoef_n(amoc1,tsi1,0,0)
ts2_amoc =   regCoef_n(amoc2,tsi2,0,0)

do i=0,28-1

ts1_y    =   ts1_amoc*amoc1(i)+onedtond(ts1_amoc@yintercept,(/nlat,nlon/))
ts2_y    =   ts2_amoc*amoc2(i)+onedtond(ts2_amoc@yintercept,(/nlat,nlon/))
tsi1(i,:,:)  =   tsi1(i,:,:)  - ts1_y
tsi2(i,:,:)  =   tsi2(i,:,:)  - ts2_y

end do


reg_ts1	=	regCoef_n(pc1,tsi1,0,0)
;reg_ts2	=	regCoef_n(pc2,tsi2,0,0)
reg_ts2	=	reg_ts1

rstd_1 =   onedtond(reg_ts1@rstd,dimsizes(reg_ts1))
rstd_2 =   onedtond(reg_ts2@rstd,dimsizes(reg_ts2))

;diff	=	reg_ts2 - reg_ts1
;prob    =   abs(diff)/sqrt((rstd_1)^2+(rstd_2)^2)

diff	=	reg_ts2
tval    =   onedtond(reg_ts2@tval,dimsizes(reg_ts2))
df  =   onedtond(reg_ts2@nptxy,dimsizes(reg_ts2))
b   =   tval
b   =   0.5
prob    =   (1-betainc(df/(df+tval^2),df/2.0,b))



reg_ts1!0	=	"lat"
reg_ts1!1	=	"lon"
reg_ts1&lat	=	lat
reg_ts1&lon	=	lon

copy_VarCoords(reg_ts1,reg_ts2)
copy_VarCoords(reg_ts1,prob)
copy_VarCoords(reg_ts1,diff)


wks        =   gsn_open_wks(drmod,"REGdiff_"+varname+"trend_WP")

res     =       True
res@cnFillOn    =       True
res@cnLinesOn   =       False
res@gsnDraw     =       False
res@gsnFrame    =       False
res@cnInfoLabelOn       =       False
res@lbLabelBarOn   =   False
res@cnLineLabelsOn      =       False
res@cnMonoFillColor = True
res@cnMonoFillPattern  =   False
res@cnLevelSelectionMode= "ExplicitLevels"
;res@cnLevels   =   (/1.64/)
res@cnLevels	=	(/0.95/)
res@cnFillPatterns =   (/-1,6/)

plot	=	gsn_csm_contour(wks,prob,res)


res1     =       True
res1@cnFillOn    =       False
res1@cnLinesOn   =       True
res1@gsnDraw     =       False
res1@gsnFrame    =       False
res1@cnInfoLabelOn       =       False
res1@lbLabelBarOn   =   False
res1@cnLineLabelsOn      =       False
res1@cnLevelSelectionMode= "ExplicitLevels"
res1@gsnContourNegLineDashPattern = 1
;res1@cnLevels	=	(/-3.*10^(-8),-2.*10^(-8),-1.*10^(-8),1.*10^(-8),2.*10^(-8),3.*10^(-8)/)
;res1@cnLevels  =   (/-30,-20,-10,10,20,30/)
res1@cnLevels  =   (/-300,-250,-200,-150,-100,-50,50,100,150,200,250,300/)
;res1@cnLevels	=	(/-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7/)

res1@cnLineThicknessF	=	2.0


;plot1   =       gsn_csm_contour(wks,reg_ts1,res1)



 res2=True
 res2@gsnDraw	=	False
 res2@gsnFrame	=	False
 res2@cnFillOn=True
 res2@cnLinesOn=False
 res2@cnLineLabelsOn=False

res2@mpFillOn    =   False
res2@cnFillOn   =       True
res2@cnLinesOn   =       False
res2@gsnPolarLabelFontHeightF   =       20.0
res2@mpGeophysicalLineThicknessF =       1.5
res2@mpGreatCircleLinesOn       =       False
 res2@gsnAddCyclic=False
 res2@mpMinLonF=100
 res2@mpMaxLonF=270
 res2@mpMinLatF=-40
 res2@mpMaxLatF=40
 res2@mpCenterLonF=370/2.
 res2@cnLevelSelectionMode =    "ManualLevels"
res2@cnMaxLevelValF    =12;1e-9;300;30;3.*10^(-8);500;2.6; 
 res2@cnMinLevelValF    =-12;-1e-9;-300;-30;-3.*10^(-8);-500;-2.6;
 res2@cnLevelSpacingF   =2;1e-10;50;10;1.*10^(-8);50;0.2
 res2@tmXBLabelFontHeightF   = .018         ; change font height
 res2@tmYLLabelFontHeightF   = .018
 res2@lbAutoManage      =       False
 res2@lbOrientation      =       "Vertical"
 res2@pmLabelBarOrthogonalPosF   =       -0.01
 res2@lbLabelFontHeightF        =       0.012

res2@gsnLeftString =   "(a) REG[OLR trend,WP index], RD-RU"
 res2@gsnLeftStringFontHeightF  =   20.0

; plot2=gsn_csm_contour_map(wks,reg_ts2(:,:),res2)
plot2=gsn_csm_contour_map(wks,diff(:,:),res2)



overlay(plot2,plot)
;overlay(plot2,plot1)
draw(plot2)
frame(wks)

end
