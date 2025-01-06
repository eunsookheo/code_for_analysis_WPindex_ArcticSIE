begin

;Draw SEP SIE (ensemble mean)

years   =   ispan(2001,2399,1)
nyear   =   399 
;years   =   ispan(2001,2280,1)
;nyear   =   280



years_ru	=	ispan(2001,2100,1)
years_rd	=	ispan(2191,2280,1)

;obs SIE
f_obs	=	addfile("/home/esheo/SRC/SIE_obs/SIE_obs.nc","r")
sie_obs	=	f_obs->sie ;1979-2021

obs_JJAS	=	(sie_obs(22:42,5)+sie_obs(22:42,6)+sie_obs(22:42,7)+sie_obs(22:42,8))/4.0

trend_obs	=	regCoef_n(ispan(2001,2021,1),obs_JJAS,0,0)
obs	=	trend_obs*ispan(2001,2021,1) + trend_obs@yintercept



;f_in	=	addfile("/home/esheo/SRC/SIE_model/SIE_40to90_JJAS.nc","r")
f_in	=	addfile("/home/esheo/SRC/SIE_model/SIE_40to90_ANN.nc","r")
;f_in	=	addfile("/home/esheo/SRC/SIE_model/SIE_40to90_ONDJFM.nc","r")

;f_in	=	addfile("SIV_JJAS_all.nc","r")
var		=	f_in->var

f_sep	=	addfile("/home/esheo/PhD_research/SIV_JASON_all.nc","r")
;var_sep	=	f_sep->sie
;var	=	f_sep->sie

trend_ru	=	regCoef_n(years_ru,var(:,0:99),0,1)
trend_rd	=	regCoef_n(years_rd,var(:,190:279),0,1)

;trend_sep	=	regCoef_n(years_ru,var_sep(:,0:99),0,1)
;cor	=	escorc(trend_ru,trend_sep)
;pr	=	rtest(cor,28,0)
;print(cor+"		"+pr)
;exit

ind_high_ru	=	maxind(trend_ru)
ind_low_ru	=	minind(trend_ru)

ind_high_rd	=	maxind(trend_rd)
ind_low_rd	=	minind(trend_rd)


high_ru	=	trend_ru(ind_high_ru)*years_ru + trend_ru@yintercept(ind_high_ru)
low_ru =   trend_ru(ind_low_ru)*years_ru + trend_ru@yintercept(ind_low_ru)

high_rd =   trend_rd(ind_high_rd)*years_rd + trend_rd@yintercept(ind_high_rd)
low_rd =   trend_rd(ind_low_rd)*years_rd + trend_rd@yintercept(ind_low_rd)


var_em  =   dim_avg_n(var,0) 

;opt =   "JJAS"
opt	=	"ANN"

var_em!0    =   "year"
var_em&year =   years

mnmx    =   new((/2,nyear/),float)
mnmx(0,:)   =  dim_min_n(var,0); var_em(:) - var_1std(:)
mnmx(1,:)   =  dim_max_n(var,0); var_em(:) + var_1std(:)

mnmx!1  =   "year"
mnmx&year   =   years

wks	 = gsn_open_wks("png","SIE_"+opt+"_40to90")
;wks     = gsn_open_wks("x11","FYarea_"+opt+"_40to90")


 res                    = True              ; plot mods desired
  res@gsnDraw            = False             ; don't draw yet
  res@gsnFrame           = False             ; don't advance frame yet

  res@vpHeightF      = 0.4               ; change aspect ratio of plot
  res@vpWidthF           = 0.7
  res@gsnYRefLine   =   0.
  res@gsnYRefLineColor  =   "grey"
  res@gsnYRefLineDashPattern = 0
 ; res@trYMaxF            = 10
  res@trYMinF   =   0
  res@trXMinF            = 2000
  res@trXMaxF   =   2281            ; set x-axis minimum
  res@gsnLeftString =   "(a)"+opt+" SIA (40N-90N), CESM1.2"
  res@gsnLeftStringFontHeightF  =   25
  res@xyLineColor       = "black" ; colors chosen
  res@xyLineThicknessF   = 3.      ; line thicknesses
  res@xyDashPattern  = 0.      ; make all lines solid

  res@tiXAxisString  = "Year"      ; add an axis title    
  res@tiYAxisString =   "SIA*10~S~6~N~ [km~S~2~N~]"

  plot1 =   gsn_csm_xy (wks,var_em&year,var_em,res)

  delete(res@xyLineColor)
  res@gsnXYFillColors = "Lightgrey"
  res@gsnXYFillOpacities = 0.5
  res@xyLineColor     = -1
  res@gsnXRefLine   =   (/2100,2191/)
  res@gsnXRefLineThicknessF =   2.
  res@gsnXRefLineColor  =   "grey"
  res@gsnXRefLineDashPatterns   =   5

  res@gsnYRefLine   =   (/1.0/)

  res@gsnYRefLineThicknessF =   2.
  res@gsnYRefLineColor  =   "grey"
  res@gsnYRefLineDashPattern   =   5
  plot2  = gsn_csm_xy (wks,mnmx&year,mnmx(:,:),res)  ; Create filled XY plot.

pmres	=	True
pmres@gsLineThicknessF	=	2.
pmres@gsLineDashPattern	=	5

dum	=	new((/4/),graphic)

pmres@gsLineColor	=	"blue"
dum(0)	=	gsn_add_polyline(wks,plot2,years_ru,high_ru,pmres)
pmres@gsLineColor   =   "red"
dum(1)  =   gsn_add_polyline(wks,plot2,years_ru,low_ru,pmres)

pmres@gsLineColor   =   "blue"
dum(2)	=	gsn_add_polyline(wks,plot2,years_rd,high_rd,pmres)
pmres@gsLineColor   =   "red"
dum(3)  =   gsn_add_polyline(wks,plot2,years_rd,low_rd,pmres)


overlay(plot2,plot1)
draw(plot2)
frame(wks)



end
