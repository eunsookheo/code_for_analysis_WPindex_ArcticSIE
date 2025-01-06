begin

;Calculate Rossby wave source for 100 year time series
;for 200hpa

expname	=	"RU_ONDJFM"
nyear	=	100

;U200/V200

f_u	=	addfile("U200_"+expname+".nc","r")
f_v	=	addfile("V200_"+expname+".nc","r")

ui	=	f_u->var
vi	=	f_v->var

lat	=	ui&lat
lon	=	ui&lon
nlat	=	dimsizes(lat)
nlon	=	dimsizes(lon)

;divergent wind/relative vorticity
;finite volume grid

dv	=	uv2dvF_Wrap(ui,vi) ;divergence
printVarSummary(dv)

ud    = new ( dimsizes(dv), typeof(dv))
vd    = new ( dimsizes(dv), typeof(dv))

dv2uvf(dv,ud,vd) ;divergent wind

vr	=	uv2vrF_Wrap(ui,vi) ;relative vorticity

f	=	coriolis_param(lat)
fc	=	conform_dims(dimsizes(dv),f,2)


avr	=	vr+fc ;absolute vorticity

print("!")
;velocity potential
;sf   = new ( dimsizes(dv), typeof(dv))
;vp    = new ( dimsizes(dv), typeof(dv))

;uv2sfvpf(ui,vi,sf,vp)


;calculate RWS

gradLatLon_vr	=	grad_latlon_cfd(avr,lat,lon,True,False) ;[model,year,lat,lon]
gradLatLon_ud   =   grad_latlon_cfd(ud,lat,lon,True,False)
gradLatLon_vd   =   grad_latlon_cfd(vd,lat,lon,True,False)

dvrdx	=	gradLatLon_vr[1]
dvrdy	=	gradLatLon_vr[0]

duddx	=	gradLatLon_ud[1]
dvddy	=	gradLatLon_vd[0]

RWS_term1	=	-ud*dvrdx-vd*dvrdy
RWS_term2	=	-avr*duddx-avr*dvddy

RWS	=	RWS_term1+RWS_term2
wind_dv	=	duddx+dvddy


;calculate trend
;[28,lat,lon]
trend_term1	=	regCoef_n(ispan(1,nyear,1),RWS_term1,0,1)
trend_term2	=	regCoef_n(ispan(1,nyear,1),RWS_term2,0,1)
trend_RWS	=	regCoef_n(ispan(1,nyear,1),RWS,0,1)
;trend_avr	=	regCoef_n(ispan(1,nyear,1),avr,0,1)
;trend_wind	=	regCoef_n(ispan(1,nyear,1),wind_dv,0,1)

;trend_vp	=	regCoef_n(ispan(1,nyear,1),vp,0,1)

trend_term2(:,0,:)	=	trend_term2@_FillValue
trend_term2(:,nlat-1,:)	=	trend_term2@_FillValue
trend_RWS(:,0,:)	=	trend_term2@_FillValue
trend_RWS(:,nlat-1,:)    =   trend_term2@_FillValue




trend_RWS!1	=	"lat"
trend_RWS!2	=	"lon"
trend_RWS&lat	=	lat
trend_RWS&lon	=	lon

copy_VarCoords(trend_RWS,trend_term1)
copy_VarCoords(trend_RWS,trend_term2)
;copy_VarCoords(trend_RWS,trend_vp)
;copy_VarCoords(trend_RWS,trend_avr)
;copy_VarCoords(trend_RWS,trend_wind)

f_outt	=	addfile("VP200_trend_"+expname+".nc","c")
f_outt->var	=	trend_vp
exit
;f_outt2  =   addfile("winddiv_trend_"+expname+".nc","c")
;f_outt2->var =   trend_wind

;exit

f_out	=	addfile("RWS_trend_"+expname+".nc","c")
f_out->term1	=	trend_term1
f_out->term2	=	trend_term2
f_out->var	=	trend_RWS




end
