; 3-D wave-activity flux derived by Takaya and Nakamura (1999, 2001)
; See (38) of Takaya and Nakamura (2001, JAS)
;
; Used data: 
;  Monthly-mean data of NCEP/NCAR reanalys 1
;   geopotential height (hgt: m) ; hgt.mon.mean.nc
;
;  Monthly climatology of NCEP/NCAR reanalys 1 
;   geopotential height (hgt: m) ; hgt.mon.mean.nc
;   air tempeature (air: degC)   ; air.mon.mean.nc
;   zonal wind (uwnd: m/s)       ; uwnd.mon.mean.nc
;   meridional wind (vwnd: m/s)  ; vwnd.mon.mean.nc
;
; The data are available at 
; http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;
; The period for calcuation is January, 1990 to December, 1990
; (you can modify this by changing fyear,fmon, lyear, lmon)
;
; The unit of level is [hPa]
;
; Basic state: monthly climatology
;
; Perturbation: deviation from climatology
;
; Calculating the flux for each month of input data.
;
; Note that 
; (i) "p" (pressure) in Eq. (5.7) is dimensionless (pressure/1000mb)


; Output files are
; x-component:  TN2001-Fx.monthly.1990.nc
; y-component:  TN2001-Fy.monthly.1990.nc
; z-component:  TN2001-Fz.monthly.1990.nc
; QG stream function anomaly: psidev.monthly.1990.nc
; Brunt Vaisala frequency: NN.monthly.1990.nc

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

expname	=	"RD_ONDJFM_90y"
ntime=	90

; monthly data
; modify directory name
zfile = addfile("Z200_"+expname+".nc","r") ; regression 

; Climatology data
btfile = addfile("T200_"+expname+".nc","r") ;t, K
;bzfile = addfile("Z500_trend_"+expname+".nc","r") ;z
bufile = addfile("U200_"+expname+".nc","r") ;u
bvfile = addfile("V200_"+expname+".nc","r") ;v

; Test a variable's type prior to reading the variable
zvarorg = zfile->var
btvar = btfile->var
;bzvar = bzfile->var
buvar = bufile->var
bvvar = bvfile->var

;[28,100,lat,lon]

btvar	=	btvar - 273.15 ;Kelvin to degC

lat = zvarorg&lat
lon = zvarorg&lon

nlat = dimsizes(lat)
nlon = dimsizes(lon)



; make arrays for climatology that are suitable for calculation
; bzvar <- bzvarorg
; btvar <- btvarorg



z_zm	=	dim_avg_n_Wrap(zvarorg,3) ;zonal mean = [28,100,lat]
do k=0,28-1
do l=0,ntime-1
do i=0,nlat-1
zvarorg(k,l,i,:)	=	zvarorg(k,l,i,:) - z_zm(k,l,i)
end do
end do
end do

czvar	=	new((/28,ntime,nlat,nlon/),float)
ctvar	=	new((/28,ntime,nlat,nlon/),float)
cuvar	=	new((/28,ntime,nlat,nlon/),float)
cvvar	=	new((/28,ntime,nlat,nlon/),float)


;height anomaly
zavar = zvarorg

do i=0,28-1
do j=0,ntime-1
czvar(i,j,:,:)	=	dim_avg_n_Wrap(zvarorg(i,:,:,:),0)
ctvar(i,j,:,:)	=	dim_avg_n_Wrap(btvar(i,:,:,:),0)
cuvar(i,j,:,:)    =   dim_avg_n_Wrap(buvar(i,:,:,:),0)
cvvar(i,j,:,:)    =   dim_avg_n_Wrap(bvvar(i,:,:,:),0)
end do
end do

;delete(bzvar)
delete(btvar)
delete(buvar)
delete(bvvar)



;  Gas constant
gc=290
;  Gravitational acceleration
ga=9.80665

;  Radius of the earth
re=6378388

; scale height
sclhgt=8000.

; pi
pi = atan(1.0)*4.

; Coriolis parameter
f =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
f!0 = "lat"
f&lat = lat
f@_FillValue = zvarorg@_FillValue

zavar	=	zavar/ga
czvar	=	czvar/ga


; missing for 10S - 10N
do ilat = 0, nlat-1
if (abs(lat(ilat) ).lt. 10. ) then
f(ilat)= f@_FillValue
end if
end do

; cosine
coslat = cos(lat(:)*pi/180.)

; 1-D -> 4-D
leveltmp = 197.908

coslattmp = conform_dims(dimsizes(zavar),coslat,2)
ftmp = conform_dims(dimsizes(zavar),f,2)


; magnitude of climatological wind
cumag = sqrt(cuvar^2 + cvvar^2)
cumag@_FillValue = cuvar@_FillValue
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function for anomaly
psidev = zavar*ga /ftmp ;[28,100,lat,lon]

;dpsidev/dlon
dpsidevdlon =  center_finite_diff_n(psidev,lon*pi/180.,True,0,3)

;ddpsidev/dlonlon
ddpsidevdlonlon =  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,3)

;dpsidev/dlat
dpsidevdlat = center_finite_diff_n(psidev, lat*pi/180., False,0,2)

;ddpsidev/dlonlat
ddpsidevdlonlat =  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,2)

;ddpsidev/dlatdlat
ddpsidevdlatlat = center_finite_diff_n(dpsidevdlat, lat*pi/180.,False,0,2)



xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
;yuterm = xvterm
yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)


; Mask out where westerlies is small or negative (less than 5 m/s).
;  by using mask

;x-component of (38)
Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( cuvar/coslattmp * xuterm + cvvar * xvterm), cuvar.lt.5.0,False)

;y-component 
Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( cuvar*xvterm + coslattmp*cvvar*yvterm),  cuvar.lt.5.0,False)


; for output
;Fx(0,:)	=	Fx@_FillValue
;Fx(nlat-1,:)	=	Fx@_FillValue

;Fy(0,:)	=	Fy@_FillValue
;Fy(nlat-1,:)	=	Fy@_FillValue

Fx!2 = "lat"
Fx&lat = lat
Fx!3 = "lon"
Fx&lon = lon

Fy!2 = "lat"
Fy&lat = lat
Fy!3 = "lon"
Fy&lon = lon

psidev!2 = "lat"
psidev&lat = lat
psidev!3 = "lon"
psidev&lon = lon


Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"

psidev@units = "m^2/s"

ncFx =     addfile("Fx_monthly_"+expname+".nc","c")
ncFy =     addfile("Fy_monthly_"+expname+".nc","c")
;ncpsidev = addfile("psidev_monthly_"+expname+".nc","c")

ncFx->Fx = Fx
ncFy->Fy = Fy
;ncpsidev->psidev = psidev



end
