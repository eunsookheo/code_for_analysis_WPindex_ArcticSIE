begin

;draw Local HC for 120E-150E


expname	=	"RU"


f_in	=	addfile("LocalHC_trend_"+expname+"_ONDJFM.nc","r")
hci	=	f_in->var

lev	=	hci&lev
lat	=	hci&lat


nlev	=	dimsizes(lev)
nlat	=	dimsizes(lat)


f_pc   =   addfile("WP_index_RU_ONDJFM.nc","r")
;f_pc   =   addfile("WP_index_RD_ONDJFM.nc","r")
pc =   f_pc->var

f_amoc =   addfile("AMOC_trend_P1.nc","r")
;f_amoc  =   addfile("AMOC_trend_P2_90y.nc","r")
amoc    =   f_amoc->var


;remove amoc

pc_amoc =   regCoef_n(amoc,pc,0,0)

do i=0,28-1
pc_y    =   pc_amoc*amoc(i)+pc_amoc@yintercept
pc(i)  =   pc(i)  - pc_y
end do

hc_amoc =    regCoef_n(amoc,hci,0,0)

do i=0,28-1

hc_y    =   hc_amoc*amoc(i)+onedtond(hc_amoc@yintercept,(/nlev,nlat/))
hci(i,:,:)  =   hci(i,:,:)  - hc_y


end do

;hc = [model,lev,lat]

reg_hc	=	regCoef_n(pc,hci,0,0)
reg_hc!0	=	"lev"
reg_hc!1	=	"lat"
reg_hc&lev	=	lev
reg_hc&lat	=	lat

tval    =   onedtond(reg_hc@tval,dimsizes(reg_hc))
df  =   onedtond(reg_hc@nptxy,dimsizes(reg_hc))
b   =   tval
b   =   0.5
prob    =   (1-betainc(df/(df+tval^2),df/2.0,b))
copy_VarCoords(reg_hc,prob)


wks =   gsn_open_wks("png","LocalHC_ONDJFM_"+expname)

res1	=	True
res1@cnFillOn    =       True
res1@cnLinesOn   =       False
res1@gsnDraw     =       False
res1@gsnFrame    =       False
res1@cnFillColor    =   "grey"
res1@cnInfoLabelOn       =       False
res1@cnLineLabelsOn     =   False
res1@lbLabelBarOn   =   False
res1@cnLineLabelsOn      =       False
res1@cnMonoFillColor = True
res1@cnMonoFillPattern  =   False
res1@cnLevelSelectionMode= "ExplicitLevels"
res1@cnLevels   =   (/0.95/)
res1@cnFillPatterns =   (/-1,12/)

plot1  =   gsn_csm_pres_hgt(wks,prob({100000:5000},:),res1)


  res2  =   True
  res2@cnLinesOn    =   False
  res2@cnFillOn =   True
  res2@tiMainString         = "(c) REG[local HC trend,WP index], "+expname        ; title
  res2@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
  res2@cnMinLevelValF       = -2e+11                  ; min level
  res2@cnMaxLevelValF       =  2e+11               ; max level
  res2@cnLineLabelsOn       = False                  ; turn on line labels
  res2@cnLevelSpacingF  =   0.4e+11
  res2@cnFillPalette        = "BlWhRe"              ; set color map
  res2@gsnDraw  =   False
  res2@gsnFrame =   False
  plot = gsn_csm_pres_hgt(wks,reg_hc({100000:5000},:),res2)


txres   =   True
txres@txFontHeightF =   0.015


unit    =   "[kgs~S~-1~N~K~S~-1~N~]"
gsn_text_ndc(wks,unit, 0.85, 0.12, txres)


overlay(plot,plot1)
draw(plot)
frame(wks)
print(lev)

end
