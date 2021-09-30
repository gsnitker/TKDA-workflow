#Measure Identified tar kilns
rm(list = ls())

library(rgdal)
library(imager)
library(rgeos)
library(sp)
library(raster)
library(tools)
library(dplyr)
library(tidyr)
library(rgdal)
library(sf)
library(spdep)
library(tidyverse)
library(reshape2)
library(geosphere)
library(ggpubr)
setwd("/Users/grantsnitker/Dropbox/Smoke/Colonial_Urban_Economy_Project/GIS/Tar_kiln_outputs/")
# #85 196 252 259 265 270
# kiln = known.kilns[270,]
# kiln$Site_num
# radius = 10
# quantile.thresh = .50
# plot = T
# Load function needed to find the trench surrounding the kiln and thus the kilns radius.
#find.trench(kiln = known.kilns[300,])
find.trench = function(kiln, radius = 10, quantile.thresh = .5, alt.method = T, plot = F){
  kiln.mound.dia.measured = kiln$TK_mnd_dia
  kiln.trench.measured = (kiln$TK_tot_dia - kiln$TK_mnd_dia) / 2
  diff = 0
  while (diff <1.5 && radius < 40 ){
    radius = radius + 5
    buffer.trench = buffer(kiln, width = radius)
    buffer.trench.crop = raster::crop(DEM, buffer.trench)
    #contour = rasterToContour(buffer.trench.crop)
    elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh,1))
    diff = as.numeric(elev.stats[2] - elev.stats[1])
    #ring = Which(buffer.trench.crop <= elev.stats[1])
    #plot(ring)
    #contour.ring = rasterToContour(ring, nlevels =1)
    #plot(buffer.trench.crop)
    #plot(contour.ring, add = T)
    #print(radius)
    #print(diff)
    }
  buffer.trench = buffer(kiln, width = radius + (radius * .35))
  buffer.trench.crop = raster::crop(DEM, buffer.trench)
  buffer.trench.crop = raster::mask(buffer.trench.crop, buffer.trench)
  
 #plot(buffer.trench.crop)
  if  (alt.method == T) {
    #contour = rasterToContour(buffer.trench.crop,nlevels =10)
    # elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh+.1,1))
    # diff = as.numeric(elev.stats[2] - elev.stats[1])
    # mound = Which(buffer.trench.crop <= elev.stats[1])
    # plot(ring)
    # contour.mound = rasterToContour(ring, nlevels =1)
    # plot(contour.mound, add = T)
    elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh+.1,1))
    diff = as.numeric(elev.stats[2] - elev.stats[1])
    trench = Which(buffer.trench.crop <= elev.stats[1])
    #plot(trench)
    contour.trench = rasterToContour(trench, nlevels =1)
    #plot(contour.trench, add = T)
    
    
    
  } else { 
    im = as.cimg(buffer.trench.crop)
    tryCatch(expr = {canny <- cannyEdges(im, alpha=1)
             px.im = as.cimg(canny)
             rast =raster(as.matrix(canny), template = buffer.trench.crop)
             contour.trench = rasterToContour(rast)}, error = function(e) { 
               elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh-.2,1))
               diff = as.numeric(elev.stats[2] - elev.stats[1])
               trench = Which(buffer.trench.crop <= elev.stats[1])
               plot(trench)
               contour.trench = rasterToContour(trench, nlevels =1)})
    }
  
  points.1 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(mean(buffer.trench@bbox[2,]),mean(buffer.trench@bbox[2,]))))
  intersect.line.1 <- as(points.1,"SpatialLines")
  intersection.points.1 = gIntersection(intersect.line.1, contour.trench, byid=T)
  if (is.null(intersection.points.1) == T) {
    line.radius.1  = NA
    trench.width.1  = NA} else {
  below.points = intersection.points.1[intersection.points.1@coords[,1] < kiln@coords[,1]]
  above.points = intersection.points.1[intersection.points.1@coords[,1] > kiln@coords[,1]]
  # find the distances between inner and outer portions of trench and center point
  below.dist = pointDistance(kiln, below.points, lonlat = F)
  below.dist = below.dist[below.dist > 8]
  below.dist = below.dist[which(below.dist < Inf)]
  below.radius = min(below.dist)
  above.dist = pointDistance(kiln, above.points, lonlat = F)
  above.dist = above.dist[above.dist > 8]
  above.dist = above.dist[which(above.dist < Inf)]
  above.radius = min(above.dist)
  # calculate kiln radius
  line.radius.1 = (below.radius + above.radius) / 2
  trench.width.above = max(above.dist) - min(above.dist)
  trench.width.below = max(below.dist) - min(below.dist)
  trench.width.1 = (trench.width.below + trench.width.above) / 2
  if (line.radius.1 == Inf) {line.radius.1 = NA}
  if (trench.width.1 == Inf) {trench.width.1 = NA}
  if (trench.width.1 == -Inf) {trench.width.1 = NA}}
  
  
  points.2 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(buffer.trench@bbox[2,1],buffer.trench@bbox[2,2])))
  intersect.line.2 <- as(points.2,"SpatialLines")
  intersection.points.2 = gIntersection(intersect.line.2, contour.trench, byid=T)
  if (is.null(intersection.points.2) == T) {
    line.radius.2  = NA
    trench.width.2  = NA} else {
    below.points = intersection.points.2[intersection.points.2@coords[,1] < kiln@coords[,1]]
    above.points = intersection.points.2[intersection.points.2@coords[,1] > kiln@coords[,1]]
    # find the distances between inner and outer portions of trench and center point
    below.dist = pointDistance(kiln, below.points, lonlat = F)
    below.dist = below.dist[below.dist > 8]
    below.dist = below.dist[which(below.dist < Inf)]
    below.radius = min(below.dist)
    above.dist = pointDistance(kiln, above.points, lonlat = F)
    above.dist = above.dist[above.dist > 8]
    above.dist = above.dist[which(above.dist < Inf)]
    above.radius = min(above.dist)
    # calculate kiln radius
    line.radius.2 = (below.radius + above.radius) / 2
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.2 = (trench.width.below + trench.width.above) / 2
    if (line.radius.2 == Inf) {line.radius.2 = NA}
    if (trench.width.2 == Inf) {trench.width.2 = NA}
    if (trench.width.2 == -Inf) {trench.width.2 = NA}}
  
  
  points.3 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(buffer.trench@bbox[2,2],buffer.trench@bbox[2,1])))
  intersect.line.3 <- as(points.3,"SpatialLines")
  intersection.points.3 = gIntersection(intersect.line.3, contour.trench, byid=T)
  if (is.null(intersection.points.3) == T) {
    line.radius.3  = NA
    trench.width.3 = NA } else {
    below.points = intersection.points.3[intersection.points.3@coords[,1] < kiln@coords[,1]]
    above.points = intersection.points.3[intersection.points.3@coords[,1] > kiln@coords[,1]]
    # find the distances between inner and outer portions of trench and center point
    below.dist = pointDistance(kiln, below.points, lonlat = F)
    below.dist = below.dist[below.dist > 8]
    below.dist = below.dist[which(below.dist < Inf)]
    below.radius = min(below.dist)
    above.dist = pointDistance(kiln, above.points, lonlat = F)
    above.dist = above.dist[above.dist > 8]
    above.dist = above.dist[which(above.dist < Inf)]
    above.radius = min(above.dist)
    # calculate kiln radius
    line.radius.3 = (below.radius + above.radius) / 2
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.3 = (trench.width.below + trench.width.above) / 2
    if (line.radius.3 == Inf) {line.radius.3 = NA}
    if (trench.width.3 == Inf) {trench.width.3 = NA}
    if (trench.width.3 == -Inf) {trench.width.3 = NA}}
  
  
  points.4 = SpatialPoints(cbind(c(mean(buffer.trench@bbox[1,]),mean(buffer.trench@bbox[1,])), c(buffer.trench@bbox[2,1],buffer.trench@bbox[2,2])))
  intersect.line.4 <- as(points.4,"SpatialLines")
  intersection.points.4 = gIntersection(intersect.line.4, contour.trench, byid=T)
  if (is.null(intersection.points.4) == T) {
    line.radius.4  = NA
    trench.width.4  = NA
    } else {
    below.points = intersection.points.4[intersection.points.4@coords[,2] < kiln@coords[,2]]
    above.points = intersection.points.4[intersection.points.4@coords[,2] > kiln@coords[,2]]
    #plot(below.points, add =T, col = "blue")
    #plot(above.points, add =T)
    # find the distances between inner and outer portions of trench and center point
    below.dist = pointDistance(kiln, below.points, lonlat = F)
    below.dist = below.dist[below.dist > 8]
    below.dist = below.dist[which(below.dist < Inf)]
    below.radius = min(below.dist)
    above.dist = pointDistance(kiln, above.points, lonlat = F)
    above.dist = above.dist[above.dist > 8]
    above.dist = above.dist[which(above.dist < Inf)]
    above.radius = min(above.dist)
    # calculate kiln radius
    line.radius.4 = (below.radius + above.radius) / 2
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.4 = (trench.width.below + trench.width.above) / 2
    if (line.radius.4 == Inf) {line.radius.4 = NA}
    if (trench.width.4 == Inf) {trench.width.4 = NA}
    if (trench.width.4 == -Inf) {trench.width.4 = NA}}

  if (plot == T){
    plot(buffer.trench.crop, main = paste0("Kiln ",kiln@data$cat))
    plot(contour.trench, add = T)
    plot(intersect.line.1, add =T)
    plot(intersection.points.1, add =T)
    plot(intersect.line.2, add =T)
    plot(intersection.points.2, add =T)
    plot(intersect.line.3, add =T)
    plot(intersection.points.3, add =T)
    plot(intersect.line.4, add =T)
    plot(intersection.points.4, add =T)
    
  }
  
  radii  = c(line.radius.1, line.radius.2, line.radius.3, line.radius.4)
  radii
  radii = na.omit(radii)
  cutoff <- quantile(radii, probs=c(.05, .95), na.rm = T)
  if (length(radii) > 2) {radii = radii[radii >= cutoff[1] & radii <= cutoff[2]]}
  radii

  trench  = c(trench.width.1, trench.width.2, trench.width.3, trench.width.4)
  trench
  trench = na.omit(trench)
  cutoff <- quantile(trench, probs=c(.1, .9), na.rm = T)
  if (length(trench) > 2) {trench = trench[trench >= cutoff[1] & trench <= cutoff[2]]}
  #trench = trench[trench >0 & trench < 20 ]
  
  trench = mean(c(trench), na.rm=T)
  line.radius = mean(c(radii), na.rm=T)
  kiln.mound.dia.calc = (line.radius * 2) / 3.281
  kiln.trench.calc = trench / 3.281
  kiln.mound.dia.calc
  kiln.mound.dia.measured
  return(list(kiln.mound.dia.calc, kiln.mound.dia.measured, kiln.trench.calc, kiln.trench.measured))
}


# Measure Tar kiln features
# after output has been cleaned in QGIS, re-import to measure
# Import GIS data
known.kilns =  shapefile("/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/shapefiles/FMNF_TKP_Validation_final_NAD83_2011_ft.shp")
FMNF = shapefile("./FMNF_boundary.shp")
DEM = raster("./FMNF_1m_Lidar_DEM.tif")

# plot(FMNF)
# plot(DEM, add =T)
# plot(known.kilns, add = T)


# measure kilns in loop
ptm <- proc.time()
options(warn=-1)
Tar.kilns.measured = data.frame("Index" =NA,"Site" = NA,"kiln.mound.dia.calc"= NA, "kiln.mound.dia.measured"= NA, "kiln.trench.calc"= NA, "kiln.trench.measured"= NA, "Kiln_Vol"=NA, "Kiln_cords_wood"=NA, "X_coord" = NA, "Y_coord" = NA )
pb = invisible(txtProgressBar(min = 0, max = length(known.kilns), initial = 0, style =3)) # create progress bar
#pb = invisible(txtProgressBar(min = 0, max = 100, initial = 0, style =3)) # create progress bar
for (i in 1:length(known.kilns)){ 
  # subset by kiln in chunk 
  #################
  kiln.metrics = find.trench(kiln = known.kilns[i,], quantile.thresh = .40, plot = F)
  kiln.mound.dia.calc = kiln.metrics[[1]]
  kiln.mound.dia.measured = kiln.metrics[[2]]

  kiln.trench.calc = kiln.metrics[[3]]
  kiln.trench.measured = kiln.metrics[[4]]
  kiln.total.dia.calc = kiln.mound.dia.calc + (kiln.trench.calc*2)
  kiln.total.dia.measured = known.kilns[i,]$TK_tot_dia
  
  
  # find volume of cone based on idealized tar kiln form identified in the literature
  #V=13πr2h
  kiln.vol = (1/3) * pi * kiln.mound.dia.calc^2 * (kiln.mound.dia.calc*2*0.55) # cubic feet; ratio of h/d (.55) is based on Greer et al. 2015; Hart 1986; Barnett 2019
  kiln.cords.wood = kiln.vol / 128
  if (i == 1){
    Tar.kilns.measured$Index = i
    Tar.kilns.measured$Site = known.kilns$Site_num[i]
    Tar.kilns.measured$kiln.mound.dia.calc = kiln.mound.dia.calc 
    Tar.kilns.measured$kiln.mound.dia.measured = kiln.mound.dia.measured 
    Tar.kilns.measured$kiln.trench.calc = kiln.trench.calc
    Tar.kilns.measured$kiln.trench.measured = kiln.trench.measured
    Tar.kilns.measured$kiln.total.dia.calc = kiln.total.dia.calc
    Tar.kilns.measured$kiln.total.dia.measured = kiln.total.dia.measured
    Tar.kilns.measured$Kiln_Vol = kiln.vol 
    Tar.kilns.measured$Kiln_cords_wood = kiln.cords.wood 
    Tar.kilns.measured$X_coord =  known.kilns[i,]@coords[,1]
    Tar.kilns.measured$Y_coord = known.kilns[i,]@coords[,2]
  } else {
    Site = known.kilns$Site_num[i]
    Tar.kilns.measured = rbind(Tar.kilns.measured,data.frame("Index" = i, "Site" = Site, "kiln.mound.dia.calc"= kiln.mound.dia.calc, "kiln.mound.dia.measured"= kiln.mound.dia.measured,"kiln.trench.calc"= kiln.trench.calc, "kiln.trench.measured"= kiln.trench.measured, "kiln.total.dia.calc" = kiln.total.dia.calc, "kiln.total.dia.measured" = kiln.total.dia.measured, "Kiln_Vol" = kiln.vol, "Kiln_cords_wood" = kiln.cords.wood, "X_coord" = known.kilns[i,]@coords[,1], "Y_coord" = known.kilns[i,]@coords[,2] ))}
  rownames(Tar.kilns.measured) <- NULL
  setTxtProgressBar(pb,i)
}
proc.time() - ptm
options(warn=0)


# 
# hist((Tar.kilns.measured$kiln.mound.dia.calc - Tar.kilns.measured$kiln.mound.dia.measured))
# 
# 
# hist((Tar.kilns.measured$kiln.mound.dia.calc/Tar.kilns.measured$kiln.mound.dia.measured))
# 
# mean((Tar.kilns.measured$kiln.mound.dia.calc/Tar.kilns.measured$kiln.mound.dia.measured))
# sd((Tar.kilns.measured$kiln.mound.dia.calc/Tar.kilns.measured$kiln.mound.dia.measured))
# 
# mean((Tar.kilns.measured$kiln.mound.dia.calc-Tar.kilns.measured$kiln.mound.dia.measured))
# sd((Tar.kilns.measured$kiln.mound.dia.calc-Tar.kilns.measured$kiln.mound.dia.measured))
# 
# plot(Tar.kilns.measured$kiln.mound.dia.calc,Tar.kilns.measured$kiln.mound.dia.measured)
# model = lm(kiln.mound.dia.measured~kiln.mound.dia.calc,data=Tar.kilns.measured)
# summary(model)
# abline(model) 
# legend("topleft",legend=paste("R2 = ", format(summary(model)$r.squared,digits=3)))


jpeg(filename = "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/figures/Tar_kiln_calc_vc_meas_mound_r.jpeg", width = 10, height = 9, units = 'in', res = 300)
p = ggscatter(Tar.kilns.measured, x ="kiln.mound.dia.calc", y = "kiln.mound.dia.measured" , 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              ylab = "Field Measured Tar Kiln Mound Dia. (m)", xlab = "Model Calculated Tar Kiln Mound Dia. (m)")
ggpar( p, ylim = c(8,24),xlim = c(8,24))
dev.off()

jpeg(filename = "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/figures/Tar_kiln_calc_vc_meas_trench_r.jpeg", width = 10, height = 9, units = 'in', res = 300)
Tar.kilns.measured$kiln.trench.measured[Tar.kilns.measured$kiln.trench.measured<=0]= NA
p  = ggscatter(Tar.kilns.measured, x ="kiln.trench.calc", y = "kiln.trench.measured" , 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               ylab = "Field Measured Tar Kiln Trench (m)", xlab = "Model Calculated Tar Kiln Trench (m)")

ggpar( p, ylim = c(0,5 ),xlim = c(0,5 ))
dev.off()



jpeg(filename = "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/figures/Tar_kiln_calc_vs_meas_tot_dia_r.jpeg", width = 10, height = 9, units = 'in', res = 300)
Tar.kilns.measured$kiln.total.dia.measured[Tar.kilns.measured$kiln.total.dia.measured<=0]= NA
p = ggscatter(Tar.kilns.measured, x ="kiln.total.dia.calc", y = "kiln.total.dia.measured" , 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              ylab = "Field Measured Tar Kiln Total Diameter (m)", xlab = "Model Calculated Tar Kiln Total Diameter (m)")
ggpar( p, ylim = c(13,30 ),xlim = c(13,30))
dev.off()

# comparison of good tar kiln performance vs poor. 

# extract kilns that are within a 10% difference between measured and calculated
# extract kilns that are over 50% difference between measured and calculated












hist(Tar.kilns.measured$Kiln_Radius_line)
Measured.Tar.Kilns = SpatialPointsDataFrame(cbind(Tar.kilns.measured$X_coord,Tar.kilns.measured$Y_coord), Tar.kilns.measured,
                                   proj4string = known.kilns@proj4string, bbox =known.kilns@bbox)
plot(Measured.Tar.Kilns)
shapefile(Measured.Tar.Kilns, "./FMNF_all_tar_kilns_measured.shp",overwrite = T)




# kiln #4 a good visual example for ppt
# error on #369
# Can use volume to estimate how many trees this might represent. See Gonzalez-Benecke et al. 2014



out = (Tar.kilns.measured$kiln.mound.dia.calc - Tar.kilns.measured$kiln.mound.dia.measured)

which(out < -4)


