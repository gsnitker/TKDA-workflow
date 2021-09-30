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
library(stats)
library(cowplot)
library(patchwork)


find.mound = function(kiln, radius = 10, quantile.thresh = .50, plot = F){
  diff = 0
  while (diff <1.5 && radius < 40 ){
    radius = radius + 5
    buffer.trench = buffer(kiln, width = radius)
    buffer.trench.crop = raster::crop(DEM, buffer.trench)
    elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh,1))
    diff = as.numeric(elev.stats[2] - elev.stats[1])
  }
  buffer.trench = buffer(kiln, width = radius + (radius * .35))
  buffer.trench.crop = raster::crop(DEM, buffer.trench)
  elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh,1))
  trench = Which(buffer.trench.crop >= elev.stats[1])
  contour.trench = rasterToContour(trench, nlevels =1)

  points.1 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(mean(buffer.trench@bbox[2,]),mean(buffer.trench@bbox[2,]))))
  intersect.line.1 <- as(points.1,"SpatialLines")
  intersection.points.1 = gIntersection(intersect.line.1, contour.trench, byid=T)
  if (is.null(intersection.points.1) == T) {
    line.radius.1  = NA
    trench.width.1  = NA
  } else {
    below.points = intersection.points.1[intersection.points.1@coords[,1] < kiln@coords[,1]]
    above.points = intersection.points.1[intersection.points.1@coords[,1] > kiln@coords[,1]]
    # find the distances between inner and outer portions of trench and center point
    below.dist = pointDistance(kiln, below.points, lonlat = F)
    below.dist = below.dist[below.dist > 10]
    below.dist = below.dist[which(below.dist < Inf)]
    below.radius = min(below.dist)
    above.dist = pointDistance(kiln, above.points, lonlat = F)
    above.dist = above.dist[above.dist > 10]
    above.dist = above.dist[which(above.dist < Inf)]
    above.radius = min(above.dist)
    # calculate kiln radius
    line.radius.1 = (below.radius + above.radius) / 2
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.1 = (trench.width.below + trench.width.above) / 2
    if (line.radius.1 == Inf) {line.radius.1 = NA}
    if (trench.width.1 == Inf) {trench.width.1 = NA}
    if (trench.width.1 == -Inf) {trench.width.1 = NA}
  }
  
  points.2 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(buffer.trench@bbox[2,1],buffer.trench@bbox[2,2])))
  intersect.line.2 <- as(points.2,"SpatialLines")
  intersection.points.2 = gIntersection(intersect.line.2, contour.trench, byid=T)
  if (is.null(intersection.points.2) == T) {
    line.radius.2  = NA
    trench.width.2  = NA
  } else {
    below.points = intersection.points.2[intersection.points.2@coords[,1] < kiln@coords[,1]]
    above.points = intersection.points.2[intersection.points.2@coords[,1] > kiln@coords[,1]]
    # find the distances between inner and outer portions of trench and center point
    below.dist = pointDistance(kiln, below.points, lonlat = F)
    below.dist = below.dist[below.dist > 10]
    below.dist = below.dist[which(below.dist < Inf)]
    below.radius = min(below.dist)
    above.dist = pointDistance(kiln, above.points, lonlat = F)
    above.dist = above.dist[above.dist > 10]
    above.dist = above.dist[which(above.dist < Inf)]
    above.radius = min(above.dist)
    # calculate kiln radius
    line.radius.2 = (below.radius + above.radius) / 2
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.2 = (trench.width.below + trench.width.above) / 2
    if (line.radius.2 == Inf) {line.radius.2 = NA}
    if (trench.width.2 == Inf) {trench.width.2 = NA}
    if (trench.width.2 == -Inf) {trench.width.2 = NA}
  }
  
  points.3 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(buffer.trench@bbox[2,2],buffer.trench@bbox[2,1])))
  intersect.line.3 <- as(points.3,"SpatialLines")
  intersection.points.3 = gIntersection(intersect.line.3, contour.trench, byid=T)
  if (is.null(intersection.points.3) == T) {
    line.radius.3  = NA
    trench.width.3 = NA 
  } else {
    below.points = intersection.points.3[intersection.points.3@coords[,1] < kiln@coords[,1]]
    above.points = intersection.points.3[intersection.points.3@coords[,1] > kiln@coords[,1]]
    # find the distances between inner and outer portions of trench and center point
    below.dist = pointDistance(kiln, below.points, lonlat = F)
    below.dist = below.dist[below.dist > 10]
    below.dist = below.dist[which(below.dist < Inf)]
    below.radius = min(below.dist)
    above.dist = pointDistance(kiln, above.points, lonlat = F)
    above.dist = above.dist[above.dist > 10]
    above.dist = above.dist[which(above.dist < Inf)]
    above.radius = min(above.dist)
    # calculate kiln radius
    line.radius.3 = (below.radius + above.radius) / 2
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.3 = (trench.width.below + trench.width.above) / 2
    if (line.radius.3 == Inf) {line.radius.3 = NA}
    if (trench.width.3 == Inf) {trench.width.3 = NA}
    if (trench.width.3 == -Inf) {trench.width.3 = NA}
  }
  
  
  points.4 = SpatialPoints(cbind(c(mean(buffer.trench@bbox[1,]),mean(buffer.trench@bbox[1,])), c(buffer.trench@bbox[2,1],buffer.trench@bbox[2,2])))
  intersect.line.4 <- as(points.4,"SpatialLines")
  intersection.points.4 = gIntersection(intersect.line.4, contour.trench, byid=T)
  if (is.null(intersection.points.4) == T) {
    line.radius.4  = NA
    trench.width.4  = NA
  } else {
    below.points = intersection.points.4[intersection.points.4@coords[,2] < kiln@coords[,2]]
    above.points = intersection.points.4[intersection.points.4@coords[,2] > kiln@coords[,2]]
    # find the distances between inner and outer portions of trench and center point
    below.dist = pointDistance(kiln, below.points, lonlat = F)
    below.dist = below.dist[below.dist > 10]
    below.dist = below.dist[which(below.dist < Inf)]
    below.radius = min(below.dist)
    above.dist = pointDistance(kiln, above.points, lonlat = F)
    above.dist = above.dist[above.dist > 10]
    above.dist = above.dist[which(above.dist < Inf)]
    above.radius = min(above.dist)
    # calculate kiln radius
    line.radius.4 = (below.radius + above.radius) / 2
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.4 = (trench.width.below + trench.width.above) / 2
    if (line.radius.4 == Inf) {line.radius.4 = NA}
    if (trench.width.4 == Inf) {trench.width.4 = NA}
    if (trench.width.4 == -Inf) {trench.width.4 = NA}
  }
  
  
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
  radii = na.omit(radii)
  cutoff <- quantile(radii, probs=c(.05, .95), na.rm = T)
  if (length(radii) > 2) {radii = radii[radii >= cutoff[1] & radii <= cutoff[2]]}
  if  (length(radii) == 0) {radii = NA}
  
  trench  = c(trench.width.1, trench.width.2, trench.width.3, trench.width.4)
  trench = na.omit(trench)
  cutoff <- quantile(trench, probs=c(.05, .95), na.rm = T)
  if (length(trench) > 2) {trench = radii[trench >= cutoff[1] & trench <= cutoff[2]]}
  if  (length(trench) == 0) {trench = NA}
  
  line.radius = mean(c(radii), na.rm=T)
  trench = mean(c(trench), na.rm=T)
  kiln.mound.dia.calc = (line.radius * 2) / 3.281
  kiln.trench.calc = trench / 3.281
  return(list(kiln.mound.dia.calc, kiln.trench.calc))
}


# Measure Tar kiln features
# after output has been cleaned in QGIS, re-import to measure
# Import GIS data
known.kilns =  shapefile("/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/FMNF_Final_dataset/FMNF_all_tar_kilns.shp")
DEM = raster("/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/R_models/GIS/Tar_kiln_outputs/FMNF_1m_Lidar_DEM.tif")

# measure kilns in loop
ptm <- proc.time()
options(warn=-1)
Tar.kilns.measured = data.frame("Index" =NA,"kiln.mound.dia.calc"= NA, "kiln.trench.calc"= NA, "kiln.total.dia.calc" = NA, "Kiln_Vol"=NA, "Kiln_cords_wood"=NA, "X_coord" = NA, "Y_coord" = NA )
pb = invisible(txtProgressBar(min = 0, max = length(known.kilns), initial = 0, style =3)) # create progress bar
for (i in 1:length(known.kilns)){ 
  #################
  kiln.metrics = find.mound(kiln = known.kilns[i,], quantile.thresh = .50, plot = F)
  kiln.mound.dia.calc = kiln.metrics[[1]]
  kiln.trench.calc = kiln.metrics[[2]]
  kiln.total.dia.calc = kiln.mound.dia.calc + (kiln.trench.calc*2)

  # find volume of cone based on idealized tar kiln form identified in the literature
  #V=13πr2h
  kiln.vol = (1/3) * pi * kiln.mound.dia.calc^2 * (kiln.mound.dia.calc*2*0.55) # cubic feet; ratio of h/d (.55) is based on Greer et al. 2015; Hart 1986; Barnett 2019
  kiln.cords.wood = kiln.vol / 128
  if (i == 1){
    Tar.kilns.measured$Index = i
    Tar.kilns.measured$kiln.mound.dia.calc = kiln.mound.dia.calc 
    Tar.kilns.measured$kiln.trench.calc = kiln.trench.calc
    Tar.kilns.measured$kiln.total.dia.calc = kiln.total.dia.calc
    Tar.kilns.measured$Kiln_Vol = kiln.vol 
    Tar.kilns.measured$Kiln_cords_wood = kiln.cords.wood 
    Tar.kilns.measured$X_coord =  known.kilns[i,]@coords[,1]
    Tar.kilns.measured$Y_coord = known.kilns[i,]@coords[,2]
  } else {
    Site = known.kilns$Site_num[i]
    Tar.kilns.measured = rbind(Tar.kilns.measured,data.frame("Index" = i, "kiln.mound.dia.calc"= kiln.mound.dia.calc,"kiln.trench.calc"= kiln.trench.calc,  "kiln.total.dia.calc" = kiln.total.dia.calc, "Kiln_Vol" = kiln.vol, "Kiln_cords_wood" = kiln.cords.wood, "X_coord" = known.kilns[i,]@coords[,1], "Y_coord" = known.kilns[i,]@coords[,2] ))}
  rownames(Tar.kilns.measured) <- NULL
  setTxtProgressBar(pb,i)
}
proc.time() - ptm
options(warn=0)

jpeg(filename = "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/figures/FMNF_all_tar_kilns_mound_hist.jpeg", width = 10, height =6 , units = 'in', res = 300)
ggplot(Tar.kilns.measured, aes(x=kiln.mound.dia.calc)) +
  geom_histogram(col = "black",binwidth = c(1)) + 
  xlab("\nMound Diameter (m)") + ylab("Frequency\n") +theme_bw(base_size = 12) + theme(axis.text.x = element_text(size = 12))
dev.off()

mean(Tar.kilns.measured$kiln.mound.dia.calc, na.rm =T)
min(Tar.kilns.measured$kiln.mound.dia.calc, na.rm =T)
max(Tar.kilns.measured$kiln.mound.dia.calc, na.rm =T)
sd(Tar.kilns.measured$kiln.mound.dia.calc, na.rm =T)


Tar.kilns.measured.shp = SpatialPoints(Tar.kilns.measured[,c(7,8)])
Tar.kilns.measured.shp = SpatialPointsDataFrame(Tar.kilns.measured.shp, Tar.kilns.measured)


shapefile(Tar.kilns.measured.shp, "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/FMNF_Final_dataset/FMNF_all_tar_kilns_measured.shp",overwrite = T)
