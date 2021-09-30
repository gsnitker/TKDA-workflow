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
setwd("/Users/grantsnitker/Dropbox/Smoke/Colonial_Urban_Economy_Project/GIS/Tar_kiln_outputs/")
# #85 196 252 259 265 270
# kiln = known.kilns[270,]
# kiln$Site_num
# radius = 10
# quantile.thresh = .50
# plot = T
# Load function needed to find the trench surrounding the kiln and thus the kilns radius.


RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

#find.trench(kiln = known.kilns[36,])
find.mound = function(kiln, radius = 10, quantile.thresh = .50, plot = F){
  #kiln.mound.dia.measured = kiln$TK_mnd_dia
 # kiln.trench.measured = kiln$TK_tch_wdt
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
 #plot(buffer.trench.crop)
    #contour = rasterToContour(buffer.trench.crop,nlevels =10)
    # elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh+.1,1))
    # diff = as.numeric(elev.stats[2] - elev.stats[1])
    # mound = Which(buffer.trench.crop <= elev.stats[1])
    # plot(ring)
    # contour.mound = rasterToContour(ring, nlevels =1)
    # plot(contour.mound, add = T)
    elev.stats = quantile(buffer.trench.crop, probs = c(quantile.thresh,1))
    #diff = as.numeric(elev.stats[2] - elev.stats[1])
    trench = Which(buffer.trench.crop >= elev.stats[1])
    #plot(buffer.trench.crop)
    contour.trench = rasterToContour(trench, nlevels =1)
    #plot(contour.trench, add = T)
  
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
    #plot(below.points, add =T, col = "blue")
    #plot(above.points, add =T)
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
  radii
  radii = na.omit(radii)
  cutoff <- quantile(radii, probs=c(.05, .95), na.rm = T)
  if (length(radii) > 2) {radii = radii[radii >= cutoff[1] & radii <= cutoff[2]]}
  if  (length(radii) == 0) {radii = NA}
  
  trench  = c(trench.width.1, trench.width.2, trench.width.3, trench.width.4)
  trench
  trench = na.omit(trench)
  cutoff <- quantile(trench, probs=c(.05, .95), na.rm = T)
  if (length(trench) > 2) {trench = radii[trench >= cutoff[1] & trench <= cutoff[2]]}
  if  (length(trench) == 0) {trench = NA}
  
  
  
  line.radius = mean(c(radii), na.rm=T)
  trench = mean(c(trench), na.rm=T)
  kiln.mound.dia.calc = (line.radius * 2) / 3.281
  kiln.trench.calc = trench / 3.281
  kiln.mound.dia.calc
  kiln.mound.dia.measured
  return(list(kiln.mound.dia.calc, kiln.mound.dia.measured, kiln.trench.calc, kiln.trench.measured))
  #return(list(kiln.mound.dia.calc, kiln.mound.dia.measured))
  }

# kiln = known.kilns[119,]
# kiln.metrics = find.mound(kiln = kiln, quantile.thresh = .50, plot = F)
# kiln.mound.dia.calc = kiln.metrics[[1]]
# kiln.mound.dia.measured = kiln.metrics[[2]]


find.trench = function(kiln){
  kiln.trench.measured = kiln$TK_tch_wdt
  radius = kiln.mound.dia.calc / 2
  buffer.trench = buffer(kiln, width = (radius * 3.281) + 20)
  buffer.trench.crop = raster::crop(DEM, buffer.trench)
  buffer.trench.mask = raster::mask(buffer.trench.crop, buffer.trench)
  
 elev.stats = quantile(buffer.trench.mask, probs = c(.4,1))
   trench = Which(buffer.trench.mask <= elev.stats[1])
   #plot(buffer.trench.mask)
   contour.trench = rasterToContour(trench, nlevels =1)
  # plot(contour.trench, add = T)

  # im = as.cimg(buffer.trench.crop)
  # canny <- cannyEdges(im, alpha=.5)
  # px.im = as.cimg(canny)
  # rast =raster(as.matrix(canny), template = buffer.trench.crop)
  # #plot(rast)
  # contour.trench = rasterToContour(rast)
  # plot(contour.trench)
  
  
  points.1 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(mean(buffer.trench@bbox[2,]),mean(buffer.trench@bbox[2,]))))
  intersect.line.1 <- as(points.1,"SpatialLines")
  intersection.points.1 = gIntersection(intersect.line.1, contour.trench, byid=T)
  if (is.null(intersection.points.1) == T) {
    trench.width.1  = NA
  } else {
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
    # calculate kiln trench
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.1 = (trench.width.below + trench.width.above) / 2
    if (trench.width.1 == Inf) {trench.width.1 = NA}
    if (trench.width.1 == -Inf) {trench.width.1 = NA}
  }
  
  
  points.2 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(buffer.trench@bbox[2,1],buffer.trench@bbox[2,2])))
  intersect.line.2 <- as(points.2,"SpatialLines")
  intersection.points.2 = gIntersection(intersect.line.2, contour.trench, byid=T)
  if (is.null(intersection.points.2) == T) {
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
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.2 = (trench.width.below + trench.width.above) / 2
    if (trench.width.2 == Inf) {trench.width.2 = NA}
    if (trench.width.2 == -Inf) {trench.width.2 = NA}
  }
  
  
  points.3 = SpatialPoints(cbind(c(buffer.trench@bbox[1,1],buffer.trench@bbox[1,2]), c(buffer.trench@bbox[2,2],buffer.trench@bbox[2,1])))
  intersect.line.3 <- as(points.3,"SpatialLines")
  intersection.points.3 = gIntersection(intersect.line.3, contour.trench, byid=T)
  if (is.null(intersection.points.3) == T) {
    trench.width.3 = NA 
  } else {
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
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.3 = (trench.width.below + trench.width.above) / 2
    if (trench.width.3 == Inf) {trench.width.3 = NA}
    if (trench.width.3 == -Inf) {trench.width.3 = NA}
  }
  
  
  points.4 = SpatialPoints(cbind(c(mean(buffer.trench@bbox[1,]),mean(buffer.trench@bbox[1,])), c(buffer.trench@bbox[2,1],buffer.trench@bbox[2,2])))
  intersect.line.4 <- as(points.4,"SpatialLines")
  intersection.points.4 = gIntersection(intersect.line.4, contour.trench, byid=T)
  if (is.null(intersection.points.4) == T) {
     trench.width.4  = NA
  } else {
    below.points = intersection.points.4[intersection.points.4@coords[,2] < kiln@coords[,2]]
    above.points = intersection.points.4[intersection.points.4@coords[,2] > kiln@coords[,2]]
   # plot(below.points, add =T, col = "blue")
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
    trench.width.above = max(above.dist) - min(above.dist)
    trench.width.below = max(below.dist) - min(below.dist)
    trench.width.4 = (trench.width.below + trench.width.above) / 2
    if (trench.width.4 == Inf) {trench.width.4 = NA}
    if (trench.width.4 == -Inf) {trench.width.4 = NA}
  }

  
  trench  = c(trench.width.1, trench.width.2, trench.width.3, trench.width.4)

  trench
  trench = na.omit(trench)
  cutoff <- quantile(trench, probs=c(.05, .95), na.rm = T)
  if (length(trench) > 2) {trench = trench[trench >= cutoff[1] & trench <= cutoff[2]]}
  #trench = trench[trench > 3 & trench < 15]
  
  if  (length(trench) == 0) {trench = NA}
  
  trench = mean(c(trench), na.rm=T)
  kiln.trench.calc = trench / 3.281
  kiln.trench.calc
  kiln.trench.measured
  return(list(kiln.trench.calc, kiln.trench.measured))
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
Tar.kilns.measured = data.frame("Index" =NA,"Site" = NA,"kiln.mound.dia.calc"= NA, "kiln.mound.dia.measured"= NA, "kiln.trench.calc"= NA, "kiln.trench.measured"= NA, "kiln.total.dia.calc" = NA, "kiln.total.dia.measured" = NA, "Kiln_Vol"=NA, "Kiln_cords_wood"=NA, "X_coord" = NA, "Y_coord" = NA )
pb = invisible(txtProgressBar(min = 0, max = length(known.kilns), initial = 0, style =3)) # create progress bar
#pb = invisible(txtProgressBar(min = 0, max = 100, initial = 0, style =3)) # create progress bar
for (i in 1:length(known.kilns)){ 
#for (i in 1:50){ 
  # subset by kiln in chunk 
  #################
  kiln.metrics = find.mound(kiln = known.kilns[i,], quantile.thresh = .50, plot = F)
  kiln.mound.dia.calc = kiln.metrics[[1]]
  kiln.mound.dia.measured = kiln.metrics[[2]]
   trench.metrics = find.trench(kiln = known.kilns[i,])
  
  kiln.trench.calc = trench.metrics[[1]]
  kiln.trench.measured = trench.metrics[[2]]
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
# 
jpeg(filename = "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/figures/Tar_kiln_calc_vc_meas_mound_r.jpeg", width = 10, height = 9, units = 'in', res = 300)
p = ggscatter(Tar.kilns.measured, x ="kiln.mound.dia.calc", y = "kiln.mound.dia.measured" ,
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "Field Measured Tar Kiln Mound Dia. (m)", xlab = "Model Calculated Tar Kiln Mound Dia. (m)")
ggpar( p, ylim = c(8,24),xlim = c(8,24))

dev.off()
cor.dia <- cor.test(Tar.kilns.measured$kiln.mound.dia.measured, Tar.kilns.measured$kiln.mound.dia.calc,
                method = "pearson")
jpeg(filename = "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/figures/Tar_kiln_calc_vc_meas_trench_r.jpeg", width = 10, height = 9, units = 'in', res = 300)
Tar.kilns.measured$kiln.trench.measured[Tar.kilns.measured$kiln.trench.measured<=0]= NA
p  = ggscatter(Tar.kilns.measured.no.na, x ="kiln.trench.calc", y = "kiln.trench.measured" , 
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
ggpar( p, ylim = c(11,30 ),xlim = c(11,30))
dev.off()


# RMSE
#mound
RMSE(Tar.kilns.measured$kiln.mound.dia.calc, Tar.kilns.measured$kiln.mound.dia.measured)
mean(Tar.kilns.measured$kiln.mound.dia.calc)
mean(Tar.kilns.measured$kiln.mound.dia.measured)

max(Tar.kilns.measured$kiln.mound.dia.calc)
max(Tar.kilns.measured$kiln.mound.dia.measured)

min(Tar.kilns.measured$kiln.mound.dia.calc)
min(Tar.kilns.measured$kiln.mound.dia.measured)

#trench
Tar.kilns.measured.no.na = na.omit(Tar.kilns.measured)
RMSE(Tar.kilns.measured.no.na$kiln.trench.calc, Tar.kilns.measured.no.na$kiln.trench.measured)
mean(Tar.kilns.measured.no.na$kiln.trench.calc)
mean(Tar.kilns.measured.no.na$kiln.trench.measured)

max(Tar.kilns.measured.no.na$kiln.trench.calc)
max(Tar.kilns.measured.no.na$kiln.trench.measured)

min(Tar.kilns.measured.no.na$kiln.trench.calc)
min(Tar.kilns.measured.no.na$kiln.trench.measured)

#overall.kiln
RMSE(Tar.kilns.measured.no.na$kiln.total.dia.calc, Tar.kilns.measured.no.na$kiln.total.dia.measured)
mean(Tar.kilns.measured.no.na$kiln.total.dia.calc)
mean(Tar.kilns.measured.no.na$kiln.total.dia.measured)

max(Tar.kilns.measured.no.na$kiln.total.dia.calc)
max(Tar.kilns.measured.no.na$kiln.total.dia.measured)

min(Tar.kilns.measured.no.na$kiln.total.dia.calc)
min(Tar.kilns.measured.no.na$kiln.total.dia.measured)


# evaluation of tar kiln trench performance
Tar.kilns.measured$trench_diff_per = abs((Tar.kilns.measured$kiln.trench.calc/Tar.kilns.measured$kiln.trench.measured) - 1)
good.trenches = subset(Tar.kilns.measured, trench_diff_per < 0.05)
bad.trenches = subset(Tar.kilns.measured, trench_diff_per > 0.50)

good.trenches$type= "Less than 5% Divergence"
bad.trenches$type= "More than 50% Divergence"
trenches.eval.dataset = rbind(good.trenches,bad.trenches )

# disturbance data collection
fun <- function(){
  disturb <- readline("Disturbance: 0 or 1?")  
  disturb <- as.numeric(unlist(strsplit(disturb, ",")))
  out <- disturb 
  return(out)
}


for (i in 1:length(trenches.eval.dataset[,1])) {
  
  kiln = known.kilns[trenches.eval.dataset$Index[i],]
  radius = trenches.eval.dataset$kiln.mound.dia.calc[i] / 2
  buffer.trench = buffer(kiln, width = (radius * 3.281) + 50) # 20
  buffer.trench.crop = raster::crop(DEM, buffer.trench)
  plot(buffer.trench.crop, col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
  disturb = fun() 
  trenches.eval.dataset$disturb[i] = disturb
}

i =1

# calculate mound prominence and surrounding terrain
for (i in 1:length(trenches.eval.dataset[,1])) {
 kiln = known.kilns[trenches.eval.dataset$Index[i],]
 radius = trenches.eval.dataset$kiln.mound.dia.calc[i] / 2
 buffer.trench = buffer(kiln, width = (radius * 3.281) + 20) # 20
 buffer.trench.crop = raster::crop(DEM, buffer.trench)
 #plot(buffer.trench.crop)
 

 buffer.kiln.shp = buffer(kiln, width = (radius * 3.281)-9) 
 buffer.kiln = raster::crop(DEM, buffer.kiln.shp)
# plot(buffer.kiln.shp , add = T)
 
 buffer.kiln.plus = buffer(kiln, width = (radius * 3.281)+9)
 buffer.trench.inverse = raster::mask(buffer.trench.crop,  buffer.kiln.plus, inverse = T)
# plot(buffer.trench.inverse)
 buffer.trench.inverse.slope = terrain(buffer.trench.inverse, opt="slope", unit="degrees", neighbors=8)
 buffer.trench.inverse.TRI = terrain(buffer.trench.inverse, opt="TRI", unit="degrees", neighbors=8)
 buffer.trench.inverse.rough = terrain(buffer.trench.inverse, opt="roughness", unit="degrees", neighbors=8)
 
 mean.kiln.height = cellStats(buffer.kiln, mean)
 mean.surrounding.height = cellStats(buffer.trench.inverse, mean)
 mean.elev = cellStats(buffer.trench.crop, mean)
 mean.TRI = cellStats(buffer.trench.inverse.TRI, mean)
 mean.slope = cellStats(buffer.trench.inverse.slope, mean)
 mean.rough = cellStats(buffer.trench.inverse.rough, mean)
 slope.range = cellStats(buffer.trench.inverse.slope, max) - cellStats(buffer.trench.inverse.slope, min)
 slope.sd = cellStats(buffer.trench.inverse.slope, sd)
 
 quant = quantile(buffer.trench.crop, probs=c(.25, .75), na.rm = T)
 quant.diff =  quant[2] - quant[1]
 
mound.prom = abs(mean.kiln.height - mean.elev)
mound.prom 
trenches.eval.dataset$mound.prom[i] = mound.prom 
trenches.eval.dataset$surrounding.slope[i] = mean.slope
trenches.eval.dataset$surrounding.TRI[i] = mean.TRI 
trenches.eval.dataset$surrounding.roughness[i] = mean.rough
trenches.eval.dataset$surrounding.slope.range[i] = slope.range
trenches.eval.dataset$surrounding.slope.sd[i] = slope.sd

trenches.eval.dataset$elev.diff[i] = quant.diff 
}

trenches.eval.dataset$disturb = disturbnace
# boxplots


# ggplot(trenches.eval.dataset, aes(x=type, y=surrounding.slope)) +
#   geom_boxplot() +
#   ggtitle("Off mound terrain flatness") +
#   xlab("") + ylab("Mean Slope (degs)")
# 
# ggplot(trenches.eval.dataset, aes(x=type, y=surrounding.TRI)) +
#   geom_boxplot() +
#   ggtitle("Off mound terrain flatness") +
#   xlab("") + ylab("Terrain Roughness Index (TRI)")
# 
# ggplot(trenches.eval.dataset, aes(x=type, y=surrounding.roughness)) +
#   geom_boxplot() +
#   ggtitle("Off mound terrain flatness") +
#   xlab("") + ylab("Roughness")

plot1 = ggplot(trenches.eval.dataset, aes(x=type, y=mound.prom)) +
  geom_boxplot() +
  ggtitle("Difference in Mound Prominence") +
  xlab("") + ylab("Meters\n")+theme_bw() +theme_bw(base_size = 12) + theme(axis.text.x = element_text(size = 12))

plot2 = ggplot(trenches.eval.dataset, aes(x=type, y=surrounding.slope.range)) +
  geom_boxplot() +
  ggtitle("Range of Off-Mound Slope Values") +
  xlab("") + ylab("Degrees\n") +theme_bw(base_size = 12) + theme(axis.text.x = element_text(size = 12))

# ggplot(trenches.eval.dataset, aes(x=type, y=elev.diff)) +
#   geom_boxplot() +
#   ggtitle("Diiference Between 25% and 75% Elevation Quantiles") +
#   xlab("") + ylab("Degrees\n") +theme_bw(base_size = 12) + theme(axis.text.x = element_text(size = 12))

plot3 = ggplot(trenches.eval.dataset, aes(x=type, y=disturb/length(disturbnace))) +
  geom_col() +
  ggtitle("Kilns with Mound Disturbance") +
  xlab("") + ylab("Proportion\n") +theme_bw(base_size = 12) + theme(axis.text.x = element_text(size = 12))

jpeg(filename = "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/figures/tar_kiln_trench_preformance_limitations.jpeg", width = 21, height = 9, units = 'in', res = 300)
p = plot1 + plot2 + plot3
p + plot_annotation(tag_levels = 'A')
dev.off()


plot_grid(plot1, plot2, plot3, col = 3, align = c("v"))
plot_grid(plot1, plot2, plot3, align = "h", axis = "b", rel_widths = c(1, 1,1))


# # t-tests
# good = subset(trenches.eval.dataset, type == "Less than 5% Divergence")
# bad = subset(trenches.eval.dataset, type == "More than 50% Divergence")
# t.test(good$mound.prom, bad$mound.prom)
# 
# good = subset(trenches.eval.dataset, type == "Less than 5% Divergence")
# bad = subset(trenches.eval.dataset, type == "More than 50% Divergence")
# t.test(good$mound.prom, bad$mound.prom)


# test visualizations
ggscatter(bad.trenches, x ="kiln.trench.calc", y = "kiln.trench.measured" , 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          ylab = "Field Measured Tar Kiln Trench (m)", xlab = "Model Calculated Tar Kiln Trench (m)")





# kiln #4 a good visual example for ppt
# error on #369
# Can use volume to estimate how many trees this might represent. See Gonzalez-Benecke et al. 2014



out = (Tar.kilns.measured$kiln.mound.dia.calc - Tar.kilns.measured$kiln.mound.dia.measured)

which(out < -4)


