#Validate tar kilns identified through model by comparison to those measure in the field
rm(list = ls())
setwd("/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/R_models/Tar_Kilns/")

# Libraries
library(imager)
library(rgeos)
library(sp)
library(raster)
library(tools)
library(dplyr)
library(tidyr)
library(rgdal)
library(sf)
library(tidyverse)
library(reshape2)
library(geosphere)
library(parallel)
library(doParallel)



# load tar kiln validation set
# need to decide if all this should be done in UTMs or the SC projection that the other layers are in
val.kilns = shapefile("/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/shapefiles/FMNF_TKP_Validation_final_NAD83_2011_ft.shp")
val.kilns = val.kilns[val.kilns$PI_Company == "ACC",]


plot(val.kilns)

# simple plots of distributions
# Total Dia
hist(val.kilns$TK_tot_dia, breaks = 10, main = "Distribution of Total Tar Kiln Diameter", xlab = "Meters")
mtext("Validation Set",3)
summary(val.kilns$TK_tot_dia)

# Mound Dia
hist(val.kilns$TK_mnd_dia, breaks = 10,main = "Distribution of Tar Kiln Mound Diameter", xlab = "Meters")
mtext("Validation Set",3)
summary(val.kilns$TK_mnd_dia)

# Trench width
trenches =((val.kilns$TK_tot_dia-val.kilns$TK_mnd_dia)/2)[((val.kilns$TK_tot_dia-val.kilns$TK_mnd_dia)/2) >= 0]
hist(trenches, breaks = 10, main = "Distribution of Tar Kiln Trench Width", xlab = "Meters")
mtext("Validation Set",3)
summary(trenches)

# Trench to mound ratio
trenches =((val.kilns$TK_tot_dia-val.kilns$TK_mnd_dia)/2)[((val.kilns$TK_tot_dia-val.kilns$TK_mnd_dia)/2) >= 0]
hist(trenches/val.kilns$TK_mnd_dia, breaks = 10, main = "Distribution of Tar Kiln Mound:Trench Ratio", xlab = "Meters")
mtext("Validation Set",3)
summary(trench.mound.ratio)


# Testing the Feature Detection Algorithm 

# Model Test 1: Detection Accuracy true positives
# load full raster of FMNF (1m res)
dem.1m = raster("/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/R_models/GIS/Tar_kiln_outputs/FMNF_1m_Lidar_DEM.tif")
kiln_names = val.kilns$Site_code

# define validation loop
# Do not run, already completed

#   pb = invisible(txtProgressBar(min = 0, max = length(val.kilns), initial = 0, style =3)) # create progress bar
# for (i in 1:length(val.kilns)){
#   focal.kiln = val.kilns[i,]
#   kiln_name = kiln_names[i]
#   file = file.exists(paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/"))
#   if (file == FALSE){
#     dir.create(paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/"))}
#   e = extent(focal.kiln) + 150# buffer around the focal kiln by 1km ( or 1312 pixels at 2.5 ft each) 
#   focal.dem = crop(dem.1m, e)
#   im = as.cimg(focal.dem)
#   plot(im)
#   px <- cannyEdges(im, alpha=.5)
#   plot(px)
#   px.im = as.cimg(px)
#   plot(px.im)
#   save.image(px.im, file = paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/", kiln_name,".jpg")) # on second thought, these need to go into their own folders so that each output can be matched up with the correct spatial information. Might need to save the rasters too.  
#   writeRaster(focal.dem, file = paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/", kiln_name,".tif"),format="GTiff", overwrite = T) # on second thought, these need to go into their own folders so that each output can be matched up with the correct spatial information. Might need to save the rasters too.  
#   setTxtProgressBar(pb,i)
#   }
#   
#   # Run FIJI portion in parallel
#   FIJI.multi.core = function(index){
#   # # FIJI 
#     kiln_name = kiln_names[index]
#     print(paste0("Starting ", kiln_name, "...", "(",index," of ", length(kiln_names),")"))
# out = system(paste("/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx -macro", paste(getwd(), ("/ImageJ_macro_tar_kilns_V6.txt"), sep = ""), paste(getwd(), ("/data/Canny_Edges_Images/Validation/"), (kiln_name),("/"),("*"), getwd(), ("/imageJ_output/"), kiln_name, sep = "")))
#   print(paste0(kiln_name, " is complete."))
#     # 
#     # # UCB Vision Sciences Plugins  - Hough Circle Transformation
#     # #Advanced mode
#     # # min radius: 8
#     # # max radius 30
#     # # radius search increment:1
#     # # default max number to be found: 65535
#     # # Hough score threshold: 0.5
#     # # transform resolution: 1000
#     # # clear neighbors radius ratio: 1.0
#   return()
#   }
#   
#   # run parallel FIJI
#   # cl = makeCluster(num.cores, type = "PSOCK")
#   # clusterExport(cl, varlist = ls())
#   # invisible(clusterEvalQ(cl, c(library(imager),library(rgeos),library(sp),
#   #                              library(raster),library(tools),library(rgdal),
#   #                              library(sf),library(tidyverse),library(reshape2),
#   #                              library(geosphere),library(parallel),library(doParallel)))) # all required libraries
#   # ptm <- proc.time() # start time
#   # output = parLapply(cl,1:length(kiln_names),FIJI.multi.core)
#   # stopCluster(cl)
#   
#   # Or single thread FIJI
#   output = lapply(1:length(kiln_names),FIJI.multi.core)


# Georectify and test outputs
for (index in 1:length(kiln_names)){ # running the loop backward to avoid headaches of the first tar kiln
    kiln_name = kiln_names[index]
  chunk.list.all = list.files(paste0("./data/Canny_Edges_Images/Validation/",kiln_name, "/"), pattern = "\\.jpg$") %>% file_path_sans_ext(.) %>% noquote(.)
  # load the results from ImageJ
  rast = raster(paste0("./data/Canny_Edges_Images/Validation/",kiln_name, "/", kiln_name,".tif"))
  imageJ.output = read.csv(paste0("./imageJ_output/", kiln_name, "_imageJ_results.csv"))
  centroid  = xyFromCell(rast,  cellFromRowCol(rast,30,30), spatial = T)
  
  if (ncol(imageJ.output)<5){
    
    ID = NA
    X..pixels.= 1
    Y..pixels. = 1
    Radius..pixels. = 10
    Score = NA 
    nCircles =NA
    Resolution = NA
    Frame..slice... = NA
    Method = "No kilns detected"
  
    imageJ.output = cbind(imageJ.output, ID, X..pixels., Y..pixels., Radius..pixels., Score ,nCircles ,Resolution ,Frame..slice..., Method)
    imageJ.output$r_units = imageJ.output$Radius..pixels. * res(rast)[1] # ft in this case; this only works if the measurements are in linear units (e.g. meters, feet)
    
    
} else {
    # replace zeros with NAs in chunk column
    imageJ.output$chunk[imageJ.output$chunk == 0] = NA
    imageJ.output = imageJ.output %>% tidyr::fill(chunk, .direction = c("down"))
    imageJ.output$r_units = imageJ.output$Radius..pixels. * res(rast)[1] # ft in this case; this only works if the measurements are in linear units (e.g. meters, feet)
    imageJ.output = imageJ.output[imageJ.output$ID != 0, ]
   borderline.radius = subset(imageJ.output, Radius..pixels. <= 6)
   imageJ.output = imageJ.output[imageJ.output$Radius..pixels. > 6, ]
   borderline.score = subset(borderline.radius, Score >= .55)
   imageJ.output  = rbind(imageJ.output, borderline.score )
}
  if (nrow(imageJ.output)<1){
    imageJ.output[1,] = imageJ.output[1,]
    imageJ.output$chunk = kiln_name
    imageJ.output$ID = NA
    imageJ.output$X..pixels.= 1
    imageJ.output$Y..pixels. = 1
    imageJ.output$Radius..pixels. = 10
    imageJ.output$ Score = NA 
    imageJ.output$ nCircles =NA
    imageJ.output$Resolution = NA
    imageJ.output$Frame..slice... = NA
    imageJ.output$Method = "No kilns detected"
    
    imageJ.output$r_units = imageJ.output$Radius..pixels. * res(rast)[1] # ft in this case; this only works if the measurements are in linear units (e.g. meters, feet)
  }
    
      chunk.list = intersect(imageJ.output$chunk, chunk.list.all)
      
      num.out = length(chunk.list)
      num.seq = c(1:num.out)
      imageJ.output.chunk = subset(imageJ.output, chunk == chunk.list[1])
      locations = cellFromRowCol(rast[[match(chunk.list[1],chunk.list.all)]], imageJ.output.chunk$Y..pixels., imageJ.output.chunk$X..pixels.)
      kiln.points = xyFromCell(rast[[match(chunk.list[1],chunk.list.all)]], locations, spatial =T)
      
      kiln.data = SpatialPointsDataFrame(kiln.points@coords, imageJ.output, bbox = kiln.points@bbox)
      kiln.data$X_coord = kiln.points@coords[,1]
      kiln.data$Y_coord = kiln.points@coords[,2]

      
      if (nrow(kiln.data)>1){
      dist.from.centroid = pointDistance(centroid, kiln.data, lonlat = F)
      keep.index = which.min(dist.from.centroid)
      
      kiln.data = kiln.data[keep.index,]}
      


      # test 1: Is the circular object slightly more prominent than the surrounding area?
      pb = invisible(txtProgressBar(min = 0, max = length(kiln.data), initial = 0, style =3)) # create progress bar
      for (i in 1:length(kiln.data)){
        kiln = kiln.data[i,]
        buffer.big = buffer(kiln, width =  kiln$r_units * 2)
        buffer.small= buffer(kiln, width = kiln$r_units * .8)
        buffer.center = buffer(kiln, width = kiln$r_units * .3)
        buffer.crop = raster::crop(rast, buffer.big)
        buffer.mask.outside = raster::mask(buffer.crop, buffer.small, inverse  = T)
        buffer.mask.inside = raster::mask(buffer.crop, buffer.small, inverse  = F)
        buffer.mask.donut = raster::mask(buffer.mask.inside, buffer.center, inverse  = T)
        outside.elev = cellStats(buffer.mask.outside,mean) + .1 #.4? half a meter adjustment...need some quantitative reasoning for this. 
        kiln.elev = cellStats(buffer.mask.donut,mean)
        #kiln.elev = cellStats(buffer.mask.inside,mean)
        kiln.data$kiln.elev[i] = kiln.elev 
        kiln.data$outside.elev[i] = outside.elev
        if(kiln.elev > outside.elev) {
          kiln.data$kiln_1[i] = 1
        } else {  kiln.data$kiln_1[i] = 0}
        
        # test 2: Is the circular object surrounded by a trench?
        buffer.trench = buffer(kiln, width = kiln$r_units * 1.8)
        buffer.trench.crop = raster::crop(rast, buffer.trench)
        elev.stats = quantile(buffer.crop, probs = seq(0, 1, 0.1))
        ring = Which(buffer.crop <= elev.stats[9])
        contour.ring = rasterToContour(ring, nlevels =1)
        #plot(contour.ring)
        
        points = SpatialPoints(cbind(c(contour.ring@bbox[1,1],contour.ring@bbox[1,2]), c(kiln@coords[2], kiln@coords[2])))
        intersect.line <- as(points,"SpatialLines")
        intersection.points = gIntersection(intersect.line, contour.ring, byid=T)
        #plot(intersection.points, add =T)
        if (is.null(intersection.points) == TRUE) {intersection.points = 0}
        
        if(length(intersection.points) >= 2) { # maybe this should 4? two edges of the trench?
          kiln.data$kiln_2[i] = 1
        } else {  kiln.data$kiln_2[i] = 0}
        
        if(kiln.data$kiln_1[i] == 1 & kiln.data$kiln_2[i] == 1) {
          kiln.data$kiln[i] = 1
        } else {kiln.data$kiln[i] = 0}
        
        setTxtProgressBar(pb,i)
      }
      
      kiln.data$Site = as.character(kiln_name)
      kiln.data = remove.duplicates(kiln.data)
      all.possible.kilns = kiln.data
      if (all.possible.kilns$Method ==  "No kilns detected"){
        all.possible.kilns$kiln_1 = 0
        all.possible.kilns$kiln_2 = 0
        all.possible.kilns$kiln = 0
      }
      tar.kilns = subset(kiln.data, kiln.data$kiln== 1)
      if (index == 1) {
        validation.dataset = all.possible.kilns
      } else {
     validation.dataset = rbind(validation.dataset,all.possible.kilns)
        }
      #plot(rast)
      #plot(known.tar.kilns, pch = 1, add = T)
      #plot(all.possible.kilns, add = T)
      #plot(tar.kilns, col ="red", add = T)
      
      #save shapefile outputs
      #shapefile(tar.kilns, paste0("./output/",lidar.name,"_tar_kilns.shp"), overwrite = T)
     # shapefile(all.possible.kilns, paste0("./output/",lidar.name,"_all_possible_tar_kilns.shp"), overwrite = T)
      #plot(rast, main = paste(kiln_name))
      #plot(tar.kilns, add = T)
      
      
      print(paste0(kiln_name, " was successfully completed."))
      #do.call(file.remove, list(list.files("./data/Canny_Edges_Images/", full.names = TRUE)))
   
}

val.data = validation.dataset@data

val.data.pos = subset(val.data, kiln == 1)
out = val.kilns[match(val.data.pos$Site, val.kilns$Site_code),]
shapefile(out, "/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/shapefiles/pos_validated_tar_kilns.shp",overwrite = T)

View(val.data.pos)
View(out@data)

circle.dect.rate = sum(val.data$Frame..slice..., na.rm = T) / length(val.data$Frame..slice...)
circle.dect.rate

prominence.rate = sum(val.data$kiln_1, na.rm = T) / length(val.data$Frame..slice...)
prominence.rate

trench.rate = sum(val.data$kiln_2, na.rm = T) / length(val.data$Frame..slice...)
trench.rate

kiln.id.rate = sum(val.data$kiln, na.rm = T) / length(val.data$Frame..slice...)
kiln.id.rate

# Metrics for Balanced Accuary #1
TP = kiln.id.rate * length(val.kilns)
FN = (1 - kiln.id.rate) * length(val.kilns)


#################################################


# Model Test 1: Detection Accuracy true negatives
# load full raster of FMNF (1m res)
rdm.val.kilns = shapefile("/Users/grantsnitker/Dropbox/Smoke/USFS_ORISE/Tar_Kilns/validation/shapefiles/Random_points_not_near_known_tar_kilns_NAD83_2011_ft.shp")
kiln_names = rdm.val.kilns$ID_code

# define validation loop
# do not run, already completed
# pb = invisible(txtProgressBar(min = 0, max = length(rdm.val.kilns), initial = 0, style =3)) # create progress bar
# for (i in 1:length(rdm.val.kilns)){
#   focal.kiln = rdm.val.kilns[i,]
#   kiln_name = kiln_names[i]
#   file = file.exists(paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/"))
#   if (file == FALSE){
#     dir.create(paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/"))}
#   e = extent(focal.kiln) + 150# buffer around the focal kiln by 1km ( or 1312 pixels at 2.5 ft each) 
#   focal.dem = crop(dem.1m, e)
#   im = as.cimg(focal.dem)
#   plot(im)
#   px <- cannyEdges(im, alpha=.5)
#   plot(px)
#   px.im = as.cimg(px)
#   plot(px.im)
#   save.image(px.im, file = paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/", kiln_name,".jpg")) # on second thought, these need to go into their own folders so that each output can be matched up with the correct spatial information. Might need to save the rasters too.  
#   writeRaster(focal.dem, file = paste0("./data/Canny_Edges_Images/Validation/",kiln_name,"/", kiln_name,".tif"),format="GTiff", overwrite = T) # on second thought, these need to go into their own folders so that each output can be matched up with the correct spatial information. Might need to save the rasters too.  
#   setTxtProgressBar(pb,i)
# }
# 
# # Run FIJI portion in parallel
# FIJI.multi.core = function(index){
#   # # FIJI 
#   kiln_name = kiln_names[index]
#   print(paste0("Starting ", kiln_name, "...", "(",index," of ", length(kiln_names),")"))
#   out = system(paste("/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx -macro", paste(getwd(), ("/ImageJ_macro_tar_kilns_V6.txt"), sep = ""), paste(getwd(), ("/data/Canny_Edges_Images/Validation/"), (kiln_name),("/"),("*"), getwd(), ("/imageJ_output/"), kiln_name, sep = "")))
#   print(paste0(kiln_name, " is complete."))
#   # 
#   # # UCB Vision Sciences Plugins  - Hough Circle Transformation
#   # #Advanced mode
#   # # min radius: 8
#   # # max radius 30
#   # # radius search increment:1
#   # # default max number to be found: 65535
#   # # Hough score threshold: 0.5
#   # # transform resolution: 1000
#   # # clear neighbors radius ratio: 1.0
#   return()
# }
# 
# # run parallel FIJI
# # cl = makeCluster(num.cores, type = "PSOCK")
# # clusterExport(cl, varlist = ls())
# # invisible(clusterEvalQ(cl, c(library(imager),library(rgeos),library(sp),
# #                              library(raster),library(tools),library(rgdal),
# #                              library(sf),library(tidyverse),library(reshape2),
# #                              library(geosphere),library(parallel),library(doParallel)))) # all required libraries
# # ptm <- proc.time() # start time
# # output = parLapply(cl,1:length(kiln_names),FIJI.multi.core)
# # stopCluster(cl)
# 
# # Or single thread FIJI
# output = lapply(1:length(kiln_names),FIJI.multi.core)


# Georectify and test outputs
for (index in 1:length(kiln_names)){ # running the loop backward to avoid headaches of the first tar kiln
  kiln_name = kiln_names[index]
  chunk.list.all = list.files(paste0("./data/Canny_Edges_Images/Validation/",kiln_name, "/"), pattern = "\\.jpg$") %>% file_path_sans_ext(.) %>% noquote(.)
  # load the results from ImageJ
  rast = raster(paste0("./data/Canny_Edges_Images/Validation/",kiln_name, "/", kiln_name,".tif"))
  imageJ.output = read.csv(paste0("./imageJ_output/", kiln_name, "_imageJ_results.csv"))
  centroid  = xyFromCell(rast,  cellFromRowCol(rast,30,30), spatial = T)
  
  if (ncol(imageJ.output)<5){
    
    ID = NA
    X..pixels.= 1
    Y..pixels. = 1
    Radius..pixels. = 10
    Score = NA 
    nCircles =NA
    Resolution = NA
    Frame..slice... = NA
    Method = "No kilns detected"
    
    imageJ.output = cbind(imageJ.output, ID, X..pixels., Y..pixels., Radius..pixels., Score ,nCircles ,Resolution ,Frame..slice..., Method)
    imageJ.output$r_units = imageJ.output$Radius..pixels. * res(rast)[1] # ft in this case; this only works if the measurements are in linear units (e.g. meters, feet)
    
    
  } else {
    # replace zeros with NAs in chunk column
    imageJ.output$chunk[imageJ.output$chunk == 0] = NA
    imageJ.output = imageJ.output %>% tidyr::fill(chunk, .direction = c("down"))
    imageJ.output$r_units = imageJ.output$Radius..pixels. * res(rast)[1] # ft in this case; this only works if the measurements are in linear units (e.g. meters, feet)
    imageJ.output = imageJ.output[imageJ.output$ID != 0, ]
    borderline.radius = subset(imageJ.output, Radius..pixels. <= 6)
    imageJ.output = imageJ.output[imageJ.output$Radius..pixels. > 6, ]
    borderline.score = subset(borderline.radius, Score >= .55)
    imageJ.output  = rbind(imageJ.output, borderline.score )
  }
  if (nrow(imageJ.output)<1){
    imageJ.output[1,] = imageJ.output[1,]
    imageJ.output$chunk = kiln_name
    imageJ.output$ID = NA
    imageJ.output$X..pixels.= 1
    imageJ.output$Y..pixels. = 1
    imageJ.output$Radius..pixels. = 10
    imageJ.output$ Score = NA 
    imageJ.output$ nCircles =NA
    imageJ.output$Resolution = NA
    imageJ.output$Frame..slice... = NA
    imageJ.output$Method = "No kilns detected"
    
    imageJ.output$r_units = imageJ.output$Radius..pixels. * res(rast)[1] # ft in this case; this only works if the measurements are in linear units (e.g. meters, feet)
  }
  
  chunk.list = intersect(imageJ.output$chunk, chunk.list.all)
  
  num.out = length(chunk.list)
  num.seq = c(1:num.out)
  imageJ.output.chunk = subset(imageJ.output, chunk == chunk.list[1])
  locations = cellFromRowCol(rast[[match(chunk.list[1],chunk.list.all)]], imageJ.output.chunk$Y..pixels., imageJ.output.chunk$X..pixels.)
  kiln.points = xyFromCell(rast[[match(chunk.list[1],chunk.list.all)]], locations, spatial =T)
  
  kiln.data = SpatialPointsDataFrame(kiln.points@coords, imageJ.output, bbox = kiln.points@bbox)
  kiln.data$X_coord = kiln.points@coords[,1]
  kiln.data$Y_coord = kiln.points@coords[,2]
  
  
  if (nrow(kiln.data)>1){
    dist.from.centroid = pointDistance(centroid, kiln.data, lonlat = F)
    keep.index = which.min(dist.from.centroid)
    
    kiln.data = kiln.data[keep.index,]}

  
  # test 1: Is the circular object slightly more prominent than the surrounding area?
  pb = invisible(txtProgressBar(min = 0, max = length(kiln.data), initial = 0, style =3)) # create progress bar
  for (i in 1:length(kiln.data)){
    kiln = kiln.data[i,]
    buffer.big = buffer(kiln, width =  kiln$r_units * 2)
    buffer.small= buffer(kiln, width = kiln$r_units * .8)
    buffer.center = buffer(kiln, width = kiln$r_units * .3)
    buffer.crop = raster::crop(rast, buffer.big)
    buffer.mask.outside = raster::mask(buffer.crop, buffer.small, inverse  = T)
    buffer.mask.inside = raster::mask(buffer.crop, buffer.small, inverse  = F)
    buffer.mask.donut = raster::mask(buffer.mask.inside, buffer.center, inverse  = T)
    outside.elev = cellStats(buffer.mask.outside,mean) + .1 #.4? half a meter adjustment...need some quantitative reasoning for this. 
    kiln.elev = cellStats(buffer.mask.donut,mean)
    #kiln.elev = cellStats(buffer.mask.inside,mean)
    kiln.data$kiln.elev[i] = kiln.elev 
    kiln.data$outside.elev[i] = outside.elev
    if(kiln.elev > outside.elev) {
      kiln.data$kiln_1[i] = 1
    } else {  kiln.data$kiln_1[i] = 0}
    
    # test 2: Is the circular object surrounded by a trench?
    buffer.trench = buffer(kiln, width = kiln$r_units * 1.8)
    buffer.trench.crop = raster::crop(rast, buffer.trench)
    elev.stats = quantile(buffer.crop, probs = seq(0, 1, 0.1))
    ring = Which(buffer.crop <= elev.stats[9])
    if (cellStats(ring, min) == cellStats(ring, max)){
      kiln.data$kiln_2[i] = 0
    } else {
    contour.ring = rasterToContour(ring, nlevels =1)
    #plot(contour.ring)
    
    points = SpatialPoints(cbind(c(contour.ring@bbox[1,1],contour.ring@bbox[1,2]), c(kiln@coords[2], kiln@coords[2])))
    intersect.line <- as(points,"SpatialLines")
    intersection.points = gIntersection(intersect.line, contour.ring, byid=T)
    #plot(intersection.points, add =T)
    if (is.null(intersection.points) == TRUE) {intersection.points = 0}
    
    if(length(intersection.points) >= 2) { # maybe this should 4? two edges of the trench?
      kiln.data$kiln_2[i] = 1
    } else {  kiln.data$kiln_2[i] = 0}
    
    }
    
    
    if(kiln.data$kiln_1[i] == 1 & kiln.data$kiln_2[i] == 1) {
      kiln.data$kiln[i] = 1
    } else {kiln.data$kiln[i] = 0}
    
    setTxtProgressBar(pb,i)
  }
  
  kiln.data$Site = as.character(kiln_name)
  kiln.data = remove.duplicates(kiln.data)
  all.possible.kilns = kiln.data
  if (all.possible.kilns$Method ==  "No kilns detected"){
    all.possible.kilns$kiln_1 = 0
    all.possible.kilns$kiln_2 = 0
    all.possible.kilns$kiln = 0
  }
  tar.kilns = subset(kiln.data, kiln.data$kiln== 1)
  if (index == 1) {
    validation.dataset = all.possible.kilns
  } else {
    validation.dataset = rbind(validation.dataset,all.possible.kilns)
  }
  #plot(rast)
  #plot(known.tar.kilns, pch = 1, add = T)
  #plot(all.possible.kilns, add = T)
  #plot(tar.kilns, col ="red", add = T)
  
  #save shapefile outputs
  #shapefile(tar.kilns, paste0("./output/",lidar.name,"_tar_kilns.shp"), overwrite = T)
  # shapefile(all.possible.kilns, paste0("./output/",lidar.name,"_all_possible_tar_kilns.shp"), overwrite = T)
  #plot(rast, main = paste(kiln_name))
  #plot(tar.kilns, add = T)
  
  
  print(paste0(kiln_name, " was successfully completed."))
  #do.call(file.remove, list(list.files("./data/Canny_Edges_Images/", full.names = TRUE)))
  
}


rdm.val.data = validation.dataset@data
View(rdm.val.data)

rdm.circle.dect.rate = sum(rdm.val.data$Frame..slice..., na.rm = T) / length(rdm.val.data$Frame..slice...)
rdm.circle.dect.rate

rdm.prominence.rate = sum(rdm.val.data$kiln_1, na.rm = T) / length(rdm.val.data$Frame..slice...)
rdm.prominence.rate

rdm.trench.rate = sum(rdm.val.data$kiln_2, na.rm = T) / length(rdm.val.data$Frame..slice...)
rdm.trench.rate

rdm.kiln.id.rate = sum(rdm.val.data$kiln, na.rm = T) / length(rdm.val.data$Frame..slice...)
rdm.kiln.id.rate

# Metrics for Balanced Accuracy #1
TN = (1 - rdm.kiln.id.rate) * length(rdm.val.kilns)
FP = (rdm.kiln.id.rate) * length(rdm.val.kilns)


# Accuracy Metrics
Sensitivity = TP/(TP+FN)
Specificity = TN/(FP+TN)

Balanced.Accuracy = (Sensitivity + Specificity) / 2

Balanced.Accuracy