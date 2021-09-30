# The Tar Kiln Feature Detection workflow (TKFD)
# September 28, 2021
# Grant Snitker - USFS Southern Research Station / ORISE
# grant.snitker@uga.edu

# Manuscript submitted to the Journal of Archaeological Science
# Please do not cite/distribute without author's permission

# Preamble

# load all libraries
rm(list = ls())
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

# set options and number of cores for parallel processing
options(warn=-1)
num.cores = detectCores()-1


##### Part 1: Load data and conduct canny edges detector
# Load lidar DEM file names 
LIDAR.list = list.files("./data/Lidar_data/") %>% file_path_sans_ext(.) %>% noquote(.)
#LIDAR.list

##### Part 2: Define function for running TKFD
l = 1
### Part 2.1: Load lidar DEM data, split into computational chunks, and run canny edges detector on raster data
run.paralell = function(l){
  lidar.name = LIDAR.list[l]
rast = raster(paste0("./data/Lidar_data/", lidar.name, ".tif")) # load lidar DEM
rast[is.na(rast[])] <- -99 # replace NA values with -99 to facilitate the canny edges detector
grid <- as(st_make_grid(rast, square = T, n =c(2,2)), 'Spatial') # create grid for split into computational chunks
split.rasters  = list() # create list to populate with split rasters

pb = invisible(txtProgressBar(min = 0, max = length(grid), initial = 0, style =3)) # create progress bar
for (i in 1:length(grid)){
   skip_to_next <- FALSE
  grid.select= raster::crop(rast, grid[i]) # crop to chunks
  im = as.cimg(grid.select) # convert raster to cimg format for canny edges
  tryCatch(px <- cannyEdges(im, alpha=.5), error = function(e) { skip_to_next <<- TRUE}) # conduct canny edges detector
  if(skip_to_next) { next }     
  px.im = as.cimg(px) # convert result to cimg format
  file = file.exists(paste0("./data/Canny_Edges_Images/",lidar.name,"/"))
  if (file == FALSE){
  dir.create(paste0("./data/Canny_Edges_Images/",lidar.name,"/"))}
  save.image(px.im, file = paste0("./data/Canny_Edges_Images/",lidar.name,"/",lidar.name,"_",(i+1000),".jpg")) # save results in file folder location so FIJI can access it
  split.rasters = append(split.rasters, grid.select)
  setTxtProgressBar(pb,i)
}

### Part 2.2: Perform Hough Circle Transformation (HCT) in FIJI 

# Hough Circle Transformation (HCT) is performed using the UCB Vision Sciences Plugin available in FIJI. For information on activating this plugin, see https://imagej.net/plugins/hough-circle-transform.

# The following parameters for Hough Circle Transformation (HCT) plugin are used in TKFD. These are adjustable by the user for other appliactions. Any adjustments to tehse values must me made in the 'FIJI_macro_TKFD.txt' macro script included in this GitHub repository.

# # Advanced mode
# # min radius: 8
# # max radius 30
# # radius search increment:1
# # default max number to be found: 65535
# # Hough score threshold: 0.5
# # transform resolution: 1000
# # clear neighbors radius ratio: 1.0

# To run FIJI headless through R, the links to paths and the system commands below must be adjusted to the specific operating system and file locations specific to each user.
# See instructions at https://stackoverflow.com/questions/28770970/running-fiji-imagej-macro-from-terminal for adapting system commands to run FIJI headless for Windows, Mac, and Linux systems.
# The code below is designed to run on a Mac OSx.

 system(paste("/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx -macro", paste(getwd(), ("/FIJI_macro_TKFD.txt"), sep = ""), paste(getwd(), ("/data/Canny_Edges_Images/"), (lidar.name),("/"),("*"), getwd(), ("/FIJI_output/"), lidar.name, sep = "")))


# Part 2.3: Process FIJI output
# get outputs from FIJI output folder
chunk.list.all = list.files(paste0("./data/Canny_Edges_Images/",lidar.name, "/")) %>% file_path_sans_ext(.) %>% noquote(.) # list results
FIJI.output = read.csv(paste0("./FIJI_output/", lidar.name, "_FIJI_results.csv")) # load the results from FIJI

if (ncol(FIJI.output)<5){out = paste0("No kilns detected in ", lidar.name)} else {
# replace zeros with NAs in chunk column
FIJI.output$chunk[FIJI.output$chunk == 0] = NA
FIJI.output = FIJI.output %>% tidyr::fill(chunk, .direction = c("down")) # fill out tabel of results
FIJI.output$r_units = FIJI.output$Radius..pixels. * res(rast)[1] # convert radius to units, ft in this case; this only works if the measurements are in linear units (e.g. meters, feet)
FIJI.output = FIJI.output[FIJI.output$ID != 0, ]
FIJI.output = FIJI.output[FIJI.output$Radius..pixels. > 6, ] # exclude outputs that have too small of a diameter
borderline.score = subset(FIJI.output, Score >= .55) # exclude outputs that are not very circular
FIJI.output  = rbind(FIJI.output, borderline.score ) # create final tabel of tar kilns

# assign spatial information to the HCT results from FIJI
if (nrow(FIJI.output)<1){out = paste0("No kilns detected in ", lidar.name)} else {
chunk.list = intersect(FIJI.output$chunk, chunk.list.all) # get the chunk number for each observation
num.out = length(chunk.list)
num.seq = c(1:num.out)
FIJI.output.chunk = subset(FIJI.output, chunk == chunk.list[num.seq[1]]) # select chunk
locations = cellFromRowCol(split.rasters[[match(chunk.list[num.seq[1]],chunk.list.all)]], FIJI.output.chunk$Y..pixels., FIJI.output.chunk$X..pixels.) # extract coordinates of tar kiln from chunk
kiln.points = xyFromCell(split.rasters[[match(chunk.list[num.seq[1]],chunk.list.all)]], locations, spatial =T) # assign spatial locations to kilns

# convert locations to a spatial points dataframe - This step repeat the previous code chuck, but for all observations beyond the first row. 
if (length(num.seq) > 1){
  for (i in num.seq[2]:length(num.seq)){
    FIJI.output.chunk = subset(FIJI.output, chunk == chunk.list[i])
    locations = cellFromRowCol(split.rasters[[match(chunk.list[i],chunk.list.all)]], FIJI.output.chunk$Y..pixels., FIJI.output.chunk$X..pixels.)
    kiln.points.add = xyFromCell(split.rasters[[match(chunk.list[i],chunk.list.all)]], locations, spatial =T)
    kiln.points = rbind(kiln.points, kiln.points.add)}
  
  kiln.data = SpatialPointsDataFrame(kiln.points@coords, FIJI.output,
                                     proj4string = kiln.points@proj4string, bbox = kiln.points@bbox)
  
  kiln.data$X_coord = kiln.points@coords[,1]
  kiln.data$Y_coord = kiln.points@coords[,2]
} else { # convert to spatial points dataframe
  kiln.data = SpatialPointsDataFrame(kiln.points@coords, FIJI.output, bbox = kiln.points@bbox)
 
  kiln.data$X_coord = kiln.points@coords[,1]
  kiln.data$Y_coord = kiln.points@coords[,2]}

# Part 2.4: Criteria tally to determine what is a tar kiln
# test 1: Is the circular object slightly more prominent than the surrounding area?
pb = invisible(txtProgressBar(min = 0, max = length(kiln.data), initial = 0, style =3)) # create progress bar
for (i in 1:length(kiln.data)){
kiln = kiln.data[i,]
buffer.big = buffer(kiln, width =  kiln$r_units * 2) # create a large buffer around possible tar kiln
buffer.small= buffer(kiln, width = kiln$r_units * .8) # create a smaller buffer around possible tar kiln
buffer.center = buffer(kiln, width = kiln$r_units * .3) # create a central buffer within possible tar kiln
buffer.crop = raster::crop(rast, buffer.big)
buffer.mask.outside = raster::mask(buffer.crop, buffer.small, inverse  = T)
buffer.mask.inside = raster::mask(buffer.crop, buffer.small, inverse  = F)
buffer.mask.donut = raster::mask(buffer.mask.inside, buffer.center, inverse  = T) # create a elevation profile in the shape of a 'donut'. This isolates the tar kiln mound
outside.elev = cellStats(buffer.mask.outside,mean) + .1 # outside of kiln elevation--0.1 a meter adjustment based on observations of tar kilns in the field and in validation dataset
kiln.elev = cellStats(buffer.mask.inside,mean) # the elevation of the kiln mound
kiln.data$kiln.elev[i] = kiln.elev 
kiln.data$outside.elev[i] = outside.elev
if(kiln.elev > outside.elev) { # compare outside and kiln elevations and tally criteria result
  kiln.data$kiln_1[i] = "yes"
} else {  kiln.data$kiln_1[i] = "no"}

# test 2: Is the circular object surrounded by a trench?
buffer.trench = buffer(kiln, width = kiln$r_units * 1.8) # create a wide buffer around the possible tar kiln
buffer.trench.crop = raster::crop(rast, buffer.trench)
elev.stats = quantile(buffer.crop, probs = seq(0, 1, 0.1)) # create elevation quantiles
ring = Which(buffer.crop <= elev.stats[9]) 
contour.ring = rasterToContour(ring, nlevels =1) # convert qunatile values ot contour lines

points = SpatialPoints(cbind(c(contour.ring@bbox[1,1],contour.ring@bbox[1,2]), c(kiln@coords[2], kiln@coords[2]))) # place spatial points on the contour lines 
intersect.line <- as(points,"SpatialLines") # darw line between points
intersection.points = gIntersection(intersect.line, contour.ring, byid=T) # determine number of intersections between points and lines
if (is.null(intersection.points) == TRUE) {intersection.points = 0}

if(length(intersection.points) >= 4) {# assign trench presenec based on number of intersections. 4 intersections suggests the presence of a low trench surrounding the kiln.
  kiln.data$kiln_2[i] = "yes"
} else {  kiln.data$kiln_2[i] = "no"}

if(kiln.data$kiln_1[i] == "yes" & kiln.data$kiln_2[i] == "yes") {
  kiln.data$kiln[i] = "yes"
} else {kiln.data$kiln[i] = "no"}

setTxtProgressBar(pb,i)
}

kiln.data$LIDAR = as.character(lidar.name) #record the lidar DEM name for reference
kiln.data = remove.duplicates(kiln.data) # remove any possible duplicates
all.possible.kilns = kiln.data # dataset of all possible kilns 
tar.kilns = subset(kiln.data, kiln.data$kiln== "yes") # dataset of tar kilns as determined by the criteria tally

# Part 2.5: Export results

#export results as shapefile 
shapefile(tar.kilns, paste0("./output/",lidar.name,"_tar_kilns.shp"), overwrite = T)
shapefile(all.possible.kilns, paste0("./output/",lidar.name,"_all_possible_tar_kilns.shp"), overwrite = T)
out = paste0(lidar.name, " was successfully completed.") # print completion statement
} 
}
print(out)
}

##### Part 3: Run TKFD
# run in parallel
cl = makeCluster(num.cores, type = "PSOCK")
clusterExport(cl, varlist = ls())
invisible(clusterEvalQ(cl, c(library(imager),library(rgeos),library(sp),
                       library(raster),library(tools),library(rgdal),
                       library(sf),library(tidyverse),library(reshape2),
                       library(geosphere),library(parallel),library(doParallel)))) # all required libraries
ptm <- proc.time() # start time
output = parLapply(cl,1:length(LIDAR.list),run.paralell)
stopCluster(cl)