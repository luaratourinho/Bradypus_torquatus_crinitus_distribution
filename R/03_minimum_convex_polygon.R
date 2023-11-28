# Credits -----------------------------------------------------------------

# Code created by Luara Tourinho (https://github.com/luaratourinho)

# Edited by
# Julia Niemeyer & Luara Tourinho

# Last update: 27 nov 2023


# Required packages

library(maps)
library(mapdata)
library(rgdal)
library(raster)
library(rgeos)
library(dplyr)


# Creating some objects to use later --------------------------------------

dir.create("./outputs")

# minimum occurrence records to run analysis
n_min <- 15

# projection WGS84
crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# South America Albers Equal Area Conic
crs.albers <-
  CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60
                  +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

# Lambert Azimuthal Equal Area (ideal for Pacific Ocean islands)
#crs.azim <-
# CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84
#               +datum=WGS84 +units=m +no_defs")


# Reading files -----------------------------------------------------------

file = "./01_occs.csv" ##enter the name of your table

sp <- read.table(file, header=TRUE, sep=",")

target_species <- unique(sp$species)

sp_names <- target_species


## Enter name of folders to read environmental layers
# I added the chosen file on env_sel directory

# Present

pres_folder = "./env_sel"
pres_files <- list.files(pres_folder, full.names = T, 'tif$|bil$')
##standardize names of the variables as in your model 
#in the order of appearance in head(pres_files)
head(pres_files)

envi <- stack(pres_files)


# Building a minimum convex polygon ---------------------------------------

for (a in 1:length(sp_names)){
  
  message("starting the analysis for ", paste0(sp_names[a]))
  
  sp.n = sp_names[[a]]
  
  # Read your table of occurrence records
  
  occurrence_records <- sp %>%
    filter(species == paste0(sp_names[a])) %>%
    select(species, lon, lat)
  
  # Number of occurrences to perform pseudoabsence sampling
 
  if (nrow(occurrence_records) < n_min){ ##Will not analyze species with less than n_min occurrences
    print('species has less than 15 records and will not be analyzed')
    target_dir = paste0("./ENM/outputs/", sp_names[a], "/")
    dir.create( target_dir )
    write(format('species has less than 15 records and will not be analyzed'), file=paste(target_dir, 'STOP.txt', sep=""))
    next
  }
  
  # Check the column names for your coordinates and place below within c("","")
  coordinates(occurrence_records) <- c("lon", "lat")
  
  # Define the original projection of your occurrence records
  proj4string(occurrence_records) <- crs.wgs84
  
  # Reproject to an equal area projection
  occurrence_records_albers <-
    spTransform(occurrence_records, crs.albers)
  
  # Minimum convex polygon (mpc) function
  mpc <- gConvexHull(occurrence_records_albers)
  
  
  # Adding a buffer to the mpc -----------------------------------------------
  
  # You can choose a fixed value for a buffer, for example
  # If you want to 500km of ratio, use b = 500000, because b is in meters
  # If you want to a proportional value as suggested by Barve et al. 2011
  # (doi:10.1016/j.ecolmodel.2011.02.011), use the function gArea.
  # The function gArea is in m2, so multiplying by 0.000001 you can get a area in km2
  # For example, if you want to draw a buffer of 20%, in km, around of the mpc
  # 20% of the polygon area is 2e-07
  # If you want to draw a buffer of 2° it is ~111km, i.e. b <- (111000)
  
  # Function for creating different size of buffer
  # This is useful when your species has a very large distribution, then
  # 20% or even 10% of its distribution area would be bigger than the world
  # in this case, the R return an error about infinite values.
  
  b20 <- (gArea(mpc)*2e-07) ##20%
  b10 <- (gArea(mpc)*1e-07) ##10%
  b5 <- (gArea(mpc)*5e-8) ##5%
  b0 <- (gArea(mpc)*1e-09) ## ~ no buffer
  b.list <- list(b20, b10, b5, b0)
  
  for (i in 1:length(b.list)) {
    
    skip_to_next <<- FALSE
    
    buffer_mpc <- gBuffer(mpc, width = b.list[[i]])
    
    tryCatch(polygon_wgs <- spTransform(buffer_mpc, crs.wgs84),
             error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next == TRUE) {
      print(paste0('buffer of place ',i, ' got error, skipping to next'))
      next
    } else {
      
      print(paste0('using buffer of place ', i, ' in the list: ', b.list[[i]], ' m2'))
      # polygon_wgs <- spTransform(buffer_mpc, crs.wgs84)
      break
    }}
  
  
  
  # Applying to the variables -----------------------------------------------
  
  # If you want to clip a worldclim environmental raster
  # (https://www.worldclim.org/), which is in WGS84
  
  # Reproject your polygon to the same projection of your environmental data
  polygon_wgs <- spTransform(buffer_mpc, crs.wgs84)
  # plot(polygon_wgs)
  
  # Cut your study area for the same extension and shape of your mpc
  present_ly <- crop(envi, polygon_wgs)
  present_ly2 <- mask(present_ly, polygon_wgs)
  
  # Plot the results
  # plot(present_ly2[[1]])
  # plot(occurrence_records, add = T)
  
  
  # Save your layers --------------------------------------------------------
  
  dir.create(paste0("./outputs/", sp_names[a]))
  
  
  dir.create(paste0("./outputs/", sp_names[a], "/Pres_env_crop/"))
  
  writeRaster(present_ly2,
    filename = paste0("./outputs/", sp_names[a], "/Pres_env_crop/",
      names(present_ly2)), bylayer = TRUE, format = "GTiff", overwrite = T)
  
}
