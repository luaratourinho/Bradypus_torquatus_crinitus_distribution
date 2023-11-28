
# Credits ---------------------------


# https://github.com/Model-R/modleR

# Elaborated by
# Diogo S. B. Rocha (https://github.com/diogosbr)
# Bruno M. Carvalho (https://github.com/###########)


# Edited by
# Luara Tourinho

# Date: 27 nov 2023


# Install packages

# remotes::install_github("Model-R/modleR", build = TRUE)
# remotes::install_github("marlonecobos/kuenm")
# remotes::install_github("mrmaxent/maxnet")

# Required packages

library(dplyr)
library(raster)
library(progress)
library(modleR)
library(foreach)
library(rgeos)
library(doParallel)
library(tidyverse)

# ENM using modleR --------------------------------------------------------

# Reading data ------------------------------------------------------------

# reading species data, only names

# reading occurrence table of all species
# in case that your occurrence table have more species than you want to project
clean_df <- read.csv("./01_occs.csv",
                     stringsAsFactors = FALSE)

target_species <- clean_df %>%
  group_by(species) %>%
  summarize(n = n())

target_species <- target_species %>%
  pull(species)

# reading climatic data from prensent conditions
wc <- list.files("./env_sel/",
                 pattern = "tif$",
                 full.names = TRUE) %>%
  stack()


# Parallel ----------------------------------------------------------------

# loops for ENM
# number of species if you have enough nodes
registerDoParallel(cores = 2) 


# Projections used in mpc -------------------------------------------------


# projections
# needed to perform the minimum convex polygon (mcp) below
crs.wgs84 <-
  CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
crs.albers <-
  CRS(
    "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs"
  ) # projected, South America Albers Equal Area Conic


# setup_sdmdata step ---------------------------------------------------------

foreach(
  sp = target_species,
  .inorder = FALSE,
  .packages = c("modleR", "sp", "raster", "rgeos", "dplyr")
) %dopar% {
  species_df <-
    clean_df[clean_df$species == sp, ] # getting occurrences
  
  # Creating the mpc buffer 
  
  coords <- species_df[, 2:3]
  coordinates(coords) <- c("lon", "lat")
  proj4string(coords) <- crs.wgs84
  coords <- spTransform(coords, crs.albers)
  mcp <- gConvexHull(coords)
  mcp_buf <- mcp %>%
    gBuffer(width = gArea(mcp) * 2e-07) %>% 
    spTransform(crs.wgs84) %>%
    SpatialPolygonsDataFrame(data = data.frame("val" = 1, row.names = "buffer"))
  
  # choosing the type of partition depending on the number of records
  partition_type <- ifelse(nrow(species_df) > 50, "crossvalidation", "bootstrap")
  
  setup_sdmdata(
    species_name = sp,
    occurrences = species_df,
    predictors = wc,
    models_dir = "./models/",
    # folder to save partitions
    buffer_type = "user",
    buffer_shape = mcp_buf,
    clean_dupl = TRUE,
    clean_nas = TRUE,
    clean_uni = TRUE,
    # remove records falling at the same pixel
    png_sdmdata = TRUE,
    n_back = nrow(species_df) * 10,
    # number of pseudoabsences
    partition_type = partition_type,
    cv_partitions = 10,
    cv_n = 1,
    boot_n = 10, 
    boot_proportion = 0.8)
}

# do_many step ------------------------------------------------------------

# partitions
foreach(sp = target_species,
        .inorder = FALSE,
        .packages = "modleR") %dopar% {
          do_many(
            species_name = sp,
            predictors = wc,
            models_dir = "./models/",
            png_partitions = FALSE,
            write_bin_cut = TRUE,
            dismo_threshold = "spec_sens",
            equalize = TRUE, # equalize presence and pseudoabsences for RF and BRT
            bioclim = TRUE,
            glm = TRUE,
            maxent = TRUE,
            rf = TRUE,
            svmk = TRUE
          )
        }


# final_model step --------------------------------------------------------

#combine partitions into one final model per algorithm
# projections path names

# final_models
foreach(sp = target_species,
        .inorder = FALSE,
        .packages = "modleR") %dopar% {
          final_model(
            species_name = sp,
            models_dir = "./models/",
            which_models = c("raw_mean", "bin_consensus"),
            consensus_level = 0.5,
            # proportion of models in the binary consensus
            proj_dir = "present",
            uncertainty = FALSE,
            png_final = FALSE,
            overwrite = TRUE
          )
          
        }
gc()


# ensemble_model step -----------------------------------------------------

# generate ensemble models, combining final models from all algorithms
# ensemble models

foreach(
  sp = target_species,
  .inorder = TRUE,
  .packages = c("modleR", "raster", "dplyr")
) %dopar% {
  species_df <- clean_df[clean_df$species == sp,]
  ensemble_model(
    species_name = sp,
    algorithms = c("bioclim", "glm", "maxent", "rf", "svmk"),
    models_dir = "./SP/birds/ENM/models/",
    performance_metric = "TSSmax",
    proj_dir = "present",
    which_ensemble = c("weighted_average", "consensus"),
    which_final = c("raw_mean", "bin_consensus"),
    ensemble_dir = "ensemble_all_algo",
    consensus_level = 0.5,
    png_ensemble = FALSE,
    uncertainty = TRUE,
    overwrite = TRUE
  )
}

gc()


# Ensemble among GCMs -----------------------------------------------------




