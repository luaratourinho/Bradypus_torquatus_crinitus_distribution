
# Credits ---------------------------

# Script created by
# Aiello-Lammens et al. (https://cran.r-project.org/web/packages/spThin/spThin.pdf)
# Aiello-Lammens et al. 2015 (doi: 10.1111/ecog.01132)

# Edited by
# Luara Tourinho (https://github.com/luaratourinho)
# and 
# Bruno Carvalho (https://github.com/brunomc-eco)

# Date: 27 nov 2023



# Geographical cleaning using spThin package ------------------------------



# Required packages

library(spThin) # Aiello-Lammens et al. 2015
library(tidyverse)
library(data.table)
library(dplyr)


# # loading clean occs and getting clean species list

clean_df <- read.csv("./species.csv",
                     stringsAsFactors = FALSE)

spp <- read.csv("./species.csv",
                           stringsAsFactors = FALSE) %>%
  pull(species)

spp <- unique(spp)


# Run to distance you are interested in

# thinning records by 5 km
thin_5 <- list()
for(i in 1:length(spp)){
  df <- clean_df %>%
    filter(species %in% spp[i])
  thinned <- thin(df,
                  lat.col = "lat",
                  long.col = "lon",
                  spec.col = "species",
                  thin.par = 1, # distance in km
                  reps = 1,
                  locs.thinned.list.return = TRUE,
                  write.files = FALSE,
                  write.log.file = FALSE)
  thin_5[[i]] <- data.frame(species = rep(spp[i], nrow(thinned[[1]])),
                            lon = thinned[[1]]$Longitude,
                            lat = thinned[[1]]$Latitude)
}
clean_df_thin_5 <- rbindlist(thin_5)


# Check thinned records
ggplot() +
 borders("world", colour="gray50", fill="gray50") +
 geom_point(data = clean_df_thin_5, aes(x = lon, y = lat),
            colour = "blue", size = 1.5) +
 # geom_point(data = clean_df_thin_10, aes(x = lon, y = lat),
 #            colour = "red", size = 1.0) +
 coord_sf(xlim = c(-160, -28), ylim = c(-60,90)) +
 theme_bw()


# counting records by species
n_5 <- clean_df_thin_5 %>%
  group_by(species) %>%
  summarize(n_thin_5 = n())


# writing outputs
write_csv(n_5, path = "./01_n_records.csv")
write_csv(clean_df_thin_5, path = "./01_occs.csv")

