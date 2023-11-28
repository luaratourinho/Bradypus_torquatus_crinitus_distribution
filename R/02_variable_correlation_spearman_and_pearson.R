
# Credits ---------------------------

# Created by
# Luara Tourinho (https://github.com/luaratourinho)

# Based on: 
# Diogo S. B. Rocha (https://github.com/diogosbr)

# Last update: 27 nov 2023


# Required packages

library(raster)
library(dplyr)
library(corrplot)
library(raster)
library(rgeos)
library(reshape)


# Cropping variables from Worldclim ---------------------------------------

# Reading rasters (biovariables 30s from https://www.worldclim.org/data/worldclim21.html)
# use pattern = '.tif$' or something else if you have multiple files in this folder
raster_files <- list.files("./env/", 
                           full.names = T, 'tif$|bil$')

head(raster_files)

envi <- stack(raster_files)

# Cropping rasters

# Choose your extension
# All America
envi.cut<-crop(envi, c(-60, -25, -50, 40))
plot(envi.cut[[1]])

# Saving rasters
dir.create(paste0("./env_cropped/", "."))
writeRaster(envi.mask, filename='./env_cropped/', format="GTiff", 
            bylayer=TRUE, suffix="names", overwrite=TRUE)



# Select variable with less collinearity ---------------------------------------

# list of files of present variables
present_list <- list.files("./env_cropped/", pattern = "tif$", full.names = T)

# object with present variables
present_ras <- stack(present_list)

# id from cells without NA
mi <- Which(present_ras[[1]], cells = TRUE)

# sampling cells to extract values
sampled <- sample(mi, 5000)

# values of selected cells from rasters of present variables
vals <- present_ras[sampled]


# An alternative using spearman ------------------------------------------------

# selecting variables to exclude with correlation 0.6 or more
# First try low values, as 0.6, if it return few variables, try higher as 0.7
exclude_vars <- caret::findCorrelation(cor(vals, method = 'spearman'), cutoff = 0.6, names = TRUE)
all_table <- as.data.frame(cor(vals, method = 'spearman'))

# selecting variables with lower correlation (<0.6)
pres_vars_sel <- present_ras[[which(!names(present_ras) %in% exclude_vars)]]
pres_vars_sel

# selecting variables with lower correlation (<0.6)
pres_vars_sel_names <- names(present_ras)[!names(present_ras) %in% exclude_vars]
pres_vars_sel_names
# [1] "X_wc2.1_30s_bio_10" "X_wc2.1_30s_bio_13" "X_wc2.1_30s_bio_17" "X_wc2.1_30s_bio_18"
# [5] "X_wc2.1_30s_bio_2"

# Ps.: I choose spearman to this analysis

# An alternative using pearson -------------------------------------------------

exclude_vars2 <- caret::findCorrelation(cor(vals, method = 'pearson'), cutoff = 0.6, names = TRUE)
all_table2 <- as.data.frame(cor(vals, method = 'pearson'))

pres_vars_sel2 <- present_ras[[which(!names(present_ras) %in% exclude_vars2)]]
pres_vars_sel2

pres_vars_sel_names2 <- names(present_ras)[!names(present_ras) %in% exclude_vars2]
pres_vars_sel_names2
# [1] "X_wc2.1_30s_bio_16" "X_wc2.1_30s_bio_17" "X_wc2.1_30s_bio_18" "X_wc2.1_30s_bio_19"
# [5] "X_wc2.1_30s_bio_2"  "X_wc2.1_30s_bio_8"  


# creating directory
dir.create(paste0("./env_sel"))
# Add the chosen variables there

write.csv(all_table, "./env_sel/5_variables_correlation_spearman.csv")
write.csv(all_table2, "./env_sel/5_variables_correlation_pearson.csv")

# Check Worldclim to chose which set of variable make more sense, if by spearman
# or pearson
# https://www.worldclim.org/data/bioclim.html

# END