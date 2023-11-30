# Credits ---------------------------

# Original routine by Vitor Cavalcante

# Edited by
# Luara Tourinho & Julia Niemeyer

# Last update: 28 nov 2023

# This script is an example for generating pseudoabsence using biomod2 


# Required packages

library(biomod2)
library(raster)
library(dplyr)


# intersect_mask <- function(x){
#   values_x <- getValues(x)
#   inter_x <- values_x %*% rep(1,nlayers(x))
#   mask <- setValues(subset(x,1),values = (inter_x>0))
#   return(mask)}


# Opening occurences ------------------------------------------------------

# Read file of occurrences
file = "./01_occs.csv" ##enter the name of your table

# minimum occurrence records to run analysis
n_min <- 15


# Reading files -----------------------------------------------------------

splist  <- read.csv("./01_n_records.csv", stringsAsFactors = FALSE)
to_filter  <- splist %>%
subset(n_thin_5 >= 15) %>%
  pull(species)

sp_occs <- read.table(file, header=TRUE, sep=",")
sp_0 <- sp_occs %>%
  subset(species %in% to_filter)

sp_names <- unique(sp_0$species)


# Loop to create pseudoabsences -------------------------------------------


for (a in 1:length(sp_names)){
# message("starting the analysis for ", paste0(sp_names[a]))

sp <- sp_0 %>%
  filter(species == paste0(sp_names[a])) %>%
  select(species, lon, lat)

message("starting the analysis for ", paste0(sp_names[a]))

if (nrow(sp) < n_min){ ##Will not analyze species with less than 15 occurences
  print('species has less than 15 records and will not be analyzed')
  next
}

#Get only lat and long
My_target_species <- sp[,2:3]


# Reading environmental variables -----------------------------------------

# Read your environmental raster selected
raster_files <-
  list.files(paste0("./outputs/", sp_names[a], "/Pres_env_crop"),
             full.names = T,
             'tif$|bil$')

head(raster_files)

environment <- stack(raster_files)

occurrence.resp <- rep(1, length(My_target_species$lon))

skip_to_next <<- FALSE
PA.list <- list(50000, 10000) 
# it will run for 50km first. For species with restricted distribution, 50km will
# return an error and the function will run for 10km

# Pseudoabsences for GLM and SVM

for (i in 1:length(PA.list)) {

  tryCatch(Mymodel <- BIOMOD_FormatingData(
  resp.var = occurrence.resp,
  expl.var = environment,
  resp.xy = My_target_species,
  resp.name = "Occurrence",
  PA.nb.rep = 1,
  PA.nb.absences = length(sp$species)*10,
  PA.strategy = "disk",
  PA.dist.min = PA.list[[i]],
  PA.dist.max = 20000000,
  na.rm = TRUE) , error = function(e) {skip_to_next <<- TRUE})

  if(skip_to_next == TRUE) {
    print(paste0('Species has restricted distribution, trying 10km'))
    next
  } else {

    print(paste0('using minimum distance of ', PA.list[[i]]))
    break
  }}
  

# Pseudoabsences for RF

  for (i in 1:length(PA.list)) {

    tryCatch(Mymodel2 <- BIOMOD_FormatingData(
  resp.var = occurrence.resp,
  expl.var = environment,
  resp.xy = My_target_species,
  resp.name = "Occurrence",
  PA.nb.rep = 1,
  PA.nb.absences = length(sp$species), ##Isso tem que entrar no loop
  PA.strategy = "disk",
  PA.dist.min = PA.list[[i]],
  PA.dist.max = 20000000,
  na.rm = TRUE) , error = function(e) {skip_to_next <<- TRUE})

    if(skip_to_next == TRUE) {
      print(paste0('Species has restricted distribution, trying 10km'))
      next
    } else {

      print(paste0('using minimum distance of ', PA.list[[i]]))
      break
    }}


# function to get PA dataset

get_PAtab <- function(bfd){
  dplyr::bind_cols(
    x = bfd@coord[, 1],
    y = bfd@coord[, 2],
    status = bfd@data.species,
    bfd@PA
  )
}

# function to get background mask

get_mask <- function(bfd){
  bfd@data.mask
}

(pres.xy <- get_PAtab(Mymodel) %>%
    filter(status == 1) %>%
    select(x, y))

(pres.xy2 <- get_PAtab(Mymodel2) %>%
    filter(status == 1) %>%
    select(x, y))

# get the coordinates of pseudoabsences
# all repetition of pseudoabsences sampling merged

(pa.all.xy <- get_PAtab(Mymodel) %>%
    filter(is.na(status)) %>%
    select(x, y)) %>%
  distinct()

(pa.all.xy2 <- get_PAtab(Mymodel2) %>%
    filter(is.na(status)) %>%
    select(x, y)) %>%
  distinct()

#plot(environment[[1]])
#points(pa.all.xy, pch = 18)
#points(pres.xy, pch = 20, col= "red")

write.csv(pa.all.xy,
          paste0("./outputs/", sp_names[a], "/pseudoabs1.csv"),
          row.names = F) 
write.csv(pa.all.xy2,
          paste0("./outputs/", sp_names[a], "/pseudoabs2.csv"),
          row.names = F) 


pseudoabs <- pa.all.xy
#head(pseudoabs)
#dim(pseudoabs)

pres = sp
# Replace using your species name in "Genus_epithet"
pres$`species` <- sub(pattern = paste0(sp_names[a]), replacement = "1", x = pres$`species`)
#tail(pres)
pseudo_0 <- rep(0,nrow(pseudoabs))
#pseudoabs$species <- NA
#pseudoabs$species <- 0
pseudoabs$species <- pseudo_0
pseudoabs <- pseudoabs[,c(3,1,2)]
names(pseudoabs) <-c("species","lon","lat")
pres_pseudo_table <- rbind(pres,pseudoabs)
#head(pres_pseudo_table)
#tail(pres_pseudo_table)
#dim(pres_pseudo_table)
names(pres_pseudo_table) <-c("pa","lon","lat")

pseudoabs2 <- pa.all.xy2
#head(pseudoabs2)
#dim(pseudoabs)

pres2 = sp
# Replace using your species name in "Genus_epithet"
pres2$`species` <- sub(pattern = paste0(sp_names[a]), replacement = "1", x = pres2$`species`)
#tail(pres)
pseudo_0 <- rep(0,nrow(pseudoabs2))
pseudoabs2$species <- pseudo_0
#pseudoabs2$species <- NA
#pseudoabs2$species <- 0
pseudoabs2 <- pseudoabs2[,c(3,1,2)]
names(pseudoabs2) <-c("species","lon","lat")
pres_pseudo_table2 <- rbind(pres2,pseudoabs2)
#head(pres_pseudo_table2)
#tail(pres_pseudo_table)
#dim(pres_pseudo_table)
names(pres_pseudo_table2) <-c("pa","lon","lat")


write.csv(
  pres_pseudo_table,
  paste0("./outputs/", sp_names[a], "/pres_pseudoabs.csv"),
  row.names = F)
write.csv(
  pres_pseudo_table2,
  paste0("./outputs/", sp_names[a], "/pres_pseudoabs2.csv"),
  row.names = F)

}

# END

