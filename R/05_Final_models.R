# Credits ---------------------------

# Script created by
# Bruno M. Carvalho (https://github.com/brunomc-eco)
#  & Luara Tourinho (https://github.com/luaratourinho)

# Edited by
# Julia Niemeyer & Luara Tourinho
# Last update: 13 Jul 2022


# Running ENM for multiple species & GCMs

# Required packages

library(parallel)
library(foreach)
library(doParallel)
require(raster)
library(rgeos)
library(rgdal)
library(dismo)
library(rJava)
library(kernlab)
library(randomForest)
library(maptools)
library(SDMTools)
library(dplyr)
library(biomod2)
library(beepr)
library(maxnet)



indexOf <- function(v,findFor) {
  i=0
  for(i2 in v) {i = i + 1
    if (i2==findFor){
      return(i)
    }
  }
  return(0)
}


# Settings ----------------------------------------------------------------

# For GLM, RandomForest and SVM
# Set your model to the biovariables you selected
# It should be exactly like the names of the environmental rasters
model <- pa ~bio_10+bio_13+bio_17+bio_18+bio_2

# partitions (4 = 75% train and 25% test; 5 = 80% and 20%; 10 = 90% and 10%)
k = 10
# number of background
bg.pt = 10000
# threshold ('spec_sens' = max. spec + sens)
t.met = 'spec_sens'
# Minimum TSS value for ensemble
tss.lim = 0.6

cont.maps = T # save continuous maps by algorithm
bin.maps = T # save binary maps by algorithm
ens.maps = T # save ensemble maps by algorithm

# Name of .csv file with species occurrences, columns: 'species', 'lon', 'lat'
file = "./01_occs.csv"

# Minimum occurrence records to run analysis
n_min <- 15

# Reading files -----------------------------------------------------------

splist  <- read.csv("./01_n_records.csv", 
                    stringsAsFactors = FALSE)

to_filter  <- splist %>%
  subset(n_thin_5 >= 15) %>%
  pull(species)

sp_occs <- read.table(file, header=TRUE, sep=",")

sp_0 <- sp_occs %>%
  subset(species %in% to_filter)

sp_0$species <-
  gsub(x = sp_0$species,
       pattern = " ",
       replacement = "_")

sp <- sp_0

sp_names <- unique(sp_0$species)


for (a in 1:length(sp_names)){
  
sp.n = sp_names[[a]]
 message("starting the analysis for ", paste0(sp.n))

presences <- sp  %>%
  filter(species == paste0(sp.n)) %>%
  select(species, lon, lat)

if (nrow(presences) < n_min){ # Will not analyze species with < 15 occurrences
  print('species has less than 15 records and will not be analyzed')
  next
}


# Running for one species
sp.names <- as.character(unique(presences$species))

# Read predictor variables (i.e. present maps cropped by mcp)
raster_files <-
  list.files(paste0("./outputs/", sp.n, "/Pres_env_crop"),
             full.names = T,
             'tif$|bil$')
#head(raster_files)
predictors <- stack(raster_files)


# ENM ---------------------------------------------------------------------


ini = Sys.time()
pb <- txtProgressBar(min = 1, max = length(sp.names)+1, style = 3)

# Number of occurrences to perform pseudoabsence sampling
lim = plyr::count(presences$species)$freq[indexOf(plyr::count(presences$species)$x, sp.n)]

started_time = Sys.time()
cat( format( started_time, "%a %b %d %X %Y"), '-', 'STARTED', '\n')
cat(format(started_time, "%a %b %d %X %Y"),
  '-','Preparing train and test datasets for',
  sp.n,'with ',lim,'lines...','\n')

# Create the directory where results from ENM will be saved
target_dir = paste(paste0("./outputs/", sp.n, "/results/", sep=""))
dir.create(target_dir)

if(file.exists(paste(target_dir, '/STARTED.txt', sep="")))
  stop("You MUST DELETE previous results folder before continue")

write(format( started_time, "%a %b %d %X %Y"), 
      file=paste(target_dir, 'STARTED.txt', sep=""))

# For using different number of pseudoabsence in each algorithm, for example:
## pseudoausencia = n*10
sp.data <- read.csv(paste0("./outputs/", 
                           sp.n,"/pres_pseudoabs.csv"), header=TRUE, sep=',')

## pseudoausencia = pres - to random forest
sp.data2 <- read.csv(paste0("./outputs/", 
                            sp.n,"/pres_pseudoabs2.csv"), header=TRUE, sep=',')


pres <- sp.data[sp.data$pa==1,2:3]
abs <- sp.data[sp.data$pa!=1,2:3]
bg <- randomPoints(predictors, bg.pt)
colnames(bg) <- c("lon", "lat")


## Crossvalidation - partition in 10 groups
set.seed(10)
foldpres <- kfold(pres, 4)
set.seed(10)
foldabs <- kfold(abs, 4)

prestrain <- list()
prestest <- list()
abstrain <- list()
abstest <- list()
for(i in 1:k){
  foldpres <- kfold(pres, 4)
  foldabs <- kfold(abs, 4)
  prestrain[[i]] <- pres[foldpres != 1,]
  prestest[[i]] <- pres[foldpres == 1,]
  abstrain[[i]] <- abs[foldabs != 1,]
  abstest[[i]] <- abs[foldabs == 1,]
}


# For SVM and GLM
train <- list()
pa_train <- list()
predtrain <- list()
testpres <- list()
testabs <- list()
data <- list()
p <- list()

for(i in 1:k){
  train[[i]] <- rbind(prestrain[[i]], abstrain[[i]])
  pa_train[[i]] <- c(rep(1, nrow(prestrain[[i]])), rep(0, nrow(abstrain[[i]])))
  predtrain[[i]] <- raster::extract(predictors, train[[i]])
  predtrain[[i]] <- na.omit(data.frame(cbind(pa=pa_train[[i]], predtrain[[i]])))
  data[[i]] <- predtrain[[i]][,-1]
  p[[i]] <- predtrain[[i]][,1]
}


# For Random Forest
pres2 <- sp.data2[sp.data2$pa==1,2:3]
abs2 <- sp.data2[sp.data2$pa!=1,2:3]
set.seed(10)
foldpres2 <- kfold(pres2, 4)
set.seed(10)
foldabs2 <- kfold(abs2, 4)

prestrain2 <- list()
prestest2 <- list()
abstrain2 <- list()
abstest2 <- list()
for(i in 1:k){
  foldpres2 <- kfold(pres2, 4)
  foldabs2 <- kfold(abs2, 4)
  prestrain2[[i]] <- pres2[foldpres2 != 1,]
  prestest2[[i]] <- pres2[foldpres2 == 1,]
  abstrain2[[i]] <- abs2[foldabs2 != 1,]
  abstest2[[i]] <- abs2[foldabs2 == 1,]
}

train2 <- list()
pa_train2 <- list()
predtrain2 <- list()
testpres2 <- list()
testabs2 <- list()

for(i in 1:k){
  train2[[i]] <- rbind(prestrain2[[i]], abstrain2[[i]])
  pa_train2[[i]] <- c(rep(1, nrow(prestrain2[[i]])), rep(0, nrow(abstrain2[[i]])))
  predtrain2[[i]] <- raster::extract(predictors, train2[[i]])
  predtrain2[[i]] <- data.frame(cbind(pa=pa_train2[[i]], predtrain2[[i]]))
}



# Running models ----------------------------------------------------------

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
     'Running Bioclim model for', sp.n, '...', '\n')
bc <- list()
evbc <- list()
bcTSS <- list()
bcAUC <- list()
bckappa <- list()
bcthres <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
       'Running Bioclim (', i, ') model for', sp.n, '...', '\n')
  bc[[i]] <- bioclim(predictors, prestrain[[i]])
  evbc[[i]] <- dismo::evaluate(prestest[[i]], abstest[[i]], bc[[i]], predictors)
  bcTSS[[i]] <- max(evbc[[i]]@TPR + evbc[[i]]@TNR)-1
  bcAUC[[i]] <- evbc[[i]]@auc
  bckappa[[i]] <- max(evbc[[i]]@kappa)
  bcthres[[i]] <- threshold(evbc[[i]], t.met)
}

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
     'Running GLM (logistic regression) for', sp.n, '...', '\n')
gm <- list()
evgm <- list()
gmTSS <- list()
gmAUC <- list()
gmkappa <- list()
gmthres <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
       'Running GLM (logistic regression) (', i, ') for', sp.n, '...', '\n')
  gm[[i]] <- glm(model, family=binomial(link="logit"), data=predtrain[[i]])
  evgm[[i]] <- dismo::evaluate(prestest[[i]], abstest[[i]], gm[[i]], predictors)
  gmTSS[[i]] <- max(evgm[[i]]@TPR + evgm[[i]]@TNR)-1
  gmAUC[[i]] <- evgm[[i]]@auc
  gmkappa[[i]] <- max(evgm[[i]]@kappa)
  gmthres[[i]] <- threshold(evgm[[i]], t.met)
}

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
     'Running Random Forest model for', sp.n, '...', '\n')
rf <- list()
evrf <- list()
rfTSS <- list()
rfAUC <- list()
rfkappa <- list()
rfthres <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
       'Running Random Forest (', i, ') model for', sp.n, '...', '\n')
  rf[[i]] <- randomForest(model, data=predtrain2[[i]], na.action=na.omit)
  evrf[[i]] <- dismo::evaluate(prestest2[[i]], abstest2[[i]], rf[[i]], predictors)
  rfTSS[[i]] <- max(evrf[[i]]@TPR + evrf[[i]]@TNR)-1
  rfAUC[[i]] <- evrf[[i]]@auc
  rfkappa[[i]] <- max(evrf[[i]]@kappa)
  rfthres[[i]] <- threshold(evrf[[i]], t.met)
}

cat(format( Sys.time(), "%a %b %d %X %Y"), '-', 
    'Running Maxent model for', sp.n, '...', '\n')
mx <- list()
evmx <- list()
mxTSS <- list()
mxAUC <- list()
mxkappa <- list()
mxthres <- list()


for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
      'Running Maxent (', i, ') model for', sp.n, '...', '\n')
  mx[[i]] <- dismo::maxent(predictors, prestrain[[i]], a=bg) #run maxent with java.
  # mx[[i]] <- maxnet(p = p[[i]], data = data[[i]]) 
  # Run maxnet without java. You can use whichever you prefer, maxent or maxnet. 
  # You must only erase # from one and add on the other.
  evmx[[i]] <- dismo::evaluate(prestest[[i]], abstest[[i]], mx[[i]], predictors)
  mxTSS[[i]] <- max(evmx[[i]]@TPR + evmx[[i]]@TNR)-1
  mxAUC[[i]] <- evmx[[i]]@auc
  mxkappa[[i]] <- max(evmx[[i]]@kappa)
  mxthres[[i]] <- threshold(evmx[[i]], t.met)
}

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Running SVM model for', sp.n, '...', '\n')
sv <- list()
evsv <- list()
svTSS <- list()
svAUC <- list()
svkappa <- list()
svthres <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
      'Running SVM (', i, ') model for', sp.n, '...', '\n')
  sv[[i]] <- ksvm(model, data=predtrain[[i]])
  evsv[[i]] <- dismo::evaluate(prestest[[i]], abstest[[i]], sv[[i]], predictors)
  svTSS[[i]] <- max(evsv[[i]]@TPR + evsv[[i]]@TNR)-1
  svAUC[[i]] <- evsv[[i]]@auc
  svkappa[[i]] <- max(evsv[[i]]@kappa)
  svthres[[i]] <- threshold(evsv[[i]], t.met)
}


# Validation table --------------------------------------------------------


cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Generating Validation table for models of', sp.n, '...', '\n')

bcTSSval <- unlist(bcTSS)
#bcTSSval <- data.frame(matrix(unlist(bcTSS), nrow=length(bcTSS), byrow=TRUE), stringsAsFactors=FALSE)
gmTSSval <- unlist(gmTSS)
#gmTSSval <- data.frame(matrix(unlist(gmTSS), nrow=length(gmTSS), byrow=TRUE), stringsAsFactors=FALSE)
rfTSSval <- unlist(rfTSS)
#rfTSSval <- data.frame(matrix(unlist(rfTSS), nrow=length(rfTSS), byrow=TRUE,), stringsAsFactors=FALSE)
mxTSSval <- unlist(mxTSS)
#mxTSSval <- data.frame(matrix(unlist(mxTSS), nrow=length(mxTSS), byrow=TRUE), stringsAsFactors=FALSE)

svTSSval <- unlist(svTSS)
bcAUCval <- unlist(bcAUC)
gmAUCval <- unlist(gmAUC)
rfAUCval <- unlist(rfAUC)
mxAUCval <- unlist(mxAUC)
svAUCval <- unlist(svAUC)
bckappaval <- unlist(bckappa)
gmkappaval <- unlist(gmkappa)
rfkappaval <- unlist(rfkappa)
mxkappaval <- unlist(mxkappa)
svkappaval <- unlist(svkappa)
mod.names <- c(rep ('bc', k), rep('gm', k), rep('rf', k), rep('mx', k), rep('sv', k))
mod.sp <- c(rep(sp.n, k*5)) #conferir a planilha

TSS <- c(bcTSSval, gmTSSval, rfTSSval, mxTSSval, svTSSval)
AUC <- c(bcAUCval, gmAUCval, rfAUCval, mxAUCval, svAUCval)
kappa <- c(bckappaval, gmkappaval, rfkappaval, mxkappaval, svkappaval)
Valid <- data.frame(mod.sp, mod.names, TSS, AUC, kappa, stringsAsFactors=FALSE)

##Pasta de resultados na pasta de cada espécie
write.csv(Valid, file = paste(target_dir, 'Valid_', sp.n, '.csv', sep=""))


# Test if at least one model is valid. 
# If no models are valid it will stop here and continue on the next species

if(all(TSS < tss.lim)){
  print('no models are valid') ##Stop this

  cat(format(Sys.time(), "%a %b %d %X %Y"),
    '-', 'No models were valid for',
    sp.n, 'with ', lim, 'lines...', '\n')
  
  #save.image("./outputs/my_analysis.rData")

  finished_time = Sys.time()
  cat( format( finished_time, "%a %b %d %X %Y"), '-', 'FINISHED', '\n')
  write(format(finished_time, "%a %b %d %X %Y"),
        file = paste0(target_dir, "FINISHED.txt", sep = ""))
  
  next
}


# Projecting Bioclim -------------------------------------------------------

#Drop bad models, project under current and future conditions
cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
     'Projecting Bioclim model of', sp.n, '...', '\n')
cur.bc <- list()
cur.bc.bin <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
      'Projecting Bioclim (', i, ') model of', sp.n, '...', '\n')
  if(bcTSS[[i]] >= tss.lim){
    cur.bc[[i]] <- predict(predictors, bc[[i]])
    cur.bc.bin[[i]] <- cur.bc[[i]] > bcthres[[i]]
    
  } else {
    cur.bc[[i]] <- NULL
    cur.bc.bin[[i]] <- NULL
      }
}

##########
cur.bc <- Filter(Negate(is.null), cur.bc) # remove null rasters
cur.bc.bin <- Filter(Negate(is.null), cur.bc.bin)

# Ensemble of binary models (Majority Rule)
cur.bc.ens <- Reduce('+', cur.bc.bin) #sum of pixels
tval <- unique(cur.bc.ens)
tval <- tval[tval != 0]
tval <- median(tval)

if(!is.null(cur.bc.ens)){
  cur.bc.ens.bin <- cur.bc.ens >= tval} else {
    cur.bc.ens.bin <- NULL
  }


##########

# Taking out the zero values and normalizing rasters

if(length(cur.bc) != 0) {
  for(z in 1:length(cur.bc)) {
    adeq = cur.bc[[z]]
    minimo <- min(adeq[], na.rm=T)
    maximo <- max(adeq[], na.rm=T)
    adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
    cur.bc[[z]] <- calc(adeq, adeq_norm)
    }

  #Ensemble of continuos models
  x <- stack(cur.bc)
  cur.bc.cont <- calc(x, fun = mean)
  
} else {
  cur.bc.cont <- NULL
}


# Projecting GLM ----------------------------------------------------------


cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Projecting GLMs of', sp.n, '...', '\n')
cur.gm <- list()
cur.gm.bin <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
       'Projecting GLMs (', i, ') of', sp.n, '...', '\n')
  if(gmTSS[[i]] >= tss.lim){
    cur.gm[[i]] <- predict(predictors, gm[[i]])
    cur.gm.bin[[i]] <- cur.gm[[i]] > gmthres[[i]]
  } else {
    cur.gm[[i]] <- NULL
    cur.gm.bin[[i]] <- NULL
 }
}

##########

cur.gm <- Filter(Negate(is.null), cur.gm)
cur.gm.bin <- Filter(Negate(is.null), cur.gm.bin)

cur.gm.ens <- Reduce('+', cur.gm.bin)
tval <- unique(cur.gm.ens)
tval <- tval[tval != 0]
tval <- median(tval)

if(!is.null(cur.gm.ens)){
cur.gm.ens.bin <- cur.gm.ens >= tval} else {
  cur.gm.ens.bin <- NULL
}

##########

if(length(cur.gm) != 0) {
for(z in 1:length(cur.gm)){
  adeq = cur.gm[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  cur.gm[[z]] <- calc(adeq, adeq_norm)
}

x <- stack(cur.gm)
cur.gm.cont <- calc(x, fun = mean)
} else {
  cur.gm.cont <- NULL
}


# Projecting RF ----------------------------------------------------------

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Projecting Random Forest models of', sp.n, '...', '\n')
cur.rf <- list()
cur.rf.bin <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
      'Projecting Random Forest (', i, ') models of', sp.n, '...', '\n')
  if(rfTSS[[i]] >= tss.lim){
    cur.rf[[i]] <- predict(predictors, rf[[i]])
    cur.rf.bin[[i]] <- cur.rf[[i]] > rfthres[[i]]
 } else {

    cur.rf[[i]] <- NULL
    cur.rf.bin[[i]] <- NULL
  }
}

##########
cur.rf <- Filter(Negate(is.null), cur.rf)
cur.rf.bin <- Filter(Negate(is.null), cur.rf.bin)

cur.rf.ens <- Reduce('+', cur.rf.bin)
tval <- unique(cur.rf.ens)
tval <- tval[tval != 0]
tval <- median(tval)

if(!is.null(cur.rf.ens)){
  cur.rf.ens.bin <- cur.rf.ens >= tval} else {
    cur.rf.ens.bin <- NULL
  }

##########

if(length(cur.rf) != 0) {
for(z in 1:length(cur.rf)){
  adeq = cur.rf[[z]]
  if(sum(adeq[], na.rm=T)!=0){
    minimo <- min(adeq[], na.rm=T)
    maximo <- max(adeq[], na.rm=T)
    adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
    cur.rf[[z]] <- calc(adeq, adeq_norm)
  }

}

cur.rf <- Filter(Negate(is.null), cur.rf)

x <- stack(cur.rf)
cur.rf.cont <- calc(x, fun = mean)

} else {
  cur.rf.cont <- NULL
}


# Projecting Maxent --------------------------------------------------------

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Projecting Maxent models of', sp.n, '...', '\n')
cur.mx <- list()
cur.mx.bin <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
      'Projecting Maxent (', i, ') models of', sp.n, '...', '\n')
  if(mxTSS[[i]] >= tss.lim){
    cur.mx[[i]] <- predict(predictors, mx[[i]])
    cur.mx.bin[[i]] <- cur.mx[[i]] > mxthres[[i]]
  } else {
    cur.mx[[i]] <- NULL
    cur.mx.bin[[i]] <- NULL
  }
}


##########
cur.mx <- Filter(Negate(is.null), cur.mx)
cur.mx.bin <- Filter(Negate(is.null), cur.mx.bin)

cur.mx.ens <- Reduce('+', cur.mx.bin)
tval <- unique(cur.mx.ens)
tval <- tval[tval != 0]
tval <- median(tval)

if(!is.null(cur.mx.ens)){
  cur.mx.ens.bin <- cur.mx.ens >= tval} else {
    cur.mx.ens.bin <- NULL
  }

##########

if(length(cur.mx) != 0) {
for(z in 1:length(cur.mx)){
  adeq = cur.mx[[z]]
  if(sum(adeq[], na.rm=T)!=0){
    minimo <- min(adeq[], na.rm=T)
    maximo <- max(adeq[], na.rm=T)
    adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
    cur.mx[[z]] <- calc(adeq, adeq_norm)
  }
}

x <- stack(cur.mx)
cur.mx.cont <- calc(x, fun = mean)

} else {
  cur.mx.cont <- NULL
}


# Projecting SVM ----------------------------------------------------------

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Projecting SVM models of', sp.n, '...', '\n')
cur.sv <- list()
cur.sv.bin <- list()

for(i in 1:k){
  cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
      'Projecting SVM (', i, ') models of', sp.n, '...', '\n')
  if(svTSS[[i]] >= tss.lim){
    cur.sv[[i]] <- predict(predictors, sv[[i]])
    cur.sv.bin[[i]] <- cur.sv[[i]] > svthres[[i]]

  } else {
    cur.sv[[i]] <- NULL
    cur.sv.bin[[i]] <- NULL
   }
}


##########
cur.sv <- Filter(Negate(is.null), cur.sv)
cur.sv.bin <- Filter(Negate(is.null), cur.sv.bin)

cur.sv.ens <- Reduce('+', cur.sv.bin)
tval <- unique(cur.sv.ens)
tval <- tval[tval != 0]
tval <- median(tval)

if(!is.null(cur.sv.ens)){
  cur.sv.ens.bin <- cur.sv.ens >= tval} else {
    cur.sv.ens.bin <- NULL
  }

##########

if(length(cur.sv) != 0) {

for(z in 1:length(cur.sv)){
  adeq = cur.sv[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  cur.sv[[z]] <- calc(adeq, adeq_norm)
}

x <- stack(cur.sv)
cur.sv.cont <- calc(x, fun = mean)
} else {
  cur.sv.cont <- NULL
}


# Ensemble ----------------------------------------------------------------
###################### ensemble current for ALL ALGORITHMS

# For binary models, by Majority Rule

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Projecting ensemble models of', sp.n, '...', '\n')

algo.cur <- list(cur.bc.ens.bin, cur.gm.ens.bin, cur.rf.ens.bin, cur.mx.ens.bin, cur.sv.ens.bin)
algo.cur <- algo.cur[lapply(algo.cur, length) > 0] #ignores NULL values already
algo.cur2 <- algo.cur[!is.na(algo.cur)]
ens.cur <- Reduce('+', algo.cur2)
tval <- unique(ens.cur)
tval <- tval[tval != 0]
tval <- median(tval)
ens.cur.bin <- ens.cur >= tval  # this one! ensemble for current binary model


# For continuous model, by TSS-weighted average

w <- c(bcTSSval, gmTSSval, rfTSSval, mxTSSval, svTSSval) #ignores NULL values already

if (length(cur.bc) != 0) {
st1 <- stack(cur.bc)
for(i in 1:length(cur.bc)){
  min_max <- range(cur.bc[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.bc", i, sp.n, ".txt",  sep="_"))
}} else {
  st1 <- NULL
  }

if (length(cur.gm) != 0){
st2 <- stack(cur.gm)
for(i in 1:length(cur.gm)){
  min_max <- range(cur.gm[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.gm", i, sp.n, ".txt",  sep="_"))
}

} else {
  st2 <- NULL
}

if (length(cur.rf) != 0){
st3 <- stack(cur.rf)
for(i in 1:length(cur.rf)){
  min_max <- range(cur.rf[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.rf", i, sp.n, ".txt",  sep="_"))
}
} else {
  st3 <- NULL
}

if (length(cur.mx) != 0){
st4 <- stack(cur.mx)
for(i in 1:length(cur.mx)){
  min_max <- range(cur.mx[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.mx", i, sp.n, ".txt",  sep="_"))
}
#rm(cur.mx)
} else {
  st4 <- NULL
}

if (length(cur.sv) != 0){
st5 <- stack(cur.sv)
for(i in 1:length(cur.sv)){
  min_max <- range(cur.sv[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.sv", i, sp.n, ".txt",  sep="_"))
}

} else {
  st5 <- NULL
}

stacks <- c(st1, st2, st3, st4, st5) 
st <- stack(stacks)

w.sem.os.nulos <- w[w >= tss.lim]
ens.cur <- weighted.mean(st, w.sem.os.nulos) # this one! ensemble for current continuous model
ens.cur.sd.w <- sum(w.sem.os.nulos * (st - ens.cur)^2) # variance
ens.cur.sd.w <- sqrt(ens.cur.sd.w) # standard deviation (i.e. uncertainty)


# Uncertainty for each algo (each partition) -------------------------------------------------------------

wbc <- c(bcTSSval)
wbc.sem.os.nulos <- wbc[wbc >= tss.lim]

if(length(wbc.sem.os.nulos) == 0) {
 print('bioclim is null')
  cur.bc.mean.w <- NULL
  cur.bc.sd.w <- NULL
  }else {
cur.bc.mean.w <- weighted.mean(st1, wbc.sem.os.nulos)
cur.bc.sd.w <- sum(wbc.sem.os.nulos * (st1 - cur.bc.mean.w)^2)
cur.bc.sd.w <- sqrt(cur.bc.sd.w) }

wgm <- c(gmTSSval)
wgm.sem.os.nulos <- wgm[wgm >= tss.lim]
if(length(wgm.sem.os.nulos)==0){
  print('glm is null')
  cur.gm.mean.w <- NULL
  cur.gm.sd.w <- NULL
  } else {
cur.gm.mean.w <- weighted.mean(st2, wgm.sem.os.nulos)
cur.gm.sd.w <- sum(wgm.sem.os.nulos * (st2 - cur.gm.mean.w)^2)
cur.gm.sd.w <- sqrt(cur.gm.sd.w)
  }

wrf <- c(rfTSSval)
wrf.sem.os.nulos <- wrf[wrf >= tss.lim]
if (length(wrf.sem.os.nulos)==0){
  print('rf is null')
  cur.rf.mean.w <- NULL
  cur.rf.sd.w <- NULL
  } else {
cur.rf.mean.w <- weighted.mean(st3, wrf.sem.os.nulos)
cur.rf.sd.w <- sum(wrf.sem.os.nulos * (st3 - cur.rf.mean.w)^2)
cur.rf.sd.w <- sqrt(cur.rf.sd.w)
  }

wmx <- c(mxTSSval)
wmx.sem.os.nulos <- wmx[wmx >= tss.lim]
if (length(wmx.sem.os.nulos)==0) {
  print('maxent is null')
  cur.mx.mean.w <- NULL
  cur.mx.sd.w <- NULL
} else {
cur.mx.mean.w <- weighted.mean(st4, wmx.sem.os.nulos)
cur.mx.sd.w <- sum(wmx.sem.os.nulos * (st4 - cur.mx.mean.w)^2)
cur.mx.sd.w <- sqrt(cur.mx.sd.w)
}

wsv <- c(svTSSval)
wsv.sem.os.nulos <- wsv[wsv >= tss.lim]
if (length(wsv.sem.os.nulos)==0){
  print('svmk is null')
  cur.sv.mean.w <- NULL
  cur.sv.sd.w <-  NULL
} else {
cur.sv.mean.w <- weighted.mean(st5, wsv.sem.os.nulos)
cur.sv.sd.w <- sum(wsv.sem.os.nulos * (st5 - cur.sv.mean.w)^2)
cur.sv.sd.w <- sqrt(cur.sv.sd.w)
}


###SAVING

cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Saving ensemble maps of', sp.n, '...', '\n')

if(!is.null(ens.cur)){
  
  writeRaster(ens.cur,
    file = paste(target_dir, '/CUR.cont_', sp.n, '.asc', sep = ""),
    overwrite = TRUE)
}

if(!is.null(ens.cur.bin)){

  writeRaster(ens.cur.bin,
    file = paste(target_dir, '/CUR.bin_', sp.n, '.asc', sep = ""),
    overwrite = TRUE)
}

#weighted.mean and weighted.sd
if(!is.null(ens.cur.sd.w)){
  
  writeRaster(ens.cur.sd.w,
    file = paste(target_dir, 'ens.cur.sd.w_', sp.n, '.asc', sep = ""),
    overwrite = TRUE)
}

if(!is.null(cur.gm.sd.w)){
writeRaster(cur.gm.sd.w, file = paste(target_dir, 'cur.gm.sd.w_', 
                                      sp.n, '.asc', sep=""),overwrite=TRUE)
}

if(!is.null(cur.rf.sd.w)){
writeRaster(cur.rf.sd.w, file = paste(target_dir, 'cur.rf.sd.w_', 
                                      sp.n, '.asc', sep=""),overwrite=TRUE)
}

if(!is.null(cur.mx.sd.w)){
writeRaster(cur.mx.sd.w, file = paste(target_dir, 'cur.mx.sd.w_', 
                                      sp.n, '.asc', sep=""),overwrite=TRUE)
}

if(!is.null(cur.sv.sd.w)){
writeRaster(cur.sv.sd.w, file = paste(target_dir, 'cur.sv.sd.w_', 
                                      sp.n, '.asc', sep=""),overwrite=TRUE)
}


cat(format(Sys.time(), "%a %b %d %X %Y"), '-', 
    'Finished train and test datasets for', sp.n, 'with ', lim, 'lines...', '\n')

#save.image("./outputs/my_analysis.rData")

finished_time = Sys.time()
cat( format( finished_time, "%a %b %d %X %Y"), '-', 'FINISHED', '\n')
write(format( finished_time, "%a %b %d %X %Y"), 
      file = paste0(target_dir, "FINISHED.txt", sep=""))

}

