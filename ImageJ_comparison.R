library(ggplot2)
library(mnsentinellakes)
library(dplyr)
library(stringr)
library(SimilarityMeasures)
library(mosaic)

DEM_hypsos <- readRDS("data/MN_bathy.rds")
DEM_hypsos$DOW <- fixlakeid(DEM_hypsos$DOW)

# Hypsos Chris did
Perc_CR <- "data/Fourth_set" #file name with the files to be checked for QC
Perc_files_CR <- list.files(Perc_CR, pattern = "csv", full.names = T)
# Make sure these are QC'ed first!
# Get DOW  
Perc_DOWs_CR <- basename(Perc_files_CR) %>% str_extract(., '[0-9]+.csv$') %>% str_extract('[0-9]+')
Perc_DOWs_CR <- fixlakeid(Perc_DOWs_CR) # 2 missing leading zeroes
DEM_hypsos_ours <- DEM_hypsos %>% filter(DOW %in% Perc_DOWs_CR)
length(unique(DEM_hypsos_ours$DOW)) #should be 100


# Hypsos Amanda did
Perc_AVP <- "data/Fourth_set_AVP" #file name with the files to be checked for QC
Perc_files_AVP <- list.files(Perc_AVP, pattern = "csv", full.names = T)
# Make sure these are QC'ed first!
# Get DOW  
Perc_DOWs_AVP <- basename(Perc_files_AVP) %>% str_extract(., '[0-9]+.csv$') %>% str_extract('[0-9]+')
Perc_DOWs_AVP <- fixlakeid(Perc_DOWs_AVP) # 2 missing leading zeroes


#extracts the depth data from the DEMs
DEM_Hypso_Perc <- NA
for (i in seq_along(Perc_DOWs_CR)) {
  dat <- DEM_hypsos_ours %>% filter(DOW == Perc_DOWs_CR[i])
  dat$depth_feet <- dat$depths / 0.3048
  dat$proportion_area <- dat$areas / max(dat$areas)
  dat <- dat[,-2] #%>% dplyr::select(DOW, depths, depth_feet, proportion_area)
  DEM_Hypso_Perc <- rbind(DEM_Hypso_Perc, dat)
}

DEM_Hypso_Perc <- DEM_Hypso_Perc %>% filter(!is.na(DOW)) # removes NA first row
DEM_Hypso_Perc$Method <- "DEM"



#extracts the depth data from chris' hypsos
Perc_Hypsos_CR <- NA
for (i in seq_along(Perc_DOWs_CR)) {
  dat <- read.csv(Perc_files_CR[i])[,1:2] %>% setNames(c('depth_feet', 'proportion_area'))
  dat$DOW <- Perc_DOWs_CR[i]
  dat$depths <- dat$depth_feet * 0.3048 # this is depth in m if we want it later
  Perc_Hypsos_CR <- rbind(Perc_Hypsos_CR, dat)
}

Perc_Hypsos_CR <- Perc_Hypsos_CR %>% filter(!is.na(DOW)) # removes NA first row
Perc_Hypsos_CR$Method <- "ImageJ_CR"



#extracts the depth data from Amanda's hypsos
Perc_Hypsos_AVP <- NA
for (i in seq_along(Perc_DOWs_AVP)) {
  dat_AVP <- read.csv(Perc_files_AVP[i])[,1:2] %>% setNames(c('depth_feet', 'proportion_area'))
  dat_AVP$DOW <- Perc_DOWs_AVP[i]
  dat_AVP$depths <- dat_AVP$depth_feet * 0.3048
  Perc_Hypsos_AVP <- rbind(Perc_Hypsos_AVP, dat_AVP)
}

Perc_Hypsos_AVP <- Perc_Hypsos_AVP %>% filter(!is.na(DOW)) # removes NA first row
Perc_Hypsos_AVP$Method <- "ImageJ_AVP"


All_Hypsos <- rbind(DEM_Hypso_Perc, Perc_Hypsos_CR, Perc_Hypsos_AVP)
#write.csv(All_Hypsos, "data/all_validation_hypsos.csv", row.names = F)


Hypsos_depth <- All_Hypsos
Hypsos_depth$proportion_depth = NA
maxdepth_lake <- data.frame("dow" = character(1), "max_depth" = numeric(1), 
                            "method" = character(1), stringsAsFactors = FALSE)

#get the maximum depth for each lake using each method
max_depth = Hypsos_depth[1,3]; dow = Hypsos_depth[1,2] 
method = "DEM"; i = 1
for (row in 1:nrow(Hypsos_depth)) {
  if (dow != Hypsos_depth[row,2]) {
    maxdepth_lake[i,1] = dow
    maxdepth_lake[i,2] = max_depth
    maxdepth_lake[i,3] = method
    i = i + 1
    dow = Hypsos_depth[row,2]
    max_depth = Hypsos_depth[row,3]
  }
  else {
    method = Hypsos_depth[row,5]
    max_depth = Hypsos_depth[row,3]
  }
}
maxdepth_lake[i,1] = dow; maxdepth_lake[i,2] = max_depth; maxdepth_lake[i,3] = method
#write.csv(maxdepth_lake, "data/max_depth_validation.csv", row.names = F)


### CALCULATE AVERAGE VERTICAL DISTANCE BETWEEN THE CURVES ####


## Difference between hypsos Chris did and DEM
avg.abs.diff <- rep(NA, length(Perc_DOWs_CR))
root.meansq.diff <- rep(NA, length(Perc_DOWs_CR))
avg.ImageJ.DEM.diff <- rep(NA, length(Perc_DOWs_CR))

for (i in seq_along(Perc_DOWs_CR)) {
  
  ImageJ.hypso <- All_Hypsos %>% filter(DOW == Perc_DOWs_CR[i], Method == "ImageJ_CR")
  DEM.hypso <- All_Hypsos %>% filter(DOW == Perc_DOWs_CR[i], Method == "DEM")
  
  max.depth <- round(max(ImageJ.hypso$depth_feet, DEM.hypso$depth_feet), digits = 1)
  depth.vals <- seq(0, max.depth, .1)
  
  ImageJ.interp <- approx(ImageJ.hypso$depth_feet, ImageJ.hypso$proportion_area, 
                          xout = depth.vals, yright = 0, yleft = 1)
  DEM.interp <- approx(DEM.hypso$depth_feet, DEM.hypso$proportion_area, 
                       xout = depth.vals, yright = 0, yleft = 1)
  
  # Avg absolute difference
  avg.abs.diff[i] <- mean(abs(ImageJ.interp$y - DEM.interp$y))
  # Root mean squared difference
  root.meansq.diff[i] <- sqrt(mean((ImageJ.interp$y - DEM.interp$y)^2))
  # Avg signed difference between ImageJ and DEM
  avg.ImageJ.DEM.diff[i] <- mean(ImageJ.interp$y - DEM.interp$y)
}

Distances_CR <- data.frame(DOW = Perc_DOWs_CR, Mean_Absolute_Diff = avg.abs.diff, 
                           Mean_ImageJ_DEM_Diff = avg.ImageJ.DEM.diff, RMSDiff = root.meansq.diff)
Distances_CR <- arrange(Distances_CR, DOW)


## Difference between hypsos Amanda did and DEM
avg.abs.diff <- rep(NA, length(Perc_DOWs_AVP))
root.meansq.diff <- rep(NA, length(Perc_DOWs_AVP))
avg.ImageJ.DEM.diff <- rep(NA, length(Perc_DOWs_AVP))

for (i in seq_along(Perc_DOWs_AVP)) {
  
  ImageJ.hypso <- All_Hypsos %>% filter(DOW == Perc_DOWs_AVP[i], Method == "ImageJ_AVP")
  DEM.hypso <- All_Hypsos %>% filter(DOW == Perc_DOWs_AVP[i], Method == "DEM")
  
  max.depth <- round(max(ImageJ.hypso$depth_feet, DEM.hypso$depth_feet), digits = 1)
  depth.vals <- seq(0, max.depth, .1)
  
  ImageJ.interp <- approx(ImageJ.hypso$depth_feet, ImageJ.hypso$proportion_area, 
                          xout = depth.vals, yright = 0, yleft = 1)
  DEM.interp <- approx(DEM.hypso$depth_feet, DEM.hypso$proportion_area, 
                       xout = depth.vals, yright = 0, yleft = 1)
  
  # Avg absolute difference
  avg.abs.diff[i] <- mean(abs(ImageJ.interp$y - DEM.interp$y))
  # Root mean squared difference
  root.meansq.diff[i] <- sqrt(mean((ImageJ.interp$y - DEM.interp$y)^2))
  # Avg signed difference between ImageJ and DEM
  avg.ImageJ.DEM.diff[i] <- mean(ImageJ.interp$y - DEM.interp$y)
}

Distances_AVP <- data.frame(DOW = Perc_DOWs_AVP, Mean_Absolute_Diff = avg.abs.diff, 
                            Mean_ImageJ_DEM_Diff = avg.ImageJ.DEM.diff, RMSDiff = root.meansq.diff)
Distances_AVP <- arrange(Distances_AVP, DOW)



## Difference between hypsos Amanda and Chris did
avg.abs.diff <- rep(NA, length(Perc_DOWs_AVP))
root.meansq.diff <- rep(NA, length(Perc_DOWs_AVP))
avg.ImageJ.DEM.diff <- rep(NA, length(Perc_DOWs_AVP))

for (i in seq_along(Perc_DOWs_AVP)) {
  
  AVP.hypso <- All_Hypsos %>% filter(DOW == Perc_DOWs_AVP[i], Method == "ImageJ_AVP")
  CR.hypso <- All_Hypsos %>% filter(DOW == Perc_DOWs_AVP[i], Method == "ImageJ_CR")
  
  max.depth <- round(max(AVP.hypso$depth_feet, CR.hypso$depth_feet), digits = 1)
  depth.vals <- seq(0, max.depth, .1)
  
  AVP.interp <- approx(AVP.hypso$depth_feet, AVP.hypso$proportion_area, 
                          xout = depth.vals, yright = 0, yleft = 1)
  CR.interp <- approx(CR.hypso$depth_feet, CR.hypso$proportion_area, 
                       xout = depth.vals, yright = 0, yleft = 1)
  
  # Avg absolute difference
  avg.abs.diff[i] <- mean(abs(AVP.interp$y - CR.interp$y))
  # Root mean squared difference
  root.meansq.diff[i] <- sqrt(mean((AVP.interp$y - CR.interp$y)^2))
  # Avg signed difference between ImageJ and DEM
  avg.ImageJ.DEM.diff[i] <- mean(AVP.interp$y - CR.interp$y)
}

Distances_AVP_CR <- data.frame(DOW = Perc_DOWs_AVP, Mean_Absolute_Diff = avg.abs.diff, 
                            Mean_ImageJ_DEM_Diff = avg.ImageJ.DEM.diff, RMSDiff = root.meansq.diff)
Distances_AVP_CR <- arrange(Distances_AVP_CR, DOW)

# Once AVP finishes write csv the merge of all the distances df


quantile(Distances_CR$Mean_Absolute_Diff, seq(0,1, .05))

mean(Distances_CR$Mean_Absolute_Diff) # .051
mean(Distances_AVP$Mean_Absolute_Diff) # 0.047
mean(Distances_AVP_CR$Mean_Absolute_Diff) #0.016

median(Distances_CR$Mean_Absolute_Diff) # .032
median(Distances_AVP$Mean_Absolute_Diff)
median(Distances_AVP_CR$Mean_Absolute_Diff)


se(Distances_CR$Mean_Absolute_Diff) # .0055 # need se function

t.test(Distances_CR$Mean_Absolute_Diff,   conf.level = 0.95)
# Mean absolute distance is .051 (.040, .061)
t.test(Distances_AVP$Mean_Absolute_Diff,   conf.level = 0.95)

# Alternative: Bootstrap mean and CI
samp.mean <- do(10000)*{mean(~Mean_Absolute_Diff, data = resample(Distances_CR))}
head(samp.mean)
(CIb <- qdata(~result, data = samp.mean, p = c(0.025, 0.975)))
mean(~result, data = samp.mean)
# .051 (.041, .062)

sum(Distances_CR$Mean_Absolute_Diff > 0.10)/100 # 14% (14) greater than 10% difference on average
sum(Distances_CR$Mean_Absolute_Diff > 0.15)/100 # 9% (9) greater than 15% difference on average

Distances.df <- arrange(Distances_CR, desc(Mean_Absolute_Diff))
head(Distances_CR, 10)
par(mfrow = c(2,1))

# Mean absolute differences
hist(Distances_CR$Mean_Absolute_Diff, breaks = 20, xlim = c(0,.3), 
     ylim = c(0,25), col = "blue", border = "white",  las = 1, xlab = "Mean absolute difference (proportion of lake area)",main = "")
text(0.25, 22, "A", col = "red", cex = 2)

# Average signed difference between ImageJ and DEM
hist(Distances_CR$Mean_ImageJ_DEM_Diff, breaks = 20, xlim = c(-0.2,0.3), 
     col = "blue", border = "white", las = 1, xlab = "Mean signed difference (proportion of lake area)", main = "")
text(0.22, 22, "B", col = "red", cex = 2)
mean(Distances_CR$Mean_ImageJ_DEM_Diff)
t.test(Distances_CR$Mean_ImageJ_DEM_Diff, conf.level = 0.95)





### Calculate depth information for all MN and WI lakes (not just validation) #####

ALL_MN <- "data/First_set" 
allMN_files <- list.files(ALL_MN, pattern = "csv", full.names = T)

ALL_WI <- "data/WI_set"
allWI_files <- list.files(ALL_WI, patter = "csv", full.names = T)

# Get DOW  
PercMN_DOWs <- basename(allMN_files) %>% str_extract(., '[0-9]+.csv$') %>% str_extract('[0-9]+')
PercMN_DOWs <- fixlakeid(PercMN_DOWs) # 2 missing leading zeroes
unique(PercMN_DOWs)

# Get WBIC
PercWI_DOWs <- basename(allWI_files) %>% str_extract(., '[0-9]+.csv$') %>% str_extract('[0-9]+')

MN_Hypsos <- NA
for (i in seq_along(PercMN_DOWs)) {
  dat <- read.csv(allMN_files[i])[,1:2] %>% setNames(c('depth_feet', 'proportion_area'))
  dat$DOW <- PercMN_DOWs[i]
  MN_Hypsos <- rbind(MN_Hypsos, dat)
}
MN_Hypsos <- MN_Hypsos %>% filter(!is.na(DOW)) # removes NA first row

WI_Hypsos <- NA
for (i in seq_along(PercWI_DOWs)) {
  dat <- read.csv(allWI_files[i])[,1:2] %>% setNames(c('depth_feet', 'proportion_area'))
  dat$DOW <- PercWI_DOWs[i]
  WI_Hypsos <- rbind(WI_Hypsos, dat)
}
WI_Hypsos <- WI_Hypsos %>% filter(!is.na(DOW)) # removes NA first row





maxdepth <- data.frame("dow" = character(1), "max_depth" = numeric(1), stringsAsFactors = FALSE)
Hypsos <- rbind(WI_Hypsos, MN_Hypsos)
unique(Hypsos$DOW)
max_depth = Hypsos[1,1]; dow = Hypsos[1,3]

i = 1
for (row in 1:nrow(Hypsos_depth)) {
  if (dow != Hypsos[row,3]) {
    maxdepth[i, 1] = dow
    maxdepth[i, 2] = max_depth
    i = i + 1
    dow = Hypsos[row,3]
    max_depth = Hypsos[row, 1]
  }
  else {
    max_depth = Hypsos[row, 1]
  }
}
maxdepth[i,1] = dow; maxdepth[i,2] = max_depth

# statistics for all lake depths
se <- function(x) sqrt(var(x)/length(x))
mean(maxdepth$max_depth)  #31.8
confint(maxdepth$max_depth)
se(maxdepth$max_depth)    #1.13
min(maxdepth$max_depth)   #4
max(maxdepth$max_depth)   #443
median(maxdepth$max_depth) #25
hist(maxdepth$max_depth, breaks = 50, xlim = c(0 , 450))
