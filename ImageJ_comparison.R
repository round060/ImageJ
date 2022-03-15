#remotes::install_github("GLEON/GLM3r")
#remotes::install_github('usgs-r/glmtools')
#remotes::install_github("jsta/glmutils")
library(glmtools)
library(ggplot2)
library(mnsentinellakes)
library(dplyr)
library(stringr)
library(SimilarityMeasures)
library(mosaic)

se <- function(x) sqrt(var(x)/length(x))
lake_areas <- read.csv("data/MN_areas.csv")

DEM_hypsos <- readRDS("data/MN_bathy.rds")
DEM_hypsos$DOW <- fixlakeid(DEM_hypsos$DOW)

# Hypsos Chris did
Perc_CR <- "data/Fourth_set_CR" #file name with the files to be checked for QC
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
maxdepth_lake <- data.frame("DOW" = character(1), "max_depth" = numeric(1), 
                            "method" = character(1), stringsAsFactors = FALSE)

#get the maximum depth for each lake using each method
max_depth = Hypsos_depth[1,3]; DOW = Hypsos_depth[1,2] 
method = "DEM"; i = 1
for (row in 1:nrow(Hypsos_depth)) {
  if (DOW != Hypsos_depth[row,2]) {
    maxdepth_lake[i,1] = DOW
    maxdepth_lake[i,2] = max_depth
    maxdepth_lake[i,3] = method
    i = i + 1
    DOW = Hypsos_depth[row,2]
    max_depth = Hypsos_depth[row,3]
  }
  else {
    method = Hypsos_depth[row,5]
    max_depth = Hypsos_depth[row,3]
  }
}
maxdepth_lake[i,1] = DOW; maxdepth_lake[i,2] = max_depth; maxdepth_lake[i,3] = method
#write.csv(maxdepth_lake, "data/max_depth_validation.csv", row.names = F)


### calulate what hypsography would be with only max-depth ###
#Use package General Lake Model#
all_hypsos <- subset(All_Hypsos, select = -c(depths))
#sets up the null models
null_hypso = data.frame(DOW = integer(), depth_feet = numeric(), proportion_area = numeric(), Method = character())
null_hypso = null_hypso[-1,]
max_depth <- maxdepth_lake %>% filter(method == "DEM")
for (i in 1:nrow(max_depth)) {
  if (max_depth[i,2] < 10) { #specifies the number of depths to interpolate depth from, more = more measurements
    nlayers = max_depth[i,2]*2
  }else{
    nlayers = round(max_depth[i,2])
  }
  res <- suppressWarnings( # does the actual generation of hypsography 
    glmutils::generate_hypsography(round(nlayers), max_depth[i,2], 1))
  res$DOW = max_depth[i,1]
  
  #Probably the worst way to do this because its running at like n^2 speed
  null_hypso = cbind(DOW = res[,4], depth_feet = res[,3], proportion_area = res[,1], Method = "ellipsoid")
  all_hypsos = rbind(all_hypsos, null_hypso)
  null_hypso = cbind(DOW = res[,4], depth_feet = res[,2], proportion_area = res[,1], Method = "linear") 
  all_hypsos = rbind(all_hypsos, null_hypso)
}

#I know im bad
all_hypsos$depth_feet = as.numeric(all_hypsos$depth_feet)
all_hypsos$proportion_area = as.numeric(all_hypsos$proportion_area)


### CALCULATE AVERAGE VERTICAL DISTANCE BETWEEN THE CURVES ####


DEM_distance <- dplyr::filter(all_hypsos, Method == "DEM")
non_DEM_distance <- dplyr::filter(all_hypsos, Method != "DEM")
methods <- unique(non_DEM_distance$Method)
Distances <- data.frame(DOW = NA, Method = NA, Mean_Absolute_Diff = NA, 
                        Mean_DEM_Diff = NA, RMS_Diff = NA)

for (j in 1:length(methods)) {
  
  temp <- dplyr::filter(all_hypsos, Method == methods[j])
  temp_dow <- unique(temp$DOW)
  avg.abs.diff <- rep(NA, length(unique(temp$DOW)))
  root.meansq.diff <- rep(NA, length(unique(temp$DOW)))
  avg.DEM.diff <- rep(NA, length(unique(temp$DOW)))
  
  for (i in seq_along(temp_dow)) {
    
    digitized.hypso <- temp %>% filter(DOW == temp_dow[i])
    DEM.hypso <- DEM_distance %>% filter(DOW == temp_dow[i])
    
    max.depth <- round(max(digitized.hypso$depth_feet, DEM.hypso$depth_feet), digits = 1)
    depth.vals <- seq(0, max.depth, .1)
    
    digitized.interp <- approx(digitized.hypso$depth_feet, digitized.hypso$proportion_area, 
                            xout = depth.vals, yright = 0, yleft = 1)
    DEM.interp <- approx(DEM.hypso$depth_feet, DEM.hypso$proportion_area, 
                         xout = depth.vals, yright = 0, yleft = 1)
    
    # Avg absolute difference
    avg.abs.diff[i] <- mean(abs(digitized.interp$y - DEM.interp$y))
    # Root mean squared difference
    root.meansq.diff[i] <- sqrt(mean((digitized.interp$y - DEM.interp$y)^2))
    # Avg signed difference between ImageJ and DEM
    avg.DEM.diff[i] <- mean(digitized.interp$y - DEM.interp$y)
  } #first for loop
  
  Distances_temp <- data.frame(DOW = temp_dow, Method = methods[j], Mean_Absolute_Diff = avg.abs.diff, 
                             Mean_DEM_Diff = avg.DEM.diff, RMS_Diff = root.meansq.diff)
  Distances_temp <- arrange(Distances_temp, DOW)
  Distances <- rbind(Distances, Distances_temp)
  
} #second for loop
Distances <- na.omit(Distances)


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


mosaic::mean(Mean_Absolute_Diff~Method, data = Distances)
mean(Distances_AVP_CR$Mean_Absolute_Diff) #0.0157

mosaic::median(Mean_Absolute_Diff~Method, data = Distances)
median(Distances_AVP_CR$Mean_Absolute_Diff) # 0.0078

mosaic::sd(Mean_Absolute_Diff~Method, data = Distances)

Distances_CR <- Distances %>%
  dplyr::filter(Method == "ImageJ_CR")
Distances_AVP <- Distances %>%
  dplyr::filter(Method == "ImageJ_AVP")

t.test(Distances_CR$Mean_Absolute_Diff,   conf.level = 0.95)
# Mean absolute distance is .051 (.040, .062)
t.test(Distances_AVP$Mean_Absolute_Diff,   conf.level = 0.95)
# Mean absolute distance is .046 (.035, .057)
t.test(Distances_AVP_CR$Mean_Absolute_Diff, conf.level = 0.95)
# Mean absolute distance is .016 (.011, .021)

sum(Distances_CR$Mean_Absolute_Diff > 0.10)/100 # 14% (14) greater than 10% difference on average
sum(Distances_CR$Mean_Absolute_Diff > 0.15)/100 # 9% (9) greater than 15% difference on average



### ggplots ###

# histogram MAD and MSD
plot_df <- Distances %>%
  dplyr::filter(Method == "ImageJ_AVP" | Method == "ImageJ_CR")
 
#change legend attributes
plot_df %>% ggplot(aes(Mean_Absolute_Diff, colour = Method)) +
  geom_freqpoly() +
  labs(y = "Frequency", x = "Mean Absolute Difference", colour = "Method")
plot_df %>% ggplot(aes(Mean_DEM_Diff, colour = Method)) +
  geom_freqpoly() +
  labs(y = "Frequency", x = "Mean Signed Difference", colour = "Method")



# good matches
good_matches <- Distances_CR[ which(Distances_CR$Mean_Absolute_Diff <= 0.0152), ] 
#arbitrary cutoff, chosen because it selects 16 lakes
good_match_data <- subset(All_Hypsos, DOW %in% good_matches$DOW)

windows()
ggplot(good_match_data, aes(x = proportion_area, y = -1 * depth_feet, group = Method)) + 
  #750 by 600 pixels works best
  geom_line(aes(linetype = Method, col = Method)) +
  geom_point(aes(shape = Method, col = Method)) + 
  facet_wrap(~DOW, scales = 'free') +
  labs(y = "Depth (feet)", x = "Proportion Area") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hypsographic Curves With High Agreement")

# poor matches
example_match_data <- subset(All_Hypsos, DOW %in% c(71004600, 69004400, 38053900))
unique(example_match_data$DOW)

windows()
ggplot(example_match_data, aes(x = proportion_area, y = -1 * depth_feet, group = Method)) + 
  #750 by 600 pixels works best
  geom_line(aes(linetype = Method, col = Method)) +
  geom_point(aes(shape = Method, col = Method)) + 
  facet_wrap(~DOW, scales = 'free') +
  labs(y = "Depth (feet)", x = "Proportion Area") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hypsographic Curves With Low Agreement")

lake_areas$DOW <- fixlakeid(lake_areas$DOW)
area_diff = merge(Distances, lake_areas)
max_depth = max_depth[,-3]
area_diff = merge(area_diff, max_depth, by = 'DOW')


#change legend attributes
area_diff %>% ggplot(aes(x = max_depth, y = Mean_Absolute_Diff, col = Method), ylab("MAD"), xlab("Max depth (ft)")) +
  geom_point() +
  geom_smooth(method = "lm") + 
  labs(y = "MAD", x = "Max depth (ft)") +
  ggtitle("MAD vs lake depth seperated by method")

area_diff %>% ggplot(aes(x = log(Area_acres), y = Mean_Absolute_Diff, col = Method), ylab("MAD"), xlab("Max depth (ft)")) +
  geom_point() +
  geom_smooth(method = "lm") + 
  labs(y = "MAD", x = "Log Lake Area (Acres)") +
  ggtitle("MAD vs lake size seperated by method")



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




maxdepth <- data.frame("DOW" = character(1), "max_depth" = numeric(1), stringsAsFactors = FALSE)
Hypsos <- rbind(WI_Hypsos, MN_Hypsos)
unique(Hypsos$DOW)
max_depth = Hypsos[1,1]; DOW = Hypsos[1,3]

i = 1
for (row in 1:nrow(Hypsos_depth)) {
  if (DOW != Hypsos[row,3]) {
    maxdepth[i, 1] = DOW
    maxdepth[i, 2] = max_depth
    i = i + 1
    DOW = Hypsos[row,3]
    max_depth = Hypsos[row, 1]
  }
  else {
    max_depth = Hypsos[row, 1]
  }
}
maxdepth[i,1] = DOW; maxdepth[i,2] = max_depth

# statistics for all lake depths
mean(maxdepth$max_depth)  #31.8
confint(maxdepth$max_depth)
se(maxdepth$max_depth)    #1.13
min(maxdepth$max_depth)   #4
max(maxdepth$max_depth)   #443
median(maxdepth$max_depth) #25
hist(maxdepth$max_depth, breaks = 50, xlim = c(0 , 450))
