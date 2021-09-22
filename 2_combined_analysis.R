#load packages
library(dynamicTreeCut)
library(WGCNA)
library(readxl)
library(data.table)
library(geosphere)
library(ggplot2)
library(ggmap)
library(tigris)

#import HI data
hi_data <- read.csv("hi_data.csv", stringsAsFactors = FALSE)

#deal with some commas in the longitude data
hi_data$longitude <- as.numeric(stringr::str_trim(gsub(", ", "", as.character(hi_data$longitude))))

#remove two NA values
hi_data <- hi_data[-c(5282, 5976),] 

#import urbanization perceptions index
upsai <- read.csv("upsai.csv", colClasses = c("character", "character", "character", "numeric", "numeric", "numeric", "character", "numeric", "numeric"))

#import CA data
ca_data <- list()
ca_files <- list.files("ca_data/130 Cali HFs FinchCatcher Outputs/.")
for(i in 1:length(ca_files)){
  ca_data[[i]] <- read_xlsx(paste0("ca_data/130 Cali HFs FinchCatcher Outputs/", ca_files[i]), col_names = FALSE)
}
ca_data <- dplyr::bind_rows(ca_data)

#import CA lat lons
ca_latlon <- read_xlsx("ca_data/California House Finch Project 2012.xlsx")
ca_latlon <- ca_latlon[-1,]
ca_data$latitude <- rep(NA, nrow(ca_data))
ca_data$longitude <- rep(NA, nrow(ca_data))

for(i in 1:nrow(ca_latlon)){
  indices <- which(ca_data$...1 == ca_latlon$`Male ID (parser)`[i])
  if(length(indices) > 0){
    ca_data$latitude[indices] <- ca_latlon$Latitude[i]
    ca_data$longitude[indices] <- ca_latlon$Longitude[i]
  }
}

#remove average slope entirely from both
hi_data <- hi_data[,-13]
ca_data <- ca_data[,-12]

#getting census tract associated with each recording location, and then getting the upsai for it
#CA <- tigris::tracts("CA")
#save(CA, file = "output/CA.RData")
load("output/CA.RData")
CA <- as(CA, "Spatial")
points <- data.frame(long = ca_data$longitude, lat = ca_data$latitude)
sp::coordinates(points) <- ~ long + lat
sp::proj4string(points) <- raster::crs(sp::proj4string(CA))
CA_over <- sp::over(points, CA)

CA_over$upsai <- sapply(1:nrow(CA_over), function(x){upsai$UPSAI_urban[which(upsai$GEOID == CA_over$GEOID[x])]})

#HI <- tigris::tracts("HI")
#save(HI, file = "output/HI.RData")
load("output/HI.RData")
HI <- as(HI, "Spatial")
points <- data.frame(long = hi_data$longitude, lat = hi_data$Latitude)
sp::coordinates(points) <- ~ long + lat
sp::proj4string(points) <- raster::crs(sp::proj4string(HI))
HI_over <- sp::over(points, HI)

HI_over$upsai <- sapply(1:nrow(HI_over), function(x){upsai$UPSAI_urban[which(upsai$GEOID == HI_over$GEOID[x])]})

#construct combined data frame of scaled data from both HI and CA
#add upsai to each data frame
hi_data$upsai <- HI_over$upsai
ca_data$upsai <- CA_over$upsai

#split first column into ind and song
hi_data <- tidyr::unite(data = tidyr::separate(data = hi_data, col = Song.ID, into = c("A", "B", "song")), col = "ind", A:B, sep = " ")

#only individuals with 20 or more songs
hi_data <- hi_data[which(hi_data$ind %in% unique(hi_data$ind)[which(sapply(1:length(unique(hi_data$ind)), function(x){length(unique(hi_data[which(hi_data$ind == unique(hi_data$ind)[x]),]$song))}) >= 20)]), ] #individuals with 20 or more songs

#restructure so columns match each other
ca_data <- ca_data[,c(1, 2, 3, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
hi_data <- hi_data[,c(1, 2, 3, 4, 5, 15, 6, 7, 8, 9, 10, 11, 12, 13, 14)]

#add binary column for whether data from HI (0) or CA (1)
hi_data$hi_ca <- rep(0, nrow(hi_data))
ca_data$hi_ca <- rep(1, nrow(ca_data))

#combine
colnames(ca_data) <- colnames(hi_data)
data <- rbind(hi_data, ca_data)

#add syllable types
load("output/ca_settings.RData")
data$type <- ca_settings
data <- data[-which(data$type == 0),]

#run basic t-tests
#t.test(data[which(data$hi_ca == 0), ]$Start._freq, data[which(data$hi_ca == 1), ]$Start._freq)
#t.test(data[which(data$hi_ca == 0), ]$End_freq, data[which(data$hi_ca == 1), ]$End_freq)
#t.test(data[which(data$hi_ca == 0), ]$Average_freq, data[which(data$hi_ca == 1), ]$Average_freq)
#t.test(data[which(data$hi_ca == 0), ]$Highest_freq, data[which(data$hi_ca == 1), ]$Highest_freq)
#t.test(data[which(data$hi_ca == 0), ]$Lowest_freq, data[which(data$hi_ca == 1), ]$Lowest_freq)
#t.test(data[which(data$hi_ca == 0), ]$Bandwidth, data[which(data$hi_ca == 1), ]$Bandwidth)
#t.test(data[which(data$hi_ca == 0), ]$Duration, data[which(data$hi_ca == 1), ]$Duration)
#t.test(data[which(data$hi_ca == 0), ]$Excursion, data[which(data$hi_ca == 1), ]$Excursion)
#t.test(data[which(data$hi_ca == 0), ]$Concavity, data[which(data$hi_ca == 1), ]$Concavity)

#convert variables to factors
data$hi_ca <- factor(data$hi_ca)
data$type <- factor(data$type)
data$ind <- factor(data$ind)

#check for variance inflation
usdm::vif(data[, c(7, 8, 9, 10, 11, 12, 13, 14, 15)])
usdm::vif(data[, c(11, 13, 14, 15)]) #best choice of variables to include

#check ICCs
performance::icc(lme4::glmer(hi_ca ~ (1|ind), data = data, family = binomial)) #ICC of 1 will lead to convergence issues
performance::icc(lme4::glmer(hi_ca ~ (1|type), data = data, family = binomial)) #ICC of 0.33 makes it a good random effect to include
 
#logistic regression for what parameters predict whether something was recorded in HI or CA
glmm_results <- rstanarm::stan_glmer(hi_ca ~ scale(Lowest_freq) + scale(Duration) + scale(Excursion) + scale(Concavity) + (1|type), data = data, family = binomial, chains = 4, iter = 5000, warmup = 1000, cores = 4)
save(glmm_results, file = "output/glmm_results.RData")
 
#logistic regression for whether unique recording locations in HI tend to be more urban than unique recording locations in CA
glmm_results_noise_1 <- rstanarm::stan_glm(hi_ca ~ scale(upsai), data = data[match(unique(data$ind), data$ind), ], family = binomial, chains = 4, iter = 5000, warmup = 1000, cores = 4)
save(glmm_results_noise_1, file = "output/glmm_results_noise_1.RData")

#check ICC of possible random effects
performance::icc(lme4::lmer(upsai ~ (1|ind), data = data)) #ICC of 1 will lead to convergence issues
performance::icc(lme4::lmer(upsai ~ (1|type), data = data)) #0.082 is relatively low but okay
performance::icc(lme4::lmer(Lowest_freq ~ (1|ind), data = data)) #0.052
performance::icc(lme4::lmer(Lowest_freq ~ (1|type), data = data)) #0.953
 
AIC(lme4::lmer(upsai ~ scale(Lowest_freq) + (1|type), data = data),
    lme4::lmer(Lowest_freq ~ scale(upsai) + (1|ind), data = data),
    lme4::lmer(Lowest_freq ~ scale(upsai) + (1|type), data = data))

#bayesian glmm for whether urbanization predicts lowest frequency of syllables, after controlling for individual and syllable type
glmm_results_noise_2 <- rstanarm::stan_glmer(upsai ~ scale(Lowest_freq) + (1|type), data = data, chains = 4, iter = 5000, warmup = 1000, cores = 4)
save(glmm_results_noise_2, file = "output/glmm_results_noise_2.RData")

#construct data frame of song types data
all_inds <- unique(data$ind)
for(i in 1:length(all_inds)){
  temp <- data[which(data$ind == all_inds[i]),]
  songs <- unique(temp$song)
  sequences <- sapply(1:length(songs), function(x){as.numeric(temp$type[which(temp$song == songs[x])])})
  if(i == 1){
    song_types_data <- data.table(ind = rep(all_inds[i], length(songs)), song = songs, sequence = sequences, hi_ca = temp$hi_ca[1], lat = temp$Latitude[1], lon = temp$longitude[1])
  }
  if(i > 1){
    song_types_data <- rbindlist(list(song_types_data, data.table(ind = rep(all_inds[i], length(songs)), song = songs, sequence = sequences, hi_ca = temp$hi_ca[1], lat = temp$Latitude[1], lon = temp$longitude[1])))
  }
  rm(list = c("temp", "songs", "sequences"))
}

#calculate distances
dist_function <- function(x, y){((stringdist::seq_dist(x, y, method = "lv") - abs(length(x)-length(y))) / min(length(x), length(y)))}
distances <- proxy::dist(song_types_data$sequence, method = dist_function)
save(distances, file = "output/distances.RData")

#convert to similarities and set threshold
similarities <- 1 - as.matrix(distances, nrow = nrow(song_types_data))
threshold <- 0.75

#get unique individuals from CA to loop through
inds <- unique(song_types_data$ind[which(song_types_data$hi_ca == 1)])

#for each individual from CA, who from HI do they share a song type with? and what song in the HI birds' repertoire?
shared_songs <- list()
for(i in 1:length(inds)){
  #get the number songs sung by individual from CA for the remainder/quotient calculation
  n_songs <- length(which(song_types_data$ind == inds[i]))
  
  #subset the similarities matrix by the HI rows and CA columns for that particular individual
  #ID the row and column of values over 0.75 using remainder/quotient and store as the HI individual and their song that matches a song of the CA bird
  shared_songs[[i]] <- data.frame(ind_ca = as.character(rep(inds[i], nrow(which(similarities[which(song_types_data$hi_ca == 0), which(song_types_data$ind == inds[i])] >= threshold, arr.ind = TRUE)))),
                                 song_ca = song_types_data$song[song_types_data$ind == inds[i]][as.data.frame(which(similarities[which(song_types_data$hi_ca == 0), which(song_types_data$ind == inds[i])] >= threshold, arr.ind = TRUE))$col],
                                 ind_hi = as.character(song_types_data$ind[as.data.frame(which(similarities[which(song_types_data$hi_ca == 0), which(song_types_data$ind == inds[i])] >= threshold, arr.ind = TRUE))$row]),
                                 song_hi = song_types_data$song[as.data.frame(which(similarities[which(song_types_data$hi_ca == 0), which(song_types_data$ind == inds[i])] >= threshold, arr.ind = TRUE))$row])
}

#save as XLSX file for manual checks
shared_songs <- as.data.frame(rbindlist(shared_songs))
writexl::write_xlsx(shared_songs, "output/shared_songs.xlsx")

#add to song types data
song_types_data$shared <- rep(0, nrow(song_types_data))
for(i in 1:nrow(shared_songs)){
  song_types_data$shared[which(song_types_data$ind == shared_songs[i,]$ind_ca & song_types_data$song == shared_songs[i,]$song_ca)] <- 1
  song_types_data$shared[which(song_types_data$ind == shared_songs[i,]$ind_hi & song_types_data$song == shared_songs[i,]$song_hi)] <- 1
}
song_types_data$shared[which(song_types_data$type %in% shared_songs)] <- 1
