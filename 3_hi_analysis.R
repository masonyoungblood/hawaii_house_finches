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
hi_data$longitude <- as.numeric(stringr::str_trim(gsub(", ", "", as.character(hi_data$longitude)))) #clean lat/lons
hi_data <- hi_data[-c(5282, 5976),] #removing two NA
hi_data_original <- hi_data #store unscaled version for later
hi_data <- hi_data[, c(5, 6, 7, 8, 9, 10, 11, 12, 14)] #reordering

#scale and center variables prior to clustering
hi_data$Start._freq <- scale(hi_data$Start._freq)
hi_data$End_freq <- scale(hi_data$End_freq)
hi_data$Average_freq <- scale(hi_data$Average_freq)
hi_data$Highest_freq <- scale(hi_data$Highest_freq)
hi_data$Lowest_freq <- scale(hi_data$Lowest_freq)
hi_data$Bandwidth <- scale(hi_data$Bandwidth)
hi_data$Duration <- scale(hi_data$Duration)
hi_data$Excursion <- scale(hi_data$Excursion)
hi_data$Concavity <- scale(hi_data$Concavity)

#set random seed
set.seed(123)

#clustering
dist_matrix <- dist(hi_data)
clustering <- fastcluster::hclust(dist_matrix, method = "average")

#dynamic tree cut with settings from CA study
ca_settings <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, maxCoreScatter = 1, minGap = 0.5, cutHeight = 2)
dunn_ca <- clValid::dunn(dist_matrix, as.numeric(ca_settings))
conn_ca <- clValid::connectivity(dist_matrix, as.numeric(ca_settings))

#dynamic tree cut with deep split of 0
hybrid_cut_deep0 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 0)
dunn_deep0 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep0))
conn_deep0 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep0))

#dynamic tree cut with deep split of 1
hybrid_cut_deep1 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 1)
dunn_deep1 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep1))
conn_deep1 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep1))

#dynamic tree cut with deep split of 2
hybrid_cut_deep2 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 2)
dunn_deep2 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep2))
conn_deep2 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep2))

#dynamic tree cut with deep split of 3
hybrid_cut_deep3 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 3)
dunn_deep3 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep3))
conn_deep3 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep3))

#dynamic tree cut with deep split of 4
hybrid_cut_deep4 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 4)
dunn_deep4 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep4))
conn_deep4 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep4))

#compare performance of different methods
cluster_indices <- data.frame(dunn_index = c(dunn_eric, dunn_deep0, dunn_deep1, dunn_deep2, dunn_deep3, dunn_deep4),
                              conn_index = c(conn_eric, conn_deep0, conn_deep1, conn_deep2, conn_deep3, conn_deep4))

#add syllables to unscaled data and overwrite
hi_data <- cbind(hi_data_original, data.frame(type = ca_settings, syl_cols = labels2colors(ca_settings)))

#import urbanization perceptions index
upsai <- read.csv("upsai.csv", colClasses = c("character", "character", "character", "numeric", "numeric", "numeric", "character", "numeric", "numeric"))

#load in HI spatial data
load("output/HI.RData")
HI <- as(HI, "Spatial")
points <- data.frame(long = hi_data$longitude, lat = hi_data$Latitude)
sp::coordinates(points) <- ~ long + lat
sp::proj4string(points) <- raster::crs(sp::proj4string(HI))
HI_over <- sp::over(points, HI)

#add upsai to HI data
HI_over$upsai <- sapply(1:nrow(HI_over), function(x){upsai$UPSAI_urban[which(upsai$GEOID == HI_over$GEOID[x])]})
hi_data$upsai <- HI_over$upsai

#split first column into ind and song
hi_data <- tidyr::unite(data = tidyr::separate(data = hi_data, col = Song.ID, into = c("A", "B", "song")), col = "ind", A:B, sep = " ")

#remove unassigned syllables (labelled with 0)
hi_data <- hi_data[-which(hi_data$type == 0), ]

#check ICC of possible random effects
performance::icc(lme4::lmer(upsai ~ (1|ind), data = hi_data)) #ICC of 1 will lead to convergence issues
performance::icc(lme4::lmer(upsai ~ (1|type), data = hi_data)) #0.092 is relatively low but okay
performance::icc(lme4::lmer(Lowest_freq ~ (1|ind), data = hi_data)) #0.054
performance::icc(lme4::lmer(Lowest_freq ~ (1|type), data = hi_data)) #0.956

AIC(lme4::lmer(upsai ~ scale(Lowest_freq) + (1|type), data = hi_data),
    lme4::lmer(Lowest_freq ~ scale(upsai) + (1|ind), data = hi_data),
    lme4::lmer(Lowest_freq ~ scale(upsai) + (1|type), data = hi_data))

#bayesian glmm for whether urbanization predicts lowest frequency of syllables, after controlling for syllable type
glmm_results_noise_2_hi <- rstanarm::stan_glmer(upsai ~ scale(Lowest_freq) + (1|type), data = hi_data, chains = 4, iter = 5000, warmup = 1000, cores = 4)
save(glmm_results_noise_2_hi, file = "output/glmm_results_noise_2_hi.RData")

#get repertoire data for each unique individual
syl_list <- list()
lat <- c()
lon <- c()
for(i in 1:length(unique(hi_data$ind))){
  syl_list[[i]] <- unique(hi_data$type[which(hi_data$ind == unique(hi_data$ind)[i])])
  lat <- c(lat, hi_data$Latitude[which(hi_data$ind == unique(hi_data$ind)[i])][1])
  lon <- c(lon, hi_data$longitude[which(hi_data$ind == unique(hi_data$ind)[i])][1])
}
reps <- data.table(ind = unique(hi_data$ind), lat = lat, lon = lon, rep = syl_list)

#get jaccard and distance matrices
jacc_matrix <- matrix(NA, nrow = nrow(reps), ncol = nrow(reps))
rownames(jacc_matrix) <- reps$ind
colnames(jacc_matrix) <- reps$ind
dist_matrix <- matrix(NA, nrow = nrow(reps), ncol = nrow(reps))
rownames(dist_matrix) <- reps$ind
colnames(dist_matrix) <- reps$ind
for(i in 1:nrow(reps)){
  for(j in (1:nrow(reps))[-i]){
    a <- intersect(reps$rep[[i]], reps$rep[[j]]) #shared syllables
    b <- reps$rep[[i]][-which(reps$rep[[i]] %in% a)] #unique syllables in first individual
    c <- reps$rep[[j]][-which(reps$rep[[j]] %in% a)] #unique syllables in second individual
    d <- abs(length(reps$rep[[i]])-length(reps$rep[[j]])) #absolutely difference in rep size
    jacc_matrix[i, j] <- length(a)/(length(a)+length(b)+length(c)-d)
    dist_matrix[i, j] <- distm(c(reps$lon[[i]], reps$lat[[i]]), c(reps$lon[[j]], reps$lat[[j]]))
    rm(list = c("a", "b", "c", "d"))
  }
}

#set diagonal to 0
diag(jacc_matrix) <- 0
diag(dist_matrix) <- 0

#store the upper triangles of both matrices in a data frame
model_data <- data.frame(x = dist_matrix[upper.tri(dist_matrix, diag = FALSE)]/1000, y = jacc_matrix[upper.tri(jacc_matrix, diag = FALSE)])

#set random seed, and split data into training and testing data
set.seed(12345)
train_indices <- sample(nrow(model_data), nrow(model_data)*0.8)
train_model_data <- model_data[train_indices,]
test_model_data <- model_data[-train_indices,]

#run different models
lin_model <- lm(y ~ x, data = train_model_data) #linear model
log_model <- lm(y ~ log(x), data = train_model_data) #logarithmic model
exp_model <- nls(y ~ SSasymp(x, yf, y0, log_alpha), data = train_model_data) #exponential
gam_model <- mgcv::gam(y ~ x, data = train_model_data) #generalized additive model

#compare model fit
model_fitting <- data.frame(AIC = AIC(lin_model, log_model, exp_model, gam_model)$AIC,
                            RMSE = c(caret::RMSE(predict(lin_model, test_model_data), test_model_data$y),
                                     caret::RMSE(predict(log_model, test_model_data), test_model_data$y),
                                     caret::RMSE(predict(exp_model, test_model_data), test_model_data$y),
                                     caret::RMSE(predict(gam_model, test_model_data), test_model_data$y)),
                            R2 = c(caret::R2(predict(lin_model, test_model_data), test_model_data$y),
                                   caret::R2(predict(log_model, test_model_data), test_model_data$y),
                                   caret::R2(predict(exp_model, test_model_data), test_model_data$y),
                                   caret::R2(predict(gam_model, test_model_data), test_model_data$y)))
rownames(model_fitting) <- c("linear", "logarithmic", "exponential", "general additive")

#re-run both linear and best fitting model on full dataset
lin_model_full <- lm(y ~ x, data = model_data) #linear model
log_model_full <- lm(y ~ log(x), data = model_data) #logarithmic model

#construct data frame of song types
all_inds <- unique(hi_data$ind)
for(i in 1:length(all_inds)){
  temp <- hi_data[which(hi_data$ind == all_inds[i]),]
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

#calculate similarities and set threshold
similarities <- 1 - as.matrix(distances, nrow = nrow(song_types_data))
threshold <- 0.75

#get rid of self-similarity of song types
diag(similarities) <- 0

#get distance and jaccard matrices
inds <- unique(song_types_data$ind)
jacc_matrix <- matrix(NA, nrow = length(inds), ncol = length(inds))
rownames(jacc_matrix) <- inds
colnames(jacc_matrix) <- inds
dist_matrix <- matrix(NA, nrow = length(inds), ncol = length(inds))
rownames(dist_matrix) <- inds
colnames(dist_matrix) <- inds
for(i in 1:length(inds)){
  for(j in 1:length(inds)){
    if(i != j){
      temp <- similarities[which(song_types_data$ind == inds[i]), which(song_types_data$ind == inds[j])]
      jacc_matrix[i, j] <- length(which(temp >= threshold))/(nrow(temp) + ncol(temp))
      dist_matrix[i, j] <- distm(c(song_types_data$lon[which(song_types_data$ind == inds[i])[1]], 
                                   song_types_data$lat[which(song_types_data$ind == inds[i])[1]]),
                                 c(song_types_data$lon[which(song_types_data$ind == inds[j])[1]],
                                   song_types_data$lat[which(song_types_data$ind == inds[j])[1]]))
      rm(temp)
    }
  }
}

#store the upper triangles of both matrices in a data frame
song_model_data <- data.frame(x = dist_matrix[upper.tri(dist_matrix, diag = FALSE)]/1000, y = jacc_matrix[upper.tri(jacc_matrix, diag = FALSE)])
