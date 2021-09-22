#import features table and subset for variables of interest
hi_data <- read.csv("hi_data.csv", stringsAsFactors = FALSE)

#import CA data
ca_data <- list()
ca_files <- list.files("ca_data/130 Cali HFs FinchCatcher Outputs/.")
for(i in 1:length(ca_files)){
  ca_data[[i]] <- readxl::read_xlsx(paste0("ca_data/130 Cali HFs FinchCatcher Outputs/", ca_files[i]), col_names = FALSE)
}
ca_data <- dplyr::bind_rows(ca_data)

#remove two NA values from HI data
hi_data <- hi_data[-c(5282, 5976),]

#construct combined data frame of scaled data from both HI and CA
#split first column into ind and song
hi_data <- tidyr::unite(data = tidyr::separate(data = hi_data, col = Song.ID, into = c("A", "B", "song")), col = "ind", A:B, sep = " ")

#only individuals with 20 or more songs
hi_data <- hi_data[which(hi_data$ind %in% unique(hi_data$ind)[which(sapply(1:length(unique(hi_data$ind)), function(x){length(unique(hi_data[which(hi_data$ind == unique(hi_data$ind)[x]),]$song))}) >= 20)]), ]

#restructure so columns match each other
hi_data <- hi_data[, c(6, 7, 8, 9, 10, 11, 12, 13, 15)]
ca_data <- as.data.frame(ca_data[, c(4, 5, 6, 7, 8, 9, 10, 11, 13)])

#add binary column for whether data from HI (0) or CA (1)
hi_data$hi_ca <- rep(0, nrow(hi_data))
ca_data$hi_ca <- rep(1, nrow(ca_data))

#combine
colnames(ca_data) <- colnames(hi_data)
data <- rbind(hi_data, ca_data)

#prepare data for clustering data
data$hi_ca <- factor(data$hi_ca)
data$Start._freq <- scale(data$Start._freq)
data$End_freq <- scale(data$End_freq)
data$Average_freq <- scale(data$Average_freq)
data$Highest_freq <- scale(data$Highest_freq)
data$Lowest_freq <- scale(data$Lowest_freq)
data$Bandwidth <- scale(data$Bandwidth)
data$Duration <- scale(data$Duration)
data$Excursion <- scale(data$Excursion)
data$Concavity <- scale(data$Concavity)

#set random seed
set.seed(123)

#construct distance matrix
dist_matrix <- stats::dist(data[,1:9])

#cluster and save
clustering <- fastcluster::hclust(dist_matrix, method = "average")
save(clustering, file = "output/clustering.RData")

#dynamic tree cut using settings from the CA study
ca_settings <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, maxCoreScatter = 1, minGap = 0.5, cutHeight = 2)
dunn_ca <- clValid::dunn(dist_matrix, as.numeric(ca_settings))
conn_ca <- clValid::connectivity(dist_matrix, as.numeric(ca_settings))
save(ca_settings, file = "output/ca_settings.RData")

#dynamic tree cut with deep split of 0
hybrid_cut_deep0 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 0)
dunn_deep0 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep0))
conn_deep0 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep0))
save(hybrid_cut_deep0, file = "output/hybrid_cut_deep0.RData")

#dynamic tree cut with deep split of 1
hybrid_cut_deep1 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 1)
dunn_deep1 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep1))
conn_deep1 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep1))
save(hybrid_cut_deep1, file = "output/hybrid_cut_deep1.RData")

#dynamic tree cut with deep split of 2
hybrid_cut_deep2 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 2)
dunn_deep2 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep2))
conn_deep2 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep2))
save(hybrid_cut_deep2, file = "output/hybrid_cut_deep2.RData")

#dynamic tree cut with deep split of 3
hybrid_cut_deep3 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 3)
dunn_deep3 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep3))
conn_deep3 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep3))
save(hybrid_cut_deep3, file = "output/hybrid_cut_deep3.RData")

#dynamic tree cut with deep split of 4
hybrid_cut_deep4 <- dynamicTreeCut::cutreeDynamic(dendro = clustering, distM = as.matrix(dist_matrix), method = "hybrid", minClusterSize = 5, deepSplit = 4)
dunn_deep4 <- clValid::dunn(dist_matrix, as.numeric(hybrid_cut_deep4))
conn_deep4 <- clValid::connectivity(dist_matrix, as.numeric(hybrid_cut_deep4))
save(hybrid_cut_deep4, file = "output/hybrid_cut_deep4.RData")

#compare performance of different methods with clustering indices
cluster_indices <- data.frame(dunn_index = c(dunn_ca, dunn_deep0, dunn_deep1, dunn_deep2, dunn_deep3, dunn_deep4),
                              conn_index = c(conn_ca, conn_deep0, conn_deep1, conn_deep2, conn_deep3, conn_deep4))
rownames(cluster_indices) <- c("ca", "deep0", "deep1", "deep2", "deep3", "deep4")
save(cluster_indices, file = "output/cluster_indices.RData")
