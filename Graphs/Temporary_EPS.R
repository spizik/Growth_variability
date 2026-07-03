


rwi.statistics <- read.table("Calculated_datasets/rwi_stats/rwi_statistics.txt", dec = ".", sep = ";", header = T)






boxplot(eps~species, rwi.statistics)

boxplot(eps~elev_zone, subset(rwi.statistics, species == "ABAL"))
boxplot(eps~elev_zone, subset(rwi.statistics, species == "PCAB"))
boxplot(eps~elev_zone, subset(rwi.statistics, species == "PISY"))
boxplot(eps~elev_zone, subset(rwi.statistics, species == "FASY"))
boxplot(eps~elev_zone, subset(rwi.statistics, species == "QUSP"))

