## Functions
## Calculations
# Computes the mean and standard deviation of relative TRW variability (variance normalized by mean TRW) for a given species.
# Prints the two summary statistics as percentages to quantify species-level variability intensity.
calculate <- function(sp){
  
  ## Testing agruments
  # sp="PCAB"
  
  sub <- subset(eco_relevance_agg, species == sp)
  print(c(round(mean(sub$TRW_variance / sub$mid_TRW)*100,1),
          round(sd(sub$TRW_variance / sub$mid_TRW)*100,1)))
}

## Figure functions
my.theme<-function(graph,legend.pos="bottom"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.x = element_line(colour="black"),
                     axis.text.x = element_text(colour="black"),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
my.theme.nox<-function(graph,legend.pos="bottom"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x=element_blank(),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks.y = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}

## Calculations ####
eco_relevance <- clim.dataset[,c("site_code", "species", "year", "cv_TRW", "mid_TRW")]
eco_relevance$TRW_variance <- eco_relevance$cv_TRW * eco_relevance$mid_TRW
eco_relevance <- eco_relevance[-which(eco_relevance$TRW_variance >5),]

eco_relevance_agg <- aggregate(mid_TRW ~ site_code + species, eco_relevance, mean) 
eco_relevance_agg <- merge(eco_relevance_agg, aggregate(TRW_variance ~ site_code + species, eco_relevance, mean)[,c("site_code","TRW_variance")],
                            by = "site_code")

eco_relevance_agg$Growth_variance <- eco_relevance_agg$TRW_variance / eco_relevance_agg$mid_TRW

res_kw_trw <- dunnTest(mid_TRW ~ species, data = eco_relevance_agg)
res_kw_trw$res$P.adj <- round(res_kw_trw$res$P.adj, 3)

res_kw_var <- dunnTest(TRW_variance ~ species, data = eco_relevance_agg)
res_kw_var$res$P.adj <- round(res_kw_var$res$P.adj, 3)

res_kw_var2 <- dunnTest(Growth_variance ~ species, data = eco_relevance_agg)
res_kw_var2$res$P.adj <- round(res_kw_var2$res$P.adj, 3)
               
## Testing ####
print("----------------------------------------------------")
print("----------------------------------------------------")
print("---------------- Differences in TRW ----------------")        
print(res_kw_trw)  
print("----------------------------------------------------")
print("----------------------------------------------------")
print("------------ Differences in TRW variabily ----------")    
print(res_kw_var)
print("----------------------------------------------------")
print("----------------------------------------------------")
print("------------------ How much is it? -----------------") 
print("ABAL")
calculate("ABAL")
print("PCAB")
calculate("PCAB")
print("PISY")
calculate("PISY")
print("FASY")
calculate("FASY")
print("QUSP")
calculate("QUSP")  
print("----------------------------------------------------")
print("----------------------------------------------------")
print("--------------- Differences, how much --------------")    
print(res_kw_var2)
print("----------------------------------------------------")
print("----------------------------------------------------")

## Figure making ####

eco_relevance_agg$x <- 0
eco_relevance_agg$x[which(eco_relevance_agg$species == "ABAL")] <- 1
eco_relevance_agg$x[which(eco_relevance_agg$species == "PCAB")] <- 2
eco_relevance_agg$x[which(eco_relevance_agg$species == "PISY")] <- 3
eco_relevance_agg$x[which(eco_relevance_agg$species == "FASY")] <- 4
eco_relevance_agg$x[which(eco_relevance_agg$species == "QUSP")] <- 5

g_trw <- ggplot(eco_relevance_agg)
g_trw <- g_trw + geom_boxplot(aes(x = x, y = mid_TRW, fill = species), alpha = 0.8)
g_trw <- g_trw + scale_fill_manual(breaks=names(cols.species), values=cols.species)
g_trw <- g_trw + scale_x_continuous("", limits = c(0.5, 5.5), breaks = c(1:5), labels = c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
g_trw <- g_trw + scale_y_continuous("mean TRW", limits = c(0, 3.5), breaks = seq(0, 10, 0.5), labels = formatC(seq(0, 10, 0.5), format = "f", digits = 1))
g_trw <- my.theme(g_trw)

g_var <- ggplot(eco_relevance_agg)
g_var <- g_var + geom_boxplot(aes(x = x, y = TRW_variance, fill = species), alpha = 0.8)
g_var <- g_var + scale_fill_manual(breaks=names(cols.species), values=cols.species)
g_var <- g_var + scale_x_continuous("", limits = c(0.5, 5.5), breaks = c(1:5), labels = c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
g_var <- g_var + scale_y_continuous("TRW variability", limits = c(0, 1.5), breaks = seq(0, 10, 0.25), labels = formatC(seq(0, 10, 0.25), format = "f", digits = 2))
g_var <- my.theme(g_var)

g_var2 <- ggplot(eco_relevance_agg)
g_var2 <- g_var2 + geom_boxplot(aes(x = x, y = Growth_variance, fill = species), alpha = 0.8)
g_var2 <- g_var2 + scale_fill_manual(breaks=names(cols.species), values=cols.species)
g_var2 <- g_var2 + scale_x_continuous("", limits = c(0.5, 5.5), breaks = c(1:5), labels = c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
g_var2 <- g_var2 + scale_y_continuous("Percentage of mean growth", limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20))
g_var2 <- my.theme(g_var2)

figure <- ggarrange(g_trw, g_var, g_var2, ncol = 1, nrow = 3, align = "hv", common.legend = T, legend = "bottom")
