
## Functions ####
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
calculated_data <- rbind(unique(abal_mod_dataset[, c("site_code", "species", "C.N_FH", "pH_L1")]),
                         unique(pcab_mod_dataset[, c("site_code", "species", "C.N_FH", "pH_L1")]),
                         unique(pisy_mod_dataset[, c("site_code", "species", "C.N_FH", "pH_L1")]),
                         unique(fasy_mod_dataset[, c("site_code", "species", "C.N_FH", "pH_L1")]),
                         unique(qusp_mod_dataset[, c("site_code", "species", "C.N_FH", "pH_L1")]))

calculated_data$x <- 0
calculated_data$x[which(calculated_data$species == "ABAL")] <- 1
calculated_data$x[which(calculated_data$species == "PCAB")] <- 2
calculated_data$x[which(calculated_data$species == "PISY")] <- 3
calculated_data$x[which(calculated_data$species == "FASY")] <- 4
calculated_data$x[which(calculated_data$species == "QUSP")] <- 5

## Figure makng ####
g_ph <- ggplot(calculated_data)
g_ph <- g_ph + geom_boxplot(aes(x = x, y = pH_L1, fill = species), alpha = 0.8)
g_ph <- g_ph + scale_fill_manual(breaks=names(cols.species), values=cols.species)
g_ph <- g_ph + scale_x_continuous("", limits = c(0.5, 5.5), breaks = c(1:5), labels = c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
g_ph <- g_ph + scale_y_continuous("pH", limits = c(3.5, 6.5), breaks = seq(0, 10, 0.5), labels = formatC(seq(0, 10, 0.5), format = "f", digits = 1))
g_ph <- my.theme(g_ph)

g_cn <- ggplot(calculated_data)
g_cn <- g_cn + geom_boxplot(aes(x = x, y = C.N_FH, fill = species), alpha = 0.8)
g_cn <- g_cn + scale_fill_manual(breaks=names(cols.species), values=cols.species)
g_cn <- g_cn + scale_x_continuous("", limits = c(0.5, 5.5), breaks = c(1:5), labels = c("ABAL", "PCAB", "PISY", "FASY", "QUSP"))
g_cn <- g_cn + scale_y_continuous("C:N", limits = c(15, 35), breaks = seq(0, 100, 2.5), labels = formatC(seq(0, 100, 2.5), format = "f", digits = 1))
g_cn <- my.theme(g_cn)

figure <- ggarrange(g_ph, g_cn, ncol = 2, nrow = 1, align = "hv", common.legend = T, legend = "bottom")
