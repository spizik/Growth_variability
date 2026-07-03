## Functions ####
## Calculations and data preparations

# Prepares a climate–growth dataset for a given species by merging CRN data with site-level climate and pollution variables across a selected period. 
# Computes yearly ranges (max–min) for all predictors, scales all variables, and returns a clean analysis-ready data frame.
prepare.data <- function(crn, dta, sp, period=c(1961:2017)){
  
  ## Testing arguments
  # crn = df.crn
  # dta = abal_mod_dataset
  # sp = "ABAL"
  # period = c(1961:2017)
  
  crn<-subset(df.crn, species==sp)
  clim<-aggregate(mean_temp~year,data=dta,mean)
  clim$cwb<-merge(clim,aggregate(mean_cwb~year,data=dta,mean),by.x="year",by.y="year")$mean_cwb
  clim$sox<-merge(clim,aggregate(sox~year,data=dta,mean),by.x="year",by.y="year")$sox
  clim$nox<-merge(clim,aggregate(nox~year,data=dta,mean),by.x="year",by.y="year")$nox
  clim$std<-merge(clim,crn,by.x="year",by.y="year")$cv.std.mid 
  
  clim <- clim[which(clim$year %in% period),]
  dta <- dta[which(dta$year %in% period),]
  
  temp<-aggregate(mean_temp~year,data=dta,max)
  clim$temp_range<-temp$mean_temp-aggregate(mean_temp~year,data=dta,min)$mean_temp
  temp<-aggregate(mean_cwb~year,data=dta,max)
  clim$cwb_range<-temp$mean_cwb-aggregate(mean_cwb~year,data=dta,min)$mean_cwb
  temp<-aggregate(sox~year,data=dta,max)
  clim$sox_range<-temp$sox-aggregate(sox~year,data=dta,min)$sox
  temp<-aggregate(nox~year,data=dta,max)
  clim$nox_range<-temp$nox-aggregate(nox~year,data=dta,min)$nox
  
  temp<-aggregate(mid_TRW~year,data=dta,max)
  clim$trw<-temp$mid_TRW-aggregate(mid_TRW~year,data=dta,min)$mid_TRW
  
  clim$mean_temp     <- scale(clim$mean_temp)
  clim$temp_range    <- scale(clim$temp_range)
  clim$cwb           <- scale(clim$cwb)
  clim$cwb_range     <- scale(clim$cwb_range)
  clim$sox           <- scale(clim$sox)
  clim$sox_range     <- scale(clim$sox_range)
  clim$nox           <- scale(clim$nox)
  clim$nox_range     <- scale(clim$nox_range)
  
  return(clim)
}

# Fits a GAM model predicting inter-site variability using tensor-product smooths of climate and pollution variables and their interactions. 
# Returns the fitted GAM object for downstream interpretation and variable-importance analysis.
calc.model <- function(input){
  
  ## Testing arguments
  # input = prepare.data(df.crn, abal_mod_dataset, "ABAL", c(1961:2017))
  
  model.gam <- gam(std ~ 
                     te(mean_temp, temp_range, k=3) +
                     te(cwb, cwb_range, k=3) +
                     te(mean_temp, cwb, k=3) +
                     te(sox, sox_range, k=3)
                   ,
                   data = input,
                   method = "REML")
  
  return(model.gam)
}

# Computes permutation-based variable importance (via iml::FeatureImp) for the fitted GAM model using MSE loss. 
# Returns a data frame of variable importances with an added column of relative contributions.
calc.importance <- function(input){
  
  ## Testing arguments
  # input = calc.model(prepare.data(df.crn, abal_mod_dataset, "ABAL", c(1961:2017)))
  
  X <- model.frame(input)[, -1] 
  y <- input$y
  
  # Predictor objekt
  pred <- Predictor$new(
    model = input,
    data = X,
    y = y,
    predict.fun = function(object, newdata) predict(object, newdata = newdata)
  )
  
  # VCalculate importance
  imp <- FeatureImp$new(pred, loss = "mse")
  
  vi_df <- as.data.frame(imp$results) 
  vi_df$perc <- vi_df$importance  / sum(vi_df$importance )
  
  return(vi_df)
}

## Graph functions
my.theme <- function(graph,legend.pos="bottom"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.x = element_line(colour="black"),
                     axis.text.x = element_text(colour="black"),
                     axis.line.y = element_line(colour="black"),
                     axis.text.y = element_text(colour="black"),
                     axis.ticks = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
barplot.primary <- function(crn, dta, sp, period){
  prepared_data <- prepare.data(crn, dta, sp, period)
  calc_model <- calc.model(prepared_data)
  
  dev_expl <- summary(calc_model)$dev.expl
  calc_imp <- calc.importance(calc_model)
  
  calc_imp$feature2 <- ""
  calc_imp$feature2[which(calc_imp$feature == "mean_temp")] <- "f"
  calc_imp$feature2[which(calc_imp$feature == "temp_range")] <- "e"
  calc_imp$feature2[which(calc_imp$feature == "cwb")] <- "d"
  calc_imp$feature2[which(calc_imp$feature == "cwb_range")] <- "c"
  calc_imp$feature2[which(calc_imp$feature == "sox")] <- "b"
  calc_imp$feature2[which(calc_imp$feature == "sox_range")] <- "a"
  
  custom_palette <- c(
    "a" = "#bdbdbd",  # sox_range – světle šedá
    "b" = "#636363",  # sox       – tmavší šedá
    "c" = "#9ecae1",  # cwb_range – světle modrá
    "d" = "#08519c",  # cwb       – tmavě modrá
    "e" = "#fcae91",  # temp_range – světle červená
    "f" = "#cb181d"   # mean_temp – tmavší červená
  )
  
  g <- ggplot(calc_imp)
  g <- g + geom_bar(aes(x = 0.5, y = perc, fill = feature2), stat = "identity")
  g <- g + scale_fill_manual(values = custom_palette, breaks = names(custom_palette))
  g <- g + coord_flip()
  g <- g + annotate("text", x = 0.5, y = 1.1, label = paste0(round(dev_expl*100,1)," %"), size = 8)
  g <- g + scale_x_continuous("", limits = c(0,1), breaks = 0.5, labels = sp) 
  g <- g + scale_y_continuous("", limits = c(0,1.25), breaks = seq(0,1,0.2), labels = formatC(seq(0,1,0.2), format = "f", digits = 1))
  g <- my.theme(g)
  g
}
barplot.secondary <- function(crn, dta, sp, period){
  prepared_data <- prepare.data(crn, dta, sp, period)
  calc_model <- calc.model(prepared_data)
  
  dev_expl <- summary(calc_model)$dev.expl
  calc_imp <- calc.importance(calc_model)
  
  calc_imp$feature2 <- ""
  calc_imp$feature2[which(calc_imp$feature == "mean_temp")] <- "f"
  calc_imp$feature2[which(calc_imp$feature == "temp_range")] <- "e"
  calc_imp$feature2[which(calc_imp$feature == "cwb")] <- "d"
  calc_imp$feature2[which(calc_imp$feature == "cwb_range")] <- "c"
  calc_imp$feature2[which(calc_imp$feature == "sox")] <- "b"
  calc_imp$feature2[which(calc_imp$feature == "sox_range")] <- "a"
  
  custom_palette <- c(
    "a" = "#bdbdbd",  # sox_range – světle šedá
    "b" = "#636363",  # sox       – tmavší šedá
    "c" = "#9ecae1",  # cwb_range – světle modrá
    "d" = "#08519c",  # cwb       – tmavě modrá
    "e" = "#fcae91",  # temp_range – světle červená
    "f" = "#cb181d"   # mean_temp – tmavší červená
  )
  
  g <- ggplot(calc_imp)
  g <- g + geom_bar(aes(x = 0.25, y = perc, fill = feature2), stat = "identity")
  g <- g + scale_fill_manual(values = custom_palette, breaks = names(custom_palette))
  g <- g + coord_flip()
  g <- g + annotate("text", x = 0.25, y = 1.1, label = paste0(round(dev_expl*100,1)," %"), size = 4)
  g <- g + scale_x_continuous("", limits = c(0,5), breaks = 0.5, labels = c("")) 
  g <- g + scale_y_continuous("", limits = c(0,1.25), breaks = seq(0,1,0.2), labels = formatC(seq(0,1,0.2), format = "f", digits = 1))
  g <- my.theme(g)
  g
}

## Calculations ####
dataset <- rbind(calc.importance(calc.model(prepare.data(df.crn, abal_mod_dataset, "ABAL", c(1961:2017)))),
                 calc.importance(calc.model(prepare.data(df.crn, pcab_mod_dataset, "PCAB", c(1961:2017)))),
                 calc.importance(calc.model(prepare.data(df.crn, pisy_mod_dataset, "PISY", c(1961:2017)))),
                 calc.importance(calc.model(prepare.data(df.crn, fasy_mod_dataset, "FASY", c(1961:2017)))),
                 calc.importance(calc.model(prepare.data(df.crn, qusp_mod_dataset, "QUSP", c(1961:2017)))))

deviances_expl <- rbind(summary(calc.model(prepare.data(df.crn, abal_mod_dataset, "ABAL", c(1961:2017))))$dev.expl,
                        summary(calc.model(prepare.data(df.crn, pcab_mod_dataset, "PCAB", c(1961:2017))))$dev.expl,
                        summary(calc.model(prepare.data(df.crn, pisy_mod_dataset, "PISY", c(1961:2017))))$dev.expl,
                        summary(calc.model(prepare.data(df.crn, fasy_mod_dataset, "FASY", c(1961:2017))))$dev.expl,
                        summary(calc.model(prepare.data(df.crn, qusp_mod_dataset, "QUSP", c(1961:2017))))$dev.expl)

dataset$species <- rep(c("ABAL", "PCAB", "PISY", "FASY", "QUSP"), each = 6)

deviances_expl <- as.data.frame(deviances_expl)
deviances_expl$species <- c("ABAL", "PCAB", "PISY", "FASY", "QUSP")

dataset$x[which(dataset$species=="ABAL")]<-5
dataset$x[which(dataset$species=="PCAB")]<-4
dataset$x[which(dataset$species=="PISY")]<-3
dataset$x[which(dataset$species=="FASY")]<-2
dataset$x[which(dataset$species=="QUSP")]<-1

deviances_expl$x[which(deviances_expl$species=="ABAL")]<-5
deviances_expl$x[which(deviances_expl$species=="PCAB")]<-4
deviances_expl$x[which(deviances_expl$species=="PISY")]<-3
deviances_expl$x[which(deviances_expl$species=="FASY")]<-2
deviances_expl$x[which(deviances_expl$species=="QUSP")]<-1

dataset$feature2 <- ""
dataset$feature2[which(dataset$feature == "mean_temp")] <- "f"
dataset$feature2[which(dataset$feature == "temp_range")] <- "e"
dataset$feature2[which(dataset$feature == "cwb")] <- "d"
dataset$feature2[which(dataset$feature == "cwb_range")] <- "c"
dataset$feature2[which(dataset$feature == "sox")] <- "b"
dataset$feature2[which(dataset$feature == "sox_range")] <- "a"

dataset$feature3 <- "pollusions"
dataset$feature3[which(dataset$feature %in% c("mean_temp", "cwb", "cwb_range", "temp_range"))] <- "climate"

print(aggregate(perc~species+feature, dataset,sum))
print(aggregate(perc~species+feature3, dataset,sum))

## Figure drawing ####
custom_palette <- c(
  "a" = "#bdbdbd",  # sox_range 
  "b" = "#636363",  # sox      
  "c" = "#9ecae1",  # cwb_range 
  "d" = "#08519c",  # cwb      
  "e" = "#fcae91",  # temp_range 
  "f" = "#cb181d"   # mean_temp 
)

g<-ggplot(data=dataset, aes(x=x, y=perc, group=feature2, fill=feature2))
g<-g+geom_bar(stat="identity", alpha = 0.8)
g<-g+scale_fill_manual(values=custom_palette,breaks=names(custom_palette))
g<-g+geom_text(data = deviances_expl, aes(x = x, y = 1.05, label = formatC(round(V1, 2), format = "f", digits = 2)),  inherit.aes=F)
g<-g+scale_y_continuous("Relative importance (%)",limits = c(0, 1.15), breaks = seq(0, 1, 0.2),labels=seq(0, 100, 20))
g<-g+scale_x_continuous("", breaks = c(1:5),labels = c("QUSP", "FASY", "PISY", "PCAB", "ABAL"))
g<-g+coord_flip()
g<-my.theme(g)
figure<-g


