##------------------------------------------------------------------------------
## Functions ####
## Calculations
# Tests for spatial autocorrelation in model residuals using Moran's I based on a 5-nearest-neighbors spatial weights matrix.
# Returns the moran.test object summarizing the strength and significance of spatial autocorrelation in the model residuals.
calculate_sp_autocorrel_in_model<-function(input_model, input_dataset){
  
  ## Testing arguments
  # input_model=mod_abal
  # input_dataset = abal_mod_dataset_scaled
  
  # Extract residuals from the model
  residuals_mod <- residuals(input_model, type = "normalized")  # Use normalized residuals
  
  # Create spatial weight matrix (k-nearest neighbors)
  coords <- cbind(input_dataset$coord_x, input_dataset$coord_y)
  # coords <- cbind(in_data$coord_x, in_data$coord_y)
  nb <- knn2nb(knearneigh(coords, k = 5))  # 5 nearest neighbors
  lw <- nb2listw(nb, style = "W")  # Convert to listw format
  
  # Moran's I test for residuals
  moran_test <- moran.test(residuals_mod, lw)
  return(moran_test)
}

# Extracts effect sizes (fixed-effect estimates) for a predefined set of predictors and their interactions from a fitted model.
# Returns a data frame with human-readable variable labels, effect sizes, and a significance flag based on p < 0.05.
calc.effect.size<-function(input){
  
  ## Testing arguments
  # input=mod_main
  
  var.names<-c("median_range",
               "median_age",
               "mid_TRW",
               
               "mean_temp",
               "mean_cwb",
               
               "sox",
               "pH_L1",
               "C.N_FH",
               
               "median_age:mean_temp",  
               "median_age:mean_cwb", 
               
               "median_range:mean_temp", 
               "median_range:mean_cwb", 
               
               "mid_TRW:pH_L1", 
               "mid_TRW:sox" , 
               
               "mean_temp:pH_L1",  
               "mean_temp:sox"
  )
  
  var.labels<-c("Age range",
                "Age",
                "TRW",
                
                "Temp.",
                "CWB",
                
                "S imisions",
                "pH",
                "C:N",
                
                "Age : Temp."    , 
                "Age : CWB"     ,
                
                "Age range : Temp." ,
                "Age range : CWB"  ,
                
                "TRW : pH" , 
                "TRW : S imisions" , 
                
                "Temp. : pH"       , 
                "Temp : S imisions" 
  )
  
  # eff_sizes<-sqrt((summary(input)[["coefficients"]][,1]^2)/(summary(input)[["coefficients"]][,4]^2+summary(input)[["coefficients"]][,3]))
  sm<-summary(input)
  eff_sizes<-as.data.frame(sm$tTable)
  eff_sizes<-eff_sizes[which(row.names(eff_sizes) %in% var.names),]
  
  pvals.out<-rep("non-significant", times = nrow(eff_sizes))
  
  pvals.out[which(eff_sizes$`p-value`<0.05)]<-"significant"
  
  eff_sizes<-data.frame(x=c(nrow(eff_sizes):1),
                        Variable = var.labels,
                        Effect_size = eff_sizes$Value,
                        pval=pvals.out)
  
  return(eff_sizes)
}

# Rescales effect sizes to emphasize relative magnitude by multiplying by 100 and applying a signed log transform.
# Shrinks very small effects toward zero, preserves sign, and returns the input data frame with transformed Effect_size values.
rescale.effects<-function(input){
  
  ## Testing arguments
  # input=calc.effect.size(mod_pcab)
  
  output<-input
  
  # Multiply by 100 to scale
  scaled_eff_size <- input$Effect_size * 100
  
  # Log-transform the absolute values
  log_transformed <- log(abs(scaled_eff_size))
  
  # Minimal effects hawing yero values
  log_transformed[which(log_transformed<0)] <- 0
  
  # Preserve the sign
  log_transformed <- log_transformed * sign(scaled_eff_size)
  
  # Replace -Inf (log(0)) with NA for clarity
  log_transformed[is.infinite(log_transformed)] <- 0
  
  # View results
  output$Effect_size <- log_transformed
  
  return(output)
}

# Extracts year-level random effects from a mixed model and fits a GLS model with AR(1) errors to estimate their temporal trend.
# Returns a data frame with year, random effect, predicted trend, and species label for plotting or further analysis.
calc.re.year<-function(in.mod,sp){
  
  ## Testing arguments
  # in.mod=mod_pcab_main
  # sp="PCAB"
  
  re_year<-as.data.frame(ranef(in.mod))
  re_year$year<-as.numeric(as.character(rownames(re_year)))
  re_year<-re_year[,c(2,1)]
  names(re_year)<-c("year","re")
  
  trend_gls <- gls(re ~ year, 
                   data = re_year, 
                   correlation = corAR1(form = ~ year))
  
  re_year$trend<-predict(trend_gls)
  re_year$species<-sp
  
  return(re_year)
  
}

# Tests for a linear temporal trend in year-level random effects using GLS with AR(1) errors.
# Returns a formatted character string containing the p-value of the year effect for the given species.
re.trend.year<-function(in.mod,sp){
  
  ## Testing arguments
  # in.mod=mod_pcab_main
  # sp="PCAB"
  
  re_year<-as.data.frame(ranef(in.mod))
  re_year$year<-as.numeric(as.character(rownames(re_year)))
  re_year<-re_year[,c(2,1)]
  names(re_year)<-c("year","re")
  
  trend_gls <- gls(re ~ year, 
                   data = re_year, 
                   correlation = corAR1(form = ~ year))
  
  sm <- summary(trend_gls)
  
  
  return(paste0("pval ", sp, " = ", round(sm$tTable[2,4], 3)))
  
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
my.theme.noy<-function(graph,legend.pos="bottom"){
  graph<-graph+theme_classic()
  graph<-graph+theme(axis.line.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y=element_blank(),
                     axis.line.x = element_line(colour="black"),
                     axis.text.x = element_text(colour="black"),
                     axis.ticks.x = element_line(color = "black"),
                     legend.position = legend.pos)
  graph
}
plot.rsqr<-function(input, sp){
  
  # input=r.squaredGLMM(mod_abal)
  # sp="ABAL"
  
  if(sp == "MAIN"){
    graph_col <- "#777777"
  } else{
    graph_col <- cols.species[sp]
  }
  
  show<-data.frame(x=1,
                   y=c(input[1],input[2]-input[1]),
                   alpha=c("b.fixed","a.random"))
  
  g<-ggplot(show)
  g<-g+geom_bar(aes(x=x,y=y,alpha=alpha),
                fill=graph_col,
                colour=graph_col,
                stat = "identity", 
                width=0.75)
  g<-g+scale_x_continuous("")
  g<-g+scale_y_continuous("", limits=c(0,1))
  g<-g+scale_alpha_manual(breaks=c("b.fixed", "a.random"), values=c(0.6,0.1))
  g<-my.theme.noy(g,"none")
  g<-g+ coord_flip()
  g
}
plot.effect.sizes<-function(input, sp){
  
  # input=mod_main
  # sp="MAIN"
  
  show<-calc.effect.size(input)
  # show$x<-c(1:nrow(show))
  
  if(sp == "MAIN"){
    graph_col <- "#777777"
  } else{
    graph_col <- cols.species[sp]
  }
  
  group_ribbons <- data.frame(
    xmin = c(13.5, 11.5, 8.5, 4.5, 2.5, 0.5),   # spodní hranice skupin podle x (menší číslo = níž v grafu)
    xmax = c(16.5, 13.5, 11.5, 8.5, 4.5, 2.5),   # horní hranice skupin
    fill = c("#E6F0DA",  # Stanoviště
             "#D0E1F9",  # Klima
             "#FCEFD9",  # Chemismus
             "#EBDCF8",  # Stanoviště & klima
             "#F8DDE0",  # Stanoviště & chemie
             "#F4E0F0"   # Klima & chemie
    ),
    label = c("Stanoviště", 
              "Klima", 
              "Chemismus",
              "Stanoviště & klima",
              "Stanoviště & chemie",
              "Klima & chemie")
  )
  group_dividers <- data.frame(
    yintercept = c(13.5, 11.5, 8.5, 4.5, 2.5)
  )
  
  g<-ggplot(show)
  g <- g + geom_rect(data = group_ribbons,
                     aes(xmin = xmin, xmax = xmax, ymin = -0.2, ymax = 0.2),
                     inherit.aes = FALSE,
                     fill = group_ribbons$fill,
                     alpha = 0.33)
  g <- g + geom_vline(data = group_dividers,
                      aes(xintercept = yintercept),
                      linetype = "dotted",
                      color = "gray40",
                      linewidth = 0.3,
                      inherit.aes=F)
  g<-g+geom_bar(aes(x=x,y=Effect_size,alpha=pval),
                fill=graph_col,
                colour=graph_col,
                stat = "identity",
                position="dodge", 
                width=0.75)
  g<-g+geom_hline(yintercept=0,colour="#000000",linewidth=0.12)
  g<-g+scale_x_continuous("Variable",
                          limits=c(0.5,nrow(show)+0.5),
                          breaks=c(1:nrow(show)),
                          labels=rev(show$names))
  g<-g+scale_y_continuous("Effect size",
                          limits=c(-0.2,0.20),
                          breaks=seq(-1.05,1.05,0.05),
                          labels=formatC(seq(-1.05,1.05,0.05), format="f", digits= 2))
  g<-g+scale_alpha_manual(breaks=c("significant", "non-significant"), values=c(0.6,0.1))
  g<-my.theme(g,"none")
  g<-g+ coord_flip()
  g
}
plot.effect.sizes.scaled.effects<-function(input, sp){
  
  # input=mod_pcab
  # sp="PCAB"
  
  show<-calc.effect.size(input)
  show<-rescale.effects(show)
  # show$x<-c(1:nrow(show))
  
  y_labels <- c(-1, -0.5, -0.3, -0.15, -0.1, -0.05, -0.02, 0, 0.02, 0.05, 0.1, 0.15, 0.3, 0.5, 1)
  
  y_ticks <- seq(-5, 5, by = 0.5)
  
  # Back-transform to the original scale
  y_ticks <- sign(y_labels) * log(abs(y_labels*100))
  
  
  g<-ggplot(show)
  g<-g+geom_hline(yintercept=0,colour="#000000",linewidth=0.12)
  g<-g+geom_bar(aes(x=x,y=Effect_size,alpha=pval),fill=cols.species[sp],stat = "identity",position="dodge")
  g<-g+scale_x_continuous("Variable",
                          limits=c(0.5,nrow(show)+0.5),
                          breaks=c(1:nrow(show)),
                          labels=rev(show$Variable))
  g<-g+scale_y_continuous("Effect size",
                          # limits=c(-5,5),
                          limits=c(-2.8,2.8),
                          breaks=y_ticks,
                          labels=formatC(y_labels,format="f",digits=2)
  )
  g<-g+scale_fill_manual(values=cols.species, breaks=names(cols.species))
  g<-g+scale_alpha_manual(breaks=c("significant", "non-significant"), values=c(1,0.25))
  g<-my.theme(g,"none")
  g<-g+ coord_flip()
  g
}
plot.panel<-function(in.mod,sp){
  
  # in.mod=mod_pcab
  # sp="PCAB"
  
  out<-ggarrange(plot.rsqr(r.squaredGLMM(in.mod),sp),
                 plot.effect.sizes(in.mod,sp),
                 align="hv",nrow=2,ncol=1,heights = c(0.2,0.8))
}
plot.panel.scaled.eff.sizes<-function(in.mod,sp){
  
  # in.mod=mod_pcab_main
  # sp="PCAB"
  
  out<-ggarrange(plot.rsqr(r.squaredGLMM(in.mod),sp),
                 plot.effect.sizes.scaled.effects(in.mod,sp),
                 align="hv",nrow=2,ncol=1,heights = c(0.15,0.85))
  return(out)
}

##------------------------------------------------------------------------------
## Dataset preparation ####
print(re.trend.year(mod_main,"MAIN"))
print(re.trend.year(mod_abal,"ABAL"))
print(re.trend.year(mod_pcab,"PCAB"))
print(re.trend.year(mod_pisy,"PISY"))
print(re.trend.year(mod_fasy,"FASY"))
print(re.trend.year(mod_qusp,"QUSP"))

dataset.year<-rbind(calc.re.year(mod_main,"MAIN"),
                    calc.re.year(mod_abal,"ABAL"),
                    calc.re.year(mod_pcab,"PCAB"),
                    calc.re.year(mod_pisy,"PISY"),
                    calc.re.year(mod_fasy,"FASY"),
                    calc.re.year(mod_qusp,"QUSP"))

##------------------------------------------------------------------------------
## Model Graphs
models<-ggarrange(plot.panel(mod_main,"MAIN"),
                  plot.panel(mod_abal,"ABAL"),
                  plot.panel(mod_pcab,"PCAB"),
                  plot.panel(mod_pisy,"PISY"),
                  plot.panel(mod_fasy,"FASY"),
                  plot.panel(mod_qusp,"QUSP"),
                  
                  nrow=1,ncol=6,labels=LETTERS[1:6])


##------------------------------------------------------------------------------
## Random effects graph ####
## random intercept - year
g<-ggplot(dataset.year)

g<-g+annotate("rect",xmin=1971,xmax=1992,ymin=-0.3,ymax=0.3,fill="#BBBBBB",alpha=0.25)
g<-g+geom_rect_pattern(
  xmin = 1992, xmax = 1998, ymin=-0.3,ymax=0.3,
  inherit.aes = FALSE,
  pattern       = "stripe",   # typ vzorku
  pattern_colour= "#BBBBBB",
  pattern_fill  = "#BBBBBB",
  pattern_angle = -45,
  pattern_density = .1,
  alpha=0.25,
  fill = NA,                 # žádná plná barva
  colour = NA
)

g<-g+geom_vline(xintercept=c(1992, 2003),linetype="dotted",size=1.5,colour="#D73027")
g<-g+geom_vline(xintercept=c(1995, 2010),linetype="dotted",size=1.5,colour="#26466D")

g<-g+geom_point(aes(x=year,y=re,colour=species), size=2.25, alpha=0.6)
g<-g+geom_line(aes(x=year,y=trend,colour=species), linewidth=1.25)

g<-g+scale_colour_manual(values=cols.species, breaks=names(cols.species))
g<-g+scale_x_continuous("",limits=c(1960,2020),breaks=seq(0,3000,10))
g<-g+scale_y_continuous("",limits=c(-0.3,0.3),breaks=seq(-1,1,0.1),labels=formatC(seq(-1,1,0.1),format="f",digits=1))

g<-my.theme(g,"bottom")
panel_a<-g


##------------------------------------------------------------------------------
## Final output ####
figure<-ggarrange(models, 
                  ggarrange(panel_a, NULL, nrow=1, ncol=2, labels=LETTERS[6], align="hv", widths=c(0.60,0.40), legend="none"),
                  nrow=2, heights=c(0.60,0.40))



 