###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

library(colorBlindness)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(ppcor)

##For this to work, the code in DAR_modelling must first be run to generate
#the preds object
data <- preds

type <- data$Type_coarse

df.isar <- data[,c("ISAR_z", "Gamma", "NI", "MeanDist", "Iso", "Temperature",  "Elev", "ISAR_c", "AreaScale", "ArchArea")]
df.ifdar <- data[,c("FDAR_pd_z", "FD_Gamma", "NI", "MeanDist", "Iso", "Temperature", "Elev", "FDAR_pd_c", "AreaScale", "ArchArea")]
df.ipdar <- data[,c("PDAR_pd_z", "PD_Gamma", "NI", "MeanDist", "Iso", "Temperature",  "Elev", "PDAR_pd_c", "AreaScale", "ArchArea")]
df.ifdar.es <- data[,c("FDAR_ES_z", "FD_Gamma", "NI", "MeanDist", "Iso", "Temperature",  "Elev", "FDAR_pd_c","AreaScale", "ArchArea")]
df.ipdar.es <- data[,c("PDAR_ES_z", "PD_Gamma", "NI", "MeanDist", "Iso", "Temperature", "Elev", "PDAR_pd_c","AreaScale", "ArchArea")]


graph_heatmap <- function(data) {
  
  type <- data$Type_coarse
  df.isar <- data[,c("ISAR_z", "Gamma", "NI", "MeanDist", "Iso", "Temperature",  "Elev", "ISAR_c", "AreaScale", "ArchArea")]
  df.ifdar <- data[,c("FDAR_pd_z", "FD_Gamma", "NI", "MeanDist", "Iso", "Temperature", "Elev", "FDAR_pd_c", "AreaScale", "ArchArea")]
  df.ipdar <- data[,c("PDAR_pd_z", "PD_Gamma", "NI", "MeanDist", "Iso", "Temperature",  "Elev", "PDAR_pd_c", "AreaScale", "ArchArea")]
  df.ifdar.es <- data[,c("FDAR_ES_z", "FD_Gamma", "NI", "MeanDist", "Iso", "Temperature",  "Elev", "FDAR_pd_c","AreaScale", "ArchArea")]
  df.ipdar.es <- data[,c("PDAR_ES_z", "PD_Gamma", "NI", "MeanDist", "Iso", "Temperature", "Elev", "PDAR_pd_c","AreaScale", "ArchArea")]
  
  df.isar[,c("Gamma", "NI", "MeanDist", "Iso", "Elev", "ISAR_c", "AreaScale", "ArchArea")] = 
    log(df.isar[,c("Gamma", "NI", "MeanDist", "Iso", "Elev", "ISAR_c", "AreaScale", "ArchArea")])
  
  df.ifdar[,c("FD_Gamma", "NI", "MeanDist", "Iso", "Elev", "FDAR_pd_c", "AreaScale", "ArchArea")] = 
    log(df.ifdar[,c("FD_Gamma", "NI", "MeanDist", "Iso", "Elev", "FDAR_pd_c", "AreaScale", "ArchArea")])
  
  df.ipdar[,c("PD_Gamma", "NI", "MeanDist", "Iso", "Elev", "PDAR_pd_c", "AreaScale", "ArchArea")] = 
    log(df.ipdar[,c("PD_Gamma", "NI", "MeanDist", "Iso", "Elev", "PDAR_pd_c", "AreaScale", "ArchArea")])
  
  df.ifdar.es[,c("FD_Gamma", "NI", "MeanDist", "Iso", "Elev", "FDAR_pd_c", "AreaScale", "ArchArea")] = 
    log(df.ifdar.es[,c("FD_Gamma", "NI", "MeanDist", "Iso", "Elev", "FDAR_pd_c", "AreaScale", "ArchArea")])
  
  df.ipdar.es[,c("PD_Gamma", "NI", "MeanDist", "Iso", "Elev", "PDAR_pd_c", "AreaScale", "ArchArea")] = 
    log(df.ipdar.es[,c("PD_Gamma", "NI", "MeanDist", "Iso", "Elev", "PDAR_pd_c", "AreaScale", "ArchArea")])
  
  
  
  fct_cor <- function(x){
    csp <- list()
    for (i in 1:(ncol(x)-1)){
      csp[[i]] <- cor.test(x[,1], x[,2:ncol(x)][,i], method = "pearson") 
    }
    
    estimate <- NA
    pv <- NA
    #df <- NA
    for (j in 1:length(csp)){
      estimate[j] <- csp[[j]]$estimate
      pv[j] <- csp[[j]]$p.value
      #df[j] <- csp[[j]]$parameter
    }
    
    df.sar <- data.frame(estimate, pv=round(pv, 3))
    rownames(df.sar) <- colnames(x[,2:ncol(x)])
    
    if (nrow(x) > 30) {df.sar[c("MeanDist", "Iso"),] =  NA}
    
    
    return(df.sar)
  }
  fct_cor.es <- function(x){
    csp <- list()
    for (i in 1:(ncol(x)-1)){
      csp[[i]] <- cor.test(x[,1], x[,2:ncol(x)][,i], method = "pearson") 
    }
    
    estimate <- NA
    pv <- NA
    #df <- NA
    for (j in 1:length(csp)){
      estimate[j] <- csp[[j]]$estimate
      pv[j] <- csp[[j]]$p.value
      #df[j] <- csp[[j]]$parameter
    }
    
    df.sar <- data.frame(estimate, pv=round(pv, 3))
    rownames(df.sar) <- colnames(x[,2:ncol(x)])
    
    df.sar[grep("_c", rownames(df.sar)),] =  NA
    if (nrow(x) > 30) {df.sar[c("MeanDist", "Iso"),] =  NA}
    
    
    return(df.sar)
  }
  fct_cor_partial <- function(x, z){
    csp <- list()
    for (i in 1:(ncol(x)-1)){
      y <- x[,2:ncol(x)][,i]
      y <- ifelse(is.na(y)==TRUE, 0, y)
      csp[[i]] <- pcor.test(x[,1], y, z, method="pearson") 
    }
    
    estimate <- NA
    pv <- NA
    df <- NA
    for (j in 1:length(csp)){
      estimate[j] <- csp[[j]]$estimate
      pv[j] <- csp[[j]]$p.value
      #df[j] <- csp[[j]]$parameter
    }
    
    df.sar <- data.frame(estimate, pv=round(pv, 3))
    rownames(df.sar) <- colnames(x[,2:ncol(x)])
    if (nrow(x) > 30) {df.sar[c("MeanDist", "Iso"),] =  NA}
    
    return(df.sar)
  }
  
  df.graph.estimate <- rbind(fct_cor(df.isar), 
                             fct_cor(df.ifdar),
                             fct_cor(df.ipdar),
                             fct_cor_partial(df.ifdar, data$ISAR_z),
                             fct_cor_partial(df.ipdar, data$ISAR_z),
                             fct_cor.es(df.ifdar.es),
                             fct_cor.es(df.ipdar.es)
  )%>%as.data.frame
  #df.graph.estimate$variables <- rep(c("Gamma", "# Islands", "MeanDist", "Iso",
  #                                     "Temperature", "Precipitation", "Elevation", 
  #                                     "c", "AreaScale", "ArchArea"), 7)
  
  df.graph.estimate$variables <- rep(c("GA", "NI", "MD", "IS",
                                       "TP",  "EL", 
                                       "C", "AS", "AA"), 7)
  
  df.graph.estimate$Types <- c(rep("ISAR", 9),rep("IFDAR", 9), rep("IPDAR", 9),
                               rep("Partial IFDAR", 9), rep("Partial IPDAR", 9),
                               rep("ES IFDAR", 9), rep("ES IPDAR", 9))
  
  
  df.graph.estimate$cat <- c(rep("Raw", 45),rep("ES", 18))
  
  df.graph.estimate$Types <- factor(df.graph.estimate$Types, 
                                    levels = c("ISAR", "IFDAR", "IPDAR", 
                                               "Partial IFDAR", "Partial IPDAR",
                                               "ES IFDAR", "ES IPDAR"))
  df.graph.estimate$cat <- factor(df.graph.estimate$cat, 
                                  levels = c("Raw", "ES"))
  
  df.graph.estimate$variables <- factor(df.graph.estimate$variables, 
                                        levels = rev(c("GA", "C", "AA", "AS", "NI", "MD", "IS", "EL", "TP")))
  
  
  #Elev2        -0.407637655 0.003        EL         IPDAR Raw
  
  #df.graph.estimate$pv <- p.adjust(df.graph.estimate$pv, "fdr")
  
  library(reshape2)
  #heatmap_com <- data.frame(design = rownames(path_select), path_select)
  df.graph.estimate$estimate <- round(df.graph.estimate$estimate, 2)
  All.nba.m <- melt(df.graph.estimate)
  All.nba.m.est <- All.nba.m[All.nba.m$variable == "estimate",]
  All.nba.m.pv <- All.nba.m[All.nba.m$variable == "pv",]
  rownames(All.nba.m.pv) <- rownames(All.nba.m.est)
  All.nba.m.lab <- All.nba.m.est[All.nba.m.pv$value<0.05,]
  All.nba.m.pv <- All.nba.m.pv[All.nba.m.pv$value<0.05,]
  All.nba.m.lab <- na.omit(All.nba.m.lab)
  All.nba.m.lab <- na.omit(All.nba.m.lab)
  #hcl.colors(30, "RdYlGn")
  g_All <- ggplot(All.nba.m.est, aes(Types, variables)) + 
    facet_grid(~ cat, scales='free_x', space="free_x") +
    geom_tile(aes(fill = value),colour = "white", 
              lwd = 0.5,
              linetype = 1) + 
    scale_fill_gradientn(colors = Blue2DarkRed18Steps,
                         limit = c(-.9,.9), 
                         name="Pearson\nCorrelation",
                         na.value="grey90")  +  theme_minimal()+
    theme(axis.text.x = element_text(angle =60, 
                                     hjust = 1, size= 8, color="black"), 
          axis.text.y = element_text(hjust = 1, size= 10, color="black"),
          strip.text = element_text(size=10, lineheight=5),
          strip.text.x = element_blank()) + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0)) + labs(x="", y="", fill="") +
    geom_point(data = All.nba.m.lab, aes(Types, variables)) +
    ggtitle("All Islands")
  g_All
  # geom_text(data = All.nba.m.lab, aes(label = value), size=2) +
  #####################################################################
  ########################## only TRUE islands ########################
  #####################################################################
  
  
  df.graph.estimate <- rbind(fct_cor(df.isar[data$Type_coarse == "Water",]), 
                             fct_cor(df.ifdar[data$Type_coarse == "Water",]),
                             fct_cor(df.ipdar[data$Type_coarse == "Water",]),
                             fct_cor_partial(df.ifdar[data$Type_coarse == "Water",], data$ISAR_z[data$Type_coarse == "Water"]),
                             fct_cor_partial(df.ipdar[data$Type_coarse == "Water",], data$ISAR_z[data$Type_coarse == "Water"]),
                             fct_cor.es(df.ifdar.es[data$Type_coarse == "Water",]),
                             fct_cor.es(df.ipdar.es[data$Type_coarse == "Water",])
  )%>%as.data.frame
  #df.graph.estimate$variables <- rep(c("Gamma", "# Islands", "MeanDist", "Iso",
  #                                     "Temperature", "Precipitation", "Elevation", 
  #                                     "c", "AreaScale", "ArchArea"), 7)
  
  df.graph.estimate$variables <- rep(c("GA", "NI", "MD", "IS",
                                       "TP", "EL", 
                                       "C", "AS", "AA"), 7)
  
  df.graph.estimate$Types <- c(rep("ISAR", 9),rep("IFDAR", 9), rep("IPDAR", 9),
                               rep("Partial IFDAR", 9), rep("Partial IPDAR", 9),
                               rep("ES IFDAR", 9), rep("ES IPDAR", 9))
  
  df.graph.estimate$cat <- c(rep("Raw", 45),rep("ES", 18))
  
  df.graph.estimate$Types <- factor(df.graph.estimate$Types, 
                                    levels = c("ISAR", "IFDAR", "IPDAR", 
                                               "Partial IFDAR", "Partial IPDAR",
                                               "ES IFDAR", "ES IPDAR"))
  df.graph.estimate$cat <- factor(df.graph.estimate$cat, 
                                  levels = c("Raw", "ES"))
  
  df.graph.estimate$variables <- factor(df.graph.estimate$variables, 
                                        levels = rev(c("GA", "C", "AA", "AS", "NI", "MD", "IS", "EL", "TP")))
  
  #df.graph.estimate$variables <- factor(df.graph.estimate$variables, levels = rev(c("Gamma", "c", "# Islands", "Elevation", "Temperature", "Precipitation")))
  
  #round(p.adjust(df.graph.estimate$pv, "BH"), 3)
  
  #df.graph.estimate$pv <- p.adjust(df.graph.estimate$pv, "BH")
  
  library(reshape2)
  #heatmap_com <- data.frame(design = rownames(path_select), path_select)
  df.graph.estimate$estimate <- round(df.graph.estimate$estimate, 2)
  True.nba.m <- melt(df.graph.estimate)
  True.nba.m.est <- True.nba.m[True.nba.m$variable == "estimate",]
  True.nba.m.pv <- True.nba.m[True.nba.m$variable == "pv",]
  rownames(True.nba.m.pv) <- rownames(True.nba.m.est)
  True.nba.m.lab <- True.nba.m.est[True.nba.m.pv$value<0.05,]
  True.nba.m.pv <- True.nba.m.pv[True.nba.m.pv$value<0.05,]
  True.nba.m.lab <- na.omit(True.nba.m.lab)
  
  g_true <- ggplot(True.nba.m.est, aes(Types, variables)) + 
    facet_grid(~ cat, scales='free_x', space="free_x") +
    geom_tile(aes(fill = value),colour = "white", 
              lwd = 0.5,
              linetype = 1) + 
    scale_fill_gradientn(colors = Blue2DarkRed18Steps,
                         limit = c(-0.9,0.9), 
                         name="Pearson\nCorrelation",
                         na.value="grey90") +  theme_minimal()+
    theme(axis.text.x = element_text(angle =60, 
                                     hjust = 1, size= 8, color="black"), 
          axis.text.y = element_text(hjust = 1, size= 10, color="black"),
          strip.text = element_text(size=10, lineheight=5),
          strip.text.x = element_blank()) + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0)) + labs(x="", y="", fill="") +
    geom_point(data = True.nba.m.lab, aes(Types, variables)) +
    ggtitle("True Islands") +
    guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "Spearman\nCorrelation"))
  g_true
  
  
  res <- list(g_All, g_true)
  
  return(res)
  
}

g1 <- graph_heatmap(data)

gPlot <- ggarrange(g1[[1]],  g1[[2]], 
                        common.legend = TRUE, legend = "top", 
                        labels = c("A", "B"), ncol=2, nrow=1) 

# jpeg("Figure6_heatmap_AllDataSets_AllSP_bcTRUE.jpeg", 
#      width = 15, height = 10, 
#      res = 300, units = "cm")
# gPlot 
# dev.off()

ggplot2::ggsave("Figure_Heat_map.pdf", gPlot, 
                width = 6, height = 4.5)




