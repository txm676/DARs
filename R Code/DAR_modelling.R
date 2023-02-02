##########################################################################
#####################################################################
######Main Script for DAR Analyses###########################
######################################################################
#########################################################################

##install packages from Github
devtools::install_github("txm676/picante", ref = "98fdbde")

library(ape)
library(picante) #NOTE THIS HAS TO BE THE VERSION ON MY GITHUB
library(dplyr)
library(sars) #needs to >= version 1.3.5
library(VGAM)
library(AICcmodavg)
library(foreach)
library(doParallel)
library(cluster)
library(ggplot2) #needs to be a recent version
#note that we noticed on some machines it was necessary to install
#a recent version of the 'vctrs' package also.

##source necessary functions
source("./R Code/DAR_SOURCE.R")

###Source the datasets, phylogeny and dendrogram
##USE BODY-SIZE CORRECTED TRAITS OR NOT (TRUE = use BS corrected traits)
body_size_CORR <- FALSE

##All species (AllSP) or just landbirds (landBird)
sp_type <- "AllSP"

#then load in the correct data given the above
source("./R Code/DAR_DATA_SOURCE.R")

####################################################################
############RUN MAIN DAR FUNCTION#####################
####################################################################

#This fits the 20 DAR models, runs the null models, and
#fits the ES-area relationship models

#Set linPow to TRUE to return the linear power z-values

#use norma = "shapiro", homoT = "cor.fitted", to add in residual
#assumptions

#check_trees runs various checks but is not essential so can be set
#to FALSE

##Note this takes many hours to run

#set up parallel back-end
cores = 10 #set depending on the machine / cluster
cl = makeCluster(cores); on.exit(stopCluster(cl))
registerDoParallel(cl)

##Run for main datasets
#main parallel for loop
rez = foreach(i=seq(from=1, to=length(ldf_all), by=1))  %dopar% {
  library(ape)
  library(sars)
  library(picante)
  library(dplyr)

  Fits <- fit_DARs(ldf_all[[i]], phy_dend = a3, phy = cons_tree, n = 0,
                   grid_start = "exhaustive",
                   grid_n = 25000,
                   null_model = "taxa.tabels", null_n = 999,
                   power_only = FALSE,
                   linPow = FALSE,
                   check_trees = FALSE)
  Fits
}

save(rez, file = "rez_all_bscFALSE.Rdata")

##Run for Alien Species
rez_alien = foreach(i=seq(from=1, to=length(ldf_alien), by=1))  %dopar% {
  library(ape)
  library(sars)
  library(picante)
  library(dplyr)

  Fits <- fit_DARs(ldf_alien[[i]], phy_dend = a3, phy = cons_tree, n = 0,
                   grid_start = "exhaustive",
                   grid_n = 25000,
                   null_model = "taxa.tabels", null_n = 999,
                   power_only = FALSE,
                   linPow = FALSE,
                   check_trees = FALSE)
  Fits
}

save(rez_alien, file = "rez_alien_bscFALSE.Rdata")

###Run for Extinct Species#######
rez_extinct = foreach(i=seq(from=1, to=length(ldf_Ex), by=1))  %dopar% {
  library(ape)
  library(sars)
  library(picante)
  library(dplyr)
  Fits <- fit_DARs(ldf_Ex[[i]], phy_dend = a3, phy = cons_tree, n = 0,
                   grid_start = "exhaustive",
                   grid_n = 25000,
                   null_model = "taxa.tabels", null_n = 999,
                   power_only = FALSE,
                   linPow = FALSE,
                   check_trees = FALSE)
  Fits
}

save(rez_extinct, file = "rez_extinct_bscFALSE.Rdata")

###################################################################
######EXTRACTING RESULTS####################################
#########################################################
 
##Set file names#####################################################
 
##create filename for figures
N1 <- "AllDataSets" #AllDataSets, true, habitat
N2 <- "AllSP" #AllSP, landBird, ResidChecks
N3 <- "bcFALSE" #bcFALSE, bcTRUE (body-size correction)
N4 <- paste0(N1,"_",N2,"_",N3)

###############################################################

##extract out different elements in rez
resM <- lapply(rez, function(x) x[[2]])#mmi model fits
resM2 <- lapply(rez, function(x) x[[3]])#EZ model fits

resA <- matrix(ncol = 34, nrow = length(rez))#main results table
for (i in 1:length(rez)){
  if (N2 == "ResidChecks"){
    #replace non-linear z values with linPow versions
    dnb <- rez[[i]][[1]]
    resA[i, ] <- linPow_replace(dnb)
  } else {
  resA[i, ] <- rez[[i]][[1]]
  }
}

resA <- as.data.frame(resA)

#There are only 48 land bird datasets as two were removed
if (N2 == "AllSP" | N2 == "ResidChecks"){
  if (length(rez) == 50 & length(lf) == 50) rownames(resA) <- lf #only add dataset name if all datasets used
} else if (N2 == "landBird"){
  if (length(rez) == 48 & length(lf) == 48) rownames(resA) <- lf #only add dataset name if all datasets used
  
}
  
colnames(resA) <- c("ISAR_best", "ISAR_shape","ISAR_asymptote","ISAR_pow_pass", 
                    "ISAR_pow_delta", "ISAR_z", "ISAR_c", "ISAR_R2", 
                    
                    "PDAR_pd_best", "PDAR_pd_shape",
                    "PDAR_pd_asymptote", "PDAR_pd_pow_pass", "PDAR_pd_pow_delta",
                    "PDAR_pd_z", "PDAR_pd_c", "PDAR_pd_R2",
                  
                    "FDAR_pd_best", "FDAR_pd_shape", "FDAR_pd_asymptote",
                    "FDAR_pd_pow_pass", "FDAR_pd_pow_delta", "FDAR_pd_z", 
                    "FDAR_pd_c", "FDAR_pd_R2",
                    
                    "PDAR_mpd_best", "PDAR_mpd_intWeight", 
                    "PDAR_mpd_linWeight", "PDAR_lin_slope", 
                    "PDAR_mpd_linR2",
                    
                    "FDAR_mpd_best", "FDAR_mpd_intWeight", 
                    "FDAR_mpd_linWeight", "FDAR_lin_slope", 
                    "FDAR_mpd_linR2")

write.csv(resA, file = paste0("DAR_results_table_",N4,".csv"))

##################################################################
######Get dataset predictor variables##############
###################################################

preds <- foreach(i=seq(from=1, to=length(ldf_all), by=1))  %dopar% {
library(ape)
library(picante)
  
  get_pred(ldf_all[[i]], n = 1, phy_dend = a3, phy = cons_tree, 
           resA = resA)
  }

preds = matrix(unlist(preds), ncol = 6, byrow = TRUE)

colnames(preds) <- c("NI", "Gamma", "PD_Gamma", "FD_Gamma", 
                     "ArchArea", "AreaScale")
rownames(preds) <- lf

preds <- as.data.frame(preds)

#get variables from resA
if (!identical(rownames(resA), lf)) stop("dataset names do not match")

preds <- preds %>%
  mutate(ISAR_c = resA$ISAR_c,
         PDAR_pd_c = resA$PDAR_pd_c,
         FDAR_pd_c = resA$FDAR_pd_c,
         ISAR_z = resA$ISAR_z,
         PDAR_pd_z = resA$PDAR_pd_z,
         FDAR_pd_z = resA$FDAR_pd_z,
         PDAR_ES_z = resA$PDAR_lin_slope,
         PDAR_ES_lin_weight = resA$PDAR_mpd_linWeight,
         FDAR_ES_z = resA$FDAR_lin_slope,
         FDAR_ES_lin_weight = resA$FDAR_mpd_linWeight,
         ISAR_R2 = resA$ISAR_R2,
         PDAR_pd_R2 = resA$PDAR_pd_R2,
         FDAR_pd_R2 = resA$FDAR_pd_R2,
         PDAR_ES_R2 = resA$PDAR_mpd_linR2,
         FDAR_ES_R2 = resA$FDAR_mpd_linR2,
  ) 

preds <- mutate_all(preds, as.numeric)

#load in elev and clim data, and island type info
ETP <- read.csv("./Data/Predictors/world_clim_all_ETP.csv")

if (N2 == "landBird"){
 ETP <- filter(ETP, Dataset %in% rownames(preds)) 
 if (nrow(ETP) != 48 | (!identical(ETP$Dataset, rownames(preds)))){
   stop("Sopranos")
 }
}

if (!identical(rownames(preds), ETP$Dataset)) stop("dataset names do not match: B")

preds <- mutate(preds, Elev = ETP$Elev_max, 
                "Temperature" = ETP$Bio1_m, "Precipitation" = ETP$Bio12_m,
                "Type_coarse" = ETP$Type_coarse,
                "Type_V_fine" = ETP$Type_V_fine,
                "Type_fine" = ETP$Type_fine,
                "Type_fine2" = ETP$Type_fine2,
                "Iso" = ETP$Iso,
                "MeanDist" = ETP$MeanDist)

write.csv(preds, file = paste0("Predictors_all_new_", N4, ".csv"))

##########################################################
####FIGURE 4##############################################
##################################################################

vars <- preds

var3 <- dplyr::select(vars, ISAR_z, FDAR_pd_z, PDAR_pd_z)
var3 <- apply(var3, 2, as.numeric) %>% as.data.frame()
var3 <- mutate(var3, Type_coarse = vars$Type_coarse)
var3 <- var3[order(var3$ISAR_z),]
var3$Time <- seq(1, nrow(var3), 1)

jpeg(paste0("z-values",N4,".jpeg"), width = 15, height = 12, res= 300, units = "cm")

ggplot(data = var3) + 
  geom_line(aes(x = Time, y = ISAR_z, colour = "ISAR")) +
  geom_line(aes(x = Time, y = FDAR_pd_z,  colour = "IFDAR")) +
  geom_line(aes(x = Time, y = PDAR_pd_z, colour = "IPDAR")) +
  geom_point(aes(x = Time, y = ISAR_z, shape = Type_coarse)) +
  geom_point(aes(x = Time, y = FDAR_pd_z, shape = Type_coarse), colour = "red") +
  geom_point(aes(x = Time, y = PDAR_pd_z, shape = Type_coarse), colour = "blue") +
  xlab("ISAR z-value rank") + ylab("DAR z-value") + theme_bw() +
  labs(shape="Island type", colour = "Legend") +
  scale_color_manual(name="DAR type",
                     breaks=c("ISAR", "IFDAR", "IPDAR"),
                     values=c("ISAR"="black", "IFDAR"="red", 
                              "IPDAR"="blue")) 


dev.off()

######################################################
######Extract ES-score & P-values############################
#####################################################

#positive Z/ES means observed > null mean, so overdispersion, and vice versa
resP <- lapply(rez, function(x) x[[4]])

Pval_mat <- matrix(nrow = length(resP), ncol = 12)
ESval_vec <- c() #a vector to store the FD ES values
FD_Zsc <- c()
PD_Zsc <- c()
ar_Zsc <- c()

#matrix to store all FD ES and PD ES values for each island
nosiz <- vapply(resP, nrow, numeric(1)) %>% sum() #number of islands
ES_FPD <- matrix(NA, ncol = 4, nrow = nosiz)
colnames(ES_FPD) <- c("f_mpd", "FD_Zsc_p", "p_mpd", "PD_Zsc_p")
k <- 1 #counter 1
k2 <- 0 #counter 2

for (i in 1:length(resP)){
  pvv <- resP[[i]]
  ESval_vec <- c(ESval_vec, pvv$f_mpd)
  ar_Zsc <- c(ar_Zsc, pvv$a)
  FD_Zsc <- c(FD_Zsc, pvv$f_mpd)
  PD_Zsc <- c(PD_Zsc, pvv$p_mpd)
  #add indiv island values to ES_FPD
  k2 <- k2 + nrow(pvv)
  ES_FPD[k:k2,] <- pvv[,c("f_mpd", "FD_Zsc_p", "p_mpd", "PD_Zsc_p")] %>% as.matrix()
  k <- k2 + 1
  
  ##FD
  if (any(pvv$FD_Zsc_p >= 0.025 & pvv$FD_Zsc_p <= 0.975)){
  f1 <- filter(pvv, FD_Zsc_p >= 0.025 & FD_Zsc_p <= 0.975)
  Pval_mat[i, 1] <- nrow(f1) #no. random FD z
  Pval_mat[i, 2] <- round((nrow(f1) / nrow(pvv)) * 100, 0) #prop. random FD z
  } else{
    Pval_mat[i, 1] <- Pval_mat[i, 2] <- 0
  }
  if (any(pvv$FD_Zsc_p < 0.025 | pvv$FD_Zsc_p > 0.975)){
    f2 <- filter(pvv, FD_Zsc_p < 0.025 | FD_Zsc_p > 0.975)
    Pval_mat[i, 3] <- length(which(f2$f_mpd > 0)) #no. overdisp FD z
    Pval_mat[i, 4] <- round((length(which(f2$f_mpd > 0)) / nrow(pvv)) * 100, 0) #prop. overdisp FD z
    Pval_mat[i, 5] <- length(which(f2$f_mpd < 0)) #no. clustering FD z
    Pval_mat[i, 6] <- round((length(which(f2$f_mpd < 0)) / nrow(pvv)) * 100, 0) #prop. clustering FD z
  } else{
    Pval_mat[i, 3] <- 0
    Pval_mat[i, 4] <- 0
    Pval_mat[i, 5] <- 0
    Pval_mat[i, 6] <- 0
  }
  
  ##PD
  if (any(pvv$PD_Zsc_p >= 0.025 & pvv$PD_Zsc_p <= 0.975)){
    f3 <- filter(pvv, PD_Zsc_p >= 0.025 & PD_Zsc_p <= 0.975)
    Pval_mat[i, 7] <- nrow(f3) #no. random FD z
    Pval_mat[i, 8] <- round((nrow(f3) / nrow(pvv)) * 100, 0) #prop. random FD z
  } else{
    Pval_mat[i, 7] <- Pval_mat[i, 8] <- 0
  }
  if (any(pvv$PD_Zsc_p < 0.025 | pvv$PD_Zsc_p > 0.975)){
    f4 <- filter(pvv, PD_Zsc_p < 0.025 | PD_Zsc_p > 0.975)
    Pval_mat[i, 9] <- length(which(f4$p_mpd > 0)) #no. overdisp FD z
    Pval_mat[i, 10] <- round((length(which(f4$p_mpd > 0)) / nrow(pvv)) * 100, 0) #prop. overdisp FD z
    Pval_mat[i, 11] <- length(which(f4$p_mpd < 0)) #no. clustering FD z
    Pval_mat[i, 12] <- round((length(which(f4$p_mpd < 0)) / nrow(pvv)) * 100, 0) #prop. clustering FD z
  } else{
    Pval_mat[i, 9] <- 0
    Pval_mat[i, 10] <- 0
    Pval_mat[i, 11] <- 0
    Pval_mat[i, 12] <- 0
  }
}

  colnames(Pval_mat) <- c("FD_No_Rand", "FD_Prop_Rand", "FD_No_Over", "FD_Prop_Over", "FD_No_clust", "FD_Prop_clust",
                        "PD_No_Rand", "PD_Prop_Rand", "PD_No_Over", "PD_Prop_Over", "PD_No_clust", "PD_Prop_clust")
  
  

  Pval_mat <- as.data.frame(Pval_mat)
  
azzume <- apply(ES_FPD, 2, mean)
if (mean(FD_Zsc) != azzume[1] | mean(PD_Zsc) != azzume[3] | 
    length(PD_Zsc) != nrow(ES_FPD)){
  stop("don't azzume")
}  

####################################################
#### Figure 6
##############################################################
ES_FPD_DF <- as.data.frame(ES_FPD)
#for each row (island) work out if FD and PD ES values are both sig or non-sig
#or one is sig and one is not etc. Codes:
#FD first and PD second, and P = plus and M = minus, and N = non-sig, so PM
# = sig positive FD and sig negative PD, while NP = non-sig FD but sig positive PD
ooc <- vector(length = nrow(ES_FPD_DF))
for (i in 1:nrow(ES_FPD_DF)){
  #for FD
  XFDV <- ES_FPD_DF[i,"FD_Zsc_p"] 
  FDV <- dplyr::case_when(
    XFDV > 0.975 ~ "P",
    XFDV < 0.025 ~ "M",
    XFDV <= 0.975 & XFDV >= 0.025 ~ "N")
  #for PD
  XPDV <- ES_FPD_DF[i,"PD_Zsc_p"]
  PDV <- dplyr::case_when(
    XPDV > 0.975 ~ "P",
    XPDV < 0.025 ~ "M",
    XPDV <= 0.975 & XPDV >= 0.025 ~ "N")
  #paste the two characters together
  ooc[i] <- paste0(FDV, PDV)
}#eo for i

##Create ES values vs each other, coloured by the codes in ooc
cbbPalette <- c("#009E73", "#56B4E9","#0072B2", "#000000", "#E69F00",  "#F0E442", 
                "#D55E00", "#CC79A7")

ES_FPD_DF$Code <- ooc

gCD <- ggplot(ES_FPD_DF, aes(x = f_mpd, y = p_mpd)) + 
  geom_point(aes(colour = Code), alpha = 0.5) +
  xlab("FD ES") + ylab("PD ES") +
  geom_smooth(method='lm', formula= y~x, colour = "black", se = FALSE) +
  geom_abline(slope = 1, intercept = 0,linetype=2, lwd = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("#E69F00","#009E73", "#0072B2", "#999999", 
                              "#56B4E9",  "#F0E442", 
                              "#D55E00")) +
  theme(legend.position = "none")


jpeg(paste0("PD_FD_",N4,".jpeg"), width = 15, height = 12, res = 300, units = "cm")
gCD
dev.off()


##FD first and PD second, and P = plus and M = minus, and N = non-sig
table(ES_FPD_DF$Code, useNA = "always")

##Work out proportion of islands clustered, overdispersed etc

#Proportion random
sum(Pval_mat$FD_No_Rand) / length(FD_Zsc)
sum(Pval_mat$PD_No_Rand) / length(PD_Zsc)

#Proportion clustered
sum(Pval_mat$FD_No_clust) / length(FD_Zsc)
sum(Pval_mat$PD_No_clust) / length(PD_Zsc)

#Proportion overdispersed
sum(Pval_mat$FD_No_Over) / length(FD_Zsc)
sum(Pval_mat$PD_No_Over) / length(PD_Zsc)

#mean values
mean(FD_Zsc)
mean(PD_Zsc)

########################################################
#######FIGURE S7
#########################################################

#first need to create a new column which states whether a point is neutral, overdisp, or clus
DF4 <- data.frame("Z" = c(FD_Zsc, PD_Zsc),
                  "DAR" = c(rep("FD", length(FD_Zsc)), rep("PD", length(PD_Zsc))))

#checked and gives same totals as in Pval_mat
DF4$Assembly <- sapply(DF4$Z, function(x){
  if (x > 1.96){
    a <- "Overdispersed"
  } else if (x < -1.96){
    a <- "Clustered"
  } else {
    a <- "Random"
  }
  a
})

e3 <- ggplot(DF4, aes(x=DAR, y=Z)) + 
  geom_point(aes(fill=DAR, color = Assembly), shape=21, alpha = 0.9,  size = 1, fill = "white",
             position=position_jitter(width=0.2, height=0)) +
  geom_boxplot(fill=NA, outlier.size = 0.2, outlier.color = "white") + labs(y="ES Value",
                                                   x="") + theme_classic() + 
  theme(axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12),
        axis.line.x  = element_line(size = 0.2),
        axis.line.y  = element_line(size = 0.2),
        panel.background = element_rect(fill = "white",size = 0.5, colour = "black")) +
  scale_colour_manual(name = "Assembly patterns", 
                      values = c("Clustered" = "#AE9C45","Random"= "#A7B8F8",
                                 "Overdispersed"= "#052955")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggtitle("a)") +
  theme(legend.position='none')


##Barplot of ES-area slopes
DF4b <- data.frame("Slope" = c(resA$FDAR_lin_slope, resA$PDAR_lin_slope),
                   "DAR" = c(rep("FD", nrow(resA)), rep("PD", nrow(resA))),
                   "Weight" = c(resA$FDAR_mpd_linWeight, resA$PDAR_mpd_linWeight))

DF4b[,c(1,3)] <- apply(DF4b[,c(1,3)], 2, as.numeric)

e4 <- ggplot(DF4b, aes(x=DAR, y=Slope)) + 
  geom_point(aes(fill=DAR, color = Weight), shape=16, alpha = 0.9,  size = 2, fill = "white",
             position=position_jitter(width=0.2, height=0)) +
  geom_boxplot(fill=NA, outlier.size = 0.2, outlier.color = "white") + labs(y="Slope",
                                                   x="") + theme_classic() + 
  theme(axis.title=element_text(size=12), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12),
        axis.line.x  = element_line(size = 0.2),
        axis.line.y  = element_line(size = 0.2),
        panel.background = element_rect(fill = "white",size = 0.5, colour = "black")) +
  ggtitle("b)") +
  theme(legend.position='none')

##plot both together
jpeg(paste0("ES_boxplots",N4,".jpeg"), width = 20, height = 10, res = 300, units = "cm")
 gridExtra::grid.arrange(e3, e4, nrow = 1)
 dev.off()

###############################################################
##Figures 2 and S2
##############################################################

#if residCheck run, some will be NA, so need to remove these first. These NAs
#are where the ISAR, IFDAR and IPDAR all had < 2 models pass the checks (a
#single NA is then returned). These NAs get removed inside sar_plots via RC =
#TRUE.
if (N2 == "ResidChecks"){
  shoop <- sapply(resM, length)
  length(which(shoop == 3))
  resM_b <- resM[which(shoop == 3)]
  resM2_b <- resM2[which(shoop == 3)]#dont care about ES value but just to keep the same length
  sar_plots(resM = resM_b, resM2 = resM2_b, N4 = N4, RC = TRUE)
} else{
  sar_plots(resM = resM, resM2 = resM2, N4 = N4, RC = FALSE) 
}
 
 ####################################################
 ########FIGURE 3########################
 #######################################################
 
 ##To plot the top 2 rows in Fig3 set the ll line to 29, for the third row
 #change it to 48
 
 #Galapagos - positive ES [lf = 29]
 #Simiakis (Aegean) - negative ES [lf = 48]

 jpeg(file = "Galapagos_DAR_plot_bcFALSE.jpeg", width = 48, height = 34, units ="cm", res = 300)
 
 ll <- resM[[29]]; ll2 <- resM2[[29]]
# ll <- resM[[48]]; ll2 <- resM2[[48]]
 
 par(mfrow = c(2, 3),mar=c(5,7,5,2) + 0.1)
 
 DAR_plot(ll, ll2, dataset = "", cx = 2.7, cxM = 3.5, p.cex = 3.1)
 
 dev.off()
 
#######################################################################
#######Number of +- Slopes for cases where linear Model Best###########
#####################################################################

##Extract all the results to do with the ES-area relationships
ES_area_results(resA)

############################################################ 
###Figure S8: boxplots of ES-area slopes by dataset type
#################################################################
jpeg(paste0("Es_slope_islType_", N4, ".jpeg"), width = 22, height = 12, 
     units = "cm", res = 300)
par(mfrow = c(1,2))
boxplot(FDAR_ES_z ~ Type_fine, data = preds,
        xlab = "Island type",
        ylab = "FD.ES-area relationship slope")
title("a)", adj = 0)
boxplot(PDAR_ES_z ~ Type_fine, data = preds,
        xlab = "Island type",
        ylab = "PD.ES-area relationship slope")
title("b)", adj = 0)
dev.off()

table(preds$Type_fine)

#############################################
####Aliens and Extinct Plots#################
############################################

alien_res <- matrix(ncol = 34, nrow = length(rez_alien))#main results table
for (i in 1:length(rez_alien)){
  if (N2 == "ResidChecks"){
    if (length(rez_alien[[i]][[2]]) == 1) next
    #replace non-linear z values with linPow versions
    dnb_al <- rez_alien[[i]][[1]]
    alien_res[i, ] <- linPow_replace(dnb_al)
  } else {
    alien_res[i, ] <- rez_alien[[i]][[1]]
  }
}

alien_res <- as.data.frame(alien_res)
rownames(alien_res) <- lfA 
colnames(alien_res) <- c("ISAR_best", "ISAR_shape","ISAR_asymptote","ISAR_pow_pass", 
                         "ISAR_pow_delta", "ISAR_z", "ISAR_c", "ISAR_R2", 
                         
                         "PDAR_pd_best", "PDAR_pd_shape",
                         "PDAR_pd_asymptote", "PDAR_pd_pow_pass", "PDAR_pd_pow_delta",
                         "PDAR_pd_z", "PDAR_pd_c", "PDAR_pd_R2",
                         
                         "FDAR_pd_best", "FDAR_pd_shape", "FDAR_pd_asymptote",
                         "FDAR_pd_pow_pass", "FDAR_pd_pow_delta", "FDAR_pd_z", 
                         "FDAR_pd_c", "FDAR_pd_R2",
                         
                         "PDAR_mpd_best", "PDAR_mpd_intWeight", 
                         "PDAR_mpd_linWeight", "PDAR_lin_slope", 
                         "PDAR_mpd_linR2",
                         
                         "FDAR_mpd_best", "FDAR_mpd_intWeight", 
                         "FDAR_mpd_linWeight", "FDAR_lin_slope", 
                         "FDAR_mpd_linR2")

#Main datasets in alien set
if (N2 == "AllSP" | N2 == "ResidChecks"){
  noAlien_nams <- c("Baiser et al (2017) Cape Verde Birds_current_noAliens.csv",     
                  "Baiser et al (2017) Cook Islands Birds_current_noAliens.csv",   
                  "Baiser et al (2017) Hawaii Birds_current_noAliens.csv",         
                  "Baiser et al (2017) Lesser Antilles Birds_current_noAliens.csv",
                  "Baiser et al (2017) Marianas Birds_current_noAliens.csv",       
                  "Baiser et al (2017) Society Islands Birds_current_noAliens.csv",
                  "Borgesetal_azores birds_noAliens.csv",                          
                  "Borgesetal_canary birds_noAliens.csv",
                  "Kubota et al. Birds Ryukus_noAliens.csv",
                  "Martin (2022) NZ_current_noAliens.csv")
} else if (N2 == "landBird"){
  noAlien_nams <- c("Baiser et al (2017) Cape Verde Birds_current_noAliens.csv",     
                    "Baiser et al (2017) Hawaii Birds_current_noAliens.csv",         
                    "Baiser et al (2017) Lesser Antilles Birds_current_noAliens.csv",
                    "Baiser et al (2017) Marianas Birds_current_noAliens.csv",
                    "Borgesetal_azores birds_noAliens.csv",                          
                    "Borgesetal_canary birds_noAliens.csv",
                    "Kubota et al. Birds Ryukus_noAliens.csv",
                    "Martin (2022) NZ_current_noAliens.csv")
}

resA_sub <- dplyr::filter(resA, rownames(resA) %in% noAlien_nams)

#check file names are in correct order
identical(sub("no.*", "", rownames(resA_sub)), 
          sub("with.*", "", rownames(alien_res)))

##########extinct results
extinct_res <- matrix(ncol = 34, nrow = length(rez_extinct))#main results table
for (i in 1:length(rez_extinct)){
  if (N2 == "ResidChecks"){
    if (length(rez_extinct[[i]][[2]]) == 1) next
    #replace non-linear z values with linPow versions
    dnb_ex <- rez_extinct[[i]][[1]]
    extinct_res[i, ] <- linPow_replace(dnb_ex)
  } else {
    extinct_res[i, ] <- rez_extinct[[i]][[1]]
  }
}

extinct_res <- as.data.frame(extinct_res)
rownames(extinct_res) <- lfEx
colnames(extinct_res) <- c("ISAR_best", "ISAR_shape","ISAR_asymptote","ISAR_pow_pass", 
                           "ISAR_pow_delta", "ISAR_z", "ISAR_c", "ISAR_R2", 
                           
                           "PDAR_pd_best", "PDAR_pd_shape",
                           "PDAR_pd_asymptote", "PDAR_pd_pow_pass", "PDAR_pd_pow_delta",
                           "PDAR_pd_z", "PDAR_pd_c", "PDAR_pd_R2",
                           
                           "FDAR_pd_best", "FDAR_pd_shape", "FDAR_pd_asymptote",
                           "FDAR_pd_pow_pass", "FDAR_pd_pow_delta", "FDAR_pd_z", 
                           "FDAR_pd_c", "FDAR_pd_R2",
                           
                           "PDAR_mpd_best", "PDAR_mpd_intWeight", 
                           "PDAR_mpd_linWeight", "PDAR_lin_slope", 
                           "PDAR_mpd_linR2",
                           
                           "FDAR_mpd_best", "FDAR_mpd_intWeight", 
                           "FDAR_mpd_linWeight", "FDAR_lin_slope", 
                           "FDAR_mpd_linR2")

##split into historic and pre-historic
#1:8 should be the historic; and 9:18 the pre-historic
if (N2 == "AllSP"| N2 == "ResidChecks"){
rownames(extinct_res)
hist_res <- extinct_res[1:8,]
pre_res <- extinct_res[9:18,]
} else if (N2 == "landBird"){
  hist_res <- extinct_res
}

##Historic

#join alien and extinct species tables together. To do this, first need to
#create columns in each with exactly same dataset names. The cells in extinct
#corresponding to datasets only in aliens are NA
#x = aliens; y = extinct
alien_res$nam_short <- sub("current.*", "", rownames(alien_res))

alien_res$nam_short <- sub("with.*", "", alien_res$nam_short)

hist_res$nam_short <- sub("historic.*", "", rownames(hist_res))

if (!all(hist_res$nam_short %in% alien_res$nam_short)) stop("all apologies")

alien_ext_res <- full_join(alien_res, hist_res, by = "nam_short", keep = TRUE)

#check file names are in correct order
identical(sub("no.*", "", rownames(resA_sub)), sub("with.*", "", rownames(alien_res)))

if (N2 == "AllSP"| N2 == "ResidChecks"){
arch <- c("Cape Verde", "Cook Isl.", "Hawaii", "L. Antilles",
          "Marianas", "Society", "Azores", "Canaries", "Ryukyu Isl.", "NZ")
} else if (N2 == "landBird"){
  arch <- c("Cape Verde", "Hawaii", "L. Antilles",
            "Marianas", "Azores", "Canaries", "Ryukyu Isl.", "NZ")
}

#x columns in alien_ext_res = aliens; y cols = extinct
exomorph_I <- data.frame("z" = c(resA_sub$ISAR_z, alien_ext_res$ISAR_z.x, 
                                 alien_ext_res$ISAR_z.y), 
                         "Type" = c(rep("B", nrow(resA_sub)),
                                    rep("C", nrow(resA_sub)),
                                    rep("A", nrow(resA_sub))),
                         "Dataset" = rep(arch, 3), 
                         "z_type" = "ISAR")

exomorph_F <- data.frame("z" = c(resA_sub$FDAR_pd_z, alien_ext_res$FDAR_pd_z.x, 
                                 alien_ext_res$FDAR_pd_z.y), 
                         "Type" = c(rep("B", nrow(resA_sub)),
                                    rep("C", nrow(resA_sub)),
                                    rep("A", nrow(resA_sub))),
                         "Dataset" = rep(arch, 3), 
                         "z_type" = "IFDAR")


exomorph_P <- data.frame("z" = c(resA_sub$PDAR_pd_z, alien_ext_res$PDAR_pd_z.x, 
                                 alien_ext_res$PDAR_pd_z.y), 
                         "Type" = c(rep("B", nrow(resA_sub)),
                                    rep("C", nrow(resA_sub)),
                                    rep("A", nrow(resA_sub))),
                         "Dataset" = rep(arch, 3), 
                         "z_type" = "IPDAR")

exomorph_ESF <- data.frame("z" = c(resA_sub$FDAR_lin_slope, alien_ext_res$FDAR_lin_slope.x, 
                                   alien_ext_res$FDAR_lin_slope.y), 
                           "Type" = c(rep("B", nrow(resA_sub)),
                                      rep("C", nrow(resA_sub)),
                                      rep("A", nrow(resA_sub))),
                           "Dataset" = rep(arch, 3), 
                           "z_type" = "FD.ES-area")

exomorph_ESP <- data.frame("z" = c(resA_sub$PDAR_lin_slope, alien_ext_res$PDAR_lin_slope.x, 
                                   alien_ext_res$PDAR_lin_slope.y), 
                           "Type" = c(rep("B", nrow(resA_sub)),
                                      rep("C", nrow(resA_sub)),
                                      rep("A", nrow(resA_sub))),
                           "Dataset" = rep(arch, 3), 
                           "z_type" = "PD.ES-area")


exomorph2 <- rbind(exomorph_I, exomorph_F, exomorph_P, exomorph_ESF, exomorph_ESP)

exomorph2$z_type <- factor(exomorph2$z_type, levels = c("ISAR", "IFDAR", "IPDAR",
                                                        "FD.ES-area", "PD.ES-area"))

exomorph2$z <- as.numeric(exomorph2$z)


if (any(is.infinite(exomorph2$z))) stop("superstar dj")

cbPalette <- c("#000000","#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#F39B7FB2")


#change Society to Society Isl
exomorph2$Dataset[which(exomorph2$Dataset == "Society")] <-
  "Society Isl."

#removes missing values (which are the NAs for Azores and Ryukyu
#not having historic dataset version)
g_alien_ext <- ggplot(exomorph2, aes(x=Type, y=z, group = Dataset)) +
  geom_point(aes(colour=Dataset), size=4.5) +
  geom_line(size=1, alpha=0.5,aes(colour=Dataset)) +
  xlab('Dataset Type') +
  ylab('z')  + 
  theme_bw() + facet_wrap(~z_type, scales = "free_y") +
  scale_colour_manual(values=cbPalette ) +
  labs(caption = "A = Historic\nB = Current - no introduced\nC = Current - with introduced") +
  theme(plot.caption.position = "plot",
        plot.caption = element_text(vjust = 8, hjust = 0, size = 11))


jpeg(file = paste0("alien_extinct_slopes_",N4,".jpeg"), width = 18, height = 16, units ="cm", res = 300)
shift_legend2(g_alien_ext)
dev.off()


##paired Wilcoxon test (with continuity correction);
#if sample size < 50, an exact P-value is calculated, but
#if there are ties or zeros (see below) exact p-values cannot be
#calculated and thus (I think) a normal approximation is used and thus
#the continuity correction is then used.

#warning for zeros means that some of the paired values are the same
#warning ties means that some of the differences between paired values are the same;
#paired wilcox works by ordering by ranked differences between pairs

wilcox_res <- matrix(ncol = 14, nrow = 5)
colnames(wilcox_res) <- c("Median_Extinct(A)", "Median_No_Alien(B)",
                          "Median_Alien(C)", "Median_difference(AB)",
                          "Median_difference(BC)","Median_difference(AC)",
                          "Numb_C>B(10)", "Numb_A>B(8)", "V(AB)", "P(AB)",
                          "V(BC)", "P(BC)","V(AC)", "P(AC)")
rownames(wilcox_res) <- unique(exomorph2$z_type)

for (i in 1:length(unique(exomorph2$z_type))){
  
  exo3 <- filter(exomorph2, z_type == unique(exomorph2$z_type)[i])
  ext3 <- filter(exo3, Type == "A")
  noalz3 <- filter(exo3, Type == "B")
  alz3 <- filter(exo3, Type == "C")
  
  #create versions with only the extinct datasets included
  wext3NA <- which(!is.na(ext3$z))
  ext3_2 <- ext3[wext3NA,]
  noalz3_2 <- noalz3[wext3NA,]
  alz3_2 <- alz3[wext3NA,]
  if (!identical(ext3_2$Dataset,alz3_2$Dataset)) stop("Polly")
  
  W_AB <- wilcox.test(ext3_2$z, noalz3_2$z, alternative = "two.sided", paired = TRUE)
  W_BC <- wilcox.test(noalz3$z, alz3$z, alternative = "two.sided", paired = TRUE)
  W_AC <- wilcox.test(ext3_2$z, alz3_2$z, alternative = "two.sided", paired = TRUE)
  
  #Note this mixes between medians of 10 datasets when A is not involved,
  #and medians of 8 datasets in anything involving A
  wilcox_res[i,] <- c(median(ext3_2$z), median(noalz3$z), median(alz3$z), 
                      median(ext3_2$z - noalz3_2$z),
                      median(noalz3$z - alz3$z),
                      median(ext3_2$z - alz3_2$z),
                      length(which(alz3$z > noalz3$z)),
                      length(which(ext3_2$z > noalz3_2$z)),
                      as.vector(W_AB$statistic),
                      as.vector(W_AB$p.value),
                      as.vector(W_BC$statistic),
                      as.vector(W_BC$p.value),
                      as.vector(W_AC$statistic),
                      as.vector(W_AC$p.value))
  
}#eo i

wilcox_res %>% round(2)

####R2 values
R2_ISAR <- cbind(alien_ext_res$ISAR_R2.y, 
                 resA_sub$ISAR_R2, alien_ext_res$ISAR_R2.x) %>%
  apply(., 2, as.numeric)
colnames(R2_ISAR) <- c("Extinct", "Cur.No.Aliens", "Cur.Aliens")

R2_IFDAR <- cbind(alien_ext_res$FDAR_pd_R2.y, resA_sub$FDAR_pd_R2, 
                  alien_ext_res$FDAR_pd_R2.x) %>%
  apply(., 2, as.numeric)
colnames(R2_IFDAR) <- c("Extinct", "Cur.No.Aliens", "Cur.Aliens")

R2_IPDAR <- cbind(alien_ext_res$PDAR_pd_R2.y, resA_sub$PDAR_pd_R2, 
                  alien_ext_res$PDAR_pd_R2.x) %>%
  apply(., 2, as.numeric)
colnames(R2_IPDAR) <- c("Extinct", "Cur.No.Aliens", "Cur.Aliens")

R2_FDARES <- cbind(alien_ext_res$FDAR_mpd_linR2.y, resA_sub$FDAR_mpd_linR2, 
                   alien_ext_res$FDAR_mpd_linR2.x) %>%
  apply(., 2, as.numeric)
colnames(R2_FDARES ) <- c("Extinct", "Cur.No.Aliens", "Cur.Aliens")

R2_PDARES <- cbind(alien_ext_res$PDAR_mpd_linR2.y, resA_sub$PDAR_mpd_linR2,
                   alien_ext_res$PDAR_mpd_linR2.x) %>%
  apply(., 2, as.numeric)
colnames(R2_PDARES) <- c("Extinct", "Cur.No.Aliens", "Cur.Aliens")

R2_all <- cbind(R2_ISAR, R2_IFDAR, R2_IPDAR, R2_FDARES, R2_PDARES)
colnames(R2_all) <- sapply(c("R2_ISAR", "R2_IFDAR", "R2_IPDAR", "R2_FDARES", "R2_PDARES"),
                           function(x) paste0(x, "_",c("Extinct", "Cur.No.Aliens", "Cur.Aliens")))

apply(R2_all, 2, function(x) round(mean(x, na.rm = TRUE), 2))

#averages across all three DARs
apply(rbind(R2_ISAR, R2_IFDAR, R2_IPDAR), 2, function(x) round(mean(x, na.rm = TRUE), 2)) 

##############################################
##Pre-historic version
##############################################################################
neomorph_I <- data.frame("z" = c(pre_res$ISAR_z), 
                         "Type" = c("Modern", "Prehistoric","Modern", "Prehistoric",
                                    "Modern", "Prehistoric", "Modern", "Prehistoric", 
                                    "Modern", "Prehistoric"),
                         "Dataset" = c("Cook", "Cook", "Haw.","Haw.", "Mar.","Mar.",
                                       "Can.", "Can.", "NZ", "NZ"), 
                         "z_type" = "ISAR")

neomorph_F <- data.frame("z" = c(pre_res$FDAR_pd_z), 
                         "Type" = c("Modern", "Prehistoric","Modern", "Prehistoric",
                                    "Modern", "Prehistoric","Modern", "Prehistoric", 
                                    "Modern", "Prehistoric"),
                         "Dataset" = c("Cook", "Cook", "Haw.","Haw.", "Mar.","Mar.",
                                       "Can.", "Can.", "NZ", "NZ"), 
                         "z_type" = "IFDAR")


neomorph_P <- data.frame("z" = c(pre_res$PDAR_pd_z), 
                         "Type" = c("Modern", "Prehistoric","Modern", "Prehistoric",
                                    "Modern", "Prehistoric", "Modern", "Prehistoric", 
                                    "Modern", "Prehistoric"),
                         "Dataset" = c("Cook", "Cook", "Haw.","Haw.", "Mar.","Mar.",
                                       "Can.", "Can.", "NZ", "NZ"),
                         "z_type" = "IPDAR")

neomorph_ESF <- data.frame("z" = c(pre_res$FDAR_lin_slope), 
                           "Type" = c("Modern", "Prehistoric","Modern", "Prehistoric",
                                      "Modern", "Prehistoric", "Modern", "Prehistoric", 
                                      "Modern", "Prehistoric"),
                           "Dataset" = c("Cook", "Cook", "Haw.","Haw.", "Mar.","Mar.",
                                         "Can.", "Can.", "NZ", "NZ"), 
                           "z_type" = "FD.ES–area")

neomorph_ESP <- data.frame("z" = c(pre_res$PDAR_lin_slope), 
                           "Type" = c("Modern", "Prehistoric","Modern", "Prehistoric",
                                      "Modern", "Prehistoric", "Modern", "Prehistoric", 
                                      "Modern", "Prehistoric"),
                           "Dataset" = c("Cook", "Cook", "Haw.","Haw.", "Mar.","Mar.",
                                         "Can.", "Can.", "NZ", "NZ"),
                           "z_type" = "PD.ES–area")


neomorph2 <- rbind(neomorph_I, neomorph_F, neomorph_P, neomorph_ESF, neomorph_ESP)

neomorph2$z_type <- factor(neomorph2$z_type, levels = c("ISAR", "IFDAR", "IPDAR",
                                                        "FD.ES–area", "PD.ES–area"))

neomorph2$z <- as.numeric(neomorph2$z)

neomorph2$Type <- factor(neomorph2$Type, levels = c("Prehistoric", "Modern"))


if (any(is.infinite(neomorph2$z))) stop("superstar dj")

#basic plot
g_pre <- ggplot(neomorph2, aes(x=Type, y=z, group = Dataset, label = Dataset)) +
  geom_point(aes(colour=Type), size=4.5) +
  geom_line(size=1, alpha=0.5) +
  xlab('Dataset Type') +
  ylab('z') +
  scale_colour_manual(values=c("#009E73", "#D55E00"), guide="none") + 
  theme_bw() + facet_wrap(~z_type, scales = "free_y") 

#add dataset names (only on the RHS: https://www.datanovia.com/en/blog/ggplot-how-to-display-the-last-value-of-each-line-as-label/)
data_ends <- neomorph2 %>% filter(Type == "Modern")
g_pre2 <- g_pre + 
  ggrepel::geom_text_repel(
    aes(label = Dataset), data = data_ends,
    fontface ="plain", color = "black", size = 3, hjust = -1)


jpeg(file = paste0("preH_slopes_",N4,".jpeg"),
     width = 18, height = 14, units ="cm", res = 300)
g_pre2 
dev.off()


#################################################################
###############REGRESSION ANALYSES######################
############################################################

library(glmmTMB)
library(dplyr)
library(tidyr)
library(car)

##test if z-values / slopes differ between DAR type

#ISAR, IFDAR, IPAR
preds$dataset <- seq(1, nrow(preds), 1)
vars_Z <- dplyr::select(preds, ISAR_z, FDAR_pd_z, PDAR_pd_z, dataset)
z_long <- tidyr::pivot_longer(vars_Z, -dataset, names_to = "Div") %>% as.data.frame
mix_z <- glmmTMB::glmmTMB(value ~ Div + (1|dataset), data = z_long, family = beta_family(),
                          REML = TRUE)
car::Anova(mix_z)


#ES
vars_Z_ES <- dplyr::select(preds, FDAR_ES_z, PDAR_ES_z, dataset)
z_long_ES <- tidyr::pivot_longer(vars_Z_ES, -dataset, names_to = "Div") %>% 
  as.data.frame
mix_z_ES <- glmmTMB::glmmTMB(value ~ Div + (1|dataset), 
                             data = z_long_ES, family = "gaussian", REML = TRUE)
car::Anova(mix_z_ES)
