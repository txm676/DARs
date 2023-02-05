##############################################################
#########Main Function####################################
####################################################################
###error messages
#Inf means either ISAR or PDAR sar_average failed (e.g. due to less than 2
#models able to be fit)

fit_DARs <- function(dat,
                     phy_dend = a3, phy = cons_tree, n = 0, 
                     null_model = "taxa.tabels", 
                     null_n = 99,
                     norma = "none", homoT = "none",
                     confInt = FALSE, ciN = 50, grid_start = "none", 
                     grid_n = NULL,
                     power_only = FALSE, 
                     linPow = FALSE, check_trees = FALSE,
                     prune_trees = FALSE){
  #phy = phylogeny, phy_dend = dendrogram
  #n is min number of species an island can have
  #null_model is the name of null model from ses.pd function
  #null_n = no. of null iterations
  #norma = sars R package residual normality test
  #homoT = sars R package residual homogeneity test
  #confInt = calculate confidence intervals around MMI SAR curve
  #ciN = number of fits to generate the CIs
  #grid_start and grid_n = whether to use grid_start in the SAR model
  #fitting (see help docs of sar_average in sars R package)
  #power_only - only fit the simple 2 par models
  #linPow = also fit the linear power model and return z for ISAR, IPDAR and IFDAR
  #check_trees = check global dendrogram PD and (S)ES values against pruned versions
  #prune_trees = prune the phylogeny and dendrogram (TRUE) or use global version
  
  #remove Extinct column if it exists
  if (any(colnames(dat) == "Extinct")){
    dat <- dat[,-(which(colnames(dat) == "Extinct"))]
  }
  
  ##format data and get species names
  dat <- rem_isl(dat, n)#remove any islands with less than n sp
  wa <- which(dat$species == "Area (ha)")
  if (length(wa) != 1) stop("Error A")
  ar <- dat[wa,2:ncol(dat)]#area values
  ar <- ar / 100 #turn into km2
  dat2 <- dat[1:(wa - 2),]
  sp <- as.character(dat2$species)
  dat3 <- dat2[ ,2:ncol(dat2)]
  if (ncol(dat3) < 7) warning("Less than seven islands")
  rownames(dat3) <- sp
  dat4 <- t(dat3)#sp as columns
  if (any(colSums(dat3) == 0) || any(rowSums(dat3) == 0)) stop("Error C")
  
  #check all species are in the phylogeny and dendrogram
  if (!all(sp %in% phy_dend$tip.label)) stop("sp missing from dendrogram")
  if (!all(sp %in% phy$tip.label)) stop("sp missing from phylogeny")
  
  #calculate PD and PD.SES using the global phylogeny
  cs <- rowSums(dat4)
  # if(!all(cs == nmF$ntaxa)) stop("Error C")
  ar <- as.vector(unlist(ar))
  
  if (prune_trees){
    #uses pruned tree
    phy <- keep.tip(phy, sp)#prune the tree to just the species in the dataset
    p_pd <- pd(dat4, phy)$PD
  } else if (!prune_trees) {
    #uses global tree
    #this is what used in the main analyses
    p_pd <- pd(dat4, phy)$PD
  }
  
  #round all PD/FD values as found that rounding
  #errors within ses.pd can result in sites with only 1 species having 
  #large SES scores (when it should be NaN / 0)
  p_pd <- round(p_pd, 3)

  null_list <- ses.pd(dat4, phy, null.model = null_model, runs = null_n,
                      check = FALSE)
  nmF <- null_list[[1]] #ses.pd result table
  nmF[,c("ES_P", "ES")] <- ES_null(null_list)
  
  #run a series of checks on the phylos and es values
  if (check_trees) {
   ct_pd <- check_trees_fun(dat4, phy, sp, nmF, null_n)
  }

  ##Create matrix to store all the ES / SES checks
  check_vals <- matrix(ncol = 3, nrow = 2)
  rownames(check_vals) <- c("PD", "FD")
  colnames(check_vals) <- c("ES0", "Pvals", "ESz")
  
  #check ES value = 0 for island that contains all species
  if (any(nmF$ntaxa == ncol(dat4))){
    check_vals[1,1] <- all(nmF$ES[which(nmF$ntaxa == ncol(dat4))] == 0)
  } else {
    check_vals[1,1] <- TRUE
  }

  #Check ES_p-values are similar to SES p-values
  check_vals[1,2] <- any(abs(round(nmF$pd.obs.p, 2) - 
                               round(nmF$ES_P, 2)) > 0.015)
  
  #as we include 1 species islands, the SES z value returns
  #a NaN for these; so change to zero to allow the correlations
  #to be run
  if (anyNA(nmF$pd.obs.z)){
    lbm <- which(is.na(nmF$pd.obs.z))
    nmF$pd.obs.z[lbm] = 0
  }
  
  #Check correlation between ES and SES z values
  check_vals[1,3] <- cor(nmF$pd.obs.z, nmF$ES)
  
  p_mpd <- nmF$ES
  
  if (!all(cs==nmF$ntaxa)) stop("sour times")

  #dont need the ES values for power_only
  if (power_only){
    df <-  data.frame("a" = ar, "s" = cs, "p_mpd" = NA, "p_pd" = p_pd)
  } else {
    df <-  data.frame("a" = ar, "s" = cs, "p_mpd" = p_mpd, "p_pd" = p_pd)
  }

  #calculate FD and add to main dataframe
  FDL <- func_mpd(d = dat4, dend = phy_dend,
                 null.FD = null_model, null_n = null_n,
                 check_trees = check_trees,
                 prune_trees = prune_trees)
  
  FD <- FDL[[1]]
  
  check_vals[2,] <- FDL[[2]]
  
  if (check_trees) ct_fd <- FDL[[3]]

  df$f_mpd <- FD$f_mpd
  df$f_pd <- FD$f_pd
  ##for null model analyses, return  P values also
  df$PD_Zsc_p <- nmF$ES_P
  df$FD_Zsc_p <- FD$ES_P

  #for Hawaii dataset don't fit 4 par models
  if (ncol(dat3) < 7){
    obj <- c("power",
             "powerR","epm1","epm2","p1","p2","loga","koba",
             "monod","negexpo","chapman","weibull3","asymp",
             "ratio","gompertz","logistic", "heleg", "linear")
  } else { #else (for all other datasets) fit all twenty
    obj <-  c("power",
              "powerR","epm1","epm2","p1","p2","loga","koba",
              "monod","negexpo","chapman","weibull3","asymp",
              "ratio","gompertz","weibull4","betap","logistic", "heleg", "linear")
  }

  if (power_only){

    ISAR <- sar_power(data = cbind(df$a, df$s),
                      grid_start = grid_start, grid_n = grid_n)

    PDAR_pd <- sar_power(data = cbind(df$a, df$p_pd),
                         grid_start = grid_start, grid_n = grid_n)

    FDAR_pd <- sar_power(data = cbind(df$a, df$f_pd),
                         grid_start = grid_start, grid_n = grid_n)

    Iz <- ISAR$par[2] %>% as.vector()
    Pz <- PDAR_pd$par[2] %>% as.vector()
    Fz <- FDAR_pd$par[2] %>% as.vector()
    return(c(Iz, Pz, Fz))

  }

  #fit twenty DAR models
  ISAR <- tryCatch(sar_average(obj = obj, data = cbind(df$a, df$s), 
                               crit = "AICc",
                               normaTest = norma, homoTest = homoT,
                               grid_start = grid_start, grid_n = grid_n, 
                               verb = FALSE,
                               confInt = confInt, ciN = ciN),
                   error = function(e) NA)

  PDAR_mpd <- tryCatch(DAR_ES(cbind(df$a, df$p_mpd)),
                       error = function(e) NA)

  PDAR_pd <- tryCatch(sar_average(obj = obj, data = cbind(df$a, df$p_pd), 
                                  crit = "AICc",
                                  normaTest = norma, homoTest = homoT,
                                  grid_start = grid_start, grid_n = grid_n, 
                                  verb = FALSE,
                                  confInt = confInt, ciN = ciN),
                      error = function(e) NA)

  FDAR_mpd <- tryCatch(DAR_ES(cbind(df$a, df$f_mpd)),
                       error = function(e) NA)

  FDAR_pd <- tryCatch(sar_average(obj = obj, data = cbind(df$a, df$f_pd), 
                                  crit = "AICc",
                                  normaTest = norma, homoTest = homoT,
                                  grid_start = grid_start, grid_n = grid_n, 
                                  verb = FALSE,
                                  confInt = confInt, ciN = ciN),
                      error = function(e) NA)

  #get linear power model results
  if (linPow){
    ISARlp <- lin_pow(data = cbind(df$a, df$s), logT = log10)
    ISARlpz <- ISARlp$Model$coefficients[2]

    PDARlp <- lin_pow(data = cbind(df$a, df$p_pd), logT = log10)
    PDARlpz <- PDARlp$Model$coefficients[2]

    FDARlp <- lin_pow(data = cbind(df$a, df$f_pd), logT = log10)
    FDARlpz <- FDARlp$Model$coefficients[2]
  }

  ###error checking
  #if using no residual checks, we want to return an Inf for
  #everything (as this could imply some sort of error). But if
  #using residual checks, then sometimes the IPDAR say will have
  # < 2 models pass checks but the ISAR will, but all model fits
  #will be given as NA. So, in this case, we only want to return NA
  #for the individual model fits. NB we aren't interested in the ES
  #relationships during the residual checks run.
  
  if (norma == "none" & homoT == "none"){

  if ((length(ISAR) == 1 & anyNA(ISAR)) ||
      (length(PDAR_pd) == 1 & anyNA(PDAR_pd)) ||
      (length(FDAR_pd) == 1 & anyNA(FDAR_pd)) ||
      (length(FDAR_mpd) == 1 & anyNA(FDAR_mpd)) ||
      (length(PDAR_mpd) == 1 & anyNA(PDAR_mpd))){
    res1 <- rep(Inf, 24)
    res2 <- rep(Inf, 10)
    RL <- NA
    RL2 <- NA
    #add in linPowresults if linPow TRUE
    if (linPow){
      RL3 <- list(c(res1,
                  "ISAR_linPow_z" = ISARlpz,
                  "PDAR_linPow_z" = PDARlpz,
                  "FDAR_linPow_z" = FDARlpz,
                  res2),
                RL, RL2, df)#also return df which includes the z-score P values
    } else{
      RL3 <- list(c(res1, res2), RL, RL2, df)
    }
    return(RL3)
  }} else { #else you are running residual checks
    #first check if any didn't work, as if not can skip
    if ((length(ISAR) == 1 & anyNA(ISAR)) ||
        (length(PDAR_pd) == 1 & anyNA(PDAR_pd)) ||
        (length(FDAR_pd) == 1 & anyNA(FDAR_pd))){
      #for the residual run, we don't use the results matrix at
      #all (only use model fit objects to make plots, and lin-pow
      #z-values), so for ease just return inf for all again.
      res1 <- rep(Inf, 24)
      res2 <- rep(Inf, 10)
      RL2 <- NA
      
      RL <- list("ISAR" = NA, "PDAR_pd" = NA,
                 "FDAR_pd" = NA)
      #ISAR
      if (!(length(ISAR) == 1 & anyNA(ISAR))) RL[[1]] <- ISAR
      #IPDAR
      if (!(length(PDAR_pd) == 1 & anyNA(PDAR_pd))) RL[[2]] <- PDAR_pd
      #IFDAR
      if (!(length(FDAR_pd) == 1 & anyNA(FDAR_pd))) RL[[3]] <- FDAR_pd
      
      #if all NA, just return the single NA.
      if (sum(sapply(RL, length)) == 3) RL <- NA
   
      #add in linPowresults if linPow TRUE
      if (linPow){
        RL3 <- list(c(res1,
                      "ISAR_linPow_z" = ISARlpz,
                      "PDAR_linPow_z" = PDARlpz,
                      "FDAR_linPow_z" = FDARlpz,
                      res2),
                    RL, RL2, df)#also return df which includes the z-score P values
      } else{
        RL3 <- list(c(res1, res2), RL, RL2, df)
      }
      return(RL3)
    }
  }

  ####deal with mmi results for ISAR, IPDAR and IFDAR
  RL <- list("ISAR" = ISAR, "PDAR_pd" = PDAR_pd,
             "FDAR_pd" = FDAR_pd)

  #get power model
  IP <- sar_power(data = cbind(df$a, df$s), grid_start = grid_start, grid_n = grid_n)
  PP_pd <- sar_power(data = cbind(df$a, df$p_pd), grid_start = grid_start, grid_n = grid_n)
  FP_pd <- sar_power(data = cbind(df$a, df$f_pd), grid_start = grid_start, grid_n = grid_n)

  #results vector
  IS <- summary(ISAR)
  PS_pd <- summary(PDAR_pd)
  FS_pd <- summary(FDAR_pd)

  IPC <- ("Power" %in% IS$Models)#is power in models that passed checks
  PPC_pd <- ("Power" %in% PS_pd$Models)#is power in models that passed checks
  FPC_pd <- ("Power" %in% FS_pd$Models)#is power in models that passed checks

  #if power in model set, get delta AICc value
  if (IPC){
    #check re-running the power model resulted in same AICc as original
    #model selection run
    if ((IP$AICc - 
        IS$Model_table[which(IS$Model_table$Model == "power"),"AICc"])> 0.1){
      stop("Kashmir")
    }
    powAICc <- IP$AICc
    minAICc <- min(IS$Model_table$AICc)
    powDeltaI <- round(powAICc - minAICc, 1)
  } else {
    powDeltaI <- NA
  }

  #same for PDAR_pd
  if (PPC_pd){
    if ((PP_pd$AICc - 
        PS_pd$Model_table[which(PS_pd$Model_table$Model == "power"),"AICc"]) > 0.1){
      stop("KashmirB")
    }
    powAICcP_pd <- PP_pd$AICc
    minAICcP_pd <- min(PS_pd$Model_table$AICc)
    powDeltaP_pd <- round(powAICcP_pd - minAICcP_pd, 1)
  } else {
    powDeltaP_pd <- NA
  }

  #same for FDAR_pd
  if (FPC_pd){
    if ((FP_pd$AICc - 
        FS_pd$Model_table[which(FS_pd$Model_table$Model == "power"),"AICc"])> 0.1){
      stop("KashmirC")
    }
    powAICcF_pd <- FP_pd$AICc
    minAICcF_pd <- min(FS_pd$Model_table$AICc)
    powDeltaF_pd <- round(powAICcF_pd - minAICcF_pd, 1)
  } else {
    powDeltaF_pd <- NA
  }

  IS_zc <- as.vector(round(IP$par, 2))
  PD_zc_pd <- as.vector(round(PP_pd$par, 2))
  FD_zc_pd <- as.vector(round(FP_pd$par, 2))

  res1 <- c("ISAR_best" = as.character(IS$Model_table$Model[1]),
            "ISAR_shape" = IS$Model_table$Shape[1],
            "ISAR_asymptote" = IS$Model_table$Asymptote[1],
            "ISAR_pow_pass" = IPC, "ISAR_pow_delta" = powDeltaI,
            "ISAR_z" = IS_zc[2],
            "ISAR_c" = IS_zc[1], "ISAR_R2" = round(IP$R2, 2),

            "PDAR_pd_best" = as.character(PS_pd$Model_table$Model[1]),
            "PDAR_pd_shape" = PS_pd$Model_table$Shape[1],
            "PDAR_pd_asymptote" = PS_pd$Model_table$Asymptote[1],
            "PDAR_pd_pow_pass" = PPC_pd, "PDAR_pd_pow_delta" = powDeltaP_pd,
            "PDAR_pd_z" = PD_zc_pd[2],
            "PDAR_pd_c" = PD_zc_pd[1], "PDAR_pd_R2" = round(PP_pd$R2, 2),

            "FDAR_pd_best" = as.character(FS_pd$Model_table$Model[1]),
            "FDAR_pd_shape" = FS_pd$Model_table$Shape[1],
            "FDAR_pd_asymptote" = FS_pd$Model_table$Asymptote[1],
            "FDAR_pd_pow_pass" = FPC_pd, "FDAR_pd_pow_delta" = powDeltaF_pd,
            "FDAR_pd_z" = FD_zc_pd[2],
            "FDAR_pd_c" = FD_zc_pd[1], "FDAR_pd_R2" = round(FP_pd$R2, 2))

  if (linPow){
    res1 <- c(res1, "ISAR_linPow_z" = ISARlpz, "PDAR_linPow_z" = PDARlpz,
              "FDAR_linPow_z" = FDARlpz)
  }

  ###deal with threshold results
  RL2 <- list("PDAR_mpd" = PDAR_mpd, "FDAR_mpd" = FDAR_mpd)

  res2 <- c(
    "PDAR_mpd_best" = rownames(PDAR_mpd[[2]])[1],
    "PDAR_mpd_intWeight" = PDAR_mpd[[2]][which(rownames(PDAR_mpd[[2]]) == "Intercept"), "AICc_weight"],
    "PDAR_mpd_linWeight" = PDAR_mpd[[2]][which(rownames(PDAR_mpd[[2]]) == "Linear"), "AICc_weight"],
    "PDAR_lin_slope" = PDAR_mpd[[2]][which(rownames(PDAR_mpd[[2]]) == "Linear"), "Slope"],
    "PDAR_mpd_linR2" = PDAR_mpd[[2]][which(rownames(PDAR_mpd[[2]]) == "Linear"), "R2"],

    "FDAR_mpd_best" = rownames(FDAR_mpd[[2]])[1],
    "FDAR_mpd_intWeight" = FDAR_mpd[[2]][which(rownames(FDAR_mpd[[2]]) == "Intercept"), "AICc_weight"],
    "FDAR_mpd_linWeight" = FDAR_mpd[[2]][which(rownames(FDAR_mpd[[2]]) == "Linear"), "AICc_weight"],
    "FDAR_lin_slope" = FDAR_mpd[[2]][which(rownames(FDAR_mpd[[2]]) == "Linear"), "Slope"],
    "FDAR_mpd_linR2" = FDAR_mpd[[2]][which(rownames(FDAR_mpd[[2]]) == "Linear"), "R2"])

  
  #if check_trees, return check_tree info
  if (check_trees){
    ct_both <- rbind(ct_pd, ct_fd)
    rownames(ct_both) <- c("pd", "fd")
    colnames(ct_both) <- c("EScor", "SEScor", "Pcor", "PDequal")
    RL3 <- list(c(res1, res2), RL, RL2, df, check_vals, ct_both)#also return df which includes the z-score P values
  } else{
  RL3 <- list(c(res1, res2), RL, RL2, df, check_vals)#also return df which includes the z-score P values
  }
  return(RL3)
}

#function to remove islands with less than n species
rem_isl <- function(dat, n){
  ws <- which(dat$species == "sp.r")
  ss <- dat[ws,][2:ncol(dat)]
  ss <- as.numeric(ss)
  if (any(ss < n)){
    wl <- which(ss < n) + 1 #+1 needed as cols are 1 to right due to species col
    dat2 <- dat[, -wl]
    #check this has not resulted in species (rows) having zero incidence
    wa <- which(dat2$species == "Area (ha)")
    dat3 <- dat2[1:(wa - 2),]
    if (any(rowSums(dat3[,2:ncol(dat3)]) == 0)){
      w0 <- as.vector(which(rowSums(dat3[,2:ncol(dat3)]) == 0))
      dat2 <- dat2[-w0,]
    }
    return(dat2)
  } else{
    return(dat)
  }
}

#This takes our output from using the global tree and compares
#it to results from the pruned trees and does a bunch of error
#checks
check_trees_fun <- function(dat4, phy, sp, nmF1, null_n){
  phy2 <- keep.tip(phy, sp)
  #check here set to TRUE as ses.pd being used inside check_trees fun
  null_list2 <- ses.pd(dat4, phy2, 
                       null.model = "taxa.labels", runs = null_n,
                       check = TRUE)
  nmF2 <- null_list2[[1]] #ses.pd result table
  nmF2[,c("ES_P", "ES")] <- ES_null(null_list2)
  
  isc <- vector(length = 4)

  #check ES values correlate
  isc[1] <- cor(nmF1$ES, nmF2$ES)
  #check SES values correlate; if any Nans convert to zero
  if (any(is.na(nmF1$pd.obs.z))){
    if (!identical(is.na(nmF1$pd.obs.z), is.na(nmF2$pd.obs.z))){
      stop("check_tree:1b")
    }
    salmon <- which(is.na(nmF1$pd.obs.z))
    nmF1$pd.obs.z[salmon] = 0
    nmF2$pd.obs.z[salmon] = 0
  }
  isc[2] <- cor(nmF1$pd.obs.z, nmF2$pd.obs.z)
  #check P-values correlate
  isc[3] <- cor(nmF1$ES_P, nmF2$ES_P)
  #check PD values have same difference
  dyl <- nmF1$pd.obs - nmF2$pd.obs
  #needs truncating rather than rounding after changing PD rounding
  # to 3 dp
  isc[4] <- all(diff(trunc(dyl*10^1)/10^1) == 0)
  return(isc)
}


###function to calculate to run FD null models
func_mpd <- function(d, dend, null.FD, null_n,
                     check_trees, prune_trees = FALSE){
  
    sp2 <- colnames(d)
    
    if (prune_trees){
      #uses pruned tree
      dend <- keep.tip(dend, sp2)
    } 
    
    #run the null models using global dendrogram and taxa TABELS
    null_list_fd <- ses.pd(d, dend, null.model = null.FD, runs = null_n,
                        check = FALSE)
    nmF_fd <- null_list_fd[[1]] #ses.pd result table
    nmF_fd[,c("ES_P", "ES")] <- ES_null(null_list_fd)
    
    if (check_trees) {
      ct_fd2 <- check_trees_fun(d, dend, sp2, nmF_fd, null_n)
    }
    
    ##Create matrix to store all the ES / SES checks
    check_vals_fd <- matrix(ncol = 3, nrow = 1)
    rownames(check_vals_fd) <- c("FD")
    colnames(check_vals_fd) <- c("ES0", "Pvals", "ESz")
    
    #check ES value = 0 for island that contains all species
    if (any(nmF_fd$ntaxa == ncol(d))){
      check_vals_fd[1,1] <- all(nmF_fd$ES[which(nmF_fd$ntaxa == ncol(d))] == 0)
    } else {
      check_vals_fd[1,1] <- TRUE
    }
    
    #Check ES_p-values are similar to SES p-values
    check_vals_fd[1,2] <- any(abs(round(nmF_fd$pd.obs.p, 2) - 
                                    round(nmF_fd$ES_P, 2)) > 0.015)
    
    #as we include 1 species islands, the SES z value returns
    #a NaN for these; so change to zero to allow the correlations
    #to be run
    if (anyNA(nmF_fd$pd.obs.z)){
      lbm_fd <- which(is.na(nmF_fd$pd.obs.z))
      nmF_fd$pd.obs.z[lbm_fd] = 0
    }
    
    #Check correlation between ES and SES z values
    check_vals_fd[1,3] <- cor(nmF_fd$pd.obs.z, nmF_fd$ES)
    
    f_pd <- nmF_fd$pd.obs
    fmpd <- nmF_fd$ES
    mat <- as.data.frame(matrix(c(fmpd, f_pd, nmF_fd$ES_P), ncol = 3))
    
  colnames(mat) <- c("f_mpd", "f_pd", "ES_P")
  
  if (check_trees){
    mat_list <- list(mat, check_vals_fd, ct_fd2)
  } else{
    mat_list <- list(mat, check_vals_fd, NA)
  }
  
  return(mat_list)
}

#use probit transformation approach to generate ES
ES_null <- function(comm){
  
  #get ES each island
  ES_IS <- matrix(nrow = nrow(comm[[1]]), ncol = 2)
  
  for (j in 1:nrow(comm[[1]])){
    #uses same code as in IUCN script
    dis <- comm[[2]][,j]
    if (mean(dis) != comm[[1]][j,"pd.rand.mean"]) stop("lonely data")
    obs <- comm[[1]][j,"pd.obs"]
    dis <- c(dis, obs)# don't forget to add the obs into the null values
    res <- (sum(dis < obs) + sum(dis == obs)/2) / length(dis)
    ES_IS[j,1] <- res
    ES_IS[j,2] <- VGAM::probitlink(res)
  }
  return(ES_IS)
}

#fit the linear model in semi-log space to SES / ES data,
#and compare with intercept only model. Rank summary table by
#AICc and generate AICc weight. Return edited summary table
DAR_ES <- function(dat){
  
  #fit threshold model just get model fits for plotting
  model_thr <- sar_threshold(dat, mod = "ContOne", non_th_models = TRUE,
                             logAxes = "area", logT = log10)
  
  model_lm <- model_thr[[1]][[2]]
  model_null <- model_thr[[1]][[3]]
  
  sum_s <- matrix(NA, ncol = 3, nrow = 2)
  colnames(sum_s) <- c("AICc", "Slope", "R2")
  rownames(sum_s) <- c("Linear", "Intercept")
  
  sum_s[1,1] <- AICcmodavg::AICc(model_lm)
  sum_s[1,2] <- as.vector(model_lm$coefficients[2])
  sum_s[1,3] <- summary(model_lm)$r.squared
  sum_s[2,1] <- AICcmodavg::AICc(model_null)
  
  if(round(sum_s[1,1],2) != summary(model_thr)$Model_table["Linear","AICc"]){
    stop("trente")
  }
  
  sum_s <- round(sum_s, 2)
  sum_s <- as.data.frame(sum_s)
  
  sum_s <- sum_s[order(sum_s$AICc),]#sort by AICc
  
  #get AICc weights
  delta_ICs <-  sum_s$AICc - min(sum_s$AICc)
  akaikesum <- sum(exp(-0.5*(delta_ICs)))
  sum_s$AICc_weight  <- round(exp(-0.5*delta_ICs) / akaikesum, 2)
  if (sum(sum_s$AICc_weight) != 1) stop("sunday morning")
  
  return(list(model_thr, sum_s))
}

##Function to move legend of a ggplot into an empty facet grid space. Normally
#ggplot places legend outside of the facet grid
#Sourced from an online forum: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  
  # now we just need a simple call to reposition the legend
  lemon::reposition_legend(p, 'center', panel=names)
}


####################################################
########DAR PLOTTING FUNCTION########################
#######################################################


DAR_plot <- function(ll, ll2, allC = TRUE, confInt = FALSE, dataset = "Test",
                     cx = 1.9, cxM = 2.4, p.cex = 1,
                     lc = c("white", "darkslategray", "grey")){
  
  ISAR <- ll$ISAR
  FDAR_pd <- ll$FDAR_pd
  PDAR_pd <- ll$PDAR_pd
  FDAR_mpd <- ll2$FDAR_mpd[[1]]
  PDAR_mpd <- ll2$PDAR_mpd[[1]]
  
  titA <- paste0(dataset, "ISAR")
  titB <- paste0(dataset, "IFDAR")
  titC <- paste0(dataset, "IPDAR")
  titD <- paste0(dataset, "FD.ES-area")
  titE <- paste0(dataset, "PD.ES-area")
  
  
  ##ISAR
  plot(ISAR, pLeg = FALSE, allCurves = allC, confInt = confInt, ModTitle = titA, pcol = "white",
       cex.axis = cx, cex.lab = cx, cex.main = cxM, mmSep = TRUE, lwd.Sep = 6, col.Sep = "black")
  points(ISAR$details$fits$power$data, cex = p.cex, pch = 16, col ="dodgerblue2")
  
  
  ##FDAR_pd
  plot(FDAR_pd, pLeg = FALSE, allCurves = allC, confInt = confInt, ModTitle = titB,
       ylab = "FD", pcol = "white",
       cex.axis = cx, cex.lab = cx, cex.main = cxM, mmSep = TRUE, lwd.Sep = 6, col.Sep = "black")
  points(FDAR_pd$details$fits$power$data, cex = p.cex, pch = 16, col ="dodgerblue2")
  
  
  ##PDAR_pd
  plot(PDAR_pd, pLeg = FALSE, allCurves = allC, confInt = confInt, ModTitle = titC,
       ylab = "PD", pcol = "white",
       cex.axis = cx, cex.lab = cx, cex.main = cxM, mmSep = TRUE, lwd.Sep = 6, col.Sep = "black")
  points(PDAR_pd$details$fits$power$data, cex = p.cex, pch = 16, col ="dodgerblue2")
  
  
  ##FDAR_mpd
  plot(FDAR_mpd, multPlot = FALSE, ModTitle = titD,
       ylab = "FD.ES", pcol = "white", xlab = expression("Log"[10]*"(Area)"),
       cex.axis = cx, cex.lab = cx, cex.main = cxM, lcol = lc, 
       mgp=c(3.5,1,0))
  points(FDAR_mpd[[4]], cex = p.cex, pch = 16, col ="dodgerblue2")
  
  ##PDAR_mpd
  plot(PDAR_mpd, multPlot = FALSE, ModTitle = titE,
       ylab = "PD.ES", pcol = "white", xlab = expression("Log"[10]*"(Area)"),
       cex.axis = cx, cex.lab = cx, cex.main = cxM, lcol = lc, 
       mgp=c(3.5,1,0))
  points(PDAR_mpd[[4]], cex = p.cex, pch = 16, col ="dodgerblue2")
  
}


################################################################################
##boxplot of mean AICc weights with error bars for each model: Figures 2 and S2
###############################################################################

#if a model is not in the set (not fitted) it gets NA rather than 0; so it is
#then (for each model) the mean weight for datasets in which it was actually
#fitted
#RC = whether this is (TRUE) or is not (FALSE) being run for 
#the residual assumption check results

sar_plots <- function(resM, resM2, cn1 = 1.8, ca1 = 1.5,
                      cm1 = 2.7, ct1 = 1.8, N4, RC = FALSE){
  
  ##FOR AICc WEIGHTS
  MM_ISAR <- get_weights(resM, type = "ISAR", calc = "weight", RCgw = RC)
  MM_FDAR_pd <- get_weights(resM, type = "FDAR_pd", calc = "weight", RCgw = RC)
  MM_PDAR_pd <- get_weights(resM, type = "PDAR_pd", calc = "weight", RCgw = RC)
  
  MM_FDAR_mpd <- get_weights(resM2, type = "FDAR_mpd", calc = "weight", RCgw = RC)
  MM_PDAR_mpd <- get_weights(resM2, type = "PDAR_mpd", calc = "weight", RCgw = RC)
  
  
  if (!all(sapply(list(rownames(MM_ISAR),
                  rownames(MM_FDAR_pd)), 
             FUN = identical, rownames(MM_PDAR_pd)) & 
        identical(rownames(MM_FDAR_mpd), rownames(MM_PDAR_mpd)))){
    stop("eggshell")
  }
  
  MM_RAW <- data.frame("Model" = rownames(MM_ISAR), 
                       "ISAR" = rowMeans(MM_ISAR, na.rm = TRUE),
                       "FDAR_pd" = rowMeans(MM_FDAR_pd, na.rm = TRUE),
                       "PDAR_pd" = rowMeans(MM_PDAR_pd, na.rm = TRUE))
  
  MM_ES <- data.frame("Model" = rownames(MM_FDAR_mpd),
                      "FDAR_mpd" = rowMeans(MM_FDAR_mpd, na.rm = TRUE),
                      "PDAR_mpd" = rowMeans(MM_PDAR_mpd, na.rm = TRUE))
  
  
  MM_ALL1 <- MM_RAW[order(MM_RAW$ISAR),]
  MM_ALL2 <- MM_RAW[order(MM_RAW$FDAR_pd),]
  MM_ALL3 <- MM_RAW[order(MM_RAW$PDAR_pd),]
  MM_ALL4 <- MM_ES[order(MM_ES$FDAR_mpd),]
  MM_ALL5 <- MM_ES[order(MM_ES$PDAR_mpd),]
  
  
  ##HAVE to plot the top and bottom row separately as if not it makes the two bars
  #in the second row really wide
  jpeg(file = paste0("mean_weight_row1",N4,".jpeg"),
       width = 36, height = 18, units ="cm", res = 300)
  
  par(mfrow=c(1,3), mar=c(6.2,7.2,4,0.7) + 0.1)
  
  barplot(MM_ALL1$ISAR, horiz=TRUE,
          names.arg= MM_ALL1$Model,las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkgreen")
  mtext(side = 1, text = "Mean AICc weight", line = 5, cex = ct1)
  title("ISAR", adj = 0, line = 0, cex.main = cm1)
  
  barplot(MM_ALL2$FDAR_pd, horiz=TRUE,
          names.arg= MM_ALL2$Model,las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkgreen")
  mtext(side = 1, text = "Mean AICc weight", line = 5, cex = ct1)
  title("IFDAR", adj = 0, line = 0, cex.main = cm1)
  
  barplot(MM_ALL3$PDAR_pd, horiz=TRUE,
          names.arg= MM_ALL3$Model,las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkgreen")
  mtext(side = 1, text = "Mean AICc weight", line = 5, cex = ct1)
  title("IPDAR", adj = 0, line = 0, cex.main = cm1)
  
  dev.off()
  
  
  
  jpeg(file = paste0("mean_weight_row2",N4,".jpeg"), 
       width = 36, height = 7, units ="cm", res = 300)
  
  par(mfrow=c(1,3), mar=c(6.2,7.2,4,0.7) + 0.1)
  
  barplot(MM_ALL4$FDAR_mpd, horiz=TRUE,
          names.arg= MM_ALL4$Model,las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkgreen")
  mtext(side = 1, text = "Mean AICc weight", line = 5, cex = ct1)
  title("FD.ES-area", adj = 0, line = 0, cex.main = cm1)
  
  barplot(MM_ALL5$PDAR_mpd, horiz=TRUE,
          names.arg= MM_ALL5$Model,las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkgreen")
  mtext(side = 1, text = "Mean AICc weight", line = 5, cex = ct1)
  title("PD.ES-area", adj = 0, line = 0, cex.main = cm1)
  
  dev.off()
  
  
  
  #####FOR BEST MODEL TALLY#############################
  
  MM_ISAR_B <- get_weights(resM, type = "ISAR", calc = "best", RCgw = RC)
  MM_FDAR_pd_B <- get_weights(resM, type = "FDAR_pd", calc = "best", RCgw = RC)
  MM_PDAR_pd_B <- get_weights(resM, type = "PDAR_pd", calc = "best", RCgw = RC)
  MM_FDAR_mpd_B <- get_weights(resM2, type = "FDAR_mpd", calc = "best", RCgw = RC)
  MM_PDAR_mpd_B <- get_weights(resM2, type = "PDAR_mpd", calc = "best", RCgw = RC)
  
  all(sapply(list(rownames(MM_ISAR_B),
                  rownames(MM_FDAR_pd_B)), 
             FUN = identical, rownames(MM_PDAR_pd_B)) & 
        identical(rownames(MM_FDAR_mpd_B), rownames(MM_PDAR_mpd_B)))
  
  
  
  MM_RAW_B <- data.frame("Model" = rownames(MM_ISAR_B), 
                         "ISAR" = MM_ISAR_B,
                         "FDAR_pd" = MM_FDAR_pd_B,
                         "PDAR_pd" = MM_PDAR_pd_B)
  
  MM_RAW_ES <- data.frame("Model" = rownames(MM_FDAR_mpd_B),
                          "FDAR_mpd" = MM_FDAR_mpd_B,
                          "PDAR_mpd" = MM_PDAR_mpd_B)
  
  MM_ALL1_B <- filter(MM_RAW_B, ISAR > 0)
  MM_ALL1_B <- MM_ALL1_B[order(MM_ALL1_B$ISAR),]
  MM_ALL2_B <- filter(MM_RAW_B, FDAR_pd > 0)
  MM_ALL2_B <- MM_ALL2_B[order(MM_ALL2_B$FDAR_pd),]
  MM_ALL3_B <- filter(MM_RAW_B, PDAR_pd > 0)
  MM_ALL3_B <- MM_ALL3_B[order(MM_ALL3_B$PDAR_pd),]
  
  MM_ALL4_B <- MM_RAW_ES
  MM_ALL4_B <- MM_ALL4_B[order(MM_ALL4_B$FDAR_mpd),]
  MM_ALL5_B <- MM_RAW_ES
  MM_ALL5_B <- MM_ALL5_B[order(MM_ALL5_B$PDAR_mpd),]
  
  
  jpeg(file = paste0("best_model_number_",N4,".jpeg"), 
       width = 36, height = 36, units ="cm", res = 300)
  
  par(mfrow=c(2,3), mar=c(6.2,7.2,4,0.7) + 0.1)
  
  
  barplot(MM_ALL1_B$ISAR, horiz=TRUE,
          names.arg= as.character(MM_ALL1_B$Model),las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkblue")
  mtext(side = 1, text = "No. best fits", line = 4, cex = ct1)
  title("ISAR", adj = 0, line = 0, cex.main = cm1)
  
  barplot(MM_ALL2_B$FDAR_pd, horiz=TRUE,
          names.arg= as.character(MM_ALL2_B$Model),las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkblue")
  mtext(side = 1, text = "No. best fits", line = 4, cex = ct1)
  title("IFDAR", adj = 0, line = 0, cex.main = cm1)
  
  barplot(MM_ALL3_B$PDAR_pd, horiz=TRUE,
          names.arg= as.character(MM_ALL3_B$Model),las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkblue")
  mtext(side = 1, text = "No. best fits", line = 4, cex = ct1)
  title("IPDAR", adj = 0, line = 0, cex.main = cm1)
  
  barplot(MM_ALL4_B$FDAR_mpd, horiz=TRUE,
          names.arg= as.character(MM_ALL4_B$Model),las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkblue")
  mtext(side = 1, text = "No. best fits", line = 4, cex = ct1)
  title("FD.ES-area", adj = 0, line = 0, cex.main = cm1)
  
  barplot(MM_ALL5_B$PDAR_mpd, horiz=TRUE,
          names.arg= as.character(MM_ALL5_B$Model),las=2, xlab = "",
          cex.names = cn1, cex.axis = ca1, col = "darkblue")
  mtext(side = 1, text = "No. best fits", line = 4, cex = ct1)
  title("PD.ES-area", adj = 0, line = 0, cex.main = cm1)
  
  
  dev.off()
  
}

###Used within sar_plots() to extract AICc weight for each model (calc = "weight"),
#or work out which was the best fitting model

#if a model is not in the set (not fitted) it gets NA rather than 0;
#so it is then (for each model) the mean weight for datasets in which it was actually fitted
get_weights <- function(resMM, type, calc = "weight", RCgw = FALSE){
  #type = "ISAR", "PDAR" etc
  #calc = "weight" for mean AICc weight, or "best" for no of times best fit
  #RCgw = whether (TRUE) or not (FALSE) resMM is the residual assumption results
  
  if (type %in% c("ISAR", "FDAR_pd", "PDAR_pd")){
    
    mods <- c("power",
              "powerR","epm1","epm2","p1","p2","loga","koba",
              "monod","negexpo","chapman","weibull3","asymp",
              "ratio","gompertz","weibull4","betap","logistic", "heleg", "linear")
    
    #if RCgw is TRUE, check whether the model fit is NA. Only doing
    #this for the residual assumption results as there should be no
    #NAs in the main results (as no models removed due to checks).
    #With exception of Hawaii, where the two 4 par models have NA as
    #not possible to calculate AICc for them due to small sample size.
    #Remove all NA fits from resMM
    if (RCgw){
      specials <- sapply(resMM, function(x) length(x[type][[1]]))
      w_species <- which(specials == 2) #these are the good fits
      resMM <- resMM[w_species] #only keep the good fits
    }
    
    if (calc == "weight"){
      
      MM <- matrix(ncol = length(resMM), nrow = length(mods))
      rownames(MM) <- mods
      
      for (i in 1:length(resMM)){
        #get names and weights from summary table
        dI <- summary(resMM[[i]][type][[1]])
        dIn <- as.character(dI$Model_table$Model)
        diW <- dI$Model_table$Weight
        #match up models with results table and input weight values
        for (j in 1:length(dIn)){
          wn <- which(rownames(MM) == dIn[j])
          MM[wn, i] <- diW[j]
        }
      }
    } else if (calc == "best"){
      
      MM <- matrix(0, ncol = 1, nrow = length(mods))
      rownames(MM) <- mods
      
      for (i in 1:length(resMM)){
        #get names and weights from summary table
        dI <- summary(resMM[[i]][type][[1]])
        dIn <- as.character(dI$Model_table$Model)
        best <- dIn[1]
        #match up models with results table and input weight values
        wn <- which(rownames(MM) == best)
        MM[wn, 1] <-  MM[wn, 1] + 1
      }
      if(!sum(MM) == length(resMM)) stop("waaaaa")
    }
  } else {
    
    #if residual assumption check results, remove the NA fits
    if (RCgw){
      specials <- sapply(resMM, length)
      w_species <- which(specials == 2) #these are the good fits
      resMM <- resMM[w_species] #only keep the good fits
    }
    
    mods <- c("Linear", "Intercept")
    
    if (calc == "weight"){
      
      MM <- matrix(ncol = length(resMM), nrow = length(mods))
      rownames(MM) <- mods
      
      for (i in 1:length(resMM)){
        #get names and weights from summary table
        dI <- resMM[[i]][type][[1]][[2]]
        dIM <- match(mods, rownames(dI))
        MM[dIM, i] <- dI$AICc_weight
      }
    } else if (calc == "best"){
      
      MM <- matrix(0, ncol = 1, nrow = length(mods))
      rownames(MM) <- mods
      
      for (i in 1:length(resMM)){
        #get names and weights from summary table
        dI <- resMM[[i]][type][[1]][[2]]
        if (which.max(dI$AICc_weight) != 1) stop("leveee")
        best <- rownames(dI)[1]
        #match up models with results table and input weight values
        MM[match(best, rownames(MM)), 1] <-  MM[match(best, rownames(MM)), 1] + 1
      }
      if(!sum(MM) == length(resMM)) stop("waaaaa")
    }
    
  }#eo if ISAR, IPDAR, IFDAR
  return(MM)
}

##################################################################
######Get dataset predictor variables##############
###################################################

##function to get total branch length of dendrogram and phylogeny for regression (FD and PD gamma)
PF_GAM <- function(dat4, phy_dend, phy){
  
  dat5 <- dat4[1, , drop = FALSE] #turn into 1 sample
  dat5[1,] <- 1 #give all species as present
  
   pd_gam <- pd(dat5, phy)$PD #work out PD for this sample (i.e. all species in the dataset)
   fd_gam <- pd(dat5, phy_dend)$PD
  
  return(c(pd_gam, fd_gam))
}


get_pred <- function(dat, n = 1, phy_dend = a3, phy = cons_tree, resA){
  dat <- rem_isl(dat, n)#remove any islands with less than n sp
  wa <- which(dat$species == "Area (ha)")
  if (length(wa) != 1) stop("Error A")
  ar <- dat[wa,2:ncol(dat)]#area values
  ar <- ar / 100 #turn into km2
  dat2 <- dat[1:(wa - 2),]
  sp <- as.character(dat2$species)
  dat3 <- dat2[ ,2:ncol(dat2)]
  rownames(dat3) <- sp
  dat4 <- t(dat3)#sp as columns
  if (any(colSums(dat3) == 0) || any(rowSums(dat3) == 0)) stop("Error C")
  PFD_GAM <- PF_GAM(dat4, phy_dend = phy_dend, phy)
  NI <- length(ar)
  GAM <- length(sp)
  AA <- sum(ar)
  AREAScale <- max(ar) / min(ar)
  
  return(c(NI, GAM, round(PFD_GAM, 2), round(AA, 2), round(AREAScale, 2)))
}

#######################################################################
#######Number of +- Slopes for cases where linear Model Best###########
#####################################################################

ES_area_results <- function(resA){
  
  cat("\nBest models for FD-ES area relationship:\n")
  print(table(resA$FDAR_mpd_best))
  
  cat("\nBest models for PD-ES area relationship:\n")
  print(table(resA$PDAR_mpd_best))
  
  cat("\nAll isls.: No. pos and neg relationships for FD-ES:\n")
  print(resA %>%
          mutate(FD_direc = ifelse(FDAR_lin_slope < 0 , "Neg", "Pos")) %>%
          filter(FDAR_mpd_best == "Linear") %>%
          group_by(FD_direc) %>%
          summarise(n()))
  
  cat("\nAll isls.: No. pos and neg relationships for PD-ES:\n")
  print(resA %>%
          mutate(PD_direc = ifelse(PDAR_lin_slope < 0 , "Neg", "Pos")) %>%
          filter(PDAR_mpd_best == "Linear") %>%
          group_by(PD_direc) %>%
          summarise(n()))
  
  #check true island datasets in correct place
  rain <- c("Azeria(2004)_noAliens.csv",                                     
  "Baiser et al (2017) Cape Verde Birds_current_noAliens.csv",     
  "Baiser et al (2017) Cook Islands Birds_current_noAliens.csv",   
  "Baiser et al (2017) Galapagos Birds_current.csv",               
  "Baiser et al (2017) Hawaii Birds_current_noAliens.csv",         
  "Baiser et al (2017) Lesser Antilles Birds_current_noAliens.csv",
  "Baiser et al (2017) Marianas Birds_current_noAliens.csv",       
  "Baiser et al (2017) Society Islands Birds_current_noAliens.csv",
  "Bengtson&Bloch(1983)BirdsBreeding.csv",                         
  "Borgesetal_azores birds_noAliens.csv",                          
  "Borgesetal_canary birds_noAliens.csv",                          
  "Bradstreet&McCracken(1978)_noAliens.csv",                       
  "Haila et al 1983_birds.csv",                                    
  "Kubota et al. Birds Ryukus_noAliens.csv",                       
  "Martin (2022) NZ_current_noAliens.csv",                         
  "Nuddsetal(1996)FathomFiveIslandsBirds.csv",                     
  "Nuddsetal(1996)GeorgianBayBirds_noAliens.csv",                  
  "OConnell et al (2020) Wakatobi.csv",                            
  "Power (1972) chanb.csv",                                        
  "Si et al (2017) TIL birds.csv",                                 
  "Simberloff & Martin (1991) Haida Gwaii.csv",                    
  "Simberloff & Martin (1991) Maddabrd_noAliens.csv",              
  "simiakis et al (2012) all islands.csv",                         
  "Sin et al (2022) Riau Lingaa.csv",                              
  "Sin et al (2022) West Sumatra.csv") 
  
  if (nrow(resA) == 50){
  if(!identical(rain, rownames(resA)[26:nrow(resA)])){
    stop ("coming down on a sunny day: A")
  }
  } else if (nrow(resA) == 48){
    #have to remove Cook and society as not in land-birds
    rain2 <- rain[-c(3,8)]
    if(!identical(rain2, rownames(resA)[26:nrow(resA)])){
      stop ("coming down on a sunny day: A")
    }
  } else{
    stop("coming down on a sunny day: C")
  }
  
  cat("\nTrue isls.: No. pos and neg relationships for FD-ES:\n")
  #true - MAKE SURE TO CHECK THESE NUMBERS
  print(resA[26:nrow(resA),] %>%
          mutate(FD_direc = ifelse(FDAR_lin_slope < 0 , "Neg", "Pos")) %>%
          filter(FDAR_mpd_best == "Linear") %>%
          group_by(FD_direc) %>%
          summarise(n()))
  
  cat("\nTrue isls.: No. pos and neg relationships for PD-ES:\n")
  #true - MAKE SURE TO CHECK THESE NUMBERS
  print(resA[26:nrow(resA),] %>%
          mutate(PD_direc = ifelse(PDAR_lin_slope < 0 , "Neg", "Pos")) %>%
          filter(PDAR_mpd_best == "Linear") %>%
          group_by(PD_direc) %>%
          summarise(n()))
  
  cat("\nAll isls.: Median slope for FD-ES:\n")
  print(median(as.numeric(resA$FDAR_lin_slope)))
  
  cat("\nTrue isls.: Median slope for FD-ES:\n")
  print(median(as.numeric(resA[26:nrow(resA),]$FDAR_lin_slope)))
  
  cat("\nHabitat isls.: Median slope for FD-ES:\n")
  print(median(as.numeric(resA[c(1:25),]$FDAR_lin_slope)))
  
  cat("\nAll isls.: Median slope for PD-ES:\n")
  print(median(as.numeric(resA$PDAR_lin_slope)))
  
  cat("\nTrue isls.: Median slope for PD-ES:\n")
  print(median(as.numeric(resA[c(26:nrow(resA)),]$PDAR_lin_slope)))
  
  cat("\nHabitat isls.: Median slope for PD-ES:\n")
  print(median(as.numeric(resA[c(1:25),]$PDAR_lin_slope)))
  
  #different oceanic dataset numbers if landbirds used
  if (nrow(resA) == 50){
    nn = 36
    if (!identical(rain[2:11], rownames(resA)[(27:nn)])){
      stop("Don't go away")
    }
  } else if (nrow(resA) == 48){
    nn = 34
    if (!identical(rain2[2:9], rownames(resA)[(27:nn)])){
      stop("Don't go away")
    }
  }
  
  cat("\nOceanic isls.: Median slope for FD-ES:\n")
  print(median(as.numeric(resA[c(27:nn),]$FDAR_lin_slope)))
  
  cat("\nOceanic isls.: Median slope for PD-ES:\n")
  print(median(as.numeric(resA[c(27:nn),]$PDAR_lin_slope)))
}


####################################################
######REPLACE NON-LINEAR z-values with linear power
#####versions#############################################
############################################################

linPow_replace <- function(dnb){

  dnb2 <- which(names(dnb) %in% c("ISAR_linPow_z", "PDAR_linPow_z",
                                  "FDAR_linPow_z"))
  #extract and remove linpow slope values
  dnbLP <- dnb[dnb2]
  dnb3 <- dnb[-dnb2]
  #replace non-linear power z values with linpow versions
  dnb4 <- which(names(dnb3) %in% c("ISAR_z", "PDAR_pd_z",
                                   "FDAR_pd_z"))
  #for the datasets where no models fit, the main function returns inf for all
  #except the lin pow results. But the other columns have no names, so the above
  #line finds nothing. For these just use the col numbers
  if (length(dnb4) == 0){
    dnb3[c(6, 14, 22)] <- dnbLP
  } else{
    dnb3[dnb4] <- dnbLP
  }
  return(dnb3)
}
