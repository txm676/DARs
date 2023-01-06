#############################################
##########Load in data#####################
######################################################

###load consensus phylogeny
cons_tree <- readRDS("Data\\CONS_TREE.rds")
  
##load in and format traits for extant species
#TRAIT_DATA: 1st element = extant species trait data; 2nd element = extinct sp
TRAIT_DATA <- readRDS("Data\\TRAIT_DATA.rds")
all_traits <- TRAIT_DATA[[1]]
  
traits <- dplyr::select(all_traits, "Beak.Length.culmen", "Beak.Length.nares", "Beak.Width", 
                        "Beak.Depth", "Tarsus.Length",
                        "Wing.Length", "Secondary1","Tail.Length", "Mass")

if (all(!rownames(traits) %in% cons_tree$tip.label)) stop("names do not match")

##KIWIs#######
#the kiwi species are extreme outliers in regard to wing length, secondary and
#tail length. What we will do is replace
#them with the mean non-kiwi observed trait values.
w3 <- which(rownames(traits) %in% c("Apteryx_australis", "Apteryx_haastii",
                                    "Apteryx_mantelli","Apteryx_owenii"))
#create version without kiwis
traits_noKiwi <- traits[-w3,]

#work out smallest trait values for remaining species (or mean)
mean_noKiwi <- apply(traits_noKiwi[,c("Wing.Length",
                                      "Secondary1","Tail.Length")], 2, mean)

#switch kiwi values for these mean values
cn_noKiwi <- c("Wing.Length", "Secondary1","Tail.Length")
for (i in 1:3){
  traits[w3,cn_noKiwi[i]] <- mean_noKiwi[i]
}

if (!identical(apply(traits[,c("Wing.Length", "Secondary1","Tail.Length")], 2, mean),
              mean_noKiwi)) stop("Township rebellion")

####Load in traits of extinct species
extinct_traits <- TRAIT_DATA[[2]]
#NB. four species are labelled in datasets as extinct but not in the extinct dataset, as they
#are in AVONET / BIRDTREE
# #[1] "Hemignathus_lucidus"      "Paroreomyza_maculata"    
# #[3] "Corvus_hawaiiensis"       "Todiramphus_cinnamominus"

#NB five species are in the dendrogram but not in the datasets (3 from
#Canaries and two from NZ). These are all marine extinct species from
#the prehistoric period, so not in the prehistoric datasets. 

traits_ex <- dplyr::select(extinct_traits, "Beak.Length.culmen", "Beak.Length.nares", 
                           "Beak.Width", 
                           "Beak.Depth", "Tarsus.Length",
                           "Wing.Length", "Kipp.s.Distance", "Tail.Length", "Mass")

#for extinct calculate Secondary as Wing - Kipps
traits_ex <- mutate(traits_ex, Secondary1 = Wing.Length - Kipp.s.Distance)
traits_ex <- dplyr::select(traits_ex, "Beak.Length.culmen", "Beak.Length.nares", "Beak.Width", 
                           "Beak.Depth", "Tarsus.Length",
                           "Wing.Length", "Secondary1", "Tail.Length", "Mass")

#round all traits to 1dp to match
traits_ex <- apply(traits_ex, 2, round, 1) %>% as.data.frame()

if(!identical(colnames(traits), colnames(traits_ex))) stop("Faint hope")

#merge the two trait datasets together
traits <- bind_rows(traits, traits_ex)

if (all(!rownames(traits) %in% cons_tree$tip.label)) stop("names do not match: B")

#format traits and build dendrogram
traits1 <- apply(traits, 2, log) %>% as.data.frame() #log-transform following Pigot et al

###########################################################
##USE BODY-SIZE CORRECTED TRAITS OR NOT
#body_size_CORR <- FALSE
#body_size_CORR <- TRUE #use BS corrected traits
###############################################################

if (body_size_CORR){
  
  ######body mass regression approach 
  #function to return residuals from trait-mass regression
  #resp = morpho trait to use as a response
  trait_residuals <- function(traitExtant, resp){
    modBlur <- lm(traitExtant[,resp] ~ traitExtant$Mass)
    residBlur <- residuals(modBlur)
    return(residBlur)
  }
  
  residTraits <- vapply(c("Beak.Length.culmen",	"Beak.Length.nares", 
                          "Beak.Width",	"Beak.Depth",
                          "Tarsus.Length", "Wing.Length", "Secondary1", 
                          "Tail.Length"), 
                        function(x) trait_residuals(traits1, resp = x),
                        FUN.VALUE = numeric(nrow(traits1)))
  
  residTraits <- cbind(residTraits, "Mass" = traits1$Mass) 
  
  
  traits2 <- apply(residTraits, 2, scale)
  rownames(traits2) <- rownames(traits)

} else {
  
  traits2 <- apply(traits1, 2, scale)
  rownames(traits2) <- rownames(traits)
}
eu.dist <- dist(traits2, method = "euclidean")
a2 <- hclust(eu.dist, method = "average") ##UPGMA as in Petchey & Gaston (2002)
a3 <- ape::as.phylo(a2)#convert into phylo object

# BAT::tree.quality(eu.dist,a3)

##############################################
####Load the datasets#####################
#############################################

if (sp_type == "AllSP"){
  #AllSP_DATA = 1 = lf; 2=ldf_null; 3=lfA;4=ldf_alien; 5 = lfEH; 6 = ldf_ExHs;
  #7 = lfEP; 8 = ldf_ExPh; 9 = lfEX; 10 = ldf_Ex
  AllSP_DATA <- readRDS("Data\\AllSP_DATA.rds")
  lf <- AllSP_DATA[[1]]#filenames
  ldf_all <- AllSP_DATA[[2]]#datasets
  
  ##check dataset names all correct and in correct order:
  rainR <- c("Azeria(2004)_noAliens.csv",                                     
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
  if (!all(rainR == lf[26:50])){
    stop("every time: A")
  }
  
  ##Alien
  lfA <- AllSP_DATA[[3]]#filenames
  ldf_alien <- AllSP_DATA[[4]]#alien datasets
  
  if(!identical(lfA,c("Baiser et al (2017) Cape Verde Birds_current_withAlien.csv",
                      "Baiser et al (2017) Cook Islands Birds_current_withAlien.csv",
                      "Baiser et al (2017) Hawaii Birds_current_withAlien.csv",
                      "Baiser et al (2017) Lesser Antilles Birds_current_withAlien.csv",
                      "Baiser et al (2017) Marianas Birds_current_withAlien.csv",
                      "Baiser et al (2017) Society Islands Birds_current_withAlien.csv",
                      "Borgesetal_azores birds_withAlien.csv",
                      "Borgesetal_canary birds_withAlien.csv",
                      "Kubota et al. Birds Ryukus_withAlien.csv",
                      "Martin (2022) NZ_current_withAlien.csv"))) stop("Wondergan")
  
  ##Historic
  lfEH <- AllSP_DATA[[5]]#filenames
  ldf_ExHs <- AllSP_DATA[[6]]#historic period datasets
  
  ##Prehistoric
  lfEP <- AllSP_DATA[[7]]#filenames
  ldf_ExPh <- AllSP_DATA[[8]]#prehistoric period datasets
  
  #join together
  lfEx <- AllSP_DATA[[9]]#filenames
  ldf_Ex <- AllSP_DATA[[10]]#historic & prehistoric period datasets
  
  if(!identical(c(lfEH, lfEP),c("Baiser et al (2017) Cape Verde Birds_historic.csv",
                                "Baiser et al (2017) Cook Islands Birds_historic.csv",
                                "Baiser et al (2017) Hawaii Birds_historic.csv",
                                "Baiser et al (2017) Lesser Antilles Birds_historic.csv",
                                "Baiser et al (2017) Marianas Birds_historic.csv",
                                "Baiser et al (2017) Society Islands Birds_historic.csv",
                                "Borgesetal_canary birds_historic.csv",
                                "Martin (2022) NZ_historic.csv",
                                "Baiser et al (2017) Cook Islands Birds_Modern_steadman.csv",
                                "Baiser et al (2017) Cook Islands Birds_Prehistoric_steadman.csv",
                                "Baiser et al (2017) Hawaii Birds_noMarine_Modern.csv",
                                "Baiser et al (2017) Hawaii Birds_noMarine_Prehistoric.csv",
                                "Baiser et al (2017) Marianas Birds_Modern_steadman.csv",
                                "Baiser et al (2017) Marianas Birds_Prehistoric_steadman.csv",
                                "Borgesetal_canary birds_noMarine_Modern.csv",
                                "Borgesetal_canary birds_noMarine_Prehistoric.csv",
                                "Martin (2022) NZ_noMarine_Modern.csv",
                                "Martin (2022) NZ_prehistoric_NOMarine.csv"))) stop("Eyen")
  
} else if (sp_type == "landBird"){

##Run for Terrestrial passerine##
  LANDBIRD_DATA <- readRDS("Data\\LANDBIRD_DATA.rds")
  lf <- LANDBIRD_DATA[[1]]#filenames
  ldf_all <- LANDBIRD_DATA[[2]]#datasets

  rainR <- c("Azeria(2004)_noAliens.csv",                                     
             "Baiser et al (2017) Cape Verde Birds_current_noAliens.csv",     
             "Baiser et al (2017) Galapagos Birds_current.csv",               
             "Baiser et al (2017) Hawaii Birds_current_noAliens.csv",         
             "Baiser et al (2017) Lesser Antilles Birds_current_noAliens.csv",
             "Baiser et al (2017) Marianas Birds_current_noAliens.csv",
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
  if (!all(rainR == lf[26:48])){
    stop("every time: B")
  }

##Terrestrial passerine alien species
lfA <- LANDBIRD_DATA[[3]]#filenames
ldf_alien <- LANDBIRD_DATA[[4]]#alien datasets

##Terrestrial passerine extinct species
lfEx <- LANDBIRD_DATA[[5]]#filenames
ldf_Ex <- LANDBIRD_DATA[[6]]#historic period datasets
}
