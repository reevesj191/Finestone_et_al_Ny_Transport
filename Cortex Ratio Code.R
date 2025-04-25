library(ggplot2)
citation("ggplot2")
library(googlesheets4)
citation("googlesheets4")
library(ggthemes)
citation("ggthemes")
library(dplyr)
citation("dplyr")

### Lithic Analysis functions
CoreSurfaceArea <- function(NYCORES, method = "Scalene"){
  
  if(method == "Scalene"){
    NYCORES$semiaxis1 <- NYCORES$LENGTH_cm/2
    NYCORES$semiaxis2 <- NYCORES$WIDTH_cm/2
    NYCORES$semiaxis3 <- NYCORES$THICKNESS_cm/2
    
    NYCORES$SA <- 4*pi*((((NYCORES$semiaxis1 ^ 1.6075)*(NYCORES$semiaxis2 ^1.6075)+(NYCORES$semiaxis1 ^1.6075)*(NYCORES$semiaxis3 ^1.6075)+(NYCORES$semiaxis2^1.6075)*(NYCORES$semiaxis3^1.6075))/3)^(1/1.6075))
    return(NYCORES$SA)
  }
}


cortex.ratio <- function(cores, detached, avg_nodule_mass = FALSE){
  require(dplyr)
  detached_group_variables <- select(detached, MASS_g, volume, cortical_SA)
  detached_group_measure <- as.data.frame(summarise(detached_group_variables,
                                      detached_mass = sum(MASS_g, na.rm = TRUE),
                                      detached_volume = sum(volume,na.rm = TRUE),
                                      detached_cortical_surface_area = sum(cortical_SA, na.rm = TRUE), ))

detached_group_measure
  cores_group_variables <- select(cores, MASS_g, volume, cortical_SA)
  cores_group_measure <- as.data.frame(summarise(cores_group_variables,
                                   cores_mass = sum(MASS_g,na.rm = TRUE),
                                   cores_volume = sum(volume,na.rm = TRUE),
                                   cores_cortical_surface_area = sum(cortical_SA,na.rm = TRUE), 
                                   nodule_volume_quartile = quantile(volume, na.rm = TRUE)[4],
                                   core_count = length(MASS_g)))

  detached_group_measure$merger <- detached_group_measure[,1]
  cores_group_measure$merger <- cores_group_measure[,1]

cortex_ratio <- merge(cores_group_measure, detached_group_measure, by = c("merger"))
  
cortex_ratio$assemblage_volume <- cortex_ratio$detached_volume + cortex_ratio$cores_volume
  
cortex_ratio$nodule_vol_core_count <- cortex_ratio$assemblage_volume / cortex_ratio$core_count
  
cortex_ratio$nodule_frequency_quartile <- cortex_ratio$assemblage_volume/cortex_ratio$nodule_volume_quartile
  
cortex_ratio$nodule_SA_Quartile <- ((((3*(cortex_ratio$nodule_volume_quartile))/(4*pi))^(2/3))*4*pi)
  
cortex_ratio$nodule_SA_core_count <- ((((3*(cortex_ratio$nodule_vol_core_count))/(4*pi))^(2/3))*4*pi)
  
cortex_ratio$observed_cortex <- cortex_ratio$detached_cortical_surface_area + cortex_ratio$cores_cortical_surface_area  

cortex_ratio$expected_cortical_SA_quartile <- cortex_ratio$nodule_SA_Quartile * cortex_ratio$nodule_frequency_quartile

cortex_ratio$expected_cortical_SA_core_count <- cortex_ratio$nodule_SA_core_count * cortex_ratio$core_count

cortex_ratio$cortex_ratio_quartile <- cortex_ratio$observed_cortex / cortex_ratio$expected_cortical_SA_quartile
  
cortex_ratio$cortex_ratio_core_count <- cortex_ratio$observed_cortex / cortex_ratio$expected_cortical_SA_core_count
  
return(cortex_ratio)
}

#cores <- group_by(EFT.Cores, area)
#detached <- group_by(EFT.Flakes, area)


volume.ratio <- function(cores, detached, avg_nodule_mass = FALSE){
  
  require(dplyr)
  detached_group_variables <- select(detached, MASS_g, volume, cortical_SA)
  detached_group_measure <- as.data.frame(summarise(detached_group_variables,
                                                    detached_mass = sum(MASS_g, na.rm = TRUE),
                                                    detached_volume = sum(volume,na.rm = TRUE),
                                                    detached_cortical_surface_area = sum(cortical_SA, na.rm = TRUE)))
  
 detached_group_measure
  
  cores_group_variables <- select(cores, MASS_g, volume, cortical_SA)
  cores_group_measure <- as.data.frame(summarise(cores_group_variables,
                                                 cores_mass = sum(MASS_g,na.rm = TRUE),
                                                 cores_volume = sum(volume,na.rm = TRUE),
                                                 cores_cortical_surface_area = sum(cortical_SA,na.rm = TRUE), 
                                                 nodule_volume_quartile = quantile(volume, na.rm = TRUE)[4],
                                                 core_count = length(MASS_g)))
  
 cores_group_measure
  
  detached_group_measure$merger <- detached_group_measure[,1]
  cores_group_measure$merger <- cores_group_measure[,1]
  
 volume_ratio <- merge(cores_group_measure, detached_group_measure, by = c("merger"))
  
  
  
##### Calculation of the Observed Volume
volume_ratio$assemblage_volume <- volume_ratio$detached_volume + volume_ratio$cores_volume
  
#### Estimation of nodules if you don't just use cortex.
volume_ratio$nodule_frequency_quartile <- volume_ratio$assemblage_volume/volume_ratio$nodule_volume_quartile
  
volume_ratio$expected_volume_core_count <-volume_ratio$nodule_volume_quartile * volume_ratio$core_count
  
volume_ratio$volume_ratio_core_count <- volume_ratio$expected_volume_core_count / volume_ratio$assemblage_volume
  
volume_ratio
return(volume_ratio)
}


### Load Nyayanga data
NYDETACHED <- read.csv("~/Desktop/NYADETACHED.csv", header = TRUE, colClasses = c("factor", "factor", "factor", "numeric", "numeric", "numeric", "numeric", "numeric"))
summary(NYDETACHED)

NYCORES <- read.csv("~/Desktop/NYACORES.csv", header = TRUE, colClasses = c("factor", "factor", "factor", "numeric", "numeric", "numeric", "numeric", "numeric"))

#convert everything to centimeters IF NECESSARY
NYCORES$LENGTH_cm <- NYCORES$LENGTH_mm * .1
NYCORES$WIDTH_cm <- NYCORES$WIDTH_mm * .1
NYCORES$THICKNESS_cm <- NYCORES$THICKNESS_mm * .1

NYDETACHED$LENGTH_cm <- NYDETACHED$LENGTH_mm * .1
NYDETACHED$WIDTH_cm <- NYDETACHED$WIDTH_mm * .1
NYDETACHED$THICKNESS_cm <- NYDETACHED$THICKNESS_mm * .1

#calculate volume
NYCORES$volume<-NYCORES$MASS_g/2.5
NYDETACHED$volume<-NYDETACHED$MASS_g/2.5

#calculatesurfaceArea
NYDETACHED$SA<-NYDETACHED$LENGTH_cm*NYDETACHED$WIDTH_cm
NYCORES$SA<-CoreSurfaceArea(NYCORES, method="Scalene")

#calculate cortical SA

NYCORES$cortical_SA<-NYCORES$SA*(NYCORES$CORTEX_PERC/100)

NYDETACHED$cortical_SA<-NYDETACHED$SA*(NYDETACHED$CORTEX_PERC/100)

#group the data according to FP and DP
NYCORES_GRP<-group_by(NYCORES,Context)
NYDETACHED_GRP<-group_by(NYDETACHED,Context)

##Calculate cortex ratio
NYCORTEXRATIO<-cortex.ratio(cores=NYCORES_GRP, detached=NYDETACHED_GRP)

NYCORTEXRATIO[,c(1,20,21)]

#Resample
n<-10000
group<-"Nonlocal"

CORES<-subset(NYCORES, Context== group)
DETACHED<-subset(NYDETACHED, Context==group)

Core.N<-round(nrow(CORES)*.75)
Detached.N<-round(nrow(DETACHED)*.75)

desired_length <- n
Resample.results <- vector(mode = "list", length = desired_length)

for (i in 1:n){
  
  core.sam<-CORES[sample(1:nrow(CORES),Core.N),]
  flake.sam<-DETACHED[sample(1:nrow(DETACHED),Detached.N),]
  core.sam<-group_by(core.sam,Context)
  flake.sam<-group_by(flake.sam,Context)
  iter<-cortex.ratio(cores=core.sam,detached=flake.sam)
  iter<-iter[,c(1,20,21)]
  Resample.results[[i]]<-(iter)
}

Resample <- do.call("rbind", Resample.results)

Resample_nonlocal<-Resample

plot1<-ggplot(Resample, aes(x = cortex_ratio_quartile)) + 
  geom_density()+ 
  labs(title = paste("Distribution of Cortex Ratio Values for", group)) +
  theme_bw()

plot1

group<-"Local"

CORES<-subset(NYCORES, Context== group)
DETACHED<-subset(NYDETACHED, Context==group)

Core.N<-round(nrow(CORES)*.75)
Detached.N<-round(nrow(DETACHED)*.75)

desired_length <- n
Resample.results <- vector(mode = "list", length = desired_length)

for (i in 1:n){
  
  core.sam<-CORES[sample(1:nrow(CORES),Core.N),]
  flake.sam<-DETACHED[sample(1:nrow(DETACHED),Detached.N),]
  core.sam<-group_by(core.sam,Context)
  flake.sam<-group_by(flake.sam,Context)
  iter<-cortex.ratio(cores=core.sam,detached=flake.sam)
  iter<-iter[,c(1,20,21)]
  Resample.results[[i]]<-(iter)
}

Resample <- do.call("rbind", Resample.results)


Resample_local<-Resample

plot2<-ggplot(Resample, aes(x = cortex_ratio_quartile)) + 
  geom_density()+ 
  labs(title = paste("Distribution of Cortex Ratio Values for", group)) +
  theme_bw()

plot2
