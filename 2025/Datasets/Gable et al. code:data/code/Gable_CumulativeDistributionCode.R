#' # Wolf Predation, Beaver Foraging, and Forest Impacts
#'
#'
#'load packages
#+ warning = FALSE, message = FALSE
library(tidyverse)
library(lme4)
library(knitr)
library(berryFunctions)
library(ggeffects)
library(ggbeeswarm)
library(ggstance)
library(nlme)
library(performance)
library(pammtools) # For step ribbon plot

#' Read in data
deploy<-read.csv(file="data/BeaverCam_DeployData.csv") #data on trail camera deployments
act<-read.csv(file="data/BeaverActivity_Data.csv") #photo data on terrestrial beaver activity
beavKills<-read.csv(file="data/BeaverKills_Dat.csv") #beaver kill site data
bha<-read.csv(file="data/Ambush_Dat.csv") #beaver hunting attempt data
forest<-read.csv(file="data/Forest_colony_data.csv") #forest area vs. increasing distance from pond
finTrails<-read.csv(file="data/TypicalTrailData_clean.csv") #data on typical trails
seas<-read.csv(file="data/FeedingTrail_SurveyData_Cleaned.csv")

#' Define cumulative distribution function
ecdf_func <- function(data, uniquevals) {
  # uniquevals are the values at which we want to estimate the ECDF
  # if not supplied, we will use all unique values in the data set
  LengthU <- length(uniquevals)
  nobs<-length(data)
  sorted <- sort(data)
  
  ECDF <- rep(0, LengthU)
  for (i in 1:LengthU) {
    ECDF[i] <- sum(sorted <= uniquevals[i]) / nobs
  }
  return(data.frame(Trail.Length = uniquevals, ECDF = ECDF))
}

# For time
ecdf_func2 <- function(timeontrail, traillength, uniquevals) {
  # uniquevals are the values at which we want to estimate the ECDF
  # if not supplied, we will use all unique values in the data set
  
  LengthU <- length(uniquevals)
  x<-order(traillength)
  
  # Sort the trail lengths and time on trail by trail length
  traillength<-traillength[x]
  timeontrail<-timeontrail[x]

  ECDF <- rep(0, LengthU)
  tottime<-sum(timeontrail) # total time on trail
  for (i in 1:LengthU) {
    ECDF[i] <- sum(timeontrail[traillength <= uniquevals[i]]) / tottime
  }
  return(data.frame(Trail.Length = uniquevals, ECDF = ECDF))
}

#' Start by characterizing the ECDF of trail lengths for wolf kills
uniquevals <- sort(unique(beavKills$feedTrailLength))
length.wolfkills <- ecdf_func(beavKills$feedTrailLength,  
                              uniquevals) 

#' Measure uncertainty with a cluster level bootstrap, resampling wolves
#' with replacement, then repeatedly calculate the ECDF
uwolves <- unique(beavKills$wolfID) # unique id for each wolf
nwolves <- length(uwolves) # number of wolves represented

# Now, bootstrap
nBoots <- 5000 # number of bootstraps
ecdf.ests <- matrix(NA, ncol=length(uniquevals), nrow=nBoots) # to hold bootstrap estimates
for(i in 1:nBoots){
  # create bootstrap
  bootIDs <- data.frame(wolfID = sample(x = uwolves, size = nwolves, replace = TRUE))
  bootDat <- merge(bootIDs, beavKills)
  
  # Estimate ecdf
  ecdf.ests[i,] <- ecdf_func(bootDat$feedTrailLength,  
                                uniquevals)[,2] 
}
length.wolfkills$LCL <-apply(ecdf.ests, 2, quantile, prob = 0.025)
length.wolfkills$UCL<-apply(ecdf.ests, 2, quantile, prob = 0.975)

  
#'  Repeat for ambush kills
#' Start by characterizing the ECDF of trail lengths for ambushes 
uniquevals2 <- sort(unique(bha$trailLength))
length.wolfambush<- ecdf_func(bha$trailLength,  
                              uniquevals = uniquevals2) 

#' Measure uncertainty with a cluster level bootstrap, resampling wolves
#' with replacement, then repeatedly calculate the ECDF
uwolves2 <- unique(bha$wolf) # unique id for each wolf
nwolves2 <- length(uwolves2) # number of wolves represented

# Now, bootstrap
ecdf.ests2 <- matrix(NA, ncol=length(uniquevals2), nrow=nBoots) # to hold bootstrap estimates
for(i in 1:nBoots){
  # create bootstrap
  bootIDs <- data.frame(wolf = sample(x = uwolves2, size = nwolves, replace = TRUE))
  bootDat <- merge(bootIDs, bha)
  
  # Estimate ecdf
  ecdf.ests2[i,] <- ecdf_func(bootDat$trailLength,  
                             uniquevals2)[,2] 
}
length.wolfambush$LCL <-apply(ecdf.ests2, 2, quantile, prob = 0.025)#, na.rm=TRUE)
length.wolfambush$UCL<-apply(ecdf.ests2, 2, quantile, prob = 0.975)#, na.rm=TRUE)


#' Next, creating ECDF for typical trails
#' in the TypicalTrailData_clean data set using the same methods as above
finTrails$type<-"Reference"
y<-data.frame(finTrails$type,finTrails$trailLength, finTrails$seasonIdent) #pulling out trail length and season
y$wolfID<-NA #creating column for combining with kill/ambush data below
y$lodgeNumber<-finTrails$lodgeNumber #adding lodge number to data 
colnames(y)<-c("type","trailLength","season", "wolfID", "lodgeNumber")


#'  Repeat for typical trails
#' Start by characterizing the ECDF of trail lengths 
uniquevals3 <- sort(unique(finTrails$trailLength))
length.typical<- ecdf_func(finTrails$trailLength,  
                              uniquevals = uniquevals3) 

#' Measure uncertainty with a cluster level bootstrap, resampling lodges
#' with replacement, then repeatedly calculate the ECDF
ulodge <- unique(finTrails$lodgeNumber) # unique id for each lodge
nlodge <- length(ulodge) # number of lodges represented

# Now, bootstrap
ecdf.ests3 <- matrix(NA, ncol=length(uniquevals3), nrow=nBoots) # to hold bootstrap estimates
for(i in 1:nBoots){
  # create bootstrap
  bootIDs <- data.frame(lodgeNumber = sample(x = ulodge, size = nlodge, replace = TRUE))
  bootDat <- merge(bootIDs, finTrails)
  
  # Estimate ecdf
  ecdf.ests3[i,] <- ecdf_func(bootDat$trailLength,  
                              uniquevals3)[,2] 
}
length.typical$LCL <-apply(ecdf.ests3, 2, quantile, prob = 0.025)#, na.rm=TRUE)
length.typical$UCL<-apply(ecdf.ests3, 2, quantile, prob = 0.975)#, na.rm=TRUE)


#' OK, and now to try to deal with the number of trips and time foraging
#' as a function of trail length...using the camera trap data 

#' Create the data from the models.
#'re-organizing the data to get the number of active trails for each lodge
t<-seas %>% group_by(lodgeNumber,season) %>% summarize(activeTrails = n())
t$season[t$season=="fallStatus"]<-"Fall" #re-labeling data categories
t$season[t$season=="springStatus"]<-"Spring"
t$season[t$season=="summerStatus"]<-"Summer"
colnames(t)[1]<-c("lodgeNumber.x") #changed name of variable to do the join below

#'combined deployment data with activity data to pair the two datasets
#'filtering out observations 1) not on feeding trails, and 2) where we did not observe beaver leave and return to water (i.e., inOutWater=="NO")
joined<-act %>% left_join(deploy, by = "camDeployID") %>% filter(ftLength>0) #removing trails that did not get measured

#'combining dataset again via a join and getting rid of trails that did not get measured
joinedALT <- joined %>% left_join(t, by=c("lodgeNumber.x","season")) %>% filter(activeTrails!="NA") #removing trails that did not get measured

#'selecting the columns we want for our analysis
cleanDat<-joinedALT %>% dplyr::select(c('beavCaptID','camDeployID','lodgeNumber.x','lodgeType.x',
                                        'inOutWater','timeOnLand','Date.Deployed','Date.Retreived','daysDeployed','season','FT_ID',
                                        'ftLength','activeTrails'))


#'summarizing data
beavers<-cleanDat %>% group_by(season,FT_ID) %>% summarize(avgTimeLand = mean(timeOnLand), #average time on land
                                                           tripsLand= n(), #number of trips on land
                                                           daysDeploy = mean(daysDeployed),
                                                           distance = mean(ftLength), #distance of the trail
                                                           lodgeID = mean(lodgeNumber.x),
                                                           numTrails=mean(activeTrails))

 

#' first model: relationship between number of foraging trips and length of feeding trail.
#'
#'kept getting error stating "model is nearly unidentifiable: large eigenvalue ratios-Rescale variables?".
#' so we centered and scaled "distance" and "numTrails" and that got rid of error.
#'
#'A quick function to center and scale data
centerScale<-function(x=""){
  mux <- mean(x, na.rm = T)
  # calculate the standard deviation of income
  sdx <- sd(x, na.rm = T)
  # centering
  tr1 <- x - mux
  # scaling
  tr2 <- tr1/sdx
  return(tr2)
}

beavers$scaledDist<-centerScale(beavers$distance) #using function to scale variable
beavers$scaledNumTrails<-centerScale(x=beavers$numTrails) #using function to scale variable

#'first model: number of foraging trips model
tripModel<-glmer.nb(tripsLand ~ scaledDist + season + scaledNumTrails + (1|lodgeID)+offset(log(daysDeploy)), beavers)

#' Consider simplified version without season, numb trails
tripModel2<-glmer.nb(tripsLand ~ distance +  (1|lodgeID)+offset(log(daysDeploy)), beavers)

AIC(tripModel, tripModel2)


#' second model: relationship between foraging trip duration and length of feeding trail
timeModel<-lmer(log(avgTimeLand) ~ distance + numTrails + season + (1|lodgeID)+(1|lodgeID:season), beavers)

#' Again, consider simplified model without numbtrails or season
timeModel2<-lmer(log(avgTimeLand) ~ distance + (1|lodgeID), beavers)
AIC(timeModel, timeModel2) #again, this looks like a good move (results in lower AIC)

#' Next steps: the full bootstrap
#' 
#' 1. Bootstrap the two models to get new coefficients 
#' 2. Resample feeding trails with replacement (using a colony-level / lodge bootstrap).
#' 3. Using the bootstrapped coefficients, simulate:
#'     a. number of trips for trails with the lenghts recorded in [2]
#'     b. the time spent on trials with lenghts recorded in [2]
#'     c. call the ecdf function to record time as a function of trail length

# Now, bootstrap
nboots <- 1000
ecdf.ests4 <- matrix(NA, ncol=length(uniquevals3), nrow=nboots) # to hold bootstrap estimates

for(i in 1:nboots){
  # create bootstrap of trail lengths
  bootIDs <- data.frame(lodgeNumber = sample(x = ulodge, size = nlodge, replace = TRUE))
  bootDat <- merge(bootIDs, finTrails)
  
  # use these trail lengths to create a data set for predicting time and number of trips
  newdat<- data.frame(distance = bootDat$trailLength, daysDeploy = 365, lodgeID = bootDat$lodgeNumber)

  # now predictions for time and number of trips for the bootstrapped trail lengths
  # use bootMer to account for uncertainty in the parameters of the trip and trail models
  boot.predtime <- exp(as.numeric(bootMer(timeModel2,
                       FUN=function(x){simulate(x, newdata=newdat, re.form=NA, allow.new.levels=T)$sim_1},
                       nsim=1)$t))
  boot.predtrips <- as.numeric(bootMer(tripModel2,
                           FUN=function(x){simulate(x, newdata=newdat, re.form=NA, allow.new.levels=T)$sim_1},
                           nsim=1)$t)
  
  # Total time = time per trip*number of trips
  boot.tottime<-boot.predtime*boot.predtrips
  
  # Estimate ecdf using second ecdf fucntion
  ecdf.ests4[i,] <- ecdf_func2(boot.tottime,  
                               newdat$distance, uniquevals = uniquevals3)[,2] 
}
# Summarize the mean and 0.025 and  0.975 percentiles of the ecdfs
time.typical<-data.frame(Trail.Length = uniquevals3) 
time.typical$ECDF<-apply(ecdf.ests4, 2, mean) 
time.typical$LCL <-apply(ecdf.ests4, 2, quantile, prob = 0.025) 
time.typical$UCL<-apply(ecdf.ests4, 2, quantile, prob = 0.975) 


#' Now, repeat the plot, but swapping out the last cumulative distribution
#' function.
#' 
#' - the population of trails associated with beaver lodges
#' - the population of trails where wolves sit and wait (i.e. ambush)
#' - the population of trails where wolves kill beavers 


#+out.width = "90%", fig.height = 5, fig.width =8
cumBeav<-ggplot(data=length.wolfkills, aes(Trail.Length, ECDF)) +
  geom_ribbon(data=time.typical, aes(ymin = LCL, ymax = UCL), fill = "#D3C225",
                  alpha = 0.2)  +
  geom_line(data=time.typical, aes(Trail.Length, ECDF, color = "Beaver foraging time"),
             size = 2) + 
  geom_ribbon(data=length.wolfkills, aes(ymin = LCL, ymax = UCL), fill = "#39568CFF",
                   alpha = 0.2)  +
  geom_line(aes(color = "Kill trails"), size = 2) +
  xlab("Feeding trail length (m)") + ylab("Proportion of total")+
  coord_cartesian(xlim=c(0,70), ylim=c(0,1.02), expand = c(0,0))+
  scale_color_manual(values = c("Beaver foraging time" = '#D3C225', "Kill trails"='#39568CFF'))+
  theme_classic()+
  theme(axis.title.y = element_text(size=24, family = "Gill Sans"), #editing font size of x-axis title
        axis.title.x = element_text(size=24, family = "Gill Sans"),
        axis.text.x = element_text(family = "Gill Sans",size=18),
        axis.text.y = element_text(family = "Gill Sans",size=18),
        legend.title = element_blank(),
        legend.position=c(0.8,0.15),
        legend.text = element_text(size=20, family="Gill Sans"))
 
cumBeav

ggsave("Fig5_ECDF_V1.png",plot=cumBeav, height=8, width=8 , dpi = 300)
   



