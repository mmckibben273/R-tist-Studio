#' # Wolf Predation, Beaver Foraging, and Forest Impacts
#'
#'
#'load packages
library(tidyverse)
library(lme4)
library(car)
library(knitr)
library(ezknitr)
library(MASS)
library(emmeans)
library(berryFunctions)
library(ggeffects)
library(ggbeeswarm)
library(ggstance)
library(nlme)
library(performance)

#' Reading in data
deploy<-read.csv(file="data/BeaverCam_DeployData.csv") #data on trail camera deployments
act<-read.csv(file="data/BeaverActivity_Data.csv") #photo data on terrestrial beaver activity
beavKills<-read.csv(file="data/BeaverKills_Dat.csv") #beaver kill site data
bha<-read.csv(file="data/Ambush_Dat.csv") #beaver hunting attempt data
forest<-read.csv(file="data/Forest_colony_data.csv") #forest area vs. increasing distance from pond
finTrails<-read.csv(file="data/TypicalTrailData_clean.csv") #data on typical trails
seas<-read.csv(file="data/FeedingTrail_SurveyData_Cleaned.csv")


#'deployment summary
deploy_sum<-deploy

#'just getting some summary deployment data
sum(deploy_sum$picsTaken, na.rm=T) #total pics taken
sum(deploy_sum$beavCapts, na.rm=T) #number of beaver events


#' number of rows <15 m and <30 m
nrow(finTrails %>% filter(trailLength<15))/nrow(finTrails)*100
nrow(finTrails %>% filter(trailLength<30))/nrow(finTrails)*100

#number of diff colonies studied
length(unique(finTrails$lodgeNumber))

#'
#'
#' ### Kill, ambush, and survey trail data comparison
#'
#'
#'organizing/pulling out data for kill trails
beavKills$type<-"Kill"
x<-data.frame(beavKills$type,beavKills$feedTrailLength, beavKills$season,beavKills$wolfID) #pulling out trail length, season, and wolf id
x$lodgeNumber<-NA #creating column for combining with beaver trail data below
colnames(x)<-c("type","trailLength", "season","wolfID", "lodgeNumber")

#'organizing/pulling out reference trail data
finTrails$type<-"Reference"
y<-data.frame(finTrails$type,finTrails$trailLength, finTrails$seasonIdent) #pulling out trail length and season
y$wolfID<-NA #creating column for combining with kill/ambush data below
y$lodgeNumber<-finTrails$lodgeNumber #adding lodge number to data 
colnames(y)<-c("type","trailLength","season", "wolfID", "lodgeNumber")

#'organizing/pulling out ambush data
bha$type<-"Ambush"
z<-data.frame(bha$type,bha$trailLength, bha$season, bha$wolf) #pulling out trail length, season, and wolf id
z$lodgeNumber<-NA #creating column for combining with beaver trail data below
colnames(z)<-c("type","trailLength", "season","wolfID","lodgeNumber")

#'combining the kill, ambush, and reference trail data
trailDat<-rbind(x,y,z)

#'reordering data for plotting
trailDat$type<-ordered(trailDat$type, levels = c( "Kill","Ambush", "Reference"))

#'
#'
#'
#' ### Visualizing raw data
#' 
#' 
#' 

#'violin plot showing distribution of the data
violinPlot<-ggplot(trailDat, aes(x=trailLength, y=type, color=type, fill=type))+
  geom_violin(scale="width", alpha=.4, width=1, adjust=.8, size=1.5)+
  geom_quasirandom(method="pseudorandom", alpha=0.5, size=3)+theme_classic()+
  theme(axis.title.x = element_text(size=23, family = "Gill Sans"),
        axis.title.y = element_text(size=23, family = "Gill Sans"), #editing font size of x-axis title
        axis.text.x = element_text(family = "Gill Sans", size=16),
        axis.text.y = element_text(family = "Gill Sans",size=16))+xlab("Feeding trail length (m)")+ylab("")+
  theme(legend.position="none")+scale_x_continuous(limits = c(0,80),expand = c(0,2))+
  stat_summary(fun="median", geom="crossbar", size=0.5, color="black", linetype="longdash", show.legend=F)+
  scale_fill_manual(values = c("Reference" = "#D3C225", "Ambush"="#20A387FF", "Kill" = "#39568CFF"),
                    guide = guide_legend(reverse = TRUE))+
  scale_color_manual(values = c("Reference" = "#D3C225", "Ambush"="#20A387FF", "Kill" = "#39568CFF"),
                     guide = guide_legend(reverse = TRUE))+
  theme(axis.title.x = element_text(size=28, family = "Gill Sans"),
        axis.title.y = element_text(size=28, family = "Gill Sans"), #editing font size of x-axis title
        axis.text.x = element_text(family = "Gill Sans", size=20),
        axis.text.y = element_text(family = "Gill Sans",size=0),
        legend.title=element_blank(),
        legend.text=element_text(size=24, family="Gill Sans"),
        legend.position = c(0.85,0.90))

violinPlot

#' exporting figure
ggsave("TrailDensityPlot_V3.png",plot=violinPlot, height=8, width=10 , dpi = 300)



#'Modeling typical feeding trail length

#' Fit reference trail model (used the finTrails df so that we could include lodgeType for covariate)
beavMod <-lmer(log(trailLength) ~ seasonIdent + lodgeType + (1|lodgeNumber)+(1 |lodgeNumber:seasonIdent), data=finTrails)
summary(beavMod)
confint(beavMod)

#'comparing seasonal values
emmeans(beavMod,specs=pairwise~seasonIdent, type="response", back.transform=F)

#' seeing if model meets assumptions
check_model(beavMod)

#'predicting values on original scale
BM<-ggemmeans(beavMod,terms=c("seasonIdent"), type="fixed",back.transform = T)
BM$group<-"Reference"

#' Fit kill/ambush trail model
wolf <- trailDat %>% filter(type !="Reference")
wolf$type<-factor(wolf$type, levels=c("Ambush","Kill"))
lmermod <-lmer(log(trailLength) ~ season*type  + (1|wolfID), data=wolf)
summary(lmermod)
confint(lmermod)

#'predicting values on original scale
KA<-ggemmeans(lmermod,terms=c("season","type"),type="fixed", back.transform = T)

#' seeing if model meets assumptions
check_model(lmermod)

#'combining wolf and beaver data for plotting
kabm<-rbind(KA,BM)
kabm$x <- ordered(kabm$x, levels = c("Spring", "Summer", "Fall")) #ordering data for plotting
kabm$group <- ordered(kabm$group, levels = c("Reference", "Ambush", "Kill")) #ordering data for plotting


#'plotting differences in trail length across seasons 
trailFig<-ggplot(kabm, aes(y=predicted, x=group, color=group, group=group))+geom_point(size=6)+geom_errorbar(aes(ymin=conf.low,ymax=conf.high), size=1.0, width=0.5)+
  facet_wrap(~x) + scale_color_manual(values = c("Reference" = "#E3D500", "Ambush"="#20A387FF", "Kill" = "#39568CFF"))+
  ylab("Feeding trail length")+theme(axis.title.y = element_text(size=24, family = "Gill Sans"), #editing font size of x-axis title
                                 axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_text(family = "Gill Sans",size=20),
                                 strip.text.x = element_text(size = 18, family = "Gill Sans"), #control facet text label
                                 legend.position=c(0.9,0.15),
                                 legend.title = element_blank(),
                                 legend.text = element_text(size=24, family="Gill Sans"))+scale_y_continuous(limits = c(5,50),expand = c(0,0))

trailFig

ggsave("FeedingTrailFig_V2.png",plot=trailFig, height=8, width=12 , dpi = 300)

#'
#'
#'
#' ## The mechanism: why are beavers killed on longer trails?
#'
#'The goal: to examine beaver activity patterns on feeding trails using remote camera data 

#'re-organizing the data to get the number of active trails for each lodge
t<-seas %>% group_by(lodgeNumber,season) %>% summarize(activeTrails = n())
t$season[t$season=="fallStatus"]<-"Fall" #re-labeling data categories
t$season[t$season=="springStatus"]<-"Spring"
t$season[t$season=="summerStatus"]<-"Summer"
colnames(t)[1]<-c("lodgeNumber.x") #changed name of variable to do the join below

#'combined deployment data with activity data to pair the two datasets
#'filtering out observations 1) not on feeding trails, 2) from island colonies, and 3) where we did not observe beaver leave and return to water (i.e., inOutWater=="NO")
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


#' ## Modeling to figure out the "mechanism"
#' 
#' 


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

#'model results
summary(tripModel) 
confint(tripModel)
check_model(tripModel) #checking model assumptions

#'using ggpredict to predict how number of foraging trips changes based on trail length and season
predDF<-ggpredict(tripModel, terms=c("scaledDist [-1.5:5, by=0.1]","season"), type="fixed")
predDF$group <- ordered(predDF$group, levels = c("Spring", "Summer", "Fall"))

#'uncentering/unscaling for plotting
mB <- mean(beavers$distance, na.rm = T) #mean
sdB <- sd(beavers$distance, na.rm = T) #sd
predDF$x2<-(predDF$x*sdB)+mB #back to original scale


#'plotting results of ggpredict() to show how umber of foraging trips changes based on trail length and season
tripFig<-ggplot(predDF, aes(x=x2,y=predicted, group=group, color=group))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=group),linetype=0, alpha=0.2)+
  geom_line(position=position_dodge(width=0.2), size=2)+
  theme_classic() + geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=group), linetype=0, alpha=0.1)+
  scale_color_manual(values = c("Spring" = "#E3D500", "Summer"="#20A387FF", "Fall" = "#39568CFF"))+
  scale_fill_manual(values = c("Spring" = "#E3D500", "Summer"="#20A387FF", "Fall" = "#39568CFF"))+
  ylab("Number of foraging trips")+xlab("Feeding trail length (m)") + 
  theme(axis.title.y = element_text(size=24, family = "Gill Sans"), #editing font size of x-axis title
        axis.title.x = element_text(size=24, family = "Gill Sans"),
        axis.text.x = element_text(family = "Gill Sans",size=20),
        axis.text.y = element_text(family = "Gill Sans",size=20),
        legend.position=c(0.2,0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=24, family="Gill Sans"))+
  coord_cartesian(xlim=c(0,83), ylim=c(0,78), expand = c(0,0))

tripFig
ggsave("TripsLandFig_V2.png",plot=tripFig, height=8, width=8 , dpi = 300)


#' second model: relationship between foraging trip duration and length of feeding trail
timeModel<-lmer(log(avgTimeLand) ~ distance + numTrails + season + (1|lodgeID)+(1|lodgeID:season), beavers)

#'model results
summary(timeModel)
confint(timeModel)
check_model(timeModel) #checking model assumptions

#'using ggpredict to predict how foraging trip duration changes based on trail length and season
predTrip<-ggpredict(timeModel, terms=c("distance","season"), type="fixed")
predTrip$group <- ordered(predTrip$group, levels = c("Spring", "Summer", "Fall"))

#'plotting results of ggpredict() to show how duration of foraging trips changes based on trail length and season
timeFig<-ggplot(predTrip, aes(x=x,y=predicted, color=group))+geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill=group),linetype=0, alpha=0.2)+
  geom_line(position=position_dodge(width=0.5), size=2)+
  theme_classic() + scale_color_manual(values = c("Spring" = "#E3D500", "Summer"="#20A387FF", "Fall" = "#39568CFF"))+
  scale_fill_manual(values = c("Spring" = "#E3D500", "Summer"="#20A387FF", "Fall" = "#39568CFF"))+
  ylab("Duration of foraging trips (min)")+xlab("Feeding trail length (m)") + 
  theme(axis.title.y = element_text(size=24, family = "Gill Sans"), #editing font size of x-axis title
        axis.title.x = element_text(size=24, family = "Gill Sans"),
        axis.text.x = element_text(family = "Gill Sans",size=20),
        axis.text.y = element_text(family = "Gill Sans",size=20),
        legend.position=c(0.85,0.1),
        legend.title = element_blank(),
        legend.text = element_text(size=24, family="Gill Sans"))+
  coord_cartesian(xlim=c(0,83),ylim=c(2,14), expand=c(0,0))


timeFig
ggsave("TimeLandFig_V2.png",plot=timeFig, height=8, width=8 , dpi = 300)


#' ## Plotting forest altered vs. distance from water
forestFig<-ggplot(forest, aes(x=buffer_distance, y=total_area))+scale_x_continuous(expand=c(0,0),limits=c(1,35))+
  scale_y_continuous(expand=c(0,0),limits=c(0,85))+theme_bw()+ylab(bquote("Forest area available "(km^2)))+xlab("Foraging distance from water (m)") + 
  theme(axis.title.y = element_text(size=24, family = "Gill Sans"), #editing font size of x-axis title
        axis.title.x = element_text(size=24, family = "Gill Sans"),
        axis.text.x = element_text(family = "Gill Sans",size=20),
        axis.text.y = element_text(family = "Gill Sans",size=20))+
  annotate("rect", xmin = 11.3, xmax = 14.6, ymin = 0, ymax = 85,alpha = .6,fill ="#E3D500")+ #adds shading for typical trails
  annotate("rect", xmin = 23.1, xmax = 28.8, ymin = 0, ymax = 85,alpha = .6,fill = "#39568CFF")+ #adds shading for kill trails
  annotate("rect", xmin = 18.4, xmax = 20.9, ymin = 0, ymax = 85,alpha = .6,fill = "#20A387FF")+ #adds shading for ambush trails
  geom_line(color="black", size=2)+geom_point(size=5, color="black")
  

forestFig
ggsave("ForestAvailableFig_V2.png",plot=forestFig, height=8, width=12 , dpi = 300)


#' ###Footer
#' 
#'spun with  ezspin(file="WolfBeaver_FT_Analysis.R",out_dir="output", keep_md=F) 
     