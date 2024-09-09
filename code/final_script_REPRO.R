#This script is for the project :
# Gobyrepro
#Here we can find data import and manipulation, plotting and statistical tests 

setwd("C://Users//ivainm//Working_Document//DYNAMAR//Manuscripts//reproductive strategy//Data")

#packages:
library(ggExtra)
library(emmeans)
library(reshape)
library(lsmeans)
library(effects)
library(multcomp)
library(adegenet)
library(ggpubr)
library(vcfR)
library(LEA)
library(xadmix)
library(genepop)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 1. BREEDING WINDOW

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Temperature data management
#Import data from the Nordkyst model
ALLmod<-read.csv("ALLmod.csv")

## Average T per month
ALLmonthly<-data.frame(Area=character(),Month=integer() ,Temperature=numeric(), TCI95=numeric() )
for (i in 1:length(levels(ALLmod$Area))) {
  for(m in 1:12) {
    ALLmonthly[m+(12*(i-1)),]<-c(levels(ALLmod$Area)[i], m,
                                 mean(ALLmod[ALLmod$Area==levels(ALLmod$Area)[i] & ALLmod$Mo==m,]$T2m),
                                 3.92*sqrt(var(ALLmod[ALLmod$Area==levels(ALLmod$Area)[i] & ALLmod$Mo==m,]$T2m)/length(ALLmod[ALLmod$Area==levels(ALLmod$Area)[i] & ALLmod$Mo==m,]$T2m)))
      }}

ALLmonthly$Monthdate<-as.Date(paste("2022",ALLmonthly$Month,"15", sep="/"))
ALLmonthly$Temperature<-as.numeric(ALLmonthly$Temperature)
ALLmonthly$TCI95<-as.numeric(ALLmonthly$TCI95)
ALLmonthly$Area<- factor(ALLmonthly$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad','Northern Limit'))
ALLmonthly$Month<-as.numeric(ALLmonthly$Month)
str(ALLmonthly)

#Fitting temperature curves to the data?
#Convert months to day unit
str(ALLmonthly)
ALLmonthly$Day<-(ALLmonthly$Month-1)*30+15
ALLmonthly$Fit<-NA
#365 days version
ALL2022<-ALLmod[ALLmod$Year=="2022",]
ALLfit<-data.frame(cbind(Area=as.character(ALL2022$Area),Date=as.character(ALL2022$datemonth),Day=NA,Temperature=NA))
ALLfit$Date<-as.Date(ALLfit$Date)
ALLfit$Day<-as.numeric(ALLfit$Day)
ALLfit$Temperature<-as.numeric(ALLfit$Temperature)

#Modst is where I store the model parameters for fitted temperature curves
MODST<-data.frame(cbind(Area=as.character(),Intercept=as.numeric(),sin=as.numeric(),
                        cos=as.numeric(),sin2=as.numeric(),cos2=as.numeric()))
for(i in 1:7){
  modx<-lm(data=ALLmonthly[ALLmonthly$Area==levels(ALLmonthly$Area)[i],], 
           Temperature ~
             sin(2*pi/365*Day)+
             cos(2*pi/365*Day)+
             sin(4*pi/365*Day)+ #adding second harmonic, needed.
             cos(4*pi/365*Day)
  );
  MODST<-rbind(MODST,c(levels(ALLmonthly$Area)[i],as.numeric(coef(modx))))
#Feed fitted values to the monthly data
  ALLmonthly[ALLmonthly$Area==levels(ALLmonthly$Area)[i],]$Fit<-fitted(modx)
  ALLfit[ALLfit$Area==levels(ALLmonthly$Area)[i],]$Day<-c(1:length(ALLfit[ALLfit$Area==levels(ALLmonthly$Area)[i],]$Area))
  ALLfit[ALLfit$Area==levels(ALLmonthly$Area)[i],]$Temperature<-predict(modx,newdata=ALLfit[ALLfit$Area==levels(ALLmonthly$Area)[i],])}
ALLfit$Area<- factor(ALLfit$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad','Northern Limit'))
#adjust year from ALLfit
year(ALLfit$Date)<-2022

#Fromatting MODST
colnames(MODST)<-c("Area","Intercept","sin","cos","sin2","cos2")
MODST$Intercept<-as.numeric(MODST$Intercept)
MODST[,3]<-as.numeric(MODST[,3])
MODST[,4]<-as.numeric(MODST[,4])
MODST[,5]<-as.numeric(MODST[,5])
MODST[,6]<-as.numeric(MODST[,6])
MODST$Area<- factor(MODST$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad','Northern Limit'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Color scheme for populations
popcols<-c("red3",  "brown1", "darkorchid4",  "mediumpurple1", "blue3",  "steelblue1", "grey15")
#Figure temperature and primary produtivity (Fig 2.a)
Aenv=ggplot(data=ALLmonthly, aes(x=Monthdate, y=Temperature, color=Area))+
  annotate(geom = "rect",xmin = as.Date("2022-04-24"),
           xmax = as.Date("2022-07-25"),
           ymin = -Inf,
           ymax = Inf,
           fill="blue", alpha =0.05) +
  geom_line(size=1, alpha=0.5)+
  geom_point(size=3)+
  scale_x_date(date_labels = "%B",breaks="month")+
  #geom_line(data=ALLfit,aes(y=Temperature,x=Date),size=2, alpha=0.15)+
  geom_line(data=ALLfit[ALLfit$Area=="Kristineberg" & ALLfit$Date>as.Date("2022-02-23") & ALLfit$Date<as.Date("2022-10-20") ,],aes(y=Temperature,x=Date),linewidth=3.5, alpha=0.35)+
  geom_line(data=ALLfit[ALLfit$Area=="Arendal" & ALLfit$Date>as.Date("2022-02-25") & ALLfit$Date<as.Date("2022-10-19") ,],aes(y=Temperature,x=Date),linewidth=3.5, alpha=0.35)+
  geom_line(data=ALLfit[ALLfit$Area=="Austevoll" & ALLfit$Date>as.Date("2022-03-08") & ALLfit$Date<as.Date("2022-10-02") ,],aes(y=Temperature,x=Date),linewidth=3.5, alpha=0.35)+
  geom_line(data=ALLfit[ALLfit$Area=="Hitra" & ALLfit$Date>as.Date("2022-03-12") & ALLfit$Date<as.Date("2022-09-06") ,],aes(y=Temperature,x=Date),linewidth=3.5, alpha=0.35)+
  geom_line(data=ALLfit[ALLfit$Area=="Helligvaer" & ALLfit$Date>as.Date("2022-03-28") & ALLfit$Date<as.Date("2022-08-16") ,],aes(y=Temperature,x=Date),linewidth=3.5, alpha=0.35)+
  geom_line(data=ALLfit[ALLfit$Area=="Ringstad" & ALLfit$Date>as.Date("2022-04-01") & ALLfit$Date<as.Date("2022-10-02") ,],aes(y=Temperature,x=Date),linewidth=3.5, alpha=0.35)+
  geom_line(data=ALLfit[ALLfit$Area=="Northern Limit" & ALLfit$Date>as.Date("2022-04-02") & ALLfit$Date<as.Date("2022-09-02") ,],aes(y=Temperature,x=Date),linewidth=3.5, alpha=0.35)+
  
  
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-02-23") & ALLfit$Area=="Kristineberg",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=24, color="red3", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-02-25") & ALLfit$Area=="Arendal",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=24, color="brown1", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-03-08") & ALLfit$Area=="Austevoll",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=24, color="darkorchid4", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-03-12") & ALLfit$Area=="Hitra",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=24, color="mediumpurple1", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-03-28") & ALLfit$Area=="Helligvaer",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=24, color="blue3", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-04-01") & ALLfit$Area=="Ringstad",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=24, color="steelblue1", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-04-02") & ALLfit$Area=="Northern Limit",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=24, color="grey15", fill="white")+
  
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-10-20") & ALLfit$Area=="Kristineberg",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=25, color="red3", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-10-19") & ALLfit$Area=="Arendal",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=25, color="brown1", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-10-02") & ALLfit$Area=="Austevoll",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=25, color="darkorchid4", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-09-06") & ALLfit$Area=="Hitra",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=25, color="mediumpurple1", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-08-16") & ALLfit$Area=="Helligvaer",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=25, color="blue3", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-10-02") & ALLfit$Area=="Ringstad",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=25, color="steelblue1", fill="white")+
  geom_point(data= ALLfit[ALLfit$Date==as.Date("2022-09-02") & ALLfit$Area=="Northern Limit",],aes(y=Temperature,x=Date),stroke=1.5, size=4, shape=25, color="grey15", fill="white")+
  
  geom_errorbar(size=1,width=2,alpha=0.6,aes(ymin = Temperature - TCI95, 
                                             ymax = Temperature + TCI95))+
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  scale_colour_manual("Population",values =popcols)+
  labs(x="Month",y=expression("Temperature " ( degree*C)),title="a.")

#Importing data for Figure 2.b and 2.c
Env2<-read.csv("env_info.txt", header = T)
Env2$Area<- factor(Env2$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad','Northern Limit'))
Coord<-read.table("Coord.txt", header=T)
Coord[Coord$Area=="Northern_limit",]$Area<-"Northern Limit"
Coord$Area<- factor(Coord$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad','Northern Limit'))
EnvD<-merge(Env2,Coord,by="Area")

#Figure 2.b
Benv<-ggplot(EnvD, aes(x=Lat))+geom_point(aes(y=Ddays, color=Area), shape=21, size=3,stroke=1.5, fill="white", alpha=0.5)+
  geom_point(aes(y=Ddaysfixed, color=Area), shape=22, size=3,stroke=1.5, fill="white")+
  scale_colour_manual("Population",values = popcols)+
  theme(plot.margin = margin(0.1,1,0.1,1, "cm")) +
  labs(x="Latitude (N)",y="Degree days ",title="b.")
#Figure 2.c
Cenv<-ggplot(EnvD, aes(x=Lat))+geom_point(aes(y=Broods, color=Area), shape=21, size=3,stroke=1.5, fill="white", alpha=0.5)+
  geom_point(aes(y=Broodsfixed, color=Area), shape=22, size=3,stroke=1.5, fill="white")+
  scale_colour_manual("Population",values = popcols)+
  theme(plot.margin = margin(0.1,1.3,0.1,1, "cm")) +
  labs(x="Latitude (N)",y="Standard broods ",title="c.")
#Arrange final figure
smalls=ggpubr::ggarrange(Benv,Cenv, nrow = 1,legend = "none", common.legend = T)
ggpubr::ggarrange(Aenv,smalls, nrow = 2,heights=c(2,1),legend="right")




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 2. Condition, density

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# BODY CONDITION
phenodata<-read.csv("phenodata.csv")
phenodata$roundness2<-as.factor(phenodata$roundness)
phenodata$area<-as.factor(phenodata$area)
phenodata$area<- factor(phenodata$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
phenodata$location<-as.factor(phenodata$location)
phenodata$timepoint<-as.factor(phenodata$timepoint)
phenodata$sex<-as.factor(phenodata$sex)
phenodata$date2<-as.Date(phenodata$date, format="%d.%m")
phenodata$Kfactor<-100*10*(phenodata$weight/phenodata$length^3) #x10 to correct for units
summary(phenodata)
#removing NAs and outliers
hist(phenodata$Kfactor)
phenodatacut<-phenodata[!is.na(phenodata$Kfactor) & !is.na(phenodata$roundness2),]
phenodatacut<-phenodatacut[phenodatacut$Kfactor<1.75 ,]
#Add GPS coordinates
phenodatacut2<-merge(phenodatacut, Coord[,1:2], by.x="area", by.y="Area")
#and timepoint as a continuous variable
phenodatacut2$timepoint2<-as.numeric(as.character(phenodatacut2$timepoint))
str(phenodatacut2)

#A subset with only males and roundness1 females
#remove unsexed,
phenodatasex<-phenodatacut2[phenodatacut2$sex!="u",]
#female roundness could be obscuring trends because of different timing of spawning across pops
phenodatasexR1<-phenodatasex[phenodatasex$roundness2=="m" | phenodatasex$roundness2=="1", ]
summary(phenodatasexR1)

# Simple model latitude -> condition -------------
modconlat<-lm(Kfactor~Lat,data=phenodatacut2)
summary(modconlat)

# Simple model nr.2 Population -> condition -------------
modconarea<-lm(Kfactor~area+sex+area:sex,data=phenodatasexR1)
summary(modconarea)
Anova(modconarea, type=3)
plot(effect("area:sex",modconarea, x.var="area"), multiline=T) 
emmeans(modconarea, list(pairwise ~ area|sex), adjust = "tukey")

# Model for sex*population*timepoint->condition -----------------------------------------------------------
modcondall<-lm(Kfactor~area+timepoint2+sex
               +area:timepoint2 +timepoint2:sex +area:sex
               +area:timepoint2:sex
               , data=phenodatasexR1) #removed non-sign interaction

Anova(modcondall, type=3)
plot(effect("area:timepoint2:sex",modcondall, x.var="timepoint2"), multiline=T) 
#compare slopes to zero for each sex and area
m.lst.cond <- lstrends(modcondall, ~area|sex, var="timepoint2",adjust = "sidak") #get slopes for each sex
summary(m.lst.cond, infer=c(TRUE,TRUE),null=0) #compare each slope to zero 

#Model for the condition of males females and unsexed only in the two northern populations where they occur----------
phenoN<-phenodatacut[phenodatacut$area=="Ringstad" | phenodatacut$area=="Helligvaer",]
modcondN<-lm(Kfactor~area*sex, data=phenoN)
#Supplementary table S3
Anova(modcondN, type=3)
plot(effect("sex",modcondN))
#Table 2
emmeans(modcondN, list(pairwise ~ sex), adjust = "tukey")


#Condition figure --------------------------------------------------------------
#Figure 5.a
fuzi<-ggplot(phenodatacut2[phenodatacut2$roundness2=="1",], aes(x=area , group=area:timepoint, y=Kfactor,fill=area))+
  geom_hline(yintercept=0.8, alpha=0.25, linewidth=1)+
  geom_point(position=position_jitterdodge(), alpha=0.1, shape=21)+
  scale_color_manual("Population",values = popcols)+
  scale_fill_manual("Population",values = popcols)+
  ylim(0.3,1.6)+
  geom_boxplot(aes(fill=area), alpha=0.2, outlier.shape = NA)+
  labs(x="Population", y="Condition factor ",title=" a. Females (roundness1)")
#Figure 5.b
muzi<-ggplot(phenodatacut2[phenodatacut2$roundness2=="m",], aes(x=area , group=area:timepoint, y=Kfactor,fill=area))+
  geom_point(position=position_jitterdodge(), alpha=0.1, shape=21)+
  geom_hline(yintercept=0.8, alpha=0.25, linewidth=1)+
  scale_color_manual("Population",values = popcols)+
  scale_fill_manual("Population",values = popcols)+
  ylim(0.3,1.6)+
  geom_boxplot(aes(fill=area), alpha=0.2, outlier.shape = NA)+
  labs(x="Population", y="Condition factor ",title=" b. Males")
muzi
#Figure 5
ggarrange(fuzi, muzi, nrow=2, common.legend = T, legend = "right")

#Supplementary Figure S4.b
ggplot(phenodatacut2[phenodatacut2$roundness2=="1" | phenodatacut2$roundness2=="m" | phenodatacut2$roundness2=="u",], aes(x=area, group=area:sex,y=Kfactor,fill=area,color=sex))+
  geom_point(position=position_jitterdodge(), alpha=0.1, shape=21)+
  scale_color_manual("Sex",values=c("black","grey60","grey90"))+
  scale_fill_manual("Population",values = popcols)+
  geom_boxplot(aes(fill=area), alpha=0.2, outlier.shape = NA)+
  labs(x="Population", y="Condition factor ",title=" ")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DENSITY
#Import data with field observations of population census and more
OSRdata<-read.csv("OSRdata.csv")
OSRdata$area<-as.factor(OSRdata$area)
OSRdata$area<- factor(OSRdata$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
#Extract the count of males and females and format in long
longdensexf<-OSRdata[,c("area", "transectlength","Nf")]
longdensexm<-OSRdata[,c("area", "transectlength","Nm")]
colnames(longdensexf)<-c("area","transectlength","N")
colnames(longdensexm)<-c("area","transectlength","N")
longdensex<-rbind(cbind(longdensexf,sex="f"),cbind(longdensexm,sex="m"))
longdensex$sex<-as.factor(longdensex$sex)
#Add latitude info
longdensex<-merge(longdensex,Coord[,1:2], by.x="area", by.y="Area")
str(OSRdata)

#Simple model of the effect of latitude on density----------
modens<-lm(N/transectlength~Lat, data=longdensex)
summary(modens)

#Simple model nr. 2 of the effect of population total adult density (both sexes together)----------
modens<-lm(Density~area+sex, data=longdensex) #non-sign interaction
summary(modens)
Anova(modens, type=2)
emmeans(modens, list(pairwise ~ area), adjust = "tukey")



#Supplementary Figure S4.a-----------------------------------
ggplot(data=longdensex,aes(x=area,y=N/transectlength,group=area:sex))+
  geom_point( aes(x=area,y=N/transectlength,fill=area),position=position_jitterdodge(), shape=21, alpha=0.3)+
  geom_boxplot(aes(group=area:sex, fill=area, color=sex), alpha=0.5, outlier.shape=NA)+
  scale_color_manual("Sex",values=c("black","grey60"))+
  scale_fill_manual("Population",values = popcols)+
  scale_y_continuous(trans="log", breaks=c(0.01,0.1,0.5,1,2,5,10))+
  labs(y="Individual density (fish/m)",x="Date", title="")

ggplot(data=OSRdata,aes(x=area,y=Density))+
  #geom_point( aes(x=area,y=N/transectlength,fill=area),position=position_jitterdodge(), shape=21, alpha=0.3)+
  geom_boxplot(aes( fill=area), alpha=0.5, outlier.shape=NA)+
  scale_color_manual("Sex",values=c("black","grey60"))+
  scale_fill_manual("Population",values = popcols)+
  scale_y_continuous(trans="log", breaks=c(0.01,0.1,0.5,1,2,5,10))+
  labs(y="Individual density (fish/m)",x="Date", title="")





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 3 FEMALE SPAWNING

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Import phenotype dataset and edit
sexrpheno<-read.csv("sexrpheno.csv")
sexrpheno$area<-as.factor(sexrpheno$area)
sexrpheno$area<- factor(sexrpheno$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
sexrpheno$date<-as.Date(sexrpheno$date)
sexrpheno$timepoint<-as.factor(sexrpheno$timepoint)
sexrpheno$Catch<-sexrpheno$f1+sexrpheno$f2+sexrpheno$f3


#Figure 3.a ()
TimeF=ggplot(data=sexrpheno, aes(x=date, y=(f2+f3)/(f1+f2+f3), color=area))+
  geom_point(aes(size=Catch, fill=area), shape=23, alpha=0.3)+
  scale_x_date(date_labels = "%B",breaks="month")+
  theme(axis.text.x = element_text(angle=0, hjust = 1))+
  geom_line(stat="smooth",method = "lm", formula = y ~ x + I(x^2),
            size = 2.5,
            alpha = 0.6)+ 
  geom_boxplot(aes(group=area:timepoint, fill=area), width=2, alpha=0.5, colour="black", outlier.shape = NA)+
  scale_colour_manual(values = popcols)+
  scale_fill_manual("population",values = popcols)+
  xlim(as.Date("2022-04-15"),as.Date("2022-08-05"))+
  ylim(0,1)+
  labs(x="Month",y="Proportion of females ready to spawn ",title="a.")

#Figure 3.b
TempF=ggplot(data=sexrpheno, aes(x=temp2022, y= (f2+f3)/(f1+f2+f3), color= area))+
  geom_point(aes( fill=Population), shape=21, alpha=0.2)+
  geom_line(stat="smooth",method = "lm", formula = y ~ x + I(x^2),
            size = 1.5,
            alpha = 0.6)+ 
  scale_colour_manual(values = popcols)+
  scale_fill_manual(values = popcols)+
  ylim(0,1)+
  labs(x="Sea temperature (ST)",y="Prop. of females ready to spawn",title="b.")


###########
# Calculating degree days and standard broods left for each sampling date/population
# Degree days is calculated by integrating temperature curves over time
# Broods left is calculated by integrating development speed as a function of temperature over time

#Import estimated development speed over time in each location
DEV<-read.table("devtemp.txt", header=T)
DEV$Area<-as.factor(DEV$Location)
DEV$Experiment<-as.factor(DEV$Experiment)

#linear regression of log development time over temperature
moddev<-lm(data=DEV[DEV$Area=="Kristineberg",],log(Dev_time)~Temperature)

#Adjust day format to get day numbers, easier for integration
sexrpheno$daynr<-difftime(sexrpheno$date,as.Date("2022-01-01") ,units="days")

##I also need the list of start to end of primary productivity dates
primprod<-read.table("primaryprod.txt", header = T)
primprod$Start<-as.Date(primprod$Start)
primprod$End<-as.Date(primprod$End)
primprod$endnr<-difftime( primprod$End,as.Date("2022-01-01"),units="days")
primprod$stanr<-difftime( primprod$Start,as.Date("2022-01-01"),units="days")
primprod$area[7]<-"Northern Limit"

#add to sexrpheno the end date of primary productivity for each location
sexrpheno$endnr<-NA
pops<-levels(sexrpheno$area)
for (bi in 1:6) { sexrpheno[sexrpheno$area==pops[bi] ,]$endnr<-primprod[primprod$area==pops[bi] ,]$endnr }

# Calculate Broods left after each samplign occurrence
#integrate from daynr to endnr-> this gives broods left from sampling day to end of breeding window
sexrpheno$broodsleft<-NA
#for each line of the dataset calculated from the column area
for (bu in 1:length(sexrpheno$area)) {
  #from the day of the sampling recorded in lowers
  lowers<-as.numeric(sexrpheno$daynr[bu])
  #the number of boods left is calculated as the integral of the development model
  sexrpheno$broodsleft[bu]<-(integrate(function(x)
    1/(exp(coef(moddev)[1]+coef(moddev)[2]*
             (as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$Intercept)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin)*sin(2*pi/365*x)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos)*cos(2*pi/365*x)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin2)*sin(4*pi/365*x)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos2)*cos(4*pi/365*x)))),lower=lowers,upper=as.numeric(sexrpheno$endnr[bu])))$value
  
}

#Add the version for fixed reproduction window where endnumber is 243
sexrpheno$broodsleftF<-NA
#for each line of the dataset calculated from the column area
for (bu in 1:length(sexrpheno$area)) {
  #from the day of the sampling recorded in lowers
  lowers<-as.numeric(sexrpheno$daynr[bu])
  #the number of boods left is calculated as the integral of the temperature model x temperature curve
  sexrpheno$broodsleftF[bu]<-(integrate(function(x)
    1/(exp(coef(moddev)[1]+coef(moddev)[2]*
             (as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$Intercept)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin)*sin(2*pi/365*x)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos)*cos(2*pi/365*x)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin2)*sin(4*pi/365*x)+
                as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos2)*cos(4*pi/365*x)))),lower=lowers,upper=243))$value
  
}


#DegreeDays
#This should come to just integrating under the temperature curve from sampling to endnr
sexrpheno$DDleft<-NA
#for each line of the dataset calculated from the column area
for (bu in 1:length(sexrpheno$area)) {
  #from the day of te sampling recorded in lowers
  lowers<-as.numeric(sexrpheno$daynr[bu])
  #the number of degree days left is calculated as the integral of the  temperature curve
  sexrpheno$DDleft[bu]<-(integrate(function(x)
    as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$Intercept)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin)*sin(2*pi/365*x)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos)*cos(2*pi/365*x)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin2)*sin(4*pi/365*x)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos2)*cos(4*pi/365*x),lower=lowers,upper=as.numeric(sexrpheno$endnr[bu])))$value
  
}

#same for fixed window
sexrpheno$DDleftF<-NA
#for each line of the dataset calculated from the column area
for (bu in 1:length(sexrpheno$area)) {
  #from the day of te sampling recorded in lowers
  lowers<-as.numeric(sexrpheno$daynr[bu])
  #the number of degree days left is calculated as the integral of the  temperature curve
  sexrpheno$DDleftF[bu]<-(integrate(function(x)
    as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$Intercept)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin)*sin(2*pi/365*x)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos)*cos(2*pi/365*x)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$sin2)*sin(4*pi/365*x)+
      as.numeric(MODST[MODST$Area==as.character(sexrpheno$area[bu]),]$cos2)*cos(4*pi/365*x),lower=lowers,upper=243))$value
  
}

#Now I can plot female spawning as a function of degree days left and dtandard broods left.
#Figure 3.d
TempFb=ggplot(data=sexrpheno, aes(x=broodsleftF, y= (f2+f3)/(f1+f2+f3), color= area))+
  geom_point(aes(fill=Population), shape=21, alpha=0.2)+
  geom_line(stat="smooth",method = "lm", formula = y ~ x + I(x^2),
            size = 1.5,
            alpha = 0.6)+ 
  scale_colour_manual(values = popcols)+
  scale_fill_manual(values = popcols)+
  scale_x_reverse()+
  labs(x="Std. broods to winter",y="",title="d.")

#Same with degree days: Figure 3.c
TempFb2=ggplot(data=sexrpheno, aes(x=DDleftF, y= (f2+f3)/(f1+f2+f3), color= area))+
  geom_point(aes(fill=Population), shape=21, alpha=0.2)+
  geom_line(stat="smooth",method = "lm", formula = y ~ x + I(x^2),
            size = 1.5,
            alpha = 0.6)+ 
  scale_colour_manual(values = popcols)+
  scale_fill_manual(values = popcols)+
  scale_x_reverse()+
  labs(x="Degree days to winter",y="",title="c.")


#Figure 3
TempFs=ggarrange(TempF,TempFb2,TempFb, nrow = 1,legend = "none")
ggarrange(TimeF,TempFs, nrow = 2,heights=c(2,1), legend="right")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 4 Adult Sex Ratio  & Nest occupancy

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#using OSRdata already imported in 2. for density
str(OSRdata)
OSRdata$timepoint2<-as.factor(OSRdata$timepoint2)
OSRdata$date<-as.Date(OSRdata$date)


#Supplementary figure S5.a
ggplot(data=OSRdata, aes(x=date, y=Nm/(Nm+Nf), color=area))+
  geom_point(aes(size=Census, fill=area), shape=22, alpha=0.3)+
  scale_x_date(date_labels = "%B",breaks="month")+
  #geom_vline(xintercept=as.Date("2022-04-24"), linetype="dashed", size=1, alpha=0.4)+
  #geom_vline(xintercept=as.Date("2022-07-25"), linetype="dashed", size=1, alpha=0.4)+
  theme(axis.text.x = element_text(angle=0, hjust = 1), legend.position = "bottom",plot.margin = margin(0.1,1,0.1,0.1, "cm"))+
  #geom_smooth (aes(color=area),alpha=0.1, size=0)+
  #stat_smooth (aes(color=area),geom="line", alpha=0.2, size=3, span=1.5) +
  #stat_smooth(aes(y=(f2+f3)/(f1+f2+f3)),method = "lm", formula = y ~ x + I(x^2),  size = 0, alpha=0.1)+ 
  geom_line(stat="smooth",method = "lm", #formula = y ~ x + I(x^2),
            size = 1.5,
            alpha = 0.6)+ 
  geom_boxplot(aes(group=area:timepoint2, fill=area), width=2, alpha=0.7, colour="black", outlier.shape = NA)+
  scale_colour_manual("Population",values = popcols)+
  scale_fill_manual("Population",values = popcols)+
  scale_size_continuous(limits = c(0, 1000), 
                        breaks = c(50, 100, 500, 1000))+
  labs(x="Month",y="Proportion of males",title="a.")

#To plot nest info I need to import and format the nest data
nestdata<-read.csv("nestdata.csv")
nestdata$Nest_size<-as.factor(nestdata$Nest_size)
nestdata$Live_eggs<-as.factor(nestdata$Live_eggs)
nestdata$Eyed<-as.factor(nestdata$Eyed)
nestdata$Timepoint<-as.factor(nestdata$Timepoint)
nestdata$male<-as.factor(nestdata$male)
nestdata$area<- factor(nestdata$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
str(nestdata)

#create proportion data for live eggs in nests
livedata<-data.frame(table(nestdata$location, nestdata$Timepoint,nestdata$Live_eggs))
colnames(livedata)<-c("location","Timepoint","Live","count")
livedata2<-reshape(livedata,idvar=c("location","Timepoint"),timevar="Live",direction="wide")
colnames(livedata2)<-c("location","Timepoint","Eggs.na","Eggs.no","Eggs.yes")
#create a frame to put the the counts from my contingency table
liveframe<-unique(nestdata[,c(2:4)])
#merge
widenest<-merge(livedata2,liveframe,by=c("location","Timepoint"))
#proportion data for male presence in nests
maledata<-data.frame(table(nestdata$location, nestdata$Timepoint,nestdata$male))
colnames(maledata)<-c("location","Timepoint","male","count")
maledata2<-reshape(maledata,idvar=c("location","Timepoint"),timevar="male",direction="wide")
colnames(maledata2)<-c("location","Timepoint","male.na","male.no","male.yes")
widenest<-merge(widenest,maledata2,by=c("location","Timepoint"))
widenest$timepoint2<-widenest$Timepoint
levels(widenest$timepoint2)<-c("Mid eason","Late season")


#Figure S5.b, nest occupancy (eggs)
NestE<-ggplot(data=widenest, aes(x=timepoint2, y=Eggs.yes/(Eggs.yes+Eggs.no)))+
  geom_boxplot(aes(x=timepoint2,group=timepoint2:area,fill=area), alpha=0.2, show.legend = F)+
  geom_boxplot(aes(x=timepoint2), alpha=0.4, outlier.shape = NA)+
  #facet_wrap(~Timepoint)+
  geom_point(aes(size=(male.yes+male.no),group=timepoint2:area,color=area),position=position_dodge(width=0.75), alpha=0.5)+
  #geom_line(aes(x=area:Timepoint, group=location, color=area), size=1,alpha=0.3)+
  scale_color_manual("Population",values =popcols, guide="none")+
  scale_fill_manual("Population",values = popcols, guide="none")+
  #theme(axis.text.x = element_text(angle=45, hjust = 1))+
  scale_size_continuous("Sample size",
                        range=c(0.1,4),
                        limits = c(0, 20), 
                        breaks = c(5, 10, 20))+
  labs(y="Proportion of nests with eggs",x="Time period", title="b.")


#Figure S5.c, nest occupancy (males)
NestM<-ggplot(data=widenest, aes(x=timepoint2, y=male.yes/(male.yes+male.no)))+
  geom_boxplot(aes(x=timepoint2,group=Timepoint:area,fill=area), alpha=0.2, show.legend = F)+
  geom_boxplot(aes(x=timepoint2), alpha=0.4, outlier.shape = NA)+
  #facet_wrap(~Timepoint)+
  geom_point(aes(size=(male.yes+male.no),group=timepoint2:area,color=area),position=position_dodge(width=0.75), alpha=0.5, show.legend = F)+
  #geom_line(aes(x=area:Timepoint, group=location, color=area), size=1,alpha=0.3)+
  scale_color_manual(values = popcols)+
  scale_fill_manual(values = popcols)+
  #theme(axis.text.x = element_text(angle=45, hjust = 1))+
  scale_size_continuous("Sample size",
                        range=c(0.1,4),
                        limits = c(0, 20), 
                        breaks = c(5, 10, 20))+
  labs(y="Proportion of nests with male",x="Time period", title="c.")

#Figure S5
Nestiz<-ggarrange(NestE,NestM,ncol=1, nrow=2,common.legend = T, legend = "bottom")
ggarrange(Malz,Nestiz,ncol=2, widths=c(1.8,1))


#---------------------------------
#Statistical models

#ASR
#Adult sex ratio compared across times and populations

#Change some factor level names for plotting convenience later
OSRdata$timepointshort<-as.factor(OSRdata$timepoint2)
levels(OSRdata$timepointshort)<-c("Early", "Mid", "Late")
OSRdata$areashort<-OSRdata$area
levels(OSRdata$areashort)<-c('KBG', 'ARD', 'AUV', 'HIT', 'HEL', 'RIG')

modASR<-glm(data=OSRdata, family="binomial", cbind(Nm,Nf)~areashort*timepointshort)
Anova(modASR, type=3)

#make a table to plot a heatmap of pariwise comparison p.values
p.val.test<-pwpm(emmeans(modASR, list(pairwise ~ areashort:timepointshort), adjust = "tukey"),means = FALSE, flip = T,reverse = T) # p-values presented compactly in matrix form
p.val.test<-sub("[<>]", "", p.val.test)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
heatmap1<-melt(p.matx)
melt(p.matx)
#heatmap
ggplot(data=heatmap1,aes(X1, forcats::fct_rev(X2), fill = value)) + geom_tile() +
  geom_text(aes(label = value), size=3)+
  scale_fill_gradientn(colours = c("limegreen","white", "darkgrey"),
                       values = c(0,0.6,0.95,1), trans="log10")+
  labs(x="",y="")+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 45, vjust = -1, hjust=0))


#
# NEST OCCUPANCY (Latitude effect)
#Add Latitude info
widenest<-merge(widenest,Coord, by.x = "area",by.y="Area")

#egg occupancy
modegg<-glm(cbind(Eggs.yes,Eggs.no)~Lat+timepoint2,data=widenest, family="binomial") #interaction exclude because non-sign
Anova(modegg, type=2)

plot(effect("Lat",modegg,x.var="Lat"), multiline=TRUE)
summary(modegg)

#male occupancy
modmale<-glm(cbind(male.yes,male.no)~Lat+timepoint2,data=widenest, family="binomial") #interaction exclude because non-sign
Anova(modmale, type=2)

plot(effect("Lat",modmale,x.var="Lat"), multiline=TRUE)
plot(effect("timepoint2",modmale), multiline=TRUE)
summary(modmale)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 5 Size distributions 

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#To plot Figure 4 I need the phenotype data, already imported for condition above...

# ...I also need info from the fish for which scales were analyzed.
#Import dataset from scale anlysis
datasett_tilde <- readxl::read_excel("C:/Users/ivainm/Working_Document/DYNAMAR/Master projects/Tilde/datasett tilde.xlsx")
#Fix a population name
datasett_tilde[datasett_tilde$loc=="Austervoll",]$loc<-"Austevoll"
datasett_tilde$loc<- factor(datasett_tilde$loc, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
datasett_tilde$class<- factor(datasett_tilde$class, levels=c('L','S','U'))
levels(datasett_tilde$class)<-c("Large", "Small", "Unsexed")
datasett_tilde$ID<-as.factor(datasett_tilde$ID)
datasett2<-datasett_tilde[datasett_tilde$sclerites!="focus",]
datasett2$sclerites<-as.numeric(datasett2$sclerites)
colnames(datasett2)<-c("ID","sclerites","scale1","scale2","scale3","loc","class","length","sex")
#I want to remove five outliers for which the Class identity is not clear (body size too close to limit value)
datasett2<-datasett2[!(datasett2$loc=="Kristineberg" & datasett2$length>48 & datasett2$length<51),]
datasett2<-datasett2[!(datasett2$loc=="Arendal" & datasett2$length>43 & datasett2$length<46),]
datasett2<-datasett2[!(datasett2$loc=="Austevoll" & datasett2$length>50 & datasett2$length<52),]
#calculate average intersclerite distance for scale analysis
datasett3<-data.frame(datasett2 %>% dplyr::rowwise() %>% dplyr::mutate(scaleav2 = mean(c(scale1, scale2, scale3), na.rm=T)) )
#extract info of the fish used for scale analysis
fishinfo<-unique(datasett3[,c("ID","loc","class","length","sex")]) 
#fix some variables
fishinfo$area<-fishinfo$loc
fishinfo$sex<-as.factor(fishinfo$sex)


#Figure 4 is focusing on early season (timepoint=1) size distributions
ggplot(phenodata[phenodata$timepoint==1,],aes(x=area,y=length,group=sex:area))+
  geom_violin(aes(group=sex:area,fill=sex), position=position_dodge(width=0.6), notch=T, width=0.7, outlier.shape = NA)+
  geom_point(aes(group=sex:area,fill=sex),position=position_jitterdodge(dodge.width=0.6), alpha=0.2, color="white")+
  #facet_grid(rows=vars(timepoint))+
  theme(text = element_text(size = 10)) +
  scale_fill_manual("Sex",values = c("f" = "goldenrod1",
                                     "m"="mediumpurple1",
                                     "u"="darkolivegreen1")) +
  scale_color_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  geom_point(data=fishinfo,aes(color=class), size=2, shape=4,position=position_dodge(width=0.6), stroke=1)+
  theme(axis.text=element_text(size = 12),axis.title=element_text(size = 11), text = element_text(size = 15),
        panel.background = element_rect(fill = 'grey75', colour = 'white')) +
  labs(x="Population",y=" Length (mm)",title=" ")


#Visual inspection of Figure 4 is used to determine for each population the limit between the two age classes.
#Kristineberg: 1 age class
#Arendal: 1 age class
#Austevoll and Hitra: Likely 2 age classes but can't be visually determined on the Figure
#Helligvaer: 52mm 
#Ringstad: 46mm


#I create a class variable for phenodata based on the criterion above
phenodata$class<-"1yo"
phenodata[phenodata$sex=="u",]$class<-"unsexed" #adults with no secondary sexual characters are set apart
phenodata[phenodata$area=="Helligvaer" & phenodata$timepoint=="1" & phenodata$sex!="u" & phenodata$length>=52,]$class<-"2yo" #T1 Helligvaer individuals >52mm are 2yo
phenodata[phenodata$area=="Ringstad" & phenodata$timepoint=="1" & phenodata$sex!="u" & phenodata$length>=46,]$class<-"2yo" #T1 Ringstad individuals >46mm are 2yo

table(phenodata[phenodata$timepoint=="1",]$area,phenodata[phenodata$timepoint=="1",]$class)

#Model comparing size of southern fish vs northern fish of the two classes
phenodataT1<-phenodata[phenodata$timepoint=="1",] #subset timepoint1
phenodataT1$popclass<-as.character(phenodataT1$area) #dummy variable area to avoid aliasing of coefficients
phenodataT1[phenodataT1$area=="Helligvaer" & phenodataT1$class=="2yo",]$popclass<-"Helligvaer_2" # variable contains population info AND class
phenodataT1[phenodataT1$area=="Helligvaer" & phenodataT1$class=="unsexed",]$popclass<-"Helligvaer_U"
phenodataT1[phenodataT1$area=="Ringstad" & phenodataT1$class=="2yo",]$popclass<-"Ringstad_2" # variable contains population info AND class
phenodataT1[phenodataT1$area=="Ringstad" & phenodataT1$class=="unsexed",]$popclass<-"Ringstad_U"

modclasssize<-lm(length~popclass,data=phenodataT1)

plot(effect("popclass",modclasssize), multiline=TRUE)
#make a table to plot a heatmap of pariwise comparison p.values
p.val.test<-pwpm(emmeans(modclasssize, list(pairwise ~ popclass), adjust = "tukey"),means = FALSE, flip = T,reverse = T) # p-values presented compactly in matrix form
p.val.test<-sub("[<>]", "", p.val.test)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
heatmap1<-melt(p.matx)
melt(p.matx)
#heatmap
ggplot(data=heatmap1,aes(X1, forcats::fct_rev(X2), fill = value)) + geom_tile() +
  geom_text(aes(label = value), size=3)+
  scale_fill_gradientn(colours = c("limegreen","white", "darkgrey"),
                       values = c(0,0.6,0.95,1), trans="log10")+
  labs(x="",y="")+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 45, vjust = -1, hjust=0))


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Side quest from Sebastian:
#Check in Northern populations if there are nest holders and round females in the Small class
str(phenodata)
#Ringstad females
ggplot(phenodata[phenodata$area=="Ringstad" & phenodata$sex=="f",],aes(x=length,group=sex))+
 geom_histogram(aes(fill=sex),alpha=0.2)+
 geom_histogram(data=phenodata[phenodata$area=="Ringstad" &
                              phenodata$roundness2=="3",],alpha=0.25, color="orange")+
  geom_histogram(data=phenodata[phenodata$area=="Ringstad" &
                                  phenodata$roundness2=="2",],alpha=0.15, color="green")+
  facet_grid(rows=vars(timepoint))+
  theme(text = element_text(size = 10)) +
  scale_fill_manual("Sex",values = c("f" = "goldenrod1",
                                     "m"="skyblue",
                                     "u"="darkolivegreen1")) +
  theme(axis.text=element_text(size = 12),axis.title=element_text(size = 11), text = element_text(size = 15),
        panel.background = element_rect(fill = 'white', colour = 'grey')) +
  labs(x="Size",y="Count",title="Females in Ringstad")

#Helligvaer females
ggplot(phenodata[phenodata$area=="Helligvaer" & phenodata$sex=="f",],aes(x=length,group=sex))+
  geom_histogram(aes(fill=sex),alpha=0.2)+
  geom_histogram(data=phenodata[phenodata$area=="Helligvaer" &
                                  phenodata$roundness2=="3",],alpha=0.25, color="orange")+
  geom_histogram(data=phenodata[phenodata$area=="Helligvaer" &
                                  phenodata$roundness2=="2",],alpha=0.15, color="green")+
  facet_grid(rows=vars(timepoint))+
  theme(text = element_text(size = 10)) +
  scale_fill_manual("Sex",values = c("f" = "goldenrod1",
                                     "m"="skyblue",
                                     "u"="darkolivegreen1")) +
  theme(axis.text=element_text(size = 12),axis.title=element_text(size = 11), text = element_text(size = 15),
        panel.background = element_rect(fill = 'white', colour = 'grey')) +
  labs(x="Size",y="Count",title="Females in Helligvaer")

#For males I need a different dataset
str(nestdata)
nestdata$timepoint<-nestdata$Timepoint
#Ringstad males
ggplot(phenodata[phenodata$area=="Ringstad" & phenodata$sex=="m",],aes(x=length))+
  geom_histogram(alpha=0.3, fill="skyblue")+
 geom_histogram(data=nestdata[nestdata$area=="Ringstad",],alpha=0.25, color="orange")+
  facet_grid(rows=vars(timepoint))+
  theme(text = element_text(size = 10)) +
  theme(axis.text=element_text(size = 12),axis.title=element_text(size = 11), text = element_text(size = 15),
        panel.background = element_rect(fill = 'white', colour = 'grey')) +
  labs(x="Size",y="Count",title="Males in Ringstad")

#Helligvaer males
ggplot(phenodata[phenodata$area=="Helligvaer" & phenodata$sex=="m",],aes(x=length))+
  geom_histogram(alpha=0.3, fill="skyblue")+
  geom_histogram(data=nestdata[nestdata$area=="Helligvaer",],alpha=0.25, color="orange")+
  facet_grid(rows=vars(timepoint))+
  theme(text = element_text(size = 10)) +
  theme(axis.text=element_text(size = 12),axis.title=element_text(size = 11), text = element_text(size = 15),
        panel.background = element_rect(fill = 'white', colour = 'grey')) +
  labs(x="Size",y="Count",title="Males in Helligvaer")

##
#Second side quest, what are the condition trends over time if we separate individuals by size class?!
#That's tricky to do in practice because the size classes are not defined the same at each time point, and actually not easy to define at T3 at all

#Model for the interaction effect of time and body size on condition in the northern populations
phenoN<-phenodatacut[phenodatacut$area=="Ringstad" | phenodatacut$area=="Helligvaer",]
str(phenoN)
modcondsize<-lm(Kfactor~length*timepoint*sex, data=phenoN)
#Supplementary table S3
Anova(modcondsize, type=3)
plot(effect("length:timepoint:sex",modcondsize,x.var="timepoint" ), multiline=T)
#Table 2
emmeans(modcondN, list(pairwise ~ sex), adjust = "tukey")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 6 Scale analysis 

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Use the datasett3 imported just above in section 5. And fishinfo, also created in 5.

#plot growth profiles from the scales. Supplementary fiugre S8.a
Growthz=ggplot(datasett3,aes(x=sclerites,y=scaleav2,group=ID,color=class))+
  geom_vline(xintercept=03, linetype="dotted", size=0.8, alpha=0.4)+
  #geom_vline(xintercept=12, linetype="dotted", size=0.8, alpha=0.4)+
  geom_vline(xintercept=21, linetype="dotted", size=0.8, alpha=0.4)+
  geom_vline(xintercept=30, linetype="dotted", size=0.8, alpha=0.4)+
  geom_point(alpha=0.2)+
  geom_line(linewidth=1, alpha=0.4)+
  scale_color_manual("Size class",values =c("brown2","navy","olivedrab4"))+
  # geom_line(aes(x=sclerites,y=scaleav2,group=class,color=class), stat="smooth",method = "loess", size = 2.5,alpha = 0.6)+ 
  facet_wrap(~loc)+
  labs(x="Sclerite number", y="Mean sclerite length",title="a.")



#We want to extract some info from these growth profiles to compare them: number of sclerites, and mean intersclerite distance

#Number of slcerite for each fish -Format data
sclemax<-plyr::dlply(datasett3, "ID", function(df)  max(df$sclerites))
maxsclerite<-plyr::ldply(sclemax)
colnames(maxsclerite)<-c("ID","max_sclerite")
maxsclerite<-merge(maxsclerite,fishinfo, by="ID")
str(maxsclerite)

#mean intersclerite distance for each fish over entire growth period-Format data
meangr<-plyr::dlply(datasett3, "ID", function(df)  mean(df$scaleav2))
meangrowth<-plyr::ldply(meangr)
colnames(meangrowth)<-c("ID","mean_growth")
maxsclerite<-merge(maxsclerite,meangrowth, by="ID")

#mean intersclerite distance for each fish Early-Late period separated
datasettscle_early<-datasett3[datasett3$sclerites<=21,] 
datasettscle_late<-datasett3[datasett3$sclerites>=21,] 
#now we also want to calculate the mean growth for each fish
means1<-plyr::dlply(datasettscle_early, "ID", function(df)  mean(df$scaleav2))
meangrowth1<-plyr::ldply(means1)
colnames(meangrowth1)<-c("ID","Meangrowth")
means2<-plyr::dlply(datasettscle_late, "ID", function(df)  mean(df$scaleav2))
meangrowth2<-plyr::ldply(means2)
colnames(meangrowth2)<-c("ID","Meangrowth")
#merge
LMdata1<-merge(fishinfo,meangrowth1, by="ID")
LMdata2<-merge(fishinfo,meangrowth2, by="ID")
#merge periods
LMdataAll<-rbind(cbind(LMdata1,Period="Early growth"),
                 cbind(LMdata2,Period="Late growth"))
LMdataAll$Period<- factor(LMdataAll$Period, levels=c('Early growth','Late growth'))
LMdataAll$sex<-as.factor(LMdataAll$sex)
summary(LMdataAll)




#Show the relationship between scale length and fish body length.  Supp. figure S8.b
Fishsize<-ggplot(maxsclerite, aes(x=mean_growth*max_sclerite/3,y=length))+geom_point(aes( fill=class),size=3,shape=21, alpha=0.8)+
  #facet_wrap(~loc)+
  geom_line(aes(), stat="smooth",method = "lm", size = 2.5,alpha = 0.4)+
  ggpubr::stat_regline_equation(label.y = 60, label.x=1.5, aes(label = ..rr.label..))+
  scale_color_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  scale_fill_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  labs(x="Scale length (mm) ",y="Body length (mm)", title="b.")
Fishsize

#Supplementary fiure S8
ggpubr::ggarrange( Growthz,Fishsize, ncol=2,nrow = 1, widths=c(1.7,1),legend="right",common.legend = T)


#Figure 6.a. Numbers of slcerites per fish class and population
Scleritz <-ggplot(maxsclerite,aes(x=loc,y=max_sclerite,group=loc:class,color=class, fill=class))+
  geom_boxplot(alpha=0.2, color="black", outlier.shape=NA)+
  facet_wrap(~class) +
  geom_point(alpha=0.5,size=2, position = position_jitterdodge())+
  scale_color_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  scale_fill_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  labs(x="", y="Sclerite number", title="a.")+
  #ylim(0.1,0.5)+
  theme(axis.text.x = element_text(angle=25, hjust = 1)) 

#Figure 6.b. mean intersclerite distance per fish class and population
Meangrowth <-ggplot(maxsclerite,aes(x=loc,y=mean_growth,group=loc:class,color=class, fill=class))+
  geom_boxplot(alpha=0.2, color="black", outlier.shape=NA)+
  facet_wrap(~class) +
  geom_point(alpha=0.5,size=2, position = position_jitterdodge())+
  scale_color_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  scale_fill_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  labs(x="", y="Mean intersclerite distance", title="b.")+
  #ylim(0.1,0.5)+
  theme(axis.text.x = element_text(angle=25, hjust = 1)) 

#Figure 6
ggarrange(Scleritz,Meangrowth, nrow=2, ncol=1, common.legend = T, legend="bottom")


#Alternative figure 6.x with the eraly and late period separated for growth
MeangrowthE <-ggplot(LMdataAll[LMdataAll$Period=="Early growth",],aes(x=loc,y=Meangrowth,color=class, fill=class))+
  geom_boxplot(alpha=0.2, color="black", outlier.shape=NA)+
  facet_wrap(~class) +
  geom_point(alpha=0.5,size=2, position = position_jitterdodge())+
  scale_color_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  scale_fill_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  labs(x="", y="Mean intersclerite distance", title="a.")+
  #ylim(0.1,0.5)+
  theme(axis.text.x = element_text(angle=25, hjust = 1))
MeangrowthE
MeangrowthL <-ggplot(LMdataAll[LMdataAll$Period=="Late growth",],aes(x=loc,y=Meangrowth,color=class, fill=class))+
  geom_boxplot(alpha=0.2, color="black", outlier.shape=NA)+
  facet_wrap(~class) +
  geom_point(alpha=0.5,size=2, position = position_jitterdodge())+
  scale_color_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  scale_fill_manual("Size class",values = c("brown2","navy","olivedrab4"))+
  labs(x="", y="Mean intersclerite distance", title="b.")+
  #ylim(0.1,0.5)+
  theme(axis.text.x = element_text(angle=25, hjust = 1))
MeangrowthL
summary(LMdataAll)



#--------------------------
#Statistical model to compare class and populations for sclerite number and intersclerite distance.

#1. exclude unsexed and do a class by population model. Sclerites
scleriteLS<-maxsclerite[maxsclerite$class!="Unsexed",]
hist (scleriteLS$max_sclerite)
modclassLS<-lm(max_sclerite~area+class,data=scleriteLS)#exclude non sign interaction
summary(modclassLS)
Anova(modclassLS, type=2)
#There is no interaction so I just look at pairwise comparison among pops
emmeans(modclassLS, list(pairwise ~ class), adjust = "tukey")
pairwizeSCL<-emmeans(modclassLS, list(pairwise ~ area), adjust = "tukey")

#make a table to plot a heatmap of pariwise comparison p.values
p.val.test<-pwpm(pairwizeSCL,means = FALSE, flip = T,reverse = T) # p-values presented compactly in matrix form
p.val.test<-sub("[<>]", "", p.val.test)
p.matx<-matrix(as.numeric((p.val.test)),nrow = length(p.val.test[,1]),ncol = length(p.val.test[,1])) #if your factor has 5 levels ncol and nrow=5
rownames(p.matx) <- colnames(p.matx) <-colnames(p.val.test)
p.matx[upper.tri(p.matx, diag=FALSE)] <- NA
heatmap1<-melt(p.matx)
melt(p.matx)
#heatmap
ggplot(data=heatmap1,aes(X1, forcats::fct_rev(X2), fill = value)) + geom_tile() +
  geom_text(aes(label = value), size=3)+
  scale_fill_gradientn(colours = c("limegreen","white", "darkgrey"),
                       values = c(0,0.6,0.95,1), trans="log10")+
  labs(x="",y="")+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 45, vjust = -1, hjust=0))



#1.2 exclude unsexed and do a class by population model. Intersclerite distance
modclassLS2<-lm(mean_growth~area*class,data=scleriteLS)
summary(modclassLS2)
Anova(modclassLS2, type=3)
emmeans(modclassLS2, list(pairwise ~ area), adjust = "tukey")
pairwizeDIS<-emmeans(modclassLS2, list(pairwise ~ area|class), adjust = "tukey", digits=9)
pairwizeDIS

#make a table to plot a heatmap of pariwise comparison p.values
p.val.test1<-pwpm(pairwizeDIS,means = FALSE, flip = T,reverse = T)[[1]] # p-values presented compactly in matrix form for Large
p.val.test2<-pwpm(pairwizeDIS,means = FALSE, flip = T,reverse = T)[[2]] # p-values presented compactly in matrix form for Small
p.val.test2<-sub("[<>]", "", p.val.test2)
p.matx1<-matrix(as.numeric((p.val.test1)),nrow = length(p.val.test1[,1]),ncol = length(p.val.test1[,1])) #if your factor has 5 levels ncol and nrow=5
p.matx2<-matrix(as.numeric((p.val.test2)),nrow = length(p.val.test2[,1]),ncol = length(p.val.test2[,1])) 
rownames(p.matx1) <- colnames(p.matx1) <-colnames(p.val.test1)
rownames(p.matx2) <- colnames(p.matx2) <-colnames(p.val.test2)
p.matx1[upper.tri(p.matx1, diag=FALSE)] <- NA
heatmap1<-melt(p.matx1)
p.matx2[upper.tri(p.matx2, diag=FALSE)] <- NA
heatmap2<-melt(p.matx2)
heatmap<-rbind(cbind(heatmap1,Class="Large"),cbind(heatmap2,Class="Small"))

#heatmap
ggplot(data=heatmap,aes(X1, forcats::fct_rev(X2), fill = value)) + geom_tile() +
  facet_wrap(~factor(Class))+
  geom_text(aes(label = value), size=3)+
  scale_fill_gradientn(colours = c("limegreen","white", "darkgrey"),
                       values = c(0,0.6,0.95,1), trans="log10")+
  labs(x="",y="")+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 45, vjust = -1, hjust=0))




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 6 Popgen plots

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Let's start by gathering the environmental varialbes
#Here I pick up temperature variables. Broods 1-3 is obsolete and should be ignored
Env<-read.table("EnvT.txt", header=T)
Env<-Env[1:6,c(1:6,10)]

#This is the updated Broods variable, plus the degreedays variables
Env2<-read.csv("env_info.txt", header = T)
Env2<-Env2[1:6,]

#This is coastal (dispersal?) distance
Dist<-read.table("Dist.txt", header=T)
Dist<-Dist[1:6,]

#gps coordinates
Coord<-read.table("Coord.txt", header=T)
Coord<-Coord[1:6,]

#I also need to standardize the variables.
STDEnv<-data.frame(STD_Month_max=scale(Env[,1]),
                   STD_Month_min=scale(Env[,2]),
                   STD_Day_max=scale(Env[,3]),
                   STD_Day_min=scale(Env[,4]),
                   STD_Month_Amp=scale(Env[,5]),
                   STD_Day_Amp=scale(Env[,6]),
                   STD_Broods=scale(Env2[,2]),
                   STD_Ddays=scale(Env2[,3]),
                   STD_Ddaysfixed=scale(Env2[,4]),
                   STD_Dist=scale(Dist[,2]),
                   STD_Lat=scale(Coord[,2]),
                   STD_Long=scale(Coord[,3]),
                   Area=Env2[,1]) 
STDEnv$Area<- factor(STDEnv$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))

##I want to create population-level distance matrices. Also for Fst 
DistM<-dist(STDEnv[,10])#Coastal distance matrix
DdaysM<-dist(STDEnv[,8])  # Degreedays matrix
DdaysfixedM<-dist(STDEnv[,9]) # Degreedays fixed matrix
TmaxM<-dist( STDEnv[,1])  #Monthly temp max
TminM<-dist( STDEnv[,2])  #etc...
TAmpM<-dist(STDEnv[,5])

#Fst data from radiator package
FstM<-matrix(c(NA,0.001649540,0.007364404,0.009187501,0.010578569,0.013108584,
               NA,NA,0.006196770,0.008337682,0.009603892,0.012435670,
               NA,NA,NA,0.001921923,0.002773006,0.006501010,
               NA,NA,NA,NA,0.000537745,0.002720528,
               NA,NA,NA,NA,NA,0.002203149,
               NA,NA,NA,NA,NA,NA),
             ncol=6, nrow=6)

FstM<-as.dist(FstM)
FstM[1:15]
#data frame and adding type of comparison: within or across southern and northern groups
MatC<-data.frame(Fst=FstM[1:15],Distance=DistM[1:15],
                 Ddays=DdaysM[1:15],Ddays_fixed=DdaysfixedM[1:15],
                 Tmax=TmaxM[1:15],Tmin=TminM[1:15],
                 Tamp=TAmpM[1:15], Comp=c("Within","Across","Across","Across","Across",
                                          "Across","Across","Across","Across",
                                          "Within","Within","Within",
                                          "Within","Within",
                                          "Within"))
MatC$Comp<-as.factor(MatC$Comp)


# Color scale for plotting
Compscale<-scale_fill_manual("Comparison",values=c("orange1","dodgerblue"))
Compshape<-scale_shape_manual("Comparison",values=c(24,23))
#---------------------
#Plot of Figure 8
A=ggplot(MatC,aes(x=Distance,y=Fst, fill=Comp, shape=Comp))+Compscale+Compshape+
  geom_point(size=3)+labs(title="a.", x="Coastal distance")#+ theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm")) 
B=ggplot(MatC,aes(x=Ddays,y=Fst, fill=Comp, shape=Comp))+Compscale+Compshape+
  geom_point(size=2)+labs(title="d.", x="DD distance")
C=ggplot(MatC,aes(x=Ddays_fixed,y=Fst, fill=Comp, shape=Comp))+Compscale+Compshape+
  geom_point(size=2)+labs(title="b.", y="",x="DDF distance")
D=ggplot(MatC,aes(x=Tmax,y=Fst, fill=Comp, shape=Comp))+Compscale+Compshape+
  geom_point(size=2)+labs(title="c.", y="",x="Max. T. distance")
E=ggplot(MatC,aes(x=Tmin,y=Fst, fill=Comp, shape=Comp))+Compscale+Compshape+
  geom_point(size=2)+labs(title="d.",y="",x="Min. T. distance")#+ theme(plot.margin = margin(0.1,0.8,0.1,0.8, "cm")) 
G=ggplot(MatC,aes(x=Tamp,y=Fst, fill=Comp, shape=Comp))+Compscale+Compshape+
  geom_point(size=2)+labs(title="e.",y="",x="T. amp. distance")#+ theme(plot.margin = margin(0.1,0.8,0.1,0.8, "cm")) 

smallz3<-ggpubr::ggarrange(C,D,E,G, common.legend=T,legend = "right")
ggpubr::ggarrange(A,smallz3, common.legend = T, legend="none", widths=c(0.8,1))

#------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#use genepop package to calculate global Fst and perform allele frequency tests

#import ped file generated from PLINK, from the neutral SNPs

#vcfR package
myvcf<-read.vcfR("plink90HWLD05mafR.vcf")
mygenind<-vcfR2genind(myvcf) #may need to edit "sep" here?
#radiator package
radiator::genomic_converter(mygenind,
                              strata = "pop_conv2.txt",
                              output = "genepop",
                              filename = "plink90HW_R_maf.txt",
                              verbose = TRUE)

#genepop
basic_info("plink90HW_R_maf.txt", outputFile = "plink90HW_R_maf.info", verbose = interactive())

#allele frequency differences
test_diff(
  "plink90HW_R_maf.txt",
  genic = F,
  pairs = FALSE, #if F, single global test
  outputFile = "testdiff2",
  settingsFile = "",
  dememorization = 10000,
  batches = 100,
  iterations = 5000,
  verbose = interactive()
)

#Fst from genepop
Fst(
  "plink90HW_R_maf.txt",
  sizes = FALSE,
  pairs = FALSE,
  outputFile = "fst1",
  dataType = "Diploid",
  verbose = interactive()
)

genfst<-read.table("genepop_fst.txt", header=T)
str(genfst)
hist(genfst$Fst, breaks=25)
n<-length(genfst$Fst)
m<-mean(genfst$Fst)
s<-sd(genfst$Fst)
#mean Fst and confidence interval calculation
margin <- qt(0.975,df=n-1)*s/sqrt(n)
c(m-margin,m,m+margin)

#------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Admixture analysis
#Admixture
#import ped file generated from PLINK, from the neutral SNPs
ped2geno("plink90HW_R_maf.ped","plink90HW_R_maf.geno")
#cross-entropy plot
obj.snmf1 = snmf("plink90HW_R_maf.geno", K = 1:8, ploidy = 2, entropy = T,
                 alpha = 100, project = "new")
plot(obj.snmf1, col = "blue4", cex = 1.4, pch = 19)

#calculate admixture for chosen number of clusters
obj.snmf <- snmf("plink90HW_R_maf.geno", K = 2, repetitions=1, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 2, run=1)

#prepare data for plotting
admix<-data.frame(qmatrix)
popmap<-read.delim("pop_map_final.txt", header=T, sep='')
outliers<-read.table("out_all.txt")
admix<-cbind(popmap[!(popmap$ID %in% outliers$V1),][,1:6],admix)
admix$Area<-factor(admix$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
admix$Timepoint<-as.factor(admix$Timepoint) 
admix<-admix[order(admix$Area),]
str(admix)
table(popmap[!(popmap$ID %in% outliers$V1),]$Area)
summary()
#plot overall admix
admix_barplot(
  admix,
  K = 7:ncol(admix),
  individuals = 1,
  sortkey = "V1",
  grouping = "Area",
  palette = "viridis",
  names = F,
  xlab = "Individuals",
  ylab = "Ancestry",
  main = "",
  noclip = FALSE
)

###
#Supplementary plot PCA
#import output from plink
pca90_maf<-read.delim("pca90_maf.eigenvec", header=F, sep= ' ')
#population infos
popmap<-read.delim("pop_map_final.txt", header=T, sep='')
#append population info to eigenvectors
pca90_maf<-cbind(pca90_maf,popmap[popmap$ID %in% pca90_maf$V1,])
pca90_maf$Area<-as.factor(pca90_maf$Area)
pca90_maf$Area<- factor(pca90_maf$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))

#PCA 
ggplot(data=pca90_maf,aes(x=V3,y=V4))+
  geom_point(aes(fill=Area), shape=21, size=3, stroke=1, alpha=0.8)+
  stat_ellipse(level = 0.95, aes(color=Area), size=1.5, alpha=0.4)+
  geom_mark
#geom_text(aes(label=ID, color=Area), nudge_y = 0.02,nudge_x = 0.01)+
scale_color_manual("Population",values = c("red3",
                                           "brown1",
                                           "darkorchid4",
                                           "mediumpurple1",
                                           "blue3",
                                           "steelblue1"))+
  scale_fill_manual("Population",values = c("red3",
                                            "brown1",
                                            "darkorchid4",
                                            "mediumpurple1",
                                            "blue3",
                                            "steelblue1"))+
  labs(x="PC 1",y="PC 2",
       title="")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)

# 6 Final discussion figure

#__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)__(=)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#here I want to compare the populations in all their dimensions
#I can start from the EnvD dataset that has 1 line per pop and some environmental variables + Latitude
EnvD
Popcomp<-EnvD[,c("Area","Broodsfixed","Lat")]

# Add mean % ancestry group 1 from admixture, from the admix dataframe, see previous section
Ancestry<-plyr::dlply(admix, "Area", function(df)  mean(df$V2))
Ancestry<-plyr::ldply(Ancestry)
Ancestry$V1<-Ancestry$V1*100
colnames(Ancestry)<-c("Area","Ancestry")
Popcomp<-merge(Popcomp,Ancestry, by="Area")

#Add early season ASR
ASR<-plyr::dlply(OSRdata[OSRdata$timepoint==1,], "area", function(df)  mean(df$ASR))
ASR<-plyr::ldply(ASR)
ASR$V1<-100*ASR$V1
colnames(ASR)<-c("Area","ASR")
Popcomp<-merge(Popcomp,ASR, by="Area")

#Add % of 2yo at T1
twoyo<-cbind(Area=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'),tyo=c(0,0,0,0,27,70))
Popcomp<-merge(Popcomp,twoyo, by="Area")
#Add % of unsexed at T1
unsex<-cbind(Area=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'),unsexed=c(0,0,0,0,1,22))
Popcomp<-merge(Popcomp,unsex, by="Area")
#Add %sapwning female at last sampling
spf<-plyr::dlply(sexrpheno[sexrpheno$timepoint=="Timepoint 3",], "area", function(df)  mean((df$f2+df$f3)/(df$f1+df$f2+df$f3)))
spf<-plyr::ldply(spf)
spf$V1<-100*spf$V1
colnames(spf)<-c("Area","spf")
Popcomp<-merge(Popcomp,spf, by="Area")
#Add condition over time for males and females (slope? sign?)
condslopef<-cbind(Area=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'),slopf=c(0,0,0.001742,0.005860,0,0.004103))
condslopem<-cbind(Area=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'),slopm=c(-0.002911,0,0.00227,0.004711,0,0.001822))
Popcomp<-merge(Popcomp,condslopef, by="Area")
Popcomp<-merge(Popcomp,condslopem, by="Area")

#Condition sexually mature adults
condmean<-plyr::dlply(phenodatasexR1, "area", function(df)  mean(df$Kfactor))
condmean<-plyr::ldply(condmean)
colnames(condmean)<-c("Area","condition")
Popcomp<-merge(Popcomp,condmean, by="Area")

#Total fish density
dens<-plyr::dlply(OSRdata, "area", function(df)  mean(df$Density))
dens<-plyr::ldply(dens)
colnames(dens)<-c("Area","density")
Popcomp<-merge(Popcomp,dens, by="Area")

#Sclerite number (mean)
str(maxsclerite)
scle<-plyr::dlply(maxsclerite, "loc", function(df)  mean(df$max_sclerite))
scle<-plyr::ldply(scle)
colnames(scle)<-c("Area","sclerites")
Popcomp<-merge(Popcomp,scle, by="Area")
#Growth speed (mean intersclerite distance)
interscle<-plyr::dlply(maxsclerite, "loc", function(df)  mean(df$mean_growth))
interscle<-plyr::ldply(interscle)
interscle$V1<-interscle$V1*1000
colnames(interscle)<-c("Area","interscelrite")
Popcomp<-merge(Popcomp,interscle, by="Area")

#finalise data
Popcomp$Area<- factor(Popcomp$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
Popcomp$unsexed<-as.numeric(Popcomp$unsexed)
Popcomp$tyo<-as.numeric(Popcomp$tyo)
Popcomp$slopf<-as.numeric(Popcomp$slopf)
Popcomp$slopm<-as.numeric(Popcomp$slopm)
Popcomp
str(Popcomp)
#plotting options
#Format data for heatmap
#scale values --OBSOLETE
Popcomp$BroodsfixedS<-scale(Popcomp$Broodsfixed)
Popcomp$AncestryS<-scale(Popcomp$Ancestry)
Popcomp$ASRS<--scale(Popcomp$ASR)
Popcomp$LatS<--scale(Popcomp$Lat)
Popcomp$tyoS<--scale(Popcomp$tyo)
Popcomp$unsexedS<--scale(Popcomp$unsexed)
Popcomp$spfS<--scale(Popcomp$spf)
Popcomp$slopfS<--scale(Popcomp$slopf)
Popcomp$slopmS<--scale(Popcomp$slopm)
Popcomp$densityS<-scale(Popcomp$density)
Popcomp$scleritesS<-scale(Popcomp$sclerites)
Popcomp$interscelriteS<--scale(Popcomp$interscelrite)
Popcomp$conditionS<-scale(Popcomp$condition)

#or scale on a 0-1 interval?
Popcomp$BroodsfixedS2<-(Popcomp$Broodsfixed-min(Popcomp$Broodsfixed))/(max(Popcomp$Broodsfixed)-min(Popcomp$Broodsfixed))
Popcomp$AncestryS2<-(Popcomp$Ancestry-min(Popcomp$Ancestry))/(max(Popcomp$Ancestry)-min(Popcomp$Ancestry))
Popcomp$ASRS2<-1-(Popcomp$ASR-min(Popcomp$ASR))/(max(Popcomp$ASR)-min(Popcomp$ASR))
Popcomp$LatS2<-1-(Popcomp$Lat-min(Popcomp$Lat))/(max(Popcomp$Lat)-min(Popcomp$Lat))
Popcomp$tyoS2<-1-(Popcomp$tyo-min(Popcomp$tyo))/(max(Popcomp$tyo)-min(Popcomp$tyo))
Popcomp$unsexedS2<-1-(Popcomp$unsexed-min(Popcomp$unsexed))/(max(Popcomp$unsexed)-min(Popcomp$unsexed))
Popcomp$spfS2<-(Popcomp$spf-min(Popcomp$spf))/(max(Popcomp$spf)-min(Popcomp$spf))
Popcomp$slopfS2<-1-(Popcomp$slopf-min(Popcomp$slopf))/(max(Popcomp$slopf)-min(Popcomp$slopf))
Popcomp$slopmS2<-1-(Popcomp$slopm-min(Popcomp$slopm))/(max(Popcomp$slopm)-min(Popcomp$slopm))
Popcomp$densityS2<-(Popcomp$density-min(Popcomp$density))/(max(Popcomp$density)-min(Popcomp$density))
Popcomp$scleritesS2<-(Popcomp$sclerites-min(Popcomp$sclerites))/(max(Popcomp$sclerites)-min(Popcomp$sclerites))
Popcomp$interscelriteS2<-1-(Popcomp$interscelrite-min(Popcomp$interscelrite))/(max(Popcomp$interscelrite)-min(Popcomp$interscelrite))
Popcomp$conditionS2<-(Popcomp$condition-min(Popcomp$condition))/(max(Popcomp$condition)-min(Popcomp$condition))

broods<-cbind(Popcomp[,c("Area","Broodsfixed","BroodsfixedS","BroodsfixedS2")],Variable="Broods (SBF)")
colnames(broods)<-c("Area","Value","Scaled","Scaled2","Variable")
lats<-cbind(Popcomp[,c("Area","Lat","LatS","LatS2")],Variable="Latitude N")
colnames(lats)<-c("Area","Value","Scaled","Scaled2","Variable")
ancestry<-cbind(Popcomp[,c("Area","Ancestry","AncestryS","AncestryS2")],Variable="% Ancestry Group 1")
colnames(ancestry)<-c("Area","Value","Scaled","Scaled2","Variable")
asr<-cbind(Popcomp[,c("Area","ASR","ASRS","ASRS2")],Variable="ASR (% males)")
colnames(asr)<-c("Area","Value","Scaled","Scaled2","Variable")
toyo<-cbind(Popcomp[,c("Area","tyo","tyoS","tyoS2")],Variable="% 2 years old")
colnames(toyo)<-c("Area","Value","Scaled","Scaled2","Variable")
unsx<-cbind(Popcomp[,c("Area","unsexed","unsexedS","unsexedS2")],Variable="% Unsexed")
colnames(unsx)<-c("Area","Value","Scaled","Scaled2","Variable")
spfs<-cbind(Popcomp[,c("Area","spf","spfS","spfS2")],Variable="% Late Spawning")
colnames(spfs)<-c("Area","Value","Scaled","Scaled2","Variable")
slopef<-cbind(Popcomp[,c("Area","slopf","slopfS","slopfS2")],Variable="F.Cond. slope")
colnames(slopef)<-c("Area","Value","Scaled","Scaled2","Variable")
slopem<-cbind(Popcomp[,c("Area","slopm","slopmS","slopmS2")],Variable="M.Cond. slope")
colnames(slopem)<-c("Area","Value","Scaled","Scaled2","Variable")
densi<-cbind(Popcomp[,c("Area","density","densityS","densityS2")],Variable="Density (fish/m)")
colnames(densi)<-c("Area","Value","Scaled","Scaled2","Variable")
sclerit<-cbind(Popcomp[,c("Area","sclerites","scleritesS","scleritesS2")],Variable="Sclerites")
colnames(sclerit)<-c("Area","Value","Scaled","Scaled2","Variable")
interscle<-cbind(Popcomp[,c("Area","interscelrite","interscelriteS","interscelriteS2")],Variable="Intersclerite dist. (microm)")
colnames(interscle)<-c("Area","Value","Scaled","Scaled2","Variable")
condi<-cbind(Popcomp[,c("Area","condition","conditionS","conditionS2")],Variable="Condition")
colnames(condi)<-c("Area","Value","Scaled","Scaled2","Variable")


Popcomplong<-rbind(broods,lats,ancestry,asr,toyo,unsx,spfs,slopef,slopem, densi, sclerit,interscle, condi )
Popcomplong$ValueR<-NA
Popcomplong[Popcomplong$Variable!="Condition" & Popcomplong$Variable!="Density (fish/m)",]$ValueR<-round(Popcomplong[Popcomplong$Variable!="Condition" & Popcomplong$Variable!="Density (fish/m)",]$Value,0)
Popcomplong[Popcomplong$Variable=="Condition",]$ValueR<-round(Popcomplong[Popcomplong$Variable=="Condition" ,]$Value,2)
Popcomplong[Popcomplong$Variable=="Density (fish/m)",]$ValueR<-round(Popcomplong[Popcomplong$Variable=="Density (fish/m)" ,]$Value,1)
#Popcomplong[Popcomplong$Variable=="F.Cond. slope",]$ValueR<-c("+/-","+","+/-","+","+/-","+")
Popcomplong[Popcomplong$Variable=="F.Cond. slope",]$ValueR<-c("--","/","--","/","--","/")
#Popcomplong[Popcomplong$Variable=="M.Cond. slope",]$ValueR<-c("+/-","+","+/-","+","-","+")
Popcomplong[Popcomplong$Variable=="M.Cond. slope",]$ValueR<-c("--","/","--","/","\\","/")
#Popcomplong[Popcomplong$Variable=="Density (fish/m)",]$ValueR<-c("A,B","B","A","A","A,B","A")
Popcomplong$Variable<-as.factor(Popcomplong$Variable)
levels(Popcomplong$Variable)
Popcomplong$Variable<-factor(Popcomplong$Variable, levels=c("% Ancestry Group 1","% Late Spawning","ASR (% males)",
                                                            "% 2 years old","% Unsexed","Sclerites",
                                                            "Intersclerite dist. (microm)",
                                                            "Density (fish/m)", "Condition",
                                                            "F.Cond. slope",
                                                            "M.Cond. slope",
                                                            "Latitude N", 
                                                            "Broods (SBF)" ))
Popcomplong
summary(Popcomplong)

#heatmap
ggplot(data=Popcomplong,aes(y=Area, x=Variable,fill=Scaled2))+
  geom_tile(aes(y=Area, x=Variable,fill=Scaled2), show.legend = F)+
  geom_text(aes(label = as.character(ValueR), color=1-Scaled2), size=4, show.legend = F)+
  scale_fill_gradient(low="steelblue1",high="red3")+
  scale_color_gradient(low="steelblue1",high="red4")+
  scale_x_discrete(position = "top")+
  labs(y="",x="")+
  theme(axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size=10,angle = 30, vjust = -1, hjust=0))+
  theme(plot.margin = unit(c(0,1.3,0,0), "cm"))






heatmap1
ggplot(data=heatmap1,aes(X1, forcats::fct_rev(X2), fill = value)) + geom_tile() +
  geom_text(aes(label = value), size=3)+
  scale_fill_gradientn(colours = c("limegreen","white", "darkgrey"),
                       values = c(0,0.6,0.95,1), trans="log10")+
  labs(x="",y="")+
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 45, vjust = -1, hjust=0))

