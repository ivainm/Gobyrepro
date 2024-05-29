#project:Gobyrepro
#Analysis 1
#Supports bioinformatic pipeline


#---------------------------------------------------------\
#=========================================================|
#---------------------------------------------------------\
#=========================================================|
#Identification of suspicious samples prior to merging of technical replicates
#In combination with output from script: nomerge.sh
#With replicates not merged to identify problems
setwd("C://Users//ivainm//Working_Document//DYNAMAR//PopGen//data_analysis//nomerge")
library(ggplot2)
library(reshape2)
library(adegenet)
library(ggpubr)
library(vcfR)
library(LEA)
library(xadmix)
library(tidyr)

#Checking missingness per individual to see if some outlier inds should be removed
missing75<-read.delim("missing_75.imiss", header=T)
#population infos
popmapnomerge<-read.delim("../vcf/popmap_SORT_nomerge.txt", header=T, sep='')
missing75<-cbind(missing75,popmapnomerge)

ggplot(data=missing75,aes(x=INDV,y=F_MISS, color=Area))+
  geom_point()+
  facet_wrap(~as.factor(Replicate))

#First problem: replicate 1 has a lot of missing data so I can basically trash all of it.
#deciding for removing individuals with more than 255 missing data. Export their names to handle in plink
write.table( missing75[missing75$F_MISS>0.75,c(1,6)], "miss_out_75.txt",sep=" ", quote=F, col.names=F, row.names=F)
write.table( missing75[missing75$F_MISS>0.75,1], "miss_out_75_vcf.txt",sep=" ", quote=F, col.names=F, row.names=F)
head(missing75)

#PCAs -With outliers excluded
#import output from plink
pca90_nomerge<-read.delim("pca90_nomerge_out.eigenvec", header=F, sep= ' ')
#add popmap info for the individuals that have passed the missingness filter
pca90_nomerge<-cbind(pca90_nomerge,popmapnomerge[!(popmapnomerge$ID %in% missing75[missing75$F_MISS>0.75,1]),])
pca90_nomerge$Area<- factor(pca90_nomerge$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
pca90_nomerge$Replicate<-as.factor(pca90_nomerge$Replicate)
pca90_nomerge$ID_rep<-as.factor(pca90_nomerge$ID_rep)

#plot
ggplot(data=pca90_nomerge,aes(x=V3,y=V4))+
  geom_point(aes(fill=Area,colour=Area), shape=21, size=3)+
  stat_ellipse(level = 0.95, aes(color=Area))+
  facet_wrap(~as.factor(Replicate))+
  labs(x="PC 1",y="PC 2",
       title="PCA stringent 90% missingness, outliers removed")

#Comparing the loadings of replicates on the different PC axes
#I need to melt and cast the data 
str(pca90_nomerge)
pca90_nomerge_min<-pca90_nomerge[,c(3,23,26)]
pca90_nomerge_wide<-spread(pca90_nomerge_min,Replicate,V3)
#population infos
popmap1<-read.delim("../vcf/pop_map_SORT.txt", header=T, sep='')
#beware of the order of the individuals when merging pop info
#pca90_nomerge_wide<-cbind(pca90_nomerge_wide,popmap1[,c(1,2,3,4,6)])
colnames(pca90_nomerge_wide)<-c("ID","Rep1","Rep2","Rep3")
str(pca90_nomerge_wide)

ggplot(data=pca90_nomerge_wide, aes(x=Rep1,y=Rep2))+
  geom_point()+
  geom_text(data=pca90_nomerge_wide[pca90_nomerge_wide$ID=="SORT-80" | pca90_nomerge_wide$ID=="SORT-82",],aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  #  geom_text(aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  labs(title="PC1 compared across replicate 1 and 2") 

ggplot(data=pca90_nomerge_wide, aes(x=Rep2,y=Rep3))+
  geom_point()+
  #geom_text(aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  geom_text(data=pca90_nomerge_wide[pca90_nomerge_wide$ID=="SORT-80" | pca90_nomerge_wide$ID=="SORT-82",],aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  labs(title="PC1 compared across replicate 2 and 3") 

ggplot(data=pca90_nomerge_wide, aes(x=Rep1,y=Rep3))+
  geom_point()+
  #geom_text(aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  geom_text(data=pca90_nomerge_wide[pca90_nomerge_wide$ID=="SORT-80" | pca90_nomerge_wide$ID=="SORT-82",],aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  labs(title="PC1 compared across replicate 1 and 3") 

#same as above with PC2
pca90_nomerge_min2<-pca90_nomerge[,c(4,23,26)]
pca90_nomerge_wide2<-spread(pca90_nomerge_min2,Replicate,V4)
colnames(pca90_nomerge_wide2)<-c("ID","Rep1","Rep2","Rep3")
str(pca90_nomerge_wide2)

ggplot(data=pca90_nomerge_wide2, aes(x=Rep1,y=Rep2))+
  geom_point()+
  #geom_text(aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  geom_text(aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  labs(title="PC2 compared across replicate 1 and 2") 

ggplot(data=pca90_nomerge_wide2, aes(x=Rep2,y=Rep3))+
  geom_point()+
  geom_text(aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  labs(title="PC2 compared across replicate 2 and 3") 

ggplot(data=pca90_nomerge_wide2, aes(x=Rep1,y=Rep3))+
  geom_point()+
  geom_text(aes(label=ID), nudge_y = 0.002,nudge_x = 0.01)+
  labs(title="PC2 compared across replicate 1 and 3") 
#The consistency of PC loadings is then used visually to exclude inconsistent replicates, or individuals that do not have at least 2/3 consistent replicates.

#---------------------------------------------------------\
#=========================================================|
#---------------------------------------------------------\
#=========================================================|

##plot depth  use SNP_DP file

#---------------------------------------------------------\
#=========================================================|
#---------------------------------------------------------\
#=========================================================|
#Plot missingness to find outliers

#-----------------------------------------------------------------------------------
#Checking missingness per individual to see if some outlier inds should be removed
missing_maf<-read.delim("missingmaf5_90.imiss", header=T)
#population infos
popmap<-read.delim("pop_map_final.txt", header=T, sep='')
missing_maf<-cbind(missing_maf,popmap)

ggplot(data=missing_maf,aes(x=INDV,y=F_MISS, color=Area))+
  geom_point()+
  facet_wrap(~as.factor(Timepoint))+
  geom_text( aes(label=ID, color=Area), nudge_y = 0.02,nudge_x = 0.01)+
  labs(x="ID",y="Misingness",  title=" With Maf 0.05 cutoff")

#4 individual really stick out with more than 40% missing data, but still less than 60%

#deciding for removing individuals with more than 40% missing data. Export their names to handle in plink
write.table( missing_maf[missing_maf$F_MISS>0.4,c(1,6)], "miss_out_40.txt",sep=" ", quote=F, col.names=F, row.names=F)
write.table( missing_maf[missing_maf$F_MISS>0.4,1], "miss_out_40_vcf.txt",sep=" ", quote=F, col.names=F, row.names=F)




#---------------------------------------------------------\
#=========================================================|
#---------------------------------------------------------\
#=========================================================|

#Hz plot
# Observed vs expected homozygosity F / clculated by plink as observed homozygosity minus expected.
#Positive is excess of Homozygosity
hom<-read.table("het90_maf.het", header=T)
popmap<-read.delim("pop_map_final.txt", header=T, sep='')
hom<-cbind(hom,popmap[!(popmap$ID %in% missing_mac[missing_mac$F_MISS>0.4,]$ID),])
hom$Area<- factor(hom$Area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
str(hom)

ggplot(data=hom, aes(x=F))+
  geom_histogram(aes(fill=Area))+
  facet_wrap(~Timepoint)+
  labs(x="F (Excess of homozygosity)",y="Count",
       title="Histogram of Excess homozygosity")

ggplot(data=hom, aes(x=Area, y=F))+
  geom_boxplot(aes(fill=Area))+
  facet_wrap(~Timepoint)+
  labs(x="Area",y="F (Excess of homozygosity)",
       title="Histogram of Excess homozygosity")



#---------------------------------------------------------\
#=========================================================|
#---------------------------------------------------------\
#=========================================================|

#fst calculations

#Pairwise FST again------------------------------------------------------------(from: final2.R)
rm(list=ls())
#install assigner
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/assigner")
#install radiator
install.packages("remotes")
remotes::install_github("thierrygosselin/radiator")
library(assigner)
library(radiator)
library(progressr)

#check that my file can be read
radiator::detect_genomic_format(data ="R90_maf.bed")

# pairwise FST calculation 
#NOTE: for fst_wc84 to work, the population (strata) has to be defined in plink using the --update-ids function.
tidy_09maf<-tidy_genomic_data(
  "R90_maf.bed",  strata = "pop_conv2.txt",  filename = NULL,  parallel.core = parallel::detectCores() - 2,  verbose = TRUE)

tidy_09maf$GT<-tidy_09maf$GT_BIN

str(tidy_09maf2)
levels(tidy_09maf$STRATA)

fst_09maf<-fst_WC84( tidy_09maf,  snprelate = FALSE,  strata = NULL,  pop.levels = c("Kristineberg","Arendal","Austevoll","Hitra","Helligvaer","Ringstad"), 
                     pairwise = T,  ci = T,  iteration.ci = 100,  quantiles.ci = c(0.025, 0.975),
                     heatmap.fst = T,  digits = 9,  filename = NULL,  parallel.core = parallel::detectCores() - 2,  verbose = T)




