
rm(list=ls())

# Top ---------------------------------------------------------------------

setwd("C:/Users/nh1087/OneDrive - USNH/Documents/NECC/Diet Data/")


library(tidyverse)
library(ggplot2)
library(readr)
library(stringr)
library(lubridate)
library(ggmap)
library(devtools)
library(grid)
library(moments)
library(vegan)
library(indicspecies)
library(usedist)
library(ggpattern)
library(RVAideMemoire)
library(remotes)
library(gridExtra)
library(RColorBrewer)
library(sf)
library(sp)
library(ggspatial)
library(geosphere)
library(rgeos)




#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}

unique(as.character(prey19$pdcomnam))

seasonPal<-c("steelblue4","goldenrod1","deeppink3","sienna2")
seasonPal2<-c("deepskyblue4","yellowgreen","goldenrod1","orange3") #BEST
seasonPal3<-c("deepskyblue4","yellowgreen","paleturquoise2","orange3")

#Data loading
load("prey19.RData")
load("pylen19.RData")
load(here("data", "shapefiles","googleMap.zoom6.eastCoast.R"))
vulScores<-read.csv(here("data", "raw", "Hareetal_climateScores.csv"))


#ALL LEVELS OF YEAR AND SEASON
ALLyearseasons<-factor(levels=c(paste(sort(rep(seq(1973,2019),4)),rep(c("Winter","Spring","Summer","Fall"),length(seq(1973,2019))))))
levels(ALLyearseasons)

#The names of my species, just to have them nice and handy
myspecies_sci<-unique(str_to_sentence(prey19$pdscinam))
myspecies_com<-unique(str_to_title(prey19$pdcomnam))
myspecies_svspp<-unique(str_pad(as.character(prey19$svspp),width=3,pad="0",side="left"))

#Size classes used by Garrison and Link
sizeClasses<-read.csv(here("data", "raw", "GarrisonLink_predSizeClasses.csv"))
#These are replicated in the sizecats column already present, except that the small and largest cats extend to the smallest and largest individuals
prey19%>%group_by(pdcomnam,sizecat)%>%summarise(min=min(pdlen),max=max(pdlen))

load("NF.prey19.RData")

#How many in new data
table(NF.prey19$pdcomnam)
n_distinct(NF.prey19$pdcomnam)-n_distinct(prey19$pdcomnam) #how many more species
NF.prey19<-NF.prey19%>%
  mutate(svspp=str_pad(svspp,width=3,side="left",pad="0"),
         id=paste0(cruise6,
                   str_pad(station,width=3,side="left",pad="0"),
                   str_pad(stratum,width=3,side="left",pad="0")),
         dietID=paste0(svspp,
                       pdsex,
                       str_pad(pdid,width=6,side="left",pad="0"),
                       str_pad(pdlen,width=3,pad="0",side="left"),
                       id))
## CSL: preyGMRI_filter not found
# n_distinct(NF.prey19$dietID)-n_distinct(preyGMRI_filter$dietID)#How many more diets
indNF.prey19<-NF.prey19%>%
  dplyr::select(-c(pynam:pyamtv))%>%
  distinct()%>%
  mutate(cruise6=as.character(cruise6),
         station=str_pad(station,width=4,pad="0",side="left"),
         pdsex=as.character(pdsex),
         pdid=str_pad(pdid,width=6,pad="0",side="left"))%>%
  filter(!duplicated(dietID)) #There is one that was both examined at sea AND preserved for lab????
dietCounts<-as.data.frame(sort(table(indNF.prey19$pdcomnam)))%>%
  mutate(predNum=order(Var1),
         predLab.1=ifelse(predNum%%2==1,str_to_sentence(as.character(Var1)),NA),
         predLab.2=ifelse(predNum%%2==0,str_to_sentence(as.character(Var1)),NA),
         new=ifelse(Var1%in%prey19$pdcomnam,"black","red"))

ggplot()+
  geom_col(data=dietCounts,aes(predNum,Freq,fill=Var1),color="black")+
  geom_text(data=dietCounts,aes(predNum,Freq+1500,label=Freq,color=new),angle=90,hjust=0.4)+
  scale_x_continuous(name="Predator Species",expand=expansion(mult=c(0.025,0.025)),
                     labels=filter(dietCounts,!is.na(predLab.1))$predLab.1,
                     breaks=seq(1,nrow(dietCounts),by=2),
                     sec.axis=dup_axis(name="",
                                       labels=filter(dietCounts,!is.na(predLab.2))$predLab.2,
                                       breaks=seq(2,nrow(dietCounts),by=2)))+
  scale_y_continuous(expand=expansion(add=c(200,2500)),name="Number of Diets")+
  scale_fill_manual(values=colorspace::rainbow_hcl(n=41,c=100))+
  scale_color_manual(values=c("black","red"))+
  theme(legend.position = "none",
        axis.text.x.bottom=element_text(angle=20,hjust=0.8,vjust=1,size=18),
        axis.text.x.top=element_text(angle=20,hjust=0.2,vjust=0,size=18),
        plot.margin=margin(t=5,r=25,b=5,l=75))

# Species SVSPP Codes -----------------------------------------------------

spec<- read.csv(here("data", "raw", "22560_UNION_FSCS_SVCAT.csv")) |>
  # read.csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv")%>%
  dplyr::select(LOGGED_SPECIES_NAME,SVSPP)%>%distinct()%>%
  filter(SVSPP>="001")%>%
  mutate(LOGGED_SPECIES_NAME=str_to_sentence(LOGGED_SPECIES_NAME))
spec<-spec[!duplicated(spec$SVSPP),]


# Important df manipulations ----------------------------------------------

prey19<-prey19%>%
  mutate(gensci=ifelse(pynam=="PRIONOTUS ALATUS"|pynam=="STERNOPTYCHIDAE","FISH",as.character(gensci)),
         analsci=ifelse(pynam=="PRIONOTUS ALATUS","TRIGLIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="PRIONOTUS ALATUS",as.character(pynam),as.character(collsci)),
         analsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(collsci)),
         year=ifelse(is.na(year),substr(cruise6,1,4),year),year=as.numeric(year))%>%
  filter(pdlen<300) #There is a Summer Flounder that was listed at 300 cm, which is unbelievable, so it will be dropped


uniquePrey19<-prey19%>%
  dplyr::select(cruise6,station,svspp,pdsex,pdid,pdcomnam,pdscinam,pdlen,pdwgt,sizecat,pdgutw,pdgutv,declat,declon,month,day,year,season,geoarea)%>%
  distinct()%>%
  mutate(cruise6=as(cruise6,"character"),pdsex=as.character(pdsex),
         station=as.character(str_sub(paste0("0000",station),-4,-1)))



# Creating my own prey categories -----------------------------------------

#The categories from Garrison and Link, 2000
gl_preycats<-read.csv(here("data","raw","GarrisonLink_preyCats.csv"))%>%
  mutate(matchingCats=gsub("p\\.","",Scientific.name),
         matchingCats=gsub("crabs","crab",matchingCats),
         matchingCats=gsub("Gammaridae","Gammaridea",matchingCats),
         matchingCats=gsub("Cnidarians","Cnidaria",matchingCats))

prey19<-prey19%>%
  mutate(INgen=ifelse(str_to_sentence(gensci) %in% gl_preycats$matchingCats,1,0),
         INanal=ifelse(str_to_sentence(analsci) %in% gl_preycats$matchingCats,1,0),
         INcoll=ifelse(str_to_sentence(collsci) %in% gl_preycats$matchingCats,1,0),
         INpy=ifelse(str_to_sentence(pynam) %in% gl_preycats$matchingCats,1,0))
prey19<-prey19%>%
  mutate(INnum=rowSums(prey19[,c("INgen","INanal","INcoll","INpy")]),
         gl_prey=ifelse(INnum==4,str_to_sentence(pynam), #STEP 1 
                  ifelse(INnum==3&INpy==1,str_to_sentence(pynam),
                   ifelse(INnum==3&INpy==0,str_to_sentence(collsci),
                    ifelse(INnum==2&INpy==1,str_to_sentence(pynam),
                     ifelse(INnum==2&INgen==1,str_to_sentence(analsci),
                      ifelse(INnum==2&INgen==0&INpy==0,str_to_sentence(collsci),
                       ifelse(INnum==1&INgen==1,str_to_sentence(gensci),
                        ifelse(INnum==1&INanal==1,str_to_sentence(analsci),
                         ifelse(INnum==1&INcoll==1,str_to_sentence(collsci),
                          ifelse(INnum==1&INpy==1,str_to_sentence(pynam),
                           ifelse(INnum==0&pynam=="EMPTY","Empty","Unobserved"))))))))))),
         gl_prey=ifelse(pynam %in% c("UROPHYCIS CHUSS","UROPHYCIS TENUIS","UROPHYCIS REGIA"),"Other hakes", #STEP 2
                  ifelse(analsci %in% c("GADIDAE","BREGMACEROTIDAE","EUCLICHTHYIDAE","LOTIDAE","MACROURIDAE",
                                        "MELANONIDAE","MERLUCCIIDAE","MORIDAE","MURAENOLEPIDIDAE","PHYCIDAE") & 
                           is.na(gl_prey),"Gadiformes", #STEP 3a
                   ifelse(analsci %in% c("PLEURONECTIDAE","PSETTODIDAE","CITHARIDAE","SCOPHTHALMIDAE","PARALICHTHYIDAE",
                                         "BOTHIDAE","PARALICHTHODIDAE","POECILOPSETTIDAE","RHOMBOSOLEIDAE",
                                         "ACHIROPSETTIDAE","SAMARIDAE","ACHIRIDAE","SOLEIDAE","CYNOGLOSSIDAE") & 
                            is.na(gl_prey), "Pleuronectiformes", #STEP 3b
                    ifelse(gensci=="FISH" & is.na(gl_prey),"Other fish", #STEP 4
                     ifelse((gensci %in% c("UROCHORDATA","BRACHIOPODA","BRYOZOA","CHAETOGNATHA","PORIFERA") | 
                              collsci %in% c("ARTHROPODA","INSECTA","HEMICHORDATA","LIMULUS POLYPHEMUS",
                                             "APHRODITIDAE","OLIGOCHAETA","HIRUDENEA","PYCNOGONIDA")) & 
                           is.na(gl_prey),"Other invertebrates", #STEP 5
                      ifelse(analsci %in% c("CEPHALOCHORDATA"),"Other", #STEP 6
                       ifelse(collsci %in% c("OSTRACODA","CUMACEA","STOMATOPODA","PENAEIDAE"),
                              "Crustacean shrimp", #STEP 7
                        ifelse(analsci %in% c("CIRRIPEDIA","COPEPODA") | 
                                 pynam %in% c("DECAPODA","DECAPODA EGGS","DECAPODA LARVAE"),"Crustacea", #STEP 8
                         ifelse(collsci %in% c("HOMARUS AMERICANUS","CALLINECTES SAPIDUS","DECAPODA LARVAE","SCYLLARIDAE") &
                                  is.na(gl_prey), "Decapoda crab", #STEP 9
                          ifelse(analsci=="EUPHAUSIACEA","Euphausiidae",gl_prey))))))))))) #STEP 10



#Checking them all
check<-unique(prey19[which(is.na(prey19$gl_prey)),c("gensci","analsci","collsci","pynam")])
nrow((prey19[which(is.na(prey19$gl_prey)),c("gensci","analsci","collsci","pynam")])) #Just the unobserved, probably should be dropped

gl_preycats[gl_preycats$matchingCats %notin% prey19$gl_prey,] #Nothing went into Zooplankton in the GL categories
unique(prey19[prey19$gl_prey %notin% gl_preycats$matchingCats,"gl_prey"]) #Empty and NA are "new" but easy to keep out so good

#Helpful to figure out #Step 1, counts how many of each level of matching the GL categories occurs
test2<-(prey19[,c("INgen","INanal","INcoll","INpy","INnum")])%>%
  group_by(INgen,INanal,INcoll,INpy,INnum)%>%
  mutate(N=n(),p=N/nrow(prey19))%>%unique()
#How many match immediately
(sum(filter(test2,INnum>0)$N)+nrow(prey19[prey19$pynam=="EMPTY",]))/nrow(prey19) 
#96.7%

amph<-filter(prey19,analsci=="AMPHIPODA")%>%dplyr::select(pynam:collsci)%>%unique()
unobs<-filter(prey19,analsci=="UNOBS")

testGroups<-prey19%>%
  mutate(totalV=sum(pyamtv))%>%
  group_by(gl_prey)%>%
  reframe(vol=sum(pyamtv),p=vol/totalV*100)%>%unique()

spGL_Groups<-prey19%>%
  group_by(pdcomnam)%>%
  mutate(totalV=sum(pyamtv))%>%
  group_by(pdcomnam,gl_prey)%>%
  reframe(vol=sum(pyamtv),p=vol/totalV*100)%>%unique()




# Summaries of data availability ------------------------------------------

summary(prey19$pdcomnam)


prey19%>%
  group_by(pdcomnam,gensci)%>%
  summarise(pynum=sum(pynum,na.rm=T))%>%
  ggplot(aes(pdcomnam,pynum,fill=gensci))+
  geom_col(color="black")+
  scale_fill_viridis_d()

#How many stomachs for each species
prey19<-prey19%>%
  mutate(uniqueID=paste(svspp,pdsex,pdid,sep="-"))

uniquePrey19<-prey19%>%
  dplyr::select(cruise6,station,svspp,pdsex,pdid,pdcomnam,pdscinam,pdlen,sizecat,pdgutw,pdgutv,declat,declon,month,day,year,season,geoarea)%>%
  distinct()%>%
  mutate(cruise6=as(cruise6,"character"),pdsex=as.character(pdsex),
         station=as.character(str_sub(paste0("0000",station),-4,-1)))




#Number of diets from each species
uniquePrey19%>%
  mutate(pdcomnam=str_to_title(gsub(" ","\n",pdcomnam)))%>%
  ggplot(aes(pdcomnam,fill=pdcomnam))+
  geom_bar(color="black",size=1)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,l=60),guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Number of Diets Collected (thous.)",
                     labels=c(0,20,40,60,80))+
  scale_x_discrete(name="Predator Species")


#Actual Season periods (THIS NEEDS SOME QAQC for sure)
uniquePrey19%>%
  mutate(doy=yday(ymd(paste(year,month,day,sep="-"))),
         season=str_to_title(season),
         season=factor(season,levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(season)%>%
  summarise(seasonStart=min(doy,na.rm=T),seasonEnd=max(doy,na.rm=T))%>%
  ggplot()+
  geom_segment(aes(x=seasonStart,xend=seasonEnd,y=1,yend=1,color=season),size=10,alpha=0.8)+
  scale_x_continuous(limits=c(0,366),expand=expansion(0),name="Date",
                     breaks=c(0,60,152,244,335),labels=c("Jan","Mar","Jun","Sep","Dec"))+
  scale_y_continuous(expand=expansion(0))+
  scale_color_manual(values=seasonPal2,name="Season")+
  theme(axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),
        panel.grid.minor=element_blank(),panel.grid.major.y=element_blank())


#Number of diets collected in each day of the year
uniquePrey19%>%
  mutate(doy=yday(ymd(paste(year,month,day,sep="-"))),
         season=str_to_title(season),
         season=factor(season,levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(season,doy)%>%
  summarise(nCatch=n())%>%
  ggplot()+
  geom_segment(aes(x=doy,xend=doy,y=-nCatch/2,yend=nCatch/2,color=season),size=1,alpha=0.8)+
  scale_x_continuous(limits=c(0,366),expand=expansion(0),name="Date",
                     breaks=c(0,60,152,244,335),labels=c("Jan","Mar","Jun","Sep","Dec"))+
  scale_y_continuous(expand=expansion(0.1),name="Number of Diets",
                     breaks=c(-2000,-1000,0,1000,2000),labels=c(4000,2000,0,2000,4000))+
  scale_color_manual(values=seasonPal2,name="Season")+
  theme(panel.grid.minor.x=element_blank())+
  guides(color=guide_legend(override.aes=list(size=10)))


#Number of diets from each species in each season
uniquePrey19%>%
  mutate(pdcomnam=str_to_title(gsub(" ","\n",pdcomnam)),
         season=factor(str_to_title(season),levels=c("Winter","Spring","Summer","Fall")))%>%
  ggplot(aes(pdcomnam,fill=season))+
  geom_bar(color="black",size=1,position=position_dodge())+
  scale_fill_manual(name="Season",values=seasonPal2)+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Number of Diets Collected")+
  scale_x_discrete(name="Predator Species")+
  theme(legend.position=c(0.2,0.69),legend.text=element_text(size=30),legend.title=element_text(size=33),
        legend.background = element_rect(color="black"))
#Calculating these for myself, to make a table
table(uniquePrey19$season,uniquePrey19$pdcomnam)



#Proportion of diets for each species in each season
uniquePrey19%>%
   mutate(pdcomnam=str_to_title(gsub(" ","\n",pdcomnam)),
         season=factor(str_to_title(season),levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(pdcomnam)%>%
  mutate(total=n())%>%
  group_by(pdcomnam,season)%>%
  reframe(prop=n()/total)%>%
  distinct()%>%
  ggplot(aes(pdcomnam,prop,fill=season))+
  geom_col(color="black",size=1)+
  scale_fill_manual(values=seasonPal2,name="Season")+
  scale_y_continuous(expand=expansion(mult=c(0)),name="Proportion of Diets Collected")+
  scale_x_discrete(name="Predator Species")


#Number of diets in each season for each year
uniquePrey19%>%
  mutate(season=factor(str_to_title(season),levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(year,season)%>%
  summarise(`Number of Diets Collected`=n())%>%
  ggplot(aes(year,`Number of Diets Collected`,color=season))+
  geom_line()+
  geom_point()+
  scale_color_manual(values=seasonPal2,guide="none")+
  facet_wrap(~season)



#How many diets over time
min<-uniquePrey19%>%
  filter(!is.na(year))%>%
  mutate(pdcomnam=str_to_title(pdcomnam),
         year=as.factor(year))%>%
  group_by(pdcomnam,year,.drop=F)%>%
  summarise(count=n())%>%distinct()%>%
  group_by(pdcomnam)%>%
  mutate(year=as.character(year))%>%
  filter(year>=min(subset(year,count>0)))%>%
  reframe(M=min(count,na.rm=T),
            Myear=subset(year,count==min(count,na.rm=T)))%>%
  mutate(Myear=ifelse(n()>1,paste(n(),"times"),Myear))%>%distinct()


uniquePrey19%>%
  filter(!is.na(year))%>%
  mutate(pdcomnam=str_to_title(pdcomnam))%>%
  group_by(pdcomnam,year)%>%
  summarise(count=n())%>%distinct()%>%
  mutate(g = c(0, cumsum(diff(year) > 1))) %>%
  ggplot()+
  geom_line(aes(year,count,color=pdcomnam,group=g),size=2)+
  geom_point(aes(year,count,color=pdcomnam,group=g),size=4)+
  geom_text(data=min,aes(2010,3333,label=paste0("Min = ",M,": ",Myear)),size=7)+
  scale_color_manual(values=colorspace::rainbow_hcl(9,c=100,l=60),guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0.035,0.1)),name="Number of Diets Collected")+
  scale_x_continuous(name="Year")+
  facet_wrap(~pdcomnam)



#How many diets over time in each season
minSeason<-uniquePrey19%>%
  filter(!is.na(year))%>%
  mutate(pdcomnam=str_to_title(pdcomnam),
         season=str_to_title(season),
         yearseason=factor(paste(year,season)))%>%
  group_by(pdcomnam,yearseason,.drop=F)%>%
  summarise(count=n())%>%distinct()%>%
  group_by(pdcomnam)%>%
  filter(as.numeric(yearseason)>=min(as.numeric(subset(yearseason,count>0))))%>%
  reframe(M=min(count,na.rm=T),Myearseason=subset(yearseason,count==min(count,na.rm=T)))%>%
  mutate(Myearseason=ifelse(n()>1,paste(n(),"times"),as.character(Myearseason)))%>%distinct()


uniquePrey19%>%
  filter(!is.na(year))%>%
  mutate(pdcomnam=str_to_title(pdcomnam),
         Season=str_to_title(season),
         Season=factor(Season,levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(pdcomnam,Season,year)%>%
  summarise(count=n())%>%distinct()%>%
  mutate(g = c(0, cumsum(diff(year) > 1)),
         sg=paste0(Season,g)) %>%
  ggplot()+
  geom_line(aes(year,count,color=Season,group=sg))+
  geom_point(aes(year,count,color=Season))+
  geom_text(data=minSeason,aes(2005,2005,label=paste0("Min = ",M,": ",Myearseason)),size=7)+
  scale_color_manual(values=seasonPal2)+
  scale_y_continuous(expand=expansion(mult=c(0.035,0.1)),name="Number of Diets Collected")+
  scale_x_continuous(name="Year")+
  facet_wrap(~pdcomnam)


#Where are the species caught
load("googleMap.zoom6.eastCoast.R")
ggmap(zoom6)+
  geom_point(data=filter(uniquePrey19,!is.na(year))%>%mutate(Species=str_to_title(pdcomnam)),
             aes(-declon,declat,fill=pdcomnam),shape=21,alpha=0.8,size=3,show.legend=F)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,end=290))+
  ylab("Latitude")+xlab("Longitude")+
  facet_wrap(~Species)



#How many diets for the different prey categories
#At the general level...
prey19%>%
  mutate(gensci=str_to_title(gsub(" ","\n",ifelse(gensci=="",as.character(pynam),as.character(gensci)))))%>%
  ggplot(aes(gensci))+
  geom_bar(color="black")+
  theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))
#Analytical level...
genCounts<-prey19%>%
  mutate(analsci=str_to_sentence(ifelse(analsci=="","Other Fish",as.character(analsci))),
         gensci=str_to_title(ifelse(gensci=="","Fish",as.character(gensci))))%>%
  group_by(gensci)%>%
  mutate(width=n_distinct(analsci),
         N=n())%>%
  dplyr::select(gensci,width,N)%>%unique()%>%
  arrange(desc(width),gensci)
NMat<-matrix(genCounts$N,ncol=4,byrow=F)

appender <- function(string, suffix = as.vector(t(NMat))[1:nrow(genCounts)]) paste(string, suffix, sep=": n=")

g<-prey19%>%
  mutate(analsci=str_to_sentence(ifelse(analsci=="","Other Fish",as.character(analsci))),
         gensci=str_to_title(ifelse(gensci=="","Fish",as.character(gensci))))%>%
  group_by(gensci)%>%
  mutate(width=n_distinct(analsci),
         N=n())%>%
  arrange((N))%>%
  ggplot(aes(analsci))+
  geom_bar(color="black")+
  theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
  facet_wrap(~fct_reorder2(gensci,N,width),scales="free",labeller=as_labeller(appender),dir="v")+
  ylab("Count")+xlab("Analytical Category")+theme(strip.text=element_text(size=30),axis.title=element_text(size=33))
g
gt = ggplot_gtable(ggplot_build(g))
gt$widths[5] = gt$widths[5]*4
grid.draw(gt)
#Collection level...
genCounts2<-prey19%>%
  mutate(collsci=str_to_sentence(ifelse(collsci=="",as.character(pynam),as.character(collsci))),
         gensci=str_to_title(ifelse(gensci=="","Fish",as.character(gensci))))%>%
  group_by(gensci)%>%
  mutate(width=n_distinct(collsci),
         N=n())%>%
  dplyr::select(gensci,width,N)%>%unique()%>%
  arrange(desc(width),gensci)
NMat<-matrix(genCounts2$N,ncol=4,byrow=F)

appender <- function(string, suffix = as.vector(t(NMat))[1:nrow(genCounts2)]) paste(string, suffix, sep=": n=")

g<-prey19%>%
  mutate(collsci=str_to_sentence(ifelse(analsci=="",as.character(pynam),as.character(collsci))),
         gensci=str_to_title(ifelse(gensci=="","Fish",as.character(gensci))))%>%
  group_by(gensci)%>%
  mutate(width=n_distinct(collsci),
         N=n())%>%
  arrange((N))%>%
  ggplot(aes(collsci))+
  geom_bar(color="black")+
  theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))+
  facet_wrap(~fct_reorder2(gensci,N,width),scales="free",labeller=as_labeller(appender),dir="v")+
  ylab("Count")+xlab("Collection Category")+theme(strip.text=element_text(size=30),axis.title=element_text(size=33))
g
gt = ggplot_gtable(ggplot_build(g))
gt$widths[5] = gt$widths[5]*10
grid.draw(gt)
#And the actual prey name
ggplot(prey19,aes(pynam))+
  geom_bar(color="black")


#How many different levels are there in all the different levels
itemDetails<-prey19%>%
  mutate(gensci=ifelse(pynam=="PRIONOTUS ALATUS"|pynam=="STERNOPTYCHIDAE","FISH",as.character(gensci)),
         analsci=ifelse(pynam=="PRIONOTUS ALATUS","TRIGLIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="PRIONOTUS ALATUS",as.character(pynam),as.character(collsci)),
         analsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(collsci)),
         total_pynam=n_distinct(pynam),
         total_gennam=n_distinct(gensci),
         total_analnam=n_distinct(analsci),
         total_colnam=n_distinct(collsci))%>%
  group_by(gensci)%>%
  mutate(gen_pynam=n_distinct(pynam),
         gen_analnam=n_distinct(analsci),
         gen_colnam=n_distinct(collsci))%>%
  ungroup()%>%
  group_by(analsci)%>%
  mutate(anal_pynam=n_distinct(pynam),
         anal_colnam=n_distinct(collsci))%>%
  ungroup()%>%
  group_by(collsci)%>%
  mutate(coll_pynam=n_distinct(pynam))%>%
  dplyr::select(pynam,gensci,analsci,collsci,total_pynam:coll_pynam)%>%unique()
itemDetails<-itemDetails[,c("gensci","analsci","collsci","pynam","total_gennam","total_analnam","total_colnam",
                            "total_pynam","gen_analnam","gen_colnam","gen_pynam","anal_colnam","anal_pynam","coll_pynam")]  

#write.csv(itemDetails,"preyDetail_counts.csv",row.names = F)

itemLists<-prey19%>%
  dplyr::select(gensci,analsci,collsci,pynam,pyamtw)%>%
  group_by(gensci)%>%
  mutate(genamtw=sum(pyamtw))%>%
  group_by(gensci,analsci)%>%
  mutate(anamtw=sum(pyamtw))%>%
  group_by(gensci,analsci,collsci)%>%
  mutate(colamtw=sum(pyamtw))%>%
  group_by(gensci,genamtw,analsci,anamtw,collsci,colamtw,pynam)%>%
  summarise(pyamtw=sum(pyamtw))%>%
  unique()
#write.csv(itemLists,"listallPrey.csv",row.names = F) 


#Plastic etc. over time
filter(prey19,collsci=="MISCELLANEOUS" & !grepl("TUBES",pynam) & pynam!="ROCK" & pynam!="SAND" & !grepl("OOD",pynam))%>%
  ggplot(aes(year,pyamtw))+geom_point()
  



#How often are the predators showing up as prey in predator diets
speciesPred<-prey19%>%
  mutate(pynam=as.character(pynam),pdscinam=as.character(pdscinam),
         selfother=ifelse(pynam==pdscinam,"Self",
                          ifelse(pynam%in%myspecies_sci,"Diff Pred","Other")))%>%
  group_by(pdcomnam)%>%
  mutate(total_items=n(),
         pdcomnam=factor(pdcomnam),selfother=factor(selfother))%>%
  count(pdcomnam,selfother,total_items,.drop=F,name="n_items")%>%
  mutate(prop_items=n_items/total_items)%>%select(-total_items)%>%
  mutate(pdcomnam=gsub(" ","\n",str_to_title(pdcomnam)))
speciesPred[is.na(speciesPred)]<-0

ggplot(speciesPred,aes(pdcomnam,prop_items,fill=selfother))+
  geom_col(color="black",size=0.5)+
  scale_x_discrete(name="Predator Species",expand=expansion(0.07))+
  scale_y_continuous(name="Proportion of Diet Items",expand=expansion(add=0.007))+
  scale_fill_viridis_d(name="Item Type")+
  theme(legend.margin = margin(0,0,0,-10))



#How often are each of the predators prey to each of the predators
speciesPrey<-prey19%>%
  group_by(pdscinam)%>%
  mutate(total_items=n())%>%
  filter(pynam%in%myspecies_sci)%>%
  mutate(pynam=factor(pynam,levels=myspecies_sci),
         pdscinam=factor(pdscinam,levels=myspecies_sci))%>%
  count(pynam,pdscinam,total_items,.drop=F,name="n_prey")%>%
  mutate(prop_prey=n_prey/total_items,
         preyLab=str_to_sentence(pynam),
         preyLab=factor(preyLab,levels=str_to_sentence(myspecies_sci)),
         pdscinam=gsub(" ","\n",str_to_sentence(pdscinam)))%>%
  group_by(pdscinam)%>%
  mutate(total_prop=sum(prop_prey,na.rm=T))
speciesPrey[is.na(speciesPrey)]<-0


ggplot(speciesPrey,aes(fct_reorder(pdscinam,total_prop),prop_prey,fill=preyLab))+
  geom_col(position = "dodge",color="black",size=0.5)+
  geom_errorbar(aes(ymin=total_prop,ymax=total_prop), color = "firebrick3") +
  #geom_segment(aes(x=pdscinam-1,xend=pdscinam,y=total_prop,yend=total_prop),size=5,shape=22,fill="firebrick2")+
  scale_x_discrete(name="Predator Species")+
  scale_y_sqrt(name="Proportion of Predator Diet",expand=expansion(add=c(0.001,0.01)))+
  scale_fill_viridis_d(name="Prey Species")+
  theme(legend.position=c(0.75,0.25),legend.background=element_rect(color="black"))+
  coord_flip()





#Who is the biggest eater of longfin squid
lsquid<-filter(prey19,grepl("LOLIGO",pynam,ignore.case = T))%>%
  mutate(Predator=str_to_title(gsub(" ","\n",pdcomnam)),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))
ggplot(lsquid)+
  geom_bar(aes(Predator,fill=geoarea),position="dodge")+
  scale_y_continuous(name="Number of Stomachs with Any Loligo sp. Squid Identified")

#Who is the biggest eater of herrings
herring<-filter(prey19,pynam%in%c("CLUPEA HARENGUS","ALOSA PSEUDOHARENGUS","ALOSA AESTIVALIS"))%>%
  mutate(Predator=str_to_title(gsub(" ","\n",pdcomnam)),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))
ggplot(herring)+
  geom_bar(aes(Predator,fill=geoarea),position="dodge")+
  scale_y_continuous(name="Number of Stomachs with Herrings Identified")
  #What about strictly juveniles/in tern range (<15 cm)
pylen19<-pylen19%>%
  dplyr::select(-pdsex)%>%
  left_join(uniquePrey19)%>%
  mutate(declon=-declon)%>%
  mutate(Species=str_to_title(gsub(" ","\n",pdcomnam)),
         species=str_to_title(pdcomnam),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))
herringLengths<-filter(pylen19,pynam%in%c("CLUPEA HARENGUS","ALOSA PSEUDOHARENGUS","ALOSA AESTIVALIS"))%>%
  filter(pylen<152)
ggplot(herringLengths)+
  geom_bar(aes(Species,fill=geoarea),position="dodge")+
  scale_y_continuous(name="Number of Stomachs with Juveniles Herrings Identified")

#Who is the biggest eater of hakes
hake<-filter(prey19,pynam%in%c("MERLUCCIUS BILINEARIS","UROPHYCIS CHUSS","UROPHYCIS TENUIS"))%>%
  mutate(Predator=str_to_title(gsub(" ","\n",pdcomnam)),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))
ggplot(hake)+
  geom_bar(aes(Predator,fill=geoarea),position="dodge")+
  scale_y_continuous(name="Number of Stomachs with Hake Identified")
  #What about strictly juveniles/in tern range (<15 cm)
hakeLengths<-filter(pylen19,pynam%in%c("MERLUCCIUS BILINEARIS","UROPHYCIS CHUSS","UROPHYCIS TENUIS"))%>%
  filter(pylen<152)
ggplot(hakeLengths)+
  geom_bar(aes(Species,fill=geoarea),position="dodge")+
  scale_y_continuous(name="Number of Stomachs with Juveniles Hakes Identified")


# Diet composition --------------------------------------------------------

#### Percent Empty ####
#Comparing species
prey19%>%
  mutate(ID=paste(svspp,cruise6,station,pdsex,pdid,pdlen,pdwgt,sep="-"))%>%
  group_by(pdcomnam)%>%
  mutate(IDs=n_distinct(ID))%>%
  group_by(gensci,pdcomnam,IDs)%>%
  reframe(N_empty=n_distinct(ID),
            p_empty=(N_empty/IDs)*100)%>%filter(gensci=="EMPTY")%>%unique()%>%ungroup()%>%
  mutate(pdcomnam=fct_reorder(as.factor(str_to_title(gsub(" ","\n",pdcomnam))),p_empty),
         p_nempty=100-p_empty)%>%
  pivot_longer(cols=c("p_empty","p_nempty"))%>%
  ggplot(aes(pdcomnam,value,value,fill=name))+
  geom_col(color="black")+
  scale_fill_brewer(palette="Set1",name="",labels=c("Empty","Not Empty"))+
  scale_x_discrete(name="Predator Species")+
  scale_y_continuous(name="Percent (%)",expand=expansion(add=1))

#comparing regions
prey19%>%
  mutate(ID=paste(svspp,cruise6,station,pdsex,pdid,pdlen,pdwgt,sep="-"))%>%
  group_by(geoarea)%>%
  mutate(IDs=n_distinct(ID))%>%
  group_by(gensci,geoarea,IDs)%>%
  reframe(N_empty=n_distinct(ID),
            p_empty=(N_empty/IDs)*100)%>%filter(gensci=="EMPTY")%>%unique()%>%ungroup()%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         p_nempty=100-p_empty)%>%
  pivot_longer(cols=c("p_empty","p_nempty"))%>%
  ggplot(aes(geoarea,value,value,fill=name))+
  geom_col(color="black")+
  scale_fill_brewer(palette="Set1",name="",labels=c("Empty","Not Empty"))+
  scale_x_discrete(name="Predator Species")+
  scale_y_continuous(name="Percent (%)",expand=expansion(add=1))


byGeo_empty<-prey19%>%
  mutate(ID=paste(svspp,cruise6,station,pdsex,pdid,pdlen,pdwgt,sep="-"))%>%
  group_by(pdcomnam,geoarea)%>%
  mutate(IDs=n_distinct(ID))%>%
  group_by(gensci,pdcomnam,geoarea,IDs)%>%
  reframe(N_empty=n_distinct(ID),
            p_empty=(N_empty/IDs)*100)%>%filter(gensci=="EMPTY")%>%unique()%>%ungroup()%>%
  mutate(pdcomnam=fct_reorder(as.factor(str_to_title(gsub(" ","\n",pdcomnam))),p_empty),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         p_nempty=100-p_empty)%>%
  pivot_longer(cols=c("p_empty","p_nempty"))

#Comparing regions with the species as replicates (for variances)
byGeo_empty%>%
  filter(name=="p_empty")%>%
  group_by(geoarea)%>%
  summarise(IDs=sum(IDs),
            meanEmpty=mean(value),se=sd(value)/sqrt(n()),
            meanNEmpty=100-meanEmpty,lower=meanNEmpty-1.96*se,upper=meanNEmpty+1.96*se)%>%
  pivot_longer(cols=c(meanEmpty,meanNEmpty))%>%
  ggplot()+
  geom_col(aes(geoarea,value,fill=name))+
  geom_errorbar(aes(geoarea,ymin=lower,ymax=upper),width=0.25)+
  geom_label(aes(geoarea,90,label=paste0("N=",IDs)),
             fill="white",color="black",alpha=0.85)+
  scale_fill_brewer(palette="Set1",name="",labels=c("Empty","Not Empty"))+
  scale_x_discrete(name="Geographical Area")+
  scale_y_continuous(name="Percent (%)",expand=expansion(add=1))


ggplot(byGeo_empty,aes(geoarea,value,fill=name))+
  geom_col(color="black")+
  geom_label(data=filter(byGeo_empty,name=="p_nempty"),aes(label=paste0("N=",IDs)),
             fill="white",color="black",alpha=0.85)+
  scale_fill_brewer(palette="Set1",name="",labels=c("Empty","Not Empty"))+
  scale_x_discrete(name="Geographic Region")+
  scale_y_continuous(name="Percent (%)",expand=expansion(add=1))+
  facet_wrap(~pdcomnam)

#Over the years for each species
s<-prey19%>%
  mutate(ID=paste(svspp,cruise6,station,pdsex,pdid,pdlen,pdwgt,sep="-"))%>%
  group_by(pdcomnam)%>%
  mutate(total=n_distinct(ID))%>%
  group_by(gensci,pdcomnam,total)%>%
  mutate(total_empty=n_distinct(ID),
         total_p_empty=(total_empty/total)*100)%>%ungroup()%>%
  group_by(pdcomnam,year)%>%
  mutate(IDs=n_distinct(ID))%>%
  group_by(gensci,pdcomnam,year,IDs,total,total_empty,total_p_empty)%>%
  reframe(N_empty=n_distinct(ID),
            p_empty=(N_empty/IDs)*100)%>%filter(gensci=="EMPTY")%>%unique()%>%ungroup()%>%
  mutate(f.pdcomnam=fct_reorder(as.factor(str_to_title(pdcomnam)),p_empty))%>%
  ggplot(aes(year,p_empty,color=pdcomnam,fill=pdcomnam))+
  geom_abline(aes(slope=0,intercept=total_p_empty),color="firebrick2",size=1.5)+
  geom_line(size=2)+
  geom_point(size=3,shape=21,stroke=1.5,color="black")+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name="Empty Stomachs (%)",limits=c(0,100))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,l=60))+
  scale_color_manual(values=colorspace::rainbow_hcl(9,c=100,l=60))+
  theme(legend.position = "none")+
  facet_wrap(~f.pdcomnam,nrow=9)
 

#Over the years for each area
g<-prey19%>%
  mutate(ID=paste(svspp,cruise6,station,pdsex,pdid,pdlen,pdwgt,sep="-"))%>%
  group_by(geoarea)%>%
  mutate(total=n_distinct(ID))%>%
  group_by(gensci,geoarea,total)%>%
  mutate(total_empty=n_distinct(ID),
         total_p_empty=(total_empty/total)*100)%>%ungroup()%>%
  group_by(geoarea,year)%>%
  mutate(IDs=n_distinct(ID))%>%
  group_by(gensci,geoarea,year,IDs,total,total_empty,total_p_empty)%>%
  summarise(N_empty=n_distinct(ID),
            p_empty=(N_empty/IDs)*100)%>%filter(gensci=="EMPTY")%>%unique()%>%ungroup()%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))%>%
  ggplot(aes(year,p_empty,shape=geoarea))+
  geom_abline(aes(slope=0,intercept=total_p_empty),color="firebrick2",size=1.5)+
  geom_line(color="black",size=2)+
  geom_point(fill="black",size=3,stroke=1.5)+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name="Percent Empty Stomachs (%)",limits=c(0,100))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_shape_manual(values=c(21:25),guide="none")+
  facet_wrap(~geoarea,nrow=1)

g
#Over the years for each species in each area
a<-prey19%>%
  mutate(ID=paste(svspp,cruise6,station,pdsex,pdid,pdlen,pdwgt,sep="-"))%>%
  group_by(pdcomnam,geoarea)%>%
  mutate(total=n_distinct(ID))%>%
  group_by(gensci,pdcomnam,geoarea,total)%>%
  mutate(total_empty=n_distinct(ID),
         total_p_empty=(total_empty/total)*100)%>%ungroup()%>%
  group_by(pdcomnam,geoarea,year)%>%
  mutate(IDs=n_distinct(ID))%>%
  group_by(gensci,pdcomnam,geoarea,year,IDs,total,total_empty,total_p_empty)%>%
  reframe(N_empty=n_distinct(ID),
            p_empty=(N_empty/IDs)*100)%>%filter(gensci=="EMPTY")%>%unique()%>%ungroup()%>%
  mutate(f.pdcomnam=fct_reorder(as.factor(str_to_title(pdcomnam)),p_empty),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))%>%
  ggplot(aes(year,p_empty,color=pdcomnam,fill=pdcomnam,shape=geoarea))+
  geom_abline(aes(slope=0,intercept=total_p_empty),color="firebrick2",size=1.5)+
  geom_line(size=2,show.legend = F)+
  geom_point(size=3,stroke=1.25,color="black")+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name="Percent Empty Stomachs (%)",limits=c(0,100))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,l=50))+
  scale_color_manual(values=colorspace::rainbow_hcl(9,c=100,l=50))+
  scale_shape_manual(values=c(21:25))+
  guides(fill="none",shape="none")+
  theme(plot.margin=unit(c(12.5,15,15,12.5),"pt"))+
  facet_wrap(~f.pdcomnam+geoarea,nrow=9,labeller = label_wrap_gen(multi_line=FALSE))

a
#Put them all together
grid.arrange(a,s,g,nrow=2,layout_matrix=rbind(c(1,1,1,2),
                                              c(1,1,1,2),
                                              c(1,1,1,2),
                                              c(1,1,1,2),
                                              c(3,3,3,NA)))




#Arrington et al. appendix data
arrington<-read.csv("../arringtonetal.2002_appendix.csv")%>%
  mutate(percentEmpty=as.numeric(gsub("%","",Percentage.with.empty.stomachs)))

myArrington<-filter(arrington,Species %in% myspecies)

ggplot(arrington,aes(Number.of.individuals.analyzed,percentEmpty,fill=Collection.location))+
  geom_point(size=3,shape=21)+
  ylim(0,100)+xlim(0,1000)


#Where did diets get examined. Up until 2004, they were probably all examined at sea
  #Up until 1991 though, they weren't flagged in any way at all
table(prey19$fhdat,prey19$year)

#### Compositional Analyses #### 
speciesPreyComp<-prey19%>%
  mutate(Species=str_to_title(gsub(" ","\n",pdcomnam)),
         gensci=str_to_title(as.character(gensci)))%>%
  filter(gensci %notin% c("Empty","Unobs","Miscellaneous"))%>% #Improves the visualization and they're not really food
  group_by(Species)%>%
  mutate(totalCount=sum(pynum,na.rm=T),totalWeight=sum(pyamtw),totalVol=sum(pyamtv))%>%
  group_by(Species,gensci)%>%
  mutate(propCount=sum(pynum,na.rm=T)/totalCount,propWeight=sum(pyamtw)/totalWeight,propVol=sum(pyamtv)/totalVol)%>%
  dplyr::select(Species,pdscinam,gensci,totalCount,propCount,totalWeight,propWeight,totalVol,propVol)%>%
  distinct()

#Count Prop
ggplot(speciesPreyComp,aes(Species,propCount,fill=gensci))+
  geom_col(position=position_stack(),color="black")+
  scale_fill_manual(values=colorspace::rainbow_hcl(21,c=100),name="Prey Category")+
  scale_y_continuous(expand=expansion(mult=c(0.001,0.001)),name="Proportion of Diet Count")

#Weight Prop
ggplot(speciesPreyComp,aes(Species,propWeight,fill=gensci))+
  geom_col(position=position_stack(),color="black")+
  scale_fill_manual(values=colorspace::rainbow_hcl(21,c=100),name="Prey Category")+
  scale_y_continuous(expand=expansion(mult=c(0.001,0.001)),name="Proportion of Diet Weight")

#Volume Prop
ggplot(speciesPreyComp,aes(Species,propVol,fill=gensci))+
  geom_col(position=position_stack(),color="black")+
  scale_fill_manual(values=colorspace::rainbow_hcl(21,c=100),name="Prey Category")+
  scale_y_continuous(expand=expansion(mult=c(0.001,0.001)),name="Proportion of Diet Volume")

#Putting the three together
pattern<-unique(speciesPreyComp$gensci)%>%
  bind_cols(pat=rep(c("Hatch","NotHatch"),11)[1:n_distinct(prey19$gensci)])
colnames(pattern)<-c("gensci","pat")
fullName<-as_labeller(c("propCount"="Proportion of Count",
                        "propWeight"="Proportion of Weight",
                        "propVol"="Proportion of Volume"))

#speciesPreyComp%>%
#  pivot_longer(cols=c(propCount,propWeight,propVol))%>%
#  left_join(pattern)%>%
#  ggplot()+
#  geom_col_pattern(color="black",aes(Species,value,fill=gensci,pattern=pat))+
#  scale_fill_manual(values=colorspace::rainbow_hcl(n_distinct(prey19$gensci),c=100,l=60),name="Prey Category")+
#  scale_y_continuous(expand=expansion(mult=c(0.001,0.001)),name="Proportion of Diet")+
#  scale_pattern_manual(values = c(Hatch = "stripe", NotHatch = "none"),guide="none") +
#  theme(legend.position="top",panel.spacing.y = unit(25,"pt"))+
#  guides(fill=guide_legend(nrow=5,override.aes=list(pattern=rep(c("stripe","none"),10)[1:n_distinct(prey19$gensci)],
#                                                    pattern_spacing=0.01)))+
#  facet_grid(rows=vars(name),labeller=fullName)


speciesPreyComp%>%
  pivot_longer(cols=c(propCount,propWeight,propVol))%>%
  left_join(pattern)%>%
  ungroup()%>%mutate(gensci=fct_reorder(gensci,value))%>%
  ggplot()+
  geom_col(color="black",aes(Species,value,fill=gensci))+
  scale_fill_manual(values=colorspace::rainbow_hcl(n_distinct(prey19$gensci),c=100,l=60),name="General\nPrey Category")+
  scale_y_continuous(expand=expansion(mult=c(0.001,0.001)),name="Proportion of Diet")+
  theme(legend.position="top",panel.spacing.y = unit(25,"pt"))+
  facet_wrap(~name,ncol=1,labeller=fullName)


#Diet Richness



# Comparing GeoRegions ----------------------------------------------------
B<-function(p) {
  a<-apply(p,2,sum)/sum(p)
  b<-(1/(sum(a^2))-1)/(length(a)-1)
  return(b)
}
#Bootstrapping function
myboot <- function(data, num) {
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data), replace=T),])
  levin <- sapply(resamples, B)
  std.dev <- sd(levin)
  list(std.dev=std.dev, levins=mean(levin))   
}
#Usage: myboot(df,999)
myboot_reps <- function(data, num) {
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data), replace=T),])
  levin <- sapply(resamples, B)
  std.dev <- sd(levin)
  return(levin)   
}


#Species abundance
prey19%>%
  filter(pynam!="EMPTY")%>%
  mutate(species=str_to_title(gsub(" ","\n",pdcomnam)),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))%>%
  group_by(species,.drop=F)%>%
  summarise(nPrey=n_distinct(collsci))%>%
  left_join(speciesLevin_coll%>%mutate(species=str_to_title(gsub(" ","\n",species))))%>%
  ggplot(aes(fct_reorder(species,levinMean),nPrey,fill=fct_reorder(species,levinMean)))+
  geom_col(position=position_dodge(),color="black")+
  scale_y_continuous(expand=expansion(add=c(0.5,20)),name="Number of Distinct Prey")+
  scale_x_discrete(name="Predator Species")+
  scale_fill_viridis_d(guide="none")


#Levin's B
### Matrix of the diets for each individual (species as columns and individuals each row)
ids<-c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat")
preyMat<-pivot_wider(prey19,id_cols=all_of(ids),
                     names_from=pynam,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
preyMat[is.na(preyMat)]<-0

#For each of the different species
speciesLevin_geo<-data.frame(species=character(),geoarea=character(),N=numeric(),levinSD=numeric(),levinMean=numeric())
species<-unique(as.character(preyMat$pdcomnam))
geos<-unique(as.character(preyMat$geoarea))
n=1

for (s in 1:length(species)) {
  for (g in 1:length(geos)) {
    subMat<-preyMat[which(preyMat$pdcomnam==species[s]&preyMat$geoarea==geos[g]),(length(ids)+1):ncol(preyMat)]
    speciesLevin_geo[n,]<-c(species[s],geos[g],nrow(subMat),
                            myboot(subMat[sample(seq(1:nrow(subMat)),
                                                 ifelse(nrow(subMat)<1000,nrow(subMat),1000),
                                                 replace=F),],
                                   999))
    n=n+1
  }
}

#Plotting
speciesLevin_geo%>%
  mutate(N=as.numeric(N),levinMean=as.numeric(levinMean),
         species=str_to_title(gsub(" ","\n",species)),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))%>%
  ggplot(aes(species,levinMean,fill=geoarea))+
  geom_col(color="black",size=1,position = position_dodge(0.9))+
  geom_errorbar(aes(species,ymin=levinMean-levinSD,ymax=levinMean+levinSD),
                width=0.25,position = position_dodge(0.9))+
  scale_fill_manual(values=colorspace::rainbow_hcl(5,c=100,l=60),name="Geographic Area")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Levin's Niche Breadth")+
  scale_x_discrete(name="Predator Species")+
  theme(legend.position=c(0.8,0.8),legend.background = element_rect(color="black"))





# Diet Metrics (Relative Consumption) -------------------------------------

#How many/percentage of diets can this be done for?
nrow(filter(uniquePrey19,!is.na(pdwgt))) #Number
nrow(filter(uniquePrey19,!is.na(pdwgt)))/nrow(uniquePrey19)*100 #percentage

#Relative consumption for each of the species
uniquePrey19%>%
  mutate(pdcomnam=str_to_title(gsub(" ","\n",pdcomnam)))%>%
  group_by(pdcomnam)%>%
  mutate(pNA=sum(is.na(pdwgt))/n()*100)%>%
  filter(pdwgt>0)%>%
  mutate(relConsump=pdgutw/pdwgt,
         relConsump_mean=mean(pdgutw/pdwgt,na.rm=T),
         relConsump_min=min(pdgutw/pdwgt,na.rm=T),
         relConsump_max=max(pdgutw/pdwgt,na.rm=T),
         relConsump_sd=sd(pdgutw/pdwgt,na.rm=T))%>%
  filter(relConsump==relConsump_max)%>% #Keep the mass of the fish with the max relative consumption
  mutate(pdwgt=mean(pdwgt))%>%
  dplyr::select(pdcomnam,pdwgt,relConsump_mean:relConsump_sd)%>%unique()%>%
  ggplot()+
  geom_col(aes(fct_reorder(pdcomnam,relConsump_mean),relConsump_mean,fill=pdcomnam),color="black",show.legend=F)+
  geom_errorbar(aes(fct_reorder(pdcomnam,relConsump_mean),
                    ymin=relConsump_mean,ymax=relConsump_max),width=0.2)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,l=60))+
  scale_x_discrete(name="Predator Species")+
  scale_y_continuous(name="Mean Relative Consumption",expand=expansion(add=c(0.001,0.01)))
  
uniquePrey19%>%
  mutate(pdcomnam=str_to_title(gsub(" ","\n",pdcomnam)))%>%
  filter(pdwgt>0 & pdgutw>0)%>%
  mutate(relConsump=pdgutw/pdwgt)%>%
  bind_rows(data.frame(cruise6=NA,station=NA,svspp=NA,pdsex=NA,pdid=NA,pdcomnam="Summer\nFlounder",
                       pdscinam=NA,pdlen=NA,pdwgt=NA,sizecat=NA,pdgutw=NA,pdgutv=NA,declat=NA,declon=NA,
                       month=NA,day=NA,year=NA,season="SUMMER",geoarea=NA,relConsump=0))%>%
  mutate(Season=factor(str_to_title(season),levels=c("Winter","Spring","Summer","Fall")))%>%
  ggplot()+
  geom_boxplot(aes(fct_reorder(pdcomnam,relConsump),relConsump),outlier.shape=NA,
               size=1.5,color="black",alpha=1)+
  geom_jitter(aes(fct_reorder(pdcomnam,relConsump),relConsump),
              shape=21,alpha=0.8,width=0.1)+
  scale_fill_manual(values=seasonPal2)+ #must be ordered w,sp,su,f
  scale_x_discrete(name="Predator Species")+
  scale_y_log10(name="Relative Consumption (g/g)",breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),
                labels=c(0.00001,0.0001,0.001,0.01,0.1,1))




# Diet Metrics (Breadth) --------------------------------------------------

### Functions to do so

#Niche breadth function
B<-function(p) {
  a<-apply(p,2,sum)/sum(p)
  b<-(1/(sum(a^2))-1)/(length(a)-1)
  return(b)
}
#Bootstrapping function
myboot <- function(data, num) {
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data), replace=T),])
  levin <- sapply(resamples, B)
  std.dev <- sd(levin)
  list(std.dev=std.dev, levins=mean(levin))   
}
#Usage: myboot(df,999)
myboot_reps <- function(data, num) {
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data), replace=T),])
  levin <- sapply(resamples, B)
  std.dev <- sd(levin)
  return(levin)   
}
 

### Matrix of the diets for each individual (species as columns and individuals each row)
preyMat<-pivot_wider(prey19,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                     names_from=pynam,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
preyMat[is.na(preyMat)]<-0
#Number of prey categories
ncol(preyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))


#Using the broadest level of resolution--general category
genpreyMat<-prey19%>%
  mutate(gensci=ifelse(gensci=="",as.character(pynam),as.character(gensci)))
genpreyMat<-pivot_wider(genpreyMat,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                     names_from=gensci,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
genpreyMat[is.na(genpreyMat)]<-0
#Number of prey categories
ncol(genpreyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))

#Levin's Niche Breadth
B(genpreyMat[,15:ncol(genpreyMat)])
#Bootstrapping, too much for the computer unless low reps
myboot(genpreyMat[,15:ncol(genpreyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesLevin_gen<-data.frame(species=character(),N=numeric(),levinSD=numeric(),levinMean=numeric())
species<-unique(as.character(genpreyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-genpreyMat[which(genpreyMat$pdcomnam==species[s]),15:ncol(genpreyMat)] #15 is the first prey item column
  speciesLevin_gen[s,]<-c(species[s],nrow(subMat),
                      myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesLevin_gen%>%
  mutate(N=as.numeric(N),levinMean=as.numeric(levinMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(species,levinMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=levinMean-levinSD,ymax=levinMean+levinSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Levin's Niche Breadth")+
  scale_x_discrete(name="Predator Species")




#More detailed prey--analytical

analpreyMat<-prey19%>%
  mutate(analsci=ifelse(analsci=="",as.character(pynam),as.character(analsci)))
analpreyMat<-pivot_wider(analpreyMat,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                        names_from=analsci,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
analpreyMat[is.na(analpreyMat)]<-0
#Number of prey categories
ncol(analpreyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))

#Levin's Niche Breadth
B(analpreyMat[,15:ncol(analpreyMat)])
#Bootstrapping, maybe too much for the computer unless low reps (test)
myboot(analpreyMat[,15:ncol(analpreyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesLevin_anal<-data.frame(species=character(),N=numeric(),levinSD=numeric(),levinMean=numeric())
species<-unique(as.character(analpreyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-analpreyMat[which(analpreyMat$pdcomnam==species[s]),15:ncol(analpreyMat)] #15 is the first prey item column
  speciesLevin_anal[s,]<-c(species[s],nrow(subMat),
                           myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesLevin_anal%>%
  mutate(N=as.numeric(N),levinMean=as.numeric(levinMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(species,levinMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=levinMean-levinSD,ymax=levinMean+levinSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Levin's Niche Breadth")+
  scale_x_discrete(name="Predator Species")



#Next most detailed prey--collection

#collpreyMat<-prey19%>%
#  mutate(collsci=ifelse(collsci=="",as.character(pynam),as.character(collsci)))
collpreyMat<-pivot_wider(prey19,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                         names_from=collsci,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
collpreyMat[is.na(collpreyMat)]<-0
#Number of prey categories
envCols<-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))
ncol(collpreyMat)-envCols

#Levin's Niche Breadth
B(collpreyMat[,(envCols+1):ncol(collpreyMat)])
#Bootstrapping, maybe too much for the computer unless low reps (test)
myboot(collpreyMat[,(envCols+1):ncol(collpreyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesLevin_coll<-data.frame(species=character(),N=numeric(),levinSD=numeric(),levinMean=numeric())
species<-unique(as.character(collpreyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-collpreyMat[which(collpreyMat$pdcomnam==species[s]),(envCols+1):ncol(collpreyMat)] #15 is the first prey item column
  speciesLevin_coll[s,]<-c(species[s],nrow(subMat),
                           myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesLevin_coll%>%
  mutate(N=as.numeric(N),levinMean=as.numeric(levinMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(fct_reorder(species,levinMean),levinMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(fct_reorder(species,levinMean),ymin=levinMean-levinSD,ymax=levinMean+levinSD),width=0.2)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,l=60),guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Levin's Niche Breadth")+
  scale_x_discrete(name="Predator Species")





#Most detailed--actual prey naming
ncol(preyMat)-envCols

#Levin's Niche Breadth
B(preyMat[,15:ncol(preyMat)])
#Bootstrapping, maybe too much for the computer unless low reps (test)
myboot(preyMat[,15:ncol(preyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesLevin_py<-data.frame(species=character(),N=numeric(),levinSD=numeric(),levinMean=numeric())
species<-unique(as.character(preyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-preyMat[which(preyMat$pdcomnam==species[s]),15:ncol(preyMat)] #15 is the first prey item column
  speciesLevin_py[s,]<-c(species[s],nrow(subMat),
                         myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesLevin_py%>%
  mutate(N=as.numeric(N),levinMean=as.numeric(levinMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(species,levinMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=levinMean-levinSD,ymax=levinMean+levinSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Levin's Niche Breadth")+
  scale_x_discrete(name="Predator Species")


#Comparing them all
speciesLevin_gen$resol<-"General Taxonomic Resolution"
speciesLevin_anal$resol<-"Analytical Taxonomic Resolution"
speciesLevin_coll$resol<-"Collection Taxonomic Resolution"
bind_rows(speciesLevin_gen,speciesLevin_anal,speciesLevin_coll)%>%
  mutate(species=str_to_title(gsub(" ","\n",species)),
         resol=factor(resol,levels=c("General Taxonomic Resolution","Analytical Taxonomic Resolution","Collection Taxonomic Resolution")))%>%
  ggplot()+
  geom_col(aes(species,levinMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=levinMean-levinSD,ymax=levinMean+levinSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Levin's Niche Breadth")+
  scale_x_discrete(name="Predator Species")+
  facet_wrap(resol~.,nrow=3)






### Using General (faster, easier, higher) over time



levinOverTime<-genpreyMat%>%
  filter(year>1)%>%
  group_by(pdcomnam,year)%>%
  summarise(N=n(),
            levinB_mean=myboot(as.matrix(across(ARTHROPODA:STERNOPTYCHIDAE)),1000)$levins,
            levinB_sd=myboot(as.matrix(across(ARTHROPODA:STERNOPTYCHIDAE)),1000)$std.dev)

levinOverTime2<-levinOverTime%>%
  filter(!(is.nan(levinB_mean)|levinB_mean==0))%>% #Dropping those times when there was only 1 fish so it couldn't calculate B
  mutate(Species=str_to_title(pdcomnam),
         levinB_lower95=levinB_mean-1.96*levinB_sd,
         levinB_upper95=levinB_mean+1.96*levinB_sd)%>%
  group_by(Species)%>%
  mutate(g = c(0, cumsum(diff(year) > 1)))%>%
  group_by(Species,g)%>%
  mutate(n=n())

ggplot(levinOverTime2)+
  geom_ribbon(aes(year,ymin=levinB_lower95,ymax=levinB_upper95,group=g),alpha=0.3,fill="grey50")+
  geom_linerange(data=filter(levinOverTime2,n==1),
                 aes(year,ymin=levinB_lower95,ymax=levinB_upper95,group=g),size=2,alpha=0.3,color="grey50")+
  geom_line(aes(year,levinB_mean,color=Species,group=g),size=1)+
  geom_point(aes(year,levinB_mean,fill=Species,size=n),color="black",shape=21)+
  scale_color_manual(values=colorspace::rainbow_hcl(9,c=100,end=280),guide="none")+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,end=280),guide="none")+
  scale_size_continuous(range=c(1,7))+
  scale_x_continuous(expand=expansion(0.1),name="Year")+
  scale_y_continuous(expand=expansion(mult=c(0.001,0.1)),name="Levin's Niche Breadth")+
  facet_wrap(~Species)+theme_bw(base_size=30)+
  theme(panel.spacing=unit(20,"pt"))



levinOverTimeS<-genpreyMat%>%
  filter(year>1)%>%
  group_by(pdcomnam,year,season)%>%
  summarise(levinB=B(as.matrix(across(ARTHROPODA:STERNOPTYCHIDAE))))

levinOverTimeS%>%
  filter(!(is.nan(levinB)|levinB==0))%>% #Dropping those times when there was only 1 fish so it could calculate B
  mutate(Species=str_to_title(pdcomnam),
         Season=str_to_title(season),
         Season=factor(Season,levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(Species,Season)%>%
  mutate(g = c(0, cumsum(diff(year) > 1)),
         sg=paste0(Season,g)) %>%
  ggplot()+
  geom_line(aes(year,levinB,color=Season,group=sg))+
  geom_point(aes(year,levinB,fill=Season),size=2,color="black",shape=21)+
  scale_color_manual(values=seasonPal2)+scale_fill_manual(values = seasonPal2)+
  scale_x_continuous(expand=expansion(0.1),name="Year")+
  scale_y_continuous(expand=expansion(mult=c(0.01,0.1)),name="Levin's Niche Breadth")+
  facet_wrap(~Species)+
  theme(panel.spacing=unit(18,"pt"))



### Shannon-Weiner 
H <- function(p) {
  a<-apply(p,2,sum)/sum(p)
  b<--1*sum(a*log(a))/log(length(a))
  return(b)
}
#Bootstrapping function
mybootH <- function(data, num) {
  resamples <- lapply(1:num, function(i) data[sample(1:nrow(data), replace=T),])
  shannon <- sapply(resamples, H)
  list(std.dev=sd(shannon,na.rm=T), shannons=mean(shannon,na.rm=T))   
}
   


### Matrix of the diets for each individual (species as columns and individuals each row)
preyMat<-pivot_wider(prey19,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                     names_from=pynam,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
preyMat[is.na(preyMat)]<-0
#Number of prey categories
ncol(preyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))


#Using the broadest level of resolution--general category
genpreyMat<-prey19%>%
  mutate(gensci=ifelse(gensci=="",as.character(pynam),as.character(gensci)))
genpreyMat<-pivot_wider(genpreyMat,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                        names_from=gensci,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
genpreyMat[is.na(genpreyMat)]<-0
#Number of prey categories
ncol(genpreyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))

#Levin's Niche Breadth
H(genpreyMat[,15:ncol(genpreyMat)])
#Bootstrapping, too much for the computer unless low reps
mybootH(genpreyMat[,15:ncol(genpreyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesSW_gen<-data.frame(species=character(),N=numeric(),shannonSD=numeric(),shannonMean=numeric())
species<-unique(as.character(genpreyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-genpreyMat[which(genpreyMat$pdcomnam==species[s]),15:ncol(genpreyMat)] #15 is the first prey item column
  speciesSW_gen[s,]<-c(species[s],nrow(subMat),
                          myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesSW_gen%>%
  mutate(N=as.numeric(N),shannonMean=as.numeric(shannonMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(species,shannonMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=shannonMean-shannonSD,ymax=shannonMean+shannonSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Shannon-Weiner Niche Breadth")+
  scale_x_discrete(name="Predator Species")




#More detailed prey--analytical

analpreyMat<-prey19%>%
  mutate(analsci=ifelse(analsci=="",as.character(pynam),as.character(analsci)))
analpreyMat<-pivot_wider(analpreyMat,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                         names_from=analsci,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
analpreyMat[is.na(analpreyMat)]<-0
#Number of prey categories
ncol(analpreyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))

#SW's Niche Breadth
H(analpreyMat[,15:ncol(analpreyMat)])
#Bootstrapping, maybe too much for the computer unless low reps (test)
mybootH(analpreyMat[,15:ncol(analpreyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesSW_anal<-data.frame(species=character(),N=numeric(),shannonSD=numeric(),shannonMean=numeric())
species<-unique(as.character(analpreyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-analpreyMat[which(analpreyMat$pdcomnam==species[s]),15:ncol(analpreyMat)] #15 is the first prey item column
  speciesSW_anal[s,]<-c(species[s],nrow(subMat),
                           myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesSW_anal%>%
  mutate(N=as.numeric(N),shannonMean=as.numeric(shannonMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(species,shannonMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=shannonMean-shannonSD,ymax=shannonMean+shannonSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean SW's Niche Hreadth")+
  scale_x_discrete(name="Predator Species")



#Next most detailed prey--collection

#collpreyMat<-prey19%>%
#  mutate(collsci=ifelse(collsci=="",as.character(pynam),as.character(collsci)))
collpreyMat<-pivot_wider(prey19,id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                         names_from=collsci,values_from=pyamtw,values_fn=sum)%>%
  dplyr::select(-EMPTY)
collpreyMat[is.na(collpreyMat)]<-0
#Number of prey categories
ncol(collpreyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))

#SW's Niche Breadth
B(collpreyMat[,15:ncol(collpreyMat)])
#Bootstrapping, maybe too much for the computer unless low reps (test)
mybootH(collpreyMat[,15:ncol(collpreyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesSW_coll<-data.frame(species=character(),N=numeric(),shannonSD=numeric(),shannonMean=numeric())
species<-unique(as.character(collpreyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-collpreyMat[which(collpreyMat$pdcomnam==species[s]),15:ncol(collpreyMat)] #15 is the first prey item column
  speciesSW_coll[s,]<-c(species[s],nrow(subMat),
                           myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesSW_coll%>%
  mutate(N=as.numeric(N),shannonMean=as.numeric(shannonMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(species,shannonMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=shannonMean-shannonSD,ymax=shannonMean+shannonSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Shannon-Weiner Niche Breadth")+
  scale_x_discrete(name="Predator Species")


#Most detailed--actual prey naming
ncol(preyMat)-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))

#SW's Niche Breadth
H(preyMat[,15:ncol(preyMat)])
#Bootstrapping, maybe too much for the computer unless low reps (test)
mybootH(preyMat[,15:ncol(preyMat)],10)


### Calculating niche breadth for each species (in total)

#For each of the different species
speciesSW_py<-data.frame(species=character(),N=numeric(),shannonSD=numeric(),shannonMean=numeric())
species<-unique(as.character(preyMat$pdcomnam))

for (s in 1:length(species)) {
  subMat<-preyMat[which(preyMat$pdcomnam==species[s]),15:ncol(preyMat)] #15 is the first prey item column
  speciesSW_py[s,]<-c(species[s],nrow(subMat),
                         myboot(subMat[sample(seq(1:nrow(subMat)),1000,replace=F),],999))
}

#Plotting
speciesSW_py%>%
  mutate(N=as.numeric(N),shannonMean=as.numeric(shannonMean),
         species=str_to_title(gsub(" ","\n",species)))%>%
  ggplot()+
  geom_col(aes(species,shannonMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=shannonMean-shannonSD,ymax=shannonMean+shannonSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean SW's Niche Breadth")+
  scale_x_discrete(name="Predator Species")


#Comparing them all
speciesSW_gen$resol<-"General Taxonomic Resolution"
speciesSW_anal$resol<-"Analytical Taxonomic Resolution"
speciesSW_coll$resol<-"Collection Taxonomic Resolution"
bind_rows(speciesSW_gen,speciesSW_anal,speciesSW_coll)%>%
  mutate(species=str_to_title(gsub(" ","\n",species)),
         resol=factor(resol,levels=c("General Taxonomic Resolution","Analytical Taxonomic Resolution","Collection Taxonomic Resolution")))%>%
  ggplot()+
  geom_col(aes(species,shannonMean,fill=species),color="black",size=1)+
  geom_errorbar(aes(species,ymin=shannonMean-shannonSD,ymax=shannonMean+shannonSD),width=0.2)+
  scale_fill_viridis_d(guide="none")+
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),name="Mean Shannon-Weiner Niche Breadth")+
  scale_x_discrete(name="Predator Species")+
  facet_wrap(resol~.,nrow=3)






### Using General (faster, easier, higher) over time


shannonOverTime<-genpreyMat%>%
  filter(year>1)%>%
  group_by(pdcomnam,year)%>%
  summarise(N=n(),
            shannonH_mean=H(as.matrix(across(ARTHROPODA:STERNOPTYCHIDAE))))

shannonOverTime2<-shannonOverTime%>%
  filter(!(is.nan(shannonH_mean)|shannonH_mean==0))%>% #Dropping those times when there was only 1 fish so it couldn't calculate B
  mutate(Species=str_to_title(pdcomnam))%>% #Add these back in if you figure out SD
#         shannonH_lower95=shannonH_mean-1.96*shannonH_sd,
#         shannonH_upper95=shannonH_mean+1.96*shannonH_sd)%>%
  group_by(Species)%>%
  mutate(g = c(0, cumsum(diff(year) > 1)))%>%
  group_by(Species,g)%>%
  mutate(n=n())

ggplot(shannonOverTime2)+
  #geom_ribbon(aes(year,ymin=shannonH_lower95,ymax=shannonH_upper95,group=g),alpha=0.3,fill="grey50")+
  #geom_linerange(data=filter(shannonOverTime2,n==1),
  #               aes(year,ymin=shannonH_lower95,ymax=shannonH_upper95,group=g),size=2,alpha=0.3,color="grey50")+
  geom_line(aes(year,shannonH_mean,color=Species,group=g),size=1)+
  geom_point(aes(year,shannonH_mean,fill=Species,size=n),color="black",shape=21)+
  scale_color_manual(values=colorspace::rainbow_hcl(9,c=100,end=280),guide="none")+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,end=280),guide="none")+
  scale_size_continuous(range=c(1,7))+
  scale_x_continuous(expand=expansion(0.1),name="Year")+
  scale_y_continuous(expand=expansion(mult=c(0.001,0.1)),name="Shannon-Weiner Niche Breadth")+
  facet_wrap(~Species)+theme_bw(base_size=30)+
  theme(panel.spacing=unit(20,"pt"))



shannonOverTimeS<-genpreyMat%>%
  filter(year>1)%>%
  group_by(pdcomnam,year,season)%>%
  summarise(shannonH=H(as.matrix(across(ARTHROPODA:STERNOPTYCHIDAE))))

shannonOverTimeS%>%
  filter(!(is.nan(shannonB)|shannonB==0))%>% #Dropping those times when there was only 1 fish so it could calculate B
  mutate(Species=str_to_title(pdcomnam),
         Season=str_to_title(season),
         Season=factor(Season,levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(Species,Season)%>%
  mutate(g = c(0, cumsum(diff(year) > 1)),
         sg=paste0(Season,g)) %>%
  ggplot()+
  geom_line(aes(year,shannonB,color=Season,group=sg))+
  geom_point(aes(year,shannonB,fill=Season),size=2,color="black",shape=21)+
  scale_color_manual(values=seasonPal2)+scale_fill_manual(values = seasonPal2)+
  scale_x_continuous(expand=expansion(0.1),name="Year")+
  scale_y_continuous(expand=expansion(mult=c(0.01,0.1)),name="Levin's Niche Breadth")+
  facet_wrap(~Species)+
  theme(panel.spacing=unit(18,"pt"))











# Diet Metrics (PSI) Prey Similarity? ------------------------------------------------------

#In Bolnick et al. 2002--wrong
P <- function(p) {
  1-0.5*sum(p-1)
}




# Dietary Overlap (Schoener) ----------------------------------------------

#This may be better served in a separate R script at some point, since this could be insular
#But for now
#This is the formula, but I don't have a way to do it fast with this
do<-function (p1,p2) { #p1 is the proportion of a prey item in predator 1, p2 is the proportion of THE SAME prey item in predator 2
  1-0.5*sum(abs(p1-p2))
}


#The full time period
#Creating a proportion matrix for the prey items for my predators
propLong<-prey19%>%
  mutate(species=str_to_title(pdcomnam))%>% #Maybe want common names later, but the myspecies is scientific so speed
  group_by(species,sizecat)%>%
  mutate(total_w=sum(pyamtw),
         total_v=sum(pyamtv))%>%
  group_by(species,sizecat,gl_prey)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()


#Following through to get the "real" overlap matrix
props1<-propLong[,-max(ncol(propLong))] #Cut off the volume prop for cleanliness
colnames(props1)<-c("species1","sizecat1","gl_prey","prop_w1")
props2<-propLong[,-max(ncol(propLong))]
colnames(props2)<-c("species2","sizecat2","gl_prey","prop_w2")

overlap_schoener<-full_join(props1,props2)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat<-pivot_wider(overlap_schoener,id_cols=c(species1,sizecat1),names_from=c(species2,sizecat2),values_from = s_do)
overlap_mat<-as.matrix(overlap_mat[,3:ncol(overlap_mat)])
rownames(overlap_mat)<-gsub(" ","\n",colnames(overlap_mat))
colnames(overlap_mat)<-gsub(" ","\n",colnames(overlap_mat))
#overlap_mat[overlap_mat==1]<-NA
for (i in 1:nrow(overlap_mat)) {
  for (j in 1:ncol(overlap_mat)) {
    overlap_mat[i,j]<-ifelse(j>i,NA,overlap_mat[i,j])
  }
}



#Bootstrapping
propMat<-propLong%>%
  pivot_wider(id_cols=c(species,sizecat),names_from = gl_prey,values_from = prop_w)%>% #Species are rows, prey are columns. Flip id and names if need opposite
  mutate(species=paste(species,sizecat,sep="_"))%>%
  select(-c("Empty","Unobserved","Miscellaneous","sizecat"))
propMat<-select(propMat,species,order(colnames(propMat)))
propMat[is.na(propMat)]<-0

resample<-function(props) {
  newProps<-matrix(nrow=nrow(props),ncol=ncol(props))
  for (i in 1:nrow(props)){
    eat<-props[i,which(props[i,]>0 & colnames(props)!="species")]
    colnames(eat)=sample(colnames(eat))
    neat=sample(props[i,which(props[i,]==0)])
    teat<-cbind(props[i,1],eat,neat)
    teat<-select(teat,species,order(colnames(teat)))
    newProps[i,]<-as.character(teat)
  }
  colnames(newProps)<-colnames(teat)
  newProps<-as.data.frame(newProps)
  newProps[,2:ncol(newProps)]<-sapply(newProps[,2:ncol(newProps)],as.numeric)
  newProps<-pivot_longer(newProps,cols=colnames(newProps[,2:ncol(newProps)]),
                         names_to="gl_prey",values_to="prop_w")
  props1<-newProps
  colnames(props1)<-c("species1","gl_prey","prop_w1")
  props2<-newProps
  colnames(props2)<-c("species2","gl_prey","prop_w2")
  
  overlap_schoener<-full_join(props1,props2)%>%
    mutate(diff=abs(prop_w1-prop_w2))%>%
    group_by(species1,species2)%>%
    summarise(ep=sum(diff),
              s_do=1-0.5*ep)
  
  overlap_mat<-pivot_wider(overlap_schoener,id_cols=species1,names_from=species2,values_from = s_do)
  overlap_mat<-as.matrix(overlap_mat[,2:ncol(overlap_mat)])
  for (i in 1:nrow(overlap_mat)) {
    for (j in 1:ncol(overlap_mat)) {
      overlap_mat[i,j]<-ifelse(j>=i,NA,overlap_mat[i,j])
    }
  }
  overlap_vect<-overlap_mat[!is.na(as.vector(overlap_mat))]
  return(overlap_vect)
}

reps=1
nreps=250
bootDiffs<-numeric()
while (reps<=nreps) {
  bootDiffs<-c(bootDiffs,resample(propMat))
  reps=reps+1
}
hist(bootDiffs)
sigGuild<-quantile(bootDiffs,probs=0.95)

#Visual of the actual matrix with significance indicated
#library(plot.matrix)
#par(mar=c(6,6,5,5.5))
#plot(overlap_mat,axis.row=list(side=2,las=1),col=viridis::viridis(n=100,option="B"),
#     polygon.key = list(border=NA), key=list(),xlab="",ylab="")

#Better visual (ggplot as always)
as.data.frame(overlap_mat)%>%
  mutate(species1=rownames(overlap_mat))%>%
  pivot_longer(cols=colnames(overlap_mat),names_to = "species2", values_to = "s_do")%>%
  mutate(species1=sort(species1,decreasing=T),species2=species2,
         sig=ifelse(s_do>sigGuild,"bold.italic","plain"))%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=s_do))+
  geom_text(aes(species1,species2,label=round(s_do,digits=2),fontface=sig,color=sig),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100],viridis::viridis(100,option="B")[sigGuild*100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5))



#Clustering these out
overlap_clust<-hclust(as.dist(1-overlap_mat),method="average")
plot(overlap_clust)

#Unsure how to better control the look of the dendrogram, probably there's a package out there that looks nice






###############################################
#The first portion, reflecting the same time period as Garrison and Link (1973-1997)
#There is one little skate over 60cm, meaning it is "large" during this time period. Remove it for ease
#Creating a proportion matrix for the prey items for my predators
propLong_gl<-prey19%>%
  filter(year<=1997)%>%
  mutate(species=str_to_title(pdcomnam))%>%
  group_by(species,sizecat)%>%
  mutate(total_w=sum(pyamtw),
         total_v=sum(pyamtv))%>%
  group_by(species,sizecat,gl_prey)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()%>%
  filter(!(species=="Little Skate" & sizecat=="L") )

#In garrison and Link these 9 species have 26 size categories
  #Spiny SML, Little SM, Silver SML, Haddock SML, Pollock SMLX, White SML, Red SML, Summer ML, Yellowtail SML
#In my version of their dataset there are 30 size categories
  #Same plus Little L, Haddock X, Summer SX
  #Little L--1 diet (63 cm)--DROPPED as a single individual
  #Haddock X--74 diets (81-88 cm)
  #Summer S--189 diets (13-20 cm)
  #Summer X--16 diets (71-300 cm)

  
#Following through to get the "real" overlap matrix
props1_gl<-propLong_gl[,-max(ncol(propLong_gl))] #Cut off the volume prop for cleanliness
colnames(props1_gl)<-c("species1","sizecat1","gl_prey","prop_w1")
props2_gl<-propLong_gl[,-max(ncol(propLong_gl))]
colnames(props2_gl)<-c("species2","sizecat2","gl_prey","prop_w2")

overlap_schoener_gl<-full_join(props1_gl,props2_gl)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat_gl<-pivot_wider(overlap_schoener_gl,id_cols=c(species1,sizecat1),names_from=c(species2,sizecat2),values_from = s_do)
overlap_mat_gl<-as.matrix(overlap_mat_gl[,3:ncol(overlap_mat_gl)])
rownames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
colnames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#overlap_mat[overlap_mat==1]<-NA
for (i in 1:nrow(overlap_mat_gl)) {
  for (j in 1:ncol(overlap_mat_gl)) {
    overlap_mat_gl[i,j]<-ifelse(j>i,NA,overlap_mat_gl[i,j])
  }
}



#Bootstrapping
propMat_gl<-propLong_gl%>%
  pivot_wider(id_cols=c(species,sizecat),names_from = gl_prey,values_from = prop_w)%>% #Species are rows, prey are columns. Flip id and names if need opposite
  mutate(species=paste(species,sizecat,sep="_"))%>%
  select(-c("Empty","Unobserved","Miscellaneous","sizecat"))
propMat_gl<-select(propMat_gl,species,order(colnames(propMat_gl)))
propMat_gl[is.na(propMat_gl)]<-0

reps=1
nreps=250
bootDiffs_gl<-numeric()
while (reps<=nreps) {
  bootDiffs_gl<-c(bootDiffs_gl,resample(propMat_gl))
  reps=reps+1
}
hist(bootDiffs_gl)
sigGuild_gl<-quantile(bootDiffs_gl,probs=0.95)

#Visual of the actual matrix with significance indicated
#library(plot.matrix)
#par(mar=c(6,6,5,5.5))
#plot(overlap_mat,axis.row=list(side=2,las=1),col=viridis::viridis(n=100,option="B"),
#     polygon.key = list(border=NA), key=list(),xlab="",ylab="")

#Better visual (ggplot as always)
as.data.frame(overlap_mat_gl)%>%
  mutate(species1=rownames(overlap_mat_gl))%>%
  pivot_longer(cols=colnames(overlap_mat_gl),names_to = "species2", values_to = "s_do")%>%
  mutate(species1=sort(species1,decreasing=T),species2=species2,
         sig=ifelse(s_do>sigGuild,"bold.italic","plain"))%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=s_do))+
  geom_text(aes(species1,species2,label=round(s_do,digits=2),fontface=sig,color=sig),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100],viridis::viridis(100,option="B")[sigGuild*100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  ggtitle("1973-1997")+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5))

#Mean overlap, how similar are the diets in this time period
mean(overlap_mat_gl[which(overlap_mat_gl<1)],na.rm=T)
sd(overlap_mat_gl[which(overlap_mat_gl<1)],na.rm=T)
range(overlap_mat_gl[which(overlap_mat_gl<1)],na.rm=T)
hist(overlap_mat_gl[which(overlap_mat_gl<1)])

#Clustering these out
overlap_clust_gl<-hclust(as.dist(1-overlap_mat_gl),method="average")
#Trying dendextend to pretty up the dendrogram
dend<-as.dendrogram(overlap_clust_gl)

#To save them out, there is some difficulty with the export losing the colorbar
png(file = "../Figures/2021 Garrison and Link Rep/fancy_1997dendrogram.png",   # The directory you want to save the file in
    width = 1250, # The width of the plot in inches
    height = 700) # The height of the plot in inches
par(mar=c(5,1,2,12))

dend %>% 
  set("branches_lwd", 4) %>%
  # Custom labels
  set("labels_cex", 1) %>%
  set("labels_col", value = viridis::viridis(14,end=0.8),h = sigGuild_gl) %>%
  set("branches_k_color", value = viridis::viridis(14,end=0.8), h = sigGuild_gl) %>%
  plot(horiz=TRUE,main="                               1973-1997 Dendrogram",axes=F,xlab="Similarity")
axis(side=1,at=c(0.6,0.5,0.4,0.3,0.2,0.1,0),labels=c(0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(v=sigGuild_gl,lty=2)
rect.dendrogram(dend, h=0.5, which=c(1:6),border="transparent",
                lty = 5, lwd = 0, col=rgb(0.2,0.2,0.2,0.15),horiz=T)

myClusts<-as.data.frame(cutree(overlap_clust_gl,h=sigGuild_gl))
myClusts$species=rownames(myClusts)
colnames(myClusts)<-c("myCluster","species")
myClusts<-myClusts%>%separate(species,into=c("species","sizecat"),sep="_")

GLclusts<-read.csv("GarrisonLink_clusters.csv")%>%
  mutate(guildCol=brewer.pal(6,"Set1")[guild],
         species=str_to_title(species))
myGLclusts<-left_join(myClusts,GLclusts)
  
colored_bars(colors = myGLclusts$guildCol, dend = dend, rowLabels = NA, horiz=T)

# Step 3: Run dev.off() to create the file!
dev.off()



#Unsure how to better control the look of the dendrogram, probably there's a package out there that looks nice





#################################################
#The latter ~half, since 1997 when Garrison and Link had data through
#Creating a proportion matrix for the prey items for my predators
propLong_new<-prey19%>%
  filter(year>1997)%>%
  mutate(species=str_to_title(pdcomnam))%>% #Maybe want common names later, but the myspecies is scientific so speed
  group_by(species,sizecat)%>%
  mutate(total_w=sum(pyamtw),
         total_v=sum(pyamtv))%>%
  group_by(species,sizecat,gl_prey)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()


#Following through to get the "real" overlap matrix
props1_new<-propLong_new[,-max(ncol(propLong_new))] #Cut off the volume prop for cleanliness
colnames(props1_new)<-c("species1","sizecat1","gl_prey","prop_w1")
props2_new<-propLong_new[,-max(ncol(propLong_new))]
colnames(props2_new)<-c("species2","sizecat2","gl_prey","prop_w2")

overlap_schoener_new<-full_join(props1_new,props2_new)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat_new<-pivot_wider(overlap_schoener_new,id_cols=c(species1,sizecat1),names_from=c(species2,sizecat2),values_from = s_do)
overlap_mat_new<-as.matrix(overlap_mat_new[,3:ncol(overlap_mat_new)])
#rownames(overlap_mat_new)<-gsub(" ","\n",colnames(overlap_mat_new))
#colnames(overlap_mat_new)<-gsub(" ","\n",colnames(overlap_mat_new))
#overlap_mat[overlap_mat==1]<-NA
for (i in 1:nrow(overlap_mat_new)) {
  for (j in 1:ncol(overlap_mat_new)) {
    overlap_mat_new[i,j]<-ifelse(j>i,NA,overlap_mat_new[i,j])
  }
}



#Bootstrapping
propMat_new<-propLong_new%>%
  pivot_wider(id_cols=c(species,sizecat),names_from = gl_prey,values_from = prop_w)%>% #Species are rows, prey are columns. Flip id and names if need opposite
  mutate(species=paste(species,sizecat,sep="_"))%>%
  select(-c("Empty","Unobserved","Miscellaneous","sizecat"))
propMat_new<-select(propMat_new,species,order(colnames(propMat)))
propMat_new[is.na(propMat_new)]<-0


reps=1
nreps=250
bootDiffs_new<-numeric()
while (reps<=nreps) {
  bootDiffs_new<-c(bootDiffs_new,resample(propMat_new))
  reps=reps+1
}
hist(bootDiffs_new)
sigGuild_new<-quantile(bootDiffs_new,probs=0.95)

#Visual of the actual matrix with significance indicated
#library(plot.matrix)
#par(mar=c(6,6,5,5.5))
#plot(overlap_mat,axis.row=list(side=2,las=1),col=viridis::viridis(n=100,option="B"),
#     polygon.key = list(border=NA), key=list(),xlab="",ylab="")

#Better visual (ggplot as always)
as.data.frame(overlap_mat_new)%>%
  mutate(species1=rownames(overlap_mat_new))%>%
  pivot_longer(cols=colnames(overlap_mat_new),names_to = "species2", values_to = "s_do")%>%
  mutate(species1=sort(species1,decreasing=T),species2=species2,
         sig=ifelse(s_do>sigGuild,"bold.italic","plain"))%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=s_do))+
  geom_text(aes(species1,species2,label=round(s_do,digits=2),fontface=sig,color=sig),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100],viridis::viridis(100,option="B")[sigGuild*100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5))

#Mean overlap, how similar are the diets in this time period
mean(overlap_mat_new[which(overlap_mat_new<1)],na.rm=T)
sd(overlap_mat_new[which(overlap_mat_new<1)],na.rm=T)
range(overlap_mat_new[which(overlap_mat_new<1)],na.rm=T)
hist(overlap_mat_new[which(overlap_mat_new<1)])

#Clustering these out
overlap_clust_new<-hclust(as.dist(1-overlap_mat_new),method="average")
#Trying dendextend to pretty up the dendrogram
dend_new<-as.dendrogram(overlap_clust_new)

#To save them out, there is some difficulty with the export losing the colorbar
png(file = "../Figures/2021 Garrison and Link Rep/fancy_2019dendrogram.png",   # The directory you want to save the file in
    width = 1250, # The width of the plot in inches
    height = 700) # The height of the plot in inches
par(mar=c(5,1,2,12))

dend_new %>% 
  set("branches_lwd", 4) %>%
  # Custom labels
  set("labels_cex", 1) %>%
  set("labels_col", value = viridis::viridis(12,end=0.8),h = sigGuild_new) %>%
  set("branches_k_color", value = viridis::viridis(12,end=0.8), h = sigGuild_new) %>%
  plot(horiz=TRUE,main="                               1998-2019 Dendrogram",axes=F,xlab="Similarity")
axis(side=1,at=c(0.6,0.5,0.4,0.3,0.2,0.1,0),labels=c(0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(v=sigGuild_gl,lty=2)
rect.dendrogram(dend_new, h=0.5, which=c(1:4),border="transparent",
                lty = 5, lwd = 0, col=rgb(0.2,0.2,0.2,0.15),horiz=T)

myClusts_new<-as.data.frame(cutree(overlap_clust_new,h=sigGuild_new))
myClusts_new$species=rownames(myClusts_new)
colnames(myClusts_new)<-c("myCluster","species")
myClusts_new<-myClusts_new%>%separate(species,into=c("species","sizecat"),sep="_")

GLclusts<-read.csv("GarrisonLink_clusters.csv")%>%
  mutate(guildCol=brewer.pal(6,"Set1")[guild],
         species=str_to_title(species))
myGLclusts_new<-left_join(myClusts_new,GLclusts)

colored_bars(colors = myGLclusts_new$guildCol, dend = dend_new, rowLabels = NA, horiz=T)

# Step 3: Run dev.off() to create the file!
dev.off()
cutree(overlap_clust_new,h=sigGuild_new)

#Unsure how to better control the look of the dendrogram, probably there's a package out there that looks nice


t.test(overlap_mat_gl[which(overlap_mat_gl<1)],overlap_mat_new[which(overlap_mat_new<1)])



# Picking out prey of interest (new to GoM, leaving GoM, endangered) -------------------------------

#Black sea bass, blue crab, lobster, northern shrimp, atlantic salmon, butterfish, others? tautog, longfin squid
speciesofinterest<-c("Centropristis striata","Callinectes","Pandalus borealis","Salmo salar","Peprilus triacanthus")
#Pulling out these species
preyInterest<-prey19%>%
  filter(grepl(paste(speciesofinterest,collapse="|"),ignore.case=T,pynam))

preyInterest%>%
  group_by(year,pynam)%>%
  summarise(mass=sum(pyamtw,na.rm=T))%>%
  ggplot(aes(year,mass,color=pynam))+
  geom_line(show.legend = F)+
  geom_point(show.legend = F)+
  facet_wrap(~pynam)




# Climate Vulnerability (Scores from Hale et al. 2016) --------------------

ggplot(vulScores,aes(Potential,Vulnerability))+
  geom_point(aes(size=Vulnerability_Certainty))

plot(vulScores[,c(2:4,6,8)])




# NMDS? -------------------------------------------------------------------

# Making a diet matrix ####
preyMat_NMDS<-filter(prey19,pynam!="EMPTY")%>%
  pivot_wider(id_cols=c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"),
                     names_from=collsci,values_from=pyamtw,values_fn=sum)
preyMat_NMDS[is.na(preyMat_NMDS)]<-0
#Number of prey categories
envCols<-length(c("svspp","pdcomnam","pdscinam","pdsex","pdid","pdlen","day","month","year","season","geoarea","pdgutw","declon","declat"))
ncol(preyMat_NMDS)-envCols

sampleMat<-preyMat_NMDS%>% #If you're subsetting down
  group_by(pdscinam)%>%
  sample_n(300)


envMat_NMDS<-sampleMat[,1:envCols]
preyMat_NMDS<-sampleMat[,(envCols+1):ncol(preyMat_NMDS)]

#Sparseness
(sum(preyMat_NMDS==0))/(nrow(preyMat_NMDS)*ncol(preyMat_NMDS))*100


#Sample a random selection of rows from the full matrix--it's too big for the NMDS
envMat_NMDS<-envMat_NMDS[rowSums(preyMat_NMDS)>0,]
preyMat_NMDS<-preyMat_NMDS[rowSums(preyMat_NMDS)>0,colSums(preyMat_NMDS)>0]
#Where were they lost?
table(as.character(envMat_NMDS$pdcomnam))
#For the all diets NMDS
original.dist<-vegdist(preyMat_NMDS)
stress_values<-numeric(6)
r2<-numeric(6)

for (n in 1:6) {
  nmds.resu <- metaMDS(preyMat_NMDS, k=n, distance = "bray", try=250, autotransform=F)
  stress_values[n]<-nmds.resu$stress
  nmds.scores<-vegan::scores(nmds.resu)
  nmds.dist<-dist(nmds.scores)
  r2[n]<-summary(lm(original.dist~nmds.dist))[[8]]
}
plot(stress_values, xlab="Number of axes", ylab="Stress",type="b")
abline(h=0.2,col="red")

View(stress_values) 

#Go back and create the output for the 2 dimensions NMDS
preyNMDS<-metaMDS(preyMat_NMDS, distance = "bray", k = 2, try=250, autotransform=F)
r2<-summary(lm(original.dist~dist(vegan::scores(preyNMDS))))[[8]]
actualStress<-preyNMDS$stress
stressplot(preyNMDS) #This is the visual of stress, the divergence of observed and ordinated distance. It's random, that's good

#Print the species scores and sample scores
NMDS_species<-as.data.frame(preyNMDS$species)
NMDS_scores<-as.data.frame(preyNMDS$points)

#PERMANOVA for the interaction of Year and fortnight 
set.seed(42)
adonis(original.dist~as.character(pdcomnam),data=envMat_NMDS,permutations=10000,method="bray")

#significant pairwise test of the species
pairwise.perm.manova(original.dist,envMat_NMDS$pdcomnam,nperm=1000)

fdisp<-betadisper(original.dist,envMat_NMDS$pdcomnam)
fdisp
permutest(fdisp)


dist_multi_centroids(original.dist,envMat_NMDS$pdcomnam)





#ISA to see what diet items might be associated with the different years that make them different from each other
spISA = multipatt(as.data.frame(preyMat_NMDS), as.character(envMat_NMDS$pdcomnam),
                          func = "IndVal.g", duleg=TRUE, control = how(nperm=4999))

#What species?
summary(spISA) #
spISA$str

#Extract them, with the year they're significant for
spIndSp<-dplyr::select(subset(spISA$sign, p.value<0.05),index,stat,p.value)
spIndSp$Species<-rownames(spIndSp)
spIndSp$fish<-colnames(spISA$B)[spIndSp$index]


#Need new axes R2 scores (from PC-ORD)
ggplot(data=NMDS_scores,aes(MDS1,MDS2))+
  geom_polygon(data=NMDS_scores%>%group_by(species=as.character(envMat_NMDS$pdcomnam))%>%slice(chull(MDS1,MDS2)),
               aes(x=MDS1,y=MDS2,fill=species,color=species),alpha=0.1,lwd=1.5,show.legend = F)+
  geom_point(size=8,aes(fill=envMat_NMDS$pdcomnam),shape=21)+
  xlab("Axis 1 (%)")+ylab("Axis 2 (%)")+
  #stat_ellipse(data=allNMDS_DF,aes(color=paste(allEnv.DF$fortnight,allEnv.DF$Year)),
  #            level=0.95,lwd=1.1,show.legend = F)+
  geom_segment(data=NMDS_species, #all species
               aes(x=0,y=0,xend=MDS1,yend=MDS2),color="grey50",lwd=1,alpha=0.5,show.legend = F)+ 
  #geom_segment(data=filter(NMDS_species,rownames(NMDS_species)%in%yIndSp$Species), #Just the indicator species
  #             aes(x=0,y=0,xend=MDS1,yend=MDS2,color=yIndSp$time),lwd=1.1,show.legend = F)+ #Color coded to the group they're indicating
  #geom_label(data=filter(NMDS_species,rownames(NMDS_species)%in%yIndSp$Species), #Just the indicator species
  #           aes(x=MDS1*1.1,y=MDS2*1.1,label=yIndSp$Species,color=yIndSp$time),size=8,show.legend = F)+ #Coded to the group they indicate
  scale_color_manual(values=colorspace::rainbow_hcl(9,c=100))+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100))+
  geom_text(aes(x=12,y=5,label=paste0("Stress = ",round(actualStress*100,digits=2)),hjust=1),size=10)+
  #geom_text(aes(x=1.5,y=1.4,label=paste0("MRPP Year p = ",round(countYear_MRPP$Pvalue,digits=4)),hjust=1),size=9)+
  theme(legend.position=c(0.125,0.75),legend.background=element_rect(color="black"),legend.margin=margin(4,8,5,8),
        legend.text=element_text(size=20),legend.title=element_text(size=22))+
  guides(fill=guide_legend(title.hjust=0,label.vjust=0,override.aes=list(shape=21)))
  coord_cartesian(xlim=c(-1.5,2),ylim=c(-1,1))


# Prey Lengths ------------------------------------------------------------

pylen19<-pylen19%>%
  left_join(uniquePrey19)%>%
  mutate(declon=-declon)%>%
  mutate(Species=str_to_title(gsub(" ","\n",pdcomnam)),
         species=str_to_title(pdcomnam))
#What amount don't have full predator info?
nrow(pylen19[is.na(pylen19$pdcomnam),])/nrow(pylen19)*100

#How many prey got measured from each species
table(pylen19$Species)

#How many prey were measured?
n_distinct(pylen19$pynam)

#The pylen are in mm
range(pylen19$pylen,na.rm=T)/10
#The pdlen are in cm
range(pylen19$pdlen,na.rm=T)

#How do the predators compare to their prey
pylen19%>%
  filter(!is.na(pdcomnam))%>%
  ggplot()+
  geom_point(aes(pdlen,pylen/10,color=Species),shape=21,size=2,fill="transparent")+
  geom_abline(aes(slope=1,intercept=0))+
  geom_smooth(aes(pdlen,pylen/10),method="lm",lty=1,se=F,color="black")+
  geom_smooth(data=pylen19%>%filter(!is.na(pdcomnam))%>%mutate(lengths=as.character(pdlen))%>%
                    group_by(species,lengths)%>%mutate(N=n())%>%filter(pylen==max(pylen)),
             aes(pdlen,pylen/10),method="lm",lty=3,se=F,color="black")+
  geom_smooth(data=pylen19%>%filter(!is.na(pdcomnam))%>%mutate(lengths=as.character(pdlen))%>%
                group_by(species,lengths)%>%mutate(N=n())%>%filter(pylen==min(pylen)),
              aes(pdlen,pylen/10),method="lm",lty=3,se=F,color="black")+
  scale_color_manual(values=colorspace::rainbow_hcl(9,c=100,l=60),guide="none")+
  scale_x_continuous(limits=c(0,130),expand=expansion(mult=0.01),name="Predator length (cm)")+
  scale_y_continuous(limits=c(0,130),expand=expansion(mult=0.01),name="Prey length (cm)")+
  theme(legend.position=c(0.33,0.75),legend.background=element_rect(color="black",fill="white"))+
  facet_wrap(~species)

#Maximum sizes for each predator size to show that max increases but min doesn't really
pylen19%>%filter(!is.na(pdcomnam))%>%mutate(lengths=as.character(pdlen))%>%
         group_by(lengths)%>%filter(pylen==max(pylen)|pylen==min(pylen))%>%
  mutate(size=ifelse(pylen==max(pylen),"Max","Min"))%>%
  ggplot()+
  geom_point(aes(pdlen,pylen/10,shape=size,fill=size),size=3,show.legend = F)+
  scale_fill_manual(values=colorspace::rainbow_hcl(2,c=100))+
  scale_shape_manual(values=c(21,22))+
  scale_x_continuous(limits=c(0,130),expand=expansion(mult=0.01),name="Predator length (cm)")+
  scale_y_continuous(limits=c(0,130),expand=expansion(mult=0.01),name="Prey length (cm)")
  
#What size distribution are the different species consuming?
pylen19%>%
  filter(!is.na(pdcomnam))%>%
  ggplot()+
  geom_point(aes(Species,pylen/10,fill=Species),shape=21,alpha=0.7,show.legend = F)+
  geom_violin(aes(Species,pylen/10,fill=Species),width=1.5,show.legend = F)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,end=290))+
  scale_y_continuous(name="Prey Length (cm)")
#And range
pylen19%>%
  filter(!is.na(pdcomnam))%>%
  group_by(Species)%>%
  summarise(minL=min(pylen,na.rm=T),maxL=max(pylen,na.rm=T))

#What is the ratio of prey to predator lengths for the different species
pylen19<-mutate(pylen19,pypdRatio=(pylen/10)/pdlen)
#How does the skewness vary by the different species
skews<-pylen19%>%
  filter(!is.na(pypdRatio))%>%
  group_by(Species)%>%
  summarise(skew=round(skewness(pypdRatio),digits=4),
            m=max(pypdRatio))
pylen19%>%
  filter(!is.na(pdcomnam))%>%
  left_join(skews)%>%
  ggplot()+
  geom_point(aes(fct_reorder(Species,skew),pypdRatio,fill=Species),shape=21,alpha=0.8,size=2,show.legend=F)+
  geom_violin(aes(fct_reorder(Species,skew),pypdRatio,fill=Species),show.legend=F)+
  #geom_text(data=skews,aes(Species,m+0.075,label=paste0("Skew=",skew)),size=6)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,l=60))+
  scale_y_continuous(name="Ratio of Prey Length:Predator Length",expand=expansion(add=c(0.035,0.12)))+
  scale_x_discrete(name="Predator Species")
#What percent of prey are below 50% of the length of their predator
nrow(filter(pylen19,pypdRatio<=0.5))/nrow(filter(pylen19,!is.na(pypdRatio)))*100
#For each of the species...
pylen19%>%
  filter(!is.na(pypdRatio))%>%
  group_by(Species)%>%
  summarise(pSmall=sum(pypdRatio<=0.5)/n())

#What prey are being eaten that are longer than the predator?
table(filter(pylen19,pypdRatio>=1)$pynam)


#How is the size of prey changing over time?
pylen19%>%
  filter(!is.na(pdcomnam))%>%
  ggplot(aes(year,pylen))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Species)
summary(preyYear<-lm(pylen~year*pdcomnam,data=filter(pylen19,!is.na(pdcomnam))))
pylen19%>%
  filter(!is.na(pdcomnam))%>%
  ggplot(aes(year,pdlen))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Species)
summary(preyYear<-lm(pdlen~year*pdcomnam,data=filter(pylen19,!is.na(pdcomnam))))
pylen19%>%
  filter(!is.na(pdcomnam))%>%
  ggplot(aes(year,pypdRatio))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Species)
summary(preyYear<-lm(pypdRatio~year*pdcomnam,data=filter(pylen19,!is.na(pdcomnam))))



#How does location influence prey size?
load("googleMap.zoom6.eastCoast.R")
ggmap(zoom6)+
  geom_point(data=filter(pylen19,!is.na(pdcomnam))%>%mutate(Species=str_to_title(pdcomnam)),
             aes(declon,declat,fill=pdcomnam),shape=21,alpha=0.8,size=3,show.legend=F)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,end=290))+
  ylab("Latitude")+xlab("Longitude")+
  facet_wrap(~Species)

#Can't really see anything, but are the size of prey influenced by space?
ggmap(zoom6)+
  geom_point(data=filter(pylen19,!is.na(pdcomnam))%>%mutate(Species=str_to_title(pdcomnam)),
             aes(declon,declat,fill=pdcomnam,size=pylen),shape=21,alpha=0.8,show.legend=F)+
  scale_fill_manual(values=colorspace::rainbow_hcl(9,c=100,end=290))+
  ylab("Latitude")+xlab("Longitude")+
  facet_wrap(~Species)


#Does prey size change with lat-long?
par(mfrow=c(1,2))
plot(filter(pylen19,!is.na(pdcomnam)&declon>-77)$declon,filter(pylen19,!is.na(pdcomnam)&declon>-77)$pylen)
abline(lm(filter(pylen19,!is.na(pdcomnam)&declon>-77)$pylen~filter(pylen19,!is.na(pdcomnam)&declon>-77)$declon),col="red",lwd=2)
summary(sizeLONG<-lm(filter(pylen19,!is.na(pdcomnam)&declon>-77)$pylen~filter(pylen19,!is.na(pdcomnam)&declon>-77)$declon))
plot(filter(pylen19,!is.na(pdcomnam)&declon>-77)$declat,filter(pylen19,!is.na(pdcomnam)&declon>-77)$pylen)
abline(lm(filter(pylen19,!is.na(pdcomnam)&declon>-77)$pylen~filter(pylen19,!is.na(pdcomnam)&declon>-77)$declat),col="red",lwd=2)
summary(sizeLAT<-lm(filter(pylen19,!is.na(pdcomnam)&declon>-77)$pylen~filter(pylen19,!is.na(pdcomnam)&declon>-77)$declat))





# Centers of Consumption --------------------------------------------------
#Loligo squid, which really means Longfin Squid
lSquid_c<-filter(prey19,grepl("loligo",pynam,ignore.case=T))%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=weighted.mean(declat,pyamtw))
nrow(lSquid_c)
ggmap(zoom6)+
  geom_point(data=lSquid_c,aes(-declon,declat,color=pyamtw),show.legend = F)+
  geom_segment(data=lSquid_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)

#American butterfish
butter_c<-filter(prey19,pynam=="PEPRILUS TRIACANTHUS")%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=weighted.mean(declat,pyamtw))
nrow(butter_c)
ggmap(zoom6)+
  geom_point(data=butter_c,aes(-declon,declat,color=pyamtw),show.legend = F)+
  geom_segment(data=butter_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)

#Black sea bass
bsb_c<-filter(prey19,pynam=="CENTROPRISTIS STRIATA")%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=weighted.mean(declat,pyamtw))
nrow(bsb_c)
ggmap(zoom6)+
  geom_point(data=bsb_c,aes(-declon,declat,color=pyamtw),show.legend = F)+
  geom_segment(data=bsb_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)

#Atlantic mackerel
mack_c<-filter(prey19,pynam=="SCOMBER SCOMBRUS")%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=weighted.mean(declat,pyamtw))
nrow(mack_c)
ggmap(zoom6)+
  geom_point(data=mack_c,aes(-declon,declat,color=pyamtw),show.legend = F)+
  geom_segment(data=mack_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)

#Silver Hake
shake_c<-filter(prey19,pynam=="MERLUCCIUS BILINEARIS")%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=weighted.mean(declat,pyamtw))
nrow(shake_c)
ggmap(zoom6)+
  geom_point(data=shake_c,aes(-declon,declat,color=pyamtw),show.legend = F)+
  geom_segment(data=shake_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)


#Together
invasives_c<-bind_rows(shake_c,mack_c,bsb_c,butter_c,lSquid_c)
ggmap(zoom6)+
  geom_point(data=invasives_c,aes(-declon,declat,color=collsci,size=pyamtw))+
  geom_segment(data=invasives_c,aes(x=-64,xend=-66,y=latMean,yend=latMean,color=collsci),lwd=2)+
  scale_color_brewer(palette="Set1",name="Species")+
  scale_size_continuous(name="Prey Mass (g)")+
  theme(legend.position="top")+
  guides(color=guide_legend(nrow=3,title.position="top"),size=guide_legend(title.position = "top"))+
  facet_wrap(~decade,nrow=1)

###
#Some that might potentially be shifting out?
###

#Lobster
lobster_c<-filter(prey19,grepl("homarus americanus",pynam,ignore.case = T))%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=mean(declat))
nrow(lobster_c)
ggmap(zoom6)+
  geom_point(data=lobster_c,aes(-declon,declat,color=year),show.legend = F)+
  geom_segment(data=lobster_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2,lty=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)

#Northern Shrimp
nshrimp_c<-filter(prey19,pynam=="PANDALUS BOREALIS")%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=mean(declat))
nrow(nshrimp_c)
ggmap(zoom6)+
  geom_point(data=nshrimp_c,aes(-declon,declat,color=year),show.legend = F)+
  geom_segment(data=nshrimp_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2,lty=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)

#Calanus
zoo_c<-filter(prey19,pynam=="CALANUS FINMARCHICUS")%>%
  mutate(decade=paste0(substr(year,1,3),"0s"))%>%
  group_by(decade)%>%
  mutate(latMean=mean(declat))
nrow(zoo_c)
ggmap(zoom6)+
  geom_point(data=zoo_c,aes(-declon,declat,color=year),show.legend = F)+
  geom_segment(data=zoo_c,aes(x=-64,xend=-66,y=latMean,yend=latMean),color="red",lwd=2,lty=2)+
  scale_color_distiller(palette="Reds")+
  facet_wrap(~decade,nrow=1)




# Centers of Distribution -------------------------------------------------

#This is all copied from the exploration script, then adjusted to keep species too
springT<-read_csv("../Trawl Data/NMFS Trawls/Spring/22561_UNION_FSCS_SVCAT.csv", 
                  col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                   STATION = col_double(), STRATUM = col_double(), 
                                   TOW = col_double()))

springT$Year<-as.numeric(substr(springT$CRUISE6,1,4))

fallT<-read_csv("../Trawl Data/NMFS Trawls/Fall/22560_UNION_FSCS_SVCAT.csv", 
                col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                 STATION = col_double(), STRATUM = col_double(), 
                                 TOW = col_double()))

fallT$Year<-as.numeric(substr(fallT$CRUISE6,1,4))

#Invasive species svs codes
  #503--Longfin Squid
  #72--Silver Hake
  #121--Atlantic Mackerel
  #131--Butterfish
  #141--Black Sea Bass
invasives<-c("503","072","121","131","141")
invasiveTrawls<-filter(bind_rows(springT,fallT),SVSPP%in%invasives)

trawlLocations_Fall<-read_csv("../Trawl Data/NMFS Trawls/Fall/22560_UNION_FSCS_SVSTA.csv", 
                              col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                               STATION = col_double(), STRATUM = col_double(), 
                                               TOW = col_double(), AREA = col_character(),
                                               EST_MONTH=col_double(),EST_DAY=col_double()))%>%
  dplyr::select(CRUISE6,ID,STATION,STRATUM,TOW,LONG=DECDEG_BEGLON,LAT=DECDEG_BEGLAT,
                Year=EST_YEAR,Month=EST_MONTH,Day=EST_DAY,Time=EST_TIME,AREA,DESSPEED,TOWDUR,SURFTEMP,BOTTEMP)

trawlLocations_Spring<-read_csv("../Trawl Data/NMFS Trawls/Spring/22561_UNION_FSCS_SVSTA.csv", 
                                col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                                 STATION = col_double(), STRATUM = col_double(), 
                                                 TOW = col_double(), AREA=col_character(),
                                                 EST_MONTH=col_double(),EST_DAY=col_double()))%>%
  dplyr::select(CRUISE6,ID,STATION,STRATUM,TOW,LONG=DECDEG_BEGLON,LAT=DECDEG_BEGLAT,
                Year=EST_YEAR,Month=EST_MONTH,Day=EST_DAY,Time=EST_TIME,AREA,DESSPEED,TOWDUR,SURFTEMP,BOTTEMP)

invasives_catchLocations<-left_join(invasiveTrawls,bind_rows(trawlLocations_Fall,trawlLocations_Spring))%>%
  mutate(Date=ymd(paste(Year,Month,Day,sep="-")),
         doy=yday(Date))


invasives_d<-filter(invasives_catchLocations,!is.na(EXPCATCHWT))%>%
  mutate(decade=paste0(substr(Year,1,3),"0s"),
                    catch=as.numeric(EXPCATCHWT))%>%
  group_by(decade,SVSPP)%>%
  mutate(latMean=weighted.mean(LAT,catch,na.rm=T))
nrow(invasives_d)
ggmap(zoom6)+
  geom_point(data=invasives_d,aes(LONG,LAT,color=SVSPP),show.legend = F)+
  geom_segment(data=invasives_d,aes(x=-64,xend=-66,y=latMean,yend=latMean,color=SVSPP),lwd=2)+
  scale_color_brewer(palette="Set1")+
  facet_wrap(~decade,nrow=1)

#Comparing with diet changes
invasives_diffD_wide<-pivot_wider(invasives_d,id_cols=SVSPP,names_from=decade,values_from=latMean,values_fn=mean)

invasives_diffD<-invasives_d%>%group_by(SVSPP,decade)%>%summarise(latMean=mean(latMean))%>%
  group_by(SVSPP)%>%mutate(diffMean=if_else(is.na(lead(latMean)),latMean-min(latMean),lead(latMean)-latMean))

ggplot(invasives_diffD,aes(decade,latMean,color=SVSPP))+
  geom_line(aes(group=SVSPP))+
  geom_point()



USCoast<-rgdal::readOGR("../Mapping Data/tl_2019_us_coastline.shp",)
eastCoast<-USCoast[USCoast@data$NAME=="Atlantic",]
eastCoast<-spTransform(eastCoast,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
eastCoast<-as(eastCoast,"SpatialPoints")

points<-SpatialPoints(invasives_catchLocations[,c("LONG","LAT")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

eastCoast.sf<-st_as_sf(eastCoast)
points.sf<-st_as_sf(points)

eastCoast.sf$id = 1:nrow(eastCoast.sf) # make sure to have unique id to trace selected features later

trawl_w_nearest_coast = st_join(points.sf, eastCoast.sf, join = st_nearest_feature)
trawl_w_nearest_coast<-trawl_w_nearest_coast%>%st_drop_geometry()%>%left_join(eastCoast.sf)%>%bind_cols(invasives_catchLocations)
trawl_w_nearest_coast<-st_as_sf(trawl_w_nearest_coast)

trawl_w_nearest_coast<-data.frame(as(trawl_w_nearest_coast,"Spatial"))%>%
  rename(trawl_long=LONG,trawl_lat=LAT,coast_long=coords.x1,coast_lat=coords.x2)%>%dplyr::select(-optional)

#Checking visually
ggmap(zoom6)+
  geom_point(data=as.data.frame(eastCoast@coords),aes(x=coords.x1,y=coords.x2),size=0.5)+
  geom_point(data=trawl_w_nearest_coast[50000:50111,],aes(coast_long,coast_lat))+
  geom_point(data=trawl_w_nearest_coast[50000:50111,],aes(x=trawl_long,y=trawl_lat),color="red")+
  geom_segment(data=trawl_w_nearest_coast[50000:50111,],aes(x=trawl_long,xend=coast_long,y=trawl_lat,yend=coast_lat))

library(smoothr)

#with gproject
USCoast<-rgdal::readOGR("../Mapping Data/tl_2019_us_coastline.shp",)
eastCoast<-USCoast[USCoast@data$NAME=="Atlantic",]
eastCoast<-spTransform(eastCoast,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
eastCoast<-raster::crop(eastCoast,raster::extent(-180,180,30.5,90)) #CROPPING TO JUST BELOW CAPE HATTERAS, WHERE TRAWLS SEEM TO STOP
eastCoast.union<-st_union(st_as_sf(eastCoast))
eastCoast.line<-st_line_merge(eastCoast.union)
eastCoast.line<-st_cast(eastCoast.line,to="LINESTRING")
trueCoast.line<-as.matrix(eastCoast.line[[1]])
plot(trueCoast.line)


#Trying unsmoothed
trueCoast.sp<-as.data.frame(trueCoast.line)%>%
  rename(long=V1,lat=V2)
trueCoast.sp<-SpatialPoints(trueCoast.sp,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
trueCoast.sp<-spTransform(trueCoast.sp,CRSobj=CRS("+proj=utm +zone=19"))
trueCoast.sp<-Line(trueCoast.sp)
trueCoast.sp<-Lines(list(trueCoast.sp),ID="coast")
trueCoast.sp<-SpatialLines(list(trueCoast.sp))

plot(trueCoast.sp)
gLength(trueCoast.sp)

points<-SpatialPoints(invasives_catchLocations[,c("LONG","LAT")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
points.sf<-spTransform(points,CRSobj = CRS("+proj=utm +zone=19"))

distAlongCoast<-gProject(trueCoast.sp,points.sf,normalize=T)
pointAlongCoast<-gInterpolate(trueCoast.sp,distAlongCoast,normalized=T)
proj4string(pointAlongCoast)<-CRS("+proj=utm +zone=19")
pointAlongCoast<-spTransform(pointAlongCoast,CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

invasives_catchLocations[,"distAlongCoast"]<-distAlongCoast
invasives_catchLocations[,c("coast_long","coast_lat")]<-pointAlongCoast@coords
invasives_catchLocations<-rename(invasives_catchLocations,trawl_long=LONG,trawl_lat=LAT)
write.csv(invasives_catchLocations,"invasiveRangeShifts_trawlCatches+truecoastDistance.csv",row.names=F)

ggmap(zoom6)+
  geom_point(data=invasives_catchLocations,aes(coast_long,coast_lat,color=distAlongCoast))+
  geom_point(data=invasives_catchLocations,aes(trawl_long,trawl_lat,fill=distAlongCoast),shape=21)+
  scale_color_distiller(palette="Reds")+
  scale_fill_distiller(palette="Blues")




#Smoothing the coast

points<-SpatialPoints(invasives_catchLocations[,c("LONG","LAT")],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
points.sf<-spTransform(points,CRSobj = CRS("+proj=utm +zone=19"))

smoothCoast<-smooth_ksmooth(trueCoast.line,smoothness=10000)
smoothCoast<-round(smoothCoast,digits=4)
smoothCoast<-smoothCoast[!duplicated(smoothCoast),]
plot(smoothCoast)

smoothCoast.sf<-as.data.frame(smoothCoast)%>%
  rename(long=V1,lat=V2)
smoothCoast.sf<-SpatialPoints(smoothCoast.sf,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
smoothCoast.sf<-spTransform(smoothCoast.sf,CRSobj=CRS("+proj=utm +zone=19"))
smoothCoast.sf<-Line(smoothCoast.sf)
smoothCoast.sf<-Lines(list(smoothCoast.sf),ID="coast")
smoothCoast.sf<-SpatialLines(list(smoothCoast.sf))

gLength(smoothCoast.sf)

distAlongCoast.s<-gProject(smoothCoast.sf,points.sf,normalize=T)
pointAlongCoast.s<-gInterpolate(smoothCoast.sf,distAlongCoast.s,normalized=T)
proj4string(pointAlongCoast.s)<-CRS("+proj=utm +zone=19")
pointAlongCoast.s<-spTransform(pointAlongCoast.s,CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

invasives_catchLocations[,"distAlongCoast.s"]<-distAlongCoast.s
invasives_catchLocations[,c("smoothCoast_long","smoothCoast_lat")]<-pointAlongCoast.s@coords
write.csv(invasives_catchLocations,"invasiveRangeShifts_trawlCatches+smoothcoastDistance.csv",row.names=F)

ggplot(invasives_catchLocations)+
  geom_point(aes(coast_long,coast_lat,color=distAlongCoast))+
  geom_point(aes(trawl_long,trawl_lat,fill=distAlongCoast),shape=21)+
  scale_color_distiller(palette="Reds")+
  scale_fill_distiller(palette="Blues")



new<-read.csv("invasiveRangeShifts_trawlCatches+truecoastDistance.csv")


ggplot(new)+
  geom_point(aes(coast_long,coast_lat,color=distAlongCoast))+
  geom_point(aes(trawl_long,trawl_lat,fill=distAlongCoast),shape=21)+
  scale_color_distiller(palette="Reds")+
  scale_fill_distiller(palette="Blues")





# Fish Size Distributions -------------------------------------------------

sizeClasses<-read.csv("GarrisonLink_predSizeClasses.csv")%>%
  left_join(spec,by=c("species_comnam"="LOGGED_SPECIES_NAME"))


fallSizes<-read_csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVLEN.csv",
                    col_types = "cccccccnnn")
springSizes<-read_csv("../Trawl Data/NMFS Trawls/2022 Redownload/Spring/22561_UNION_FSCS_SVLEN.csv",
                      col_types = "cccccccnnn")


allSizes<-bind_rows(fallSizes,springSizes)%>%
  left_join(sizeClasses)%>%
  filter(SVSPP%in%myspecies_svspp)%>%
  group_by(SVSPP)%>%
  mutate(sizeClass=case_when(between(LENGTH,min(small_min),min(small_max))~"S",
                             between(LENGTH,min(medium_min),min(medium_max))~"M",
                             between(LENGTH,min(large_min),min(large_max))~"L",
                             between(LENGTH,min(xlarge_min),min(xlarge_max))~"XL",
                             TRUE ~ "S")) #Because the only ones left are less than the small measure, but they're actually labeled small

sizeClassProps<-allSizes%>%
  group_by(ID,species_scinam,species_comnam,SVSPP)%>%
  mutate(N=sum(EXPNUMLEN))%>%
  group_by(ID,pdscinam=species_scinam,species_comnam,svspp=SVSPP,sizecat=sizeClass)%>%
  summarise(p=sum(EXPNUMLEN)/N)%>%distinct()

allGMRItrawls_clean<-read_csv("../Trawl Data/NMFS Trawls/Complete/NMFS_survdat_gmri_tidy.csv",
                              col_types="cccccccnnncnnTnnnnnnnnnnccncn")%>%
  mutate(stratum_full=str_pad(stratum, width = 5, pad = "0", side = "left"),
         tow_full=str_pad(tow, width = 3, pad = "0", side = "left"),
         station_full=str_pad(station, width = 4, pad = "0", side = "left"),
         ID=paste0(cruise6,stratum_full,tow_full,station_full))%>%
  left_join(sizeClasses,by=c("svspp"="SVSPP"))%>%
  filter(svspp%in%myspecies_svspp)%>%
  group_by(svspp)%>%
  mutate(sizecat=case_when(between(length_cm,min(small_min),min(small_max))~"S",
                           between(length_cm,min(medium_min),min(medium_max))~"M",
                           between(length_cm,min(large_min),min(large_max))~"L",
                           between(length_cm,min(xlarge_min),min(xlarge_max))~"XL",
                           TRUE ~ "S"),
         merge="yes") #Because the only ones left are less than the small measure, but they're actually labeled small
GMRIabundances<-allGMRItrawls_clean%>%
  dplyr::select(id,svspp,catchsex,abundance)%>%distinct()%>%
  group_by(id,svspp)%>%
  summarise(trawlabundance=sum(abundance))

sumGMRItrawls<-allGMRItrawls_clean%>%
  left_join(GMRIabundances)%>%
  dplyr::group_by(id,svspp,comname,sizecat,trawlabundance,merge)%>%
  summarise(N_adj=sum(numlen_adj),p=N_adj/trawlabundance,
            sizeabundance=p*trawlabundance)%>%distinct()%>%ungroup()

prey19<-prey19%>%
  mutate(svspp=str_pad(svspp,width=3,side="left",pad="0"),
         id=paste0(cruise6,
                   str_pad(station,width=3,side="left",pad="0"),
                   str_pad(stratum,width=3,side="left",pad="0")),
         dietID=paste0(svspp,
                       pdsex,
                       str_pad(pdid,width=6,side="left",pad="0"),
                       str_pad(pdlen,width=3,pad="0",side="left"),
                       id))%>%
  group_by(id,svspp)%>%
  mutate(nDiets=n_distinct(dietID))%>%
  left_join(sizeClasses,by=c("svspp"="SVSPP"))%>%
  group_by(svspp)%>%
  mutate(sizecat2=case_when(between(pdlen,min(small_min),min(small_max))~"S",
                            between(pdlen,min(medium_min),min(medium_max))~"M",
                            between(pdlen,min(large_min),min(large_max))~"L",
                            between(pdlen,min(xlarge_min),min(xlarge_max))~"XL",
                            TRUE ~ "S"))%>%
  select(-c(species_scinam:xlarge_max))
  


strata_key <- list(
  "Georges Bank"          = as.character(13:23),
  "Gulf of Maine"         = as.character(24:40),
  "Southern New England"  = stringr::str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
  "Mid-Atlantic Bight"    = as.character(61:76))
preyGMRI<-left_join(prey19,sumGMRItrawls,by=c("svspp","sizecat2"="sizecat","id"))%>% #merge by sizecat2 because this is defined from Garrison and Link in the same way as sizecat is defined in the GMRI, clearly the diet df is defining sizecat some other way (especially if not entirely for spiny dogfish, sexual dimorphism maybe??)
  mutate(strata = str_pad(stratum, width = 5, pad = "0", side = "left"),
         strat_num = str_sub(strata, 3,4),
         survey_area =  dplyr::case_when(
           strat_num %in% strata_key$`Georges Bank`         ~ "GB",
           strat_num %in% strata_key$`Gulf of Maine`        ~ "GoM",
           strat_num %in% strata_key$`Southern New England` ~ "SNE",
           strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
           TRUE                                             ~ "stratum not in key"))

#Just cut out the diets that don't have trawl data (whether because GMRI didn't provide it or it's missing entirely)
preyGMRI_filter<-filter(preyGMRI,id%in%allGMRItrawls_clean$id)%>% #We know that GMRI didn't have the full trawls, but we only want to use the diets from their trawl list (keep things even) most of what's missing are the really southern and ScS diets
  filter(!is.na(p))%>%
  mutate(comname=str_to_sentence(comname),
         season=factor(str_to_sentence(season),levels=c("Spring","Fall")),
         geoarea=droplevels(geoarea))%>%
  ungroup()

#Each individual gets a single value, instead of by prey item
indGMRI_filter<-preyGMRI_filter%>%
  dplyr::select(-c(pynam:gl_prey))%>%
  distinct()%>%
  mutate(cruise6=as.character(cruise6),
         station=str_pad(station,width=4,pad="0",side="left"),
         pdsex=as.character(pdsex),
         pdid=str_pad(pdid,width=6,pad="0",side="left"))

allGeoCom<-dplyr::select(preyGMRI_filter,geoarea,comname)%>%
  distinct()%>%expand(geoarea,comname)


# Weighting Diets by Trawl Capture ----------------------------------------

clusterSD <- function(n,tM,M,pq,FW) { #Function to calculate the SD for cluster means, from Buckel et al. 1999
  ssq<- sum(M^2*(pq-FW)^2,na.rm=T)
  var<- (1/(n*(tM/n)^2))*ssq/(n-1)
  sd<-sqrt(var)
  return(sd)
}
#Where:
#n is the number of samples (trawls) in the GROUP
#tM is the total abundance of the fish in the GROUP, not equal to the sum of M because need to include those caught in trawls where a diet item may not have shown up
#M is the abundance of fish in a sample (trawl)
#pq is the p frequency or q relative abundance of prey item in that sample (trawl)
#FW is the mean F frequency or W relative abundance of the prey item in the GROUP
#Note that tM/n is equal to the mean number of fish caught in each sample (trawl) for that GROUP


#Need to have whole diet totals, and whole trawl measures
trawlDiets<-preyGMRI_filter%>%
  group_by(id,comname)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv))%>%
  group_by(id,pdscinam,gl_prey)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets))%>%
  dplyr::select(geoarea,year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,trawlabundance,pdscinam,comname,gl_prey,totalwt,totalv,nDiets,qikw,qikv,pik)%>% 
  distinct()%>% #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.
  ungroup()

nDiets<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(geoarea,year,season,comname,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(geoarea,year,season,comname)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(geoarea,year,season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(geoarea,year,season,comname)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(geoarea,year,season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(geoarea,year,season,comname)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietSum<-trawlDiets%>%
  group_by(geoarea,year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(geoarea,year,season,pdscinam,comname,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       trawlabundance,pik, Fk))%>%distinct()
trawlDietSum[is.na(trawlDietSum)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum%>%group_by(geoarea,year,season,pdscinam)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0


nDiets_ny<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(geoarea,season,comname,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(geoarea,season,comname)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun_ny<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(geoarea,season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(geoarea,season,comname)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_ny<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(geoarea,season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(geoarea,season,comname)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietSum_ny<-trawlDiets%>%
  group_by(geoarea,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_ny)%>%
  left_join(sumAbun_ny)%>%
  left_join(sumAbun_nEmpty_ny)%>%
  group_by(geoarea,season,pdscinam,comname,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       trawlabundance,pik, Fk))%>%distinct()
trawlDietSum_ny[is.na(trawlDietSum_ny)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum_ny%>%group_by(geoarea,season,pdscinam)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0



#Means over time just for the species
nDiets_spy<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(year,season,comname,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(year,season,comname)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun_spy<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(year,season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(year,season,comname)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_spy<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(year,season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(year,season,comname)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietSum_spy<-trawlDiets%>%
  group_by(year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_spy)%>%
  left_join(sumAbun_spy)%>%
  left_join(sumAbun_nEmpty_spy)%>%
  group_by(year,season,comname,pdscinam,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       trawlabundance,pik, Fk))%>%distinct()
trawlDietSum_spy[is.na(trawlDietSum_spy)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum_spy%>%group_by(year,season,comname)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 (or they're all empties)


#Means for a whole species, regardless of the geoarea
nDiets_sp<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(season,comname,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(season,comname)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun_sp<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(season,comname)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_sp<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(season,comname,id,trawlabundance)%>%distinct()%>%
  group_by(season,comname)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietSum_sp<-trawlDiets%>%
  group_by(season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_sp)%>%
  left_join(sumAbun_sp)%>%
  left_join(sumAbun_nEmpty_sp)%>%
  group_by(season,comname,pdscinam,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       trawlabundance,pik, Fk))%>%distinct()
trawlDietSum_sp[is.na(trawlDietSum_sp)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum_sp%>%group_by(season,comname)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 (or they're all empties)





#Need to have whole diet totals, and whole trawl measures
#need to adjust the trawl numbers first
trawlN_geo<-preyGMRI_filter%>%
  dplyr::select(id,svspp,trawlabundance)%>%distinct()%>%
  group_by(id)%>%
  summarise(wholetrawlabundance=sum(trawlabundance))
trawlDiets_geo<-preyGMRI_filter%>%
  left_join(trawlN_geo)%>%
  group_by(id,geoarea)%>%
  mutate(totalwt=sum(pyamtw),
         totalv =sum(pyamtv),
         nDiets =n_distinct(dietID))%>%
  group_by(id,geoarea,gl_prey)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets))%>%
  dplyr::select(geoarea,year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,trawlabundance=wholetrawlabundance,gl_prey,totalwt,totalv,nDiets,qikw,qikv,pik)%>% 
  distinct() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.


#Means for a whole geoarea
nDiets_geoy<-trawlDiets_geo%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(year,season,geoarea,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(year,season,geoarea)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun_geoy<-trawlDiets_geo%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(year,season,geoarea,id,trawlabundance)%>%distinct()%>%
  group_by(year,season,geoarea)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_geoy<-trawlDiets_geo%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(year,season,geoarea,id,trawlabundance)%>%distinct()%>%
  group_by(year,season,geoarea)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietSum_geoy<-trawlDiets_geo%>%
  group_by(year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_geoy)%>%
  left_join(sumAbun_geoy)%>%
  left_join(sumAbun_nEmpty_geoy)%>%
  group_by(year,season,geoarea,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       trawlabundance,pik, Fk))%>%distinct()
trawlDietSum_geoy[is.na(trawlDietSum_geoy)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum_geoy%>%group_by(year,season,geoarea)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 (or they're all empties)


#Means for a whole geoarea
nDiets_geo<-trawlDiets_geo%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(season,geoarea,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(season,geoarea)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun_geo<-trawlDiets_geo%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(season,geoarea,id,trawlabundance)%>%distinct()%>%
  group_by(season,geoarea)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_geo<-trawlDiets_geo%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(season,geoarea,id,trawlabundance)%>%distinct()%>%
  group_by(season,geoarea)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietSum_geo<-trawlDiets_geo%>%
  group_by(season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_geo)%>%
  left_join(sumAbun_geo)%>%
  left_join(sumAbun_nEmpty_geo)%>%
  group_by(season,geoarea,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       trawlabundance,pik, Fk))%>%distinct()
trawlDietSum_geo[is.na(trawlDietSum_geo)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum_geo%>%group_by(season,geoarea)%>%summarise(N=sum(Wk))
table(check$N) #When this is 0 or 1, then you have all your means correct because they sum to 1 or they're all empties








#Need to have whole diet totals, and whole trawl measures
trawlDiets_sizes<-preyGMRI_filter%>%
  group_by(id,comname)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv))%>%
  group_by(id,pdscinam,gl_prey)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets))%>%
  dplyr::select(geoarea,year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,sizeabundance,pdscinam,comname,sizecat=sizecat2,gl_prey,totalwt,totalv,nDiets,qikw,qikv,pik)%>% 
  distinct() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.


nDiets_sizes<-trawlDiets_sizes%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(geoarea,year,season,comname,sizecat,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(geoarea,year,season,comname,sizecat)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun_sizes<-trawlDiets_sizes%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(geoarea,year,season,comname,sizecat,id,sizeabundance)%>%distinct()%>%
  group_by(geoarea,year,season,comname,sizecat)%>%summarise(sumAbun=sum(sizeabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_sizes<-trawlDiets_sizes%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(geoarea,year,season,comname,sizecat,id,sizeabundance)%>%distinct()%>%
  group_by(geoarea,year,season,comname,sizecat)%>%summarise(sumAbun_nEmpty=sum(sizeabundance))


trawlDietSum_sizes<-trawlDiets_sizes%>%
  group_by(geoarea,year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(geoarea,year,season,pdscinam,comname,sizecat,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(sizeabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(sizeabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(sizeabundance*pik) /sumAbun, #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizeabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizeabundance,qikw,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       sizeabundance,pik, Fk))%>%distinct()
trawlDietSum_sizes[is.na(trawlDietSum_sizes)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum_sizes%>%group_by(geoarea,year,season,pdscinam)%>%summarise(N=sum(Wk))
table(check$N) #you have all your means correct because they sum to 1 (or they're all empties)

#Saving these, so they can be better used within other analyses
#write.csv(trawlDiets,"trawl_speciesClusterCompositions.csv",row.names=F)
#write.csv(trawlDietSum,"geoareayearseason_speciesClusterCompositions.csv",row.names=F)





# Emptiness over time -----------------------------------------------------

empty<-filter(trawlDietSum,gl_prey=="Empty")%>%
  full_join(allGeoCom)%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))
empty_ny<-filter(trawlDietSum_ny,gl_prey=="Empty")%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=levels(empty$comname)))
empty_sp<-filter(trawlDietSum_sp, gl_prey=="Empty")%>%
  mutate(comname=factor(comname,levels=levels(empty$comname)))
empty_spy<-filter(trawlDietSum_spy,gl_prey=="Empty")%>%
  mutate(comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))
  #mutate(comname=factor(comname,levels=levels(empty$comname)))
empty_geo<-filter(trawlDietSum_geo,gl_prey=="Empty")%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))
empty_geoy<-filter(trawlDietSum_geoy,gl_prey=="Empty")%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))


#Over the years for each species
s<-ggplot()+
  #geom_rect(data=empty_sp,
  #            aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_sp,
  #            aes(slope=0,intercept=Fk*100,color=season,lty=season),size=2)+
  geom_line(data=empty_spy,
            aes(year,Fk,lty=season),size=2)+
  geom_errorbar(data=empty_spy,
                aes(year,ymin=Fk-FkSD,ymax=Fk+FkSD,color=season),width=0.2,size=1.1)+
  geom_point(data=empty_spy,
             aes(year,Fk,fill=comname,color=season),size=4,shape=21,stroke=0.75)+
  geom_smooth(data=empty_spy,
              aes(year,Fk,color=season,lty=season),size=3,alpha=0.5,method="lm",show.legend = F)+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Proportion Empty Stomachs")+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_viridis_d(option="C",end=0.8)+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  guides(shape="none",color="none",fill="none",
         linetype=guide_legend(override.aes=list(size=5,lty=c(1,3),color=seasonPal2[c(2,4)])))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=30),
        legend.key.width=unit(66,"pt"),
        axis.title=element_blank(),axis.text.y=element_blank(),
        plot.margin=unit(c(12.5,15,15,12.5),"pt"))+
  facet_grid(comname~.)
s


#Calculating the drop for each species:season
summary(s_eLM<-lm(Fk~year*season*comname,data=empty_spy))
preds_s_e<-data.frame(empty_spy,pred=predict(s_eLM))%>%
  group_by(season,comname)%>%
  summarise(drop=(max(pred)-min(pred))*100,
            direction=ifelse(lag(pred)<pred,"+","-"))%>%
  filter(!is.na(direction))%>%distinct()



#Over the years for each area
g<-ggplot()+
  #geom_rect(data=empty_geo,
  #            aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_geo,
  #            aes(slope=0,intercept=Fk*100,lty=season,color=season),size=2)+
  geom_line(data=empty_geoy,
            aes(year,Fk,lty=season),color="black",size=2)+
  geom_errorbar(data=empty_geoy,
              aes(year,ymin=Fk-FkSD,ymax=Fk+FkSD,color=season),size=1.1)+
  geom_point(data=empty_geoy,
             aes(year,Fk,shape=geoarea,color=season),fill="black",size=4,stroke=0.75)+
  geom_smooth(data=empty_geoy,
              aes(year,Fk,lty=season,color=season),method="lm",alpha=0.5,size=3)+
  scale_x_continuous(name="Year",breaks=c(1980,2000,2020),limits=c(1973,2019))+
  scale_y_continuous(name="",limits=c(0,1))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_manual(values=seasonPal2[c(2,4)])+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  scale_shape_manual(values=c(21:25))+
  theme(axis.title.y=element_text(size=40),
        plot.margin=unit(c(12.5,15,15,12.5),"pt"),
        legend.position="none")+
  facet_wrap(~geoarea,nrow=1,strip.position = "bottom")
g


#Calculating the drop for each species:season
summary(g_eLM<-lm(Fk~year*season*geoarea,data=empty_geoy))
preds_g_e<-data.frame(empty_geoy,pred=predict(g_eLM))%>%
  group_by(season,geoarea)%>%
  summarise(drop=(max(pred)-min(pred))*100,
            direction=ifelse(lag(pred)<pred,"+","-"))%>%
  filter(!is.na(direction))%>%distinct()





#Over the years for each species in each area
a<-ggplot()+
  #geom_rect(data=empty_ny,
  #          aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_ny,
  #            aes(slope=0,intercept=Fk*100,lty=season,color=season),size=2)+
  geom_line(data=empty,
            aes(year,Fk,lty=season),size=2,show.legend = F)+
  geom_errorbar(data=empty,
                aes(year,ymin=Fk-FkSD,ymax=Fk+FkSD,color=season),size=1.1)+
  geom_point(data=empty,
             aes(year,Fk,fill=comname,color=season,shape=geoarea),size=4,stroke=0.75)+
  geom_smooth(data=empty,
              aes(year,Fk,lty=season,color=season),method="lm",alpha=0.5,size=3)+
  scale_x_continuous(name="Year",breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Proportion Empty Stomachs",limits=c(0,1))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_viridis_d(option="C",end=0.8)+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  scale_shape_manual(values=c(21:25))+
  theme(plot.margin=unit(c(12.5,15,15,12.5),"pt"),legend.position = "none",
        axis.title.y=element_text(size=40),axis.title.x=element_blank(),
        strip.text=element_blank())+
  facet_wrap(~comname+geoarea,nrow=9,labeller = label_wrap_gen(multi_line=FALSE))
a

#Put them all together
grid.arrange(a,s,g,nrow=2,layout_matrix=rbind(c(1,1,1,1,1,2,2),
                                              c(1,1,1,1,1,2,2),
                                              c(1,1,1,1,1,2,2),
                                              c(1,1,1,1,1,2,2),
                                              c(3,3,3,3,3,NA,NA)))



# Composition -----------------------------------------------------------

#A prettier order of the gl_preycats so that they can be better understood in a colorscheme
unique(gl_preycats$matchingCats)
#The fish:
  #Ammodytes sp, Clupeidae, Clupea harengus, Cottidae, Engraulidae, Fish eggs, Fish larvae, Gadiformes, Illex sp, 
  #Lepophidium profundorum, Loligo sp, Macrozoarces americanus, Merluccius bilinearis, Other hakes, Other fish, Peprilus triacanthus, 
  #Pleuronectiformes, Rajiformes, Scombridae, Unidentified fish
#Benthic:
  #Anthozoa, Bivalvia, Cancridae, Crangonidae, Crustacean shrimp, Crustacea, Decapoda crab, Decapoda shrimp, Echinodermata,
  #Gastropoda, Holothuroidea, Hydrozoa, Isopoda, Mollusca, Ophiuroidea, Paguroidea, Polychaeta, Worms
#Pelagic:
  #Amphipoda, Cephalopoda, Cnidaria, Ctenophora, Euphausiidae, Gammaridea, Mysidacea, Other invertebrates, Pandalidae, Zooplankton
#Other:
  #Animal remains, Miscellaneous, Other
trawlDietSum_sp<-trawlDietSum_sp%>%
  mutate(gl_prey=factor(gl_prey,levels=c("Clupea harengus","Lepophidium profundorum","Macrozoarces americanus","Merluccius bilinearis", 
                                         "Peprilus triacanthus","Ammodytes sp","Clupeidae", "Cottidae", "Engraulidae",
                                         "Scombridae","Pleuronectiformes", "Rajiformes" ,"Gadiformes","Other hakes",
                                         "Other fish","Fish eggs", "Fish larvae", "Unidentified fish", "Illex sp","Loligo sp",
                                         "Cephalopoda", "Cnidaria", "Ctenophora", "Euphausiidae","Gammaridea","Hyperiidae","Pandalidae", 
                                         "Mysidacea", "Cancridae", "Crangonidae","Crustacean shrimp","Crustacea", "Decapoda crab", "Decapoda shrimp",
                                         "Amphipoda","Anthozoa", "Bivalvia","Echinodermata", "Gastropoda","Holothuroidea", 
                                         "Hydrozoa","Isopoda","Mollusca", "Ophiuroidea", "Paguroidea", "Polychaeta", "Worms",
                                         "Zooplankton", "Other invertebrates", "Animal remains", "Miscellaneous", "Other",
                                         "Unobserved","Empty")),
         comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))
  

ggplot(data=filter(trawlDietSum_ny,gl_prey!="Empty"))+
  geom_col(aes(comname,Wk,fill=gl_prey),color="black",position="stack")+
  theme(axis.text.x=element_text(angle=30,hjust=0.8,vjust=1,size=15),
        legend.position="top")+
  scale_x_discrete(name="Species Common Name")+
  scale_y_continuous(name="Mean Proportion of Diet Mass",
                     expand=expansion(add=c(0.01,0.01)))+
  scale_fill_manual(values=colorspace::rainbow_hcl(n=n_distinct(filter(trawlDietSum_ny,gl_prey!="Empty")$gl_prey),l=69,c=100),
                    name="Prey Category")+
  facet_grid(season~factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))

#Short and long for poster
ggplot(data=filter(trawlDietSum_sp,gl_prey%notin%c("Empty")))+
  geom_col(aes(comname,Wk,fill=gl_prey),color="black",position="stack")+
  theme(axis.text.x=element_text(angle=20,hjust=0.9,vjust=1,size=30),
        legend.position="right")+
  guides(fill=guide_legend(nrow=25))+
  scale_x_discrete(name="Species Common Name")+
  scale_y_continuous(name="Mean Proportion of Diet Mass",
                     expand=expansion(add=c(0.01,0.01)))+
  scale_fill_manual(values=colorspace::rainbow_hcl(n=n_distinct(filter(trawlDietSum_ny,gl_prey!="Empty")$gl_prey),l=69,c=100),
                    name="Prey Category")+
  facet_wrap(season~.,nrow=1)
#Second version for poster
ggplot(data=filter(trawlDietSum_sp,gl_prey%notin%c("Empty")))+
  geom_col(aes(paste(comname,season),Wk,fill=gl_prey),color="black",position="stack")+
  theme(axis.text.x=element_text(angle=20,hjust=0.9,vjust=1,size=30),
        legend.position="right")+
  guides(fill=guide_legend(nrow=25))+
  scale_x_discrete(name="Species Common Name")+
  scale_y_continuous(name="Mean Proportion of Diet Mass",
                     expand=expansion(add=c(0.01,0.01)))+
  scale_fill_manual(values=colorspace::rainbow_hcl(n=n_distinct(filter(trawlDietSum_ny,gl_prey!="Empty")$gl_prey),l=69,c=100),
                    name="Prey Category")

#Simplified, remove anyone who has less than a total of 1% averaged across all panels
rarePrey_ny<-trawlDietSum_ny%>%
  group_by(gl_prey)%>%
  summarise(meanWk=mean(Wk))%>%
  filter(meanWk>0.01)
rareTrawlDietSum_ny<-filter(trawlDietSum_ny,gl_prey%in%rarePrey_ny$gl_prey)

ggplot(data=rareTrawlDietSum_ny)+
  geom_col(aes(comname,Wk,fill=gl_prey),color="black",position="stack")+
  theme(axis.text.x=element_text(angle=30,hjust=0.8,vjust=1,size=15),
        legend.position="top")+
  scale_x_discrete(name="Species Common Name")+
  scale_y_continuous(name="Mean Proportion of Diet Mass",
                     expand=expansion(add=c(0.01,0.01)))+
  scale_fill_manual(values=colorspace::rainbow_hcl(n=n_distinct(rareTrawlDietSum_ny$gl_prey),l=69,c=100),
                    name="Prey Category")+
  facet_grid(season~factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))





#And with Volume, probably looks the same
ggplot(data=filter(trawlDietSum_ny,gl_prey!="Empty"))+
  geom_col(aes(comname,Vk,fill=gl_prey),color="black",position="stack")+
  theme(axis.text.x=element_text(angle=30,hjust=0.8,vjust=1,size=15),
        legend.position="top")+
  scale_x_discrete(name="Species Common Name")+
  scale_y_continuous(name="Mean Proportion of Diet Volume",
                     expand=expansion(add=c(0.01,0.01)))+
  scale_fill_manual(values=colorspace::rainbow_hcl(n=n_distinct(filter(trawlDietSum_ny,gl_prey!="Empty")$gl_prey),l=69,c=100),
                    name="Prey Category")+
  facet_grid(season~factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))


#Simplified, remove anyone who has less than a total of 1% averaged across all panels
rarePrey_ny<-trawlDietSum_ny%>%
  group_by(gl_prey)%>%
  summarise(meanVk=mean(Vk))%>%
  filter(meanVk>0.01)
rareTrawlDietSum_ny<-filter(trawlDietSum_ny,gl_prey%in%rarePrey_ny$gl_prey)

ggplot(data=rareTrawlDietSum_ny)+
  geom_col(aes(comname,Vk,fill=gl_prey),color="black",position="stack")+
  theme(axis.text.x=element_text(angle=30,hjust=0.8,vjust=1,size=15),
        legend.position="top")+
  scale_x_discrete(name="Species Common Name")+
  scale_y_continuous(name="Mean Proportion of Diet Volume",
                     expand=expansion(add=c(0.01,0.01)))+
  scale_fill_manual(values=colorspace::rainbow_hcl(n=n_distinct(rareTrawlDietSum_ny$gl_prey),l=69,c=100),
                    name="Prey Category")+
  facet_grid(season~factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))





# Levin's Breadth ---------------------------------------------------------


#Calculated from the Wk in each GROUP
#Can't have any variance, it will be a single value


matDietSum_spy_Wk<-trawlDietSum_spy%>%
  mutate(gl_prey=paste("Prey",gl_prey,sep="."),
         Wk2=Wk^2)%>%
  pivot_wider(id_cols=c(year,season,comname),
              names_from="gl_prey",values_from="Wk2")%>%
  ungroup()%>%
  dplyr::select(-Prey.Empty)%>%
  mutate(levin=1/select(.,starts_with("Prey"))%>%rowSums(na.rm=T),
         levin=ifelse(is.infinite(levin),NA,levin),
         levinStd=(levin-1)/(n_distinct(trawlDietSum$gl_prey)-2))%>%
  mutate(comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))
  #mutate(comname=fct_reorder(factor(comname,levels=unique(trawlDietSum$comname)),levin,na.rm=T))

s_l<-ggplot()+
  #geom_rect(data=empty_sp,
  #          aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_sp,
  #            aes(slope=0,intercept=Fk*100,color=season,lty=season),size=2)+
  geom_line(data=matDietSum_spy_Wk,
            aes(year,levinStd,lty=season),size=2)+
  geom_point(data=matDietSum_spy_Wk,
             aes(year,levinStd,fill=comname,color=season),size=4,shape=21,stroke=0.75)+
  geom_smooth(data=matDietSum_spy_Wk,
              aes(year,levinStd,lty=season,color=season),method="lm",size=3,alpha=0.5,show.legend=F)+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Levin's Standardized Breadth Index",limits=c(0,0.28))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_viridis_d(option="C",end=0.8)+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  guides(shape="none",color="none",fill="none",
         linetype=guide_legend(override.aes=list(size=5,lty=c(1,3),color=seasonPal2[c(2,4)])))+
  theme(legend.title=element_blank(),plot.margin=unit(c(12.5,15,15,12.5),"pt"),
        legend.text=element_text(size=30),legend.key.width=unit(66,"pt"),
        axis.title=element_blank())+
  facet_grid(comname~.)
s_l

#Calculating the drop for each species:season
summary(s_lLM<-lm(levinStd~year*season*comname,data=matDietSum_spy_Wk))
preds_s_l<-data.frame(filter(matDietSum_spy_Wk,!is.na(levinStd)),pred=predict(s_lLM))%>%
  group_by(season,comname)%>%
  summarise(drop=(max(pred)-min(pred)),
            direction=ifelse(lag(pred)<pred,"+","-"))%>%
  filter(!is.na(direction))%>%distinct()



matDietSum_geo_Wk<-trawlDietSum_geoy%>%
  mutate(gl_prey=paste("Prey",gl_prey,sep="."),
         Wk2=Wk^2)%>%
  pivot_wider(id_cols=c(geoarea,year,season),
              names_from="gl_prey",values_from="Wk2")%>%
  ungroup()%>%
  dplyr::select(-Prey.Empty)%>%
  mutate(levin=1/select(.,starts_with("Prey"))%>%rowSums(na.rm=T),
         levin=ifelse(is.infinite(levin),NA,levin),
         levinStd=(levin-1)/(n_distinct(trawlDietSum$gl_prey)-2),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))


#Over the years for each area
g_l<-ggplot()+
  #geom_rect(data=empty_geo,
  #          aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_geo,
  #            aes(slope=0,intercept=Fk*100,lty=season,color=season),size=2)+
  #geom_ribbon(data=empty_geoy,
  #            aes(year,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.4,size=2)+
  geom_line(data=matDietSum_geo_Wk,
            aes(year,levinStd,lty=season),color="black",size=2)+
  geom_point(data=matDietSum_geo_Wk,
             aes(year,levinStd,shape=geoarea,color=season),fill="black",size=4,stroke=0.75)+
  geom_smooth(data=matDietSum_geo_Wk,
              aes(year,levinStd,lty=season,color=season),method="lm",size=3,alpha=0.5,show.legend=F)+
  scale_x_continuous(name="Year",breaks=c(1980,2000,2020))+
  scale_y_continuous(name="",limits=c(0,0.28))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_manual(values=seasonPal2[c(2,4)])+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  scale_shape_manual(values=c(21:25))+
  theme(legend.position="none",plot.margin=unit(c(12.5,15,15,12.5),"pt"),
        axis.title.y=element_text(size=40))+
  facet_wrap(~geoarea,nrow=1,strip.position = "bottom")
g_l

#Calculating the drop for each geoarea:season
summary(g_lLM<-lm(levinStd~year*season*geoarea,data=matDietSum_geo_Wk))
preds_g_l<-data.frame(matDietSum_geo_Wk,pred=predict(g_lLM))%>%
  group_by(season,geoarea)%>%
  summarise(drop=(max(pred)-min(pred)))




matDietSum_Wk<-trawlDietSum%>%
  mutate(gl_prey=paste("Prey",gl_prey,sep="."),
         Wk2=Wk^2)%>%
  pivot_wider(id_cols=c(geoarea,year,season,comname),
                           names_from="gl_prey",values_from="Wk2")%>%
  ungroup()%>%
  dplyr::select(-Prey.Empty)%>%
  mutate(levin=1/select(.,starts_with("Prey"))%>%rowSums(na.rm=T),
         levin=ifelse(is.infinite(levin),NA,levin),
         levinStd=(levin-1)/(n_distinct(trawlDietSum$gl_prey)-2))%>%
  full_join(allGeoCom)%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=levels(matDietSum_spy_Wk$comname)))


a_l<-ggplot()+
  #geom_rect(data=empty_ny,
  #          aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_ny,
  #            aes(slope=0,intercept=Fk*100,lty=season,color=season),size=2)+
  #geom_ribbon(data=empty,
  #            aes(year,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.4,size=2)+
  #scale_fill_manual(values=seasonPal2[c(2,4)])+
  #ggnewscale::new_scale_fill()+
  geom_line(data=matDietSum_Wk,
            aes(year,levinStd,lty=season),size=2,show.legend = F)+
  geom_point(data=matDietSum_Wk,
             aes(year,levinStd,fill=comname,color=season,shape=geoarea),size=4,stroke=0.75)+
  geom_smooth(data=matDietSum_Wk,
              aes(year,levinStd,lty=season,color=season),method="lm",size=3,alpha=0.5,show.legend=F)+
  scale_x_continuous(name="Year",breaks=c(1980,2000,2020),labels=c("","",""))+
  scale_y_continuous(name="Levin's Standardized Breadth Index",limits=c(0,0.28))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_viridis_d(option="C",end=0.8)+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  scale_shape_manual(values=c(21:25))+
  theme(plot.margin=unit(c(12.5,15,15,12.5),"pt"),legend.position = "none",
        axis.title.y=element_text(size=40),axis.title.x=element_blank(),
        strip.text=element_blank())+
  facet_wrap(~comname+geoarea,nrow=9,labeller = label_wrap_gen(multi_line=FALSE))
a_l


#Put them all together
grid.arrange(a_l,s_l,g_l,nrow=2,layout_matrix=rbind(c(1,1,1,1,1,2,2),
                                                    c(1,1,1,1,1,2,2),
                                                    c(1,1,1,1,1,2,2),
                                                    c(1,1,1,1,1,2,2),
                                                    c(3,3,3,3,3,NA,NA)))



#Very simple lm, just to check feasibility
summary(levinLM<-lm(levinStd~year+comname+geoarea+season,data=matDietSum_Wk))


# Relative Consumption ----------------------------------------------------

#Has to be a single value for each individual, so remove the different prey
indGMRI_filter #That's in here from above

#Not sure how much I trust pdwgt, but it actually looks pretty good with only a couple possible errorenous
ggplot(indGMRI_filter)+
  geom_point(aes(pdlen,pdwgt))+
  geom_smooth(aes(pdlen,pdwgt))+
  ylab("Mass (g)")+
  xlab("Length (cm)")+
  facet_wrap(~comname)

trawlDiets_ind<-indGMRI_filter%>%
  group_by(id,comname)%>%
  mutate(relConsump=pdgutw/pdwgt,
         relConsump=ifelse(!is.finite(relConsump),NA,relConsump),
         meanC=mean(relConsump,na.rm=T),
         meanWt=mean(pdgutw),
         meanV =mean(pdgutv))%>%
  dplyr::select(geoarea,year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,trawlabundance,comname,meanWt,meanV,meanC,nDiets)%>% 
  distinct()%>%ungroup() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.

#c means that's only the stomachs with something in them (c=consumed)
trawlDiets_cind<-indGMRI_filter%>%
  filter(pdgutw>0&pdwgt>0)%>%
  group_by(id,comname)%>%
  mutate(relConsump=pdgutw/pdwgt,
         meanC=mean(relConsump,na.rm=T),
         meanWt=mean(pdgutw),
         meanV =mean(pdgutv),
         nDiets=n_distinct(dietID))%>%
  dplyr::select(geoarea,year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,trawlabundance,comname,meanWt,meanV,meanC,nDiets)%>% 
  distinct()%>%ungroup() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.


#The abundance totals are the same for the full diet analysis above

trawlDietSum_ind<-trawlDiets_ind%>%
  group_by(geoarea,year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  group_by(year,season,geoarea,comname,nTrawls,tnDiets,sumAbun)%>%
  summarise(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*meanV)/sumAbun, #Same
            Ck=sum(trawlabundance*meanC,na.rm=T)/sumAbun, #This is relative consumption, so there might be NAs if a fish wasn't weighed
            WkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanWt,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanV, Vk),
            CkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanC, Ck))%>%distinct()%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         season=factor(season,levels=c("Spring","Fall")))%>%ungroup()%>%
  full_join(allGeoCom)%>%
  mutate(comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))

#Averages for the whole species
ggplot(trawlDietSum_ind,aes(comname,Ck*100))+
  geom_boxplot()+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),labels=c("0.001","0.01","0.1","1","10"))

trawlDietSum_cind<-trawlDiets_cind%>%
  group_by(geoarea,year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  group_by(year,season,geoarea,comname,nTrawls,tnDiets,sumAbun)%>%
  summarise(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*meanV)/sumAbun, #Same
            Ck=sum(trawlabundance*meanC,na.rm=T)/sumAbun, #This is relative consumption, so there might be NAs if a fish wasn't weighed
            WkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanWt,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanV, Vk),
            CkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanC, Ck))%>%distinct()%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         season=factor(season,levels=c("Spring","Fall")))%>%ungroup()%>%
  full_join(allGeoCom)%>%
  mutate(comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))


trawlDietSum_cind_spy<-trawlDiets_cind%>%
  group_by(year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_spy)%>%
  left_join(sumAbun_spy)%>%
  group_by(year,season,comname,nTrawls,tnDiets,sumAbun)%>%
  summarise(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*meanV)/sumAbun, #Same
            Ck=sum(trawlabundance*meanC,na.rm=T)/sumAbun, #This is relative consumption, so there might be NAs if a fish wasn't weighed
            WkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanWt,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanV, Vk),
            CkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanC, Ck))%>%distinct()%>%
  mutate(season=factor(season,levels=c("Spring","Fall")))%>%ungroup()%>%
  mutate(comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))
  #mutate(comname=fct_reorder(factor(comname,levels=unique(trawlDietSum_ind$comname)),Ck,na.rm=T))


trawlDietSum_cind_geoy<-trawlDiets_cind%>%
  group_by(geoarea,year,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_geoy)%>%
  left_join(sumAbun_geoy)%>%
  group_by(year,season,geoarea,nTrawls,tnDiets,sumAbun)%>%
  summarise(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*meanV)/sumAbun, #Same
            Ck=sum(trawlabundance*meanC,na.rm=T)/sumAbun, #This is relative consumption, so there might be NAs if a fish wasn't weighed
            WkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanWt,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanV, Vk),
            CkSD=clusterSD(nTrawls,sumAbun,trawlabundance,meanC, Ck))%>%distinct()%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         season=factor(season,levels=c("Spring","Fall")))%>%ungroup()


a_c<-ggplot()+
  #geom_rect(data=empty_ny,
  #          aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_ny,
  #            aes(slope=0,intercept=Fk*100,lty=season,color=season),size=2)+
  geom_line(data=trawlDietSum_cind,
            aes(year,Ck,lty=season),size=2,show.legend = F)+
  geom_errorbar(data=trawlDietSum_cind,
                aes(year,ymin=Ck-CkSD,ymax=Ck+CkSD,color=season),size=1.1)+
  geom_point(data=trawlDietSum_cind,
             aes(year,Ck,fill=comname,color=season,shape=geoarea),size=4,stroke=0.75)+
  geom_smooth(data=trawlDietSum_cind,
              aes(year,Ck,lty=season,color=season),method="lm",size=3,alpha=0.5,show.legend=F)+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020),labels=c("","",""))+
  scale_y_log10(name="Relative Consumption (g/g)",
                breaks=c(0.001,0.01,0.1),
                labels=c("0.001","0.01","0.1"),
                limits=c(0.001,0.1))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_viridis_d(option="C",end=0.8)+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  scale_shape_manual(values=c(21:25))+
  theme(plot.margin=unit(c(12.5,15,15,12.5),"pt"),legend.position = "none",
        axis.title.y=element_text(size=40),axis.title.x=element_blank(),
        strip.text=element_blank())+
  facet_wrap(~comname+geoarea,nrow=9,labeller = label_wrap_gen(multi_line=FALSE))




g_c<-ggplot()+
  #geom_rect(data=empty_geo,
  #          aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_geo,
  #            aes(slope=0,intercept=Fk*100,lty=season,color=season),size=2)+
  #geom_ribbon(data=empty_geoy,
  #            aes(year,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.4,size=2)+
  geom_line(data=trawlDietSum_cind_geoy,
            aes(year,Ck,lty=season),color="black",size=2)+
  geom_errorbar(data=trawlDietSum_cind_geoy,
                aes(year,ymin=Ck-CkSD,ymax=Ck+CkSD,color=season),size=1.1)+
  geom_point(data=trawlDietSum_cind_geoy,
             aes(year,Ck,shape=geoarea,color=season),fill="black",size=4,stroke=0.75)+
  geom_smooth(data=trawlDietSum_cind_geoy,
              aes(year,Ck,lty=season,color=season),method="lm",size=3,alpha=0.5,show.legend=F)+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_log10(name="",
                breaks=c(0.001,0.01,0.1),
                labels=c("0.001","0.01","0.1"),
                limits=c(0.001,0.1))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_manual(values=seasonPal2[c(2,4)])+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  scale_shape_manual(values=c(21:25))+
  theme(legend.position="none",axis.title=element_text(size=40),
        plot.margin=unit(c(12.5,15,15,12.5),"pt"))+
  facet_wrap(~geoarea,nrow=1,strip.position="bottom")
g_c

#Calculating the drop for each geoarea:season
summary(g_cLM<-lm(Ck~year*season*geoarea,data=trawlDietSum_cind_geoy))
preds_g_c<-data.frame(trawlDietSum_cind_geoy,pred=predict(g_cLM))%>%
  group_by(season,geoarea)%>%
  summarise(drop=(max(pred)-min(pred))*100)



s_c<-ggplot()+
  #geom_rect(data=empty_sp,
  #          aes(xmin=-Inf,xmax=Inf,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.2,size=2)+
  #geom_abline(data=empty_sp,
  #            aes(slope=0,intercept=Fk*100,color=season,lty=season),size=2)+
  #geom_ribbon(data=empty_spy,
  #            aes(year,ymin=(Fk-1.96*FkSD)*100,ymax=(Fk+1.96*FkSD)*100,fill=season),alpha=0.4,size=2)+
  #scale_fill_manual(values=seasonPal2[c(2,4)])+
  #ggnewscale::new_scale_fill()+
  geom_line(data=trawlDietSum_cind_spy,
            aes(year,Ck,lty=season),size=2)+
  geom_point(data=trawlDietSum_cind_spy,
             aes(year,Ck,fill=comname,color=season),size=4,shape=21,stroke=0.75)+
  geom_smooth(data=trawlDietSum_cind_spy,
              aes(year,Ck,lty=season,color=season),method="lm",size=3,alpha=0.5,show.legend=F)+
  scale_x_continuous(name=" ",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_log10(name="Relative Consumption (g/g)",
                breaks=c(0.001,0.01,0.1),
                labels=c("0.001","0.01","0.1"),
                limits=c(0.001,0.1))+
  #scale_size_continuous(name="Number of\nDiets",range=c(1.5,7))+
  scale_fill_viridis_d(option="C",end=0.8)+
  scale_color_manual(values=seasonPal2[c(2,4)])+
  guides(shape="none",color="none",fill="none",
         linetype=guide_legend(override.aes=list(size=5,lty=c(1,3),color=c(seasonPal2[c(2,4)]))))+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=30),legend.key.width=unit(66,"pt"),
        axis.title=element_blank(),axis.text.y=element_blank(),
        plot.margin=unit(c(12.5,15,15,12.5),"pt"))+
  facet_grid(comname~.)
s_c

#Calculating the drop for each species:season
summary(s_cLM<-lm(Ck~year*season*comname,data=trawlDietSum_cind_spy))
preds_s_c<-data.frame(trawlDietSum_cind_spy,pred=predict(s_cLM))%>%
  group_by(season,comname)%>%
  summarise(drop=(max(pred)-min(pred))*100)



#Put them all together
grid.arrange(a_c,s_c,g_c,nrow=2,layout_matrix=rbind(c(1,1,1,1,1,2,2),
                                                    c(1,1,1,1,1,2,2),
                                                    c(1,1,1,1,1,2,2),
                                                    c(1,1,1,1,1,2,2),
                                                    c(3,3,3,3,3,NA,NA)))





#Poster panels
grid.arrange(s,s_l,s_c,nrow=1) #Without legend, don't want unbalanced matrix
grid.arrange(g,g_l,g_c,nrow=3)

# Prey Length -------------------------------------------------------------

pylenTrawl<-pylen19%>%
  mutate(svspp=str_pad(svspp,width=3,pad="0",side="left"),
         pdid=str_pad(pdid,width=6,pad="0",side="left"),
         pdsex=as.character(pdsex))%>%
  left_join(indGMRI_filter)%>%
  mutate(pypdRatio=(pylen/10)/pdlen,
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))%>%
  filter(!is.na(trawlabundance))

pylenTrawlmean<-pylenTrawl%>%
  group_by(id,comname)%>%
  mutate(meanPyPd=mean(pypdRatio,na.rm=T))%>%
  dplyr::select(geoarea,season,
                id,trawlabundance,comname,meanPyPd)%>%
  distinct()

pylenSum<-pylenTrawlmean%>%
  group_by(geoarea,season)%>%
  mutate(nTrawls=n_distinct(id))%>%
  group_by(geoarea,season,comname)%>%
  summarise(PyPd=sum(trawlabundance*meanPyPd)/sum(trawlabundance),
            skew=skewness(meanPyPd),
            m=max(meanPyPd))%>%distinct()

ggplot()+
  geom_violin(data=pylenTrawlmean,aes(geoarea,meanPyPd,fill=season),
              position=position_dodge(width=0.9))+
  geom_jitter(data=pylenTrawlmean,aes(geoarea,meanPyPd,fill=season),
              shape=21,position=position_jitterdodge(dodge.width=0.9,jitter.width=0.1),alpha=0.6)+
  geom_errorbar(data=pylenSum,aes(geoarea,ymin=PyPd,ymax=PyPd,fill=season),
                size=1.25,position=position_dodge(width=0.9),width=0.5)+
  geom_text(data=pylenSum,aes(geoarea,m+0.1,label=round(skew,digit=2),fill=season),
            position=position_dodge(width=0.9))+
  scale_color_viridis_d(option="C",end=0.95)+
  scale_fill_viridis_d(option="C",end=0.95)+
  scale_y_continuous(name="Prey:Predator Length Ratio")+
  scale_x_discrete(name="Geographic Area")+
  theme(legend.position="none",
        axis.text.x=element_text(angle=30,hjust=1,vjust=1.1))+
  facet_wrap(.~comname)


#All lms by species
pylen_cm<-pylenTrawl$pylen/10
pdlen<-pylenTrawl$pdlen
comname<-as.character(pylenTrawl$comname)
Xmat<-model.matrix(~pdlen*comname-1-pdlen)

summary(pylenLM<-lm(pylen_cm~Xmat))
LMcoefs<-data.frame(slopes=pylenLM$coefficients[11:19])
LMcoefs$comname=factor(gsub("Xmatpdlen:comname","",rownames(LMcoefs)),
                       levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                "Haddock","Pollock","Yellowtail flounder","Summer flounder"))

#How do the predators compare to their prey
ggplot(pylenTrawl)+
  geom_point(aes(pdlen,pylen/10,color=comname),shape=21,size=2,fill="transparent",alpha=0.7)+
  geom_abline(aes(slope=1,intercept=0))+
  geom_smooth(aes(pdlen,pylen/10),method="lm",lty=1,se=F,color="black",size=1.25)+
  geom_smooth(data=pylenTrawl%>%mutate(lengths=as.character(pdlen))%>%
                group_by(comname,lengths)%>%mutate(N=n())%>%filter(pylen==max(pylen)),
              aes(pdlen,pylen/10),method="lm",lty=3,se=F,color="black",size=1.5)+
  geom_smooth(data=pylenTrawl%>%mutate(lengths=as.character(pdlen))%>%
                group_by(comname,lengths)%>%mutate(N=n())%>%filter(pylen==min(pylen)),
              aes(pdlen,pylen/10),method="lm",lty=3,se=F,color="black",size=1.5)+
  geom_text(data=LMcoefs,aes(x=25,y=100,label=paste("Slope =",round(slopes,3)),color=comname),vjust=0)+
  scale_color_viridis_d(option="C",end=0.8,guide="none")+
  scale_x_continuous(limits=c(0,130),expand=expansion(mult=0.01),name="Predator length (cm)")+
  scale_y_continuous(limits=c(0,130),expand=expansion(mult=0.01),name="Prey length (cm)")+
  theme(legend.position=c(0.33,0.75),legend.background=element_rect(color="black",fill="white"))+
  facet_wrap(~comname)




# extras ------------------------------------------------------------------




#Looking at some filtering
preyTrawls_filtered<-filter(preyGMRI,
                            stratum >= 01010,
                            stratum <= 01760,
                            stratum != 1310,
                            stratum != 1320,
                            stratum != 1330,
                            stratum != 1350,
                            stratum != 1410,
                            stratum != 1420,
                            stratum != 1490)%>%
  mutate(
    strata = str_pad(stratum, width = 5, pad = "0", side = "left"),
    strat_num = str_sub(strata, 3,4),
    survey_area =  dplyr::case_when(
      strat_num %in% strata_key$`Georges Bank`         ~ "GB",
      strat_num %in% strata_key$`Gulf of Maine`        ~ "GoM",
      strat_num %in% strata_key$`Southern New England` ~ "SNE",
      strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
      TRUE                                             ~ "stratum not in key"))%>%
  filter(season%in%c("SPRING","FALL") & survey_area!="stratum not in key")

test<-filter(preyTrawls_filtered, id%notin%preyGMRI_filter$id)
n_distinct(test$dietID)
check<-filter(preyGMRI_filter,is.na(p))
check2<-filter(check,sizecat!=sizecat2)%>%
  select(svspp:sizecat2)%>%
  left_join(sumGMRItrawls)

n_distinct(prey19$dietID)






extras<-filter(preyGMRI,nDiets>trawlabundance)%>%
  dplyr::select(id,dietID,nDiets,trawlabundance)%>%distinct()
n_distinct(extras$id) #number of trawls that have more diets than abundance
n_distinct(extras$id)/n_distinct(prey19$id)*100 #percent

check<-prey19%>%filter(sizecat!=sizecat2)%>%
  select(svspp,pdsex,dietID,id,sizecat,sizecat2)


trawlDiets<-read_csv("SUMMARIZED/merged_prey_GMRItrawls.csv",
                     col_types = cols(ID = col_character()))%>%
  left_join(sumGMRItrawls)%>%
  group_by(ID,svspp,catchsex,sizecat)%>%
  mutate(nNAs=sum(is.na(p)))%>%
  group_by(ID,svspp)%>%
  mutate(tNAs=sum(is.na(p)),
         p=ifelse(is.na(p),nNAs/abundance,p),
         p=ifelse(!is.na(merge)&tNAs>0,p-(tNAs/abundance),p),
         abundance_sizecat=p*abundance)%>%ungroup()
nrow(trawlDiets[is.na(trawlDiets$p),])/nrow(trawlDiets)

test<-trawlDiets%>%
  group_by(ID,comname,abundance)%>%
  summarise(dietAbundance=n(),
            propDiet=dietAbundance/abundance)%>%distinct()
nrow(test[test$propDiet==1,])/nrow(test)*100
summary(test$propDiet)



