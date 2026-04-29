


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

library(R2jags)




#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}


seasonPal2<-c("deepskyblue4","yellowgreen","goldenrod1","orange3") #BEST

#Data loading
load("prey19.RData")
load("pylen19.RData")
load("googleMap.zoom6.eastCoast.R")
vulScores<-read.csv("../Hareetal_climateScores.csv")


#ALL LEVELS OF YEAR AND SEASON
ALLyearseasons<-factor(levels=c(paste(sort(rep(seq(1973,2019),4)),rep(c("Winter","Spring","Summer","Fall"),length(seq(1973,2019))))))
#levels(ALLyearseasons)

#The names of my species, just to have them nice and handy
keyspecies_sci<-unique(str_to_sentence(prey19$pdscinam))
keyspecies_com<-unique(str_to_title(prey19$pdcomnam))
keyspecies_svspp<-unique(str_pad(as.character(prey19$svspp),width=3,pad="0",side="left"))

#Size classes used by Garrison and Link
sizeClasses<-read.csv("GarrisonLink_predSizeClasses.csv") 
#These are replicated in the sizecats column already present, except that the small and largest cats extend to the smallest and largest individuals
prey19%>%group_by(pdcomnam,sizecat)%>%summarise(min=min(pdlen),max=max(pdlen))




# Species SVSPP Codes -----------------------------------------------------

spec<-read.csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv")%>%
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
load("NF.prey19")
uniquePrey19.all<-NF.prey19%>%
  dplyr::select(cruise6,station,svspp,pdsex,pdid,pdcomnam,pdscinam,pdlen,pdwgt,sizecat,pdgutw,pdgutv,declat,declon,month,day,year,season,geoarea)%>%
  distinct()%>%
  mutate(cruise6=as(cruise6,"character"),pdsex=as.character(pdsex),
         station=as.character(str_sub(paste0("0000",station),-4,-1)))



#The categories from Garrison and Link, 2000
gl_preycats<-read.csv("GarrisonLink_preyCats.csv")%>%
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



sizeClasses<-read.csv("GarrisonLink_predSizeClasses.csv")%>%
  left_join(spec,by=c("species_comnam"="LOGGED_SPECIES_NAME"))


allGMRItrawls_clean<-read_csv("../Trawl Data/NMFS Trawls/Complete/NMFS_survdat_gmri_tidy.csv",
                              col_types="cccccccnnncnnTnnnnnnnnnnccncn")%>%
  mutate(stratum_full=str_pad(stratum, width = 5, pad = "0", side = "left"),
         tow_full=str_pad(tow, width = 3, pad = "0", side = "left"),
         station_full=str_pad(station, width = 4, pad = "0", side = "left"),
         ID=paste0(cruise6,stratum_full,tow_full,station_full))%>%
  left_join(sizeClasses,by=c("svspp"="SVSPP"))%>%
  filter(svspp%in%keyspecies_svspp)%>%
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
  reframe(N_adj=sum(numlen_adj),p=N_adj/trawlabundance,
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
  distinct()# %>%expand(geoarea,comname)


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


empty<-filter(trawlDiets,gl_prey=="Empty")%>%
  right_join(dplyr::select(trawlDiets,geoarea:comname,nDiets)%>%distinct())%>%
  full_join(allGeoCom)%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")),
         pik=ifelse(is.na(pik),0,pik),
         freq=pik*nDiets,
         fullFreq=pik*trawlabundance)
empty2<-filter(empty,!is.na(year))

#c means that's only the stomachs with something in them (c=consumed)
trawlDiets_cind<-indGMRI_filter%>%
  filter(pdwgt>0 & pdgutw>0)%>% #Needs a mass, and only looking at those fish that ate (since emptiness is modeled elsewhere)
  group_by(id,comname)%>%
  mutate(relConsump=pdgutw/pdwgt,
         meanC=mean(relConsump,na.rm=T),
         meanWt=mean(pdgutw),
         meanV =mean(pdgutv),
         meanMass=mean(pdwgt),
         meanTL=mean(pdlen),
         nDiets=n_distinct(dietID),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))%>%
  dplyr::select(geoarea,year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,trawlabundance,comname,meanWt,meanV,meanC,meanMass,meanTL,nDiets)%>% 
  distinct()%>%ungroup() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.

hist(log(trawlDiets_cind$meanC))

matDietSum_Wk<-trawlDiets%>%
  #filter(gl_prey!="Empty")%>%
  mutate(prop2=qikw^2,
         gl_prey=paste0("Prey.",gl_prey))%>%
  pivot_wider(id_cols=c(geoarea,year,season,trawlabundance,totalwt,id,comname),
              names_from="gl_prey",values_from="prop2")%>%
  mutate(richness=rowSums(select(.,starts_with("Prey"))>0,na.rm=T),
         simpson=1-select(.,starts_with("Prey"))%>%rowSums(na.rm=T),
         shannon=-1*rowSums(select(.,starts_with("Prey"))*log(select(.,starts_with("Prey"))),na.rm=T),
         levin=1/select(.,starts_with("Prey"))%>%rowSums(na.rm=T),
         levin=ifelse(is.infinite(levin),NA,levin),
         levinStd=(levin-1)/(n_distinct(prey19$gl_prey)-2))%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")))
#mutate(comname=fct_reorder(factor(comname,levels=unique(trawlDietSum$comname)),levin,na.rm=T))
trawlDiets_Levin<-dplyr::select(matDietSum_Wk,geoarea:comname,richness:levinStd)%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")),
         comname=factor(comname,levels=c("Spiny dogfish","Little skate","White hake","Red hake","Silver hake",
                                         "Haddock","Pollock","Yellowtail flounder","Summer flounder")),
         simpson=ifelse(richness==0,NA,simpson),
         shannon=ifelse(richness==0,NA,shannon))

ggplot(trawlDiets_Levin,aes(levinStd,simpson))+
  geom_point()
simpson<-ggplot(trawlDiets_Levin,aes(comname,simpson))+
  geom_boxplot()+
  geom_jitter(width=0.2)
levin<-ggplot(trawlDiets_Levin,aes(comname,levinStd))+
  geom_boxplot()+
  geom_jitter(width=0.2)
grid.arrange(levin,simpson)

hist(trawlDiets_Levin$levinStd,breaks=100)

# Start of Bayesian Work --------------------------------------------------

# JAGS for Empty Stomachs -------------------------------------------------

sink("speciesEmpty_LM.txt")
cat("
    model {
    # PRIORS
    for (i in 1:nspecies) {
      beta.species[i] ~ dnorm(int.s.mean, int.s.pres)
      for (j in 1:ngeo) {
         beta.year[i,j] ~ dnorm(slope.mean,slope.pres)
      }
    }
    for (i in 1:ngeo) {
      beta.geo[i] ~ dnorm(int.r.mean, int.s.pres)
    }
    
    # Hyperpriors
      #Diffuse priors for geographic region random intercepts
    int.r.mean ~ dnorm(0,0.001)
    int.r.var ~ dunif(0,10)
    int.r.pres <- 1/int.r.var^2
    
      #Diffuse priors for species random intercepts
    int.s.mean ~ dnorm(0,0.001)
    int.s.var ~ dunif(0,10)
    int.s.pres <- 1/int.s.var^2
    
      #Diffuse priors for slope parameters
    slope.mean ~ dnorm(0,0.001)
    slope.var ~ dunif(0,10)
    slope.pres <- 1/slope.var^2
      
      #Diffuse prior for overdispersion term
    tau <- 1/sigma^2
    sigma ~ dunif(0,100)
    
    # LIKELIHOOD
    for (i in 1:ntrawls) {
      
      # looping through each trawl
      y[i] ~ dbin(p[i],N[i]) #Match distribution to data, in this case binomial
      
      logit(p[i])<-beta.geo[geoarea[i]]
                    +beta.species[species[i]]
                    +beta.year[species[i],geoarea[i]] * year[i]
                    +eps[i] #Linear predictor of mean, two intercepts plus slope by year+overdispersion random term
      eps[i] ~ dnorm(0,tau) #Distribution for overdispersion term
      
      # RESIDUALS 
      presid[i] <- (y[i]-(N[i]*p[i]))/(sqrt(p[i]*(1-p[i])*N[i])+0.0001)          #Pearson are the difference between data and mean/variance (for binom that's sqrt(p*(1-p)), but find for each)
      p2[i] <- presid[i]^2
      #Generate new data and calculate its residuals
        new[i] ~ dbin(p[i],N[i])                                        #Use same distribution as above
        presid.new[i] <- (new[i]-(p[i]*N[i]))/(sqrt(p[i]*(1-p[i])*N[i])+0.0001)  #Again, use mean and variance for that distribution
        p2.new[i] <- presid.new[i]^2
      
    } #end of ntrawls
      
    # CHECK MODEL FIT - BAYESIAN P-VALUE
    psum <- sum(p2[])  
    psum.new <- sum(p2.new[]) 
    Bayes.P <- step(psum.new/psum - 1) 
    
    } #end of the model
",fill=T)
sink()

win.data<-list(ntrawls=nrow(empty2),
               y=empty2$freq,
               N=empty2$nDiets,
               nspecies=length(keyspecies_com),
               species=as.integer(empty2$comname),
               year=as.integer(empty2$year-min(empty2$year,na.rm=T)),
               geoarea=as.integer(factor(empty2$geoarea,levels=c("MAB","SNE","GB","GoM","ScS"))),
               ngeo=length(unique(empty2$geoarea)))

inits<-function()list()

params<-c("beta.geo","beta.species","beta.year","psum","psum.new","Bayes.P")

nc<-3; nt<-1; ni<-5000; nb<-1000

outSpecies_LM <- jags(win.data,inits,params,"speciesEmpty_LM.txt",
                      n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb,
                      working.directory=getwd())

print(outSpecies_LM,digits=3)
traceplot(outSpecies_LM)

oo <- outSpecies_LM$BUGSoutput
#save(oo,file="EmptyBugs.RData")


# JAGS for Levin Breadth --------------------------------------------------
  #Have to first squeeze the values to be (0,1)
trawlDiets_Levin<-trawlDiets_Levin%>%
  mutate(squeezeSimpson=(simpson*(n_distinct(trawlDiets_Levin$id)-1)+0.5)/n_distinct(trawlDiets_Levin$id),
         squeezeLevin=(levinStd*(n_distinct(trawlDiets_Levin$id)-1)+0.5)/n_distinct(trawlDiets_Levin$id))

sink("dietsLevin.txt")
cat("
    model {
    # PRIORS
    for (i in 1:nspecies) {
      beta.species[i] ~ dnorm(int.s.mean, int.s.pres)
      for (j in 1:ngeo) {
         beta.year[i,j] ~ dnorm(slope.mean,slope.pres)
         phi[i,j] ~ dgamma(0.1,0.1)
      }
    }
    for (i in 1:ngeo) {
      beta.geo[i] ~ dnorm(int.r.mean, int.s.pres)
    }
    
    # Hyperpriors
      #Diffuse priors for geographic region intercepts
    int.r.mean ~ dnorm(0,0.001)
    int.r.var ~ dunif(0,10)
    int.r.pres <- 1/int.r.var^2
    
      #Diffuse priors for species intercepts
    int.s.mean ~ dnorm(0,0.001)
    int.s.var ~ dunif(0,10)
    int.s.pres <- 1/int.s.var^2
    
      #Diffuse priors for slope parameters
    slope.mean ~ dnorm(0,0.001)
    slope.var ~ dunif(0,10)
    slope.pres <- 1/slope.var^2
    
      #Diffuse prior for overdispersion
    tau <- 1/sigma^2
    sigma ~ dunif(0,100)
    
    # LIKELIHOOD
    for (i in 1:ntrawls) {
      
      # looping through each trawl
      y[i] ~ dbeta(alpha[i], beta[i])
      alpha[i] <- mu[i] * phi[species[i],geoarea[i]]
      beta[i]  <- (1-mu[i]) * phi[species[i],geoarea[i]]
      
      logit(mu[i])<-beta.geo[geoarea[i]]
                    +beta.species[species[i]]
                    +beta.year[species[i],geoarea[i]] * year[i]
                    +eps[i] #Linear predictor of mean, two intercepts and the mean+overdispersion random term
      eps[i] ~ dnorm(0,tau) #Overdispersion term distribution is mean 0 and variance tau
      
      # RESIDUALS 
      presid[i] <- (y[i]-mu[i])/sqrt((alpha[i]*beta[i])/((alpha[i]+beta[i]+1)*(alpha[i]+beta[i])^2))        #Pearson are the difference between data and mean/variance (for binom that's sqrt(p*(1-p)), but find for each)
      p2[i] <- presid[i]^2
      #Generate new data and calculate its residuals
        new[i] ~ dbeta(alpha[i], beta[i])  #Use same distribution as above
        presid.new[i] <- (new[i]-mu[i])/sqrt((alpha[i]*beta[i])/((alpha[i]+beta[i]+1)*(alpha[i]+beta[i])^2))  #Again, use mean and variance for that distribution
        p2.new[i] <- presid.new[i]^2
      
    } #end of ntrawls
      
    # CHECK MODEL FIT - BAYESIAN P-VALUE
    psum <- sum(p2[])  
    psum.new <- sum(p2.new[]) 
    Bayes.P <- step(psum.new/psum - 1) 
    
    } #end of the model
",fill=T)
sink()

win.data<-list(ntrawls=nrow(trawlDiets_Levin),
               y=trawlDiets_Levin$squeezeSimpson,
               nspecies=length(keyspecies_com),
               species=as.integer(trawlDiets_Levin$comname),
               year=as.integer(trawlDiets_Levin$year-min(trawlDiets_Levin$year,na.rm=T)),
               geoarea=as.integer(factor(trawlDiets_Levin$geoarea,levels=c("MAB","SNE","GB","GoM","ScS"))),
               ngeo=length(unique(trawlDiets_Levin$geoarea)))

inits<-function()list()

params<-c("phi","beta.geo","beta.species","beta.year","psum","psum.new","Bayes.P")

nc<-3; nt<-1; ni<-5000; nb<-1000

outLevin<-jags(win.data,inits,params,"dietsLevin.txt",
               n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb,
               working.directory=getwd())
print(outLevin,digits=3)
traceplot(outLevin)
ooo <- outLevin$BUGSoutput
#save(ooo,file="levinBUGS.RData")


# JAGS for Consumption ----------------------------------------------------
sink("dietsConsumption.txt")
cat("
    model {
    # PRIORS
    for (i in 1:nspecies) {
      beta.species[i] ~ dnorm(int.s.mean, int.s.pres)
      for (j in 1:ngeo) {
         beta.year[i,j] ~ dnorm(slope.mean,slope.pres)
      }
    }
    for (i in 1:ngeo) {
      beta.geo[i] ~ dnorm(int.r.mean, int.s.pres)
    }
    
    # Hyperpriors
      #For the regional intercept, diffuse priors
    int.r.mean ~ dnorm(0,0.001) 
    int.r.var ~ dunif(0,10)
    int.r.pres <- 1/int.r.var^2
    
      #For the species intercept, diffuse priors
    int.s.mean ~ dnorm(0,0.001)
    int.s.var ~ dunif(0,10)
    int.s.pres <- 1/int.s.var^2
    
      #For the slope, diffuse priors
    slope.mean ~ dnorm(0,0.001)
    slope.var ~ dunif(0,10)
    slope.pres <- 1/slope.var^2
    
      #For the variance in data
    tau <- 1/sigma^2
    sigma ~ dunif(0,100)
    
    # LIKELIHOOD
    for (i in 1:ntrawls) {
      
      # looping through each trawl
      y[i] ~ dlnorm(mu[i],tau)   #Match distribution to data, log-normal here
      
      mu[i]<-beta.geo[geoarea[i]]
              +beta.species[species[i]]
              +beta.year[species[i],geoarea[i]] * year[i] #Linear predictor of two intercepts and slope
      
      # RESIDUALS 
      presid[i] <- (y[i]-mu[i])/sigma       #Pearson are the difference between data and mean/variance (for binom that's sqrt(p*(1-p)), but find for each)
      p2[i] <- presid[i]^2
      #Generate new data and calculate its residuals
        new[i] ~ dlnorm(mu[i],tau)                                      #Use same distribution as above
        presid.new[i] <- (new[i]-mu[i])/sigma  #Again, use mean and variance for that distribution
        p2.new[i] <- presid.new[i]^2
      
    } #end of ntrawls
      
    # CHECK MODEL FIT - BAYESIAN P-VALUE
    psum <- sum(p2[])  
    psum.new <- sum(p2.new[]) 
    Bayes.P <- step(psum.new/psum - 1) 
    
    } #end of the model
",fill=T)
sink()

win.data<-list(ntrawls=nrow(trawlDiets_cind),
               y=trawlDiets_cind$meanC,
               nspecies=length(keyspecies_com),
               species=as.integer(trawlDiets_cind$comname),
               year=as.integer(trawlDiets_cind$year-min(trawlDiets_cind$year,na.rm=T)),
               geoarea=as.integer(trawlDiets_cind$geoarea),
               ngeo=length(unique(trawlDiets_cind$geoarea)))

inits<-function()list()

params<-c("beta.geo","beta.year","beta.species","tau","psum","psum.new","Bayes.P")

nc<-3; nt<-1; ni<-5000; nb<-1000

outConsump<-jags(win.data,inits,params,"dietsConsumption.txt",
                 n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb,
                 working.directory=getwd())
print(outConsump,digits=3)
traceplot(outConsump)


o <- outConsump$BUGSoutput
#save(o,file="ConsumpBUGS.RData")

# Summarizing Output ------------------------------------------------------

#Percent Empty Output
empty.geo.ints <- oo$mean$beta.geo
empty.species.ints <- oo$mean$beta.species
empty.slopes <- as.numeric(oo$mean$beta.year)
predEmptyData <- data.frame(geo.ints=rep(empty.geo.ints,each=9), 
                              species.ints=empty.species.ints,
                              slopes=empty.slopes,
                              comname=sort(str_to_sentence(keyspecies_com)),
                              geoarea=factor(c(rep("MAB",9),rep("SNE",9),rep("GB",9),rep("GoM",9),rep("ScS",9)),
                                             levels=c("MAB","SNE","GB","GoM","ScS")))
year_values <- seq(min(empty2$year)-min(empty2$year),
                   max(empty2$year)-min(empty2$year),by=1)

predEmptyData <- cbind(predEmptyData,years=sort(rep(year_values,nrow(predEmptyData))))%>%
  mutate(LOGITpreds=geo.ints+species.ints+slopes*years,
         preds=plogis(LOGITpreds))

ggplot()+
  geom_point(data=empty,aes(year,pik,color=comname,shape=geoarea),alpha=0.5,show.legend=F)+
  geom_line(data=filter(predEmptyData,!(comname=="Summer flounder"&geoarea=="ScS")),aes(years+min(empty2$year),preds),size=2)+
  scale_x_continuous(limits=c(min(trawlDiets$year),max(trawlDiets$year)),name="Year")+
  scale_y_continuous(name="Proportion of Stomachs Empty")+
  theme(axis.text.x=element_text(size=12,angle=30,hjust=1,vjust=1))+
  facet_grid(factor(geoarea,levels=c("ScS","GoM","GB","SNE","MAB"))~comname)

predEmptyData%>%group_by(geoarea)%>%filter(years==0)%>%summarise(mean=mean(preds),range=range(preds))

#Relative Consumption Output
consump.geo.ints <- o$mean$beta.geo
consump.species.ints <- o$mean$beta.species
consump.slopes <- as.numeric(o$mean$beta.year)
predConsumpData <- data.frame(geo.ints=rep(consump.geo.ints,each=9), 
                              species.ints=consump.species.ints,
                              slopes=consump.slopes,
                              comname=sort(str_to_sentence(keyspecies_com)),
                              geoarea=factor(c(rep("MAB",9),rep("SNE",9),rep("GB",9),rep("GoM",9),rep("ScS",9)),
                                             levels=c("MAB","SNE","GB","GoM","ScS")))
year_values <- seq(min(trawlDiets_cind$year)-min(trawlDiets_cind$year),
                   max(trawlDiets_cind$year)-min(trawlDiets_cind$year),by=1)

predConsumpData <- cbind(predConsumpData,years=sort(rep(year_values,nrow(predConsumpData))))%>%
  mutate(LOGpreds=geo.ints+species.ints+slopes*years,
         preds=exp(LOGpreds))

ggplot()+
  geom_point(data=trawlDiets_cind,aes(year,meanC,color=comname,shape=geoarea),alpha=0.5,show.legend=F)+
  geom_line(data=predConsumpData,aes(years+min(trawlDiets_cind$year),preds),size=2)+
  scale_x_continuous(limits=c(min(trawlDiets$year),max(trawlDiets$year)),name="Year")+
  scale_y_log10(name="Relative Consumption (g/g)")+
  theme(axis.text.x=element_text(size=12,angle=30,hjust=1,vjust=1))+
  facet_grid(factor(geoarea,levels=c("ScS","GoM","GB","SNE","MAB"))~comname)

predConsumpData%>%filter(years==0 & comname%in%c("Haddock","Little skate","Spiny dogfish","Summer flounder","Yellowtail flounder"))%>%summarise(mean(preds))

#Breadth Output
breadth.geo.ints <- ooo$mean$beta.geo
breadth.species.ints <- ooo$mean$beta.species
breadth.slopes <- as.numeric(ooo$mean$beta.year)
predBreadthData <- data.frame(geo.ints=rep(breadth.geo.ints,each=9), 
                            species.ints=breadth.species.ints,
                            slopes=breadth.slopes,
                            comname=sort(str_to_sentence(keyspecies_com)),
                            geoarea=factor(c(rep("MAB",9),rep("SNE",9),rep("GB",9),rep("GoM",9),rep("ScS",9)),
                                           levels=c("MAB","SNE","GB","GoM","ScS")))
year_values <- seq(min(trawlDiets_Levin$year)-min(trawlDiets_Levin$year),
                   max(trawlDiets_Levin$year)-min(trawlDiets_Levin$year),by=1)

predBreadthData <- cbind(predBreadthData,years=sort(rep(year_values,nrow(predBreadthData))))%>%
  mutate(LOGITpreds=geo.ints+species.ints+slopes*years,
         preds=plogis(LOGITpreds))

ggplot()+
  geom_point(data=trawlDiets_Levin,aes(year,squeezeSimpson,color=comname,shape=geoarea),alpha=0.5,show.legend=F)+
  geom_line(data=filter(predBreadthData,!(comname=="Summer flounder"&geoarea=="ScS")),
            aes(years+min(trawlDiets_Levin$year),preds),size=2)+
  scale_x_continuous(limits=c(min(trawlDiets$year),max(trawlDiets$year)),name="Year")+
  scale_y_continuous(name="Simpson's Diversity Index")+
  theme(axis.text.x=element_text(size=12,angle=30,hjust=1,vjust=1))+
  facet_grid(factor(geoarea,levels=c("ScS","GoM","GB","SNE","MAB"))~comname)

predBreadthData%>%filter(years==0)%>%summarise(mean(preds))
predBreadthData%>%group_by(comname,geoarea)%>%summarise(decrease=max(preds)-min(preds))%>%ungroup()%>%summarise(mean(decrease))
0.142/0.2817

predBreadthData%>%filter(years==0)%>%group_by(geoarea)%>%summarise(mean(preds))
predBreadthData%>%filter(years==0)%>%group_by(comname)%>%summarise(mean(preds))


#differences in residuals
(o$mean$psum.new-o$mean$psum)/o$mean$psum
(oo$mean$psum.new-oo$mean$psum)/oo$mean$psum
(ooo$mean$psum.new-ooo$mean$psum)/ooo$mean$psum
