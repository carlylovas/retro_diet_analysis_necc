
rm(list=ls())

# Top ---------------------------------------------------------------------


setwd("/Volumes/My Passport for Mac/NECC/Diet Data")

library(cowplot)
library(tidyr)
library(tibble)
library(dplyr)
library(forcats)
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
library(dendextend)
library(car)
library(ggpubr)
library(MuMIn)
library(lme4)
library(ggspatial)
library(sf)
library(colorspace)
library(rnaturalearth) 
library(rnaturalearthdata)
library(rgeos)
library(caret)
library(Hmisc)
library(lsmeans)
library(ggrepel)
library(stringi)
library(MCMCglmm)


#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}




# Data Loading ----------------------------------
load("NF.prey19.RData")
allprey<-NF.prey19
rm(NF.prey19)

#ALL LEVELS OF YEAR AND SEASON
ALLyearseasons<-factor(levels=c(paste(sort(rep(seq(1973,2019),4)),rep(c("Winter","Spring","Summer","Fall"),length(seq(1973,2019))))))

#The names of my species, just to have them nice and handy
species_sci<-unique(str_to_sentence(allprey$pdscinam))
species_com<-unique(str_to_title(allprey$pdcomnam))
#helping me decode them
species_key<-unique(dplyr::select(allprey,svspp,pdscinam,pdcomnam))%>%
  mutate(svspp = str_pad(svspp,3,side="left"))

#These are all the species in the trawls
spec<-read.csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv") %>%
  dplyr::select(LOGGED_SPECIES_NAME,SVSPP)%>%distinct()%>%
  filter(SVSPP>="001")%>%
  mutate(LOGGED_SPECIES_NAME=str_to_sentence(LOGGED_SPECIES_NAME))
spec<-spec[!duplicated(spec$SVSPP),]

#Size classes used by Garrison and Link
sizeClasses<-read.csv("GarrisonLink_predSizeClasses.csv")
sizeClasses <- sizeClasses %>%
  left_join(spec,by=c("species_comnam"="LOGGED_SPECIES_NAME"))%>%
  rename(svspp=SVSPP)



# Important df manipulations ----------------------------------------------

allprey<-allprey%>%
  mutate(cruise6=str_pad(cruise6,width=6,side="left",pad="0"),
         svspp=str_pad(svspp,width=3,side="left",pad="0"),
         station=str_pad(station,width=3,side="left",pad="0"),
         stratum=str_pad(stratum,width=3,side="left",pad="0"),
         id=paste0(cruise6,
                   station,
                   stratum),
         pdid=str_pad(pdid,width=6,side="left",pad="0"),
         dietID=paste0(svspp,
                       pdsex,
                       pdid,
                       str_pad(pdlen,width=3,pad="0",side="left"),
                       id),
         season=str_to_sentence(season))%>%
  group_by(id,svspp)%>%
  mutate(nDiets=n_distinct(dietID))%>%
  left_join(sizeClasses,by=c("svspp"))%>%
  group_by(svspp)%>%
  mutate(sizecat2=case_when(between(pdlen,min(small_min),min(small_max))~"S",
                            between(pdlen,min(medium_min),min(medium_max))~"M",
                            between(pdlen,min(large_min),min(large_max))~"L",
                            between(pdlen,min(xlarge_min),min(xlarge_max))~"XL",
                            TRUE ~ "S"),
         year=ifelse(is.na(year),substr(cruise6,1,4),year),year=as.numeric(year))%>%
  select(-c(species_scinam:xlarge_max))%>%
  #Instances where there is a pynam but no gensci, analsci, OR collsci
  mutate(gensci=ifelse(pynam %in% c("PRIONOTUS ALATUS", #spiny searobin
                                    "STERNOPTYCHIDAE", #hatchetfishes
                                    "EPIGONUS PANDIONIS", #bieye
                                    "MANTA BIROSTRIS", #giant manta ray
                                    "DACTYLOPTERUS VOLITANS", #Flying gurnard
                                    "SCYLIORHINUS RETIFER", #Chain catshark
                                    "SELENE SETAPINNIS", #Atlantic moonfish
                                    "OGCOCEPHALIDAE", #Batfishes
                                    "SYNAGROPS BELLUS", #blackmouth bass
                                    "LEUCORAJA GARMANI", #rosette skate
                                    "PARASUDIS TRUCULENTA", #longnose greeneye
                                    "MONOLENE SESSILICAUDA"), #deepwater flounder
                       "FISH",as.character(gensci)),
         analsci=case_when(pynam %notin% c("PRIONOTUS ALATUS", #spiny searobin
                                           "STERNOPTYCHIDAE", #hatchetfishes
                                           "EPIGONUS PANDIONIS", #bieye
                                           "MANTA BIROSTRIS", #giant manta ray
                                           "DACTYLOPTERUS VOLITANS", #Flying gurnard
                                           "SCYLIORHINUS RETIFER", #Chain catshark
                                           "SELENE SETAPINNIS", #Atlantic moonfish
                                           "OGCOCEPHALIDAE", #Batfishes
                                           "SYNAGROPS BELLUS", #blackmouth bass
                                           "LEUCORAJA GARMANI", #rosette skate
                                           "PARASUDIS TRUCULENTA", #longnose greeneye
                                           "MONOLENE SESSILICAUDA")~as.character(analsci),
                           pynam=="PRIONOTUS ALATUS"~"TRIGLIDAE",
                           pynam=="STERNOPTYCHIDAE"~"STERNOPTYCHIDAE",
                           pynam=="EPIGONUS PANDIONIS"~"EPIGONIDAE",
                           pynam=="MANTA BIROSTRIS"~"MOBULIDAE",
                           pynam=="DACTYLOPTERUS VOLITANS"~"DACTYLOPTERIDAE",
                           pynam=="SCYLIORHINUS RETIFER"~"SCYLIORHINIDAE",
                           pynam=="SELENE SETAPINNIS"~"CARANGIDAE",
                           pynam=="OGCOCEPHALIDAE"~"OGCOCEPHALIDAE",
                           pynam=="SYNAGROPS BELLUS"~"SYNAGROPIDAE",
                           pynam=="LEUCORAJA GARMANI"~"RAJIFORMES",
                           pynam=="PARASUDIS TRUCULENTA"~"CHLOROPHTHALMIDAE",
                           pynam=="MONOLENE SESSILICAUDA"~"BOTHIDAE"),
         collsci=ifelse(pynam %in% c("PRIONOTUS ALATUS", #spiny searobin
                                     "STERNOPTYCHIDAE", #hatchetfishes
                                     "EPIGONUS PANDIONIS", #bieye
                                     "MANTA BIROSTRIS", #giant manta ray
                                     "DACTYLOPTERUS VOLITANS", #Flying gurnard
                                     "SCYLIORHINUS RETIFER", #Chain catshark
                                     "SELENE SETAPINNIS", #Atlantic moonfish
                                     "OGCOCEPHALIDAE", #Batfishes
                                     "SYNAGROPS BELLUS", #blackmouth bass
                                     "LEUCORAJA GARMANI", #rosette skate
                                     "PARASUDIS TRUCULENTA", #longnose greeneye
                                     "MONOLENE SESSILICAUDA"),
                        as.character(pynam),as.character(collsci)),
         #pynam2.0=gsub("DECAPODA CRAB ?[A-z]+","DECAPODA CRAB",pynam),
         #pynam2.0=gsub("DECAPODA SHRIMP ?[A-z]+","DECAPODA SHRIMP",pynam2.0),
         pynam2.0=gsub("LOLIGO [A-z]+","LOLIGO SP",pynam),
         pynam2.0=gsub("ILLEX [A-z]+","ILLEX SP",pynam2.0),
         pynam2.0=gsub("AMMODYTES [A-z]+","AMMODYTES SP",pynam2.0),
         pynam2.0=ifelse(gensci=="FISH",gsub(" EGGS$","",pynam2.0),pynam2.0),
         pynam2.0=ifelse(gensci=="FISH",gsub(" LARVAE$","",pynam2.0),pynam2.0))


# Creating prey categories that mirror GL -----------------------------------------

#The categories from Garrison and Link, 2000
gl_preycats<-read.csv(here("data", "raw","GarrisonLink_preyCats.csv")) %>%
  mutate(matchingCats=gsub("p\\.","",Scientific.name),
         matchingCats=gsub("crabs","crab",matchingCats),
         matchingCats=gsub("Gammaridae","Gammaridea",matchingCats),
         matchingCats=gsub("Cnidarians","Cnidaria",matchingCats))

allprey<-allprey%>%
  mutate(INgen=ifelse(str_to_sentence(gensci) %in% gl_preycats$matchingCats,1,0),
         INanal=ifelse(str_to_sentence(analsci) %in% gl_preycats$matchingCats,1,0),
         INcoll=ifelse(str_to_sentence(collsci) %in% gl_preycats$matchingCats,1,0),
         INpy=ifelse(str_to_sentence(pynam2.0) %in% gl_preycats$matchingCats,1,0))
allprey<-allprey%>%
  ungroup()%>%
  mutate(INnum=rowSums(allprey[,c("INgen","INanal","INcoll","INpy")]),
         gl_prey=case_when(INnum==4~str_to_sentence(pynam), #STEP 1 
                           INnum==3&INpy==1~str_to_sentence(pynam2.0),
                           INnum==3&INpy==0~str_to_sentence(collsci),
                           INnum==2&INpy==1~str_to_sentence(pynam2.0),
                           INnum==2&INgen==1~str_to_sentence(analsci),
                           INnum==2&INgen==0&INpy==0~str_to_sentence(collsci),
                           INnum==1&INgen==1~str_to_sentence(gensci),
                           INnum==1&INanal==1~str_to_sentence(analsci),
                           INnum==1&INcoll==1~str_to_sentence(collsci),
                           INnum==1&INpy==1~str_to_sentence(pynam2.0),
                           INnum==0&pynam=="EMPTY"~"Empty"),
         gl_prey=case_when(pynam %in% c("UROPHYCIS CHUSS","UROPHYCIS TENUIS",
                                        "UROPHYCIS REGIA")~"Other hakes", #STEP 2
                           (analsci %in% c("GADIDAE","BREGMACEROTIDAE",
                                           "EUCLICHTHYIDAE","LOTIDAE","MACROURIDAE",
                                           "MELANONIDAE","MERLUCCIIDAE","MORIDAE",
                                           "MURAENOLEPIDIDAE","PHYCIDAE")
                            | pynam2.0 %in% c("GADIDAE","BREGMACEROTIDAE",
                                              "EUCLICHTHYIDAE","LOTIDAE","MACROURIDAE",
                                              "MELANONIDAE","MERLUCCIIDAE","MORIDAE",
                                              "MURAENOLEPIDIDAE","PHYCIDAE"))
                           & (is.na(gl_prey) 
                              | gl_prey %in% c("Fish larvae",
                                               "Fish eggs"))~"Gadiformes", #STEP 3a
                           (analsci %in% c("PLEURONECTIDAE","PSETTODIDAE",
                                           "CITHARIDAE","SCOPHTHALMIDAE",
                                           "PARALICHTHYIDAE","BOTHIDAE",
                                           "PARALICHTHODIDAE","POECILOPSETTIDAE",
                                           "RHOMBOSOLEIDAE","ACHIROPSETTIDAE",
                                           "SAMARIDAE","ACHIRIDAE",
                                           "SOLEIDAE","CYNOGLOSSIDAE")
                            | pynam2.0 %in% c("PLEURONECTIDAE","PSETTODIDAE",
                                              "CITHARIDAE","SCOPHTHALMIDAE",
                                              "PARALICHTHYIDAE","BOTHIDAE",
                                              "PARALICHTHODIDAE","POECILOPSETTIDAE",
                                              "RHOMBOSOLEIDAE","ACHIROPSETTIDAE",
                                              "SAMARIDAE","ACHIRIDAE",
                                              "SOLEIDAE","CYNOGLOSSIDAE",
                                              "ETROPUS SP","GLYPTOCEPHALUS CYNOGLOSSUS",
                                              "SCOPHTHALMUS AQUOSUS"))
                           & (is.na(gl_prey)
                              | gl_prey %in% c("Fish larvae",
                                               "Fish eggs"))~"Pleuronectiformes", #STEP 3b
                           grepl("MYOXOCEPHALUS",pynam2.0)~"Cottidae", #STEP 4
                           pynam %in% c("FISH SCALES",
                                        "FISH OTOLITHS",
                                        "FISH")~"Unidentified fish", #STEP 5
                           gensci=="FISH" 
                           & is.na(gl_prey)~"Other fish", #STEP 6
                           #gensci=="FISH"
                           #   & grepl("EGGS",pynam)~"Fish eggs", #STEP 6a
                           #gensci=="FISH"
                           #   & grepl("LARVAE",pynam)~"Fish larvae", #STEP 6b
                           gensci %in% c("CHAETOGNATHA")
                           | analsci %in% c("COPEPODA")
                           | collsci=="OSTRACODA"
                           | pynam=="PLANKTON" 
                           | grepl("megalop",pynam2.0,ignore.case=T)
                           | grepl("zoea",pynam2.0,ignore.case=T)
                           | (gensci=="ARTHROPODA" 
                              & grepl("larvae",pynam2.0,ignore.case=T))~"Zooplankton", #STEP 7
                           collsci %in% c("OLIGOCHAETA",
                                          "HIRUDENEA")~"Worms", #STEP 8
                           analsci %in% c("CEPHALOCHORDATA")~"Other", #STEP 9
                           collsci %in% c("CUMACEA",
                                          "STOMATOPODA")~"Crustacean shrimp", #STEP 10
                           collsci %in% c("PENAEIDAE",
                                          "HOMARUS AMERICANUS",
                                          "SCYLLARIDAE")
                           | grepl("PALINURA",pynam2.0)~"Decapoda shrimp",#STEP 11a
                           collsci %in% c("CALLINECTES SAPIDUS")~"Decapoda crab", #STEP 11b
                           analsci %in% c("CIRRIPEDIA") 
                           | collsci %in% c("DECAPODA","DECAPODA EGGS",
                                            "DECAPODA LARVAE") 
                           | pynam %in% c("DECAPODA","DECAPODA EGGS",
                                          "DECAPODA LARVAE")~"Crustacea", #STEP 12
                           analsci=="EUPHAUSIACEA"~"Euphausiidae", #STEP 13
                           collsci %in% c("APHRODITIDAE")~"Polychaeta", #STEP 14
                           gensci %in% c("UROCHORDATA","BRACHIOPODA",
                                         "BRYOZOA",
                                         "PORIFERA")
                           | collsci %in% c("ARTHROPODA","INSECTA",
                                            "HEMICHORDATA","LIMULUS POLYPHEMUS",
                                            "PYCNOGONIDA","HALACARIDAE")
                           | pynam=="INVERTEBRATA"~"Other invertebrates", #STEP 15
                           !is.na(gl_prey)~gl_prey), #Keep gl_prey from above
         gl_prey=factor(gl_prey,levels=c(gl_preycats$matchingCats,"Empty")))%>%
  left_join(gl_preycats,by=c("gl_prey"="matchingCats"))%>%
  mutate(GL_pyscinam=factor(ifelse(gl_prey=="Empty","Empty",Scientific.name),
                            levels=c(gl_preycats$Scientific.name,"Empty")),
         GL_pycomnam=factor(ifelse(gl_prey=="Empty","Empty",Common.name),
                            levels=c(gl_preycats$Common.name,"Empty")))%>%
  dplyr::select(-c(Scientific.name,Common.name))


#check<-filter(allprey,is.na(GL_pycomnam))[,c("gensci","analsci","collsci","pynam","INnum","INpy")]
#sort(table(allprey$GL_pyscinam))
#printOut<-allprey%>%
#  group_by(gensci,analsci,collsci,pynam,GL_pyscinam)%>%
#  summarise(N=n())%>%
#  arrange(GL_pyscinam)
#write.csv(printOut,"allprey_glpreycats.v15.csv",row.names = F)

uniqueDiets<-allprey%>%
  mutate(Empty=ifelse(pynam=="EMPTY","Y","N"))%>%
  dplyr::select(cruise6,station,svspp,pdsex,pdid,pdcomnam,
                pdscinam,dietID,pdlen,pdwgt,sizecat2,Empty,pdgutw,pdgutv,
                declat,declon,month,day,year,season,geoarea,id,nDiets)%>%
  distinct()


nrow(filter(uniqueDiets, Empty=="Y"))/nrow(uniqueDiets)


#Making an order list of the GL cats, so I can set factor levels to this 
#First the fishes, sp-gen-fam-other-unid
gl_fish<-c("Clupea harengus","Clupeidae","Peprilus triacanthus","Ammodytes spp.",
           "Lepophidium profundorum","Macrozoarces americanus","Merluccius bilinearis",
           "Other hakes","Gadiformes","Cottidae","Pleuronectiformes",
           "Rajiformes","Scombridae","Engraulidae",
           "Other fish","Unidentified fish","Fish larvae","Fish eggs")
gl_inverts<-c("Loligo spp.","Illex spp.","Cephalopoda",
              "Crangonidae","Euphausiidae","Pandalidae","Mysidacea",
              "Crustacean shrimp","Crustacea","Zooplankton",
              "Cancridae","Decapoda crabs","Decapoda shrimp","Paguroidea",
              "Hyperiidae","Gammaridae","Amphipoda","Isopoda",
              "Bivalvia","Gastropoda","Mollusca","Echinodermata","Ophiuroidea",
              "Holothuroidea","Hydrozoa","Anthozoa","Cnidarians","Ctenophora",
              "Polychaeta","Worms","Other invertebrates")
gl_other<-c("Animal remains","Other","Miscellaneous","Empty",NA)
gl_levels<-c(gl_fish,gl_inverts,gl_other)

gl_comlevels<-c("Atlantic herring","Herrings","Butterfish","Sand lance","Fawn cusk-eel","Ocean pout",
                "Silver hake","Red, white, and spotted hake","Gadid fish","Sculpins","Flatfish","Rays, skates",
                "Mackerels","Anchovies","Other fish","Unidentified fish remains","Fish larvae","Fish eggs",
                "Loligo squid","Illex squid","Unclassified cephalopods","Crangonid shrimps","Euphausiid shrimp, krill",
                "Pandalid shrimps","Mysid shrimps","Unclassified shrimp","Unclassified crustaceans","Zooplankton",
                "Cancer crabs","Decapod crabs","Decapod shrimp","Hermit crabs","Hyperiid amphipods","Gammarid amphipods",
                "Unclassified amphipods","Isopods","Bivalves","Snails","Unclassified molluscs","Unclassified echinoderms",
                "Brittle stars","Sea cucumbers","Hydroids","Anemones","Jellies and hydroids","Comb-jellies",
                "Polychaete worms","Assorted worms","Other invertebrates","Unidentified animal remains",
                "Unidentified organic material","Inorganic material","Empty",NA)

allprey<-allprey%>%
  mutate(cruise6=as.character(cruise6),
         station=as.character(station),
         stratum=as.character(stratum))


# Adding in the trawl data (from GMRI) ------------------------------------

#reading in prey trawls df from "connecting_food+trawls.R" where prey19 is merged with GMRI's clean data for trawls
trawls<-read_csv(here("data","processed","NMFS_survdat_gmri_tidy.csv"),
                 col_types=c("cccccccnnncnnTnnnnnnnnnnccccn"))
TRAWLS<-trawls%>%
  filter(svspp %in% c(species_key$svspp,"401"))%>%
  dplyr::select(svspp,comname,catchsex,cruise6,stratum,station,tow,id,decdeg_beglon,decdeg_beglat,est_year,season,biomass_kg,abundance)%>%distinct()
colnames(TRAWLS)<-toupper(colnames(TRAWLS))
colnames(TRAWLS)<-gsub("DECDEG_BEG","",colnames(TRAWLS))
colnames(TRAWLS)<-gsub("EST_","",colnames(TRAWLS))
colnames(TRAWLS)<-gsub("_KG","",colnames(TRAWLS))

abundance<-trawls%>%
  ungroup()%>%
  dplyr::select(year=est_year,season,id,pdcomnam=comname,catchsex,abundance,biomass_kg)%>%distinct()%>%
  group_by(year,season,id,pdcomnam)%>%
  summarise(abundance=sum(abundance),
            biomass_kg=sum(biomass_kg))

trawls<-trawls%>%
  left_join(sizeClasses,by=c("svspp"))%>%
  filter(svspp %in% allprey$svspp)%>%
  group_by(svspp)%>%
  mutate(sizecat2=case_when(dplyr::between(length_cm,min(small_min,na.rm=T), min(small_max,na.rm=T))~"S",
                            dplyr::between(length_cm,min(medium_min,na.rm=T),min(medium_max,na.rm=T))~"M",
                            dplyr::between(length_cm,min(large_min,na.rm=T), min(large_max,na.rm=T))~"L",
                            dplyr::between(length_cm,min(xlarge_min,na.rm=T),min(xlarge_max,na.rm=T))~"XL",
                            TRUE~"S"))%>%
  group_by(svspp,id,sizecat2)%>%
  mutate(sizecat_abundance=sum(numlen_adj))%>%
  dplyr::select(-c(abundance,biomass_kg,catchsex,length_cm:numlen_adj,n_len_class,small_min:xlarge_max))%>%
  distinct()%>%
  left_join(abundance,by=c("id"="id","est_year"="year","season"="season","comname"="pdcomnam"))
#Want to use just one row for each trawl, so combined the sexes and removed the length classes
preyTrawls_nSize<-left_join(allprey,trawls%>%ungroup()%>%dplyr::select(-c(sizecat2,sizecat_abundance))%>%distinct())%>% #Remove the sizecat because there are some diets from fish in sizecats that weren't caught
  mutate(abundance=ifelse(is.na(abundance),nDiets,abundance),
         abundance=ifelse(nDiets>abundance,nDiets,abundance)) #Sometimes there are more diets collected than fish caught (like kinda a lot of the time)

preyTrawls<-left_join(allprey,trawls%>%distinct())%>% #Remove the sizecat because there are some diets from fish in sizecats that weren't caught
  mutate(abundance=ifelse(is.na(abundance),nDiets,abundance),
         abundance=ifelse(nDiets>abundance,nDiets,abundance)) #Sometimes there are more diets collected than fish caught (like kinda a lot of the time)


uniqueDietTrawls<-left_join(uniqueDiets,trawls%>%distinct())

rm(trawls)





#Need to have whole diet totals, and whole trawl measures
trawlDiets<-preyTrawls%>%
  group_by(id,pdscinam,sizecat2)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv),
         totalN=n_distinct(dietID))%>%
  group_by(id,pdscinam,sizecat2,GL_pyscinam)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(totalN==0,0,n_distinct(dietID)/totalN),
         pdcomnam=str_to_title(pdcomnam),
         pdscinam=str_to_sentence(pdscinam))%>%
  dplyr::select(year,season,geoarea,id,sizecat_abundance,abundance,
                #geoarea,season,    #Any covariates to include?
                pdscinam,pdcomnam,sizecat2,GL_pyscinam,
                totalwt,totalv,totalN,qikw,qikv,pik)%>% 
  distinct() 

#if grouping, add those groups in here both at select and group_by
nDiets<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case the whole pred sizecat)
  mutate(period=ifelse(as.numeric(substr(id,1,4))<=1997,"GL","New"))%>%
  ungroup()%>%
  dplyr::select(period,pdscinam,sizecat2,id,totalN)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(period,pdscinam,sizecat2)%>%summarise(nDiets=sum(totalN))#Calculate the number of diets
sumAbun<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  mutate(period=ifelse(as.numeric(substr(id,1,4))<=1997,"GL","New"))%>%
  ungroup()%>%
  dplyr::select(period,pdscinam,sizecat2,id,sizecat_abundance)%>%distinct()%>%
  group_by(period,pdscinam,sizecat2)%>%summarise(sumAbun=sum(sizecat_abundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  mutate(period=ifelse(as.numeric(substr(id,1,4))<=1997,"GL","New"))%>%
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(period,pdscinam,sizecat2,id,sizecat_abundance)%>%distinct()%>%
  group_by(period,pdscinam,sizecat2)%>%summarise(sumAbun_nEmpty=sum(sizecat_abundance))


trawlDietSum<-trawlDiets%>%
  mutate(period=ifelse(as.numeric(substr(id,1,4))<=1997,"GL","New"))%>%
  ungroup()%>%
  #group_by(geoarea,est_year,season)%>% #NOT grouping, but add if you are
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(period,pdscinam,pdcomnam,sizecat2,nTrawls,nDiets,sumAbun,sumAbun_nEmpty,GL_pyscinam)%>%
  reframe(Wk=sum(sizecat_abundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(sizecat_abundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(sizecat_abundance*pik)/sumAbun #Frequency of occurrence, doesn't need to remove empties
  )%>%distinct() #removed the SD, weren't using anyway


check<-trawlDietSum%>%group_by(period,pdscinam,sizecat2)%>%summarise(N=sum(Wk,na.rm=T))
sum(check$N<0&check$N>0) #When this is 0, then you have all your means correct because they sum to 1 (or they're all empties)


# Clustering -----------------------
#### Comparing the two time periods ####

#Only want to use those species-sizecats that are in both time periods
GLclusts<-read.csv(here("data","raw", "GarrisonLink_clusters.csv"))%>%
  mutate(guildCol=brewer.pal(6,"Set1")[guild],
         species=str_to_title(species),
         guildName=case_when(guild==1~"Crab Eaters",
                             guild==2~"Planktivores",
                             guild==3~"Amphipod/Shrimp Eaters",
                             guild==4~"Shrimp/Small Fish Eaters",
                             guild==5~"Benthivores",
                             guild==6~"Piscivores"))

species_sizecats_inGLtime<-preyTrawls%>%
  filter(year<=1997
         & !is.na(GL_pyscinam))%>%
  mutate(species=str_to_title(pdcomnam))%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  reframe(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_sizecats_inNEWtime<-preyTrawls%>%
  filter(year>1997
         & !is.na(GL_pyscinam))%>%
  mutate(species=str_to_title(pdcomnam))%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  reframe(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_sizecats_inBOTH<-inner_join(dplyr::select(species_sizecats_inGLtime,species,sizecat2)%>%
                                      distinct(),
                                    dplyr::select(species_sizecats_inNEWtime,species,sizecat2)%>%
                                      distinct())


###### Cropping down the Original time, and clustering again #####
shared_props_gl<-species_sizecats_inGLtime%>%
  filter(paste(species,sizecat2) %in% paste(species_sizecats_inBOTH$species,species_sizecats_inBOTH$sizecat2))

shared_props_gl%>%
  filter(prop_w>0)%>%
  group_by(GL_pyscinam)%>%
  reframe(prop=max(prop_w))%>%
  print(n=52)

#Following through to get the "real" overlap matrix
props1_gl<-shared_props_gl[,-max(ncol(shared_props_gl))] #Cut off the volume prop for cleanliness
colnames(props1_gl)<-c("species1","sizecat1","GL_pyscinam","prop_w1")
props1_gl<-complete(props1_gl,species1,sizecat1,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species1,sizecat1)%>%mutate(t=sum(prop_w1,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props1_gl[is.na(props1_gl)]<-0
props2_gl<-shared_props_gl[,-max(ncol(shared_props_gl))]
colnames(props2_gl)<-c("species2","sizecat2","GL_pyscinam","prop_w2")
props2_gl<-complete(props2_gl,species2,sizecat2,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species2,sizecat2)%>%mutate(t=sum(prop_w2,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props2_gl[is.na(props2_gl)]<-0

overlap_shared_gl<-full_join(props1_gl,props2_gl,relationship="many-to-many")%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  reframe(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat_shared_gl<-pivot_wider(overlap_shared_gl,
                                   id_cols=c(species1,sizecat1),
                                   names_from=c(species2,sizecat2),
                                   values_from = s_do)
overlap_mat_shared_gl<-as.matrix(overlap_mat_shared_gl[,3:ncol(overlap_mat_shared_gl)])
#rownames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#colnames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#overlap_mat[overlap_mat==1]<-NA
for (i in 1:nrow(overlap_mat_shared_gl)) {
  for (j in 1:ncol(overlap_mat_shared_gl)) {
    overlap_mat_shared_gl[i,j]<-ifelse(j>i,NA,overlap_mat_shared_gl[i,j])
  }
}


#Bootstrapping
propMat_gl<-shared_props_gl%>%
  pivot_wider(id_cols=c(species,sizecat2),names_from = GL_pyscinam,values_from = prop_w)%>% #Species are rows, prey are columns. Flip id and names if need opposite
  mutate(species=paste(species,sizecat2,sep="_"))%>%
  select(-c("Empty","sizecat2"))
propMat_gl<-select(propMat_gl,species,order(colnames(propMat_gl)))
propMat_gl[is.na(propMat_gl)]<-0




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
bootDiffs_gl<-numeric()
while (reps<=nreps) {
  bootDiffs_gl<-c(bootDiffs_gl,resample(propMat_gl)) #resampling across all the zeroes too is a lower significance
  reps
  reps=reps+1
}
hist(bootDiffs_gl)
sigGuild_gl<-quantile(bootDiffs_gl,probs=0.95)
abline(v=sigGuild_gl,col="red",lwd=2)


#Visual of the actual matrix with significance indicated
#library(plot.matrix)
#par(mar=c(6,6,5,5.5))
#plot(overlap_mat,axis.row=list(side=2,las=1),col=viridis::viridis(n=100,option="B"),
#     polygon.key = list(border=NA), key=list(),xlab="",ylab="")

#Better visual (ggplot as always)
as.data.frame(overlap_mat_shared_gl)%>%
  mutate(species1=colnames(overlap_mat_shared_gl))%>%
  pivot_longer(cols=colnames(overlap_mat_shared_gl),names_to = "species2", values_to = "s_do")%>%
  mutate(species1=sort(species1,decreasing=T),species2=species2,
         sig=ifelse(s_do>sigGuild_gl,"bold.italic","plain"))%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=s_do))+
  geom_text(aes(species1,species2,label=round(s_do,digits=2),fontface=sig,color=sig),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100],
                              viridis::viridis(100,option="B")[sigGuild_gl*100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  ggtitle("1973-1997")+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

#Mean overlap, how similar are the diets in this time period
mean(overlap_mat_shared_gl[which(overlap_mat_shared_gl<1)],na.rm=T)
sd(overlap_mat_shared_gl[which(overlap_mat_shared_gl<1)],na.rm=T)
range(overlap_mat_shared_gl[which(overlap_mat_shared_gl<1)],na.rm=T)
hist(overlap_mat_shared_gl[which(overlap_mat_shared_gl<1)],main="All species 1973-1997",xlab="Dietary Overlap")
abline(v=sigGuild_gl,col="red",lwd=2)

#Clustering these out
overlap_shared_clust_gl<-hclust(as.dist(1-overlap_mat_shared_gl),method="average")
#Trying dendextend to pretty up the dendrogram
dend<-as.dendrogram(overlap_shared_clust_gl)


##### Cropping down the new time and clustering again #####
shared_props_new<-species_sizecats_inNEWtime%>%
  filter(paste(species,sizecat2) %in% paste(species_sizecats_inBOTH$species,species_sizecats_inBOTH$sizecat2))

shared_props_new%>%
  filter(prop_w>0)%>%
  group_by(GL_pyscinam)%>%
  reframe(prop=max(prop_w))%>%
  print(n=52)

#Following through to get the "real" overlap matrix
props1_new<-shared_props_new[,-max(ncol(shared_props_new))] #Cut off the volume prop for cleanliness
colnames(props1_new)<-c("species1","sizecat1","GL_pyscinam","prop_w1")
props1_new<-complete(props1_new,species1,sizecat1,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species1,sizecat1)%>%mutate(t=sum(prop_w1,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props1_new[is.na(props1_new)]<-0
props2_new<-shared_props_new[,-max(ncol(shared_props_new))]
colnames(props2_new)<-c("species2","sizecat2","GL_pyscinam","prop_w2")
props2_new<-complete(props2_new,species2,sizecat2,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species2,sizecat2)%>%mutate(t=sum(prop_w2,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props2_new[is.na(props2_new)]<-0

overlap_shared_new<-full_join(props1_new,props2_new)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat_shared_new<-pivot_wider(overlap_shared_new,
                                    id_cols=c(species1,sizecat1),
                                    names_from=c(species2,sizecat2),
                                    values_from = s_do)
overlap_mat_shared_new<-as.matrix(overlap_mat_shared_new[,3:ncol(overlap_mat_shared_new)])
#rownames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#colnames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#overlap_mat[overlap_mat==1]<-NA
for (i in 1:nrow(overlap_mat_shared_new)) {
  for (j in 1:ncol(overlap_mat_shared_new)) {
    overlap_mat_shared_new[i,j]<-ifelse(j>i,NA,overlap_mat_shared_new[i,j])
  }
}


#Bootstrapping
propMat_new<-shared_props_new%>%
  pivot_wider(id_cols=c(species,sizecat2),names_from = GL_pyscinam,values_from = prop_w)%>% #Species are rows, prey are columns. Flip id and names if need opposite
  mutate(species=paste(species,sizecat2,sep="_"))%>%
  select(-c("Empty","Miscellaneous","sizecat2"))
propMat_new<-select(propMat_new,species,order(colnames(propMat_new)))
propMat_new[is.na(propMat_new)]<-0


reps=1
nreps=250
bootDiffs_new<-numeric()
while (reps<=nreps) {
  bootDiffs_new<-c(bootDiffs_new,resample(propMat_new))
  print(reps)
  reps=reps+1
}
hist(bootDiffs_new)
sigGuild_new<-quantile(bootDiffs_new,probs=0.95)
abline(v=sigGuild_new,col="red",lwd=2)

#Visual of the actual matrix with significance indicated
#library(plot.matrix)
#par(mar=c(6,6,5,5.5))
#plot(overlap_mat,axis.row=list(side=2,las=1),col=viridis::viridis(n=100,option="B"),
#     polygon.key = list(border=NA), key=list(),xlab="",ylab="")

#Better visual (ggplot as always)
as.data.frame(overlap_mat_shared_new)%>%
  mutate(species1=colnames(overlap_mat_shared_new))%>%
  pivot_longer(cols=colnames(overlap_mat_shared_new),names_to = "species2", values_to = "s_do")%>%
  mutate(species1=sort(species1,decreasing=T),species2=species2,
         sig=ifelse(s_do>sigGuild_new,"bold.italic","plain"))%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=s_do))+
  geom_text(aes(species1,species2,label=round(s_do,digits=2),fontface=sig,color=sig),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100],viridis::viridis(100,option="B")[sigGuild_new*100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

#Mean overlap, how similar are the diets in this time period
mean(overlap_mat_shared_new[which(overlap_mat_shared_new<1)],na.rm=T)
sd(overlap_mat_shared_new[which(overlap_mat_shared_new<1)],na.rm=T)
range(overlap_mat_shared_new[which(overlap_mat_shared_new<1)],na.rm=T)
hist(overlap_mat_shared_new[which(overlap_mat_shared_new<1)],main="All Species 1998-2019",xlab="Dietary Overlap")
abline(v=sigGuild_new,col="red",lwd=2)
par(mar=c(5,5,5,5))
#Clustering these out
overlap_clust_shared_new<-hclust(as.dist(1-overlap_mat_shared_new),method="average")
#Trying dendextend to pretty up the dendrogram
dend_new<-as.dendrogram(overlap_clust_shared_new)


##### Guilds in the two times #####

#Cutting at the significance level for "clusters"
myClusts<-as.data.frame(cutree(overlap_shared_clust_gl,h=1-sigGuild_gl))
myClusts$species=rownames(myClusts)
colnames(myClusts)<-c("myCluster_1997","species")
myClusts<-myClusts%>%separate(species,into=c("species","sizecat"),sep="_")
#Cutting at a reasonable level to group those a bit higher for "guilds"
myGuilds<-as.data.frame(cutree(overlap_shared_clust_gl,k=5))
myGuilds$species=rownames(myGuilds)
colnames(myGuilds)<-c("myGuild_1997","species")
myGuilds<-myGuilds%>%separate(species,into=c("species","sizecat"),sep="_")
#Pasting those together in a way that keeps "cluster" as a subset within "guild"
myGuildClust<-full_join(myGuilds,myClusts)%>%
  arrange(myGuild_1997,myCluster_1997)
myGuildClust$clust<-1
for (i in 2:nrow(myGuildClust)) {
  if(myGuildClust$myGuild_1997[i-1]!=myGuildClust$myGuild_1997[i]) 
  {myGuildClust$clust[i]<-1}
  else if(myGuildClust$myCluster_1997[i-1]!=myGuildClust$myCluster_1997[i])
  {myGuildClust$clust[i]<-myGuildClust$clust[i-1]+1}
  else 
  {myGuildClust$clust[i]<-myGuildClust$clust[i-1]}
}
myGuildClust<-myGuildClust%>%mutate(clust=letters[clust],
                                    guildclust_1997=paste0(myGuild_1997,clust))%>%
  dplyr::select(-clust)
#Combining with GL actual clusters and guilds
myGLclusts<-full_join(myGuildClust,GLclusts)
#Simplifying for labels on the dendrograms
myGLclust_dend<-filter(myGLclusts,paste(species,sizecat,sep="_") %in% c(dend%>%labels))
myGLclust_dend<-myGLclust_dend[match(c(dend%>%labels), paste(myGLclust_dend$species,myGLclust_dend$sizecat,sep="_")),]


#Cutting at the significance level for "clusters"
myClusts_new<-as.data.frame(cutree(overlap_clust_shared_new,h=1-sigGuild_new))
myClusts_new$species=rownames(myClusts_new)
colnames(myClusts_new)<-c("myCluster_2019","species")
myClusts_new<-myClusts_new%>%separate(species,
                                      into=c("species","sizecat"),sep="_")
#Cutting at a reasonable level to group those a bit higher for "guilds"
myGuilds_new<-as.data.frame(cutree(overlap_clust_shared_new,k=13)) #This is at a sigLevel of ~0.45
myGuilds_new$species=rownames(myGuilds_new)
colnames(myGuilds_new)<-c("myGuild_2019","species")
myGuilds_new<-myGuilds_new%>%separate(species,
                                      into=c("species","sizecat"),sep="_")
#Pasting those together in a way that keeps "cluster" as a subset within "guild"
myGuilds_new<-full_join(myGuilds_new,myClusts_new)
myGuilds_new<-arrange(myGuilds_new,myGuild_2019,myCluster_2019)
myGuilds_new$clust<-1
for (i in 2:nrow(myGuilds_new)) {
  if(myGuilds_new$myGuild_2019[i-1]!=myGuilds_new$myGuild_2019[i]) 
  {myGuilds_new$clust[i]<-1}
  else if(myGuilds_new$myCluster_2019[i-1]!=myGuilds_new$myCluster_2019[i])
  {myGuilds_new$clust[i]<-myGuilds_new$clust[i-1]+1}
  else 
  {myGuilds_new$clust[i]<-myGuilds_new$clust[i-1]}
}
myGuilds_new<-myGuilds_new%>%mutate(clust=letters[clust],
                                    guildclust_2019=paste0(myGuild_2019,clust))%>%
  dplyr::select(-clust)
#Combining with GL actual clusters and guilds
myGLclusts_new<-full_join(myGuilds_new,GLclusts)
#Simplifying for labels on the dendrograms
myGLclust_new_dend<-filter(myGLclusts_new,paste(species,sizecat,sep="_") %in% c(dend_new%>%labels))
myGLclust_new_dend<-myGLclust_new_dend[match(c(dend_new%>%labels), paste(myGLclust_new_dend$species,myGLclust_new_dend$sizecat,sep="_")),]

#Just down to those that I have data for
compareClusts<-full_join(myGLclusts,myGLclusts_new)
sharedClusts<-filter(compareClusts,!is.na(myCluster_1997)
                     & !is.na(myCluster_2019))
#Naming
sharedClusts<-sharedClusts%>%
  mutate(nameClust_1997=case_when(guildclust_1997=="1a"~"Shrimp eaters",
                                  guildclust_1997=="1b"~"Polychaete/Amphipod eaters",
                                  guildclust_1997=="1c"~"Shrimp eaters",
                                  guildclust_1997=="1d"~"Polychaete/Amphipod eaters",
                                  guildclust_1997=="1e"~"Polychaete/Amphipod eaters",
                                  guildclust_1997=="2a"~"Benthivores",
                                  guildclust_1997=="2b"~"Benthivores",
                                  guildclust_1997=="3a"~"Piscivores",
                                  guildclust_1997=="3b"~"Crab eaters",
                                  guildclust_1997=="3c"~"Piscivores",
                                  guildclust_1997=="3d"~"Piscivores",
                                  guildclust_1997=="4a"~"Piscivores", #Really focused on squid, but compromise
                                  guildclust_1997=="5a"~"Piscivores",
                                  TRUE ~ "Unnamed"),
         nameClust_2019=case_when(myCluster_2019=="1"~"Shrimp eaters",
                                  myCluster_2019=="2"~"Polychaete/Amphipod eaters",
                                  myCluster_2019=="3"~"Piscivores",
                                  myCluster_2019=="4"~"Crab eaters",
                                  myCluster_2019=="5"~"Benthivores",
                                  TRUE ~ "Unnamed"))





##### Alluvial Plot showing changes in Cluster groups #####
# install.packages("ggalluvial")
library(ggalluvial)


sharedClusts%>%
  group_by(nameClust_1997,nameClust_2019)%>%
  summarise(freq=n())%>%
  mutate(nameClust_2019=as.character(nameClust_2019),
         nameClust_1997=factor(nameClust_1997,levels=c("Piscivores","Shrimp eaters","Crab eaters",
                                                       "Polychaete/Amphipod eaters","Benthivores","Not Assigned")),
         nameClust_2019=factor(nameClust_2019,levels=c("Piscivores","Shrimp eaters","Crab eaters",
                                                       "Polychaete/Amphipod eaters","Benthivores","Not Assigned"))) %>%
  ggplot(aes(axis1 = nameClust_1997, axis2 = nameClust_2019, y = freq)) +
  geom_alluvium(aes(fill = nameClust_1997),alpha=1)+
  geom_stratum(alpha=0.2,size=1.2)+
  geom_text(stat = "stratum",
            aes(label = gsub(" ","\n",after_stat(stratum)),color=after_stat(stratum)),fontface="bold",size=8.5) +
  scale_x_continuous(breaks=c(1,2),labels = c("1973-1997", "1998-2019"),
                     expand = c(0.15, 0.05))+
  scale_fill_viridis_d(option="D",end=1)+
  scale_color_manual(values=c("black","black","black","black","black"))+
  theme_void()+
  theme(legend.position = "none",axis.text.x=element_text(size=20,vjust=5),
        plot.title=element_text(size=30,hjust=0.5))+
  ggtitle(label="Guild Assignments in Each Time Period")


##Tanglegram ####
abs(order.dendrogram(dend_new)-rev(order.dendrogram(dend_new)))^1

#Doing the more robust comparisons (traditonally used) requires they have identical groups
dend_cor<-dendlist(dend,dend_new)%>%
  cor_cophenetic()
dend_tangle<-dendlist(dend,dend_new)%>%
  untangle(method="step2side")%>%
  entanglement()
dendlist(dend,dend_new)%>%
  untangle(method="step2side")%>%
  plot(common_subtrees_color_branches=T,highlight_distinct_edges=F,highlight_branches_lwd=F,
       edge.lwd=2,margin_inner=12,
       main=paste0("1973-1997       |       1998-2019"))


# NMDS but with guilds separated out --------------------------------------

#Total amount of catch of all species in each trawl
guildtrawlHauls<-preyTrawls%>%
  mutate(species=str_to_title(pdcomnam))%>%
  left_join(sharedClusts,by=c("sizecat2"="sizecat","species"="species"))%>%
  ungroup()%>%
  dplyr::select(id,species,nameClust_1997,sizecat_abundance)%>%distinct()%>%
  group_by(id,nameClust_1997)%>%
  reframe(trawlabundance=sum(sizecat_abundance,na.rm=T),
          trawlabundance=ifelse(trawlabundance==0,1,trawlabundance),
          nameClust_1997=fct_explicit_na(nameClust_1997,na_level="Not Assigned"))%>%distinct()


#Need to have whole diet totals, and whole trawl measures
NMDSguildtrawlDiets<-preyTrawls%>%
  mutate(species=str_to_title(pdcomnam))%>%
  left_join(sharedClusts,by=c("sizecat2"="sizecat","species"="species"))%>%
  group_by(id,nameClust_1997)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv),
         totalN=n_distinct(dietID),
         nameClust_1997=fct_explicit_na(nameClust_1997,na_level="Not Assigned"))%>%
  left_join(guildtrawlHauls)%>%
  group_by(id,GL_pyscinam)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(totalN==0,0,n_distinct(dietID)/totalN),
         doy=yday(est_towdate))%>%
  dplyr::select(year,id,nameClust_1997,trawlabundance,
                GL_pyscinam,totalwt,totalv,totalN,qikw,qikv,pik,decdeg_beglat:botsalin,doy,-est_towdate)%>%
  group_by(id)%>%
  fill(decdeg_beglat:doy,.direction="downup")%>%
  distinct()



#if grouping, add those groups in here both at select and group_by
nDiets_NMDS<-NMDSguildtrawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case the whole pred sizecat)
  ungroup()%>%
  dplyr::select(year,id,nameClust_1997,totalN)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(year,nameClust_1997)%>%summarise(nDiets=sum(totalN))#Calculate the number of diets
sumAbun_NMDS<-NMDSguildtrawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(year,id,nameClust_1997,trawlabundance)%>%distinct()%>%
  group_by(year,nameClust_1997)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_NMDS<-NMDSguildtrawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(year,id,nameClust_1997,trawlabundance)%>%distinct()%>%
  group_by(year,nameClust_1997)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))

NMDSguildyearsum<-NMDSguildtrawlDiets%>%
  group_by(year,nameClust_1997)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets_NMDS)%>%
  left_join(sumAbun_NMDS)%>%
  left_join(sumAbun_nEmpty_NMDS)%>%
  mutate(surftemp=weighted.mean(surftemp,trawlabundance,na.rm=T),
         bottemp=weighted.mean(bottemp,trawlabundance,na.rm=T),
         avgdepth=weighted.mean(avgdepth,trawlabundance,na.rm=T),
         surfsalin=weighted.mean(surfsalin,trawlabundance,na.rm=T),
         botsalin=weighted.mean(botsalin,trawlabundance,na.rm=T),
         doy=weighted.mean(doy,trawlabundance,na.rm=T))%>%
  group_by(year,nameClust_1997,nTrawls,nDiets,sumAbun,sumAbun_nEmpty,surftemp,bottemp,avgdepth,surfsalin,botsalin,doy,GL_pyscinam)%>%
  reframe(sumAbun_nEmpty=ifelse(is.na(sumAbun_nEmpty),nDiets,sumAbun_nEmpty),
            Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik)/sumAbun #Frequency of occurrence, doesn't need to remove empties
            #remove the SDs they don't work anyway
            )%>%
  distinct()

check<-NMDSguildyearsum%>%group_by(year,nameClust_1997)%>%summarise(N=sum(Wk))
sum(check$N<0&check$N>0) #When this is 0, then you have all your means correct because they sum to 1 (or they're empties)


idCols<-c("year","nameClust_1997","surftemp","bottemp","avgdepth","surfsalin","botsalin","doy")
NMDSguildmat<-pivot_wider(NMDSguildyearsum,id_cols=all_of(idCols),
                          names_from="GL_pyscinam",values_from="Wk")%>%
  dplyr::select(-c(Empty,"NA"))

library(vegan)
library(RVAideMemoire)
library(usedist)
library(indicspecies)

NMDSguildmat.diet<-(NMDSguildmat[,(length(idCols)+1):ncol(NMDSguildmat)])
#Percent zeroes--number of zeroes times dim of matrix
(sum(is.na(NMDSguildmat.diet))/(nrow(NMDSguildmat.diet)*ncol(NMDSguildmat.diet)))*100

#Filter out rows with no data
NMDSguildmat$sumRow<-rowSums(NMDSguildmat.diet,na.rm=T)
NMDSguildmat<-filter(NMDSguildmat,sumRow>0)

NMDSguildmat.diet<-(NMDSguildmat[,(length(idCols)+1):(ncol(NMDSguildmat)-1)])
NMDSguildmat.diet[is.na(NMDSguildmat.diet)]<-0
NMDSguildmat.env<-(NMDSguildmat[,1:length(idCols)])

#Using the whole dataset--but can use the one without minor species (Jellyfish, Sea Angel)
set.seed(42)


#removing the not assigned...
NMDSguildmat.wna<-filter(NMDSguildmat,nameClust_1997%notin%c("Not Assigned","Unnamed"))

NMDSguildmat.diet<-(NMDSguildmat.wna[,(length(idCols)+1):(ncol(NMDSguildmat.wna)-1)])
NMDSguildmat.diet[is.na(NMDSguildmat.diet)]<-0
NMDSguildmat.env<-(NMDSguildmat.wna[,1:length(idCols)])

#For the all diets NMDS
original.dist<-vegdist(NMDSguildmat.diet)
stress_values<-numeric(6)


for (n in 1:6) {
  nmds.resu <- metaMDS(NMDSguildmat.diet, k=n, distance = "bray", try=50, autotransform=F)
  stress_values[n]<-nmds.resu$stress
}
plot(stress_values, xlab="Number of axes", ylab="Stress",type="b")
abline(h=0.2,col="red")
stress_values
#So 3 axes is the first below 0.2
stress_differences<-as.numeric()
for (i in 1:length(stress_values)-1) {
  stress_differences[i]<-stress_values[i]-stress_values[i+1]
}
stress_differences
#So the drop in stress from 3 to 4 axes is the first less than 0.05
#So final 3 axes is best


#Go back and create the output for the 3 dimensions NMDS
set.seed(42)
guildNMDS<-metaMDS(NMDSguildmat.diet, distance = "bray", k = 3, try=250, autotransform=F ,na.rm = T)
r2<-summary(lm(as.vector(original.dist)~as.vector(dist(vegan::scores(guildNMDS)$sites))))[[8]]
actualStress<-guildNMDS$stress
stressplot(guildNMDS) #This is the visual of stress, the divergence of observed and ordinated distance. It's random, that's good

#Writing it out to use in PC-ORD for comparison and for the axes R2 values
#temp<-cbind(NMDSguildmat.env,guildNMDS)
#write.csv(temp,"/Users/nh1087/Documents/NMDS_Matrix.csv",row.names = F)



#### Print the species scores and sample scores ####
NMDS_guild_species<-as.data.frame(guildNMDS$species)%>%
  mutate(GL_comnam=rownames(.))
NMDS_guild_scores<-as.data.frame(guildNMDS$points)%>%
  bind_cols(NMDSguildmat.env)%>%
  mutate(nameClust_1997=factor(nameClust_1997,levels=c("Benthivores","Crab eaters",
                                                       "Piscivores","Polychaete/Amphipod eaters",
                                                       "Shrimp eaters","Not Assigned")),
         decade=paste0("'",substr(year,3,3),"0s"))
NMDS_guild_scores[is.na(NMDS_guild_scores)]<-NA

guildCentroids <- NMDS_guild_scores%>%
  group_by(nameClust_1997)%>%
  summarise(MDS1=mean(MDS1),
            MDS2=mean(MDS2),
            MDS3=mean(MDS3))%>%
  mutate(guildLabels=case_when(nameClust_1997=="Piscivores"~"P",
                               nameClust_1997=="Polychaete/Amphipod eaters"~"P/A",
                               nameClust_1997=="Shrimp eaters"~"S",
                               nameClust_1997=="Not Assigned"~"NA",
                               nameClust_1997=="Benthivores"~"B",
                               nameClust_1997=="Crab eaters"~"C",
                               TRUE~"NA"))


#### Analysing post-hoc ####
set.seed(42)
adonis2(original.dist~year*nameClust_1997,data=NMDS_guild_scores,permutations=1000,method="bray")

summary(NMDS_lm<-lm(cbind(MDS1,MDS2,MDS3)~year*nameClust_1997,data=NMDS_guild_scores))

predLine<-cbind(NMDS_guild_scores[,c(4,5)],data.frame(predict(NMDS_lm)))%>%
  group_by(nameClust_1997)%>%
  filter(year == min(year) | year == max(year))%>%
  mutate(time=ifelse(year==min(year),"start","stop"))%>%
  pivot_wider(id_cols=nameClust_1997,values_from=starts_with("MDS"),names_from=time)%>%
  mutate(guildLabels=case_when(nameClust_1997=="Piscivores"~"P",
                               nameClust_1997=="Polychaete/Amphipod eaters"~"P/A",
                               nameClust_1997=="Shrimp eaters"~"S",
                               nameClust_1997=="Not Assigned"~"NA",
                               nameClust_1997=="Benthivores"~"B",
                               nameClust_1997=="Crab eaters"~"C",
                               TRUE~"NA"),
         MDS1_label=MDS1_start-((MDS1_stop-MDS1_start)/10),
         MDS2_label=MDS2_start-((MDS2_stop-MDS2_start)/10),
         MDS3_label=MDS3_start-((MDS3_stop-MDS3_start)/10))


##### Vegan vector fitting ####
fit.guild<-envfit(guildNMDS~NMDS_guild_scores$nameClust_1997:NMDS_guild_scores$year,choices=1:3)
fit.guild
fit.guild.df<-data.frame(fit.guild[["vectors"]]$arrows)%>%
  rownames_to_column(var="var")

#Dispersion in different years
ydisp<-betadisper(original.dist,NMDS_guild_scores$year)
year_disp_guild<-data.frame(distance=ydisp[["distances"]],NMDSguildmat.wna[,1:2])
year_disp_guild$nameClust_1997<-factor(year_disp_guild$nameClust_1997,levels=c("Benthivores","Crab eaters","Piscivores",
                                                                               "Polychaete/Amphipod eaters","Shrimp eaters",
                                                                               "Not Assigned"))
year_dispersion<-data.frame(distance_to_median=ydisp$distances)%>%
  mutate(year=as.numeric(as.character(ydisp$group)))

ggplot()+
  geom_point(data=year_disp_guild,aes(year,distance,color=nameClust_1997),alpha=0.8,size=3)+
  geom_smooth(data=year_disp_guild,aes(year,distance,color=nameClust_1997,fill=nameClust_1997),method="lm",size=1)+
  geom_point(data=year_dispersion,aes(year,distance_to_median),size=5)+
  geom_smooth(data=year_dispersion,aes(year,distance_to_median),method="lm",color="black",size=2)+
  scale_color_manual(values=c(viridis::viridis(n=5),"grey"),name="Adult Feeding Guild",labels=c("Benthivores","Crab eaters","Piscivores",
                                                                                                "Polychaete/\nAmphipod eaters","Shrimp eaters",
                                                                                                "Not Assigned"))+
  scale_fill_manual(values=c(viridis::viridis(n=5),"grey"),name="Adult Feeding Guild",labels=c("Benthivores","Crab eaters","Piscivores",
                                                                                               "Polychaete/\nAmphipod eaters","Shrimp eaters",
                                                                                               "Not Assigned"))+
  xlab("Year")+
  scale_y_continuous(limits=c(0.16,0.8),name="Distance to Annual Centroid")+
  guides(fill=guide_legend(override.aes = list(alpha=0)),
         color=guide_legend(override.aes = list(linewidth=15),title.position="top",ncol=3))+
  theme(legend.key.spacing = unit(0,"pt"),legend.title=element_text(size=30,hjust=0.5),
        legend.position="top",legend.key = element_rect(fill="transparent"),
        legend.box.margin = margin(-10,0,-20,0))

summary(ydisp_lm<-lm(distance_to_median~year,data=year_dispersion))
summary(ygdisp_lm<-lm(distance~year*nameClust_1997,data=year_disp_guild))
(max(predict(ydisp_lm))-min(predict(ydisp_lm)))/max(predict(ydisp_lm))
gdisps<-bind_cols(year_disp_guild,smoothDist=predict(ygdisp_lm))%>%
  group_by(year)%>%
  mutate(mDist=mean(smoothDist))%>%
  group_by(nameClust_1997)%>%
  reframe(max=max(smoothDist),min=min(smoothDist),
          tmax=max(mDist),tmin=min(mDist))%>%
  mutate(d=(max-min)/max,
         td=(tmax-tmin)/tmax)


#Comparisons for each guild pairwise over time
pairwise_dist<-as.matrix(original.dist)
rownames(pairwise_dist)<-paste(NMDS_guild_scores$year,NMDS_guild_scores$nameClust_1997,sep="_")
colnames(pairwise_dist)<-paste(NMDS_guild_scores$year,NMDS_guild_scores$nameClust_1997,sep="_")
pairwise_dist<-as.data.frame(pairwise_dist)%>%
  pivot_longer(cols=matches("[0-9]"),names_to="guild_2",values_to="distance")%>%
  bind_cols(as.vector(sapply(paste(NMDS_guild_scores$year,NMDS_guild_scores$nameClust_1997,sep="_"), 
                             function (x) rep(x,nrow(NMDS_guild_scores)))))%>%
  mutate(guild_1=...3)%>%
  dplyr::select(-"...3")%>%
  filter(guild_1!=guild_2)%>%
  separate(guild_1,into=c("year_1","guild_1"),sep="_")%>%
  separate(guild_2,into=c("year_2","guild_2"),sep="_")%>%
  mutate(pairwise=paste(guild_1,guild_2,sep="--"),
         pair_rev=paste(guild_2,guild_1,sep="--"),
         year=as.numeric(year_1))%>%
  filter(year_1==year_2)
pairwise_dist<-pairwise_dist[!duplicated(lapply(as.data.frame(t(pairwise_dist), stringsAsFactors=FALSE), sort)),]
n_distinct(pairwise_dist$pairwise)
ggplot(data=pairwise_dist,aes(year,distance))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~pairwise,nrow=5)
summary(pairDist<-lm(distance~pairwise*year,data=pairwise_dist))



#### Plotting NMDS ####
axes12<-ggplot()+
  geom_point(data=NMDS_guild_scores,
             aes(MDS1,MDS2,color=nameClust_1997,shape=nameClust_1997),size=3,stroke=2,alpha=1)+
  #geom_segment(data=fit.guild.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
  #             color="blue",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  #geom_label(data=fit.guild.df,aes(x=NMDS1*1.1,y=NMDS2*1.1,label=var))+
  #geom_polygon(data=NMDS_guild_scores%>%group_by(nameClust_1997)%>%slice(chull(MDS1,MDS2)),
  #             aes(x=MDS1,y=MDS2,fill=nameClust_1997,color=nameClust_1997),alpha=0.1,lwd=1.5,show.legend=F)+
  #geom_segment(data=predLine,aes(x=MDS1_start,xend=MDS1_stop,y=MDS2_start,yend=MDS2_stop,group=nameClust_1997),
  #             size=4.5,arrow=arrow(length=unit(0.1,"inches"),angle=45),lineend="round",linejoin="round",show.legend=F)+
  #geom_segment(data=predLine,aes(x=MDS1_start,xend=MDS1_stop,y=MDS2_start,yend=MDS2_stop,color=nameClust_1997),
  #             size=3,arrow=arrow(length=unit(0.1,"inches"),angle=45),lineend="round",linejoin="round",show.legend=F)+
  #geom_point(data=guildCentroids,aes(MDS1,MDS2,color=nameClust_1997),
  #           size=12,stroke=2,show.legend=F,shape=21)+
#geom_label(data=predLine,aes(MDS1_label,MDS2_label,label=guildLabels,color=nameClust_1997),
#           size=7,fontface="bold",show.legend=F)+
#geom_text(aes(x=0.25,y=1.99,label=paste0("Stress == ",round(actualStress,3))),hjust=0,size=8,parse=T)+
#geom_text(aes(x=0.25,y=1.55,label=paste0("R^{2} == ",round(r2,3))),hjust=0,size=8,parse=T)+
#geom_segment(data=NMDS_guild_species,aes(x=0,xend=MDS1,y=0,yend=MDS2))+
#geom_text(data=NMDS_guild_species,aes(x=MDS1,y=MDS2,label=GL_comnam))+
scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c("-2","-1","0","1","2"),limits=c(-1.5,1.5))+
  scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c("-2","-1","0","1","2"),limits=c(-2.5,2.5))+
  scale_fill_manual(values=viridis::viridis_pal()(5),name="Adult Feeding Guild",guide="none")+
  scale_color_manual(values=viridis::viridis_pal()(5),name="Adult Feeding Guild")+
  scale_shape_manual(values=c(0:2,5,6),guide="none")+
  guides(color=guide_legend(override.aes = list(shape=c(0:2,5,6),size=6,alpha=1,stroke=2),
                            ncol=3,title.position="top"))+
  theme(axis.text=element_text(size=28),axis.title=element_text(size=33),
        legend.text=element_text(size=23),legend.title=element_text(size=33))+
  xlab("Axis 1")+ylab("Axis 2")


axes13<-ggplot()+
  geom_point(data=NMDS_guild_scores,
             aes(MDS1,MDS3,color=nameClust_1997,shape=nameClust_1997),size=3,stroke=2,alpha=1)+
  #geom_segment(data=fit.guild.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS3),
  #             color="blue",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  #geom_label(data=fit.guild.df,aes(x=NMDS1*1.1,y=NMDS3*1.1,label=var))+
  #geom_polygon(data=NMDS_guild_scores%>%group_by(nameClust_1997)%>%slice(chull(MDS1,MDS3)),
  #             aes(x=MDS1,y=MDS3,fill=nameClust_1997,color=nameClust_1997),alpha=0.1,lwd=1.5,show.legend=F)+
  #geom_segment(data=predLine,aes(x=MDS1_start,xend=MDS1_stop,y=MDS3_start,yend=MDS3_stop,group=nameClust_1997),
  #             size=4.5,arrow=arrow(length=unit(0.1,"inches"),angle=45),lineend="round",linejoin="round",show.legend=F)+
  #geom_segment(data=predLine,aes(x=MDS1_start,xend=MDS1_stop,y=MDS3_start,yend=MDS3_stop,color=nameClust_1997),
  #             size=3,arrow=arrow(length=unit(0.1,"inches"),angle=45),lineend="round",linejoin="round",show.legend=F)+
  #geom_point(data=guildCentroids,aes(MDS1,MDS3,color=nameClust_1997),
  #           size=12,stroke=2,show.legend=F,shape=21)+
#geom_label(data=predLine,aes(MDS1_label,MDS3_label,label=guildLabels,color=nameClust_1997),
#           size=7,fontface="bold",show.legend=F)+
#geom_segment(data=NMDS_guild_species,aes(x=0,xend=MDS1,y=0,yend=MDS3))+
#geom_text(data=NMDS_guild_species,aes(x=MDS1,y=MDS3,label=GL_comnam))+
scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c("-2","-1","0","1","2"),limits=c(-1.5,1.5))+
  scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c("-2","-1","0","1","2"),limits=c(-1.2,1.55))+
  scale_fill_manual(values=viridis::viridis_pal()(5),name="Adult Feeding Guild",guide="none")+
  scale_color_manual(values=viridis::viridis_pal()(5),name="Adult Feeding Guild")+
  scale_shape_manual(values=c(0:2,5,6),guide="none")+
  guides(color=guide_legend(override.aes = list(shape=c(0:2,5,6),size=6,alpha=1,stroke=2),
                            title.position="top",ncol=3))+
  theme(axis.text=element_text(size=28),axis.title=element_text(size=33),
        legend.text=element_text(size=23),legend.title=element_text(size=33))+
  xlab("Axis 1")+ylab("Axis 3")


axes23<-ggplot()+
  geom_point(data=NMDS_guild_scores,
             aes(MDS2,MDS3,color=nameClust_1997,shape=nameClust_1997),size=3,stroke=2,alpha=1)+
geom_segment(data=fit.guild.df,aes(x=0,xend=NMDS2,y=0,yend=NMDS3),
            color="blue",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
# geom_label(data=fit.guild.df,aes(x=NMDS2*1.1,y=NMDS2*1.1,label=var))+
geom_polygon(data=NMDS_guild_scores%>%group_by(nameClust_1997)%>%slice(chull(MDS2,MDS3)),
            aes(x=MDS2,y=MDS3,fill=nameClust_1997,color=nameClust_1997),alpha=0.1,lwd=1.5,show.legend=F)+
geom_segment(data=predLine,aes(x=MDS2_start,xend=MDS2_stop,y=MDS3_start,yend=MDS3_stop,group=nameClust_1997),
            size=4.5,arrow=arrow(length=unit(0.1,"inches"),angle=45),lineend="round",linejoin="round",show.legend="none")+
geom_segment(data=predLine,aes(x=MDS2_start,xend=MDS2_stop,y=MDS3_start,yend=MDS3_stop,color=nameClust_1997),
            size=3,arrow=arrow(length=unit(0.1,"inches"),angle=45),lineend="round",linejoin="round",show.legend="none")+
geom_point(data=guildCentroids,aes(MDS2,MDS3,color=nameClust_1997),
          size=12,stroke=2,alpha=0.7,show.legend=F,shape=21)+
# geom_label(data=predLine,aes(MDS2_label,MDS3_label,label=guildLabels,color=nameClust_1997),
          # size=7,fontface="bold",show.legend=F)+
# geom_segment(data=NMDS_guild_species,aes(x=0,xend=MDS2,y=0,yend=MDS3))+
# geom_text(data=NMDS_guild_species,aes(x=MDS2,y=MDS3,label=GL_comnam))+
scale_x_continuous(breaks=c(-2,-1,0,1,2),labels=c("-2","-1","0","1","2"),limits=c(-2.5,2.5))+
  scale_y_continuous(breaks=c(-2,-1,0,1,2),labels=c("-2","-1","0","1","2"),limits=c(-1.2,1.55))+
  scale_fill_manual(values=viridis::viridis_pal()(5),name="Adult Feeding Guild",guide="none")+
  scale_color_manual(values=viridis::viridis_pal()(5),name="Adult Feeding Guild")+
  scale_shape_manual(values=c(0:2,5,6),guide="none")+
  guides(color=guide_legend(override.aes = list(shape=c(0:2,5,6),size=6,alpha=1,stroke=2),
                            title.position="top",ncol=3))+
  theme(axis.text=element_text(size=28),axis.title=element_text(size=33),
        legend.text=element_text(size=23),legend.title=element_text(size=33))+
  xlab("Axis 2")+ylab("Axis 3")

png(filename="../Figures/Presentations/c3_nmds_guild_yeararrows_firstpoints.png",height=500,width=1400)
ggarrange(axes12,axes13,axes23,common.legend = T,legend="top",nrow=1)
dev.off()


### species accumulation curves, Testing which species should be included -----------
idTrawls_wide<-preyTrawls%>%
  filter(GL_pycomnam!="Empty")%>%
  pivot_wider(id_cols=c(dietID,pdcomnam,year,season),
              names_from=GL_pycomnam,values_from=pyamtw,values_fn=sum)
idTrawls_wide[is.na(idTrawls_wide)] <- 0
species_accum<-subset(species_com,toupper(species_com)%in%idTrawls_wide$pdcomnam)
species_accumulation_rates<-data.frame(species=species_accum,
                                       totalDiets=numeric(length=length(species_accum)),
                                       max_richness=numeric(length=length(species_accum)),
                                       nDiets_10=numeric(length=length(species_accum)),
                                       nDiets_50=numeric(length=length(species_accum)),
                                       nDiets_75=numeric(length=length(species_accum)),
                                       nDiets_80=numeric(length=length(species_accum)),
                                       nDiets_90=numeric(length=length(species_accum)),
                                       nDiets_99=numeric(length=length(species_accum)))
for (i in 1:length(species_accum)) {
  s<-specaccum(filter(idTrawls_wide,pdcomnam==toupper(species_accum)[i])[,5:ncol(idTrawls_wide)])
  species_accumulation_rates$totalDiets[i]<-length(s[["richness"]])
  species_accumulation_rates$max_richness[i]<-max(s[["richness"]])
  species_accumulation_rates$nDiets_10[i]<-Position(function(x) x > 0.10*max(s[["richness"]]), s[["richness"]])
  species_accumulation_rates$nDiets_50[i]<-Position(function(x) x > 0.50*max(s[["richness"]]), s[["richness"]])
  species_accumulation_rates$nDiets_75[i]<-Position(function(x) x > 0.75*max(s[["richness"]]), s[["richness"]])
  species_accumulation_rates$nDiets_80[i]<-Position(function(x) x > 0.80*max(s[["richness"]]), s[["richness"]])
  species_accumulation_rates$nDiets_90[i]<-Position(function(x) x > 0.90*max(s[["richness"]]), s[["richness"]])
  species_accumulation_rates$nDiets_99[i]<-Position(function(x) x > 0.99*max(s[["richness"]]), s[["richness"]])
}
ggplot(species_accumulation_rates)+
  geom_point(aes(x=totalDiets,y=max_richness),size=5)+
  ylab("Maximum Diet Richness")+xlab("Total Diets Collected")

species_accumulation_rates%>%
  pivot_longer(cols=starts_with("nDiets_"),names_to="Level",values_to="richness")%>%
  filter(Level%in%c("nDiets_10","nDiets_50","nDiets_90"))%>%
  ggplot()+
  geom_point(aes(x=species,y=richness,color=Level),size=5)+
  scale_color_discrete(labels=c("10%","50%","90%"))+
  ylab("Number of Diets")+xlab("Species")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.33),legend.position=c(0.8,0.8))


species_accum_curves<-pivot_longer(species_accumulation_rates,cols=starts_with("nDiets"),
                                   names_to="n_prey",values_to="n_diets")%>%
  mutate(n_prey=case_when(n_prey=="nDiets_10" ~ 0.1 *max_richness,
                          n_prey=="nDiets_50" ~ 0.5 *max_richness,
                          n_prey=="nDiets_75" ~ 0.75*max_richness,
                          n_prey=="nDiets_80" ~ 0.8 *max_richness,
                          n_prey=="nDiets_90" ~ 0.9 *max_richness,
                          n_prey=="nDiets_99" ~ 0.99*max_richness))

ggplot(species_accum_curves)+
  geom_line(aes(x=n_diets/1000,y=n_prey))+
  geom_point(aes(x=totalDiets/1000,y=max_richness))+
  scale_x_continuous(name="Diets (thous)",labels=scales::number_format(accuracy=1))+
  scale_y_continuous(name="Prey Categories")+
  facet_wrap(~species,scales="fixed",nrow=8)



# Individual Species Analyses ---------------------------------------------

#Need to have whole diet totals, and whole trawl measures
trawlDiets<-preyTrawls_nSize%>%
  ungroup()%>%
  filter(GL_pyscinam!="Empty")%>%
  dplyr::select(id,pdcomnam,pdscinam,pyamtw,pyamtv,GL_pyscinam,dietID,year,season,nDiets,abundance)%>%
  group_by(id,pdcomnam)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv))%>%
  group_by(id,pdcomnam,GL_pyscinam)%>%
  mutate(preyTotal=sum(pyamtw),
         qikw=ifelse(totalwt==0,0,preyTotal/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets))%>%
  dplyr::select(year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,abundance,pdcomnam,pdscinam,GL_pyscinam,totalwt,totalv,nDiets,preyTotal,qikw,qikv,pik)%>% 
  distinct()%>% #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.
  group_by(id,pdcomnam)%>%
  mutate(test=sum(qikw))%>%
  ungroup()


#What's going on with trawl #1980072491130????
empty<-preyTrawls_nSize%>%
  group_by(id,pdcomnam)%>%
  mutate(nEmpty=sum(pynam=="EMPTY"),
         pik=nEmpty/nDiets,
         fullFreq=pik*abundance)%>%
  dplyr::select(year,season,id,pdcomnam,pdscinam,abundance,nDiets,nEmpty,pik,fullFreq)%>%distinct()

#empty<-filter(preyTrawls,GL_pyscinam=="Empty")%>%
#  group_by(id,pdcomnam)%>%
#  mutate(pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets))%>%
#  dplyr::select(id,pdcomnam,pdscinam,year,season,nDiets,pik)%>%distinct()%>%
#  right_join(abundance)%>%
#  mutate(pik=ifelse(is.na(pik),0,pik),
#         freq=pik*nDiets,
#         fullFreq=pik*abundance)


#c means that's only the stomachs with something in them (c=consumed)
trawlDiets_cind<-uniqueDiets%>%
  filter(pdwgt>0 & pdgutw>0)%>% #Needs a mass, and only looking at those fish that ate (since emptiness is modeled elsewhere)
  filter(dietID!="02720000030112011022571280")%>%
  left_join(abundance%>%mutate(pdcomnam=toupper(pdcomnam)))%>%
  group_by(id,pdcomnam,pdscinam)%>%
  mutate(relConsump=pdgutw/pdwgt,
         meanC=mean(relConsump,na.rm=T),
         meanWt=mean(pdgutw),
         meanV =mean(pdgutv),
         meanMass=mean(pdwgt),
         meanTL=mean(pdlen),
         abundance=ifelse(is.na(abundance),nDiets,abundance))%>%
  dplyr::select(year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                id,abundance,pdcomnam,pdscinam,meanWt,meanV,meanC,meanMass,meanTL,nDiets)%>% 
  distinct()%>%ungroup() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.



matDietSum_Wk<-trawlDiets%>%
  mutate(prop2=qikw^2,
         GL_pyscinam=paste0("Prey.",GL_pyscinam))%>%
  pivot_wider(id_cols=c(year,season,abundance,nDiets,totalwt,id,pdcomnam,pdscinam),
              names_from="GL_pyscinam",values_from="prop2")%>%
  right_join(dplyr::select(preyTrawls_nSize,year,season,abundance,nDiets,id,pdcomnam)%>%distinct())%>%
  mutate(richness=rowSums(select(.,starts_with("Prey"))>0,na.rm=T),
         sumP2=select(.,starts_with("Prey"))%>%rowSums(na.rm=T),
         simpson=1-sumP2,
         levin=1/sumP2,
         levin=ifelse(is.infinite(levin),NA,levin),
         levinStd=(levin-1)/(ncol(select(.,starts_with("Prey")))-1))


#mutate(comname=fct_reorder(factor(comname,levels=unique(trawlDietSum$comname)),levin,na.rm=T))
trawlDiets_breadth<-dplyr::select(matDietSum_Wk,year:pdscinam,richness:levinStd)%>%
  mutate(simpson=ifelse(richness==0,NA,simpson))



#Hare et al. 2016 expert rankings on species ecology
hare_qualScores<-read_csv(here("data","raw","S1Dataset.csv"))
colnames(hare_qualScores)<-gsub(" ","\\.",colnames(hare_qualScores))  
hare_qualScores<-hare_qualScores%>%
  mutate(Species=str_to_title(Species),
         Species=gsub("Sand Lances","Northern Sand Lance",Species),
         Species=gsub("Monkfish","Goosefish",Species),
         Species=gsub("Inshore ","",Species),
         totalScore=1*Low+2*Moderate+3*High+4*Very.High,
         meanScore=ifelse(Attribute.Category=="Sensitivity Attribute",totalScore/25,totalScore/20))
hare_directScores<-read_csv(here("data","raw","S2Dataset.csv"))
colnames(hare_directScores)<-gsub(" ","\\.",colnames(hare_directScores))
hare_directScores<-hare_directScores%>%
  mutate(Positive=as.numeric(gsub(",","",Positive)),
         Species=str_to_title(Species),
         Species=gsub("Sand Lances","Northern Sand Lance",Species),
         Species=gsub("Monkfish","Goosefish",Species),
         Species=gsub("Inshore ","",Species),
         Species=gsub("Menahden","Menhaden",Species),
         Species=gsub("Cleanrnose","Clearnose",Species),
         totalScore=-1*Negative+0*Neutral+1*Positive,
         meanScore=totalScore/12,
         Attribute="Climate direction effects",
         Attribute.Category="Climate direction")

hare_allScores<-bind_rows(hare_qualScores,hare_directScores)

unique(filter(allprey, str_to_title(pdcomnam) %notin% hare_qualScores$Species)$pdcomnam) #Almost all my fish are included, so that's good



#Getting their overall summary scores (the big main table)
hare_overallScores<-hare_allScores%>%
  group_by(Species,Attribute.Category)%>%
  summarise(overallScore=mean(meanScore),
            nVery.High=sum(meanScore>=3.5),
            nHigh=sum(meanScore>=3),
            nModerate=sum(meanScore>=2.5))%>%
  mutate(scoreClass=case_when(Attribute.Category=="Climate direction" & overallScore<=-0.333 ~ "Negative",
                              Attribute.Category=="Climate direction" & overallScore>= 0.333 ~ "Positive",
                              Attribute.Category=="Climate direction" & overallScore>-0.333 & overallScore<0.333 ~ "Neutral",
                              Attribute.Category%in%c("Sensitivity Attribute","Exposure Factor") & nVery.High>=3 ~ "Very High",
                              Attribute.Category%in%c("Sensitivity Attribute","Exposure Factor") & nHigh>=2 & nVery.High<3 ~ "High",
                              Attribute.Category%in%c("Sensitivity Attribute","Exposure Factor") & nModerate>=2 & nHigh<2 ~ "Moderate",
                              TRUE ~ "Low"),
         scoreClass_num=case_when(scoreClass=="Low"~1,
                                  scoreClass=="Moderate"~2,
                                  scoreClass=="High"~3,
                                  scoreClass=="Very High"~4,
                                  TRUE~0),
         scoreClass_num=ifelse(scoreClass_num==0,NA,scoreClass_num))%>%
  group_by(Species)%>%
  mutate(vulnerability=prod(scoreClass_num,na.rm=T),
         vulClass=case_when(vulnerability>=12~"Very High",
                            vulnerability>=8 & vulnerability<12~"High",
                            vulnerability>=4 & vulnerability<8~"Moderate",
                            vulnerability<4~"Low"))%>%dplyr::select(-c(nVery.High,nHigh,nModerate))


#Calculating the overall change potential for each species
hare_changePotential<-hare_allScores%>%
  filter(Attribute%in%c("Adult Mobility","Dispersal of Early Life Stages","Habitat Specificity","Sensitivity to Temperature"))%>%
  mutate(totalScore=ifelse(Attribute=="Sensitivity to Temperature",totalScore,4*Low+3*Moderate+2*High+1*Very.High),
         meanScore=totalScore/25)%>%
  group_by(Species)%>%
  summarise(overallScore=mean(meanScore),
            nVery.High=sum(meanScore>=3.5),
            nHigh=sum(meanScore>=3),
            nModerate=sum(meanScore>=2.5))%>%
  mutate(scoreClass=case_when(nVery.High>=3 ~ "Very High",
                              nHigh>=2 & nVery.High<3 ~ "High",
                              nModerate>=2 & nHigh<2 ~ "Moderate",
                              TRUE ~ "Low"),
         scoreClass_num=case_when(scoreClass=="Low"~1,
                                  scoreClass=="Moderate"~2,
                                  scoreClass=="High"~3,
                                  scoreClass=="Very High"~4,
                                  TRUE~0),
         scoreClass_num=ifelse(scoreClass_num==0,NA,scoreClass_num),
         scoreClass=factor(scoreClass,levels=c("Low","Moderate","High","Very High")))%>%
  dplyr::select(-c(nVery.High,nHigh,nModerate))


hare_everything<-hare_overallScores%>%
  pivot_wider(id_cols=Species,names_from=Attribute.Category,values_from=scoreClass)%>%
  left_join(dplyr::select(hare_overallScores,Species,Climate.Vulnerability=vulClass)%>%distinct())%>%
  left_join(dplyr::select(hare_changePotential,Species,Change.Potential=scoreClass))%>%
  mutate(Climate.Direction=factor(`Climate direction`,levels=c("Negative","Neutral","Positive")),
         Exposure.Factor=factor(`Exposure Factor`,levels=c("Low","Moderate","High","Very High")),
         Sensitivity.Attribute=factor(`Sensitivity Attribute`,levels=c("Low","Moderate","High","Very High")),
         Climate.Vulnerability=factor(Climate.Vulnerability,levels=c("Low","Moderate","High","Very High")))%>%
  dplyr::select(-c(`Climate direction`,`Sensitivity Attribute`,`Exposure Factor`))

hare_everything_FHD<-filter(hare_everything,Species %in% str_to_title(allprey$pdcomnam))

indResponses<-left_join(empty,trawlDiets_breadth)%>%
  left_join(trawlDiets_cind)%>%
  mutate(Species=str_to_title(pdcomnam),
         simpson2=simpson+0.0000001)%>% #to make the range (0,1) instead of [0,1]
  left_join(hare_everything_FHD)%>%
  mutate(Climate.Vulnerability=fct_explicit_na(Climate.Vulnerability,"Undetermined"),
         Climate.Direction=fct_explicit_na(Climate.Direction,"Undetermined"),
         Change.Potential=fct_explicit_na(Change.Potential,"Undetermined"))%>%
  left_join(preyTrawls_nSize%>%ungroup()%>%
              dplyr::select(id,lat=decdeg_beglat,lon=decdeg_beglon,geoarea,
                            est_towdate,month,day,avgdepth,surftemp,bottemp,surfsalin)%>%
              distinct(id,.keep_all=T)%>%
              mutate(DOY=yday(mdy(paste(month,day,substr(id,1,4),sep="-"))),
                     TOD=hour(est_towdate)))%>%
  mutate(season=factor(season,levels=c("Spring","Fall")),
         geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))

goodSpecies<-indResponses%>%
  group_by(Species,year,season)%>%
  summarise(nTotal=sum(nDiets))%>%
  group_by(Species,season)%>%
  mutate(nnTotals=(nTotal>20),
         nnTotals=sum(nnTotals))%>%
  filter(nnTotals>20)
n_distinct(paste(goodSpecies$Species,goodSpecies$season))
n_distinct(goodSpecies$Species)

indResponses.good<-filter(indResponses, paste(Species,season) %in% paste(goodSpecies$Species,goodSpecies$season))




#  Single summary plot of metrics over time -----------------------------------------------------------------


species_key<-arrange(species_key,pdcomnam)
species_key_good<-filter(species_key,pdcomnam %in% indResponses.good$pdcomnam)


sp.pEmpty.y<-glm(cbind(nEmpty,nDiets)~year,
                 data=indResponses.good,family=binomial(),weight=sqrt(abundance)+nDiets)
summary(sp.pEmpty.y)
pik.pred<-data.frame(year=filter(indResponses.good,!is.na(year))$year,
                     predict=predict.glm(sp.pEmpty.y,type="response"),
                     se=predict.glm(sp.pEmpty.y,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
pe.change<-(max(pik.pred$predict)-min(pik.pred$predict))*sign(sp.pEmpty.y$coefficients[2])

#Truly for all the species
sp.pEmpty.std<-glm(cbind(nEmpty,nDiets)~Species*year*season+surftemp+DOY,
                   data=filter(indResponses.good,!is.na(surftemp)),family=binomial())
summary(sp.pEmpty.std)

sp.pEmpty.s<-glm(cbind(nEmpty,nDiets)~Species*year*season,data=indResponses.good,family=binomial(),weight=sqrt(abundance)+nDiets)
summary(sp.pEmpty.s)
#Anova(sp.pEmpty.s,sp.pEmpty.std)

sp.plotResponses.pe<-data.frame(filter(indResponses.good,!is.na(abundance)),
                                pred=predict(sp.pEmpty.s,type="response",se=T))%>%
  mutate(se.lower=pred.fit-1.96*pred.se.fit,
         se.upper=pred.fit+1.96*pred.se.fit)%>%
  dplyr::select(Species,year,season,pred.fit,se.lower,se.upper)%>%distinct()
sp.responseCats.pe<-sp.plotResponses.pe%>%
  group_by(Species,season)%>%
  filter(year==min(year) | year==max(year))%>%
  mutate(year=ifelse(year==min(year),"Min","Max"))%>%
  pivot_wider(id_cols=c("Species","season"),names_from="year",values_from="pred.fit")%>%
  mutate(slope=Max*100-100*Min,
         pChange=Max/Min,
         responseCat=case_when(slope>=5 ~"Strong Positive",
                               slope<5 & slope>=1.666~"Weak Positive",
                               slope<=-1.666 & slope>-5~"Weak Negative",
                               slope<=-5~"Strong Negative",
                               TRUE~"Neutral"),
         responseCat2=case_when(pChange>=1.3~"Strong Positive",
                                pChange>=1.1 & pChange<1.3~"Weak Positive",
                                pChange<=0.9 & pChange>0.7~"Weak Negative",
                                pChange<=0.7 ~ "Strong Negative",
                                TRUE~"Neutral"),
         responseCat=factor(responseCat,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")),
         responseCat2=factor(responseCat2,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")))
table(sp.responseCats.pe$responseCat2)
n_distinct(filter(sp.responseCats.pe,responseCat2=="Strong Positive")$Species)
#Summary of change
median(sp.responseCats.pe$pChange-1)
range(sp.responseCats.pe$pChange-1)


#Pulling more info from the predicted GLM output
pe.lst <- lstrends(sp.pEmpty.s, c("Species","season"), var="year")
pe.lst<-as.data.frame(pe.lst)%>%
  mutate(psig=ifelse(sign(asymp.LCL)==sign(asymp.UCL),"*",""))%>%
  rename(species=Species)
pik.sp.pred<-data.frame(year=filter(indResponses.good,!is.na(year))$year,
                        season=filter(indResponses.good,!is.na(year))$season,
                        species=filter(indResponses.good,!is.na(year))$Species,
                        predict=predict.glm(sp.pEmpty.s,type="response"),
                        se=predict.glm(sp.pEmpty.s,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
pe.sp.change<-pik.sp.pred%>%
  arrange(species,season,year)%>%
  group_by(species,season)%>%
  slice(c(1,n()))%>%
  reframe(diff=predict-lag(predict))%>%
  filter(!is.na(diff))
lm.sum.pEmpty<-as.data.frame(summary(sp.pEmpty.s)$coefficients)%>%
  rename(p_val=`Pr(>|z|)`)%>%
  mutate(coef=rownames(summary(sp.pEmpty.s)$coefficients))%>%
  filter(grepl("[A-z ]+:[A-z ]+:[A-z ]",coef) | grepl("Species[A-z ]+:year",coef) | coef=="year" | coef=="year:seasonFall")%>%
  separate(coef,into=c("species","p_sig","season"),sep=":")%>%
  mutate(season=ifelse(grepl("season",species) | grepl("season",season),"Fall","Spring"),
         species=ifelse(grepl("Species",species),gsub("Species","",species),"Acadian Redfish"),
         season=ifelse(species %in% c("Smooth Skate","Bluefish"),"Fall",season),
         p_sig=case_when(p_val>0.1~"",
                         between(p_val,0.05,0.1)~".",
                         between(p_val,0.01,0.05)~"*",
                         between(p_val,0.001,0.01)~"**",
                         p_val<0.001~"***"))%>%
  left_join(filter(pik.sp.pred,year==max(year)))%>%
  mutate(year=year+1)%>%
  left_join(pe.sp.change)%>%
  left_join(pe.lst)



sp.dBreadth.y<-glm(simpson~year,data=filter(indResponses.good,!is.na(simpson)),family=quasibinomial(logit),weight=sqrt(abundance)+nDiets)
summary(sp.dBreadth.y)
db.pred<-data.frame(year=filter(indResponses.good,!is.na(simpson))$year,
                    predict=predict.glm(sp.dBreadth.y,type="response"),
                    se=predict.glm(sp.dBreadth.y,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
db.change<-(max(db.pred$predict)-min(db.pred$predict))*sign(sp.dBreadth.y$coefficients[2])

#Diet Breadth
sp.dBreadth<-glm(simpson~Species*year*season,data=indResponses.good,family=quasibinomial(logit),weight=sqrt(abundance)+nDiets)
summary(sp.dBreath)

sp.plotResponses.db<-data.frame(filter(indResponses.good,!is.na(simpson)&!is.na(abundance)),
                                pred=predict(sp.dBreadth,type="response",se=T))%>%
  mutate(se.lower=pred.fit-1.96*pred.se.fit,
         se.upper=pred.fit+1.96*pred.se.fit)%>%
  dplyr::select(Species,year,season,pred.fit,se.lower,se.upper)%>%distinct()
sp.responseCats.db<-sp.plotResponses.db%>%
  group_by(Species,season)%>%
  filter(year==min(year) | year==max(year))%>%
  mutate(year=ifelse(year==min(year),"Min","Max"))%>%
  pivot_wider(id_cols=c("Species","season"),names_from="year",values_from="pred.fit")%>%
  mutate(slope=Max*100-100*Min,
         pChange=Max/Min,
         responseCat=case_when(slope>=5 ~"Strong Positive",
                               slope<5 & slope>=1.666~"Weak Positive",
                               slope<=-1.666 & slope>-5~"Weak Negative",
                               slope<=-5~"Strong Negative",
                               TRUE~"Neutral"),
         responseCat2=case_when(pChange>=1.3~"Strong Positive",
                                pChange>=1.1 & pChange<1.3~"Weak Positive",
                                pChange<=0.9 & pChange>0.7~"Weak Negative",
                                pChange<=0.7 ~ "Strong Negative",
                                TRUE~"Neutral"),
         responseCat=factor(responseCat,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")),
         responseCat2=factor(responseCat2,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")))
table(sp.responseCats.db$responseCat2)
n_distinct(filter(sp.responseCats.db,responseCat2=="Strong Negative")$Species)
#Summary of change
range(sp.responseCats.db$pChange-1)
median(sp.responseCats.db$pChange-1)



db.lst <- lstrends(sp.dBreadth, c("Species","season"), var="year")
db.lst<-as.data.frame(db.lst)%>%
  mutate(psig=ifelse(sign(asymp.LCL)==sign(asymp.UCL),"*",""))%>%
  rename(species=Species)
db.sp.pred<-data.frame(year=filter(indResponses.good,!is.na(simpson))$year,
                       season=filter(indResponses.good,!is.na(simpson))$season,
                       species=filter(indResponses.good,!is.na(simpson))$Species,
                       predict=predict.glm(sp.dBreadth,type="response"),
                       se=predict.glm(sp.dBreadth,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
db.sp.change<-db.sp.pred%>%
  arrange(species,season,year)%>%
  group_by(species,season)%>%
  slice(c(1,n()))%>%
  reframe(diff=predict-lag(predict))%>%
  filter(!is.na(diff))
lm.sum.dBreadth<-as.data.frame(summary(sp.dBreadth)$coefficients)%>%
  rename(p_val=`Pr(>|t|)`)%>%
  mutate(coef=rownames(summary(sp.dBreadth)$coefficients))%>%
  filter(grepl("[A-z ]+:[A-z ]+:[A-z ]",coef) | grepl("Species[A-z ]+:year",coef) | coef=="year" | coef=="year:seasonFall")%>%
  separate(coef,into=c("species","p_sig","season"),sep=":")%>%
  mutate(season=ifelse(grepl("season",species) | grepl("season",season),"Fall","Spring"),
         species=ifelse(grepl("Species",species),gsub("Species","",species),"Acadian Redfish"),
         season=ifelse(species %in% c("Smooth Skate","Bluefish"),"Fall",season),
         p_sig=case_when(p_val>0.1~"",
                         between(p_val,0.05,0.1)~".",
                         between(p_val,0.01,0.05)~"*",
                         between(p_val,0.001,0.01)~"**",
                         p_val<0.001~"***"))%>%
  left_join(filter(db.sp.pred,year==max(year)))%>%
  mutate(year=year+1)%>%
  left_join(db.sp.change)%>%
  left_join(db.lst)


sp.rConsump.y<-glm(meanC~year,data=filter(indResponses.good,!is.na(meanC)),family=Gamma(link="log"),weight=sqrt(abundance)+nDiets)
summary(sp.rConsump.y)
rc.pred<-data.frame(year=filter(indResponses.good,!is.na(meanC))$year,
                    predict=predict.glm(sp.rConsump.y,type="response"),
                    se=predict.glm(sp.rConsump.y,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
rc.change<-(max(rc.pred$predict)-min(rc.pred$predict))*sign(sp.rConsump.y$coefficients[2])


#Relative Consumption
sp.rConsump<-glm(meanC~Species*year*season,data=indResponses.good,family=Gamma(link="log"),weight=sqrt(abundance)+nDiets)
summary(sp.rConsump)

sp.plotResponses.rc<-data.frame(filter(indResponses.good,!is.na(meanC)&!is.na(abundance)),
                                pred=predict(sp.rConsump,type="response",se=T))%>%
  mutate(se.lower=pred.fit-1.96*pred.se.fit,
         se.upper=pred.fit+1.96*pred.se.fit)%>%
  dplyr::select(Species,year,season,pred.fit,se.lower,se.upper)%>%distinct()
sp.responseCats.rc<-sp.plotResponses.rc%>%
  group_by(Species,season)%>%
  filter(year==min(year) | year==max(year))%>%
  mutate(year=ifelse(year==min(year),"Min","Max"))%>%
  pivot_wider(id_cols=c("Species","season"),names_from="year",values_from="pred.fit")%>%
  mutate(slope=Max*100-100*Min,
         pChange=Max/Min,
         responseCat=case_when(slope>=5 ~"Strong Positive",
                               slope<5 & slope>=1.666~"Weak Positive",
                               slope<=-1.666 & slope>-5~"Weak Negative",
                               slope<=-5~"Strong Negative",
                               TRUE~"Neutral"),
         responseCat2=case_when(pChange>=1.3~"Strong Positive",
                                pChange>=1.1 & pChange<1.3~"Weak Positive",
                                pChange<=0.9 & pChange>0.7~"Weak Negative",
                                pChange<=0.7 ~ "Strong Negative",
                                TRUE~"Neutral"),
         responseCat=factor(responseCat,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")),
         responseCat2=factor(responseCat2,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")))
table(sp.responseCats.rc$responseCat2)
n_distinct(filter(sp.responseCats.rc,responseCat2=="Strong Negative")$Species)
#Summary of change
median(sp.responseCats.rc$pChange-1)
range(sp.responseCats.rc$pChange-1)



rc.lst <- lstrends(sp.rConsump, c("Species","season"), var="year")
rc.lst<-as.data.frame(rc.lst)%>%
  mutate(psig=ifelse(sign(lower.CL)==sign(upper.CL),"*",""))%>%
  rename(species=Species)
rc.sp.pred<-data.frame(year=filter(indResponses.good,!is.na(meanC))$year,
                       season=filter(indResponses.good,!is.na(meanC))$season,
                       species=filter(indResponses.good,!is.na(meanC))$Species,
                       predict=predict.glm(sp.rConsump,type="response"),
                       se=predict.glm(sp.rConsump,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
rc.sp.change<-rc.sp.pred%>%
  arrange(species,season,year)%>%
  group_by(species,season)%>%
  slice(c(1,n()))%>%
  reframe(diff=predict-lag(predict))%>%
  filter(!is.na(diff))
lm.sum.relC<-as.data.frame(summary(sp.rConsump)$coefficients)%>%
  rename(p_val=`Pr(>|t|)`)%>%
  mutate(coef=rownames(summary(sp.rConsump)$coefficients))%>%
  filter(grepl("[A-z ]+:[A-z ]+:[A-z ]",coef) | grepl("Species[A-z ]+:year",coef) | coef=="year" | coef=="year:seasonFall")%>%
  separate(coef,into=c("species","p_sig","season"),sep=":")%>%
  mutate(season=ifelse(grepl("season",species) | grepl("season",season),"Fall","Spring"),
         species=ifelse(grepl("Species",species),gsub("Species","",species),"Acadian Redfish"),
         season=ifelse(species %in% c("Smooth Skate","Bluefish"),"Fall",season),
         p_sig=case_when(p_val>0.1~"",
                         between(p_val,0.05,0.1)~".",
                         between(p_val,0.01,0.05)~"*",
                         between(p_val,0.001,0.01)~"**",
                         p_val<0.001~"***"))%>%
  left_join(filter(rc.sp.pred,year==max(year)))%>%
  mutate(year=year+1)%>%
  left_join(rc.sp.change)%>%
  left_join(rc.lst)





#Plots
PIK<-ggplot()+
  #geom_point(data=indResponses.good,aes(x=year,y=pik,fill=Species,shape=season,size=abundance),alpha=0.2)+
  #geom_ribbon(data=pik.sp.pred,aes(x=year,ymin=min,ymax=max,fill=species,group=paste(season,species)),alpha=0.3)+
  geom_line(data=pik.sp.pred,aes(x=year,y=predict,color=species,group=paste(season,species),linetype=season),
            size=2,alpha=0.5)+
  geom_ribbon(data=pik.pred,aes(x=year,ymin=min,ymax=max),alpha=0.7,show.legend = F)+
  geom_line(data=pik.pred,aes(x=year,y=predict),size=3.5,show.legend=F)+
  #geom_text(data=filter(lm.sum.pEmpty,psig=="*"),aes(x=year,y=predict,color=species,label=round(year.trend,digits=2)),size=3)+
  #geom_text(aes(x=2019.75,y=filter(pik.pred,year==max(year))$predict,label=round(pe.change,digits=2)),size=6,hjust=0)+
  scale_fill_viridis_d(option="A")+
  scale_color_viridis_d(option="A")+
  scale_shape_manual(values=c(21,24))+
  scale_size_continuous(range=c(1,6))+
  scale_y_continuous(name="Empty Stomachs",labels=scales::percent_format(),
                     expand=expansion(add=0.035))+
  scale_x_continuous(name="Year",breaks=c(1975,1985,1995,2005,2015),expand=expansion(add=c(0.5,0.5)))+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),
        plot.margin=unit(c(1,0.5,0.25,1),units="cm"),legend.position = "none")+
  guides(fill=guide_legend(override.aes = list(fill="transparent",color="transparent")))

DB<-ggplot()+
  #  geom_point(data=indResponses.good,aes(x=year,y=simpson,fill=Species,shape=season,size=abundance),alpha=0.2)+
  #geom_ribbon(data=db.sp.pred,aes(x=year,ymin=min,ymax=max,fill=species,group=paste(species,season)),alpha=0.3)+
  geom_line(data=db.sp.pred,aes(x=year,y=predict,color=species,linetype=season,group=paste(species,season)),
            size=2,alpha=0.5)+
  geom_ribbon(data=db.pred,aes(x=year,ymin=min,ymax=max),fill="black",alpha=0.7)+
  geom_line(data=db.pred,aes(x=year,y=predict),color="black",size=3.5)+
  #geom_text(data=filter(lm.sum.dBreadth,psig=="*"),aes(x=year,y=predict,color=species,label=round(year.trend,digits=2)),size=4)+
  #geom_text(aes(x=2019.25,y=filter(db.pred,year==max(year))$predict,label=round(db.change,digits=2)),size=6,hjust=0)+
  scale_fill_viridis_d(option="A")+
  scale_color_viridis_d(option="A")+
  scale_shape_manual(values=c(21,24))+
  scale_size_continuous(range=c(1,6))+
  scale_y_continuous("Diet Breadth",expand=expansion(add=0.035))+
  scale_x_continuous(name="Year",breaks=c(1975,1985,1995,2005,2015),expand=expansion(add=c(0.5,0.5)))+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),strip.text=element_blank(),
        plot.margin=unit(c(1,0.5,0.25,1),units="cm"),legend.position="none")

RC<-ggplot()+
  #geom_point(data=indResponses.good,aes(x=year,y=meanC,fill=Species,shape=season,size=abundance),alpha=0.2)+
  #geom_ribbon(data=rc.sp.pred,aes(x=year,ymin=min,ymax=max,fill=species,group=paste(species,season)),alpha=0.3)+
  geom_line(data=rc.sp.pred,aes(x=year,y=predict,color=species,linetype=season,group=paste(species,season)),
            size=2,alpha=0.5)+
  geom_ribbon(data=rc.pred,aes(x=year,ymin=min,ymax=max),fill="black",alpha=0.7)+
  geom_line(data=rc.pred,aes(x=year,y=predict),color="black",size=3.5)+
  #geom_text(data=filter(lm.sum.relC,psig=="*"),aes(x=year,y=predict,color=species,label=round(year.trend,digits=2)),size=4)+
  #geom_text(aes(x=2019.25,y=filter(rc.pred,year==max(year))$predict,label=round(rc.change,digits=2)),size=6,hjust=0)+
  scale_fill_viridis_d(option="A")+
  scale_color_viridis_d(option="A")+
  scale_shape_manual(values=c(21,24))+
  scale_size_continuous(range=c(1,6))+
  scale_y_continuous(name="Relative Consumption (g/g)")+
  scale_x_continuous(name="Year",breaks=c(1975,1985,1995,2005,2015),expand=expansion(add=c(0.5,0.5)),
                     limits=c(1973,2019))+
  theme(strip.text=element_blank(),
        plot.margin=unit(c(1,0.5,1,1),units="cm"),legend.position="none")+
  guides(fill=guide_legend(override.aes = list(fill="transparent",color="transparent")))

png(filename="../Figures/Dissertation/c2_posthoc_metrics_overtime_Allspecies_customWeight3.png",width=750,height=1111)
align_plots(PIK,DB,RC)
dev.off()


pik.sp.pred<-pik.sp.pred%>%complete(species,season)
db.sp.pred<-db.sp.pred%>%complete(species,season)
rc.sp.pred<-rc.sp.pred%>%complete(species,season)
noDataCombos<-filter(pik.sp.pred,is.na(year))%>%select(Species=species,season)
indResponses.good<-bind_rows(noDataCombos,indResponses.good)%>%
  mutate(pdcomnam=toupper(Species))




# Supplemental Comparing to method changes over time ------------------------------------------

#before
sp.ind.pEmpty.preserved<-glm(cbind(nEmpty,nDiets)~Species*season*year,
                             data=filter(indResponses.good,year<1981 | 
                                           (year<=1984 & Species%in%c("Atlantic Cod","Haddock",
                                                                      "Silver Hake","Yellowtail Flounder",
                                                                      "Winter Flounder","Atlantic Herring",
                                                                      "Atlantic Mackerel"))),
                             family=binomial())
pe.preserved <- as.data.frame(lstrends(sp.ind.pEmpty.preserved, c("Species","season"), var="year"))%>%
  mutate(year.trend.preserved=year.trend)
#after
sp.ind.pEmpty.sea<-glm(cbind(nEmpty,nDiets)~Species*season*year,
                       data=filter(indResponses.good, year>=1985 | 
                                     (year>=1981 & Species%notin%c("Atlantic Cod","Haddock",
                                                                   "Silver Hake","Yellowtail Flounder",
                                                                   "Winter Flounder","Atlantic Herring",
                                                                   "Atlantic Mackerel"))),family=binomial())
pe.sea <- as.data.frame(lstrends(sp.ind.pEmpty.sea, c("Species","season"), var="year"))%>%
  mutate(year.trend.sea=year.trend)

#Combining
ba.pe<-left_join(pe.preserved%>%select(Species,season,year.trend.preserved),
                 pe.sea%>%select(Species,season,year.trend.sea))%>%
  left_join(pe.lst%>%select(Species=species,season,year.trend))%>%
  mutate(sampling=ifelse(year.trend.preserved>year.trend.sea,"preserved","at sea"))
table(ba.pe$sampling)

pe.sp.pred.sea<-data.frame(year=filter(indResponses.good, year>=1985 | 
                                         (year>=1981 & Species%notin%c("Atlantic Cod","Haddock",
                                                                       "Silver Hake","Yellowtail Flounder",
                                                                       "Winter Flounder","Atlantic Herring",
                                                                       "Atlantic Mackerel")))$year,
                           season=filter(indResponses.good, year>=1985 | 
                                           (year>=1981 & Species%notin%c("Atlantic Cod","Haddock",
                                                                         "Silver Hake","Yellowtail Flounder",
                                                                         "Winter Flounder","Atlantic Herring",
                                                                         "Atlantic Mackerel")))$season,
                           species=filter(indResponses.good, year>=1985 | 
                                            (year>=1981 & Species%notin%c("Atlantic Cod","Haddock",
                                                                          "Silver Hake","Yellowtail Flounder",
                                                                          "Winter Flounder","Atlantic Herring",
                                                                          "Atlantic Mackerel")))$Species,
                           predict=predict.glm(sp.ind.pEmpty.sea,type="response"),
                           se=predict.glm(sp.ind.pEmpty.sea,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
pe.sp.pred.preserved<-data.frame(year=filter(indResponses.good,year<1981 | 
                                               (year<=1984 & Species%in%c("Atlantic Cod","Haddock",
                                                                          "Silver Hake","Yellowtail Flounder",
                                                                          "Winter Flounder","Atlantic Herring",
                                                                          "Atlantic Mackerel")))$year,
                                 season=filter(indResponses.good,year<1981 | 
                                                 (year<=1984 & Species%in%c("Atlantic Cod","Haddock",
                                                                            "Silver Hake","Yellowtail Flounder",
                                                                            "Winter Flounder","Atlantic Herring",
                                                                            "Atlantic Mackerel")))$season,
                                 species=filter(indResponses.good,year<1981 | 
                                                  (year<=1984 & Species%in%c("Atlantic Cod","Haddock",
                                                                             "Silver Hake","Yellowtail Flounder",
                                                                             "Winter Flounder","Atlantic Herring",
                                                                             "Atlantic Mackerel")))$Species,
                                 predict=predict.glm(sp.ind.pEmpty.preserved,type="response"),
                                 se=predict.glm(sp.ind.pEmpty.preserved,type="response",se.fit=T)$se.fit)%>%
  distinct()%>%
  mutate(min=predict-1.96*se,max=predict+1.96*se)
ggplot()+
  geom_line(data=pik.sp.pred,
            aes(x=year,y=predict,color=species,group=paste(species,season),linetype=season),
            color="black",linewidth=3,alpha=0.5,show.legend=F)+
  geom_line(data=pe.sp.pred.sea,
            aes(x=year,y=predict,group=paste(species,season),linetype=season),
            color="red",linewidth=1.5)+
  geom_line(data=pe.sp.pred.preserved,
            aes(x=year,y=predict,group=paste(species,season),linetype=season),
            color="blue",linewidth=1.5)+
  geom_vline(xintercept=1981,linetype="dashed")+
  geom_vline(xintercept=1985,linetype="dashed")+
  scale_y_continuous(name="Empty Stomachs",labels=scales::percent_format(),
                     expand=expansion(add=0.01))+
  scale_x_continuous(name="Year",breaks=c(1975,1995,2015),labels=c("  1975","1995","2015  "),expand=expansion(add=c(0.5,0.5)))+
  theme(plot.margin=unit(c(1,0.5,0.25,1),units="cm"),legend.position=c(0.7,0.05),legend.direction="horizontal",
        strip.text=element_text(size=15),legend.key.size=unit(0.6,units="in"))+
  facet_wrap(~species)

ba.pe%>%
  pivot_longer(cols=c(year.trend,year.trend.preserved,year.trend.sea))%>%
  filter(name!="year.trend.preserved")%>%
  ggplot()+
  geom_col(aes(x=paste(Species,season),y=value,fill=name),position=position_dodge())+
  scale_fill_discrete(label=c("Total","Only At Sea"))+
  ylab("Slope of Change")+
  theme(axis.text.x=element_blank(),legend.title=element_blank(),axis.title.x=element_blank(),
        legend.position=c(0.66,0.8))

mean(ba.pe$year.trend,na.rm=T)
mean(ba.pe$year.trend.sea,na.rm=T)

sp.responseCats.pe.sea<-pe.sp.pred.sea%>%
  group_by(species,season)%>%
  filter(year==min(year) | year==max(year))%>%
  mutate(year=ifelse(year==min(year),"Min","Max"))%>%
  pivot_wider(id_cols=c("species","season"),names_from="year",values_from="predict")%>%
  mutate(slope=Max*100-100*Min,
         pChange=Max/Min,
         responseCat=case_when(slope>=5 ~"Strong Positive",
                               slope<5 & slope>=1.666~"Weak Positive",
                               slope<=-1.666 & slope>-5~"Weak Negative",
                               slope<=-5~"Strong Negative",
                               TRUE~"Neutral"),
         responseCat2=case_when(pChange>=1.3~"Strong Positive",
                                pChange>=1.1 & pChange<1.3~"Weak Positive",
                                pChange<=0.9 & pChange>0.7~"Weak Negative",
                                pChange<=0.7 ~ "Strong Negative",
                                TRUE~"Neutral"),
         responseCat=factor(responseCat,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")),
         responseCat2=factor(responseCat2,levels=c("Strong Negative","Weak Negative","Neutral","Weak Positive","Strong Positive")))
table(filter(sp.responseCats.pe,!is.na(Species))$responseCat2)
table(sp.responseCats.pe.sea$responseCat2)
nrow(filter(sp.responseCats.pe.sea,grepl("Negative",responseCat2)))/59
nrow(filter(sp.responseCats.pe,grepl("Negative",responseCat2)))/59
nrow(filter(sp.responseCats.pe.sea,responseCat2=="Neutral"))/59

ggplot()+
  geom_bar(data=sp.responseCats.pe.sea,aes(x=1,fill=responseCat2),position=position_stack(),color="black")+
  geom_bar(data=sp.responseCats.pe,aes(x=2,fill=responseCat2),position=position_stack(),color="black")+
  scale_y_continuous(breaks=c(0,0.25*59,0.5*59,0.75*59,59),
                     labels=c("0%","25%","50%","75%","100%"),
                     expand=expansion(add=0),
                     name="Percent of Species-Season")+
  scale_x_discrete(name="At Sea Only        Total     \nPercent Empty Stomachs",
                   expand=expansion(add=0.1))+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T),
                    name="Response\nCategory",labels=gsub(" ","\n",levels(sp.responseCats.pe$responseCat2)))+
  theme(legend.box.margin=margin(0,0,0,-22),legend.key.height=unit(2,"cm"))

#how many changed less than total
seaChange<-left_join(sp.responseCats.pe.sea%>%rename(Species=species),sp.responseCats.pe,by=c("Species","season"))%>%
  select(Species,season,responseCat2.sea=responseCat2.x,responseCat2=responseCat2.y)%>%
  mutate(sea_score=case_when(responseCat2.sea=="Strong Negative"~2,
                             responseCat2.sea=="Weak Negative"~1,
                             responseCat2.sea=="Neutral"~0,
                             responseCat2.sea=="Weak Positive"~-1,
                             responseCat2.sea=="Strong Positive"~-2),
         total_score=case_when(responseCat2=="Strong Negative"~2,
                               responseCat2=="Weak Negative"~1,
                               responseCat2=="Neutral"~0,
                               responseCat2=="Weak Positive"~-1,
                               responseCat2=="Strong Positive"~-2),
         change=total_score-sea_score)
mean(seaChange$change)

#How many are higher the first survey after switching to at sea
atSeaSwitch1<-filter(indResponses.good, (year==1984 | year==1985) & 
                                          Species%in%c("Atlantic Cod","Haddock",
                                                       "Silver Hake","Yellowtail Flounder",
                                                       "Winter Flounder","Atlantic Herring",
                                                       "Atlantic Mackerel"))
atSeaSwitch2<-filter(indResponses.good, (year==1980 | year==1981) & 
                       Species%notin%c("Atlantic Cod","Haddock",
                                    "Silver Hake","Yellowtail Flounder",
                                    "Winter Flounder","Atlantic Herring",
                                    "Atlantic Mackerel"))
atSeaSwitch<-bind_rows(atSeaSwitch1,atSeaSwitch2)%>%
  group_by(Species,year)%>%
  reframe(pEmpty=sum(nEmpty)/sum(nDiets))%>%
  group_by(Species)%>%
  filter(n()>1)%>%
  mutate(year=ifelse(year==1980|year==1984,"before","after"))%>%
  pivot_wider(id_cols=Species,names_from=year,values_from=pEmpty)%>%
  mutate(diff=after-before)
t.test(atSeaSwitch$before,atSeaSwitch$after,paired=T)


#### Individual plots #####


for (i in 1:nrow(species_key_good)) {
  
  empty<-ggplot()+
    geom_point(data=filter(indResponses.good,pdcomnam %in% species_key_good[i,3]),
               aes(x=year,y=pik,fill=Species,shape=season,size=abundance),alpha=0.2)+
    geom_ribbon(data=filter(pik.sp.pred,species==str_to_title(species_key_good[i,3])),
                aes(x=year,ymin=min,ymax=max,fill=species,group=paste(season,species)),alpha=0.5)+
    geom_line(data=filter(pik.sp.pred,species==str_to_title(species_key_good[i,3])),
              aes(x=year,y=predict,color=species,linetype=season,group=paste(season,species)),size=2)+
    geom_ribbon(data=pik.pred,aes(x=year,ymin=min,ymax=max),alpha=0.7,show.legend = F)+
    geom_line(data=pik.pred,aes(x=year,y=predict),size=2,show.legend=F)+
    geom_text_repel(data=filter(lm.sum.pEmpty,species==str_to_title(species_key_good[i,3])),
                    aes(x=ifelse(sign(round(diff,digits=2))>0,2019.75,2019.25),y=predict,color=species,label=round(diff,digits=2)),
                    size=6,hjust=0,direction="y",box.padding=0)+
    scale_fill_manual(values=viridis::viridis(nrow(species_key_good),option="A")[i])+
    scale_color_manual(values=viridis::viridis(nrow(species_key_good),option="A")[i])+
    scale_shape_manual(values=c(21,24))+
    scale_size_continuous(range=c(1,6))+
    scale_y_continuous(name="Empty Stomachs",labels=scales::percent_format(),
                       limits=c(0,1),expand=expansion(add=0.01))+
    scale_x_continuous(name="Year",breaks=c(1975,1985,1995,2005,2015),expand=expansion(add=c(0.5,3.5)))+
    theme(axis.text.x=element_blank(),axis.title.x=element_blank(),plot.title.position = "plot",
          plot.margin=unit(c(1,0.5,0.25,1),units="cm"),legend.position = "none",plot.title=element_text(hjust=0.5))+
    guides(fill=guide_legend(override.aes = list(fill="transparent",color="transparent")))+
    ggtitle(str_to_title(species_key_good[i,3]))
  
  breadth<-ggplot()+
    geom_point(data=filter(indResponses.good,pdcomnam%in%species_key_good[i,3]),
               aes(x=year,y=simpson,fill=Species,shape=season,size=abundance),alpha=0.2)+
    geom_ribbon(data=filter(db.sp.pred,species==str_to_title(species_key_good[i,3])),
                aes(x=year,ymin=min,ymax=max,fill=species,group=paste(species,season)),alpha=0.5)+
    geom_line(data=filter(db.sp.pred,species==str_to_title(species_key_good[i,3])),
              aes(x=year,y=predict,color=species,linetype=season,group=paste(species,season)),size=2)+
    geom_ribbon(data=db.pred,aes(x=year,ymin=min,ymax=max),fill="black",alpha=0.7)+
    geom_line(data=db.pred,aes(x=year,y=predict),color="black",size=2)+
    geom_text_repel(data=filter(lm.sum.dBreadth,species==str_to_title(species_key_good[i,3])),
                    aes(x=ifelse(sign(round(diff,digits=2))>0,2019.75,2019.25),y=predict,color=species,label=round(diff,digits=2)),
                    size=6,hjust=0,direction="y",box.padding=0)+
    scale_fill_manual(values=viridis::viridis(nrow(species_key_good),option="A")[i])+
    scale_color_manual(values=viridis::viridis(nrow(species_key_good),option="A")[i])+
    scale_shape_manual(values=c(21,24))+
    scale_size_continuous(range=c(1,6))+
    scale_y_continuous("Diet Breadth",limits=c(0,1),expand=expansion(add=0.01))+
    scale_x_continuous(name="Year",breaks=c(1975,1985,1995,2005,2015),expand=expansion(add=c(0.5,3.5)))+
    theme(axis.text.x=element_blank(),axis.title.x=element_blank(),strip.text=element_blank(),
          plot.margin=unit(c(1,0.5,0.25,1),units="cm"),legend.position="none")
  
  consump<-ggplot()+
    geom_point(data=filter(indResponses.good,pdcomnam%in%species_key_good[i,3]),
               aes(x=year,y=meanC,fill=Species,shape=season,size=abundance),alpha=0.2)+
    geom_ribbon(data=filter(rc.sp.pred,species==str_to_title(species_key_good[i,3])),
                aes(x=year,ymin=min,ymax=max,fill=species,group=paste(species,season)),alpha=0.5)+
    geom_line(data=filter(rc.sp.pred,species==str_to_title(species_key_good[i,3])),
              aes(x=year,y=predict,color=species,linetype=season,group=paste(species,season)),size=2)+
    geom_ribbon(data=rc.pred,aes(x=year,ymin=min,ymax=max),fill="black",alpha=0.7)+
    geom_line(data=rc.pred,aes(x=year,y=predict),color="black",size=2)+
    geom_text_repel(data=filter(lm.sum.relC,species==str_to_title(species_key_good[i,3])),
                    aes(x=ifelse(sign(round(diff,digits=2))>0,2019.75,2019.25),y=predict,color=species,label=round(diff,digits=2)),
                    size=6,hjust=0,direction="y",box.padding=0)+
    scale_fill_manual(values=viridis::viridis(nrow(species_key_good),option="A")[i])+
    scale_color_manual(values=viridis::viridis(nrow(species_key_good),option="A")[i])+
    scale_shape_manual(values=c(21,24))+
    scale_size_continuous(range=c(1,6))+
    scale_y_log10(name="Relative Consumption (g/g)",breaks=c(0.00001,0.0001,0.001,0.01,0.1,1),
                  labels=c("0.00001","0.0001","0.001","0.01","0.1","1"),
                  minor_breaks=log10_minor_break())+
    scale_x_continuous(name="Year",breaks=c(1975,1985,1995,2005,2015),expand=expansion(add=c(0.5,3.5)))+
    theme(strip.text=element_blank(),
          plot.margin=unit(c(1,0.5,1,1),units="cm"),legend.position="none")+
    guides(fill=guide_legend(override.aes = list(fill="transparent",color="transparent")))
  
  png(filename=paste0("../Figures/2024 Spring/Manuscript Test Figures/Manuscript 8.5/indMetrics_overtime/",species_key_good[i,3],".png"),width=750,height=1111)
  align_plots(empty,breadth,consump)
  dev.off()
  
}

# Total Species Metric Summaries ---------------------------
#Total frequency of empty
sum(indResponses$nEmpty)/sum(indResponses$nDiets)
nrow(uniqueDietTrawls%>%filter(Empty=="Y"))/nrow(uniqueDietTrawls)
pe<-uniqueDietTrawls%>%
  group_by(comname)%>%
  summarise(n=n(),
            e=sum(Empty=="Y"),
            pe=sum(Empty=="Y")/n())
pe.ind<-indResponses.good%>%
  group_by(Species,year)%>%
  reframe(pe=mean(pik))
spMean.pe<-indResponses%>%
  group_by(pdcomnam)%>%
  summarise(mean.pe=mean(nEmpty/nDiets,na.rm=T))
#Total simpson's diversity index
range(indResponses.good$simpson,na.rm=T)
spMean.db<-indResponses%>%
  group_by(pdcomnam)%>%
  summarise(mean.db=mean(simpson,na.rm=T))
median(indResponses.good$simpson,na.rm=T)


#Rel consumption
spMean.rc<-indResponses.good%>%
  group_by(pdcomnam)%>%
  summarise(mean.rc=mean(meanC,na.rm=T))
median(indResponses.good$meanC,na.rm=T)

#Absolute energy
## CSL: trawlEDs_good not found ##
summary(trawlEDs_good$meankJ)

#ED
summary(trawlEDs_good$meanED)

#Rel energy
summary(filter(trawlEDs_good,!is.infinite(meanRE))$meanRE)
filter(guildEDs,!is.infinite(meanRE))%>%
  group_by(guild)%>%
  reframe(re=median(meanRE,na.rm=T))


#Correlation between diet breadth and percentage empty and relative consumption
cor.test(spMean.db$mean.db,spMean.pe$mean.pe,method="pearson")
cor.test(spMean.rc$mean.rc,spMean.pe$mean.pe,method="pearson")

cor.test(indResponses.good$meanC,indResponses.good$pik,method="pearson")


left_join(spMean.db,spMean.pe)%>%
  ggplot()+
  geom_point(aes(mean.db,mean.pe,color=pdcomnam))
left_join(spMean.rc,spMean.pe)%>%
  ggplot()+
  geom_point(aes(mean.rc,mean.pe,color=pdcomnam))


#Total relative consumption
median(indResponses.good$meanC,na.rm=T)
#How often was more than 10% achieved for each species
indResponses.good%>%
  mutate(binge=ifelse(meanC>0.1,"Y","N"))%>%
  group_by(pdcomnam)%>%
  summarise(Nbinge=sum(binge=="Y",na.rm=T))%>%
  arrange(-Nbinge)%>%
  mutate(Species_comnam=str_to_title(pdcomnam))%>%
  left_join(dietCountandModel)%>%
  mutate(pbinge=Nbinge/nDiets)%>%
  arrange(-pbinge)
filter(indResponses.good,pdcomnam=="WINTER FLOUNDER")%>%arrange(-meanC)




# random draws of diet breadth
set.seed(42)
species_simpsons<-filter(indResponses,!is.na(simpson))%>%
  ungroup()%>%
  select(id,Species,simpson)%>%arrange(Species)%>%
  pivot_wider(names_from=Species,values_from=simpson,values_fill=NA)
key_wo_shark<-filter(species_key,pdcomnam!="ATLANTIC SHARPNOSE SHARK")
randomDB_draws<-data.frame(matrix(ncol=nrow(key_wo_shark),nrow=100))
for (j in 1:nrow(key_wo_shark)) {
  i=1
  i_simpsons<-subset(species_simpsons[,j+1],!is.na(species_simpsons[,j+1]))
  while (i<101) {
    randomDB_draws[i,j]<-sapply(sample_n(i_simpsons,100,replace=T),mean)
    i=i+1
  }
}
colnames(randomDB_draws)<-key_wo_shark$pdcomnam
pivot_longer(randomDB_draws,cols=everything(),
             names_to="Species",values_to="Simpson")%>%
  ggplot()+
  geom_boxplot(aes(x=fct_reorder(str_to_title(Species),Simpson),y=Simpson),
               outlier.shape=NA,linewidth=2)+
  geom_jitter(aes(x=fct_reorder(str_to_title(Species),Simpson),y=Simpson),
              width=0.2,alpha=0.8)+
  geom_errorbar(data=spMean.db%>%filter(simpson>0),
                aes(x=Species,ymin=simpson,ymax=simpson),
                linewidth=1,color="firebrick2",alpha=0.5)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.33))+
  xlab("Species")

