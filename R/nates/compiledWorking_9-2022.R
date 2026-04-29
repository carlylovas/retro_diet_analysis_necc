
rm(list=ls())

# Top ---------------------------------------------------------------------

# GL repeat ---------------------------------------------------------------


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
library(dendextend)
library(car)
library(ggpubr)
library(MuMIn)
library(lme4)

#Personalization
theme_set(theme_bw(base_size=25))
options(dplyr.summarise.inform = FALSE)
'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}

seasonPal<-c("steelblue4","goldenrod","deeppink3","sienna2")
seasonPal2<-c("deepskyblue4","yellowgreen","lightgoldenrod","orange3") #BEST
seasonPal3<-c("deepskyblue4","yellowgreen","paleturquoise2","orange3")

#Data loading
load("NF.prey19.RData")
allprey<-NF.prey19
rm(NF.prey19)

#ALL LEVELS OF YEAR AND SEASON
ALLyearseasons<-factor(levels=c(paste(sort(rep(seq(1973,2019),4)),rep(c("Winter","Spring","Summer","Fall"),length(seq(1973,2019))))))

#The names of my species, just to have them nice and handy
species_sci<-unique(str_to_sentence(allprey$pdscinam))
species_com<-unique(str_to_title(allprey$pdcomnam))
#helping me decode them
species_key<-unique(dplyr::select(allprey,pdscinam,pdcomnam))

#These are all the species in the trawls
spec<-read.csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv")%>%
  dplyr::select(LOGGED_SPECIES_NAME,SVSPP)%>%distinct()%>%
  filter(SVSPP>="001")%>%
  mutate(LOGGED_SPECIES_NAME=str_to_sentence(LOGGED_SPECIES_NAME))
spec<-spec[!duplicated(spec$SVSPP),]


#Size classes used by Garrison and Link
sizeClasses<-read.csv("GarrisonLink_predSizeClasses.csv")%>%
  left_join(spec,by=c("species_comnam"="LOGGED_SPECIES_NAME"))%>%
  rename(svspp=SVSPP)
#These are replicated in the sizecats column already present, except that the small and largest cats extend to the smallest and largest individuals
allprey%>%group_by(pdcomnam,sizecat)%>%summarise(min=min(pdlen),max=max(pdlen))


#Species with SDMs (from Allyn et al. 2020 PLOS One paper)
SDMs<-read.csv("COCA2020_ModeledSpeciesTable.csv")
SDMs$SPECIES<-gsub(" $","",SDMs$SPECIES)
table(SDMs$COCA2020)

#Species with SDMs and diet data
filter(SDMs, str_to_sentence(SDMs$SPECIES)%in%species_sci & COCA2020=="Y")


#Hare et al. 2016 expert rankings on species ecology
hare_qualScores<-read_csv("S1Dataset.csv")
colnames(hare_qualScores)<-gsub(" ","\\.",colnames(hare_qualScores))  
hare_qualScores<-hare_qualScores%>%
  mutate(Species=str_to_title(Species),
         Species=gsub("Sand Lances","Northern Sand Lance",Species),
         Species=gsub("Monkfish","Goosefish",Species),
         Species=gsub("Inshore ","",Species),
         totalScore=1*Low+2*Moderate+3*High+4*Very.High,
         meanScore=ifelse(Attribute.Category=="Sensitivity Attribute",totalScore/25,totalScore/20))
hare_directScores<-read_csv("S2Dataset.csv")
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
  left_join(sizeClasses)%>%
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
gl_preycats<-read.csv("GarrisonLink_preyCats.csv")%>%
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


check<-filter(allprey,is.na(GL_pycomnam))[,c("gensci","analsci","collsci","pynam","INnum","INpy")]
sort(table(allprey$GL_pyscinam))
printOut<-allprey%>%
  group_by(gensci,analsci,collsci,pynam,GL_pyscinam)%>%
  summarise(N=n())%>%
  arrange(GL_pyscinam)
#write.csv(printOut,"allprey_glpreycats.v15.csv",row.names = F)

uniqueDiets<-allprey%>%
  dplyr::select(cruise6,station,svspp,pdsex,pdid,pdcomnam,
                pdscinam,dietID,pdlen,pdwgt,sizecat2,pdgutw,pdgutv,
                declat,declon,month,day,year,season,geoarea,id,nDiets)%>%
  distinct()


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

allprey<-mutate(allprey,GL_pyscinam=factor(GL_pyscinam,levels=gl_levels))


# Adding in the trawl data (from GMRI) ------------------------------------

#reading in prey trawls df from "connecting_food+trawls.R" where prey19 is merged with GMRI's clean data for trawls
trawls<-read_csv("../Trawl Data/NMFS Trawls/Complete/NMFS_survdat_gmri_tidy.csv",
                 col_types=c("cccccccnnncnnTnnnnnnnnnnccccn"))

n_distinct(trawls$id)

abundance<-trawls%>%
  ungroup()%>%
  dplyr::select(year=est_year,season,id,pdcomnam=comname,catchsex,abundance,biomass_kg)%>%distinct()%>%
  group_by(year,season,id,pdcomnam)%>%
  summarise(abundance=sum(abundance),
            biomass_kg=sum(biomass_kg))

trawls<-trawls%>%
  left_join(sizeClasses)%>%
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
  summarise(Wk=sum(sizecat_abundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(sizecat_abundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(sizecat_abundance*pik)/sumAbun, #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizecat_abundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizecat_abundance,qikw,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,sizecat_abundance,pik,Fk))%>%distinct()


check<-trawlDietSum%>%group_by(period,pdscinam,sizecat2)%>%summarise(N=sum(Wk,na.rm=T))
sum(check$N<0&check$N>0) #When this is 0, then you have all your means correct because they sum to 1 (or they're all empties)

#Saving these, so they can be better used within other analyses
#write.csv(trawlDiets,"trawl_speciesClusterCompositions.csv",row.names=F)
#write.csv(trawlDietSum,"geoareayearseason_speciesClusterCompositions.csv",row.names=F)








# Data Summaries ----------------------------------------------------------

#Number of diets for each predator
dietCounts<-as.data.frame(sort(table(uniqueDiets$pdcomnam)))%>%
  mutate(species=str_to_title(Var1))%>%
  left_join(filter(SDMs, str_to_sentence(SDMs$SPECIES)%in%species_sci & COCA2020=="Y")%>%
              mutate(species=str_to_title(COMNAME)))%>%
  ungroup()%>%
  mutate(species=fct_reorder(species,Freq),
         predNum=order(species),
         predLab.1=ifelse(predNum%%2==1,str_to_sentence(as.character(species)),NA),
         predLab.2=ifelse(predNum%%2==0,str_to_sentence(as.character(species)),NA),
         color=ifelse(is.na(COCA2020),"grey50","pink"))

ggplot(data=dietCounts)+
  geom_col(aes(x=predNum,y=Freq,fill=COCA2020),color="black")+
  geom_text(aes(x=predNum,Freq+1500,label=Freq,color=COCA2020),angle=90,hjust=0.4,show.legend=F)+
  scale_x_continuous(name="Predator Species",expand=expansion(mult=c(0.025,0.025)),
                     labels=filter(dietCounts,!is.na(predLab.1))$predLab.1,
                     breaks=seq(1,nrow(dietCounts),by=2),
                     sec.axis=dup_axis(name="",
                                       labels=filter(dietCounts,!is.na(predLab.2))$predLab.2,
                                       breaks=seq(2,nrow(dietCounts),by=2)))+
  scale_y_continuous(expand=expansion(add=c(200,2500)),name="Number of Diets")+
  scale_fill_manual(values=c("lightpink2","grey50"),name="SDM Fit",labels=c("Yes","No"))+
  theme(legend.position = c(0.2,0.8),legend.background=element_rect(color="black"),
        axis.text.x.bottom=element_text(angle=20,hjust=0.8,vjust=1,size=18),
        axis.text.x.top=element_text(angle=20,hjust=0.2,vjust=0,size=18),
        plot.margin=margin(t=5,r=25,b=5,l=75))

#and combining with who has a model (from later, so will not work if starting from scratch)
dietCountandModel<-left_join(dietCounts%>%ungroup()%>%mutate(Species_comnam=species),
                             left_join(species_key,indResponses.good%>%ungroup()%>%mutate(inModel="Y")%>%
                               dplyr::select(pdcomnam,pdscinam,inModel)%>%distinct())%>%
                               mutate(Species_comnam=str_to_title(pdcomnam),Species_scinam=str_to_sentence(pdscinam)))%>%
  dplyr::select(Species_comnam,Species_scinam,nDiets=Freq,COCA2020,NECC=inModel)%>%
  mutate(COCA2020=ifelse(is.na(COCA2020),"N","Y"),
         NECC=ifelse(is.na(NECC),"N","Y"))%>%arrange(Species_comnam)
#write.csv(dietCountandModel,"speciesList_inNECC.csv",row.names=F)
#Hare et al. climate scores
hare_allScores%>%
  filter(Species %in% str_to_title(allprey$pdcomnam))%>%
  ggplot()+
  geom_col(aes(Species,meanScore,fill=Attribute.Category,group=paste(Attribute.Category,Attribute)),color="black",position="dodge")+
  scale_fill_viridis_d()+
  theme(axis.text.x=element_text(angle=30,vjust=1.05,hjust=1))

hare_allScores%>%
  filter(Species %in% str_to_title(allprey$pdcomnam))%>%
  group_by(Species,Attribute.Category)%>%
  summarise(meanScore=mean(meanScore))%>%
  ggplot()+
  geom_col(aes(Species,meanScore,fill=Attribute.Category),color="black",position="dodge")+
  scale_fill_viridis_d()+
  theme(axis.text.x=element_text(angle=30,vjust=1.05,hjust=1))


#Tables for exposure scores, sensitivity, and climate direction
hare_overallScores%>%
  filter(Species %in% str_to_title(allprey$pdcomnam))%>%
  mutate(scoreClass=factor(scoreClass,levels=c("Low","Moderate","High","Very High")))%>%
  pivot_wider(id_cols=Species,names_from=Attribute.Category,values_from=scoreClass)%>%
  dplyr::select(-c(`Climate direction`))%>%
  ggplot()+
  geom_text(aes(x=`Exposure Factor`,y=`Sensitivity Attribute`,label=Species))

table(hare_overallScores%>%ungroup()%>%
        filter(Species %in% str_to_title(allprey$pdcomnam) & Attribute.Category!="Climate direction")%>%
        mutate(scoreClass=factor(scoreClass,levels=c("Low","Moderate","High","Very High")))%>%dplyr::select(Attribute.Category,scoreClass))

table(hare_changePotential%>%
        filter(Species %in% str_to_title(allprey$pdcomnam))%>%dplyr::select(scoreClass))


#Are the Hare variables dependent or related to each other
cVulcDir<-table(hare_everything$Climate.Direction,hare_everything$Climate.Vulnerability)
chisq.test(cVulcDir)


cVulcPot<-table(hare_everything$Climate.Vulnerability,hare_everything$Change.Potential)
chisq.test(cVulcPot)


cDircPot<-table(hare_everything$Climate.Direction,hare_everything$Change.Potential)
chisq.test(cDircPot)

plotSpacing<-c(mean(c(-1.233,0.1)),mean(c(0.1,1.433)),mean(c(1.433,2.766)),mean(c(2.766,4.1)))
hare_everything_FHD_plot<-hare_everything_FHD%>%
  arrange(Climate.Vulnerability,desc(Change.Potential),Climate.Direction,Species)%>%
  mutate(Change.Potential.num=case_when(Change.Potential=="Very High"~4,
                                        Change.Potential=="High"~2.667,
                                        Change.Potential=="Moderate"~1.333,
                                        Change.Potential=="Low"~0),
         Climate.Vulnerability.num=case_when(Climate.Vulnerability=="Very High"~plotSpacing[4],
                                             Climate.Vulnerability=="High"~plotSpacing[3],
                                             Climate.Vulnerability=="Moderate"~plotSpacing[2],
                                             Climate.Vulnerability=="Low"~plotSpacing[1]))
for (i in 2:nrow(hare_everything_FHD_plot)) {
  hare_everything_FHD_plot$Change.Potential.num[i]<-ifelse(hare_everything_FHD_plot$Change.Potential[i]==hare_everything_FHD_plot$Change.Potential[i-1] &
                                                             hare_everything_FHD_plot$Climate.Vulnerability[i]==hare_everything_FHD_plot$Climate.Vulnerability[i-1],
                                                           hare_everything_FHD_plot$Change.Potential.num[i-1]-0.1,
                                                           hare_everything_FHD_plot$Change.Potential.num[i])
}


hare_everything_plot<-hare_everything%>%
  arrange(Climate.Vulnerability,desc(Change.Potential),Climate.Direction,Species)%>%
  mutate(Change.Potential.num=case_when(Change.Potential=="Very High"~4,
                                        Change.Potential=="High"~2.667,
                                        Change.Potential=="Moderate"~1.333,
                                        Change.Potential=="Low"~0),
         Climate.Vulnerability.num=case_when(Climate.Vulnerability=="Very High"~plotSpacing[4],
                                             Climate.Vulnerability=="High"~plotSpacing[3],
                                             Climate.Vulnerability=="Moderate"~plotSpacing[2],
                                             Climate.Vulnerability=="Low"~plotSpacing[1]))
for (i in 2:nrow(hare_everything_plot)) {
  hare_everything_plot$Change.Potential.num[i]<-ifelse(hare_everything_plot$Change.Potential[i]==hare_everything_plot$Change.Potential[i-1] &
                                                             hare_everything_plot$Climate.Vulnerability[i]==hare_everything_plot$Climate.Vulnerability[i-1],
                                                           hare_everything_plot$Change.Potential.num[i-1]-0.1,
                                                           hare_everything_plot$Change.Potential.num[i])
}



library(png)
allFish.files<-list.files(path="../Figures/Silhouettes",pattern="silhouette*")

silhouettePNGs<-function(fileName) {
  rasterGrob(readPNG(paste0("../Figures/Silhouettes/",fileName)),interpolate=T)
}
allFish.pngs<-lapply(allFish.files,silhouettePNGs)


ggplot(hare_everything_FHD_plot)+
  annotation_custom(allFish.pngs[[1]], xmin=1.05,xmax=1.38,ymin=2.56,ymax=2.76)+    #Acadian Redfish
  annotation_custom(allFish.pngs[[2]], xmin=3.66,xmax=4,ymin=-0.1,ymax=0.1)+        #Alewife
  annotation_custom(allFish.pngs[[3]], xmin=-0.3,xmax=0.03,ymin=2.56,ymax=2.76)+    #American Plaice
  annotation_custom(allFish.pngs[[4]], xmin=0.16,xmax=0.5,ymin=2.46,ymax=2.66)+     #Atlantic Cod
  annotation_custom(allFish.pngs[[5]], xmin=0.16,xmax=0.5,ymin=1.66,ymax=1.86)+     #Atlantic Croaker
  annotation_custom(allFish.pngs[[6]], xmin=-1.2,xmax=-0.866,ymin=2.46,ymax=2.66)+  #Atlantic Herring
  annotation_custom(allFish.pngs[[7]], xmin=1.05,xmax=1.38,ymin=2.36,ymax=2.56)+    #Atlantic Mackerel
  #annotation_custom(allFish.pngs[[8]], xmin=0,xmax=1,ymin=1,ymax=1.1)+             #Atlantic Sharpnose Shark
  annotation_custom(allFish.pngs[[9]], xmin=2.33,xmax=2.66,ymin=2.36,ymax=2.56)+    #Black Sea Bass
  annotation_custom(allFish.pngs[[10]],xmin=-0.3,xmax=0.03,ymin=1.56,ymax=1.76)+    #Bluefish
  annotation_custom(allFish.pngs[[11]],xmin=-0.3,xmax=0.03,ymin=3.7,ymax=3.9)+      #Butterfish
  #annotation_custom(allFish.pngs[[12]],xmin=0,xmax=1,ymin=1,ymax=1.1)+             #Fawn Cusk-eel
  #annotation_custom(allFish.pngs[[13]],xmin=0,xmax=1,ymin=1,ymax=1.1)+             #Fourspot Flounder
  #annotation_custom(allFish.pngs[[14]],xmin=0,xmax=1,ymin=1,ymax=1.1)+             #Gulfstream Flounder
  annotation_custom(allFish.pngs[[15]],xmin=-0.3,xmax=0.03,ymin=2.36,ymax=2.56)+    #Haddock
  annotation_custom(allFish.pngs[[16]],xmin=-1.2,xmax=-0.866,ymin=2.26,ymax=2.46)+  #Little Skate
  annotation_custom(allFish.pngs[[17]],xmin=-1.2,xmax=-0.866,ymin=1.46,ymax=1.66)+  #Longfin Squid
  #annotation_custom(allFish.pngs[[18]],xmin=0,xmax=1,ymin=1,ymax=1.1)+             #Longhorn Sculpin
  annotation_custom(allFish.pngs[[19]],xmin=-1.2,xmax=-0.866,ymin=1.86,ymax=2.06)+  #Monkfish
  annotation_custom(allFish.pngs[[20]],xmin=0.16,xmax=0.5,ymin=2.26,ymax=2.46)+     #Northern Sand Lance
  annotation_custom(allFish.pngs[[21]],xmin=2.33,xmax=2.66,ymin=1.23,ymax=1.43)+    #Ocean Pout
  annotation_custom(allFish.pngs[[22]],xmin=1.05,xmax=1.38,ymin=2.16,ymax=2.36)+    #Pollock
  annotation_custom(allFish.pngs[[23]],xmin=-0.3,xmax=0.03,ymin=1.76,ymax=1.96)+    #Red Hake
  annotation_custom(allFish.pngs[[24]],xmin=1.05,xmax=1.38,ymin=1.56,ymax=1.76)+    #Scup
  #annotation_custom(allFish.pngs[[25]],xmin=0,xmax=1,ymin=1,ymax=1.1)+             #Sea Raven
  annotation_custom(allFish.pngs[[26]],xmin=-1.2,xmax=-0.866,ymin=3.6,ymax=3.8)+    #Shortfin Squid
  annotation_custom(allFish.pngs[[27]],xmin=-0.3,xmax=0.03,ymin=2.16,ymax=2.36)+    #Silver Hake
  annotation_custom(allFish.pngs[[28]],xmin=-0.3,xmax=0.03,ymin=3.9,ymax=4.1)+      #Smooth Dogfish
  annotation_custom(allFish.pngs[[29]],xmin=0.16,xmax=0.5,ymin=2.06,ymax=2.26)+     #Smooth Skate
  annotation_custom(allFish.pngs[[30]],xmin=-1.2,xmax=-0.866,ymin=3.8,ymax=4)+      #Spiny Dogfish
  #annotation_custom(allFish.pngs[[31]],xmin=0,xmax=1,ymin=1,ymax=1.1)+             #Spotted Hake
  annotation_custom(allFish.pngs[[32]],xmin=0.16,xmax=0.5,ymin=1.86,ymax=2.06)+     #Summer Flounder
  annotation_custom(allFish.pngs[[33]],xmin=3.66,xmax=4,ymin=1.2,ymax=1.4)+         #Tautog
  annotation_custom(allFish.pngs[[34]],xmin=2.33,xmax=2.66,ymin=2.56,ymax=2.76)+    #Thorny Skate
  annotation_custom(allFish.pngs[[35]],xmin=1.05,xmax=1.38,ymin=1.76,ymax=1.96)+    #Weakfish
  annotation_custom(allFish.pngs[[36]],xmin=1.05,xmax=1.38,ymin=1.96,ymax=2.16)+    #White Hake
  annotation_custom(allFish.pngs[[37]],xmin=-1.2,xmax=-0.866,ymin=1.66,ymax=1.86)+  #Windowpane
  annotation_custom(allFish.pngs[[38]],xmin=3.66,xmax=4,ymin=2.56,ymax=2.76)+       #Winter Flounder
  annotation_custom(allFish.pngs[[39]],xmin=-1.2,xmax=-0.866,ymin=2.06,ymax=2.26)+  #Winter Skate
  annotation_custom(allFish.pngs[[40]],xmin=1.53,xmax=1.86,ymin=2.46,ymax=2.66)+    #Witch Flounder
  annotation_custom(allFish.pngs[[41]],xmin=-0.3,xmax=0.03,ymin=1.96,ymax=2.16)+    #Yellowtail Flounder
  geom_text(aes(x=Climate.Vulnerability.num,y=Change.Potential.num,color=Climate.Direction,label=gsub("Northern ","",Species)),
            size=5,fontface="bold",key_glyph="point")+
  scale_y_continuous(minor_breaks=c(-1.233,0.1,1.433,2.766,4.1),limits=c(-1.233,4.1),expand=expansion(0),
                     breaks=c(plotSpacing),labels=c("Low","Moderate","High","Very High"),name="Change Potential")+
  scale_x_continuous(minor_breaks=c(-1.233,0.1,1.433,2.766,4.1),limits=c(-1.233,4.1),expand=expansion(0),
                     breaks=c(plotSpacing),labels=c("Low","Moderate","High","Very High"),name="Climate Vulnerability")+
  scale_color_manual(values=c(colorspace::diverge_hcl(n=3,palette = "Berlin",rev=T,c=100,l=60)),
                     name="Climate Direction")+
  guides(color=guide_legend(override.aes = list(size=8)))+
  theme(panel.grid.minor = element_line(size=1.5,color="black"),axis.ticks=element_blank(),
        panel.background=element_rect(color="black",size=3),
        panel.grid.major=element_blank(),legend.position=c(plotSpacing[4]/4,plotSpacing[4]/4))





# Dietary Overlap (Schoener) ----------------------------------------------

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

#Another version where the zeroes can also be shuffled into the mix
resampleAllPrey<-function(props) {
  newProps<-matrix(nrow=nrow(props),ncol=ncol(props))
  for (i in 1:nrow(props)){
    eat<-props[i,which(props[i,]>=0 & colnames(props)!="species")]
    colnames(eat)=sample(colnames(eat))
    teat<-cbind(props[i,1],eat)
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


# Using cluster sample means ----------------------------------------------


clusterSD <- function(n,tM,M,pq,FW) { #Function to calculate the SD for cluster means, from Buckel et al. 1999
  var<- (1/(n*(tM/n)^2))*((sum(M^2*(pq-FW)^2))/n-1)
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



# FOO on all FHD species --------------------------------------------------

n_distinct(filter(preyTrawls,is.na(decdeg_beglat))$id)/n_distinct(preyTrawls$id)
#print<-dplyr::select(preyTrawls,-c(pdscinam,pdcomnam,gl_prey,comname))
#write.csv(print,file="SUMMARIZED/foodHabits_trawlData.csv",row.names = F)
save(print,file="SUMMARIZED/foodHabits_trawlData.RData")

colnames(print)

#The skates have had some genus changes since the data started, will just change the ones that are in FHD
preyTrawls<-preyTrawls%>%
  mutate(collsci=case_when(collsci=="RAJA ERINACEA"~"LEUCORAJA ERINACEA",
                           collsci=="RAJA OCELLATA"~"LEUCORAJA OCELLATA",
                           collsci=="RAJA RADIATA" ~"AMBLYRAJA RADIATA",
                           collsci=="RAJA SENTA"   ~"MALACORAJA SENTA",
                           TRUE~collsci))

sum(str_to_sentence(unique(preyTrawls$collsci)) %in% species_sci)
subset(species_sci, species_sci %notin% str_to_sentence(preyTrawls$collsci))
#These 6 species were never IDed in the diets of any FHD predators (offshore anyway)

nrow(filter(preyTrawls,str_to_sentence(collsci) %in% species_sci)) #How many times are FHD species found as prey
FHDprey<-preyTrawls%>%
  filter(str_to_sentence(collsci) %in% species_sci)%>%
  group_by(collsci)%>%
  summarise(prey_N=n(),
            prey_prop=prey_N/nrow(preyTrawls),
            prey_perc=prey_prop*100)%>%
  bind_rows(data.frame(collsci=subset(species_sci,species_sci %notin% str_to_sentence(preyTrawls$collsci)),
                       prey_N=0,prey_prop=0,prey_perc=0))%>%
  mutate(collsci=fct_reorder(as.factor(str_to_sentence(collsci)),prey_N))

ggplot(FHDprey)+
  geom_col(aes(collsci,prey_prop),fill="seagreen4",color="black")+
  geom_text(aes(collsci,prey_prop+0.0006,label=paste0("N=",prey_N)),angle=90,size=5)+
  scale_x_discrete(name="Prey Species")+
  scale_y_continuous(name="Percent of Diets",labels=scales::percent_format(),expand=expansion(add=c(0,0.0006)))+
  theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))

test<-filter(preyTrawls,is.na(species_scinam))

#How often are each of the FHD preds eating other FHD preds
FHDpred<-preyTrawls%>%
  mutate(FHDpred=ifelse(str_to_sentence(collsci) %in% species_sci, "pred","no"))%>%
  group_by(pdscinam,pdcomnam,dietID)%>%
  summarise(FHDpred=ifelse(any(FHDpred=="pred"),"pred","no"))%>%
  group_by(pdscinam,pdcomnam)%>%
  summarise(pred_N=sum(FHDpred=="pred"),
            total_N=n_distinct(dietID),
            pred_prop=pred_N/total_N,
            pred_perc=pred_prop*100)%>%
  ungroup()%>%
  mutate(pdscinam=fct_reorder(as.factor(str_to_sentence(pdscinam)),pred_perc))

sum(FHDpred$pred_N) #Not quite as many as times when FHD species are prey, so some must have eaten two different FHD species in same "meal"

ggplot(FHDpred)+
  geom_col(aes(pdscinam,pred_prop),fill="seagreen4",color="black")+
  geom_text(aes(pdscinam,pred_prop+0.02,label=paste0("N=",pred_N)),angle=90,size=5)+
  scale_x_discrete(name="Predator Species")+
  scale_y_continuous(name="Percent of Diets",labels=scales::percent_format(),expand=expansion(add=c(0,0.02)))+
  theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))



# Garrison-Link Repeat -----------------------------------------------------------

#### Original Time ####

#The first portion, reflecting the same time period as Garrison and Link (1973-1997)
#There is one little skate over 60cm, meaning it is "large" during this time period. Remove it for ease
#Creating a proportion matrix for the prey items for my predators
propLong_gl<-preyTrawls%>%
  filter(year<=1997
         & !is.na(GL_pyscinam))%>%
  mutate(species=str_to_title(pdcomnam))%>%
  group_by(species,sizecat2)%>%
  mutate(nDiets=n_distinct(dietID),
         total_w=sum(pyamtw),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()


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
colnames(props1_gl)<-c("species1","sizecat1","GL_pyscinam","prop_w1")
props1_gl<-complete(props1_gl,species1,sizecat1,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species1,sizecat1)%>%mutate(t=sum(prop_w1,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props1_gl[is.na(props1_gl)]<-0
props2_gl<-propLong_gl[,-max(ncol(propLong_gl))]
colnames(props2_gl)<-c("species2","sizecat2","GL_pyscinam","prop_w2")
props2_gl<-complete(props2_gl,species2,sizecat2,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species2,sizecat2)%>%mutate(t=sum(prop_w2,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props2_gl[is.na(props2_gl)]<-0

overlap_schoener_gl<-full_join(props1_gl,props2_gl)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat_gl<-pivot_wider(overlap_schoener_gl,id_cols=c(species1,sizecat1),names_from=c(species2,sizecat2),values_from = s_do)
overlap_mat_gl<-as.matrix(overlap_mat_gl[,3:ncol(overlap_mat_gl)])
#rownames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#colnames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#overlap_mat[overlap_mat==1]<-NA
for (i in 1:nrow(overlap_mat_gl)) {
  for (j in 1:ncol(overlap_mat_gl)) {
    overlap_mat_gl[i,j]<-ifelse(j>i,NA,overlap_mat_gl[i,j])
  }
}





#### Updated Time ####


#The latter ~half, since 1997 when Garrison and Link had data through
#Creating a proportion matrix for the prey items for my predators
propLong_new<-preyTrawls%>%
  filter(year>1997
         & !is.na(GL_pyscinam))%>%
  mutate(species=str_to_title(pdcomnam))%>% #Maybe want common names later, but the myspecies is scientific so speed
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()


#Following through to get the "real" overlap matrix
props1_new<-propLong_new[,-max(ncol(propLong_new))] #Cut off the volume prop for cleanliness
colnames(props1_new)<-c("species1","sizecat1","GL_pyscinam","prop_w1")
props1_new<-complete(props1_new,species1,sizecat1,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species1,sizecat1)%>%mutate(t=sum(prop_w1,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props1_new[is.na(props1_new)]<-0
props2_new<-propLong_new[,-max(ncol(propLong_new))]
colnames(props2_new)<-c("species2","sizecat2","GL_pyscinam","prop_w2")
props2_new<-complete(props2_new,species2,sizecat2,GL_pyscinam)%>%
  filter(GL_pyscinam!="Empty")%>%
  group_by(species2,sizecat2)%>%mutate(t=sum(prop_w2,na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
props2_new[is.na(props2_new)]<-0

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







#### Comparing the two time periods ####

#Only want to use those species-sizecats that are in both time periods
GLclusts<-read.csv("GarrisonLink_clusters.csv")%>%
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
  summarise(prop_w=sum(pyamtw)/total_w,
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
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_sizecats_inBOTH<-inner_join(dplyr::select(species_sizecats_inGLtime,species,sizecat2)%>%
                                     distinct(),
                                   dplyr::select(species_sizecats_inNEWtime,species,sizecat2)%>%
                                     distinct())


###### Cropping down the Original time, and clustering again #####
shared_props_gl<-species_sizecats_inGLtime%>%
  filter(paste(species,sizecat2) %in% paste(species_sizecats_inBOTH$species,species_sizecats_inBOTH$sizecat2))

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

overlap_shared_gl<-full_join(props1_gl,props2_gl)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
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

reps=1
nreps=250
bootDiffs_gl<-numeric()
while (reps<=nreps) {
  bootDiffs_gl<-c(bootDiffs_gl,resample(propMat_gl)) #resampling across all the zeroes too is a lower significance
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




#To save them out, there is some difficulty with the export losing the colorbar
png(file = "../Figures/2022 Fall/guilds/GLtime_shared_dendrogram.png",   # The directory you want to save the file in
    width = 1500, # The width of the plot in inches
    height = 1500) # The height of the plot in inches
par(mar=c(5,1,2,12))

dend %>% 
  set("branches_lwd", 4) %>%
  # Custom labels
  set("labels_cex", 1) %>%
  #set("labels_col", value = viridis::viridis(14,end=0.8),h = sigGuild_gl) %>%
  #set("branches_k_color", value = viridis::viridis(14,end=0.8), h = sigGuild_gl) %>%
  plot(horiz=TRUE,main="                               1973-1997 Dendrogram",axes=F,xlab="Similarity")
axis(side=1,at=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(v=1-sigGuild_gl,lty=2)
#rect.dendrogram(dend, h=0.5, which=c(1:6),border="transparent",
#                lty = 5, lwd = 0, col=rgb(0.2,0.2,0.2,0.15),horiz=T) #This is a way to add boxes indicating the guilds
legend(0.94,92,legend=ifelse(is.na(unique(myGLclust_dend$guildName)),"Not Included",unique(myGLclust_dend$guildName)),
       fill=unique(myGLclust_dend$guildCol),title="Garrison and Link Feeding Guild",cex=1.7)

colored_bars(colors = myGLclust_dend$guildCol, dend = dend,sort_by_labels_order = F, rowLabels = NA, horiz=T)

# Step 3: Run dev.off() to create the file!
dev.off()




##### Cropping down the new time and clustering again #####
shared_props_new<-species_sizecats_inNEWtime%>%
  filter(paste(species,sizecat2) %in% paste(species_sizecats_inBOTH$species,species_sizecats_inBOTH$sizecat2))

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
propMat_new<-propLong_new%>%
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

#To save them out, there is some difficulty with the export losing the colorbar
png(file = "../Figures/2022 Fall/guilds/newTime_shared_dendrogram.png",   # The directory you want to save the file in
    width = 1500, # The width of the plot in inches
    height = 1500) # The height of the plot in inches
par(mar=c(5,1,2,12))

dend_new %>% 
  set("branches_lwd", 4) %>%
  # Custom labels
  set("labels_cex", 1) %>%
  #set("labels_col", value = viridis::viridis(12,end=0.8),h = sigGuild_new) %>%
  #set("branches_k_color", value = viridis::viridis(12,end=0.8), h = sigGuild_new) %>%
  plot(horiz=T,main="                               1998-2019 Dendrogram",axes=F,xlab="Similarity")
axis(side=1,at=c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0),labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(v=1-sigGuild_new,lty=2)
#rect.dendrogram(dend_new, h=0.5, which=c(1:4),border="transparent",
#                lty = 5, lwd = 0, col=rgb(0.2,0.2,0.2,0.15),horiz=T)
#legend(0.9,90,legend=ifelse(is.na(unique(myGLclust_new_dend$guildName)),"Not Included",unique(myGLclust_new_dend$guildName)),
#       fill=unique(myGLclust_new_dend$guildCol),title="Garrison and Link Feeding Guild",cex=1.9)

colored_bars(colors = myGLclust_new_dend$guildCol, dend = dend_new,sort_by_labels_order = F,
             rowLabels = NA,horiz=T)

# Step 3: Run dev.off() to create the file!
dev.off()


#Simply how much overlap is there in the early vs the late period--t-test of which has a higher average overlap (2019)
  #This suggests that the diets in 2019 are more similar, like they've COMPRESSED
t.test(overlap_mat_shared_gl[which(overlap_mat_shared_gl<1)],overlap_mat_shared_new[which(overlap_mat_shared_new<1)])


#Create a dendlist, which has both the timelines together
dend_compare<-dendlist(dend,dend_new)
par(mar=c(5,1,2,12))
dend_diff(dend,dend_new)



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



sharedClusts%>%
  mutate(group=paste(species,sizecat))%>%
  dplyr::select(group,myCluster_1997,myCluster_2019)%>%
  pivot_longer(cols=c("myCluster_1997","myCluster_2019"),names_to="year",values_to="cluster")%>%
  mutate(year=as.numeric(gsub("myCluster_","",year)))%>%
  ggplot()+
  geom_point(aes(year,cluster,color=group),show.legend = F)+
  geom_line(aes(year,cluster,color=group),show.legend = F)+
  scale_color_viridis_d()

##### Alluvial Plot showing changes in Cluster groups #####
# install.packages("ggalluvial")
library(ggalluvial)

sharedClusts%>%
  group_by(guildclust_1997,myCluster_2019)%>%
  summarise(freq=n())%>%
  mutate(myCluster_2019=as.character(myCluster_2019))%>%
  ggplot(aes(axis1 = guildclust_1997, axis2 = myCluster_2019, y = freq)) +
  geom_alluvium(aes(fill = myCluster_2019))+
  geom_stratum()+
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_continuous(breaks=c(1,2),labels = c("1973-1997", "1998-2019"),
                   expand = c(0.15, 0.05))+
  theme_void()+
  theme(legend.position = "none",axis.text.x=element_text(size=20,vjust=8),
        plot.title=element_text(size=30,hjust=0.5))+
  ggtitle(label="Guild Assignments in Each Time Period")

sharedClusts%>%
  group_by(nameClust_1997,nameClust_2019)%>%
  summarise(freq=n())%>%
  mutate(nameClust_2019=as.character(nameClust_2019))%>%
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



#Some simple summaries about the guilds
guildSummaries_gl<-sharedClusts%>%
  mutate(species_sizecat=paste(species,sizecat,sep="_"))%>%
  left_join(uniqueDiets%>%filter(year<=1997)%>%
              group_by(species_sizecat=paste(str_to_title(pdcomnam),sizecat2,sep="_"))%>%
              summarise(N=n()))%>%
  group_by(guildclust_1997)%>%
  summarise(N_spsize=n(),
            N_species=n_distinct(species),
            N_diets=sum(N,na.rm=T))
guildSummaries_new<-sharedClusts%>%
  mutate(species_sizecat=paste(species,sizecat,sep="_"))%>%
  left_join(uniqueDiets%>%filter(year>1997)%>%
              group_by(species_sizecat=paste(str_to_title(pdcomnam),sizecat2,sep="_"))%>%
              summarise(N=n()))%>%
  group_by(myCluster_2019)%>%
  summarise(N_spsize=n(),
            N_species=n_distinct(species),
            N_diets=sum(N,na.rm=T))


##### Switching Guilds #####
switchingClusts<-sharedClusts%>%
  group_by(guildclust_1997)%>%
  mutate(N_total=n())%>%
  group_by(guildclust_1997,myCluster_2019)%>%
  mutate(N=n(),
         prop=N/N_total,
         sizecat=factor(sizecat,levels=c("S","M","L","XL")))%>%
  dplyr::select(species,sizecat,guildclust_1997,myCluster_2019,N,prop)
switchers<-filter(switchingClusts,prop<0.5)

##### Categorizing Guilds and Switches #####
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

#Keep this to show that the FOURSPOT FLOUNDER L that were their own cluster in <1998 didn't really change (still squids)
fourspotL<-filter(preyTrawls,pdcomnam=="FOURSPOT FLOUNDER" & sizecat2=="L" & year>1997)%>%
  mutate(d=n_distinct(dietID))%>%group_by(GL_pyscinam)%>%
  summarise(N=n(),M=sum(pyamtw,na.rm=T),d=d)%>%distinct()
#Same happening with the American Plaice, they're still really eating a lot of Brittle Stars
amplaiceML<-filter(preyTrawls,pdcomnam=="AMERICAN PLAICE" & sizecat2%in%c("L","M") & year>1997)%>%
  mutate(d=n_distinct(dietID))%>%group_by(GL_pyscinam)%>%
  summarise(N=n(),M=sum(pyamtw,na.rm=T),d=d)%>%distinct()

#Who switched by names
nameSwitchers<-sharedClusts%>%
  filter(nameClust_1997!=nameClust_2019)%>%
  dplyr::select(species,sizecat,guildclust_1997,nameClust_1997,myCluster_2019,nameClust_2019)

trueSwitchers<-inner_join(switchers,nameSwitchers)
trueSwitchers_dietSum_1997<-left_join(trueSwitchers%>%mutate(species=paste(species,sizecat,sep="_"))%>%ungroup()%>%
                                   dplyr::select(species,nameClust=nameClust_1997),
                                 propMat_gl)%>%
  pivot_longer(cols=matches("^[A-Z]",ignore.case=F),names_to="GL_pyscinam",values_to="prop")%>%
  arrange(species,nameClust,desc(prop))%>%
  group_by(species,nameClust)%>%
  mutate(cumProp=cumsum(prop),
         nItems_to90=sum(cumProp<0.90)+1,
         GL_pyscinam_labels=str_pad(GL_pyscinam,max(str_length(preyTrawls$GL_pyscinam),na.rm=T),side="left"))%>%
  slice(1:nItems_to90)%>%
  mutate(time="<1998")

trueSwitchers_dietSum_2019<-left_join(trueSwitchers%>%mutate(species=paste(species,sizecat,sep="_"))%>%ungroup()%>%
                                        dplyr::select(species,nameClust=nameClust_2019),
                                      propMat_new)%>%
  pivot_longer(cols=matches("^[A-Z]",ignore.case=F),names_to="GL_pyscinam",values_to="prop")%>%
  arrange(species,nameClust,desc(prop))%>%
  group_by(species,nameClust)%>%
  mutate(cumProp=cumsum(prop),
         nItems_to90=sum(cumProp<0.90)+1,
         GL_pyscinam_labels=str_pad(GL_pyscinam,max(str_length(preyTrawls$GL_pyscinam),na.rm=T),side="left"))%>%
  slice(1:nItems_to90)%>%
  mutate(time=">1997")

bind_rows(trueSwitchers_dietSum_1997,trueSwitchers_dietSum_2019)%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,prop,fill=time),position=position_dodge(),color="black")+
  geom_text(data=trueSwitchers_dietSum_1997,aes(x=0.5,y=0.28,label=paste("<1998:",nameClust)),hjust=0,color=viridis::viridis(2,end=0.7,option="B")[1])+
  geom_text(data=trueSwitchers_dietSum_2019,aes(x=0.5,y=0.25,label=paste(">1997:",nameClust)),hjust=0,color=viridis::viridis(2,end=0.7,option="B")[2])+
  scale_fill_viridis_d(name="Time Period",guide="none",end=0.7,option="B")+
  facet_wrap(~species,scales="free")+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.1))+
  xlab("Prey Item")+ylab("Proportion of Diet")

#Proportion true switchers
nrow(trueSwitchers)/nrow(sharedClusts)

###### FOO in each time Guild ######
#GL Time#
propMat_gl_guild<-left_join(sharedClusts%>%mutate(species=paste(species,sizecat,sep="_"))%>%
                              dplyr::select(species,guild=myGuild_1997,guildclust=guildclust_1997,
                                            cluster=myCluster_1997,namedClust=nameClust_1997),
                            propMat_gl)
propSum_gl_guild<-pivot_longer(propMat_gl_guild,cols=matches("[A-Z]",ignore.case=F),names_to="GL_pyscinam",values_to="prop")%>%
  group_by(guild,GL_pyscinam)%>%
  mutate(guildprop=mean(prop))%>%
  group_by(guild,guildclust,GL_pyscinam,guildprop)%>%
  summarise(clustprop=mean(prop))%>%
  arrange(guild,guildclust,desc(guildprop))%>%
  group_by(guild)%>%
  mutate(cumGuildProp=cumsum(guildprop))
propSum_gl_cluster<-pivot_longer(propMat_gl_guild,cols=matches("^[A-Z]",ignore.case=F),
                                 names_to="GL_pyscinam",values_to="prop")%>%
  group_by(cluster,guildclust,namedClust,GL_pyscinam)%>%
  summarise(clusterprop=mean(prop))%>%
  arrange(cluster,desc(clusterprop))%>%
  group_by(cluster)%>%
  mutate(cumClusterProp=cumsum(clusterprop),
         nItems_to90=sum(cumClusterProp<0.90)+1,
         GL_pyscinam_labels=str_pad(GL_pyscinam,max(str_length(preyTrawls$GL_pyscinam),na.rm=T),side="left"))

majorPropSum_gl<-propSum_gl_guild%>%
  mutate(cumGuildProp=gsub("[1-9]\\.","0\\.",cumGuildProp))%>%
  filter(cumGuildProp<=0.9)
majorClusterItems_gl<-propSum_gl_cluster%>%
  group_by(cluster)%>%
  slice(1:nItems_to90)

guild1<-filter(majorPropSum_gl,guild==1)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild2<-filter(majorPropSum_gl,guild==2)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild3<-filter(majorPropSum_gl,guild==3)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild4<-filter(majorPropSum_gl,guild==4)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild5<-filter(majorPropSum_gl,guild==5)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild6<-filter(majorPropSum_gl,guild==6)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild7<-filter(majorPropSum_gl,guild==7)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild8<-filter(majorPropSum_gl,guild==8)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")

grid.arrange(guild1,guild2,guild3,guild4,guild5,nrow=5)


#Doing each cluster on its own

cluster1<-filter(majorClusterItems_gl,guildclust=="1a")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="1a")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="1a")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster2<-filter(majorClusterItems_gl,guildclust=="1b")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="1b")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="1b")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster3<-filter(majorClusterItems_gl,guildclust=="1c")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="1c")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="1c")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster4<-filter(majorClusterItems_gl,guildclust=="1d")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="1d")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="1d")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster5<-filter(majorClusterItems_gl,guildclust=="1e")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="1e")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="1e")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster6<-filter(majorClusterItems_gl,guildclust=="2a")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="2a")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="2a")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster7<-filter(majorClusterItems_gl,guildclust=="2b")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="2b")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="2b")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster8<-filter(majorClusterItems_gl,guildclust=="3a")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="3a")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="3a")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster9<-filter(majorClusterItems_gl,guildclust=="3b")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="3b")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="3b")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster10<-filter(majorClusterItems_gl,guildclust=="3c")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="3c")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="3c")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster11<-filter(majorClusterItems_gl,guildclust=="3d")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="3d")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="3d")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster12<-filter(majorClusterItems_gl,guildclust=="4a")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="4a")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="4a")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster13<-filter(majorClusterItems_gl,guildclust=="5a")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",guildclust,", ",namedClust,"; Species=",filter(guildSummaries_gl,guildclust_1997=="5a")$N_species,
                     "; Diets=",filter(guildSummaries_gl,guildclust_1997=="5a")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")

grid.arrange(cluster1,cluster2,cluster3,cluster4,cluster5,cluster6,cluster7,
             cluster8,cluster9,cluster10,cluster11,cluster12,cluster13,ncol=3)


#NEW Time#
propMat_new_guild<-left_join(sharedClusts%>%mutate(species=paste(species,sizecat,sep="_"))%>%
                               dplyr::select(species,guild=myGuild_2019,guildclust=guildclust_2019,
                                             cluster=myCluster_2019,namedClust=nameClust_2019),
                             propMat_new)
propSum_new_guild<-pivot_longer(propMat_new_guild,cols=matches("^[A-Z]",ignore.case=F),
                                names_to="GL_pyscinam",values_to="prop")%>%
  group_by(guild,GL_pyscinam)%>%
  mutate(guildprop=mean(prop))%>%
  group_by(guild,guildclust,GL_pyscinam,guildprop)%>%
  summarise(clustprop=mean(prop))%>%
  arrange(guild,guildclust,desc(guildprop))%>%
  group_by(guild)%>%
  mutate(cumGuildProp=cumsum(guildprop))
propSum_new_cluster<-pivot_longer(propMat_new_guild,cols=matches("^[A-Z]",ignore.case=F),
                                  names_to="GL_pyscinam",values_to="prop")%>%
  group_by(cluster,namedClust,GL_pyscinam)%>%
  summarise(clusterprop=mean(prop))%>%
  arrange(cluster,desc(clusterprop))%>%
  group_by(cluster)%>%
  mutate(cumClusterProp=cumsum(clusterprop),
         nItems_to90=sum(cumClusterProp<0.90)+1,
         GL_pyscinam_labels=str_pad(GL_pyscinam,max(str_length(preyTrawls$GL_pyscinam),na.rm=T),side="left"))

majorPropSum_new<-propSum_new_guild%>%
  mutate(cumGuildProp=gsub("[1-9]\\.","0\\.",cumGuildProp))%>%
  filter(cumGuildProp<=0.9)
majorClusterItems_new<-propSum_new_cluster%>%
  group_by(cluster)%>%
  slice(1:nItems_to90)


guild1<-filter(majorPropSum_new,guild==1)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild2<-filter(majorPropSum_new,guild==2)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild3<-filter(majorPropSum_new,guild==3)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild4<-filter(majorPropSum_new,guild==4)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild5<-filter(majorPropSum_new,guild==5)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild6<-filter(majorPropSum_new,guild==6)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild7<-filter(majorPropSum_new,guild==7)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")
guild8<-filter(majorPropSum_new,guild==8)%>%
  mutate(GL_pyscinam=fct_reorder(GL_pyscinam,guildprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam,clustprop,fill=guildclust),color="black",position="dodge")+
  facet_wrap(~paste0("Guild=",guild))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none")

grid.arrange(guild1,guild2,guild3,guild4,guild5,guild6,guild7,guild8,nrow=4)


cluster1<-filter(majorClusterItems_new,cluster=="1")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",cluster,", ",namedClust,"; Species=",filter(guildSummaries_new,myCluster_2019=="1")$N_species,
                     "; Diets=",filter(guildSummaries_new,myCluster_2019=="1")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster2<-filter(majorClusterItems_new,cluster=="2")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",cluster,", ",namedClust,"; Species=",filter(guildSummaries_new,myCluster_2019=="2")$N_species,
                     "; Diets=",filter(guildSummaries_new,myCluster_2019=="2")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster3<-filter(majorClusterItems_new,cluster=="3")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",cluster,", ",namedClust,"; Species=",filter(guildSummaries_new,myCluster_2019=="3")$N_species,
                     "; Diets=",filter(guildSummaries_new,myCluster_2019=="3")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster4<-filter(majorClusterItems_new,cluster=="4")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",cluster,", ",namedClust,"; Species=",filter(guildSummaries_new,myCluster_2019=="4")$N_species,
                     "; Diets=",filter(guildSummaries_new,myCluster_2019=="4")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")
cluster5<-filter(majorClusterItems_new,cluster=="5")%>%
  mutate(GL_pyscinam_labels=fct_reorder(GL_pyscinam_labels,clusterprop,.desc=T))%>%
  ggplot()+
  geom_col(aes(GL_pyscinam_labels,clusterprop),color="black",position="dodge")+
  facet_wrap(~paste0("Cluster=",cluster,", ",namedClust,"; Species=",filter(guildSummaries_new,myCluster_2019=="5")$N_species,
                     "; Diets=",filter(guildSummaries_new,myCluster_2019=="5")$N_diets))+
  theme(axis.text.x=element_text(angle=30,hjust=1,vjust=1.05),
        legend.position="none",axis.title.x=element_blank())+
  ylab("Proportion of Diet")

grid.arrange(cluster1,cluster2,cluster3,cluster4,cluster5,nrow=5)



##### Tanglegram #####


abs(order.dendrogram(dend_shared_new)-rev(order.dendrogram(dend_shared_new)))^1
#Doing the more robust comparisons (traditonally used) requires they have identical groups
dend_cor<-dendlist(dend_shared_gl,dend_shared_new)%>%
  cor_cophenetic()
dend_cor2<-dendlist(dend_shared_gl,dend_shared_new)%>%
  cor_bakers_gamma()
dend_tangle<-dendlist(dend_shared_gl,dend_shared_new)%>%
  untangle(method="step2side")%>%
  entanglement()
dendlist(dend_shared_gl,dend_shared_new)%>%
  untangle(method="step2side")%>%
  plot(common_subtrees_color_branches=T,highlight_distinct_edges=F,highlight_branches_lwd=F,
       edge.lwd=2,margin_inner=12,
       main=paste0("1973-1997       |       1998-2019"))

cophenetic(dend_shared_gl)

##### Overlap from one time to another #####
crossTime<-full_join(props1_gl,props2_new,
                     by=c("species1"="species2",
                          "sizecat1"="sizecat2",
                          "GL_pyscinam"))
crossTime[is.na(crossTime)]<-0

crossTime<-crossTime%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

ggplot(crossTime,aes(species1,s_do))+
  geom_boxplot(outlier.shape=NA,color="black")+
  geom_point(aes(fill=factor(sizecat1,levels=c("S","M","L","XL")),
                 shape=factor(sizecat1,levels=c("S","M","L","XL"))),
             size=3.5)+
  scale_shape_manual(name="Size Class",
                     values=c(21,22,23,24))+
  scale_fill_viridis_d(name="Size Class")+
  scale_y_continuous(limits=c(0,1),
                     name="Overlap")+
  scale_x_discrete(name="Species")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))





# NMDS work ---------------------------------------------------------------

#Total amount of catch of all species in each trawl
trawlHauls<-preyTrawls_nSize%>%
  ungroup()%>%
  dplyr::select(id,pdcomnam,abundance)%>%distinct()%>%
  group_by(id)%>%
  summarise(trawlabundance=sum(abundance,na.rm=T),
            trawlabundance=ifelse(trawlabundance==0,1,trawlabundance))

#Need to have whole diet totals, and whole trawl measures
NMDStrawlDiets<-preyTrawls_nSize%>%
  group_by(id)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv),
         totalN=n_distinct(dietID))%>%
  left_join(trawlHauls)%>%
  group_by(id,GL_pyscinam)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(totalN==0,0,n_distinct(dietID)/totalN),
         doy=yday(est_towdate))%>%
  dplyr::select(year,season,geoarea,id,trawlabundance,
                GL_pyscinam,totalwt,totalv,totalN,qikw,qikv,pik,decdeg_beglat:botsalin,doy,-est_towdate)%>%
  group_by(id)%>%
  fill(decdeg_beglat:doy,.direction="downup")%>%
  distinct()



#if grouping, add those groups in here both at select and group_by
nDiets_NMDS<-NMDStrawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case the whole pred sizecat)
  ungroup()%>%
  dplyr::select(year,geoarea,id,totalN)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(year,geoarea)%>%summarise(nDiets=sum(totalN))#Calculate the number of diets
sumAbun_NMDS<-NMDStrawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(year,geoarea,id,trawlabundance)%>%distinct()%>%
  group_by(year,geoarea)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty_NMDS<-NMDStrawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(year,geoarea,id,trawlabundance)%>%distinct()%>%
  group_by(year,geoarea)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))

NMDSyearsum<-NMDStrawlDiets%>%
  group_by(year,geoarea,season)%>%
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
  group_by(year,geoarea,season,nTrawls,nDiets,sumAbun,sumAbun_nEmpty,surftemp,bottemp,avgdepth,surfsalin,botsalin,doy,GL_pyscinam)%>%
  summarise(sumAbun_nEmpty=ifelse(is.na(sumAbun_nEmpty),nDiets,sumAbun_nEmpty),
            Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik)/sumAbun, #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,trawlabundance,pik,Fk))%>%
  distinct()

check<-NMDSyearsum%>%group_by(geoarea,year,season)%>%summarise(N=sum(Wk))
sum(check$N<0&check$N>0) #When this is 0, then you have all your means correct because they sum to 1 (or they're empties)


idCols<-c("year","geoarea","season","surftemp","bottemp","avgdepth","surfsalin","botsalin","doy")
NMDSmat<-pivot_wider(NMDSyearsum,id_cols=all_of(idCols),
                                 names_from="GL_pyscinam",values_from="Wk")%>%
  dplyr::select(-c(Empty,"NA"))

library(vegan)
library(RVAideMemoire)
library(usedist)
library(indicspecies)

NMDSmat.diet<-(NMDSmat[,(length(idCols)+1):ncol(NMDSmat)])
NMDSmat.env<-(NMDSmat[,1:length(idCols)])

#Percent zeroes--number of zeroes times dim of matrix
(sum(is.na(NMDSmat.diet))/(nrow(NMDSmat.diet)*ncol(NMDSmat.diet)))*100

NMDSmat.diet<-NMDSmat.diet%>%
  mutate(sumRow=rowSums(NMDSmat.diet,na.rm=T))%>%
  filter(sumRow>0)%>%
  dplyr::select(-sumRow)%>%
  as.matrix()
NMDSmat.diet[is.na(NMDSmat.diet)]<-0

#Using the whole dataset--but can use the one without minor species (Jellyfish, Sea Angel)
set.seed(42)

#For the all diets NMDS
original.dist<-vegdist(NMDSmat.diet)
stress_values<-numeric(6)
r2<-numeric(6)

for (n in 1:6) {
  nmds.resu <- metaMDS(NMDSmat.diet, k=n, distance = "bray", try=250, autotransform=F)
  stress_values[n]<-nmds.resu$stress
  nmds.scores<-vegan::scores(nmds.resu)
  nmds.dist<-dist(nmds.scores)
  r2[n]<-summary(lm(original.dist~nmds.dist))[[8]]
}
plot(stress_values, xlab="Number of axes", ylab="Stress",type="b")
abline(h=0.2,col="red")

View(stress_values) 

#Go back and create the output for the 2 dimensions NMDS
seasonNMDS<-metaMDS(NMDSmat.diet, distance = "bray", k = 3, try=250, autotransform=F,na.rm=T)
r2<-summary(lm(original.dist~dist(vegan::scores(count_NMDS))))[[8]]
actualStress<-count_NMDS$stress
stressplot(count_NMDS) #This is the visual of stress, the divergence of observed and ordinated distance. It's random, that's good

#Writing it out to use in PC-ORD for comparison and for the axes R2 values
#temp<-cbind(allEnv.Mat[,1],allCount.Mat)
#write.csv(temp,"/Users/nh1087/Documents/NMDS_Matrix.csv",row.names = F)



#Print the species scores and sample scores
NMDS_species<-as.data.frame(count_NMDS$species)
NMDS_scores<-as.data.frame(count_NMDS$points)%>%
  bind_cols(NMDSmat.env)%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))
NMDS_scores[is.na(NMDS_scores)]<-NA

#With seasons
NMDS_species<-as.data.frame(seasonNMDS$species)
NMDS_scores<-as.data.frame(seasonNMDS$points)%>%
  bind_cols(NMDSmat.env)%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))
NMDS_scores[is.na(NMDS_scores)]<-NA



basinCentroids <- NMDS_scores%>%
  group_by(geoarea)%>%
  summarise(MDS1=mean(MDS1),
            MDS2=mean(MDS2),
            MDS3=mean(MDS3))

set.seed(42)
adonis(original.dist~year*geoarea,data=NMDS_scores,permutations=1000,method="bray")

summary(NMDS_lm<-lm(cbind(MDS1,MDS2,MDS3)~year*geoarea,data=NMDS_scores))
library(car)
Anova(NMDS_lm)

predLine<-cbind(NMDS_scores[,c(4,5)],data.frame(predict(NMDS_lm)))
predWide<-predLine%>%
  filter(year==2019 | year==1973)%>%
  pivot_wider(id_cols=geoarea,names_from=year,values_from=starts_with("MDS"))



##### Vegan vector fitting ####
fit<-envfit(seasonNMDS,NMDS_scores[,c("year","geoarea","season","bottemp","avgdepth","doy")],choices=1:3,na.rm=T)
fit.df<-data.frame(fit[["vectors"]]$arrows)%>%
  rownames_to_column(var="var")


##### Plotting NMDS, Highlighting yearly changes, but keeping basins involved #####
axes12<-ggplot(data=NMDS_scores,aes(MDS1,MDS2))+
  geom_point(data=NMDS_scores,aes(MDS1,MDS2,fill=year),shape=21,color="black",size=6)+
  geom_segment(data=fit.df,aes(x=0,xend=NMDS1,
                               y=0,yend=NMDS2),
               color="blue",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  geom_label(data=fit.df,aes(x=NMDS1*1.1,y=NMDS2*1.1,label=var))+
  geom_point(data=basinCentroids,aes(MDS1,MDS2),
             color="red",size=12,stroke=2,alpha=0.7,show.legend=F,shape=21)+
  geom_label(data=basinCentroids,aes(MDS1,MDS2,label=geoarea),
             color="red",size=3,alpha=0.7,show.legend=F)+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"))+
  scale_y_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab("Axis 1")+ylab("Axis 2")+
  facet_wrap(~geoarea,nrow=1)


axes13<-ggplot(data=NMDS_scores,aes(MDS1,MDS3))+
  geom_point(data=NMDS_scores,aes(MDS1,MDS3,fill=year),shape=21,color="black",size=6)+
  geom_segment(data=fit.df,aes(x=0,xend=NMDS1,
                               y=0,yend=NMDS3),
               color="blue",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  geom_label(data=fit.df,aes(x=NMDS1*1.1,y=NMDS3*1.1,label=var))+
  geom_point(data=basinCentroids,aes(MDS1,MDS3),
             color="red",size=12,stroke=2,alpha=0.7,show.legend=F,shape=21)+
  geom_label(data=basinCentroids,aes(MDS1,MDS3,label=geoarea),
             color="red",size=3,alpha=0.7,show.legend=F)+
  #geom_segment(data=NMDS_species,aes(x=0,xend=MDS1,y=0,yend=MDS2))+
  #geom_text(data=NMDS_species,aes(x=MDS1,y=MDS2,label=row.names(NMDS_species)))+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"))+
  scale_y_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab("Axis 1")+ylab("Axis 3")+
  #xlim(-2.4,2.4)+ylim(-2.4,2.4)+
  facet_wrap(~geoarea,nrow=1)


axes23<-ggplot(data=NMDS_scores,aes(MDS2,MDS3))+
  geom_point(data=NMDS_scores,aes(MDS2,MDS3,fill=year),shape=21,color="black",size=6)+
  geom_segment(data=fit.df,aes(x=0,xend=NMDS2,
                               y=0,yend=NMDS3),
               color="blue",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  geom_label(data=fit.df,aes(x=NMDS2*1.1,y=NMDS3*1.1,label=var))+
  geom_point(data=basinCentroids,aes(MDS2,MDS3),
             color="red",size=12,stroke=2,alpha=0.7,show.legend=F,shape=21)+
  geom_label(data=basinCentroids,aes(MDS2,MDS3,label=geoarea),
             color="red",size=3,alpha=0.7,show.legend=F)+
  #geom_segment(data=NMDS_species,aes(x=0,xend=MDS1,y=0,yend=MDS2))+
  #geom_text(data=NMDS_species,aes(x=MDS1,y=MDS2,label=row.names(NMDS_species)))+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"))+
  scale_y_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab("Axis 2")+ylab("Axis 3")+
  #xlim(-2.4,2.4)+ylim(-2.4,2.4)+
  facet_wrap(~geoarea,nrow=1)


ggarrange(axes12,axes13,axes23,common.legend = T,legend="top",ncol=1)

actualStress
r2
max(NMDS_scores$MDS3)

#Plotting for the presentation (one year only, first with just the one survey then two)
axes12.1s<-ggplot(data=filter(NMDS_scores,year=="1997" & season=="Spring"),aes(MDS1,MDS2))+
  geom_point(data=filter(NMDS_scores,year=="1997" & season=="Spring"),aes(MDS1,MDS2),color="black",size=6)+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"),limits=c(-1.21,1.02))+
  scale_y_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"),limits=c(-1.21,2.39))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab("Axis 1")+ylab("Axis 2")+
  facet_wrap(~geoarea,nrow=1)


axes13.1s<-ggplot(data=filter(NMDS_scores,year=="1997" & season=="Spring"),aes(MDS1,MDS3))+
  geom_point(data=filter(NMDS_scores,year=="1997" & season=="Spring"),aes(MDS1,MDS3),color="black",size=6)+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"),limits=c(-1.21,1.02))+
  scale_y_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"),limits=c(-1.1,1.17))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab("Axis 1")+ylab("Axis 3")+
  facet_wrap(~geoarea,nrow=1)


axes23.1s<-ggplot(data=filter(NMDS_scores,year=="1997" & season=="Spring"),aes(MDS2,MDS3))+
  geom_point(data=filter(NMDS_scores,year=="1997" & season=="Spring"),aes(MDS2,MDS3),color="black",size=6)+
  scale_x_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"),limits=c(-1.21,2.39))+
  scale_y_continuous(breaks=c(-1,0,1,2),labels=c("-1","0","1","2"),limits=c(-1.1,1.17))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab("Axis 2")+ylab("Axis 3")+
  facet_wrap(~geoarea,nrow=1)


ggarrange(axes12.1s,axes13.1s,axes23.1s,common.legend = T,legend="top",ncol=1)



# Same thing with PCA (easier to explain and more common) -----------------

PCA<-princomp(NMDSmat.diet)

PCA_vars<-PCA$sdev^2 / sum(PCA$sdev^2)

PCA_species<-as.data.frame(PCA$loadings[,1:3])
PCA_scores<-as.data.frame(PCA$scores[,1:3])%>%
  bind_cols(NMDSmat.env)%>%
  mutate(geoarea=factor(geoarea,levels=c("MAB","SNE","GB","GoM","ScS")))


basinPCA_Centroids <- PCA_scores%>%
  group_by(geoarea)%>%
  summarise(Comp.1=mean(Comp.1),
            Comp.2=mean(Comp.2),
            Comp.3=mean(Comp.3))

summary(PCA_lm<-lm(cbind(Comp.1,Comp.2,Comp.3)~year*geoarea,data=PCA_scores))

Anova(PCA_lm)

PCA_predLine<-cbind(PCA_scores[,c(4,5)],data.frame(predict(PCA_lm)))
PCA_predWide<-PCA_predLine%>%
  filter(year==2019 | year==1973)%>%
  pivot_wider(id_cols=geoarea,names_from=year,values_from=starts_with("Comp"))

pca12<-ggplot()+
  geom_point(data=PCA_scores,aes(Comp.1,Comp.2,fill=year),shape=21,color="black",size=6)+
  geom_segment(data=PCA_predWide,aes(x=Comp.1_1973,xend=Comp.1_2019,
                                 y=Comp.2_1973,yend=Comp.2_2019),
               color="red",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  geom_point(data=basinPCA_Centroids,aes(Comp.1,Comp.2),
             color="red",size=12,stroke=2,alpha=0.7,show.legend=F,shape=21)+
  geom_label(data=basinPCA_Centroids,aes(Comp.1,Comp.2,label=geoarea),
             color="red",size=3,alpha=0.7,show.legend=F)+
  #geom_segment(data=NMDS_species,aes(x=0,xend=MDS1,y=0,yend=MDS2))+
  #geom_text(data=NMDS_species,aes(x=MDS1,y=MDS2,label=row.names(NMDS_species)))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank(),legend.position="top")+
  xlab(paste0("Axis 1 (",round(PCA_vars[1],3)*100,"%)"))+ylab(paste0("Axis 2 (",round(PCA_vars[2],3)*100,"%)"))+
  xlim(-.5,.36)+ylim(-.5,.36)+
  facet_wrap(~geoarea,nrow=1)
pca12

pca13<-ggplot()+
  geom_point(data=PCA_scores,aes(Comp.1,Comp.3,fill=year),shape=21,color="black",size=6)+
  geom_segment(data=PCA_predWide,aes(x=Comp.1_1973,xend=Comp.1_2019,
                                 y=Comp.3_1973,yend=Comp.3_2019),
               color="red",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  geom_point(data=basinPCA_Centroids,aes(Comp.1,Comp.3),
             color="red",size=12,stroke=2,alpha=0.7,show.legend=F,shape=21)+
  geom_label(data=basinPCA_Centroids,aes(Comp.1,Comp.3,label=geoarea),
             color="red",size=3,alpha=0.7,show.legend=F)+
  #geom_segment(data=NMDS_species,aes(x=0,xend=MDS1,y=0,yend=MDS2))+
  #geom_text(data=NMDS_species,aes(x=MDS1,y=MDS2,label=row.names(NMDS_species)))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab(paste0("Axis 1 (",round(PCA_vars[1],3)*100,"%)"))+ylab(paste0("Axis 3 (",round(PCA_vars[3],3)*100,"%)"))+
  xlim(-.5,.36)+ylim(-.5,.36)+
  facet_wrap(~geoarea,nrow=1)


pca23<-ggplot()+
  geom_point(data=PCA_scores,aes(Comp.2,Comp.3,fill=year),shape=21,color="black",size=6)+
  geom_segment(data=PCA_predWide,aes(x=Comp.2_1973,xend=Comp.2_2019,
                                     y=Comp.3_1973,yend=Comp.3_2019),
               color="red",alpha=0.7,lwd=2,arrow=arrow(length=unit(0.1,"inches")),lineend="round",linejoin="round")+
  geom_point(data=basinPCA_Centroids,aes(Comp.2,Comp.3),
             color="red",size=12,stroke=2,alpha=0.7,show.legend=F,shape=21)+
  geom_label(data=basinPCA_Centroids,aes(Comp.2,Comp.3,label=geoarea),
             color="red",size=3,alpha=0.7,show.legend=F)+
  #geom_segment(data=NMDS_species,aes(x=0,xend=MDS1,y=0,yend=MDS2))+
  #geom_text(data=NMDS_species,aes(x=MDS1,y=MDS2,label=row.names(NMDS_species)))+
  scale_fill_viridis_c(option="E",name="Year")+
  theme(legend.key.width=unit(2.5,"cm"),strip.text=element_blank())+
  xlab(paste0("Axis 2 (",round(PCA_vars[2],3)*100,"%)"))+ylab(paste0("Axis 3 (",round(PCA_vars[3],3)*100,"%)"))+
  xlim(-.5,.36)+ylim(-.5,.36)+
  facet_wrap(~geoarea,nrow=1)


ggarrange(pca12,pca13,pca23,common.legend = T,legend="top",ncol=1)




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

indResponses.good<-filter(indResponses, paste(Species,season) %in% paste(goodSpecies$Species,goodSpecies$season))






# GLMs for individual responses -------------------------------------------
library(lme4)
#Truly for all the species
sp.pEmpty.std<-glm(cbind(nEmpty,nDiets)~Species*year*season+surftemp+DOY,
                   data=filter(indResponses.good,!is.na(surftemp)),family=binomial())
summary(sp.pEmpty.std)
sp.pEmpty.s<-glm(cbind(nEmpty,nDiets)~Species*year*season,data=indResponses.good,family=binomial())
summary(sp.pEmpty.s)
Anova(sp.pEmpty.s,sp.pEmpty.std)

sp.plotResponses.pe<-data.frame(indResponses.good,pred=predict(sp.pEmpty.s,type="response",se=T))%>%
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



#Diet Breadth
sp.dBreath<-glm(simpson~Species*year*season,data=indResponses.good,family=quasibinomial(logit))
summary(sp.dBreath)

sp.plotResponses.db<-data.frame(filter(indResponses.good,!is.na(simpson)),pred=predict(sp.dBreath,type="response",se=T))%>%
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



#Relative Consumption
sp.rConsump<-glm(meanC~Species*year*season,data=indResponses.good,family=Gamma(link="log"))
summary(sp.rConsump)

sp.plotResponses.rc<-data.frame(filter(indResponses.good,!is.na(meanC)),pred=predict(sp.rConsump,type="response",se=T))%>%
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



# Post-hoc Analysis -------------------------------------------------------

#Factors to include
  #Predator guild
  #Predator climate sensitivity
  #Seasonal migration
  #Range centroid shift
  #Foraging specificity (specialist v generalist)

reducedClusts<-sharedClusts%>%
  group_by(species,nameClust_1997)%>%
  mutate(N_sizes=n(),
         Species=species)%>%
  group_by(Species)%>%
  filter(N_sizes==max(N_sizes))%>%
  dplyr::select(Species,guild=nameClust_1997)%>%distinct()
doubleClusts<-bind_rows(reducedClusts[duplicated(reducedClusts$Species),],
                        as.data.frame(t(rev(as.data.frame(t(reducedClusts)))))[duplicated(rev(reducedClusts$Species)),])
#I want to keep the piscivores for all these and amphi for Little Skate
reducedClusts<-reducedClusts%>%
  group_by(Species)%>%
  mutate(N=n(),
         guild=case_when(N==1~guild,
                         N==2 & Species=="Little Skate"~"Polychaete/Amphipod eaters",
                         N==2 & Species!="Little Skate"~"Piscivores",
                         T~"NA"))%>%distinct()%>%dplyr::select(-N)
finalClusts<-sharedClusts%>%
  mutate(sizecat2=as.numeric(factor(sizecat,levels=c("S","M","L","XL"))))%>%
  group_by(species)%>%
  filter(sizecat2==max(sizecat2))%>%
  dplyr::select(Species=species,finalguild=nameClust_1997)




#Percent Empty
pChange.drivers.pe<-left_join(sp.responseCats.pe,hare_everything_FHD)%>%
  left_join(reducedClusts)%>%
  left_join(finalClusts)%>%
  left_join(dplyr::select(filter(hare_qualScores,Attribute=="Prey Specificity"),
                          Species,Attribute,meanScore))%>%
  ungroup()%>%
  mutate(hare_cVul=fct_explicit_na(Climate.Vulnerability,na_level="Unassessed"))
  

#Simple ANOVA concept
summary(posthoc.aov.pe<-aov(pChange~hare_cVul*finalguild+Error(Species),data=pChange.drivers.pe))
#One up, since prey specificity is continuous

summary(posthoc.lm.pe<-lmer(pChange~meanScore*finalguild+(1|Species),
                            data=pChange.drivers.pe))
#Plotting
pChange.drivers.pe%>%
  group_by(guild,hare_cVul,.drop=F)%>%
  summarise(mean=mean(pChange-1), #Do minus 1 to adjust to the anomaly (if positive, increase and negative is decrease)
            lower=mean-(sd(pChange)/sqrt(n())),
            upper=mean+(sd(pChange)/sqrt(n())))%>%
  ggplot()+
  geom_col(aes(guild,mean,fill=hare_cVul),color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(guild,ymin=lower,ymax=upper,group=hare_cVul),position=position_dodge(0.9),width=0.2)
pChange.drivers.pe%>%
  mutate(guild=fct_reorder(as.factor(finalguild),pChange))%>%
  complete(guild,hare_cVul)%>%
  mutate(pChange2=ifelse(is.na(pChange),0,pChange-1))%>%
  ggplot()+
  geom_boxplot(aes(guild,pChange2*100,fill=hare_cVul),outlier.shape=NA)+
  geom_jitter(aes(guild,(pChange-1)*100,fill=hare_cVul),shape=21,position=position_jitterdodge(0.2),size=5)+
  scale_x_discrete(drop=F,name="Adult Feeding Guild")+
  scale_fill_discrete(drop=F,name="Climate Vulnerability\n(Hare et al. 2016)")+
  scale_y_continuous(name="Percent Change in Response (%Empty)",breaks=c(-25,0,50,100,150,200),
                     labels=c("-25%","0%","50%","100%","150%","200%"))+
  theme(legend.position=c(0.15,0.8),legend.background = element_rect(color="black"))

pChange.drivers.pe%>%
  mutate(guild=fct_reorder(as.factor(finalguild),pChange))%>%
  complete(guild,hare_cVul)%>%
  mutate(pChange2=ifelse(is.na(pChange),0,pChange-1))%>%
  ggplot()+
  geom_point(aes(meanScore,(pChange-1)*100,fill=guild),size=5,shape=21)+
  geom_smooth(aes(meanScore,(pChange-1)*100,color=guild,group=guild),method="lm",alpha=0.3,size=2,show.legend=F)+
  scale_x_continuous(name="Prey Specificity (Hare et al. 2016)",breaks=c(1,2),
                     labels=c("             Generalist","Specialist           "))+
  scale_y_continuous(name="Percent Change in Response (%Empty)",breaks=c(-50,0,50,100,150,200),
                     labels=c("-50%","0%","50%","100%","150%","200%"))+
  scale_fill_viridis_d(end=0.94,name="Adult Feeding Guild")+
  scale_color_viridis_d(end=0.94,name="Adult Feeding Guild")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank(),
        legend.position=c(0.15,0.8),legend.background = element_rect(color="black"))


#Relative Consumption
pChange.drivers.rc<-left_join(sp.responseCats.rc,hare_everything_FHD)%>%
  left_join(reducedClusts)%>%
  left_join(finalClusts)%>%
  left_join(dplyr::select(filter(hare_qualScores,Attribute=="Prey Specificity"),
                          Species,Attribute,meanScore))%>%
  ungroup()%>%
  mutate(hare_cVul=fct_explicit_na(Climate.Vulnerability,na_level="Unassessed"))


#Simple ANOVA concept
summary(posthoc.aov.rc<-aov(pChange~hare_cVul*finalguild+Error(Species),data=pChange.drivers.rc))
#One up, since prey specificity is continuous

summary(posthoc.lm.rc<-lmer(pChange~meanScore+(1|Species),
                            data=pChange.drivers.rc))
#Plotting
pChange.drivers.rc%>%
  group_by(guild,hare_cVul,.drop=F)%>%
  summarise(mean=mean(pChange-1), #Do minus 1 to adjust to the anomaly (if positive, increase and negative is decrease)
            lower=mean-(sd(pChange)/sqrt(n())),
            upper=mean+(sd(pChange)/sqrt(n())))%>%
  ggplot()+
  geom_col(aes(guild,mean,fill=hare_cVul),color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(guild,ymin=lower,ymax=upper,group=hare_cVul),position=position_dodge(0.9),width=0.2)
pChange.drivers.rc%>%
  mutate(guild=fct_reorder(as.factor(finalguild),pChange))%>%
  complete(guild,hare_cVul)%>%
  mutate(pChange2=ifelse(is.na(pChange),0,pChange-1))%>%
  ggplot()+
  geom_boxplot(aes(guild,pChange2*100,fill=hare_cVul),outlier.shape=NA)+
  geom_jitter(aes(guild,(pChange-1)*100,fill=hare_cVul),shape=21,position=position_jitterdodge(0.2),size=5)+
  scale_x_discrete(drop=F,name="Adult Feeding Guild")+
  scale_fill_discrete(drop=F,name="Climate Vulnerability\n(Hare et al. 2016)")+
  scale_y_continuous(name="Percent Change in Response (Relative Consumption)",breaks=c(-100,-75,-50,-25,0,25),
                     labels=c("-100%","-75%","-50%","-25%","0%","25%"))+
  theme(legend.position=c(0.15,0.81),legend.background = element_rect(color="black"))

pChange.drivers.rc%>%
  mutate(guild=fct_reorder(as.factor(finalguild),pChange))%>%
  complete(guild,hare_cVul)%>%
  mutate(pChange2=ifelse(is.na(pChange),0,pChange-1))%>%
  ggplot()+
  geom_point(aes(meanScore,(pChange-1)*100,fill=guild),size=5,shape=21)+
  geom_smooth(aes(meanScore,(pChange-1)*100,color=guild,group=guild),method="lm",alpha=0.3,size=2,show.legend=F)+
  scale_x_continuous(name="Prey Specificity (Hare et al. 2016)",breaks=c(1,2),
                     labels=c("             Generalist","Specialist           "))+
  scale_y_continuous(name="Percent Change in Response (Relative Consumption)",breaks=c(-100,-75,-50,-25,0,25),
                     labels=c("-100%","-75%","-50%","-25%","0%","25%"))+
  scale_fill_viridis_d(end=0.94,name="Adult Feeding Guild")+
  scale_color_viridis_d(end=0.94,name="Adult Feeding Guild")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank(),
        legend.position=c(0.15,0.18),legend.background = element_rect(color="black"))



#Diet Breadth
pChange.drivers.db<-left_join(sp.responseCats.db,hare_everything_FHD)%>%
  left_join(reducedClusts)%>%
  left_join(finalClusts)%>%
  left_join(dplyr::select(filter(hare_qualScores,Attribute=="Prey Specificity"),
                          Species,Attribute,meanScore))%>%
  ungroup()%>%
  mutate(hare_cVul=fct_explicit_na(Climate.Vulnerability,na_level="Unassessed"))


#Simple ANOVA concept
summary(posthoc.aov.db<-aov(pChange~hare_cVul*finalguild+Error(Species),data=pChange.drivers.db))
#One up, since prey specificity is continuous

summary(posthoc.lm.db<-lm(log(pChange)~meanScore,
                              data=pChange.drivers.db))
(exp(coef(posthoc.lm.db)["meanScore"]) - 1) * 100
test<-data.frame(filter(pChange.drivers.db,!is.na(meanScore)),pred=predict(posthoc.lm.db,se.fit=T))%>%
  mutate(ci.lower=pred.fit-1.96*pred.se.fit,
         ci.upper=pred.fit+1.96*pred.se.fit,
         pred.Nat=exp(pred.fit),
         lower.Nat=exp(ci.lower),
         upper.Nat=exp(ci.upper))

#Plotting
pChange.drivers.db%>%
  group_by(guild,hare_cVul,.drop=F)%>%
  summarise(mean=mean(pChange-1), #Do minus 1 to adjust to the anomaly (if positive, increase and negative is decrease)
            lower=mean-(sd(pChange)/sqrt(n())),
            upper=mean+(sd(pChange)/sqrt(n())))%>%
  ggplot()+
  geom_col(aes(guild,mean,fill=hare_cVul),color="black",position=position_dodge(0.9))+
  geom_errorbar(aes(guild,ymin=lower,ymax=upper,group=hare_cVul),position=position_dodge(0.9),width=0.2)
pChange.drivers.db%>%
  mutate(guild=fct_reorder(as.factor(finalguild),pChange))%>%
  complete(guild,hare_cVul)%>%
  mutate(pChange2=ifelse(is.na(pChange),0,pChange-1))%>%
  ggplot()+
  geom_boxplot(aes(guild,pChange2*100,fill=hare_cVul),outlier.shape=NA)+
  geom_jitter(aes(guild,(pChange-1)*100,fill=hare_cVul),shape=21,position=position_jitterdodge(0.2),size=5)+
  scale_x_discrete(drop=F,name="Adult Feeding Guild")+
  scale_fill_discrete(drop=F,name="Climate Vulnerability\n(Hare et al. 2016)")+
  scale_y_continuous(name="Percent Change in Response (Diet Breadth)",breaks=c(-100,-75,-50,-25,0,25,50,75,100),
                     labels=c("-100%","-75%","-50%","-25%","0%","25%","50%","75%","100%"))+
  theme(legend.position=c(0.8,0.81),legend.background = element_rect(color="black"))

pChange.drivers.db%>%
  mutate(guild=fct_reorder(as.factor(finalguild),pChange))%>%
  complete(guild,hare_cVul)%>%
  mutate(pChange2=ifelse(is.na(pChange),0,pChange-1))%>%
  ggplot()+
  geom_point(aes(meanScore,(pChange-1)*100),size=5,shape=19)+
  geom_smooth(aes(meanScore,(pChange-1)*100),method="lm",alpha=0.3,size=2,show.legend=F)+
  geom_hline(aes(yintercept=0),lty=2,size=2)+
  scale_x_continuous(name="Prey Specificity (Hare et al. 2016)",breaks=c(1,2),
                     labels=c("             Generalist","Specialist           "))+
  scale_y_continuous(name="Percent Change in Diet Breadth",breaks=c(-100,-75,-50,-25,0,25,50,75,100),
                     labels=c("-100%","-75%","-50%","-25%","0%","25%","50%","75%","100%"))+
  scale_fill_viridis_d(end=0.94,name="Adult Feeding Guild")+
  scale_color_viridis_d(end=0.94,name="Adult Feeding Guild")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank(),
        legend.position=c(0.8,0.8),legend.background = element_rect(color="black"))

hist(log(pChange.drivers.db$pChange))
ggplot(test)+
  geom_point(aes(meanScore,(pChange-1)*100),size=5,shape=19)+
  geom_ribbon(aes(ymin=(lower.Nat-1)*100,ymax=(upper.Nat-1)*100,x=meanScore),alpha=0.5,fill="grey50")+
  geom_line(aes(meanScore,(pred.Nat-1)*100),size=2,show.legend=F,color="blue")+
  geom_hline(aes(yintercept=0),lty=2,size=2)+
  scale_x_continuous(name="Prey Specificity (Hare et al. 2016)",breaks=c(1,2),
                     labels=c("             Generalist","Specialist             "))+
  scale_y_continuous(name="Percent Change in Diet Breadth",breaks=c(-100,-75,-50,-25,0,25,50,75,100),
                     labels=c("-100%","-75%","-50%","-25%","0%","25%","50%","75%","100%"))+
  scale_fill_viridis_d(end=0.94,name="Adult Feeding Guild")+
  scale_color_viridis_d(end=0.94,name="Adult Feeding Guild")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank(),
        legend.position=c(0.8,0.8),legend.background = element_rect(color="black"))

summary(posthoc.lm.db2<-lm(log(meanScore)~pChange,
                          data=pChange.drivers.db))
test2<-data.frame(filter(pChange.drivers.db,!is.na(meanScore)),pred=predict(posthoc.lm.db2,se.fit=T))%>%
  mutate(ci.lower=pred.fit-1.96*pred.se.fit,
         ci.upper=pred.fit+1.96*pred.se.fit,
         pred.Nat=exp(pred.fit),
         lower.Nat=exp(ci.lower),
         upper.Nat=exp(ci.upper))
ggplot(test2)+
  geom_point(aes((pChange-1)*100,meanScore),size=5,shape=19)+
  geom_ribbon(aes(ymin=lower.Nat,ymax=upper.Nat,x=(pChange-1)*100),alpha=0.5,fill="grey50")+
  geom_line(aes((pChange-1)*100,pred.Nat),size=2,show.legend=F,color="blue")+
  geom_vline(aes(xintercept=0),lty=2,size=2)+
  scale_y_continuous(name="Prey Specificity (Hare et al. 2016)",breaks=c(1,2),
                     labels=c("Generalist","Specialist"))+
  scale_x_continuous(name="Percent Change in Diet Breadth",breaks=c(-100,-75,-50,-25,0,25,50,75,100),
                     labels=c("-100%","-75%","-50%","-25%","0%","25%","50%","75%","100%"))+
  scale_fill_viridis_d(end=0.94,name="Adult Feeding Guild")+
  scale_color_viridis_d(end=0.94,name="Adult Feeding Guild")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank(),
        legend.position=c(0.8,0.8),legend.background = element_rect(color="black"))


#Combining all guilds for presentation
bind_rows(pChange.drivers.db,pChange.drivers.pe,pChange.drivers.rc)%>%
  ggplot()+
  geom_boxplot(aes(guild,(pChange-1)*100,fill=type),outlier.shape=NA,size=1.5)+
  geom_jitter(aes(guild,(pChange-1)*100,fill=type),shape=21,position=position_jitterdodge(0.15),size=5)+
  scale_x_discrete(drop=F,name="Adult Feeding Guild")+
  scale_fill_viridis_d(option="C")+
  scale_y_continuous(name="Percent Change in Emptiness",breaks=c(-100,-50,0,50,100,150,200),
                     labels=c("-100","-50%","0%","50%","100%","150%","200%"))+
  theme(legend.position=c(0.15,0.8),legend.background = element_rect(color="black"))

summary(db.guild<-aov(pChange~finalguild,data=pChange.drivers.db))
summary(rc.guild<-aov(pChange~finalguild,data=pChange.drivers.rc))
TukeyHSD(rc.guild)
summary(pe.guild<-aov(pChange~finalguild,data=pChange.drivers.pe))
TukeyHSD(pe.guild)
bind_rows(pChange.drivers.db,pChange.drivers.pe,pChange.drivers.rc)%>%
  mutate(finalguild=fct_reorder(finalguild,pChange))%>%
  ggplot()+
  geom_hline(aes(yintercept=0),lty=2,size=2)+
  geom_boxplot(aes(type,(pChange-1)*100,fill=finalguild),outlier.shape=NA,size=1.5,alpha=0.7)+
  geom_jitter(aes(type,(pChange-1)*100,fill=finalguild),shape=21,position=position_jitterdodge(0.15),size=5)+
  geom_text(aes(x=2,y=225,label="***"),size=9,inherit.aes=F)+
  geom_text(data=data.frame(x=c(1.7,1.85,2,2.15,2.3),y=c(210,210,210,210,210),label=c("ab","b","a","a","ab")),
            aes(x,y,label=label),size=7)+
  geom_text(aes(x=3,y=80,label="***"),size=9,inherit.aes=F)+
  geom_text(data=data.frame(x=c(2.7,2.85,3,3.15,3.3),y=c(65,65,65,65,65),label=c("ab","a","abc","bc","c")),
            aes(x,y,label=label),inherit.aes=F,size=7)+
  scale_x_discrete(drop=F,name="Feeding Response",labels=c("Diet\nBreadth","Percent Empty\nStomachs","Relative\nConsumption"))+
  scale_fill_viridis_d(option="C",name="Adult Feeding Guild")+
  scale_y_continuous(name="Percent Change in Response",breaks=c(-100,-50,0,50,100,150,200),
                     labels=c("-100","-50%","0%","50%","100%","150%","200%"))+
  theme(legend.position=c(0.2,0.85),legend.background = element_rect(color="black"))


#Loaded GLM approach (with dredge)

hareAttributes<-c("Adult Mobility","Complexity in Reproductive Strategy","Dispersal of Early Life Stages", 
                  "Habitat Specificity","Population Growth Rate","Prey Specificity","Sensitivity to Temperature",
                  "Spawning Cycle","Stock Size/Status")


testingClusts_Hare<-left_join(data.frame(Species=arrange(sp.responseCats.pe,Species,season)$Species,
                                         season=arrange(sp.responseCats.pe,Species,season)$season,
                                         PE=arrange(sp.responseCats.pe,Species,season)$pChange,
                                         DB=arrange(sp.responseCats.db,Species,season)$pChange,
                                         RC=arrange(sp.responseCats.rc,Species,season)$pChange),
                              finalClusts)%>%
  left_join(hare_qualScores%>%filter(Attribute%in%hareAttributes)%>%
              pivot_wider(id_cols=Species,names_from=Attribute,values_from=meanScore))
testingClusts_Hare<-filter(testingClusts_Hare,!is.na(`Adult Mobility`))

hist(testingClusts_Hare$PE)
hist(testingClusts_Hare$RC)
hist(testingClusts_Hare$DB)

options(na.action="na.fail")
summary(ph.lmer.pe<-lm(PE~.-Species-season-DB-RC,data=testingClusts_Hare))
summary(ph.lmer.rc<-lm(RC~.-Species-season-DB-PE,data=testingClusts_Hare))
summary(ph.lmer.db<-lm(DB~.-Species-season-PE-RC,data=testingClusts_Hare))

dredge(ph.lmer.pe)
dredge(ph.lmer.rc)
dredge(ph.lmer.db)


##### Multivariate approaches #####
metricMat<-as.matrix(data.frame(PE=arrange(sp.responseCats.pe,Species,season)$pChange,
                                DB=arrange(sp.responseCats.db,Species,season)$pChange,
                                RC=arrange(sp.responseCats.rc,Species,season)$pChange))
dataMat<-as.matrix(arrange(left_join(indResponses.good%>%ungroup()%>%dplyr::select(Species,season)%>%distinct(),
                                     hare_qualScores%>%filter(Attribute%in%hareAttributes)%>%
                                                      pivot_wider(id_cols=Species,names_from=Attribute,values_from=meanScore)),
                           Species,season))
metricMat<-metricMat[-which(is.na(dataMat[,3])),]
dataMat<-dataMat[-which(is.na(dataMat[,3])),]
#fullMat<-cbind(dataMat,metricMat)

dim(metricMat)
dim(dataMat)
spMat<-dataMat[,1:2]
dataMat<-dataMat[,3:ncol(dataMat)]
dataMat[,1:length(hareAttributes)]<-apply(dataMat[,1:length(hareAttributes)], 2, as.numeric)
dim(dataMat)

#### PLSR ####
library(pls)

pls_fit.pe<-plsr(metricMat[,1]~dataMat,scale=TRUE,validation="CV") #Doesn't work with categorical variables
summary(pls_fit.pe)
pls_fit.db<-plsr(metricMat[,2]~dataMat,scale=TRUE,validation="CV")
summary(pls_fit.db)
pls_fit.rc<-plsr(metricMat[,3]~dataMat,scale=TRUE,validation="CV")
summary(pls_fit.rc)


#### CCA (Attempting to see how a matrix of descriptive variables correlate with the matrix of responses in each "site") ####
  #So I want to see how species features correlate with their feeding responses. Need features to be continuous
library(CCA)

#Visualize correlations
plot(as.data.frame(metricMat))
####The matcor function looks for all possible correlations both within variablies in the first and second matrix, as well as between the two matrices
matcor(metricMat, dataMat)


###here is the cca, I put in here scale=TRUE for the data to be standardized
cc1 <- cca(metricMat,dataMat)
###quick plot of cca 
plot(cc1)

###summary will give all the scores, eigenvalues, and other output from the cca
summary(cc1)

#Looking at the cc1 data
envScore<-as.data.frame(cc1$CCA$biplot) #Env. characteristics
spScore<-cc1$CCA$wa #Species scores, where every species goes based on its combo of feeding metrics
metScore<-as.data.frame(cc1$CCA$v) #Scores in feeding metrics

spScore<-as.data.frame(cbind(spMat,spScore))
spScore[,3:4]<-apply(spScore[,3:4],2,as.numeric)

ggplot()+
  geom_text(data=spScore,aes(CCA1,CCA2,label=Species,color=season))+
  geom_text(data=metScore,aes(CCA1,CCA2,label=rownames(metScore)))+
  geom_text(data=envScore,aes(CCA1,CCA2,label=rownames(envScore)),color="blue",size=6)


#### PerMANOVA ####
testingClusts<-left_join(data.frame(Species=arrange(sp.responseCats.pe,Species,season)$Species,
                                    season=arrange(sp.responseCats.pe,Species,season)$season,
                                    PE=arrange(sp.responseCats.pe,Species,season)$pChange,
                                    DB=arrange(sp.responseCats.db,Species,season)$pChange,
                                    RC=arrange(sp.responseCats.rc,Species,season)$pChange),
                         finalClusts)
adonis(cbind(PE,DB,RC)~finalguild,data=testingClusts,permutations=10000,method="bray")

testingClusts_Hare<-left_join(data.frame(Species=arrange(sp.responseCats.pe,Species,season)$Species,
                                    season=arrange(sp.responseCats.pe,Species,season)$season,
                                    PE=arrange(sp.responseCats.pe,Species,season)$pChange,
                                    DB=arrange(sp.responseCats.db,Species,season)$pChange,
                                    RC=arrange(sp.responseCats.rc,Species,season)$pChange),
                         finalClusts)%>%
  left_join(hare_qualScores%>%filter(Attribute%in%hareAttributes)%>%
              pivot_wider(id_cols=Species,names_from=Attribute,values_from=meanScore))
testingClusts_Hare<-filter(testingClusts_Hare,!is.na(`Adult Mobility`))
adonis(cbind(PE,DB,RC)~finalguild+`Adult Mobility`,data=testingClusts_Hare,permutations=10000,method="bray")



# Plotting Summarized GLM results for Species (Pretty Tables) -------------

#Percent Empty
sp.responseCats.pe<-sp.responseCats.pe%>%
  arrange(responseCat2,Species,season)%>%
  group_by(responseCat2)%>%
  mutate(N=n(),y=N)

for (i in 2:nrow(sp.responseCats.pe)) {
  sp.responseCats.pe$y[i]<-ifelse(sp.responseCats.pe$responseCat2[i]==sp.responseCats.pe$responseCat2[i-1],
                                  sp.responseCats.pe$y[i-1]-1,
                                  sp.responseCats.pe$y[i])
  sp.responseCats.pe$N[i]<-ifelse(sp.responseCats.pe$responseCat2[i]==sp.responseCats.pe$responseCat2[i-1],
                                  0,
                                  sp.responseCats.pe$N[i])
}


ggplot(sp.responseCats.pe)+
  geom_col(aes(responseCat2,N,fill=responseCat2),color="black")+
  geom_text(aes(x=responseCat2,y-0.33,label=paste0(Species," (",season,")"),color=season),fontface="bold",angle=270)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T)[2:5],guide="none")+
  scale_y_continuous(expand=expansion(add=c(0,2)),limits=c(0,50))+
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  xlab("Response in Percent Empty Stomachs")+
  coord_flip()
#Simplifieid for presenttion
ggplot(sp.responseCats.pe)+
  geom_col(aes(responseCat2,N,fill=responseCat2),color="black")+
  geom_text(aes(x=responseCat2,y-0.33,label=Species,color=season),fontface="bold",angle=270,size=5.5)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T)[2:5],guide="none")+
  scale_y_continuous(expand=expansion(add=c(0,2)),limits=c(0,50))+
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  xlab("Response in Percent Empty Stomachs")+
  coord_flip()


#Diet Breadth
sp.responseCats.db<-sp.responseCats.db%>%
  arrange(responseCat2,Species,season)%>%
  group_by(responseCat2)%>%
  mutate(N=n(),y=N)

for (i in 2:nrow(sp.responseCats.db)) {
  sp.responseCats.db$y[i]<-ifelse(sp.responseCats.db$responseCat2[i]==sp.responseCats.db$responseCat2[i-1],
                                                       sp.responseCats.db$y[i-1]-1,
                                                       sp.responseCats.db$y[i])
  sp.responseCats.db$N[i]<-ifelse(sp.responseCats.db$responseCat2[i]==sp.responseCats.db$responseCat2[i-1],
                                  0,
                                  sp.responseCats.db$N[i])
}



ggplot(sp.responseCats.db)+
  geom_col(aes(responseCat2,N,fill=responseCat2),color="black")+
  geom_text(aes(x=responseCat2,y-0.33,label=paste0(Species," (",season,")"),color=season),fontface="bold",angle=270)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T),guide="none")+
  scale_y_continuous(expand=expansion(add=c(0,2)),limits=c(0,50))+
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  xlab("Response in Diet Breadth")+
  coord_flip()
#Simplified for presentation
ggplot(sp.responseCats.db)+
  geom_col(aes(responseCat2,N,fill=responseCat2),color="black")+
  geom_text(aes(x=responseCat2,y-0.33,label=Species,color=season),fontface="bold",angle=270,size=5.5)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T),guide="none")+
  scale_y_continuous(expand=expansion(add=c(0,2)),limits=c(0,50))+
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  xlab("Response in Diet Breadth")+
  coord_flip()


#Relative Consumption
sp.responseCats.rc<-sp.responseCats.rc%>%
  arrange(responseCat2,Species,season)%>%
  group_by(responseCat2)%>%
  mutate(N=n(),y=N)

for (i in 2:nrow(sp.responseCats.rc)) {
  sp.responseCats.rc$y[i]<-ifelse(sp.responseCats.rc$responseCat2[i]==sp.responseCats.rc$responseCat2[i-1],
                                  sp.responseCats.rc$y[i-1]-1,
                                  sp.responseCats.rc$y[i])
  sp.responseCats.rc$N[i]<-ifelse(sp.responseCats.rc$responseCat2[i]==sp.responseCats.rc$responseCat2[i-1],
                                  0,
                                  sp.responseCats.rc$N[i])
}



ggplot(sp.responseCats.rc)+
  geom_col(aes(responseCat2,N,fill=responseCat2),color="black")+
  geom_text(aes(x=responseCat2,y-0.33,label=paste0(Species," (",season,")"),color=season),fontface="bold",angle=270)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T),guide="none")+
  scale_y_continuous(expand=expansion(add=c(0,2)),limits=c(0,50))+
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  xlab("Response in Relative Consumption")+
  coord_flip()
#Simplified, maybe easier for a presentation
ggplot(sp.responseCats.rc)+
  geom_col(aes(responseCat2,N,fill=responseCat2),color="black")+
  geom_text(aes(x=responseCat2,y-0.33,label=Species,color=season),fontface="bold",angle=270,size=5.5)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T),guide="none")+
  scale_y_continuous(expand=expansion(add=c(0,2)),limits=c(0,50))+
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  xlab("Response in Relative Consumption")+
  coord_flip()


#Proportions "Bad" over time
bad.rc<-sp.responseCats.rc%>%filter(grepl("Negative",responseCat2))
nrow(bad.rc)/nrow(sp.responseCats.rc)
propBad.rc<-n_distinct(bad.rc$Species)/n_distinct(sp.responseCats.rc$Species)
bad.db<-sp.responseCats.db%>%filter(grepl("Negative",responseCat2))
nrow(bad.db)/nrow(sp.responseCats.db)
propBad.db<-n_distinct(bad.db$Species)/n_distinct(sp.responseCats.db$Species)
bad.pe<-sp.responseCats.pe%>%filter(grepl("Positive",responseCat2))
nrow(bad.pe)/nrow(sp.responseCats.pe)
propBad.pe<-n_distinct(bad.pe$Species)/n_distinct(sp.responseCats.pe$Species)

allBad<-bind_rows(bad.rc%>%mutate(type="RC"),bad.db%>%mutate(type="DB"),bad.pe%>%mutate(type="PE"))%>%
  mutate(responseCat2=gsub(" Positive","",responseCat2),
         responseCat2=gsub(" Negative","",responseCat2))%>%
  pivot_wider(id_cols=c("Species","season"),names_from="type",values_from=responseCat2)%>%
  mutate(RC.num=ifelse(is.na(RC),0,ifelse(RC=="Strong",2,1)),
         DB.num=ifelse(is.na(DB),0,ifelse(DB=="Strong",2,1)),
         PE.num=ifelse(is.na(PE),0,ifelse(PE=="Strong",2,1)),
         badSeason=RC.num+DB.num+PE.num)%>%
  group_by(Species)%>%
  mutate(badTotal=sum(badSeason))
n_distinct(allBad$Species)
filter(sp.responseCats.db,Species%notin%allBad$Species)$Species
n_distinct(totalBad$Species)


#Total scores: strong bad =-2, weak bad =-1, neutral =0, weak good =+1, strong good =+2
sp.responseCats.db<-sp.responseCats.db%>%
  mutate(direction=ifelse(grepl("Negative",responseCat2),-1,1),
         magnitude=ifelse(grepl("Strong",responseCat2),2,ifelse(grepl("Weak",responseCat2),1,0)),
         score=direction*magnitude,
         type="DB")
sp.responseCats.rc<-sp.responseCats.rc%>%
  mutate(direction=ifelse(grepl("Negative",responseCat2),-1,1),
         magnitude=ifelse(grepl("Strong",responseCat2),2,ifelse(grepl("Weak",responseCat2),1,0)),
         score=direction*magnitude,
         type="RC")
sp.responseCats.pe<-sp.responseCats.pe%>%
  mutate(direction=ifelse(grepl("Positive",responseCat2),-1,1),
         magnitude=ifelse(grepl("Strong",responseCat2),2,ifelse(grepl("Weak",responseCat2),1,0)),
         score=direction*magnitude,
         type="PE")
sp.responseScores<-bind_rows(sp.responseCats.db,sp.responseCats.rc,sp.responseCats.pe)%>%
  pivot_wider(id_cols=c("Species","season"),names_from="type",values_from="score")%>%
  mutate(seasonScore=RC+PE+DB)%>%
  group_by(Species)%>%
  mutate(totalScore=sum(seasonScore),
         meanscore=mean(seasonScore))
totalBad<-filter(sp.responseScores,RC<0 & DB<0 & PE<0)


cor(sp.responseScores$RC,sp.responseScores$PE)

#Blank plot for presentation

ggplot(sp.responseCats.rc,aes(responseCat2,N,fill=responseCat2),color="black")+
  #geom_col()+
  #geom_text(aes(x=responseCat2,y-0.33,label=paste0(Species," (",season,")"),color=season),fontface="bold")+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T),guide="none")+
  scale_y_continuous(expand=expansion(add=c(0,2)),limits=c(0,50))+
  scale_x_discrete(drop=F)+
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  xlab("Response in [   ]")+
  coord_flip()


#Example neutral relationship
ggplot()+
  geom_point(aes(x=rep(1:100),y=rnorm(100,sd=0.1)))+
  geom_segment(aes(x=1,xend=100,y=0,yend=0),size=2)+
  ylim(-1,1)+
  theme_void()

ggplot()+
  geom_point(aes(x=seq(1:100),y=1.5*seq(1:100)+rnorm(100,mean=0,sd=50)))+
  geom_segment(aes(x=1,xend=100,y=1,yend=150),size=2)+
  ylim(-200,200)+
  theme_void()


#To have them all on one plot--kinda ugly
sp.responseCats.all<-bind_rows(sp.responseCats.db%>%mutate(type="Diet Breadth")%>%group_by(type,responseCat2,.drop=F)%>%summarise(N=max(N)),
                               sp.responseCats.rc%>%mutate(type="Relative Consumption")%>%group_by(type,responseCat2,.drop=F)%>%summarise(N=max(N)))%>%
  bind_rows(sp.responseCats.pe%>%mutate(type="Percent Empty")%>%group_by(type,responseCat2,.drop=F)%>%summarise(N=max(N)))%>%
  mutate(N=ifelse(N<0,0,N))

ggplot(sp.responseCats.all)+
  geom_col(aes(responseCat2,N,fill=responseCat2,group=type,color=type),position="dodge",size=1)+
  scale_fill_manual(values=colorspace::diverge_hcl(n=5,palette="Berlin",l=c(20,50),rev=T),guide="none")+
  scale_color_manual(values=c("red","blue","green"))




# Confounding Impacts ---------------------------------------------

#Temperature should raise %Empty
#DOY is a complex relationship, but if it hits at spawning times for example, then it might raise %Empty
testingPE<-glm(cbind(nEmpty,nDiets)~DOY+surftemp+Species:DOY+Species:surftemp,
               data=indResponses,family=binomial())
summary(testingPE)

#But first we need to know how the DOY has changed with year
#Have the survey always taken the same amount of time? Yes
indResponses%>%
  group_by(year,season)%>%
  summarise(`Duration of Survey (d)`=max(DOY)-min(DOY))%>%
  ggplot(aes(year,`Duration of Survey (d)`,fill=season,color=season))+
  geom_col(color="black",position="dodge")+
  geom_smooth(method="lm",size=2)

testingDOY<-lm(DOY~year*season,data=indResponses)
testingDOY.pred<-data.frame(indResponses,pred=predict(testingDOY))
testingDOY.sum<-testingDOY.pred%>%
  group_by(season)%>%
  summarise(diff=max(pred)-min(pred),
            labelY=210-diff*2)
ggplot()+
  geom_boxplot(data=testingDOY.pred,aes(year,DOY,group=year),size=1.3,outlier.shape=NA)+
  geom_point(data=testingDOY.pred,aes(year,DOY,fill=geoarea),shape=21,size=3,alpha=0.5)+
  geom_line(data=testingDOY.pred,aes(year,pred),size=2,color="black")+
  geom_text(data=testingDOY.sum,aes(2001,labelY,label=paste0("Survey Days Later=",round(diff,1))),size=8,hjust=0)+
  scale_fill_manual(values=colorspace::qualitative_hcl(5,palette="Harmonic"),name="Geographic\nArea")+
  ylab("Day of Year")+xlab("Year")+
  theme(legend.position=c(0.83,0.28),legend.direction="horizontal",
        legend.background=element_rect(fill="white",color="black"),
        strip.text=element_text(size=30))+
  guides(fill=guide_legend(override.aes=list(size=10,alpha=1),nrow=2,byrow=T))+
  facet_wrap(~season)

#Looking at spawning times
spawnTimes<-read_csv("spawningTimes.csv",col_types = list(Spawn_Start=col_date(format="%m/%d/%Y"),
                                                          Spawn_End=col_date(format="%m/%d/%Y")))%>%
  mutate(Spawn_Start_doy=yday(Spawn_Start),
         Spawn_End_doy=yday(Spawn_End),
         plot_start1=ifelse(Spawn_Start_doy<Spawn_End_doy,Spawn_Start_doy,1),
         plot_end1=ifelse(Spawn_Start_doy<Spawn_End_doy,Spawn_End_doy,Spawn_End_doy),
         plot_start2=ifelse(Spawn_Start_doy<Spawn_End_doy,NA,Spawn_Start_doy),
         plot_end2=ifelse(Spawn_Start_doy<Spawn_End_doy,NA,366))

ggplot(spawnTimes)+
  geom_rect(aes(xmin=plot_start1,xmax=plot_end1,ymin=0,ymax=1),fill="grey50",alpha=0.5)+
  geom_rect(aes(xmin=plot_start2,xmax=plot_end2,ymin=0,ymax=1),fill="grey50",alpha=0.5)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  scale_x_continuous(expand=expansion(0))+
  scale_y_continuous(expand=expansion(0))+
  facet_wrap(~Species)

#Overlapping spawn times with trawl times
ggplot()+
  geom_bar(data=uniqueDiets%>%mutate(Species=str_to_title(pdcomnam)),aes(x=yday(ymd(paste(year,month,day,sep="-")))),fill="black")+
  geom_rect(data=spawnTimes,aes(xmin=plot_start1,xmax=plot_end1,ymin=0,ymax=Inf),fill="grey50",alpha=0.5)+
  geom_rect(data=spawnTimes,aes(xmin=plot_start2,xmax=plot_end2,ymin=0,ymax=Inf),fill="grey50",alpha=0.5)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  scale_x_continuous(expand=expansion(0),name="Day of Year")+
  scale_y_continuous(expand=expansion(0),name="Number of Diets")+
  facet_wrap(~Species,scales="free_y")

#How are PIK and RC values in and out of spawning season
ggplot()+
  geom_col(data=indResponses%>%group_by(Species,DOY)%>%
             summarise(nDiets=sum(nDiets,na.rm=T),mean=sum(nEmpty,na.rm=T)/sum(nDiets,na.rm=T),sd=sd(pik,na.rm=T))%>%
             filter(nDiets>10),aes(x=DOY,y=mean),fill="black")+
  geom_rect(data=spawnTimes,aes(xmin=plot_start1,xmax=plot_end1,ymin=0,ymax=Inf),fill="grey50",alpha=0.5)+
  geom_rect(data=spawnTimes,aes(xmin=plot_start2,xmax=plot_end2,ymin=0,ymax=Inf),fill="grey50",alpha=0.5)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  scale_x_continuous(expand=expansion(0),name="Day of Year")+
  scale_y_continuous(expand=expansion(0),name="Mean Percent of Empty Stomachs")+
  facet_wrap(~Species)

ggplot()+
  geom_col(data=indResponses%>%group_by(Species,DOY)%>%
             summarise(tDiets=sum(nDiets,na.rm=T),mean=weighted.mean(meanC,nDiets,na.rm=T),sd=sd(meanC,na.rm=T))%>%
             filter(tDiets>10),aes(x=DOY,y=mean),fill="black")+
  geom_rect(data=spawnTimes,aes(xmin=plot_start1,xmax=plot_end1,ymin=0,ymax=Inf),fill="grey50",alpha=0.5)+
  geom_rect(data=spawnTimes,aes(xmin=plot_start2,xmax=plot_end2,ymin=0,ymax=Inf),fill="grey50",alpha=0.5)+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  scale_x_continuous(expand=expansion(0),name="Day of Year")+
  scale_y_continuous(expand=expansion(0),name="Mean Relative Consumption")+
  facet_wrap(~Species,scales="free_y")




#And also temp with year
testingTEMP<-lm(surftemp~year*season*geoarea,data=indResponses)
summary(testingTEMP)
testingTEMP.pred<-data.frame(filter(indResponses,!is.na(surftemp)),pred=predict(testingTEMP))
testingTEMP.sum<-testingTEMP.pred%>%
  group_by(season,geoarea)%>%
  summarise(diff=max(pred)-min(pred),
            nTrawls=n_distinct(id))%>%
  group_by(season)%>%
  summarise(diff=weighted.mean(diff,nTrawls),
            labelY=51-diff^4)
ggplot()+
  geom_boxplot(data=testingTEMP.pred,aes(year,surftemp,group=year),size=1.3,outlier.shape=NA)+
  geom_point(data=testingTEMP.pred,aes(year,surftemp,fill=geoarea),shape=21,size=3,alpha=0.5)+
  geom_line(data=testingTEMP.pred,aes(year,pred,group=geoarea,color=geoarea),size=2,show.legend=F)+
  geom_text(data=testingTEMP.sum,aes(2002,labelY,label=paste0("Temperature Rise=",round(diff,1))),size=8,hjust=0)+
  scale_fill_manual(values=colorspace::qualitative_hcl(5,palette="Harmonic"),name="Geographic\nArea")+
  scale_color_manual(values=colorspace::qualitative_hcl(5,palette="Harmonic"))+
  ylab("Surface Temperature")+xlab("Year")+
  theme(legend.position=c(0.16,0.82),legend.direction="horizontal",
        legend.background=element_rect(fill="white",color="black"),
        strip.text=element_text(size=30))+
  guides(fill=guide_legend(override.aes=list(size=10,alpha=1),nrow=2,byrow=T))+
  facet_wrap(~season)


#Time of day as an effect of the amount in stomachs/pik
indResponses%>%
  group_by(TOD,Species)%>%
  summarise(PIK=sum(nEmpty,na.rm=T)/sum(nDiets,na.rm=T),
            relC=weighted.mean(meanC,abundance,na.rm=T))%>%
  ggplot()+
  geom_col(aes(x=TOD,y=PIK))+
  coord_polar()+
  facet_wrap(~Species)
indResponses%>%
  group_by(TOD,Species)%>%
  summarise(PIK=sum(nEmpty,na.rm=T)/sum(nDiets,na.rm=T),
            relC=weighted.mean(meanC,abundance,na.rm=T))%>%
  ggplot()+
  geom_col(aes(x=TOD,y=relC))+
  scale_y_sqrt()+
  coord_polar()+
  facet_wrap(~Species)



#If the exponent for the gastric evacuation temperature coeffieicnt is 0.115 as recommended by Durbin et al (1983)
#The equation is R = ae^bT where a would be a prey-species constant
#Don't know what a should be, but it can be a constant
#For this, should use the bottemp estimates, even if they are a bit more limited
testingBOTTEMP<-lm(bottemp~year*season,data=indResponses)
summary(testingBOTTEMP)
testingBOTTEMP.pred<-data.frame(filter(indResponses,!is.na(bottemp)),pred=predict(testingBOTTEMP))
testingBOTTEMP.sum<-testingBOTTEMP.pred%>%
  group_by(season)%>%
  summarise(diff=max(pred)-min(pred),
            labelY=51-diff^7.5)

ggplot()+
  geom_boxplot(data=testingBOTTEMP.pred,aes(year,bottemp,group=year),size=1.3,outlier.shape=NA)+
  geom_point(data=testingBOTTEMP.pred,aes(year,bottemp,fill=geoarea),shape=21,size=3,alpha=0.5)+
  geom_line(data=testingBOTTEMP.pred,aes(year,pred),color="black",size=2)+
  geom_text(data=testingBOTTEMP.sum,aes(2002,labelY,label=paste0("Temperature Rise=",round(diff,2))),size=8,hjust=0)+
  scale_fill_manual(values=colorspace::qualitative_hcl(5,palette="Harmonic"),name="Geographic\nArea")+
  scale_color_manual(values=colorspace::qualitative_hcl(5,palette="Harmonic"))+
  ylab("Bottom Temperature")+xlab("Year")+
  theme(legend.position=c(0.16,0.87),legend.direction="horizontal",
        legend.background=element_rect(fill="white",color="black"),
        strip.text=element_text(size=30))+
  guides(fill=guide_legend(override.aes=list(size=10,alpha=1),nrow=2,byrow=T))+
  facet_wrap(~season)


###### Seeing the final "answer" of temp rises resulting in faster gastric evac ######
testingBOTTEMP.detailSum<-testingBOTTEMP.pred%>%
  group_by(season,year)%>%
  summarise(Temp=mean(pred))
gastricEvac.temp<-testingBOTTEMP.detailSum%>%
  mutate(R.temp=exp(0.111*Temp))
gastricEvac.temp.sum<-gastricEvac.temp%>%
  group_by(season)%>%
  summarise(diff=max(R.temp)-min(R.temp),
            pChange=diff/mean(R.temp)*100,
            pChange2=max(R.temp)/min(R.temp))
mean(gastricEvac.temp.sum$pChange2)
1-mean(sp.responseCats.rc$pChange) 
nrow(filter(sp.responseCats.rc,pChange-1<(1-mean(gastricEvac.temp.sum$pChange2))))/nrow(sp.responseCats.rc) #Proportion when the %decline in rc is more than %decline from ge

ggplot()+
  geom_line(data=gastricEvac.temp,aes(Temp,R.temp,color=season),size=2)+
  geom_text(aes(x=11,y=2.1,label="*Assuming universal b=0.111\n(Durbin et al. 1983)"),hjust=1)+
  scale_x_continuous(name="Temperature (°C)")+
  scale_y_continuous(name="Relative Gastric Evacuation Rate*")+
  scale_color_viridis_d(option="D",begin=0.1,end=0.98,name="Season When\nTemperatures were Observed")+
  theme(legend.position=c(0.3,0.7),legend.background=element_rect(color="black"))


(3.5*0.01*24*0.01)/((3.5*mean(gastricEvac.temp.sum$pChange2))*0.01*24*(0.01/(1+mean(abs(sp.responseCats.rc$pChange-1)))))
#The initial here is 11% higher DR than the later
#Note that it DOES NOT MATTER what the initial R or S might have been for a species, it's relative changes in both

#What does that mean for daily rations over time
#Important Assumptions!!!!
#R is same for all species/prey, especially the a coefficient
#Body Size and Temperature are both averaged out
#Assuming that the b coefficient is constant across species
#Assuming that the temperature increase recorded by trawls reflects temperature change experienced by any and all fish



#But how many more stomachs would be empty based on temperature alone???
#Trying a simulation study
a<-c(0.002,0.004,0.011,0.028,0.041)
pik<-data.frame(Temp=filter(testingBOTTEMP.detailSum,season=="Fall")$Temp)
temp<-numeric(length=indResponses%>%group_by(Species)%>%summarise(N=n_distinct(id))%>%ungroup()%>%summarise(N=round(mean(N),0))%>%as.numeric())

for (j in 1:length(a)) {
  for (i in 1:nrow(filter(testingBOTTEMP.detailSum,season=="Fall"))) {
    for (k in 1:length(temp)) {
      b=0.111
      Temp=filter(testingBOTTEMP.detailSum,season=="Fall")$Temp
      R=a[j]*exp(b*Temp)
      
      time=seq(0,2000)
      Wp=exp(-1*R[1]*time)
      tEvac=min(which(Wp<0.01))
      Wp<-exp(-1*R[i]*time[1:(1*tEvac)])
      diets=sample(Wp,indResponses%>%group_by(id)%>%summarise(N=mean(sum(nDiets)))%>%ungroup()%>%summarise(N=round(mean(N),0))%>%as.numeric(),replace=T)
      temp[k]=sum(diets<0.01)/length(diets)
      pik[i,j+1]=mean(temp)
      colnames(pik)<-c(colnames(pik)[1:j],paste("a",a[j],sep="_"))
    }
  }
}

#Trying a simpler simulation, since the above takes a long time
a<-c(0.002,0.004,0.011,0.028,0.041) #Series from Stehlik et al. 2021
pik<-data.frame(Temp=filter(testingBOTTEMP.detailSum,season=="Fall")$Temp,
                startingMult=sort(rep(seq(1,2,by=0.1),47)))
#what should count as "empty" in simulations?
minRC<-filter(uniqueDiets, pdwgt>0)%>%
  mutate(relConsump=pdgutw/pdwgt)%>%ungroup()%>%filter(relConsump>0)%>%summarise(min(relConsump))%>%as.numeric()


for (j in 1:length(a)) {
  for (i in 1:nrow(pik)) {
    
      b=0.111
      R=a[j]*exp(b*pik[i,1])
      R1=a[j]*exp(b*pik[1,1])
    
      time=seq(0,10000)
      Wp=exp(-1*R1*time)
      tEvac=min(which(Wp<minRC))
      Wp<-exp(-1*R*time[1:(pik[i,2]*tEvac)])
      diets=sample(Wp,indResponses%>%group_by(Species,year)%>%summarise(N=mean(sum(nDiets)))%>%ungroup()%>%summarise(N=round(mean(N),0))%>%as.numeric(),replace=T)
      pik[i,j+2]=sum(diets<minRC)/length(diets)
      colnames(pik)<-c(colnames(pik)[1:(j+1)],paste("a",a[j],sep="_"))
   
  }
}
1-minRC

#Visual to explain the process
R1.plot=a[1]*exp(b*min(filter(testingBOTTEMP.detailSum,season=="Fall")$Temp))
R2.plot=a[1]*exp(b*max(filter(testingBOTTEMP.detailSum,season=="Fall")$Temp))

Wp1.plot=exp(-1*R1.plot*time)
Wp2.plot=exp(-1*R2.plot*time)
tEvac=min(which(Wp1.plot<0.01))
tEvac2=min(which(Wp2.plot<0.01))
difftEvac<-tEvac-tEvac2
Wp1.plot=Wp1.plot[1:(2*tEvac)]
Wp2.plot=Wp2.plot[1:(2*tEvac)]
Wp.plot<-data.frame(Wp=c(Wp1.plot,Wp2.plot),
                    Temp=as.character(c(rep(round(min(filter(testingBOTTEMP.detailSum,season=="Fall")$Temp),1),tEvac*2),
                                        rep(round(max(filter(testingBOTTEMP.detailSum,season=="Fall")$Temp),1),tEvac*2))),
                    Hours=rep(seq(1,2*tEvac),2))
test<-sample(Wp2.plot[1:tEvac],256,replace=T)
sum(test<0.01)
50/256
132/782
ggplot()+
  geom_point(data=Wp.plot,aes(Hours,Wp,color=Temp))+
  geom_hline(aes(yintercept=0.01),lty=2)+
  geom_segment(aes(x=tEvac+10,xend=tEvac,y=0.2,yend=0.02),lineend="round",arrow=arrow(length=unit(0.1,"in")),size=2)+
  geom_label(aes(x=tEvac+15,y=0.22,label="tEvac"),color="firebrick2")+
  geom_segment(aes(x=tEvac*1.1,xend=tEvac*1.1,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.1,y=0.075,label="1.1"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.2,xend=tEvac*1.2,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.2,y=0.075,label="1.2"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.3,xend=tEvac*1.3,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.3,y=0.075,label="1.3"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.4,xend=tEvac*1.4,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.4,y=0.075,label="1.4"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.5,xend=tEvac*1.5,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.5,y=0.075,label="1.5"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.6,xend=tEvac*1.6,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.6,y=0.075,label="1.6"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.7,xend=tEvac*1.7,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.7,y=0.075,label="1.7"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.8,xend=tEvac*1.8,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.8,y=0.075,label="1.8"),hjust=0.2)+
  geom_segment(aes(x=tEvac*1.9,xend=tEvac*1.9,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*1.9,y=0.075,label="1.9"),hjust=0.2)+
  geom_segment(aes(x=tEvac*2,xend=tEvac*2,y=0,yend=0.05))+
  geom_text(aes(x=tEvac*2,y=0.075,label="2"),hjust=0.2)+
  theme(legend.position=c(0.8,0.7),legend.background=element_rect(color="black"))+
  geom_text(aes(x=45,y=0.45,label=paste("Time Diff to Evac =",difftEvac,"hours")),hjust=0)

pik%>%
  pivot_longer(cols=starts_with("a_"),names_to="a_Coefficient",values_to="PIK")%>%
  mutate(a_Coefficient=gsub("a_","",a_Coefficient))%>%
  ggplot()+
  geom_line(aes(x=Temp,y=PIK,color=a_Coefficient),size=2)+
  geom_smooth(aes(x=Temp,y=PIK),method="lm",size=3)+
  scale_color_viridis_d(option="B",end=0.8,name="a Coefficient")+
  scale_y_continuous(labels=scales::percent_format(),name="Empty Stomachs (%)",limits=c(0,1))+
  scale_x_continuous(name="Temperature (°C)",limits=c(9.7,11.55))+
  theme(legend.position=c(0.88,0.15),legend.background = element_rect(color="black"))+
  facet_wrap(~startingMult)

pEmpty_temp<-pik%>%
  pivot_longer(cols=starts_with("a_"),names_to="a_Coefficient",values_to="PIK")%>%
  group_by(Temp,startingMult)%>%
  summarise(meanPIK=mean(PIK))
pEmpty_temp_eachFrequency<-pEmpty_temp%>%
  ungroup()%>%
  filter(Temp==min(Temp) | Temp==max(Temp))%>%
  pivot_wider(id_cols="startingMult",names_from="Temp",values_from="meanPIK")
colnames(pEmpty_temp_eachFrequency)<-c("startingMult","lowTemp","highTemp")
pEmpty_temp_eachFrequency<-pEmpty_temp_eachFrequency%>%
  mutate(diff=highTemp-lowTemp,
         pChange=highTemp/lowTemp)
median(sp.responseCats.pe$pChange_fromTemp-1) #A typical increase across all the species
sp.responseCats.pe<-sp.responseCats.pe%>%
  rowwise()%>%
  mutate(whichMult=which(abs(Min-pEmpty_temp_eachFrequency$lowTemp)==min(abs(Min-pEmpty_temp_eachFrequency$lowTemp))),
         initial_pEmpty_fromTemp=pEmpty_temp_eachFrequency$lowTemp[whichMult],
         final_pEmpty_fromTemp=pEmpty_temp_eachFrequency$highTemp[whichMult],
         pChange_fromTemp=pEmpty_temp_eachFrequency$pChange[whichMult])
ggplot(sp.responseCats.pe)+
  geom_point(aes(Min,initial_pEmpty_fromTemp),color="grey50",alpha=0.5,size=2)+
  geom_segment(aes(x=Min,xend=Max,y=initial_pEmpty_fromTemp,yend=final_pEmpty_fromTemp),color="grey50",alpha=0.5)+
  geom_point(aes(Max,final_pEmpty_fromTemp),size=4)+
  geom_abline(slope=1,intercept=0)+
  scale_x_continuous(name="Observed %Empty",breaks=c(0,0.1,0.2,0.3,0.4,0.5),limits=c(0,0.55))+
  scale_y_continuous(name="Simulated %Empty",breaks=c(0,0.1,0.2,0.3,0.4,0.5),limits=c(0,0.55))

median(sp.responseCats.pe$pChange_fromTemp-1) #A typical increase across all the species (mean is really skewed by the 730% of smooth dogfish)
nrow(filter(sp.responseCats.pe,pChange>pChange_fromTemp))/nrow(sp.responseCats.pe)

table(sp.responseCats.pe$whichMult)

# Gastric Evacuation Rates for everyone -----------------------------------

max.col()
#Can I just map a species' relConsump's in order and get a gastric evac plot???
rcTOge<-filter(uniqueDiets, pdwgt>0)%>%
  mutate(relConsump=pdgutw/pdwgt)%>%
  arrange(pdcomnam,desc(relConsump))%>%
  group_by(pdcomnam)%>%
  mutate(x=order(relConsump,decreasing=T))%>%
  left_join(dplyr::select(preyTrawls,id,bottemp))
ggplot(rcTOge)+
  geom_point(aes(x,relConsump))+
  facet_wrap(~pdcomnam,scales="free")

#Simple exponential model for this data (starting with one species)
rcTOge.test<-filter(rcTOge,pdcomnam=="ATLANTIC COD"&relConsump>0)
summary(GEexp.lm.test<-lm(log(relConsump)~x,data=rcTOge.test))
summary(GEexp.glm.test<-glm((relConsump)~x,data=rcTOge.test,family=Gamma(link="log")))

rcTOge.test<-data.frame(rcTOge.test,pred=exp(predict(GEexp.lm.test)),pred.glm=predict.glm(GEexp.glm.test,type="response"))
ggplot(rcTOge.test)+
  geom_point(aes(x,relConsump))+
  geom_line(aes(x,pred))+
  geom_line(aes(x,pred.glm),color="blue")

#Getting a and b parameters
summary(GEexp.ab.test<-glm(pred.glm~bottemp,data=rcTOge.test,family=Gamma(link="log")))
#how does this compare to most changes in PIK?
mean(sp.responseCats.pe$pChange)
range(sp.responseCats.pe$pChange)



#Visualizing the average PIK for each DOY
indResponses%>%
  group_by(DOY,Species)%>%
  summarise(pikMean=weighted.mean(pik,abundance,na.rm=T),
            PIK=sum(nEmpty,na.rm=T)/sum(nDiets,na.rm=T))%>%
  ggplot()+
  geom_col(aes(DOY,PIK))+
  facet_wrap(~Species)




# Relating Foraging Species Responses to other Variables ------------------

envForaging<-left_join(sp.responseScores,hare_everything_FHD)%>%
  dplyr::select(Species,totalScore,meanscore,Climate.Vulnerability:Sensitivity.Attribute)%>%
  distinct()

ggplot(envForaging,aes(Climate.Direction,meanscore))+
  geom_boxplot()+
  geom_jitter(height = 0)+
  geom_text(aes(label=Species))

climateAOV<-aov(meanscore~Climate.Direction,data=envForaging)
summary(climateAOV)


#Extracting the Feeding Guilds from the GL Repeat




#More summarized using Hare et al. categories
library(lme4)
summary(percEmpty<-glm(cbind(freq,nDiets)~year:season+year:season:Climate.Vulnerability+year:season:Change.Potential+year:season:Climate.Direction,
                       data=indResponses,family=binomial(logit)))
summary(percEmpty<-glm(cbind(freq,nDiets)~year*season*Climate.Vulnerability, data=indResponses,family=binomial(logit)))

plotResponses.pe<-data.frame(indResponses,pred=predict(percEmpty,type="response",se=T))%>%
  mutate(se.lower=pred.fit-1.96*pred.se.fit,
         se.upper=pred.fit+1.96*pred.se.fit)%>%
  dplyr::select(year,season,Climate.Vulnerability,pred.fit,se.lower,se.upper)%>%distinct()

summary(percEmpty<-glmer(cbind(freq,nDiets)~year*season*Climate.Vulnerability+(year|Species), data=indResponses,family=binomial(link="logit")))




#Diet Breadth

library(betareg)
summary(dbreadth.beta<-betareg(simpson2~year*season*Climate.Vulnerability,data=indResponses))
summary(dbreadth.glm<-glm(simpson~year*season*Climate.Vulnerability,data=indResponses,family=quasibinomial(logit)))

plotResponses.db<-data.frame(filter(indResponses,!is.na(simpson)),pred=predict(dbreadth.glm,type="response",se=T))%>%
  mutate(se.lower=pred.fit-1.96*pred.se.fit,
         se.upper=pred.fit+1.96*pred.se.fit)%>%
  dplyr::select(year,season,Climate.Vulnerability,pred.fit,se.lower,se.upper)%>%distinct()



summary(relConsump<-glm(meanC~year:season+year:season:Climate.Vulnerability+
                          year:season:Change.Potential+year:season:Climate.Direction,data=indResponses,family=Gamma(link="log")))
summary(relConsump<-glm(meanC~year*season*Climate.Vulnerability,data=indResponses,family=Gamma(link="log")))

plotResponses.rc<-data.frame(filter(indResponses,!is.na(meanC)),pred=predict(relConsump,type="response",se=T))%>%
  mutate(se.lower=pred.fit-1.96*pred.se.fit,
         se.upper=pred.fit+1.96*pred.se.fit)%>%
  dplyr::select(year,season,Climate.Vulnerability,pred.fit,se.lower,se.upper)%>%distinct()


#Cheap-o versions for rushing the abstract
summary(forAbstract<-glm(meanC~year*Species,data=indResponses,family=Gamma(link="log")))
pred<-data.frame(filter(indResponses,!is.na(meanC)),pred=predict.glm(forAbstract,type="response"))%>%
  filter(year=="1992" | year=="2019")%>%
  pivot_wider(id_cols=Species,names_from=year,names_prefix = "year",values_from=pred,values_fn=mean)%>%
  mutate(diff=year1992-year2019)


summary(forAbstract2<-glm(cbind(freq,nDiets)~year*Species,data=indResponses,family=binomial(logit)))
pred2<-data.frame(indResponses,pred=predict.glm(forAbstract2,type="response"))%>%
  filter(year=="1992" | year=="2019")%>%
  pivot_wider(id_cols=Species,names_from=year,names_prefix = "year",values_from=pred,values_fn=mean)%>%
  mutate(diff=year1992-year2019)

summary(forAbstract3<-glm(simpson2~year*Species,data=indResponses,family=binomial(logit)))
pred3<-data.frame(filter(indResponses,!is.na(simpson)),pred=predict.glm(forAbstract3,type="response"))%>%
  filter(year=="1992" | year=="2019")%>%
  pivot_wider(id_cols=Species,names_from=year,names_prefix = "year",values_from=pred,values_fn=mean)%>%
  mutate(diff=year1992-year2019)





#### Plotting those GLMs  ####

#Percent Empty: summarized
indResponses%>%
  group_by(year,season,Climate.Vulnerability)%>%
  summarise(nDiets=sum(nDiets),
            freq=sum(freq),
            pik=freq/nDiets)%>%
  ggplot()+
  geom_point(aes(x=year,y=pik,size=nDiets,fill=season),shape=21)+
  facet_grid(season~Climate.Vulnerability)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Proportion Empty Stomachs")
#...Raw data
ggplot()+
  geom_point(data=indResponses,aes(x=year,y=pik,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=plotResponses.pe,aes(year,ymin=se.lower,ymax=se.upper,fill=season),alpha=0.5)+
  geom_line(data=plotResponses.pe,aes(year,pred.fit,color=season),size=2)+
  facet_grid(season~Climate.Vulnerability)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Proportion Empty Stomachs")
#...Individually for species
ggplot()+
  geom_point(data=indResponses.good,aes(x=year,y=pik,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=sp.plotResponses.pe,aes(year,ymin=se.lower,ymax=se.upper,fill=season,color=season),alpha=0.5)+
  geom_line(data=sp.plotResponses.pe,aes(year,pred.fit,color=season),size=2)+
  facet_wrap(~Species,nrow=4)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Proportion Empty Stomachs")
#...Individually Summarized
indResponses.good%>%
  group_by(Species,year,season)%>%
  summarise(nDiets=sum(nDiets),
            freq=sum(freq),
            pik=freq/nDiets)%>%
  ggplot()+
  geom_point(aes(x=year,y=pik,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=sp.plotResponses.pe,aes(year,ymin=se.lower,ymax=se.upper,fill=season,color=season),alpha=0.5)+
  geom_line(data=sp.plotResponses.pe,aes(year,pred.fit,color=season),size=2)+
  facet_wrap(~Species,nrow=4)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Proportion Empty Stomachs")


#Simpson diet breadth: summarized
indResponses%>%
  group_by(year,season,Climate.Vulnerability)%>%
  summarise(breadth=mean(simpson,na.rm=T),
            nDiets=sum(nDiets))%>%
  ggplot()+
  geom_point(aes(x=year,y=breadth,size=nDiets,fill=season),shape=21)+
  facet_grid(season~Climate.Vulnerability)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Simpson's Diet Breadth",limits=c(0,1))
#...Raw data
ggplot()+
  geom_point(data=indResponses,aes(x=year,y=simpson,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=plotResponses.db,aes(year,ymin=se.lower,ymax=se.upper,fill=season),alpha=0.5)+
  geom_line(data=plotResponses.db,aes(year,pred.fit,color=season),size=2)+
  facet_grid(season~Climate.Vulnerability)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Simpson's Diet Breadth")
#...Individually for species
ggplot()+
  geom_point(data=indResponses.good,aes(x=year,y=simpson,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=sp.plotResponses.db,aes(year,ymin=se.lower,ymax=se.upper,fill=season,color=season),alpha=0.5)+
  geom_line(data=sp.plotResponses.db,aes(year,pred.fit,color=season),size=2)+
  facet_wrap(~Species,nrow=4)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Diet Breadth (Simpson's)")
#...Individually Summarized
indResponses.good%>%
  group_by(year,season,Species)%>%
  summarise(breadth=mean(simpson,na.rm=T),
            nDiets=sum(nDiets))%>%
  ggplot()+
  geom_point(aes(x=year,y=breadth,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=sp.plotResponses.db,aes(year,ymin=se.lower,ymax=se.upper,fill=season,color=season),alpha=0.5)+
  geom_line(data=sp.plotResponses.db,aes(year,pred.fit,color=season),size=2)+
  facet_wrap(~Species,nrow=4)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Diet Breadth (Simpson's)")



#Relative Consumption: summarized
indResponses%>%
  group_by(year,season,Climate.Vulnerability)%>%
  summarise(relConsump=mean(meanC,na.rm=T),
            nDiets=sum(nDiets))%>%
  ggplot()+
  geom_point(aes(x=year,y=relConsump,fill=season,size=nDiets),shape=21)+
  facet_grid(season~Climate.Vulnerability)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_continuous(name="Relative Consumption")
#...Raw data
ggplot()+
  geom_point(data=indResponses,aes(x=year,y=meanC,fill=season,size=nDiets),shape=21)+
  geom_ribbon(data=plotResponses.rc,aes(year,ymin=se.lower,ymax=se.upper,fill=season),alpha=0.5)+
  geom_line(data=plotResponses.rc,aes(year,pred.fit,color=season),size=2)+
  facet_grid(season~Climate.Vulnerability)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_log10(name="Relative Consumption")
#...Individually for species
ggplot()+
  geom_point(data=indResponses.good,aes(x=year,y=meanC,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=sp.plotResponses.rc,aes(year,ymin=se.lower,ymax=se.upper,fill=season,color=season),alpha=0.5)+
  geom_line(data=sp.plotResponses.rc,aes(year,pred.fit,color=season),size=2)+
  facet_wrap(~Species,nrow=4)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_log10(name="Relative Consumption")
#...Individually Summarized
indResponses.good%>%
  group_by(year,season,Species)%>%
  summarise(rc=mean(meanC,na.rm=T),
            nDiets=sum(nDiets))%>%
  ggplot()+
  geom_point(aes(x=year,y=rc,size=nDiets,fill=season),shape=21)+
  geom_ribbon(data=sp.plotResponses.rc,aes(year,ymin=se.lower,ymax=se.upper,fill=season,color=season),alpha=0.5)+
  geom_line(data=sp.plotResponses.rc,aes(year,pred.fit,color=season),size=2)+
  facet_wrap(~Species,nrow=4)+
  scale_color_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_fill_manual(values=seasonPal2[c(4,2)],guide="none")+
  scale_x_continuous(name="Year",limits=c(1973,2019),breaks=c(1980,2000,2020))+
  scale_y_log10(name="Relative Consumption")





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




#### Levin's Breadth ####


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


#### Relative Consumption ####

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










# Comparing the species across the decades --------------------------------

#The 1970s
props_70s<-preyTrawls%>%
  filter(year<1980
         & !is.na(GL_pyscinam)
         & GL_pyscinam!="Empty")%>%
  mutate(species=str_to_title(pdcomnam),
         decade="1970s")%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(decade,species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_70s<-unique(paste(props_70s$species,props_70s$sizecat2,sep="_"))
#The 1980s
props_80s<-preyTrawls%>%
  filter(year>=1980 & year<1990
         & !is.na(GL_pyscinam)
         & GL_pyscinam!="Empty")%>%
  mutate(species=str_to_title(pdcomnam),
         decade="1980s")%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(decade,species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_80s<-unique(paste(props_80s$species,props_80s$sizecat2,sep="_"))
#The 1990s
props_90s<-preyTrawls%>%
  filter(year>=1990 & year<2000
         & !is.na(GL_pyscinam)
         & GL_pyscinam!="Empty")%>%
  mutate(species=str_to_title(pdcomnam),
         decade="1990s")%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(decade,species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_90s<-unique(paste(props_90s$species,props_90s$sizecat2,sep="_"))
#The 2000s
props_00s<-preyTrawls%>%
  filter(year>=2000 & year<2010
         & !is.na(GL_pyscinam)
         & GL_pyscinam!="Empty")%>%
  mutate(species=str_to_title(pdcomnam),
         decade="2000s")%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(decade,species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_00s<-unique(paste(props_00s$species,props_00s$sizecat2,sep="_"))
#The 2010s
props_10s<-preyTrawls%>%
  filter(year>=2010 & year<2020
         & !is.na(GL_pyscinam)
         & GL_pyscinam!="Empty")%>%
  mutate(species=str_to_title(pdcomnam),
         decade="2010s")%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  filter(nDiets>20)%>%
  group_by(decade,species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()
species_10s<-unique(paste(props_10s$species,props_10s$sizecat2,sep="_"))

#Keep only those species_sizecats that are in all the decades
species_allDecades<-Reduce(intersect, list(species_70s,species_80s,species_90s,species_00s,species_10s))
#N groups and N species
length(species_allDecades)
length(unique(gsub("_[A-Z]+","",species_allDecades)))

#Calculate the proportions of all diet items and the amount shared across decades
list<-list(props_70s,props_80s,props_90s,props_00s,props_10s)
props_list<-list()
for (i in 1:length(list)) {
  list[[i]]<-filter(list[[i]],paste(species,sizecat2,sep="_") %in% species_allDecades)
  props_list[[i]]<-list[[i]][,-max(ncol(list[[i]]))] #Cut off the volume prop for cleanliness
  colnames(props_list[[i]])<-c("decade","species1","sizecat1","GL_pyscinam",paste0("prop_w_",unique(props_list[[i]][,"decade"])))
  props_list[[i]]<-complete(props_list[[i]],decade,species1,sizecat1,GL_pyscinam)%>%
    filter(GL_pyscinam!="Empty")%>%
    group_by(decade,species1,sizecat1)%>%mutate(t=sum(across(matches("prop_w")),na.rm=T))%>%filter(t>0)%>%dplyr::select(-t)
  props_list[[i]][is.na(props_list[[i]])]<-0
}

props_df<-bind_cols(props_list)%>%
  dplyr::select(species=species1...2,sizecat=sizecat1...3,GL_pyscinam=GL_pyscinam...4,
                matches("prop_w"))
overlap_df<-props_df%>%
  mutate(diff_1970s_1980s=abs(prop_w_1970s-prop_w_1980s),
         diff_1970s_1990s=abs(prop_w_1970s-prop_w_1990s),
         diff_1970s_2000s=abs(prop_w_1970s-prop_w_2000s),
         diff_1970s_2010s=abs(prop_w_1970s-prop_w_2010s),
         diff_1980s_1990s=abs(prop_w_1980s-prop_w_1990s),
         diff_1980s_2000s=abs(prop_w_1980s-prop_w_2000s),
         diff_1980s_2010s=abs(prop_w_1980s-prop_w_2010s),
         diff_1990s_2000s=abs(prop_w_1990s-prop_w_2000s),
         diff_1990s_2010s=abs(prop_w_1990s-prop_w_2010s),
         diff_2000s_2010s=abs(prop_w_2000s-prop_w_2010s))%>%
  group_by(species,sizecat)%>%
  summarise(ep_1970s_1980s=sum(diff_1970s_1980s), s_do_1970s_1980s=1-0.5*ep_1970s_1980s,
            ep_1970s_1990s=sum(diff_1970s_1990s), s_do_1970s_1990s=1-0.5*ep_1970s_1990s,
            ep_1970s_2000s=sum(diff_1970s_2000s), s_do_1970s_2000s=1-0.5*ep_1970s_2000s,
            ep_1970s_2010s=sum(diff_1970s_2010s), s_do_1970s_2010s=1-0.5*ep_1970s_2010s,
            ep_1980s_1990s=sum(diff_1980s_1990s), s_do_1980s_1990s=1-0.5*ep_1980s_1990s,
            ep_1980s_2000s=sum(diff_1980s_2000s), s_do_1980s_2000s=1-0.5*ep_1980s_2000s,
            ep_1980s_2010s=sum(diff_1980s_2010s), s_do_1980s_2010s=1-0.5*ep_1980s_2010s,
            ep_1990s_2000s=sum(diff_1990s_2000s), s_do_1990s_2000s=1-0.5*ep_1990s_2000s,
            ep_1990s_2010s=sum(diff_1990s_2010s), s_do_1990s_2010s=1-0.5*ep_1990s_2010s,
            ep_2000s_2010s=sum(diff_2000s_2010s), s_do_2000s_2010s=1-0.5*ep_2000s_2010s)
long_overlap<-pivot_longer(overlap_df,cols=matches("s_do"),
                           names_to="decades",values_to="s_do")%>%
  dplyr::select(-matches("^ep"))%>%
  mutate(decades=gsub("s_do_","",decades),
         species1=species,sizecat1=sizecat)%>%
  separate(decades,into=c("decade1","decade2"))


overlap_wide<-pivot_wider(long_overlap,
                          id_cols=c(species,sizecat,decade1),
                          names_from=c(species1,sizecat1,decade2),
                          values_from = s_do)
overlap_mat<-as.matrix(overlap_wide[,4:ncol(overlap_wide)])
rownames(overlap_mat)<-paste(overlap_wide$species,overlap_wide$sizecat,overlap_wide$decade1,sep="_")

#Visualize
ggplot(long_overlap)+
  geom_tile(aes(paste(species1,sizecat1,decade1),paste(species,sizecat,decade1),fill=s_do))+
  geom_text(aes(paste(species1,sizecat1,decade1),paste(species,sizecat,decade1),
                label=round(s_do,digits=2)),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

#or individually
species_mats<-list()
for (i in 1:length(species_allDecades)) {
  #overlap matrix
  df<-filter(long_overlap,paste(species,sizecat,sep="_")==species_allDecades[i])
  species_mats[[i]]<-ggplot(df)+
    geom_tile(aes(paste(species1,sizecat1,decade1),paste(species,sizecat,decade2),fill=s_do))+
    geom_text(aes(paste(species1,sizecat1,decade1),paste(species,sizecat,decade2),label=round(s_do,digits=2)),size=4)+
    scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white",
                         limits = c(0, 1), oob = scales::squish)+
    scale_color_manual(values=c(viridis::viridis(100,option="B")[100]),guide=NULL)+
    scale_x_discrete(name="",expand=expansion(0),labels=unique(df$decade1))+
    scale_y_discrete(name="",expand=expansion(0),labels=unique(df$decade2))+
    ggtitle(paste(df$species,df$sizecat))+
    theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
          plot.title.position="plot",plot.title=element_text(hjust=0.5))
  print(species_mats[[i]])
  ggsave(species_mats[[i]],file=paste0("../Figures/2022 GL Re-up/correct/species/",species_allDecades[i],"_overlapMat.png"))
  #clustering
  dist<-t(1-overlap_mat[which(grepl(species_allDecades[i],rownames(overlap_mat))),
                        which(grepl(species_allDecades[i],rownames(overlap_mat)))])
  dist<-rbind(c(NA,NA,NA,NA),dist)
  dist<-cbind(dist,c(NA,NA,NA,NA,NA))
  rownames(dist)[1]<-colnames(dist)[1]
  rownames(dist)<-str_extract(rownames(dist),"[0-9]+")
  dist<-as.dist(dist,diag=T)
  #dendrogram of clusters
  dend<-as.dendrogram(hclust(dist,method="average"))
  png(filename=paste0("../Figures/2022 GL Re-up/correct/species/",species_allDecades[i],"_dend.png"),
      width=650,height=650)
  dend %>% 
    set("branches_lwd", 4) %>%
    # Custom labels
    set("labels_cex", 1) %>%
    #set("labels_col", value = viridis::viridis(14,end=0.8),h = sigGuild_gl) %>%
    #set("branches_k_color", value = viridis::viridis(14,end=0.8), h = sigGuild_gl) %>%
    plot(horiz=TRUE,main=species_allDecades[i],axes=T,xlab="")
  dev.off()
}



#How the patterns look across the different species
decadePatterns<-read_csv("species_acrossDecades_Patterns.v2.csv")

ggplot(decadePatterns)+
  geom_bar(aes(Pattern_desc,fill=Pattern),show.legend = F)
table(decadePatterns$Pattern)
26/58*100
8/58*100
19/58*100
5/58*100




# Trying the first decade to the second -----------------------------------




#Cropping down the Original time, and clustering again
shared_props_new<-preyTrawls%>%
  filter(year<=1983
         & !is.na(GL_pyscinam)
         & GL_pyscinam!="Empty")%>%
  mutate(species=str_to_title(pdcomnam))%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  group_by(species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()

#Cropping down the new time and clustering again
shared_props_d2<-preyTrawls%>%
  filter(year>1983 & year<=1993
         & !is.na(GL_pyscinam)
         & GL_pyscinam!="Empty")%>%
  mutate(species=str_to_title(pdcomnam))%>%
  group_by(species,sizecat2)%>%
  mutate(total_w=sum(pyamtw),
         nDiets=n_distinct(dietID),
         total_v=sum(pyamtv))%>%
  group_by(species,sizecat2,GL_pyscinam)%>% #Just doing it for the general cats for now, trying with Garrison and Link categories
  summarise(prop_w=sum(pyamtw)/total_w,
            prop_v=sum(pyamtv)/total_v)%>%unique()%>%ungroup()

newd2Shared<-inner_join(dplyr::select(shared_props_new,species,sizecat2)%>%distinct(),
                        dplyr::select(shared_props_d2,species,sizecat2)%>%distinct())

shared_props_new<-filter(shared_props_new,paste(species,sizecat2) %in% paste(newd2Shared$species,newd2Shared$sizecat2))
shared_props_d2<-filter(shared_props_d2,paste(species,sizecat2) %in% paste(newd2Shared$species,newd2Shared$sizecat2))


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


overlap_clust_shared_new<-hclust(as.dist(1-overlap_mat_shared_new),method="average")
#Trying dendextend to pretty up the dendrogram
dend_shared_new<-as.dendrogram(overlap_clust_shared_new)



#Following through to get the "real" overlap matrix
props1_d2<-shared_props_d2[,-max(ncol(shared_props_d2))] #Cut off the volume prop for cleanliness
colnames(props1_d2)<-c("species1","sizecat1","GL_pyscinam","prop_w1")
props2_d2<-shared_props_d2[,-max(ncol(shared_props_d2))]
colnames(props2_d2)<-c("species2","sizecat2","GL_pyscinam","prop_w2")

overlap_shared_d2<-full_join(props1_d2,props2_d2)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat_shared_d2<-pivot_wider(overlap_shared_d2,
                                   id_cols=c(species1,sizecat1),
                                   names_from=c(species2,sizecat2),
                                   values_from = s_do)
overlap_mat_shared_d2<-as.matrix(overlap_mat_shared_d2[,3:ncol(overlap_mat_shared_d2)])
#rownames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#colnames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#overlap_mat[overlap_mat==1]<-NA
for (i in 1:nrow(overlap_mat_shared_d2)) {
  for (j in 1:ncol(overlap_mat_shared_d2)) {
    overlap_mat_shared_d2[i,j]<-ifelse(j>i,NA,overlap_mat_shared_d2[i,j])
  }
}


overlap_clust_shared_d2<-hclust(as.dist(1-overlap_mat_shared_d2),method="average")
#Trying dendextend to pretty up the dendrogram
dend_shared_d2<-as.dendrogram(overlap_clust_shared_d2)





####Ideas####
#Create a tanglegram of the dendrograms
#Different categories exist in the two (because of L Little Skate)
#How many times is the predator_sizecat in a new cluster, for each predator
#Direction of moves between guilds? i.e., are piscivores becoming benthivores?
#Mean guild number, if the numbers can be considered ordinal
#Guild membership size, what's growing and what's decreasing





# Linking Variables to Species Foraging Success ---------------------------







# OLD ---------------------------------------------------------------------




#SHOULD THESE BE MEANS ACROSS SIZECAT, THE ABUNDANCES AREN'T THEY'RE JUST AT SPECIES




#Following through to get the "real" overlap matrix
props1_gl<-trawlDietSum[,c("pdcomnam","sizecat2","GL_pyscinam","Wk")] #Cut off the volume prop for cleanliness
colnames(props1_gl)<-c("species1","sizecat1","gl_prey","prop_w1")
props2_gl<-trawlDietSum[,c("pdcomnam","sizecat2","GL_pyscinam","Wk")]
colnames(props2_gl)<-c("species2","sizecat2","gl_prey","prop_w2")

overlap_schoener_gl<-full_join(props1_gl,props2_gl)%>%
  mutate(diff=abs(prop_w1-prop_w2))%>%
  group_by(species1,sizecat1,species2,sizecat2)%>%
  summarise(ep=sum(diff),
            s_do=1-0.5*ep)

overlap_mat_gl<-pivot_wider(overlap_schoener_gl,id_cols=c(species1,sizecat1),names_from=c(species2,sizecat2),values_from = s_do)
overlap_mat_gl<-as.matrix(overlap_mat_gl[,3:ncol(overlap_mat_gl)])
#rownames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
#colnames(overlap_mat_gl)<-gsub(" ","\n",colnames(overlap_mat_gl))
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
  select(-c("Empty","sizecat"))
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
  mutate(species1=colnames(overlap_mat_gl))%>%
  pivot_longer(cols=colnames(overlap_mat_gl),names_to = "species2", values_to = "s_do")%>%
  mutate(sig=ifelse(s_do>sigGuild_gl,"bold.italic","plain"))%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=s_do))+
  geom_text(aes(species1,species2,label=round(s_do,digits=2),fontface=sig,color=sig),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100],viridis::viridis(100,option="B")[sigGuild_gl*100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  ggtitle("All Years, cluster means")+
  theme(panel.border = element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

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
png(file = "../Figures/2021 Garrison and Link Rep/clusterMeans_allYearsDend.png",   # The directory you want to save the file in
    width = 1250, # The width of the plot in inches
    height = 700) # The height of the plot in inches
par(mar=c(5,1,2,12))

dend %>% 
  set("branches_lwd", 4) %>%
  # Custom labels
  set("labels_cex", 1) %>%
  #set("labels_col", value = viridis::viridis(14,end=0.8),h = sigGuild_gl) %>%
  #set("branches_k_color", value = viridis::viridis(14,end=0.8), h = sigGuild_gl) %>%
  plot(horiz=TRUE,main="                               1973-1997 Dendrogram",axes=F,xlab="Similarity")
axis(side=1,at=c(0.6,0.5,0.4,0.3,0.2,0.1,0),labels=c(0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(v=sigGuild_gl,lty=2)
#rect.dendrogram(dend, h=0.5, which=c(1:6),border="transparent",
#                lty = 5, lwd = 0, col=rgb(0.2,0.2,0.2,0.15),horiz=T) #This is a way to add boxes indicating the guilds

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






# connecting --------------------------------------------------------------



#Blending the food habits and the trawls (cleaned from GMRI, Adam Kettering)


rm(list=ls())

# Top ---------------------------------------------------------------------



setwd("C:/Users/nh1087/OneDrive - USNH/Documents/NECC/Trawl Data")


library(readr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(lubridate)
library(ggmap)


load("googleMap.zoom6.eastCoast.R")


allGMRItrawls_clean<-read_csv("NMFS Trawls/Complete/NMFS_survdat_gmri_tidy.csv",col_types=list(id=col_character()))
GMRItrawls_totals<-dplyr::select(allGMRItrawls_clean,-c(numlen,numlen_adj,length_cm,catchsex,n_len_class))%>%
  mutate(est_day=day(est_towdate))%>%
  group_by(svspp,cruise6,station,stratum)%>%
  mutate(abundance=sum(abundance),biomass_kg=sum(biomass_kg))%>%
  distinct()%>%
  mutate(stratum_full=str_pad(stratum, width = 5, pad = "0", side = "left"),
         tow_full=str_pad(tow, width = 3, pad = "0", side = "left"),
         station_full=str_pad(station, width = 4, pad = "0", side = "left"),
         ID=paste0(cruise6,stratum_full,tow_full,station_full))

load("../Diet Data/prey19.RData")
prey19<-prey19%>%
  mutate(svspp=ifelse(svspp<10,paste0("00",svspp),ifelse(svspp<100,paste0("0",svspp),as.character(svspp))),
         station=ifelse(station<10,paste0("00",station),ifelse(station<100,paste0("0",station),as.character(station))),
         season=str_to_title(season),
         geoarea=as.character(geoarea),
         pdcomnam=tolower(pdcomnam),
         pdscinam=str_to_sentence(pdscinam),
         declon=-1*declon)%>%
  rename(est_year=year,est_month=month,est_day=day,comname=pdcomnam)


preyTrawls<-left_join(prey19,GMRItrawls_totals)%>%
  distinct()

strata_key <- list(
  "Georges Bank"          = as.character(13:23),
  "Gulf of Maine"         = as.character(24:40),
  "Southern New England"  = stringr::str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
  "Mid-Atlantic Bight"    = as.character(61:76))


preyTrawls_filtered<-filter(preyTrawls,
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
  filter(season%in%c("Spring","Fall") & survey_area!="stratum not in key")
#proportion of diets without trawl info after filtering out accounted for trawls (filtered from GMRI or Summer-Winter)
nrow(filter(preyTrawls_filtered,is.na(tow)))/nrow(preyTrawls_filtered) 
#How many "good" diets are lost? None
nrow(filter(preyTrawls_filtered,!is.na(tow)))/nrow(filter(preyTrawls,!is.na(tow)))

#Just cut out the diets that don't have trawl data (whether because GMRI didn't provide it or it's missing entirely)
preyTrawls_final<-filter(preyTrawls_filtered,!is.na(tow))
#Proportion of diets that remain
nrow(preyTrawls_final)/nrow(preyTrawls)

write.csv(preyTrawls_final,"../Diet Data/SUMMARIZED/merged_prey_GMRItrawls.csv",row.names = F)


# Weighting Diets by Trawl Capture ----------------------------------------

clusterSD <- function(n,tM,M,pq,FW) { #Function to calculate the SD for cluster means, from Buckel et al. 1999
  var<- (1/(n*(tM/n)^2))*((sum(M^2*(pq-FW)^2))/n-1)
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
trawlDiets<-preyTrawls_final%>%
  mutate(pdid=str_pad(pdid,width=6,pad="0",side="left"),
         dietID=paste(svspp,cruise6,station,pdsex,pdid,str_pad(pdlen,width=3,pad="0",side="left"),sep="0"))%>%
  group_by(ID,pdscinam)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv),
         totalN=n_distinct(dietID))%>%
  group_by(ID,pdscinam,collsci)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(totalN==0,0,n_distinct(dietID)/totalN))%>%
  dplyr::select(geoarea,est_year,season, #Can add in other variables here if useful, like lat-long or temps, they're for grouping later
                ID,abundance,biomass_kg,pdscinam,comname,collsci,totalwt,totalv,totalN,qikw,qikv,pik)%>% 
  distinct() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.


nDiets<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(geoarea,est_year,season,pdscinam,ID,totalN)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(geoarea,est_year,season,pdscinam)%>%summarise(nDiets=sum(totalN))#Calculate the number of diets
sumAbun<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(geoarea,est_year,season,pdscinam,ID,abundance)%>%distinct()%>%
  group_by(geoarea,est_year,season,pdscinam)%>%summarise(sumAbun=sum(abundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(geoarea,est_year,season,pdscinam,ID,abundance)%>%distinct()%>%
  group_by(geoarea,est_year,season,pdscinam)%>%summarise(sumAbun_nEmpty=sum(abundance))


trawlDietSum<-trawlDiets%>%
  group_by(geoarea,est_year,season)%>%
  mutate(nTrawls=n_distinct(ID))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(geoarea,est_year,season,pdscinam,comname,nTrawls,nDiets,sumAbun,sumAbun_nEmpty,collsci)%>%
  summarise(Wk=sum(abundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(abundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(abundance*pik)/sumAbun, #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,abundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,abundance,qikw,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,abundance,pik,Fk))%>%distinct()
trawlDietSum[is.na(trawlDietSum)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum%>%group_by(geoarea,est_year,season,pdscinam)%>%summarise(N=sum(Wk))
sum(check$N<0&check$N>0) #When this is 0, then you have all your means correct because they sum to 1 (or they're all empties)

#Saving these, so they can be better used within other analyses
#write.csv(trawlDiets,"trawl_speciesClusterCompositions.csv",row.names=F)
#write.csv(trawlDietSum,"geoareayearseason_speciesClusterCompositions.csv",row.names=F)


# Investigating the Missing Trawls ----------------------------------------


stillMissing<-filter(preyTrawls_filtered,is.na(tow))%>%select(cruise6,station,stratum)%>%distinct()


#just curious when the survey area doesnt match the given geoarea
curious<-filter(strata_sf_regions,geoarea!=survey_area)


oldNMFS<-bind_rows(springT,fallT)
colnames(oldNMFS)<-tolower(colnames(oldNMFS))
oldNMFS<-mutate(oldNMFS,cruise6=as.double(cruise6),station=str_pad(as.character(station), width = 3, pad = "0", side = "left"))%>%
  dplyr::select(1:2,4:5)%>%
  distinct()
fullCheck<-left_join(strata_sf_regions,oldNMFS)%>%
  arrange(id)%>%
  rename("id_fromNMFSrawData"="id")

#write.csv(fullCheck,"possiblyMissingTrawls.v2.csv",row.names=F)


#Checking in the other NMFS data I have
check<-filter(targetspecies_catchLocations,CRUISE6=="197411"&STATION=="270")
check1<-filter(allGMRItrawls_clean,station=="294"&cruise6=="197303")
check2<-filter(prey19,station=="583"&cruise6=="197910")



numFoodTrawls<-prey19%>%
  ungroup()%>%
  dplyr::select(cruise6,station,stratum,est_month,est_day,est_year)%>%
  distinct()%>%nrow()
numTrawls<-allGMRItrawls_clean%>%
  ungroup()%>%
  dplyr::select(cruise6,station,stratum,est_month,est_day,est_year)%>%
  distinct()%>%nrow()



# NMFS Redownload ---------------------------------------------------------


spring<-read_csv("NMFS Trawls/2022 Redownload/Spring/22561_UNION_FSCS_SVCAT.csv",col_types="ccccccccnn")
fall<-read_csv("NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv",col_types="ccccccccnn")
summer<-read_csv("NMFS Trawls/2022 Redownload/Summer/22562_UNION_FSCS_SVCAT.csv",col_types="ccccccccnn")
winter<-read_csv("NMFS Trawls/2022 Redownload/Winter/22563_UNION_FSCS_SVCAT.csv",col_types="ccccccccnn")

allTrawls<-bind_rows(spring,fall,summer,winter)
spfa<-bind_rows(spring,fall)

allGMRItrawls_clean<-read_csv("NMFS Trawls/Complete/NMFS_survdat_gmri_tidy.csv",
                              col_types="cccccccnnncnnTnnnnnnnnnnccncn")%>%
  mutate(stratum_full=str_pad(stratum, width = 5, pad = "0", side = "left"),
         tow_full=str_pad(tow, width = 3, pad = "0", side = "left"),
         station_full=str_pad(station, width = 4, pad = "0", side = "left"),
         ID=paste0(cruise6,stratum_full,tow_full,station_full))

compare<-right_join(allGMRItrawls_clean,spfa)

compareTest<-filter(compare,is.na(id))%>%
  mutate(STRATUM=as.numeric(STRATUM))%>%
  filter(STRATUM >= 01010,
         STRATUM <= 01760,
         STRATUM != 01310,
         STRATUM != 01320,
         STRATUM != 01330,
         STRATUM != 01350,
         STRATUM != 01410,
         STRATUM != 01420,
         STRATUM != 01490)%>%select(33:37)%>%distinct()
table(compareTest$CRUISE6) #Missing mostly from early years


#From above, the still missing trawls.
stillMissing
#are any of them in the redownload?
spfa$simpleID<-paste0(spfa$CRUISE6,spfa$STATION,spfa$STRATUM)
stillMissing$simpleID<-paste(stillMissing$cruise6,stillMissing$station,stillMissing$stratum,sep="0")
found<-stillMissing[which(stillMissing$simpleID%in%spfa$simpleID),]






springB<-read_csv("NMFS Trawls/2022 Redownload/Spring/22561_UNION_FSCS_SVLEN.csv",col_types="ccccccccnnccnn")
fallB<-read_csv("NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVLEN.csv",col_types="ccccccccnnccnn")
summerB<-read_csv("NMFS Trawls/2022 Redownload/Summer/22562_UNION_FSCS_SVLEN.csv",col_types="ccccccccnnccnn")
winterB<-read_csv("NMFS Trawls/2022 Redownload/Winter/22563_UNION_FSCS_SVLEN.csv",col_types="ccccccccnnccnn")

allTrawlsB<-bind_rows(springB,fallB,summerB,winterB)




# trawl exploration -------------------------------------------------------



rm(list=ls())

# Top ---------------------------------------------------------------------

setwd("C:/Users/nh1087/OneDrive - USNH/Documents/NECC/Trawl Data")


library(tidyverse)
library(ggplot2)
library(readr)
library(stringr)
library(lubridate)
library(ggmap)
library(rphylopic)
library(devtools)
library(ggpattern)



#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}

seasonPal2<-c("deepskyblue4","yellowgreen","lightgoldenrod1","orange3")

#Data loading
load("googleMap.zoom6.eastCoast.R")

load("NF.prey19.RData")
dietSeasons<-prey19%>%
  mutate(doy=yday(ymd(paste(year,month,day,sep="-"))),
         season=str_to_title(season),
         season=factor(season,levels=c("Winter","Spring","Summer","Fall")))%>%
  group_by(season)%>%
  summarise(seasonStart=min(doy,na.rm=T),seasonEnd=max(doy,na.rm=T))


springNMFS<-read_csv("NMFS Trawls/2022 Redownload/22561_FSCSTables/22561_UNION_FSCS_SVSTA.csv")
springNMFS<-springNMFS%>%
  mutate(day=as.numeric(yday(mdy_hms(BEGIN_EST_TOWDATE))),
         year=as.numeric(year(mdy_hms(BEGIN_EST_TOWDATE))),
         season="SPRING",
         data="Bottom Trawl")%>%
  select(year,day,season,lat=DECDEG_BEGLAT,lon=DECDEG_BEGLON,data)%>%
  filter(day<170) #Drop the one outlier trawl date from late June

fallNMFS<-read_csv("NMFS Trawls/2022 Redownload/22560_FSCSTables/22560_UNION_FSCS_SVSTA.csv")
fallNMFS<-fallNMFS%>%
  mutate(day=as.numeric(yday(mdy_hms(BEGIN_EST_TOWDATE))),
         year=as.numeric(year(mdy_hms(BEGIN_EST_TOWDATE))),
         season="FALL",
         data="Bottom Trawl")%>%
  select(year,day,season,lat=DECDEG_BEGLAT,lon=DECDEG_BEGLON,data)


nhmeT<-read_csv("NH-ME Trawls/MaineDMR_Trawl_Survey_Tow_Data_2021-01-20.csv")
nhmeT<-nhmeT%>%
  mutate(day=yday(Start_Date),
         data="Inshore Trawl")%>%
  dplyr::select(Year,day,lat=Start_Latitude,lon=Start_Longitude,data)%>%
  filter(day>100)%>%
  mutate(season=ifelse(day<200,"SPRING","FALL"))
nhmeFT<-filter(nhmeT,season=="FALL")
nhmeST<-filter(nhmeT,season=="SPRING")

ichthy<-read_csv("Ichthyoplankton/Plankton.csv")
ichthy<-ichthy%>%
  mutate(day=yday(date),
         data="Ichthyoplankton")%>%
  select(year,day,lat,lon,data)


seine<-read_csv("NHFG Seines/NHFG Seine Survey Data_FUREY_1-19-21_samples.csv")
seine<-seine%>%
  mutate(day=yday(ymd(paste(Year,Month,Day,sep="-"))),
         data="Beach Seine")%>%
  select(year=Year,day,data)



ternSil<-image_data("f164783d-3bab-45f3-9885-c1b382202369",size=1024)[[1]]
herringSil<-image_data("0b89df58-7eae-40c5-9676-d35e3449afb2",size=1024)[[1]]
hakeSil<-image_data("8d92b454-3131-4bbd-ac9c-e1df14c2fc5a",size=1024)[[1]]
amphSil<-image_data("fd6af059-2365-4a6e-806a-ce22bb537103",size=1024)[[1]]
zooSil<-image_data("3517b09b-f068-4d04-8f19-72c3dc66d85c",size=1024)[[1]]
silversideSil<-image_data("4e31c12d-1c23-4a64-bcc6-0730ab6d0617",size=1024)[[1]]
killiSil<-image_data("ffc40581-a3b2-4296-a901-9ff9d8cd4805",size=1024)[[1]]

ggplot()+
  geom_segment(data=ichthy,aes(x=min(day,na.rm=T),xend=max(day,na.rm=T),y=1,yend=1),
               lwd=15,color="palegreen2")+
  geom_segment(data=seine,aes(x=min(day,na.rm=T),xend=max(day,na.rm=T),y=2,yend=2),
               lwd=15,color="gold1")+
  geom_segment(data=nhmeST,aes(x=min(day,na.rm=T),xend=max(day,na.rm=T),y=3,yend=3),
               lwd=15,color="darkorange2")+
  geom_segment(data=nhmeFT,aes(x=min(day,na.rm=T),xend=max(day,na.rm=T),y=3,yend=3),
               lwd=15,color="turquoise2")+
  geom_segment(data=fallNMFS,aes(x=min(day,na.rm=T),xend=max(day,na.rm=T),y=4,yend=4),
               lwd=15,color="dodgerblue2")+
  geom_segment(data=springNMFS,aes(x=min(day,na.rm=T),xend=max(day,na.rm=T),y=4,yend=4),
               lwd=15,color="firebrick2")+
  geom_segment(data=dietSeasons,aes(x=seasonStart,xend=seasonEnd,y=5,yend=5,color=season,alpha=season),
               lwd=15,show.legend = F)+
  #Adding silhouettes of the different groups
  add_phylopic(herringSil,alpha=1,x=35,y=4.17,ysize=9)+
  add_phylopic(hakeSil,alpha=1,x=30,y=3.84,ysize=10)+
  add_phylopic(herringSil,alpha=1,x=100,y=2.86,ysize=8)+
  add_phylopic(hakeSil,alpha=1,x=97,y=3.15,ysize=8)+
  add_phylopic(silversideSil,alpha=1,x=140,y=2.11,ysize=5)+
  add_phylopic(killiSil,alpha=1,x=141,y=1.9,ysize=4)+
  add_phylopic(amphSil,alpha=1,x=-3,y=1.1,ysize=7)+
  add_phylopic(zooSil,alpha=1,x=-4,y=0.94,ysize=10)+
  scale_x_continuous(name="Day of Year",limits=c(-20,380),expand=expansion(add=0),
                     breaks=c(1,92,183,275,366),
                     labels=c("Jan 1","Apr 1","Jul 1","Oct 1","Dec 31"))+
  scale_y_continuous(name="",expand=expansion(mult=0.1),
                     breaks=c(1,2,3,4,5),
                     labels=c("NMFS Ichthyoplankton","NHFG Beach Seine",
                              "NH-ME Inshore Trawl","NMFS Bottom Trawl","NMFS Diet Collections"))+
  scale_color_manual(values=seasonPal2)+
  scale_alpha_manual(values=c(0.8,0.8,1,1))+
  theme(panel.grid.minor.y=element_blank())



# Unique trawl locations and areas ----------------------------------------

allTrawls<-bind_rows(fallNMFS,springNMFS)%>%
  dplyr::select(STRATUM,STATION,AREA,DECDEG_BEGLAT,DECDEG_BEGLON)%>%
  distinct()

#write.csv(allTrawls,"allTrawls.csv",row.names = F)

# Actual Exploration of the NMFS Trawls -----------------------------------

#This is all copied from the exploration script, then adjusted to keep species too
springT<-read_csv("NMFS Trawls/Spring/22561_UNION_FSCS_SVCAT.csv", 
                  col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                   STATION = col_double(), STRATUM = col_double(), 
                                   TOW = col_double()))

springT$Year<-as.numeric(substr(springT$CRUISE6,1,4))

fallT<-read_csv("NMFS Trawls/Fall/22560_UNION_FSCS_SVCAT.csv", 
                col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                 STATION = col_double(), STRATUM = col_double(), 
                                 TOW = col_double()))

fallT$Year<-as.numeric(substr(fallT$CRUISE6,1,4))

#Extracting the species to keep, should ROCKLING be added? They were never aged...
speciesList<-fallT%>%
  filter(grepl(paste(c("^(?=.*silver)(?=hake).+","^(?=.*little)(?=.*skate).+","^(?=.*yellowtail)(?=.*flounder).+",
                       "haddock","^(?=.*white)(?=.*hake).+","^(?=.*spiny)(?=.*dogfish).+",
                       "^(?=.*summer)(?=.*flounder).+","^(?=.*red)(?=.*hake).+","pollock",
                       "Squalus acanthias","Merluccius bilinearis","Urophycis chuss","Paralichthys dentatus",
                       "Leucoraja erinacea","Limanda ferruginea","Melanogrammus aeglefinus ",
                       "Pollachius virens","Urophycis tenuis"),
                     collapse="|"),perl=T,ignore.case=T,LOGGED_SPECIES_NAME))%>%
  dplyr::select(LOGGED_SPECIES_NAME)%>%unique()

targetSpecies_trawls<-fallT%>%
  mutate(season="Fall")%>%
  bind_rows(mutate(springT,season="Spring",CATCHSEX=as.character(CATCHSEX)))%>%
  filter(LOGGED_SPECIES_NAME %in% speciesList$LOGGED_SPECIES_NAME)


trawlLocations_Fall<-read_csv("NMFS Trawls/Fall/22560_UNION_FSCS_SVSTA.csv", 
                              col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                               STATION = col_double(), STRATUM = col_double(), 
                                               TOW = col_double(), AREA = col_character(),
                                               EST_MONTH=col_double(),EST_DAY=col_double()))%>%
  dplyr::select(CRUISE6,ID,STATION,STRATUM,TOW,LONG=DECDEG_BEGLON,LAT=DECDEG_BEGLAT,
                Year=EST_YEAR,Month=EST_MONTH,Day=EST_DAY,Time=EST_TIME,AREA,DESSPEED,TOWDUR,SURFTEMP,BOTTEMP)

trawlLocations_Spring<-read_csv("NMFS Trawls/Spring/22561_UNION_FSCS_SVSTA.csv", 
                                col_types = cols(CRUISE6 = col_character(), ID = col_character(), 
                                                 STATION = col_double(), STRATUM = col_double(), 
                                                 TOW = col_double(), AREA=col_character(),
                                                 EST_MONTH=col_double(),EST_DAY=col_double()))%>%
  dplyr::select(CRUISE6,ID,STATION,STRATUM,TOW,LONG=DECDEG_BEGLON,LAT=DECDEG_BEGLAT,
                Year=EST_YEAR,Month=EST_MONTH,Day=EST_DAY,Time=EST_TIME,AREA,DESSPEED,TOWDUR,SURFTEMP,BOTTEMP)

targetspecies_catchLocations<-left_join(targetSpecies_trawls,bind_rows(trawlLocations_Fall,trawlLocations_Spring))%>%
  mutate(Date=ymd(paste(Year,Month,Day,sep="-")),
         doy=yday(Date))






#TRYING TO GET THE GEO AREA DEFINITIONS FROM THE PREY DF, BUT THEY AREN"T PERFECT
prey19%>%
  group_by(geoarea)%>%
  summarise(minSTRAT=min(stratum),maxSTRAT=max(stratum),
            minSTATION=min(station),maxSTATION=max(station))

fake<-prey19
colnames(fake)<-str_to_upper(colnames(fake))
test<-targetspecies_catchLocations%>%
  mutate(CRUISE6=as.integer(CRUISE6),
         SVSPP=as.integer(SVSPP),
         DECLON=-1*LONG,
         DECLAT=LAT)%>%
  dplyr::select(CRUISE6,AREA,DECLON,DECLAT)%>%
  left_join(dplyr::select(fake,CRUISE6,GEOAREA,DECLON,DECLAT))%>%unique()%>%
  mutate(GEOAREA=ifelse(is.na(GEOAREA),"UnId",as.character(GEOAREA)))%>%
  arrange(desc(GEOAREA))
table(test$GEOAREA)

ggmap(zoom6)+
  geom_point(data=test,aes(-DECLON,DECLAT,fill=GEOAREA),shape=21,alpha=0.6,size=3)+
  geom_polygon(data=test%>%group_by(GEOAREA)%>%slice(chull(-DECLON,DECLAT))%>%filter(GEOAREA!="UnId"),
               aes(x=-DECLON,y=DECLAT,fill=GEOAREA),color="black",alpha=0.1,lwd=1.5,show.legend = F)+
  geom_polygon(data=test%>%mutate(as.numeric(AREA))%>%filter(!is.na(AREA))%>%
                 group_by(AREA)%>%slice(chull(-DECLON,DECLAT)),
               aes(x=-DECLON,y=DECLAT,group=AREA),fill="white",color="black",alpha=0.1,lwd=1,show.legend = F)+
  scale_fill_manual(values=colorspace::rainbow_hcl(6,c=100,end=290))+
  ylab("Latitude")+xlab("Longitude")





nhmeC<-read_csv("NH-ME Trawls/MaineDMR_Trawl_Survey_Catch_Data_2021-01-20.csv")
unique(nhmeC$Common_Name)
seabass<-filter(nhmeC,Common_Name=="Sea Bass Black")




# what are sections -------------------------------------------------------

allTrawls<-bind_rows(fallNMFS,springNMFS)

stations<-select(allTrawls,CRUISE6,STATION,lat=DECDEG_BEGLAT,lon=DECDEG_BEGLON)%>%
  distinct()

nStations_cruise<-allTrawls%>%
  group_by(substr(CRUISE6,1,4))%>%
  summarise(nStation=n_distinct(ID))

nCruise_station<-allTrawls%>%
  group_by(STATION)%>%
  summarise(nCruises=n())

station1<-filter(allTrawls,STATION=="0001")
station2<-filter(allTrawls,STATION=="0002")

survAreas<-st_read("../Mapping Data/BTS_Strata.shp")
survAreas<-st_transform(survAreas,4326)

#Using shapefiles to build a clean map
world<-map_data("world")
states<-map_data("state")

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_spatial_point(data=station1,
                     aes(DECDEG_BEGLON,DECDEG_BEGLAT),
                     crs=4326,size=5)+
  geom_spatial_point(data=station2,
                     aes(DECDEG_BEGLON,DECDEG_BEGLAT),
                     crs=4326,size=5,color="red")+
  geom_sf(data=survAreas,fill="transparent",col="black")+
  scale_color_discrete_sequential(palette="Hawaii",name="Geographic\nRegion",l2=80)+
  guides(color=guide_legend(override.aes = list(size=9)))+
  annotation_scale(location="br",width_hint=0.4,height=unit(0.125,"in"),
                   pad_x=unit(0.25,"in"),pad_y=unit(0.4,"in"),text_cex=1.2)+
  geom_rect(aes(xmin=-63.35,xmax=-62.5,ymin=35.8,ymax=36.2),fill="azure1")+
  annotation_north_arrow(location="br",pad_x=unit(0.5,"in"),pad_y=unit(0.75,"in"),
                         style = north_arrow_fancy_orienteering(text_size=15,line_width=1.25),
                         height=unit(1,"in"),width=unit(0.9,"in"))+
  coord_sf(xlim=c(-77.8,st_bbox(survAreas)$xmax),
           ylim=c(35.2,st_bbox(survAreas)$ymax))+
  theme(legend.position=c(0.225,0.8),legend.background=element_rect(color="black"),
        panel.grid.major = element_line(color=rgb(0.5,0.5,0.5,alpha=0.7),linetype="dashed",size=0.5),
        panel.background = element_rect(fill="azure1"))+
  xlab("Longitude")+ylab("Latitude")

length(unique(allTrawls$CRUISE6))



# food habits exploration -------------------------------------------------


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
load("googleMap.zoom6.eastCoast.R")
vulScores<-read.csv("../Haleetal_climateScores.csv")


#ALL LEVELS OF YEAR AND SEASON
ALLyearseasons<-factor(levels=c(paste(sort(rep(seq(1973,2019),4)),rep(c("Winter","Spring","Summer","Fall"),length(seq(1973,2019))))))
levels(ALLyearseasons)

#The names of my species, just to have them nice and handy
myspecies_sci<-unique(str_to_sentence(prey19$pdscinam))
myspecies_com<-unique(str_to_title(prey19$pdcomnam))
myspecies_svspp<-unique(str_pad(as.character(prey19$svspp),width=3,pad="0",side="left"))

#Size classes used by Garrison and Link
sizeClasses<-read.csv("GarrisonLink_predSizeClasses.csv") 
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
n_distinct(NF.prey19$dietID)-n_distinct(preyGMRI_filter$dietID)#How many more diets
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



# Creating my own prey categories -----------------------------------------

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
  summarise(vol=sum(pyamtv),p=vol/totalV*100)%>%unique()

spGL_Groups<-prey19%>%
  group_by(pdcomnam)%>%
  mutate(totalV=sum(pyamtv))%>%
  group_by(pdcomnam,gl_prey)%>%
  summarise(vol=sum(pyamtv),p=vol/totalV*100)%>%unique()




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
uniqueDiets%>%
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
uniqueDiets%>%
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
  summarise(prop=n()/total)%>%distinct()%>%
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
  summarise(M=min(count,na.rm=T),
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
  summarise(M=min(count,na.rm=T),Myearseason=subset(yearseason,count==min(count,na.rm=T)))%>%
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
  summarise(N_empty=n_distinct(ID),
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
  summarise(N_empty=n_distinct(ID),
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
  summarise(N_empty=n_distinct(ID),
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
  summarise(N_empty=n_distinct(ID),
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
  summarise(N_empty=n_distinct(ID),
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




# Bayesian for class ------------------------------------------------------




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



# Bioenergetics -----------------------------------------------------------



setwd("C:/Users/nh1087/OneDrive - USNH/Documents/NECC/")


library(tidyverse)
library(ggplot2)
library(readr)
library(ggnewscale)



#Editing the table from Steimle et al. 1985, never have to do again
#preyEDs<-read.csv("Steimle1985_ED.csv")%>%
#  separate(col="Taxa",into=c("Taxa","N_combustions","Ash_perc","H2O_perc",
#                             "Dry_KJg_mean","Dry_KJg_SD","Ash_KJg","Wet_KJg","Shell_perc"),
#           sep=" ")
#
#write.csv(file="Steimle1985_ED_2.0.csv",preyEDs)


preyEDs<-read.csv("Diet Data/Steimle1988_ED_2.0.csv")





# Bioenergetics--Wuenschel Concept ----------------------------------------

#Starting with all the Deslauriers parameters
params<-read_csv("Bioenergetics/jim-breck-FB4-1e6be54/Parameters_official.csv")

#Copied from Deslauriers code
Cf2T <- function(df,Temperature) { ### Temperature function equation 2 (Hanson et al. 1997; equation from Kitchell et al. 1977)
  if (Temperature < df$CTM) { 
    V <- (df$CTM - Temperature) / (df$CTM - df$CTO)
    CY <- log(df$CQ) * (df$CTM - df$CTO + 2)
    CZ <- log(df$CQ) * (df$CTM - df$CTO)
    CX <- (CZ^2 * (1+(1+40/CY)^0.5)^2)/400
    ft <- V^CX * exp(CX * (1 - V))
  } else if (Temperature >= df$CTM) {ft  <-  0}  
  
  if(ft < 0) {ft  <-  0}  ## prevent negative values
  return(ft)
}


Cf3T <- function(df,Temperature) { ### Temperature function equation 3 (Hanson et al. 1997; equation from Thornton and Lessem 1978)
  CG1 = (1/(df$CTO - df$CQ))*log((0.98*(1-df$CK1))/(df$CK1*0.02))
  CG2 = (1/(df$CTL - df$CTM))*log((0.98*(1-df$CK4))/(df$CK4*0.02))
  L1 <- exp(CG1*(Temperature-df$CQ))
  KA <- (df$CK1*L1) / (1 + df$CK1*(L1-1))
  L2 <- exp(CG2*(df$CTL-Temperature))
  KB <- (df$CK4*L2) / (1 + df$CK4*(L2-1))
  ft <- KA * KB
  return(ft)
}

cTemps<-data.frame(Species=as.character(),Temp=as.numeric(),propCmax=as.numeric())
for (i in 1:nrow(params)) {
  if (is.na(params[i,"CEQ"])) {
    cTemps<-bind_rows(cTemps,data.frame(Species=params[i,"Species"],Temp=NA,propCmax=NA))
  } else if (params[i,"CEQ"]==2) {
    cTemps<-bind_rows(cTemps,data.frame(Species=params[i,"Species"],Temp=c(1:30),propCmax=Cf2T(params[i,],c(1:30))))
  } else if (params[i,"CEQ"]==3) {
    cTemps<-bind_rows(cTemps,data.frame(Species=params[i,"Species"],Temp=c(1:30),propCmax=Cf3T(params[i,],c(1:30))))
  } else {
    print("NA")
  }
}
ggplot(cTemps,aes(Temp,propCmax,color=Species))+
  geom_line()
#Just the cod and herring
ggplot(filter(cTemps,Species %in% c("Atlantic cod (juvenile & adult)","Baltic herring (YOY)")),aes(Temp,propCmax,color=Species))+
  geom_line(size=3)


consumption <- function(df,Temperature, W, p) { ### Consumption function
  Cmax <- df$CA * W ^ df$CB 
  if(is.na(df$CEQ)) {ft = NA   # reformatted to minimize if-tests; JEB
  } else if(df$CEQ == 2) {ft = Cf2T(df,Temperature)
  } else if(df$CEQ == 3) {ft = Cf3T(df,Temperature)
  } else {ft = NA}
  
  return(C=Cmax * p * ft)
}

plot(consumption(params[1,],c(1:30),100,1))
codherringMasses<-filter(fishMass,species_comnam %in% c("Atlantic cod","Atlantic herring"))


cTotal<-data.frame(Species=as.character(),Temp=as.numeric(),Weight=as.numeric(),C=as.numeric())
for (i in 1:nrow(params)) {
  cTotal<-bind_rows(cTotal,data.frame(Species=params[i,"Species"],
                                      Temp=sort(rep(0:30,6)),
                                      Weight=c(codherringMasses$Weight,
                                               codherringMasses$minWeight,
                                               codherringMasses$maxWeight),
                                      C=consumption(params[i,],sort(rep(0:30,6)),
                                                    c(codherringMasses$Weight,
                                                      codherringMasses$minWeight,
                                                      codherringMasses$maxWeight),1)))
}
ggplot(cTotal,aes(Temp,C,color=Species))+
  geom_line()
#Just the cod and herring
ggplot(filter(cTotal,(Species=="Atlantic cod (juvenile & adult)"&Weight==max(Weight)) |
                (Species=="Baltic herring (YOY)" & Weight==min(Weight))),
       aes(Temp,C,color=Species))+
  geom_line(size=3)


Rf1T <- function(df,Temperature) { ### Temperature function equation 1 (Hanson et al. 1997; Stewart et al. 1983)
  ft <- exp(df$RQ*Temperature)
  return(ft)
}

RACTf1T <- function(df,W,Temperature) { ### Temperature function equation 1 with activity component (Hanson et al. 1997; Stewart et al. 1983)
  if(Temperature <= df$RTL) {VEL <- df$ACT * W ^ df$RK4 * exp(df$BACT * Temperature)
  } else if(Temperature >  df$RTL) {VEL <- df$RK1 * W ^ df$RK4 * exp(df$RK5 * Temperature)}  
  ACTIVITY <- exp(df$RTO * VEL)
  return(ACTIVITY)
}

Rf2T <- function(df,Temperature) { ### Temperature function equation 2 (Hanson et al. 1997; Kitchell et al. 1977)
  if (Temperature< df$RTM) {
    RY <- log(df$RQ) * (df$RTM - df$RTO + 2)
    RZ <- log(df$RQ) * (df$RTM - df$RTO)
    RX <- (RZ^2 * (1+(1+40/RY)^0.5)^2)/400
    V <- (df$RTM - Temperature) / (df$RTM - df$RTO)
    ft <- V^RX * exp(RX * (1 - V))
  } else if (Temperature>=df$RTM) {ft <- 0.000001}
  
  if(ft < 0) {ft <- 0.000001}  
  return(ft)
}


respiration <- function(df,Temperature, W) { ### Respiration function
  Rmax <- df$RA * W ^ df$RB  
  if(df$REQ == 1) {
    ft <- Rf1T(df,Temperature)  
    ACTIVITY <- RACTf1T(df,W,Temperature) 
  } else if(df$REQ == 2) {
    ft <- Rf2T(df,Temperature)
    ACTIVITY <- df$ACT
  }
  R <- (Rmax * ft * ACTIVITY) 
  return(R)
}
plot(consumption(params[5,],c(-20:25),100,1),type="l")
lines(respiration(params[5,],c(-20:25),100))

rTotal<-data.frame(Species=as.character(),Temp=as.numeric(),Weight=as.numeric(),R=as.numeric())
for (i in 1:nrow(params)) {
  rTotal<-bind_rows(rTotal,data.frame(Species=params[i,"Species"],
                                      Temp=sort(rep(0:30,6)),
                                      Weight=c(codherringMasses$Weight,
                                               codherringMasses$minWeight,
                                               codherringMasses$maxWeight),
                                      R=respiration(params[i,],sort(rep(0:30,6)),
                                                    c(codherringMasses$Weight,
                                                      codherringMasses$minWeight,
                                                      codherringMasses$maxWeight))))
}
rTotal$R<-ifelse(rTotal$R>1,1,rTotal$R)
ggplot(rTotal,aes(Temp,R,color=Species))+
  geom_line()
#Just the cod and herring
ggplot()+
  geom_line(data=filter(rTotal,(Species=="Atlantic cod (juvenile & adult)"&Weight==1000) |
                          (Species=="Baltic herring (YOY)" & Weight==10)),
            aes(Temp,R,color=Species),size=3)+
  geom_line(data=filter(cTotal,(Species=="Atlantic cod (juvenile & adult)"&Weight==1000) |
                          (Species=="Baltic herring (YOY)" & Weight==10)),
            aes(Temp,C,color=Species))+
  facet_wrap(~Species,scales="free")

egestion1 <- function(df,C) {
  # egestion
  Eg = df$FA * C
  return(Eg)
}

cTotal$Eg<-egestion1(params[5,],cTotal$C)

SpDynAct <- function(df,C,Eg) { ### Specific dynamic action function (Hanson et al. 1997)
  S <- df$SDA *(C-Eg)
  return(S)
}

cTotal$SDA<-SpDynAct(params[5,],cTotal$C,cTotal$Eg)

crTotal<-left_join(cTotal%>%mutate(Weight=round(Weight,0)),
                   rTotal%>%mutate(Weight=round(Weight,0)))%>%distinct()

codherringCRtotal<-bind_rows(filter(crTotal, Species=="Atlantic cod (juvenile & adult)")[c(T,F),],
                             filter(crTotal, Species=="Baltic herring (YOY)")[!c(T,F),])%>%
  group_by(Species)%>%
  mutate(weightGroup=case_when(Weight==min(Weight)~"Min",
                               Weight==max(Weight)~"Max",
                               TRUE~"Mean"),
         Species=case_when(grepl("herring",Species)~"Atlantic herring",
                           grepl("cod",Species)~"Atlantic cod"))%>%
  pivot_wider(id_cols=c("Species","Temp"),
              names_from="weightGroup",values_from=c("C","Eg","SDA","R"))

ggplot(codherringCRtotal)+
  geom_ribbon(aes(Temp,ymin=R_Min+SDA_Min,ymax=R_Max+SDA_Max,
                  fill="R+SDA",color="R+SDA"),alpha=0.1,size=0.1)+
  geom_ribbon(aes(Temp,ymin=C_Max-Eg_Max,ymax=C_Min-Eg_Min,
                  fill="C",color="C"),alpha=0.1,size=0.1)+
  geom_line(aes(Temp,R_Mean+SDA_Mean,color="R+SDA"),size=2)+
  geom_line(aes(Temp,C_Mean-Eg_Mean,color="C"),size=2)+
  facet_wrap(~Species,scales="free")+
  scale_fill_manual(values=c("dodgerblue2","firebrick2"),name="")+
  scale_color_manual(values=c("dodgerblue2","firebrick2"),name="")+
  theme_bw()+theme(legend.position="top")+
  ylab("Energy and materials (g/g/day)")+xlab("Temperature")




# Actual Data on Temp and Consumption from Trawls -------------------------------------------------------------

#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

seasonPal2<-c("deepskyblue4","yellowgreen","lightgoldenronew","orange3") #BEST

#Data loading
load("Diet Data/NF.prey19.RData")
allprey<-NF.prey19
rm(NF.prey19)

#ALL LEVELS OF YEAR AND SEASON
ALLyearseasons<-factor(levels=c(paste(sort(rep(seq(1973,2019),4)),rep(c("Winter","Spring","Summer","Fall"),length(seq(1973,2019))))))

#The names of my species, just to have them nice and handy
species_sci<-unique(str_to_sentence(allprey$pdscinam))
species_com<-unique(str_to_title(allprey$pdcomnam))

#These are all the species in the trawls
spec<-read.csv("Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv")%>%
  dplyr::select(LOGGED_SPECIES_NAME,SVSPP)%>%distinct()%>%
  filter(SVSPP>="001")%>%
  mutate(LOGGED_SPECIES_NAME=str_to_sentence(LOGGED_SPECIES_NAME))
spec<-spec[!duplicated(spec$SVSPP),]


#Size classes used by Garrison and Link
sizeClasses<-read.csv("Diet Data/GarrisonLink_predSizeClasses.csv")%>%
  left_join(spec,by=c("species_comnam"="LOGGED_SPECIES_NAME"))%>%
  rename(svspp=SVSPP)
#These are replicated in the sizecats column already present, except that the small and largest cats extend to the smallest and largest individuals
allprey%>%group_by(pdcomnam,sizecat)%>%summarise(min=min(pdlen),max=max(pdlen))

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
  left_join(sizeClasses)%>%
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


#The categories from Garrison and Link, 2000
gl_preycats<-read.csv("Diet Data/GarrisonLink_preyCats.csv")%>%
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


check<-filter(allprey,is.na(GL_pycomnam))[,c("gensci","analsci","collsci","pynam","INnum","INpy")]
sort(table(allprey$GL_pyscinam))
printOut<-allprey%>%
  group_by(gensci,analsci,collsci,pynam,GL_pyscinam)%>%
  summarise(N=n())%>%
  arrange(GL_pyscinam)
#write.csv(printOut,"allprey_glpreycats.v15.csv",row.names = F)

uniqueDiets<-allprey%>%
  dplyr::select(cruise6,station,svspp,pdsex,pdid,pdcomnam,
                pdscinam,dietID,pdlen,pdwgt,sizecat,pdgutw,pdgutv,
                declat,declon,month,day,year,season,geoarea,id)%>%
  distinct()


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

allprey<-mutate(allprey,GL_pyscinam=factor(GL_pyscinam,levels=gl_levels))


#reading in prey trawls df from "connecting_food+trawls.R" where prey19 is merged with GMRI's clean data for trawls
fullTrawls<-read_csv("Trawl Data/NMFS Trawls/Complete/NMFS_survdat_gmri_tidy.csv",
                     col_types=c("cccccccnnncnnTnnnnnnnnnnccccn"))
trawls<-fullTrawls%>%
  left_join(sizeClasses)%>%
  filter(svspp %in% allprey$svspp)%>%
  group_by(svspp)%>%
  mutate(sizecat2=case_when(dplyr::between(length_cm,min(small_min,na.rm=T), min(small_max,na.rm=T))~"S",
                            dplyr::between(length_cm,min(medium_min,na.rm=T),min(medium_max,na.rm=T))~"M",
                            dplyr::between(length_cm,min(large_min,na.rm=T), min(large_max,na.rm=T))~"L",
                            dplyr::between(length_cm,min(xlarge_min,na.rm=T),min(xlarge_max,na.rm=T))~"XL",
                            TRUE~"S"))%>%
  group_by(svspp,id)%>%
  mutate(species_abundance=sum(numlen_adj))%>%
  group_by(svspp,id,sizecat2)%>%
  mutate(sizecat_abundance=sum(numlen_adj))%>%
  dplyr::select(-c(catchsex,abundance,biomass_kg,length_cm:numlen_adj,n_len_class,small_min:xlarge_max))%>%
  distinct()
#Want to use just one row for each trawl, so combined the sexes and removed the length classes
preyTrawls<-left_join(allprey,trawls)%>%
  group_by(svspp,id,sizecat2)%>%
  mutate(sizecat_abundance=ifelse(is.na(sizecat_abundance),n_distinct(dietID),sizecat_abundance))

uniqueDietTrawls<-preyTrawls[!duplicated(preyTrawls$dietID),]
#There is one diet that's in here twice as it was both examined at sea and preserved, will just drop it
uniqueDietTrawls<-filter(uniqueDietTrawls,dietID != "19720000010222009021591080")

bioenergDietSpecies<-uniqueDietTrawls%>%
  ungroup()%>%rowwise()%>%
  mutate(Temp=sum(surftemp,bottemp,na.rm=T)/2,
         C=ifelse(pdgutw==0,0,pdgutw/pdwgt))%>%
  group_by(id,species_comnam)%>%
  mutate(C_mean=mean(C,na.rm=T),
         T_mean=mean(Temp,na.rm=T),
         bT_mean=mean(bottemp,na.rm=T),
         C_mean=ifelse(is.nan(C_mean),NA,C_mean),
         T_mean=ifelse(is.nan(T_mean),NA,T_mean),
         bT_mean=ifelse(is.nan(bT_mean),NA,bT_mean))%>%
  group_by(year,species_comnam)%>%
  summarise(year_abundance=sum(species_abundance),
            C_simplemean=mean(C_mean),
            C_clustermean=sum(C_mean*species_abundance)/year_abundance,
            C_SD=sd(C_mean),
            T_clustermean=sum(T_mean*species_abundance)/year_abundance,
            T_SD=sd(T_mean),
            bT_clustermean=sum(bT_mean*species_abundance)/year_abundance)


#More detailed, kept for a point at every trawl
bioenergDietTrawls<-uniqueDietTrawls%>%
  ungroup()%>%rowwise()%>%
  mutate(Temp=sum(surftemp,bottemp,na.rm=T)/2,
         C=ifelse(pdgutw==0,0,pdgutw/pdwgt))%>%
  group_by(id,species_comnam)%>%
  mutate(C_mean=mean(C,na.rm=T),
         T_mean=mean(Temp,na.rm=T),
         bT_mean=mean(bottemp,na.rm=T),
         C_mean=ifelse(is.nan(C_mean),NA,C_mean),
         T_mean=ifelse(is.nan(T_mean),NA,T_mean),
         bT_mean=ifelse(is.nan(bT_mean),NA,bT_mean))


#Mean weight for each species in the FHD
fishMass<-uniqueDietTrawls%>%
  group_by(species_comnam,id)%>%
  mutate(m=mean(pdwgt,na.rm=T))%>%
  group_by(species_comnam)%>%
  summarise(Weight_simple=mean(pdwgt,na.rm=T),
            Weight=sum(m*species_abundance,na.rm=T)/sum(species_abundance,na.rm=T),
            minWeight=min(pdwgt,na.rm=T)+1,
            maxWeight=max(pdwgt,na.rm=T))



ggplot(filter(bioenergDietTrawls,species_comnam %in% c("Atlantic cod","Atlantic herring")))+
  geom_point(aes(T_mean,C_mean,color=species_comnam,size=year))




# Plotting together -------------------------------------------------------


codherring_bioenerg<-filter(crTotal,Species %in% c("Atlantic cod (juvenile & adult)","Baltic herring (YOY)"))%>%
  mutate(species_comnam=ifelse(Species=="Atlantic cod (juvenile & adult)","Atlantic cod","Atlantic herring"))%>%
  filter((species_comnam=="Atlantic cod" & Weight==1000) |
           (species_comnam=="Atlantic herring" & Weight==10))

codherring_dietsTrawl<-filter(bioenergDietTrawls,
                              species_comnam %in% c("Atlantic cod", "Atlantic herring") &
                                year>=1990)%>%
  mutate(Species=species_comnam,
         decade=factor(paste0(substr(year,1,3),"0s"),levels=c("1990s","2000s","2010s")))
codherring_dietsYears<-filter(bioenergDietSpecies,
                              species_comnam %in% c("Atlantic cod", "Atlantic herring") &
                                year>=1990)%>%
  mutate(Species=species_comnam,
         decade=factor(paste0(substr(year,1,3),"0s"),levels=c("1990s","2000s","2010s")))


ggplot(codherringCRtotal)+
  geom_ribbon(aes(Temp,ymin=R_Min+SDA_Min,ymax=R_Max+SDA_Max,
                  fill="R+SDA",color="R+SDA"),alpha=0.1,size=0.1)+
  geom_ribbon(aes(Temp,ymin=C_Max-Eg_Max,ymax=C_Min-Eg_Min,
                  fill="C",color="C"),alpha=0.1,size=0.1)+
  geom_line(aes(Temp,R_Mean+SDA_Mean,color="R+SDA"),size=2)+
  geom_line(aes(Temp,C_Mean-Eg_Mean,color="C"),size=2)+
  geom_errorbar(data=codherring_dietsYears,aes(x=T_clustermean,
                                               ymin=C_clustermean-C_SD,
                                               ymax=C_clustermean+C_SD),
                alpha=0.5)+
  geom_errorbarh(data=codherring_dietsYears,aes(xmin=T_clustermean-T_SD,
                                                xmax=T_clustermean+T_SD,
                                                y=C_clustermean),
                 alpha=0.5)+
  geom_point(data=codherring_dietsYears,aes(T_clustermean,C_clustermean,size=year))+
  facet_wrap(~Species,scales="free")+
  scale_fill_manual(values=c("dodgerblue2","firebrick2"),name="")+
  scale_color_manual(values=c("dodgerblue2","firebrick2"),name="")+
  theme_bw()+theme(legend.position="top")+
  ylab("Energy and materials (g/g/day)")+xlab("Temperature")

dummy <- codherring_dietsTrawl%>%
  group_by(Species)%>%
  summarise(Energy=range(C_mean,na.rm=T))%>%
  mutate(Temp=range(codherringCRtotal$Temp))

ggplot(data=codherringCRtotal)+
  geom_ribbon(aes(Temp,ymin=R_Min+SDA_Min,ymax=R_Max+SDA_Max,
                  fill="R+SDA",color="R+SDA"),alpha=0.1,size=0.1)+
  geom_ribbon(aes(Temp,ymin=C_Max-Eg_Max,ymax=C_Min-Eg_Min,
                  fill="C",color="C"),alpha=0.1,size=0.1)+
  geom_line(aes(Temp,R_Mean+SDA_Mean,color="R+SDA"),size=2)+
  geom_line(aes(Temp,C_Mean-Eg_Mean,color="C"),size=2)+
  scale_fill_manual(values=c("dodgerblue2","firebrick2"),name="")+
  scale_color_manual(values=c("dodgerblue2","firebrick2"),name="")+
  theme_bw()+theme(legend.position="top")+
  new_scale_color()+
  geom_point(data=codherring_dietsTrawl,aes(T_mean,C_mean,color=year),size=1.25,alpha=0.5)+
  geom_blank(data=dummy,aes(Temp,Energy))+
  scale_color_viridis_c(name="Year")+
  facet_wrap(Species~decade,scales="free")+
  ylab("Energy and materials (g/g/day)")+xlab("Temperature")




# Prey EDs ----------------------------------------------------------------

preyEDs<-read_csv("Diet Data/preyEDs_litReview.csv")





# Fleet Data --------------------------------------------------------------



rm(list=ls())

setwd("C:/Users/nh1087/OneDrive - USNH/Documents/NECC/Fleet Data/")


library(tidyverse)
library(ggplot2)
library(readr)
library(stringr)
library(lubridate)
library(sf)
library(grid)
library(gridExtra)

theme_set(theme_bw(base_size=25))




# Data --------------------------------------------------------------------

world<-map_data("world")
states<-map_data("state")

statAreas<-st_read("../Mapping Data/Statistical_Areas_2010_withNames.shp",quiet=T)
statAreas<-st_transform(statAreas,4326)

#US EEZ (and other boundaries)
eez<-st_read("../Mapping Data/USMaritimeLimitsNBoundaries.shp",quiet=T)
ATLeez<-filter(eez,REGION=="Atlantic Coast and Gulf of Mexico")
#plot(eez[,3])

#World EEZ
eez.world <- st_read("../Mapping Data/eez_v11.shp",quiet=T)
eez.usa2 <- eez.world[eez.world$TERRITORY1 == "United States", ]
#dat.eez.usa2 <- fortify.shape(dat.eez.usa2) # a 10298x30 dataframe


#My area of interest
statArea_crop<-function(df) {
  inLONG<-filter(df,between(cell_ll_lon,st_bbox(statAreas)$xmin,st_bbox(statAreas)$xmax))
  inLAT<-filter(inLONG,between(cell_ll_lat,st_bbox(statAreas)$ymin,st_bbox(statAreas)$ymax))
  return(inLAT)
}
st_bbox(statAreas)

fleetData<-read_csv("fishing-vessels-v2.csv")
USfleetData<-filter(fleetData,flag_registry=="USA")

files <- list.files(path = "mmsi-daily-csvs-10-v2-2020/", pattern = "*.csv", full.names = T)
fleet_daily_2020 <- sapply(files, read_csv, simplify=FALSE) %>% 
  bind_rows(.id = "id")
USfleet_daily_2020 <- right_join(fleet_daily_2020,USfleetData)


czmp<-st_read("../Mapping Data/CZMP_counties_2009.shp",quiet=T)
AtlanticStates<-c("Alabama","Connecticut","Delaware","DC","Florida","Georgia","Louisiana","Maine","Maryland","Massachusetts","Mississippi",
                  "New Hampshire","New Jersey","New York","North Carolina","Pennsylvania","Rhode Island","South Carolina","Texas","Virginia")
AtlanticStatesFIPS<-c("01","09","10","11","12","13","22","23","24","25","28","33","34","36","37","42","44","45","48","51")
ATLczmp<-filter(czmp,STATEFP %in% AtlanticStatesFIPS)

plot(ATLczmp[,1])


# Exploratory Analyses ----------------------------------------------------


fleet_01_01_2020<-read_csv("mmsi-daily-csvs-10-v2-2020/2020-01-01.csv")%>%
  mutate(fishing=ifelse(fishing_hours>0,"Fishing","Not Fishing"))

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_point(data=fleet_01_01_2020,aes(cell_ll_lon,cell_ll_lat,fill=fishing),shape=21)+
  scale_fill_viridis_d(option="C",name="",direction=-1)+
  theme(legend.position=c(0.85,0.15),legend.background = element_rect(fill="transparent"))+
  coord_sf(xlim=c(-77.8,st_bbox(statAreas)$xmax),
           ylim=c(35.2,st_bbox(statAreas)$ymax))+
  ggtitle("Fleet Positions: 01-01-2020")


fleet_01_07_2020<-read_csv("mmsi-daily-csvs-10-v2-2020/2020-07-01.csv")%>%
  mutate(fishing=ifelse(fishing_hours>0,"Fishing","Not Fishing"))

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_point(data=fleet_01_07_2020,aes(cell_ll_lon,cell_ll_lat,fill=fishing),shape=21)+
  scale_fill_viridis_d(option="C",name="",direction=-1)+
  theme(legend.position=c(0.85,0.15),legend.background = element_rect(fill="transparent"))+
  coord_sf(xlim=c(-77.8,st_bbox(statAreas)$xmax),
           ylim=c(35.2,st_bbox(statAreas)$ymax))+
  ggtitle("Fleet Positions: 01-07-2020")



fishingFleet_01_07_2020<-left_join(fleet_01_07_2020,fleetData)
fishingFleet_01_01_2020<-left_join(fleet_01_01_2020,fleetData)


table(fleetData$flag_registry)

boat_daily_2020<-filter(fleet_daily_2020,mmsi=="303672000")
boat_daily_2020$fishing<-ifelse(boat_daily_2020$fishing_hours>0,"Fishing","Not Fishing")

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite",color="black")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_point(data=boat_daily_2020,aes(cell_ll_lon,cell_ll_lat,fill=fishing),shape=21)+
  geom_point(data=filter(boat_daily_2020,fishing=="Fishing"),aes(cell_ll_lon,cell_ll_lat,fill=fishing),shape=21)+
  scale_fill_viridis_d(option="C",name="",direction=-1)+
  theme(legend.position=c(0.5,0.15),legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_line(color=rgb(0.5,0.5,0.5,alpha=0.7),linetype="dashed",size=0.5),
        panel.background = element_rect(fill="azure1"))+
  coord_sf(xlim=c(-120,25),
           ylim=c(17.5,60))+
  ggtitle("MMSI 303672000 Positions: 2020")



#All fishing points
fleet_daily_2020$fishing<-ifelse(fleet_daily_2020$fishing_hours>0,"Fishing","Not Fishing")


fleet_effort_2020<-fleet_daily_2020%>%
  group_by(cell_ll_lat,cell_ll_lon)%>%
  summarise(hours=sum(hours),
            fishing_hours=sum(fishing_hours))

crop_fleet_effort_2020<-statArea_crop(fleet_effort_2020)
facet_fleet_effort_2020<-pivot_longer(crop_fleet_effort_2020,cols=c("hours","fishing_hours"),
                                      names_to="Activity",values_to = "Hours")

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_point(data=crop_fleet_effort_2020,aes(cell_ll_lon,cell_ll_lat,color=hours),alpha=0.5,shape=15)+
  #geom_point(data=crop_fleet_effort_2020,aes(cell_ll_lon,cell_ll_lat,color=fishing_hours),alpha=0.95,shape=15)+
  scale_color_viridis_c(option="C",name="",trans="log",breaks=c(0.1,10,1000,100000),limits=c(0.01,500000))+
  theme(legend.position="none")+
  coord_sf(xlim=c(-77.8,st_bbox(statAreas)$xmax),
           ylim=c(35.2,st_bbox(statAreas)$ymax))+
  ggtitle("All fleet positions: 2020")+
  ylab("Latitude")+xlab("Longitude")

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  #geom_point(data=crop_fleet_effort_2020,aes(cell_ll_lon,cell_ll_lat,color=hours),alpha=0.5,shape=15)+
  geom_point(data=crop_fleet_effort_2020,aes(cell_ll_lon,cell_ll_lat,color=fishing_hours),alpha=0.95,shape=15)+
  scale_color_viridis_c(option="C",name="Hours",trans="log",
                        breaks=c(0.1,10,1000,100000),labels=c("0.1","10","1,000","100,000"),limits=c(0.01,500000))+
  coord_sf(xlim=c(-77.8,st_bbox(statAreas)$xmax),
           ylim=c(35.2,st_bbox(statAreas)$ymax))+
  ggtitle("All fleet fishing positions: 2020")

activity.labs<-c("Positions","Fishing Positions")
names(activity.labs)<-c("hours","fishing_hours")
ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_point(data=facet_fleet_effort_2020,
             aes(cell_ll_lon,cell_ll_lat,color=Hours),alpha=0.95,shape=15)+
  scale_color_viridis_c(option="C",name="Hours",trans="log",
                        breaks=c(0.1,10,1000,100000),labels=c("0.1","10","1,000","100,000"),limits=c(0.01,500000))+
  coord_sf(xlim=c(-77.8,st_bbox(statAreas)$xmax),
           ylim=c(35.2,st_bbox(statAreas)$ymax))+
  xlab("Longitude")+ylab("Latitude")+
  facet_wrap(~Activity,nrow=1,labeller=labeller(Activity=activity.labs))+
  theme(strip.background = element_blank(), strip.placement = "outside")

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite",color="grey50",size=0.25)+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50",size=0.25)+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_sf(data=ATLczmp)+
  geom_sf(data=filter(ATLeez,EEZ=="1")[,3])+
  geom_point(data=USfleet_daily_2020,
             aes(cell_ll_lon,cell_ll_lat,color=fishing_hours),alpha=0.95,shape=15)+
  scale_color_viridis_c(option="C",name="Hours",trans="log",
                        breaks=c(0.1,10,1000,100000),labels=c("0.1","10","1,000","100,000"),limits=c(0.01,500000))+
  coord_sf(xlim=c(-97.8,-42.5),
           ylim=c(21,51.5))+
  xlab("Longitude")+ylab("Latitude")
#facet_wrap(~Activity,nrow=1,labeller=labeller(Activity=activity.labs))+
#theme(strip.background = element_blank(), strip.placement = "outside")


unique(eez$CZ)

# Pulling out all the individuals and finding ports -----------------------

#Points
#what fishing is in the EEZ and what's in the state waters
USsummary<-USfleet_daily_2020%>%
  group_by(cell_ll_lon,cell_ll_lat)%>%
  summarise(Hours=sum(hours),
            Fishing_Hours=sum(fishing_hours))%>%ungroup()%>%
  filter(!is.na(cell_ll_lat))%>%
  st_as_sf(coords=c("cell_ll_lon","cell_ll_lat"))
st_crs(USsummary)<-st_crs(eez.usa2)
ATLczmp<-st_transform(ATLczmp,crs=st_crs(USsummary))

fishing_in_eez <- st_join(USsummary, eez.usa2, join = st_within)
fishing_in_states <- st_join(USsummary,ATLczmp,join = st_within)

fishing_in<-st_join(fishing_in_eez,fishing_in_states)%>%
  mutate(location=paste(gsub("[0-9]+","State Waters",STATEFP),TERRITORY1,sep="-"),
         location=case_when(grepl("State Waters",location)~"State Waters",
                            location=="NA-United States"~"Federal Waters",
                            location=="NA-NA"~"International Waters"))
unique(fishing_in$location)
ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite",color="grey50",size=0.25)+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50",size=0.25)+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_sf(data=fishing_in,aes(color=location))+
  geom_sf(data=eez.usa2,fill="transparent")+
  #geom_point(data=USsummary,
  #           aes(cell_ll_lon,cell_ll_lat,color=Fishing_Hours),alpha=0.95,shape=15)+
  #scale_color_viridis_c(option="C",name="Hours",trans="log",
  #                      breaks=c(0.1,10,1000,100000),labels=c("0.1","10","1,000","100,000"),limits=c(0.01,500000))+
  coord_sf(xlim=c(-97.8,-42.5),
           ylim=c(21,51.5))+
  xlab("Longitude")+ylab("Latitude")


#State waters
stateWater<-fortify(ATLczmp)
#US EEZ
class(eez.usa2)

#what fishing is in the EEZ and what's in the state waters
USsummary<-USfleet_daily_2020%>%
  group_by(cell_ll_lon,cell_ll_lat)%>%
  summarise(Hours=sum(hours),
            Fishing_Hours=sum(fishing_hours),
            stateWaters=point.in.polygon(cell_ll_lon,cell_ll_lat))
ATLeezPoly<-st_cast(ATLeez$geometry[[1]],"MULTIPOLYGON")

ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite",color="grey50",size=0.25)+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50",size=0.25)+
  geom_sf(data=statAreas,fill="transparent",col="black")+
  geom_sf(data=ATLczmp)+
  geom_sf(data=eez.usa2)+
  #geom_point(data=USsummary,
  #           aes(cell_ll_lon,cell_ll_lat,color=Fishing_Hours),alpha=0.95,shape=15)+
  #scale_color_viridis_c(option="C",name="Hours",trans="log",
  #                      breaks=c(0.1,10,1000,100000),labels=c("0.1","10","1,000","100,000"),limits=c(0.01,500000))+
  coord_sf(xlim=c(-97.8,-42.5),
           ylim=c(21,51.5))+
  xlab("Longitude")+ylab("Latitude")
#facet_wrap(~Activity,nrow=1,labeller=labeller(Activity=activity.labs))+
#theme(strip.background = element_blank(), strip.placement = "outside")





fleetDays<-fleet_daily_2020%>%
  dplyr::select(date,mmsi)%>%
  distinct()

check<-filter(fleet_daily_2020,is.na(date))

eezPolygon<-raster::bind(states,eez)


USfleet_daily_2020<-fleet_daily_2020%>%
  filter()



# Landings ----------------------------------------------------------------

landings<-read_csv("foss_landings.csv")



# Using Global Fishing Watch API ------------------------------------------

devtools::install_github("GlobalFishingWatch/gfwr") #not working



# Bottom ------------------------------------------------------------------




