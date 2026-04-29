
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
library(colorspace)
library(R2jags)
library(cowplot)




#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}


seasonPal2<-c("deepskyblue4","yellowgreen","goldenrod1","orange3") #BEST

#Data loading
load("NF.prey19.RData")
allprey<-NF.prey19
rm(NF.prey19)
load("prey19.RData")
load("googleMap.zoom6.eastCoast.R")
vulScores<-read.csv("../Hareetal_climateScores.csv")


#ALL LEVELS OF YEAR AND SEASON
ALLyearseasons<-factor(levels=c(paste(sort(rep(seq(1973,2019),4)),rep(c("Winter","Spring","Summer","Fall"),length(seq(1973,2019))))))
#levels(ALLyearseasons)

#The names of key 9 species, just to have them nice and handy
keyspecies_sci<-unique(str_to_sentence(prey19$pdscinam))
keyspecies_com<-unique(str_to_title(prey19$pdcomnam))
keyspecies_svspp<-unique(str_pad(as.character(prey19$svspp),width=3,pad="0",side="left"))

#And then all the species in the food habits database
species_sci<-unique(str_to_sentence(allprey$pdscinam))
species_com<-unique(str_to_title(allprey$pdcomnam))
species_svspp<-unique(str_pad(as.character(allprey$svspp),width=3,pad="0",side="left"))


# Species SVSPP Codes -----------------------------------------------------

spec<-read.csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv")%>%
  dplyr::select(LOGGED_SPECIES_NAME,SVSPP)%>%distinct()%>%
  filter(SVSPP>="001")%>%
  mutate(LOGGED_SPECIES_NAME=str_to_sentence(LOGGED_SPECIES_NAME))
spec<-spec[!duplicated(spec$SVSPP),]


# Important df manipulations ----------------------------------------------

allprey<-allprey%>%
  mutate(gensci=ifelse(pynam=="PRIONOTUS ALATUS"|pynam=="STERNOPTYCHIDAE","FISH",as.character(gensci)),
         analsci=ifelse(pynam=="PRIONOTUS ALATUS","TRIGLIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="PRIONOTUS ALATUS",as.character(pynam),as.character(collsci)),
         analsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(collsci)),
         year=ifelse(is.na(year),substr(cruise6,1,4),year),year=as.numeric(year))%>%
  filter(pdlen<300) #There is a Summer Flounder that was listed at 300 cm, which is unbelievable, so it will be dropped



uniquePrey.all<-allprey%>%
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

allprey<-allprey%>%
  mutate(INgen=ifelse(str_to_sentence(gensci) %in% gl_preycats$matchingCats,1,0),
         INanal=ifelse(str_to_sentence(analsci) %in% gl_preycats$matchingCats,1,0),
         INcoll=ifelse(str_to_sentence(collsci) %in% gl_preycats$matchingCats,1,0),
         INpy=ifelse(str_to_sentence(pynam) %in% gl_preycats$matchingCats,1,0))
allprey<-allprey%>%
  mutate(INnum=rowSums(allprey[,c("INgen","INanal","INcoll","INpy")]),
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
  filter(svspp%in%species_svspp)%>%
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

allprey<-allprey%>%
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
                            TRUE ~ "S"),
         geoarea=factor(geoarea,levels=c("MAB",
                                         "SNE",
                                         "GB",
                                         "GoM",
                                         "ScS")))%>%
  select(-c(species_scinam:xlarge_max))


strata_key <- list(
  "Georges Bank"          = as.character(13:23),
  "Gulf of Maine"         = as.character(24:40),
  "Southern New England"  = stringr::str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
  "Mid-Atlantic Bight"    = as.character(61:76))
preyGMRI<-left_join(allprey,sumGMRItrawls,by=c("svspp","sizecat2"="sizecat","id"))%>% #merge by sizecat2 because this is defined from Garrison and Link in the same way as sizecat is defined in the GMRI, clearly the diet df is defining sizecat some other way (especially if not entirely for spiny dogfish, sexual dimorphism maybe??)
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


allTrawls<-read_csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVSTA.csv")
statAreas<-st_read("../Mapping Data/Statistical_Areas_2010_withNames.shp",quiet=T)
statAreas<-filter(statAreas,Id%in%allTrawls$AREA)
statAreas<-st_transform(statAreas,4326)

survAreas<-st_read("../Mapping Data/BTS_Strata.shp")
survAreas<-st_transform(survAreas,4326)

#Using shapefiles to build a clean map
world<-map_data("world")
states<-map_data("state")


dMap<-ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="antiquewhite")+
  geom_polygon(data=states,aes(long,lat,group=group),fill="transparent",color="grey50")+
  geom_spatial_point(data=preyGMRI_filter,
                     aes(-1*declon,declat,color=geoarea),
                     crs=4326,size=0.5)+
  geom_sf(data=survAreas,fill="transparent",col="black")+
  scale_color_discrete_sequential(palette="Hawaii",name="Geographic\nRegion",l2=80)+
  guides(color=guide_legend(override.aes = list(size=9)))+
  annotation_scale(location="br",width_hint=0.4,height=unit(0.125,"in"),
                   pad_x=unit(0.25,"in"),pad_y=unit(0.4,"in"),text_cex=1.2)+
  geom_rect(aes(xmin=-63.35,xmax=-62.5,ymin=35.8,ymax=36.2),fill="azure1")+
  annotation_north_arrow(location="br",pad_x=unit(0.5,"in"),pad_y=unit(0.75,"in"),
                         style = north_arrow_fancy_orienteering(text_size=15,line_width=1.25),
                         height=unit(1,"in"),width=unit(0.9,"in"))+
  coord_sf(xlim=c(-77.8,st_bbox(statAreas)$xmax),
           ylim=c(35.2,st_bbox(statAreas)$ymax))+
  theme(legend.position=c(0.225,0.8),legend.background=element_rect(color="black"),
        panel.grid.major = element_line(color=rgb(0.5,0.5,0.5,alpha=0.7),linetype="dashed",size=0.5),
        panel.background = element_rect(fill="azure1"))+
  xlab("Longitude")+ylab("Latitude")

iMap<-ggplot()+
  geom_polygon(data=world,aes(long,lat,group=group),fill="grey25",color="grey25")+
  geom_sf(data=st_as_sfc(st_bbox(statAreas)),color="red",fill="transparent",size=2)+
  coord_sf(xlim=c(-180,0),ylim=c(0,70))+
  theme_void()+
  theme(panel.border = element_rect(color="black",fill="transparent"))

ggdraw()+
  draw_plot(dMap)+
  draw_plot(iMap,x=0.725,y=0.275,width=0.225,height=0.225)




