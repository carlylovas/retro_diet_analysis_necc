### CSL: This was a test run of this section to see if it generates any of the manuscript figures. It does not. 


#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}


seasonPal2<-c("deepskyblue4","yellowgreen","goldenrod1","orange3") #BEST

#Data loading
load(here("data","processed","prey19.RData"))
load(here("data","processed","pylen19.RData"))
load(here("data","shapefiles","googleMap.zoom6.eastCoast.R"))
vulScores<-read.csv(here("data","raw","Hareetal_climateScores.csv"))

# Species SVSPP Codes -----------------------------------------------------

spec<-read.csv(here("data","raw","22560_UNION_FSCS_SVCAT.csv")) %>%
  dplyr::select(LOGGED_SPECIES_NAME,SVSPP)%>%distinct()%>%
  filter(SVSPP>="001")%>%
  mutate(LOGGED_SPECIES_NAME=str_to_sentence(LOGGED_SPECIES_NAME))
spec<-spec[!duplicated(spec$SVSPP),]


# Important df manipulations ----------------------------------------------
# 
# prey19<-prey19%>%
#   mutate(gensci=ifelse(pynam=="PRIONOTUS ALATUS"|pynam=="STERNOPTYCHIDAE","FISH",as.character(gensci)),
#          analsci=ifelse(pynam=="PRIONOTUS ALATUS","TRIGLIDAE",as.character(analsci)),
#          collsci=ifelse(pynam=="PRIONOTUS ALATUS",as.character(pynam),as.character(collsci)),
#          analsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(analsci)),
#          collsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(collsci)),
#          year=ifelse(is.na(year),substr(cruise6,1,4),year),year=as.numeric(year))%>%
#   filter(pdlen<300)


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



#The names of my species, just to have them nice and handy
myspecies_sci<-unique(str_to_sentence(prey19$pdscinam))
myspecies_com<-unique(str_to_title(prey19$pdcomnam))
myspecies_svspp<-unique(str_pad(as.character(prey19$svspp),width=3,pad="0",side="left"))

# Fish Size Distributions -------------------------------------------------
sizeClasses<-read.csv(here("data","raw","GarrisonLink_predSizeClasses.csv"))%>%
  left_join(spec,by=c("species_comnam"="LOGGED_SPECIES_NAME"))

fallSizes<-read_csv(here("data","raw","22560_UNION_FSCS_SVLEN.csv"),
                    col_types = "cccccccnnn")
springSizes<-read_csv(here("data","raw","22561_UNION_FSCS_SVLEN.csv"),
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
  reframe(p=sum(EXPNUMLEN)/N)%>%distinct()

allGMRItrawls_clean<-read_csv(here("data","processed","NMFS_survdat_gmri_tidy.csv"),
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
  distinct() # %>%expand(geoarea,comname)


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
  reframe(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(sizeabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(drop=(max(pred)-min(pred))*100,
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
  reframe(drop=(max(pred)-min(pred))*100,
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
  reframe(drop=(max(pred)-min(pred)),
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
  reframe(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(Wk=sum(trawlabundance*meanWt)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
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
  reframe(PyPd=sum(trawlabundance*meanPyPd)/sum(trawlabundance),
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
  reframe(dietAbundance=n(),
            propDiet=dietAbundance/abundance)%>%distinct()
nrow(test[test$propDiet==1,])/nrow(test)*100
summary(test$propDiet)



