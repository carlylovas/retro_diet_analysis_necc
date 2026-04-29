

setwd("C:/Users/nh1087/OneDrive - USNH/Documents/NECC/Diet Data/")


library(tidyverse)
library(ggplot2)
library(readr)
library(stringr)
library(lubridate)

theme_set(theme_bw(base_size=25))

#Data loading
load("prey19.RData")
prey19<-prey19%>%
  mutate(gensci=ifelse(pynam=="PRIONOTUS ALATUS"|pynam=="STERNOPTYCHIDAE","FISH",as.character(gensci)),
         analsci=ifelse(pynam=="PRIONOTUS ALATUS","TRIGLIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="PRIONOTUS ALATUS",as.character(pynam),as.character(collsci)),
         analsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(collsci)),
         year=ifelse(is.na(year),substr(cruise6,1,4),year),year=as.numeric(year))%>%
  filter(pdlen<300) #There is a Summer Flounder that was listed at 300 cm, which is unbelievable, so it will be dropped


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
  dplyr::select(-c(pynam:dietID))%>%
  distinct()%>%
  mutate(cruise6=as.character(cruise6),
         station=str_pad(station,width=4,pad="0",side="left"),
         pdsex=as.character(pdsex),
         pdid=str_pad(pdid,width=6,pad="0",side="left"))

nrow(indGMRI_filter)/n_distinct(preyGMRI$dietID)*100

allGeoCom<-dplyr::select(preyGMRI_filter,geoarea,comname)%>%
  distinct()%>%expand(geoarea,comname)



# Connecting with the full trawl data -------------------------------------

fallT<-read_csv("../Trawl Data/NMFS Trawls/2022 Redownload/Fall/22560_UNION_FSCS_SVCAT.csv",col_types = "ccccccccnn")
springT<-read_csv("../Trawl Data/NMFS Trawls/2022 Redownload/Spring/22561_UNION_FSCS_SVCAT.csv",col_types = "ccccccccnn")
fullT<-bind_rows(fallT,springT)%>%
  mutate(svspp=str_pad(SVSPP,width=3,side="left",pad="0"),
         station=substr(STATION,2,4),
         stratum=substr(STRATUM,2,5),
         id=paste0(CRUISE6,
                   station,
                   stratum))%>%
  dplyr::select(id,svspp,EXPCATCHNUM,EXPCATCHWT)%>%
  group_by(id,svspp)%>%
  summarise(abundance=sum(EXPCATCHNUM),
            wt=sum(EXPCATCHWT))


preyFull<-left_join(prey19,fullT)%>%
  filter(!is.na(abundance))


smithTest<-read_csv("SmithTest/final2.csv")%>%
  mutate(Wk_percent=Wk*100,
         diff_percents=Wk_percent-relmsw,
         relative_error=abs(Wk_percent-relmsw)/relmsw)

Ntrawls<-preyGMRI_filter%>%
  group_by(comname,collsci)%>%
  summarise(N=n_distinct(id))

smithTest<-left_join(smithTest,Ntrawls)

check<-filter(smithTest,relative_error>0.5)

ggplot(smithTest,aes(relmsw,relative_error))+geom_point()+
  ylim(0,700)+ylab("Relative Error")+xlab("Relative Mean Mass")
ggplot(smithTest,aes(N,relative_error))+geom_point()+
  ylim(0,700)+xlim(0,3100)+ylab("Relative Error")+xlab("Number of Diets with Item")
ggplot(smithTest,aes(cv,relative_error))+geom_point()+
  ylim(0,700)+xlim(0,1.01)+ylab("Relative Error")+xlab("CV")

check<-smithTest%>%group_by(comname)%>%summarise(N=sum(relmsw,na.rm=T))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0
smithCheck<-smithTest%>%
  group_by(comname)%>%
  summarise(total=sum(meansw,na.rm=T),
            max=max(meansw,na.rm=T))


#What I think he might actually be doing (just the mean mass weighted by trawl)
smithDietSum<-preyFull%>%
  group_by(pdcomnam)%>%
  mutate(sumabundance=sum(abundance))%>%
  group_by(pdcomnam,collsci)%>%
  summarise(meansw=sum(abundance*pyamtw)/sumabundance)%>%distinct()
  
check<-smithDietSum%>%group_by(pdcomnam)%>%summarise(N=sum(meansw,na.rm=T))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0
smithCheck2<-smithDietSum%>%
  group_by(pdcomnam)%>%
  summarise(total=sum(meansw,na.rm=T),
            max=max(meansw,na.rm=T))



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
  group_by(id,pdscinam,collsci)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets))%>%
  dplyr::select(id,trawlabundance,pdscinam,comname,collsci,totalwt,totalv,nDiets,qikw,qikv,pik)%>% 
  distinct()%>% #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.
  ungroup()

nDiets<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(comname,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(comname)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(comname,id,trawlabundance)%>%distinct()%>%
  group_by(comname)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(comname,id,trawlabundance)%>%distinct()%>%
  group_by(comname)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietSum<-trawlDiets%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(pdscinam,comname,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,collsci)%>%
  summarise(Wk=sum(trawlabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(trawlabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(trawlabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,trawlabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       trawlabundance,pik, Fk))%>%distinct()

check<-trawlDietSum%>%group_by(pdscinam)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0



#Size categories and years
#Need to have whole diet totals, and whole trawl measures
trawlDiets_size<-preyGMRI_filter%>%
  group_by(id,comname,sizecat2)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv))%>%
  group_by(id,pdscinam,sizecat2,collsci)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets),
         year=factor(year),
         comname=factor(comname),
         sizecat=factor(sizecat2,levels=c("XS","S","M","L","XL")))%>%
  ungroup()%>%
  dplyr::select(id,sizeabundance,year,pdscinam,comname,sizecat,collsci,totalwt,totalv,nDiets,qikw,qikv,pik)%>% 
  distinct() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.


nDiets<-trawlDiets_size%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(year,comname,sizecat,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(year,comname,sizecat,.drop=F)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun<-trawlDiets_size%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(year,comname,sizecat,id,sizeabundance)%>%distinct()%>%
  group_by(year,comname,sizecat,.drop=F)%>%summarise(sumAbun=sum(sizeabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDiets_size%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(year,comname,sizecat,id,sizeabundance)%>%distinct()%>%
  group_by(year,comname,sizecat,.drop=F)%>%summarise(sumAbun_nEmpty=sum(sizeabundance))


trawlDietSum_size<-trawlDiets_size%>%
  group_by(year)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(year,pdscinam,comname,sizecat,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,collsci)%>%
  summarise(Wk=sum(sizeabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(sizeabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(sizeabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizeabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizeabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       sizeabundance,pik, Fk))%>%distinct()

check<-trawlDietSum_size%>%group_by(year,pdscinam,sizecat)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0

write.csv(trawlDietSum,"eachPredator_summary.csv",row.names=F)
write.csv(trawlDietSum_size,"PredSizeYear_summary.csv",row.names=F)




# Absolute Mass--What Brian Really Did ------------------------------------



#Need to have whole diet totals, and whole trawl measures
trawlDietsAbs<-preyGMRI_filter%>%
  group_by(id,comname,collsci)%>%
  mutate(qikw=mean(pyamtw))%>%
  dplyr::select(id,trawlabundance,comname,comname,collsci,qikw)%>% 
  distinct()%>% #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.
  ungroup()

sumAbun<-trawlDietsAbs%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(comname,id,trawlabundance)%>%distinct()%>%
  group_by(comname)%>%summarise(sumAbun=sum(trawlabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDietsAbs%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(qikw!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(comname,id,trawlabundance)%>%distinct()%>%
  group_by(comname)%>%summarise(sumAbun_nEmpty=sum(trawlabundance))


trawlDietAbsSum<-preyGMRI_filter%>%
  group_by(comname)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(comname,nTrawls,sumAbun,sumAbun_nEmpty,collsci)%>%
  summarise(Wk=sum(trawlabundance*pyamtw)/sumAbun, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            WkSD=clusterSD(nTrawls,sumAbun,trawlabundance,pyamtw,Wk))%>%distinct()

check<-trawlDietAbsSum%>%group_by(comname)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0



#Size categories and years
#Need to have whole diet totals, and whole trawl measures
trawlDiets_size<-preyGMRI_filter%>%
  group_by(id,comname,sizecat2)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv))%>%
  group_by(id,pdscinam,sizecat2,collsci)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(nDiets==0,0,n_distinct(dietID)/nDiets),
         year=factor(year),
         comname=factor(comname),
         sizecat=factor(sizecat2,levels=c("XS","S","M","L","XL")))%>%
  ungroup()%>%
  dplyr::select(id,sizeabundance,year,pdscinam,comname,sizecat,collsci,totalwt,totalv,nDiets,qikw,qikv,pik)%>% 
  distinct() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.


nDiets<-trawlDiets_size%>% #Want to know how many diets were analyzed for a GROUP (in this case each region for a season (and year))
  ungroup()%>%
  dplyr::select(year,comname,sizecat,id,nDiets)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(year,comname,sizecat,.drop=F)%>%summarise(tnDiets=sum(nDiets))#Calculate the number of diets
sumAbun<-trawlDiets_size%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(year,comname,sizecat,id,sizeabundance)%>%distinct()%>%
  group_by(year,comname,sizecat,.drop=F)%>%summarise(sumAbun=sum(sizeabundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDiets_size%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(year,comname,sizecat,id,sizeabundance)%>%distinct()%>%
  group_by(year,comname,sizecat,.drop=F)%>%summarise(sumAbun_nEmpty=sum(sizeabundance))


trawlDietSum_size<-trawlDiets_size%>%
  group_by(year)%>%
  mutate(nTrawls=n_distinct(id))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(year,pdscinam,comname,sizecat,nTrawls,tnDiets,sumAbun,sumAbun_nEmpty,collsci)%>%
  summarise(Wk=sum(sizeabundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(sizeabundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(sizeabundance*pik) /sumAbun,        #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizeabundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,sizeabundance,qikv,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,       sizeabundance,pik, Fk))%>%distinct()

check<-trawlDietSum_size%>%group_by(year,pdscinam,sizecat)%>%summarise(N=sum(Wk))
table(check$N) #You have all your means correct because they sum to 1 or they're all 0

write.csv(trawlDietSum,"eachPredator_summary.csv",row.names=F)
write.csv(trawlDietSum_size,"PredSizeYear_summary.csv",row.names=F)




#Diet priorities
test<-read.table("dietPriorities.txt",sep=" ")
colnames(test)<-c("Genus","Species",paste(floor(seq(1973,2008,by=0.5)),c("Sp","Fa"),sep="_"))
write.csv(test,"dietPriorities.csv",row.names=F)

