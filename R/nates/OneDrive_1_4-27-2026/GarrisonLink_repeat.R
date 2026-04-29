
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
library(dendextend)


#Personalization
theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')

yDate<-function(x) {
  print(paste0("Non-leap year: ",month(as.Date("2019-01-01")+x-1,label=T),"-",day(as.Date("2019-01-01")+x-1)))
  print(paste0("Leap year: ",month(as.Date("2020-01-01")+x-1,label=T),"-",day(as.Date("2020-01-01")+x-1)))
}

unique(as.character(prey19$pdcomnam))

seasonPal<-c("steelblue4","goldenrod1","deeppink3","sienna2")
seasonPal2<-c("deepskyblue4","yellowgreen","lightgoldenrod1","orange3") #BEST
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

#Size classes used by Garrison and Link
sizeClasses<-read.csv("GarrisonLink_predSizeClasses.csv") 
#These are replicated in the sizecats column already present, except that the small and largest cats extend to the smallest and largest individuals
prey19%>%group_by(pdcomnam,sizecat)%>%summarise(min=min(pdlen),max=max(pdlen))




# Important df manipulations ----------------------------------------------

prey19<-prey19%>%
  mutate(gensci=ifelse(pynam=="PRIONOTUS ALATUS"|pynam=="STERNOPTYCHIDAE","FISH",as.character(gensci)),
         analsci=ifelse(pynam=="PRIONOTUS ALATUS","TRIGLIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="PRIONOTUS ALATUS",as.character(pynam),as.character(collsci)),
         analsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(collsci)),
         year=ifelse(is.na(year),substr(cruise6,1,4),year),year=as.numeric(year))%>%
  filter(pdlen<300) #There is a Summer Flounder that was listed at 300 cm, which is unbelievable, so it will be dropped


check<-filter(test,is.na(sizeClass))
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
                                                ifelse(INnum==0&pynam=="EMPTY","Empty",NA))))))))))),
         gl_prey=ifelse(pynam %in% c("UROPHYCIS CHUSS","UROPHYCIS TENUIS","UROPHYCIS REGIA"),"Other hakes", #STEP 2
                    ifelse(analsci %in% c("GADIDAE","BREGMACEROTIDAE","EUCLICHTHYIDAE","LOTIDAE","MACROURIDAE",
                                          "MELANONIDAE","MERLUCCIIDAE","MORIDAE","MURAENOLEPIDIDAE","PHYCIDAE") & 
                                 is.na(gl_prey),"Gadiformes", #STEP 3a
                      ifelse(analsci %in% c("PLEURONECTIDAE","PSETTODIDAE","CITHARIDAE","SCOPHTHALMIDAE","PARALICHTHYIDAE",
                                            "BOTHIDAE","PARALICHTHODIDAE","POECILOPSETTIDAE","RHOMBOSOLEIDAE",
                                            "ACHIROPSETTIDAE","SAMARIDAE","ACHIRIDAE","SOLEIDAE","CYNOGLOSSIDAE") & 
                                   is.na(gl_prey), "Pleuronectiformes", #STEP 3b
                          ifelse((gensci=="FISH" & is.na(gl_prey))|pynam=="FISH","Unidentified fish", #STEP 4
                             ifelse((gensci %in% c("UROCHORDATA","BRACHIOPODA","BRYOZOA","CHAETOGNATHA","PORIFERA") | 
                                     collsci %in% c("ARTHROPODA","INSECTA","HEMICHORDATA","LIMULUS POLYPHEMUS",
                                                    "APHRODITIDAE","OLIGOCHAETA","HIRUDENEA","PYCNOGONIDA")) & 
                                     is.na(gl_prey),"Other invertebrates", #STEP 5
                                ifelse(analsci %in% c("CEPHALOCHORDATA"),"Other", #STEP 6
                                   ifelse(collsci %in% c("OSTRACODA","CUMACEA","STOMATOPODA","PENAEIDAE"),
                                          "Crustacean shrimp", #STEP 7
                                      ifelse(analsci %in% c("CIRRIPEDIA","COPEPODA") | 
                                             pynam %in% c("DECAPODA","DECAPODA EGGS","DECAPODA LARVAE"),"Crustacea", #STEP 8
                                        ifelse(collsci %in% c("HOMARUS AMERICANUS","CALLINECTES SAPIDUS",
                                                              "DECAPODA LARVAE","SCYLLARIDAE") &
                                               is.na(gl_prey), "Decapoda crab", #STEP 9
                                           ifelse(analsci=="EUPHAUSIACEA","Euphausiidae",gl_prey)))))))))), #STEP 10
         gl_prey=factor(gl_prey,levels=c(gl_preycats$matchingCats,"Empty")))


#Future work with gl_prey, to make things easier
ref_glprey<-prey19%>%
  dplyr::select(gl_prey,pynam,collsci,analsci,gensci)%>%
  distinct()


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

#reading in prey trawls df from "connecting_food+trawls.R" where prey19 is merged with GMRI's clean data for trawls
preyTrawls<-read_csv("SUMMARIZED/merged_prey_GMRItrawls.csv")
preyTrawls<-preyTrawls%>%
  mutate(gensci=ifelse(pynam=="PRIONOTUS ALATUS"|pynam=="STERNOPTYCHIDAE","FISH",as.character(gensci)),
         analsci=ifelse(pynam=="PRIONOTUS ALATUS","TRIGLIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="PRIONOTUS ALATUS",as.character(pynam),as.character(collsci)),
         analsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(analsci)),
         collsci=ifelse(pynam=="STERNOPTYCHIDAE","STERNOPTYCHIDAE",as.character(collsci)))%>%
  filter(pdlen<300)%>% #There is a Summer Flounder that was listed at 300 cm, which is unbelievable, so it will be dropped
  left_join(dplyr::select(ref_glprey,pynam,gl_prey)%>%distinct()) #add in the gl prey cats


#Need to have whole diet totals, and whole trawl measures
trawlDiets<-preyTrawls%>%
  mutate(pdid=str_pad(pdid,width=6,pad="0",side="left"),
         dietID=paste(svspp,cruise6,station,pdsex,pdid,str_pad(pdlen,width=3,pad="0",side="left"),sep="0"))%>%
  group_by(ID,pdscinam,sizecat)%>%
  mutate(totalwt=sum(pyamtw),
         totalv=sum(pyamtv),
         totalN=n_distinct(dietID))%>%
  group_by(ID,pdscinam,sizecat,gl_prey)%>%
  mutate(qikw=ifelse(totalwt==0,0,sum(pyamtw)/totalwt),
         qikv=ifelse(totalv==0,0,sum(pyamtv)/totalv),
         pik=ifelse(totalN==0,0,n_distinct(dietID)/totalN))%>%
  dplyr::select(#Can add in variables here if useful for grouping, like lat-long or temps, they're for grouping later
                ID,abundance,biomass_kg,pdscinam,comname,sizecat,gl_prey,totalwt,totalv,totalN,qikw,qikv,pik)%>% 
  distinct() #Doing "summarise" the old fashioned way, makes it easier to edit later when I want to do different vars. Just have to drop individual columns, like dietID, pdlen, etc.

#if grouping, add those groups in here both at select and group_by
nDiets<-trawlDiets%>% #Want to know how many diets were analyzed for a GROUP (in this case the whole pred sizecat)
  ungroup()%>%
  dplyr::select(pdscinam,sizecat,ID,totalN)%>%
  distinct()%>% #Reduce replicated N values to the unique trawl counts, only want to count each trawl once
  group_by(pdscinam,sizecat)%>%summarise(nDiets=sum(totalN))#Calculate the number of diets
sumAbun<-trawlDiets%>% #Same as above, except for the abundance for that GROUP
  ungroup()%>%
  dplyr::select(pdscinam,sizecat,ID,abundance)%>%distinct()%>%
  group_by(pdscinam,sizecat)%>%summarise(sumAbun=sum(abundance)) #Total abundance for a GROUP from each of the trawl IDs
sumAbun_nEmpty<-trawlDiets%>% #Same as above except only for those trawls where at least some diets had mass
  ungroup()%>%
  filter(totalwt!=0)%>% #Whether a true empty or appearing empty, they need to be removed in calculating the total abundance to divide by because they are not contributing mass so they aren't considered in mean
  dplyr::select(pdscinam,sizecat,ID,abundance)%>%distinct()%>%
  group_by(pdscinam,sizecat)%>%summarise(sumAbun_nEmpty=sum(abundance))


trawlDietSum<-trawlDiets%>%
  ungroup()%>%
  #group_by(geoarea,est_year,season)%>% #NOT grouping, but add if you are
  mutate(nTrawls=n_distinct(ID))%>%
  left_join(nDiets)%>%
  left_join(sumAbun)%>%
  left_join(sumAbun_nEmpty)%>%
  group_by(pdscinam,comname,sizecat,nTrawls,nDiets,sumAbun,sumAbun_nEmpty,gl_prey)%>%
  summarise(Wk=sum(abundance*qikw)/sumAbun_nEmpty, #Need to exclude the trawl abundance counts when the whole trawl was "empty" or else you're "reserving space" in diet prop
            Vk=sum(abundance*qikv)/sumAbun_nEmpty, #Same
            Fk=sum(abundance*pik)/sumAbun, #Frequency of occurrence, doesn't need to remove empties
            WkSD=clusterSD(nTrawls,sumAbun_nEmpty,abundance,qikw,Wk), #Calculate the SD for the cluster means, see above for function
            VkSD=clusterSD(nTrawls,sumAbun_nEmpty,abundance,qikw,Vk),
            FkSD=clusterSD(nTrawls,sumAbun,abundance,pik,Fk))%>%distinct()
levels(trawlDietSum$gl_prey)<-c(levels(trawlDietSum$gl_prey),"Unobserved")
trawlDietSum[is.na(trawlDietSum$gl_prey),"gl_prey"]<-"Unobserved"
trawlDietSum[is.na(trawlDietSum)]<-0 #Convert NAs to 0, because the contribution of empties is 0

check<-trawlDietSum%>%group_by(pdscinam,sizecat)%>%summarise(N=sum(Wk))
sum(check$N<0&check$N>0) #When this is 0, then you have all your means correct because they sum to 1 (or they're all empties)

#Saving these, so they can be better used within other analyses
#write.csv(trawlDiets,"trawl_speciesClusterCompositions.csv",row.names=F)
#write.csv(trawlDietSum,"geoareayearseason_speciesClusterCompositions.csv",row.names=F)





#####
#SHOULD THESE BE MEANS ACROSS SIZECAT, THE ABUNDANCES AREN'T THEY'RE JUST AT SPECIES
########



#Following through to get the "real" overlap matrix
props1_gl<-trawlDietSum[,c("comname","sizecat","gl_prey","Wk")] #Cut off the volume prop for cleanliness
colnames(props1_gl)<-c("species1","sizecat1","gl_prey","prop_w1")
props2_gl<-trawlDietSum[,c("comname","sizecat","gl_prey","Wk")]
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




# Original Time -----------------------------------------------------------


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
  mutate(species1=sort(species1,decreasing=T),species2=species2,
         sig=ifelse(s_do>sigGuild_gl,"bold.italic","plain"))%>%
  ggplot()+
  geom_tile(aes(species1,species2,fill=s_do))+
  geom_text(aes(species1,species2,label=round(s_do,digits=2),fontface=sig,color=sig),size=4)+
  scale_fill_viridis_c(option="B",name="Schoener's\nDietary\nOverlap",na.value="white")+
  scale_color_manual(values=c(viridis::viridis(100,option="B")[100],viridis::viridis(100,option="B")[sigGuild_gl*100]),guide=NULL)+
  scale_x_discrete(name="",expand=expansion(0))+
  scale_y_discrete(name="",expand=expansion(0))+
  ggtitle("1973-1997")+
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
png(file = "../Figures/2021 Garrison and Link Rep/fixed_1997dendrogram.png",   # The directory you want to save the file in
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









# Updated Time ------------------------------------------------------------


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
  select(-c("Empty","Miscellaneous","sizecat"))
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

#Visual of the actual matrix with significance indicated
#library(plot.matrix)
#par(mar=c(6,6,5,5.5))
#plot(overlap_mat,axis.row=list(side=2,las=1),col=viridis::viridis(n=100,option="B"),
#     polygon.key = list(border=NA), key=list(),xlab="",ylab="")

#Better visual (ggplot as always)
as.data.frame(overlap_mat_new)%>%
  mutate(species1=colnames(overlap_mat_new))%>%
  pivot_longer(cols=colnames(overlap_mat_new),names_to = "species2", values_to = "s_do")%>%
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
mean(overlap_mat_new[which(overlap_mat_new<1)],na.rm=T)
sd(overlap_mat_new[which(overlap_mat_new<1)],na.rm=T)
range(overlap_mat_new[which(overlap_mat_new<1)],na.rm=T)
hist(overlap_mat_new[which(overlap_mat_new<1)])

#Clustering these out
overlap_clust_new<-hclust(as.dist(1-overlap_mat_new),method="average")
#Trying dendextend to pretty up the dendrogram
dend_new<-as.dendrogram(overlap_clust_new)

#To save them out, there is some difficulty with the export losing the colorbar
png(file = "../Figures/2021 Garrison and Link Rep/fixed_2019dendrogram.png",   # The directory you want to save the file in
    width = 1250, # The width of the plot in inches
    height = 700) # The height of the plot in inches
par(mar=c(5,1,2,12))

dend_new %>% 
  set("branches_lwd", 4) %>%
  # Custom labels
  set("labels_cex", 1) %>%
  #set("labels_col", value = viridis::viridis(12,end=0.8),h = sigGuild_new) %>%
  #set("branches_k_color", value = viridis::viridis(12,end=0.8), h = sigGuild_new) %>%
  plot(horiz=TRUE,main="                               1998-2019 Dendrogram",axes=F,xlab="Similarity")
axis(side=1,at=c(0.6,0.5,0.4,0.3,0.2,0.1,0),labels=c(0.4,0.5,0.6,0.7,0.8,0.9,1))
abline(v=sigGuild_new,lty=2)
#rect.dendrogram(dend_new, h=0.5, which=c(1:4),border="transparent",
#                lty = 5, lwd = 0, col=rgb(0.2,0.2,0.2,0.15),horiz=T)

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




#Comparing the two time periods
t.test(overlap_mat_gl[which(overlap_mat_gl<1)],overlap_mat_new[which(overlap_mat_new<1)])


# Tanglegram --------------------------------------------------------------

#Create a dendlist, which has both the timelines together
dend_compare<-dendlist(dend,dend_new)
par(mar=c(5,1,2,12))
dend_diff(dend,dend_new)




#Doing the more robust comparisons (traditonally used) requires they have identical groups
#Redo the new without Large Little Skate
#Clustering these out
lls<-which(colnames(overlap_mat_new)=="Little Skate_L")
overlap_clust_new2<-hclust(as.dist(1-overlap_mat_new[-lls,-lls]),method="average")
#Trying dendextend to pretty up the dendrogram
dend_new<-as.dendrogram(overlap_clust_new2)

dend_cor<-dendlist(dend,dend_new)%>%
  cor_cophenetic()
dend_tangle<-dendlist(dend,dend_new)%>%
  untangle(method="step2side")%>%
  entanglement()
dendlist(dend,dend_new)%>%
  untangle(method="step2side")%>%
  plot(common_subtrees_color_branches=T,
       main=paste0("1973-1997  |  Entanglement=",round(dend_tangle,2),", Correlation=",round(dend_cor,2),"  |  1998-2019"))



####Ideas####
#Create a tanglegram of the dendrograms
  #Different categories exist in the two (because of L Little Skate)
#How many times is the predator_sizecat in a new cluster, for each predator
#Direction of moves between guilds? i.e., are piscivores becoming benthivores?
  #Mean guild number, if the numbers can be considered ordinal
  #Guild membership size, what's growing and what's decreasing


