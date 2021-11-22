#01.11.2021. ArcticKelp analysis
#updated the R file, so no more issues with the library. 


#.libPaths(c("~/R-packages", .libPaths()))
#.libPaths('C:/Users/A22171/Documents/R/R-3.5.1/library')

library(tidyverse)
library(ggridges)
library(ggpubr)
library(RColorBrewer)
library(writexl)
library(tidyr)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(plotrix)


setwd("C:/Users/A22171/Dropbox/ArcticKelp Canada/Arctic kelp extent ms/Analysis/ArcticKelpCanada/data/open source data")

adf_summary<-adf_summary%>%filter(!(family=="kelp.cover"))
colnames(study_sites)


fig3<-adf_summary %>% 
  ggplot(aes(x = reorder(Site, desc(Site)), y = mean_cover, fill = family)) +
  facet_grid(. ~ factor(fDepth))+
  geom_col(position=position_stack(reverse = TRUE),  col='black') + theme_bw() +
  labs(x = '', y = "Kelp cover (%)", fill = "family") + ylim(0,100)+
  scale_fill_manual(name = "Kelp", breaks = c("Agarum", "Alaria","S.latissima","L.solidung","L.digitata","Saccorhiza"),values= c('darkred','red','orange','gold','skyblue4','skyblue1'), 
                  labels= c("A. clathratum", "A. esculenta","S. latissima", "L. solidungula","H. nigripes/L. digitata","S. dermatodes")) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),legend.title=element_blank(), strip.background = element_blank(), 
        legend.position = 'top',text = element_text(size = 14),legend.text = element_text(size = 10, face="italic")) +
  coord_flip() 

fig3 



#tiff("figures/kelp_cover_vs_sites_V2.tiff", width = 25, height = 23, units = 'cm', res = 300)
#fig3 # Make plot
#dev.off()

######################################################################################################
#red algae

#get summary data on red algae cover
adf_summary_red<-read.csv2('other seaweed cover summary by site.csv')


#averages of other species.
#mean statistics palmeria
desmarestia<-adf_summary_red%>% filter(family=='Desmarestia')
mean(desmarestia$mean_cover); range(desmarestia$mean_cover)#mean
sd(desmarestia$mean_cover)/sqrt(nrow(desmarestia)) #SE

fucus<-adf_summary_red%>% filter(family=='Fucus')
mean(fucus$mean_cover); range(fucus$mean_cover)#mean
sd(fucus$mean_cover)/sqrt(nrow(desmarestia)) #SE

#mean statistics palmeria
palmeria<-adf_summary_red%>% filter(family=='Red.algae.fleshy')
mean(palmeria$mean_cover); range(palmeria$mean_cover)#mean
sd(palmeria$mean_cover)/sqrt(nrow(palmeria)) #SE

#mean statistics filamentous reds
filamentous<-adf_summary_red%>% filter(family=='Filamentous')
mean(filamentous$mean_cover); range(filamentous$mean_cover)#mean
sd(filamentous$mean_cover)/sqrt(nrow(filamentous)) #SE

reds<-adf_summary_red%>% filter(family=='Alga.Roja.fol')
mean(reds$mean_cover); range(reds$mean_cover)#mean
sd(reds$mean_cover)/sqrt(nrow(reds)) #SE

greens<-adf_summary_red%>% filter(family=='Green')
mean(greens$mean_cover); range(greens$mean_cover)#mean
sd(greens$mean_cover)/sqrt(nrow(greens)) #SE

#plot other seaweeds

figred<-adf_summary_red %>% 
  ggplot(aes(x = reorder(Site, desc(Site)), y = mean_cover, fill = family)) +
  facet_grid(. ~ factor(fDepth))+
  geom_col(position=position_stack(reverse = TRUE),  col='black') +
  theme_bw() +
  # geom_errorbar(aes(ymin=mean_cover,ymax=mean_cover+se), col='black')+
  labs(x = '', y = "other macroalgal cover (%)", fill = "family") +
  ylim(0,100)+
  #scale_fill_brewer(palette = 'Reds')+
  scale_fill_manual(name = "Red algae", breaks = c('Fucus', 'Saccorhiza','Desmarestia', 'Red.algae.fleshy', "Filamentous","Alga.Roja.fol","Green"),
                    values= c(brewer.pal(3,"YlOrRd")[1:2], brewer.pal(4,"Reds"), 'skyblue2'),#"#FEE5D9", "#FC9272", "#EF3B2C",'#99000D','coral4'
                   labels= c(expression(paste(italic("Fucus "),'spp.', sep=' ')), expression(italic("S. dermatodea")),expression(paste(italic("Desmarestia "),'spp.', sep=' ')), expression(italic('P. palmata')),"red/brown filamentous",'reds',"greens" )) +
  
 # scale_fill_manual(name = "Kelp", breaks = c("Agarum", "Alaria","S.latissima","L.solidung","L.digitata","Saccorhiza"),values= c('darkred','red','orange','gold','skyblue4','skyblue1'), 
#                    labels= c("A. clathratum", "A. esculenta","S. latissima", "L. solidungula","H. nigripes/L. digitata","S. dermatodes")) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),legend.title=element_blank(), strip.background = element_blank(), 
        legend.position = 'top',text = element_text(size = 14),legend.text = element_text(size = 10), legend.text.align = 0) +
  coord_flip() 

figred 

#dev.off()
#tiff("figures/red_cover_vs_sites.tiff", width = 25, height = 23, units = 'cm', res = 300)
#figred # Make plot
#dev.off()

