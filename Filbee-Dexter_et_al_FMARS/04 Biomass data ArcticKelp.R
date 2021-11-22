
###########################################################################################
library(here)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(plotrix)
library(vegan)
library(randomForest)
library(purrr)
library(carData)
library(car)
library(agricolae)
library(coin)
library(survival)
library(ggridges)
library(ggpubr)
library(RColorBrewer)

####biomass by depth
biomass_depth_long= read.csv2(file= 'kelpbiomassbydepth.csv')

unique(biomass_depth_long$Species)

depthfig=ggplot(biomass_depth_long, aes(y = Biomass_kg_m2, x = factor(Species), fill=factor(Depth_m))) +
  geom_boxplot(notch = F)+
  geom_point(data=depth_summary, aes(y = Mean_biomass_m2, x=Species), col='red',size=3, shape=18, position=position_dodge(width=0.75))+
  geom_errorbar(data=depth_summary,
                aes(x=Species, ymin=Mean_biomass_m2-SE, ymax=Mean_biomass_m2+SE, fill=factor(Depth_m)),
                position = position_dodge(width=0.75), colour='red', width=0,inherit.aes = FALSE)+
  scale_fill_manual(values=c('#7EC0EE','#3A5FCD','#191970'))+
  theme_classic(base_size = 14) +
  labs(y = ("Biomass" ~ (kg ~ m^{-2})),fill = "Depth (m)", x = '')+scale_y_continuous(limits = c(0,12), expand = c(0, 0.3)) +
  theme(legend.position = 'top',  axis.text.x = element_text(angle = 45, vjust=0.5))+ # face = "italic",
  scale_x_discrete(labels=c(expression(italic("A. esculenta")),expression(italic("A. clathratum")),expression(paste(italic("Fucus"),'spp.', sep=' ')),
                            "digitated spp.", expression(italic("L. solidungula")),   expression(italic("S. latissima")), "all kelp"))

depthfig


################################################@###



#setup our environment
#get WW master data
setwd("C:/Users/A22171/Dropbox/ArcticKelp Canada/Arctic kelp extent ms/Analysis/Biomass analysis/Camille's tests")
biomass<-read.csv('Biomass_spreadsheet_V5.csv', dec=',', sep = ';') %>%
  mutate(#Depth_m=as.factor(Depth_m),
    Quadrat_id=as.factor(Quadrat_id),
    Tot_quad_region_10_15=as.numeric(as.character(Tot_quad_region_10_15)),
    Tot_quadrat_depth=as.numeric(as.character(Tot_quadrat_depth)),
    Tot_depth=as.numeric(as.character(Tot_depth)),
    Weight_kg = as.numeric(Weight_g)*0.001)
#View(biomass)




#Summary per site per depth
biomass_depth_reg<-biomass %>%filter(!(Quadrat_id%in%'Transect')) %>% filter(!(Species%in%'Drift')) %>% filter(!(Species%in%'None')) %>%
  filter(!(Species%in%'Filamentous')) %>%  filter(!(Species%in%'Unknown')) %>% filter(!(Species%in%'epiphytes on A. esculenta')) %>%
  filter(!(Species%in%'Red algae')) %>% filter(!(Species%in%'Desmarestia')) %>% filter(!(Species%in%'S. polyschides'))%>% filter(!(Species%in%'Sachoriza'))%>%
  dplyr::group_by(Region,Site, Depth_m,Species)%>% 
  dplyr::summarise(Biomass_kg_m2 = (sum(Weight_kg, na.rm=T)/mean(Tot_quadrat_depth)*4))
#View(biomass_depth_reg)

#Pivot it to wide to build the total kelp column, remove spaces, add zeros instead of NA's
biomass_depth<-biomass_depth_reg %>%
  pivot_wider(names_from = Species, values_from = Biomass_kg_m2)%>%
  rename_all(function(x) gsub(" ", "", x))%>%
  replace(is.na(.), 0) %>%
  mutate(Total.Kelp = Agarum+S.latissima+A.esculenta+L.solidungula+L.digitata+Fucus)
#View(biomass_depth)

#Pivot it back to long to prepare for the figure
biomass_depth_long<-biomass_depth %>% 
  pivot_longer(-c(Region,Site,Depth_m),
  names_to = "Species", values_to = "Biomass_kg_m2")
#View(biomass_depth_long)

#Make summary with mean and std error
depth_summary<-biomass_depth_long%>% 
  dplyr::group_by(Depth_m,Species)%>% 
  dplyr::summarise(Mean_biomass_m2 = mean(Biomass_kg_m2),
            sd_biomass = sd(Biomass_kg_m2),
            count = n())

depth_summary=depth_summary %>%mutate(SE=sd_biomass/sqrt(count))
View(depth_summary)

#stats
m1=lmer(Biomass_kg_m2~Species*Depth_m+(1 | Region/Site), data=biomass_depth_long)
summary(m1)
anova(m1)



#################################################################################################################

#At the thallus level, so mean per species per depth
biomass_reg2<-biomass %>% filter(!(Species%in%'None')) %>% filter(!(Species%in%'Drift')) %>%
  filter(!(Species%in%'Filamentous')) %>%  filter(!(Species%in%'Unknown')) %>% filter(!(Species%in%'epiphytes on A. esculenta')) %>%
  filter(!(Species%in%'Red algae')) %>% filter(!(Species%in%'Desmarestia')) %>% filter(!(Species%in%'S. polyschides')) %>%filter(!(Species%in%'Sachoriza'))%>%
  group_by(Depth_m, Species)%>% summarise_each(funs(mean,std.error),Weight_kg)



#Fig.3 Visualize with line graph based on means and SE, = Plant mean biomass
ggplot(biomass_reg2, aes(x=Depth_m, y=mean, color=Species)) + 
  geom_line(aes(linetype=Species),size=1)+
  geom_point()+
  geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.8,
                position=position_dodge(0.05))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values= c('darkred','red','orange','gold','skyblue1',
                               'skyblue4')) +
  labs(y =expression("Thallus mean biomass" ~ (kg)), fill = "Species",x = "Depth (m)") +
  scale_x_continuous(breaks=seq(0,20,5))+
  theme(legend.position = 'top', legend.title = element_blank())



