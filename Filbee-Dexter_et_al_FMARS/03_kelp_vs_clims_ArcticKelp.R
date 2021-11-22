# analyses/5_kelp_vs_clims.R
# The purpose of this script is to explore relationships between kelp and 
# environmental data
# it outputs kelpforMDS, the main sheet used in PRIMER

library(mgcv)
library(ape)
library(ggplot2)

# Source data ----------------------------------------------------------

#load the environmental data
load("study_site_BO2.RData")

levels(study_site_BO$Region)[levels(study_site_BO$Region)=="Nunatsiavut"] <- "Davis Strait"
levels(study_site_BO$Region)[levels(study_site_BO$Region)=="Ellesmere"] <- "Ellesmere I"


#################################################################################
adf_summary2=read.csv2('kelp environmental variables 10 m.csv')
adf_summary3=read.csv2('kelp environmental variables.csv')
######


#spatial autocorrelation analysis
#generate distance matrix


kelp.dists <- as.matrix(dist(cbind(adf_summary2$lon, adf_summary2$lat)))
kelp.dists.inv <- 1/kelp.dists
kelp.dists.inv[is.infinite(kelp.dists.inv)] <- 0
diag(kelp.dists.inv) <- 0
#kelp.dists.inv[1:5, 1:5] #inspect matrix

#calculate Moran\s I
Moran.I(adf_summary2$kelp.cover, kelp.dists.inv, scaled=T)


#substrate type 
mean(adf_summary3$SandP)
sd(adf_summary3$SandP)/sqrt(nrow(adf_summary3))

kelp_formds <- adf_summary3 %>% 
  dplyr::select(Region,site,depth, Agarum, Alaria, 
                S.latissima, Desmarestia, Fucus, Saccorhiza, 
                L.solidung, L.digitata, Green, Red.algae.fleshy, 
                Filamentous, 'Coralline algae %', BedrockP, BoulderP, 
                CobbleP, PebblesP, SandP, BO_parmean,BO2_phosphatemean_bdmin, 
                BO2_nitratemean_bdmin, BO2_tempmean_bdmin, BO2_salinitymean_bdmin,
                              BO2_icecovermean_ss, BO2_icethickmean_ss) 


#this file  goes into PRIMER. 


################################################################################
#exploring specific relationships highlighted in multivariate analyses 
#temperature, sea ice, substrate and nutrients

range(adf_summary3$BO2_icecovermean_ss)
range(adf_summary3$BO2_icethickrange_ss)
range(adf_summary3$BO2_salinitymean_bdmin)
range(adf_summary3$BO2_phosphatemean_bdmin)
range(adf_summary3$BO2_nitratemean_bdmin)


#order the regions for the plots
adf_summary2$Region=factor(adf_summary2$Region, 
                           levels = c( "Ellesmere I", "North Baffin Island", "Baffin Bay","Davis Strait" ,                
                                       "Foxe Basin"  , "Roes Welcome Sound",  "Hudson Strait" , "Labrador Sea" ,  "Hudson Bay" ))

adf_summary3$Region=factor(adf_summary3$Region, 
                           levels = c( "Ellesmere I", "North Baffin Island", "Baffin Bay","Davis Strait" ,                
                                       "Foxe Basin"  , "Roes Welcome Sound",  "Hudson Strait" , "Labrador Sea" ,  "Hudson Bay" ))

#average across depths
adf_summary_site <- adf_summary3 %>% 
  group_by(site, Region) %>% 
  summarise_if(is.numeric, mean, na.rm = F)


#sea ice and S. latissima cover
m1=gam(S.latissima ~ s(BO2_icethickrange_ss, k = 3,  bs = "cs"), data = adf_summary_site)
m2=gam(S.latissima ~ s(BO2_icethickrange_ss, k = 4,  bs = "cs"), data = adf_summary_site)
m3=gam(S.latissima ~ s(BO2_icethickrange_ss, k = 5,  bs = "cs"), data = adf_summary_site)
anova(m1,m2,m3)

#explore model fit
summary(m2)
res=resid(m2)
plot(adf_summary_site$BO2_icethickrange_ss, res,  ylab="Residuals", xlab="Sea ice") 

fig3r=ggplot(adf_summary_site, aes(x = BO2_icethickrange_ss, y=S.latissima)) +
  geom_smooth(col='blue', size=1, method = "gam", formula = y ~ s(x, bs = "cs", k=4))+
  geom_point(aes(col = Region), size=2)+ 
  xlab('sea ice thickness (m)')+
  scale_color_manual(values=c('black','blue3','cornflowerblue', 'mediumpurple1',
                              'darkcyan','cyan3', 'gold', 'red',
                              'salmon1'), name='')+ ylim(-10,100)+
  scale_size(range = c(1.5, 7), name="Kelp (%)")+ylab(expression(paste(italic("S. latissima"), " (%)", sep=' ')))+
  theme(legend.position = 'right',text = element_text(size = 20))+theme_classic()
fig3r

#only 10 m depth
m1=gam(S.latissima ~ s(BO2_icethickrange_ss, k = 3,  bs = "cs"), data = adf_summary_site)
m2=gam(S.latissima ~ s(BO2_icethickrange_ss, k = 4,  bs = "cs"), data = adf_summary_site)
m3=gam(S.latissima ~ s(BO2_icethickrange_ss, k = 5,  bs = "cs"), data = adf_summary_site)
anova(m1,m2,m3)

#explore model fit
summary(m3)
plot(m3)
res=resid(m3)
plot(adf_summary2$BO2_icethickrange_ss, res,  ylab="Residuals", xlab="sea ice",    main="response") 

fig3r10m=ggplot(adf_summary2, aes(x = BO2_icethickrange_ss, y=S.latissima)) +
  #geom_smooth(col='blue', size=1)
geom_smooth(col='blue', size=1, method = "gam", formula = y ~ s(x, bs = "cs", k=5))+
  geom_point(aes(col = Region), size=2)+ 
xlab('sea ice thickness (m)')+
  scale_color_manual(values=c('black','blue3','cornflowerblue', 'mediumpurple1',
                              'darkcyan','cyan3', 'gold', 'red',
                              'salmon1'), name='')+ ylim(-10,100)+
  scale_size(range = c(1.5, 7), name="Kelp (%)")+ylab(expression(paste(italic("S. latissima"), " (%)", sep=' ')))+
  theme(legend.position = 'right',text = element_text(size = 20))+theme_classic()
fig3r10m


#GAMs for rock and Agarum 
#select k
m1=gam(Agarum ~ s(rock, k = 3,  bs = "cs"), data = adf_summary3)
m2=gam(Agarum ~ s(rock, k = 4,  bs = "cs"), data = adf_summary3)
m3=gam(Agarum ~ s(rock, k = 5,  bs = "cs"), data = adf_summary3)
anova(m1,m2,m3)

#explore model fit
summary(m1)
plot(m1)
res=resid(m1)
plot(adf_summary3$rock, res,  ylab="Residuals", xlab="Rock",    main="response") 

fig3s=ggplot(adf_summary3, aes(x = rock, y=Agarum)) +
  geom_smooth(col='blue', size=1, method = "gam", formula = y ~ s(x, bs = "cs", k=3))+
  geom_point(aes(col = Region), size=2)+ 
  # geom_smooth(col='blue', size=1)+
  xlab('rock (%)')+
  # scale_colour_continuous(name = "Open water\nperiod (%)", low='gold', high='red')+
  scale_color_manual(values=c('black','blue3','cornflowerblue', 'mediumpurple1',
                              'darkcyan','cyan3', 'gold', 'red',
                              'salmon1'),  name='')+
  scale_size(range = c(1.5, 7), name="Kelp %")+ylab(expression(paste(italic("A. clathratum"), " (%)", sep=' ')))+
  theme(legend.position = 'top',text = element_text(size = 20))+theme_classic()
fig3s

#10 m only
#select k
m1=gam(Agarum ~ s(rock, k = 3,  bs = "cs"), data = adf_summary2)
m2=gam(Agarum ~ s(rock, k = 4,  bs = "cs"), data = adf_summary2)
m3=gam(Agarum ~ s(rock, k = 5,  bs = "cs"), data = adf_summary2)
anova(m1,m2,m3)

#explore model fit
summary(m1)
plot(m1)
res=resid(m1)
plot(adf_summary2$rock, res,  ylab="Residuals", xlab="Rock",    main="response") 


fig3s10m=ggplot(adf_summary2, aes(x = rock, y=Agarum)) +
  geom_smooth(col='blue', size=1, method = "gam", formula = y ~ s(x, bs = "cs", k=3))+
  geom_point(aes(col = Region), size=2)+ 
  # geom_smooth(col='blue', size=1)+
  xlab('rock (%)')+
  # scale_colour_continuous(name = "Open water\nperiod (%)", low='gold', high='red')+
  scale_color_manual(values=c('black','blue3','cornflowerblue', 'mediumpurple1',
                              'darkcyan','cyan3', 'gold', 'red',
                              'salmon1'),  name='')+
  scale_size(range = c(1.5, 7), name="Kelp %")+ylab(expression(paste(italic("A. clathratum"), " (%)", sep=' ')))+
  theme(legend.position = 'top',text = element_text(size = 20))+theme_classic()
fig3s10m


fig4=ggarrange(fig3r10m, fig3s10m, nrow=1, common.legend = T, legend = 'right')

#tiff("figures/kelp_cover_vs_sand and ice.v2.tiff", width = 20, height = 11.9, units = "cm", res = 300)
#fig4
#dev.off()

###################################################################

#relationship with other environmental variables
# Visuals -----------------------------------------------------------------
figice=ggplot(adf_summary2, aes(y = BO2_icecovermean_ss, x=lat_BO, 
                                col = Region, size=kelp.cover)) +
  geom_point()+ylab('Sea ice cover (%)')+
  xlab(expression(paste("Latitude (",degree,")")))+theme_classic()+
  guides(col=guide_legend(title=NULL), size = guide_legend("Kelp %"))
 

figtemp=ggplot(adf_summary2, aes(y = BO2_tempmean_bdmin, x=lat_BO, col = Region, size=kelp.cover)) +
  geom_point()+ylab(expression(paste("Sea temperature (",degree,"C)")))+
  xlab(expression(paste("Latitude (",degree,")")))+  
  theme_classic()+  guides(col=guide_legend(title=NULL), size = guide_legend("Kelp %"))

ggarrange( figtemp,figice, legend = 'right', common.legend = T)