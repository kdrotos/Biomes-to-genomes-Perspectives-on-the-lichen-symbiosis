# KDrotos PhD Chapter 4 Analyses and Figures

# author: Katherine Drotos
# last updated 2023-02-27

# Setup -----
rm(list=ls()) # clears out anything that exists in the environment

setwd("C:/Users/Katherine/Documents/Guelph/PhD/THESIS/Chapter 4 - Nuclear content data/Analyses")

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(ggthemes)
library(data.table)
library(ggridges)
library(ggtree)
library(ggtreeExtra)
library(forcats) # for help with re-ordering
library(ggpubr)


library(picante)
library(phytools)
library(vegan)

library(qqplotr)

# __theme_Publication install -----

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# Comparison of KOH treatments (figures included in Chapter 3)----

KOH_dilutions <- read.csv("FIAD_KOH_dilutions.csv")

# need to transpose data to long form

KOH_dilutions_tr <- as.data.frame(t(KOH_dilutions)) #unnecessary with the tidyr gather
KOH_dilutions_tr_L <- gather(KOH_dilutions)

KOH_p1 <- ggplot(KOH_dilutions_tr_L, aes(x=key, y=value)) + geom_violin()
KOH_p1

KOH_p2 <- ggplot(KOH_dilutions_tr_L, aes(x=key, y=value)) +
  geom_violin() + ylab(expression(paste("Genome size (pg)"))) + theme_classic() + coord_flip()

# boxplot looks cleaner

KOH_p3 <- ggplot(KOH_dilutions_tr_L, aes(x=key, y=value)) +
  geom_boxplot() + ylab(expression(paste("Genome size (pg)"))) + theme_classic() + coord_flip()

#__Graphing within species KOH treatment comparisons ----

#_____Xanthomendoza fallax ----

Xfallax <- subset(KOH_dilutions_tr_L, key == 'Xfallax_water'| key == 'Xfallax_0.1KOH' | key == 'Xfallax_1KOH' | key == 'Xfallax_10KOH')

Xfallax$key <- factor(Xfallax$key, levels=c('Xfallax_water', 'Xfallax_0.1KOH', 'Xfallax_1KOH', 'Xfallax_10KOH'))

Xf_p1 <- ggplot(Xfallax, aes(x=key, y=value)) +
  geom_boxplot() + ylab(expression(paste("Genome size (pg)"))) + theme_classic() + coord_flip()

Xf_p2 <- ggplot(Xfallax, aes(x=key, y=value, fill=key)) +
  geom_boxplot() + ylab(expression(paste("Genome size (pg)")))+ theme_classic() + coord_flip() 

Xf_xticks <- c("Water", "0.1% KOH", "1% KOH", "10% KOH") # create vector for renaming x ticks

Xf_p3 <- ggplot(Xfallax, aes(x=key, y=value, fill=key)) +
  geom_boxplot() + scale_fill_brewer(palette = "Blues") + 
  scale_x_discrete(labels = Xf_xticks) + 
  ylab(expression(paste("Genome size (pg)")))+ 
  theme_classic() + 
  xlab(expression(paste("Preparation treatment")))+ 
  theme_Publication() + 
  theme(legend.position="none")

#_____Physcia stellaris ----

Pstellaris <- subset(KOH_dilutions_tr_L, key == 'Physcia_stellaris_water' | key == 'Physcia_stellaris_0.1KOH' | key == 'Physcia_stellaris_1KOH' | key == 'Physcia_stellaris_10KOH')

Ps_xticks <- c("Water", "0.1% KOH", "1% KOH", "10% KOH")

Pstellaris$key <- factor(Pstellaris$key, levels=c('Physcia_stellaris_water', 'Physcia_stellaris_0.1KOH', 'Physcia_stellaris_1KOH', 'Physcia_stellaris_10KOH'))

Ps_p1 <- ggplot(Pstellaris, aes(x=key, y=value, fill=key)) +
  geom_boxplot() + 
  scale_fill_brewer(palette = "Reds") + 
  scale_x_discrete(labels = Xf_xticks) + 
  ylab(expression(paste("Genome size (pg)")))+ 
  theme_classic() + 
  xlab(expression(paste("Preparation treatment"))) + 
  theme_Publication() + 
  theme(legend.position="none")


# Genome size data by species -----

LFF_gs <- read.csv("FIAD_data_by_species.csv")

LFF_gs_tr <- gather(LFF_gs)

LFF_gs_p1 <- ggplot(LFF_gs_tr, aes(x=key, y=value)) + geom_boxplot() + coord_flip()

# rainbow colouring

LFF_gs_p2 <- ggplot(LFF_gs_tr, aes(x=key, y=value, fill=key)) +
  geom_boxplot() + 
  scale_fill_discrete() + 
  ylab(expression(paste("Genome size (pg)")))+ 
  theme_classic() + theme_Publication() + 
  theme(legend.position="none") + 
  theme(axis.title.y=element_blank()) + 
  coord_flip()

# grey with jitters

LFF_gs_p3 <- ggplot(LFF_gs_tr, aes(x=key, y=value, fill=key)) +
  geom_boxplot(fill = "grey") + 
  geom_jitter(width = 0.3, size = 0.6, colour=4) +
  ylab(expression(paste("Genome size (pg)")))+ 
  theme_classic() + theme_Publication() + 
  theme(legend.position="none") + 
  theme(axis.title.y=element_blank()) + 
  coord_flip()  


# General stats ----
#_____Mean of each species ----

setDT(LFF_gs_tr)

LFF_gs_tr[,list(mean=mean(value, na.rm=TRUE)), by=key]

write.table(LFF_gs_tr[,list(mean=mean(value, na.rm=TRUE)), by=key])

#_____Count of each species -----

countLFF_gs <- LFF_gs_tr %>% group_by(key) %>% summarise(total_non_na = sum(!is.na(value)))

print.data.frame(countLFF_gs)

#_____Standard deviation of each species -----

sd_LFF_gs <- LFF_gs_tr[ ,list(sd=sd(value, na.rm=TRUE)), by = key]

sd_LFF_gs[order(sd_LFF_gs$key),]

# Comparison of lineages ----
#__Bringing in trees

testtree <- ape::read.tree(file="Lecanoromycetes_Nelsen2020_edited.txt")
plot(testtree)

plot(testtree, type="fan", show.tip.label = FALSE)

#okay that works, but it is very messy. I think I'm going to build from the backbone down
# will test each stage to make sure it works

Fphyla_tree <- ape::read.tree(text='(Other_fungi,(Basidiomycota,Ascomycota));')

# Ascomycota big tree -----
# makin' trees to graft from the deeper branches forward

Asco_tree <- ape::read.tree(text='(Taphrinomycotina,(Saccharomycotina,Pezizomycotina));')

Pezizo <- ape::read.tree(text= '((Pezizomycetes,(((Arthoniomycetes, Dothidiomycetes),(Xylonomycetes,(Eurotiomycetes,Lecanoromycetes))),(Leotiomycetes,Sordariomycetes))));')

p1 <- ape::read.tree(text='(Other_fungi,(Basidiomycota,(Taphrinomycotina,(Saccharomycotina,Pezizomycotina))));')

Sordariomycetes <- ape::read.tree(text= '(Xylariales,(((Microascales,(Hypocreales,Glomerellales)),(((((Sordariales,Boliniales),Chaetosphaerales), Coniochaetales)),((Togniniales,Diaporthales),(Magnaporthales,Ophiostomatales))))));')

Leotiomycetes <- ape::read.tree(text='(Rhytismatales,(Helotiales,Erysiphales));')

Dothideomycetes <- ape::read.tree(text= '(((Hysteriales,Phaeotrichales),(Pleosporales,(Botryosphaeriales, Trypetheliales))),((Capnodiales,(Dothideales,Myriangales))));')

Eurotiomycetes <- ape::read.tree(text='(Mycocalicilaes,(Chaetothyriales,(Ascosphaerales,(Eurotiales,Onygenales))))')

Lecanoromycetes <- ape::read.tree(text= '((Umbilicariales,(Pertusariales, Ostropales)),(Lecideales,(Caliciales,(Peltigerales,(Teloschistales,Lecanorales)))));')

Taphr1 <- ape::read.tree(text= '(Neolectales,((Taphrinales,(Pneumocystidales,Schizosaccharomycetales))));')

# Big, total, working tree! ----

FullAscoTree <- ape::read.tree(text= '(Other_fungi,(Basidiomycota,((Taphrinales,(Pneumocystidales,Schizosaccharomycetales)),(Saccharomycetales,((Pezizales,((Orbiliales,((((((Candelariales,((Mytilinidiales,((Hysteriales,Phaeotrichales),(Pleosporales,(Botryosphaeriales, Trypetheliales)))),((Capnodiales,Dothideales)))))),((Xylonales,Symbiotaphrinales),((Mycocaliciales,(Chaetothyriales,(Ascosphaerales,(Eurotiales,Onygenales)))),((Umbilicariales,(Pertusariales, Ostropales)),(Lecideales,(Caliciales,(Peltigerales,(Teloschistales,Lecanorales))))))))),(((Rhytismatales,(Helotiales,Erysiphales))),(Xylariales,(((Microascales,(Hypocreales,Glomerellales)),(((((Sordariales,Boliniales),Chaetosphaeriales), Coniochaetales)),(Diaporthales,(Magnaporthales,Ophiostomatales))))))))))))))));')

# plotted nicely

FA_Tree <- ggtree(FullAscoTree, mapping=NULL) + geom_tiplab(mapping=NULL)

# Lichen tree only ----
# do each family that we have more than 2 data points for

Cladoniaceae <- ape::read.tree(text='(Cladonia_grayi,(Cladonia_cristatella, Cladonia_rangiferina));')

Lobariaceae <- ape::read.tree(text='(Lobaria_pulmonaria, Lobaria_oregana);')

Parmeliaceae <- ape::read.tree(text='(Letharia_columbiana,((Pseudevernia_consocians,(Cetrariella_delisei,(Evernia_mesomorpha, Evernia_perfragilis))),(Platismatia_tuckermanii,(Parmelia_sulcata,((Xanthoparmelia_cumberlandia, Xanthoparmelia_conspersa),(Punctelia_rudecta,(Parmotrema_crinitum,(Flavoparmelia_caperata, Flavoparmelia_baltimorensis))))))));')

Physciaceae <- ape::read.tree(text='(Physcia_stellaris,(Physconia_muscigena, Anaptychia_crinalis));')
  
Teloschistaceae <- ape::read.tree(text='(Xanthomendoza_fallax,(Xanthoria_elegans, Xanthoria_parietina));')

# plotting just the groups from my data

Lecanorales1 <- ape::read.tree(text='(Pilocarpaceae,(Parmeliaceae, Cladoniaceae));')

Lichen_orders1 <- ape::read.tree(text='(Trypetheliales,(Candelariales,(Umbilicariales,(Pertusariales, Ostropales),(Lecideales,(Caliciales,(Peltigerales,(Teloschistales,Lecanorales)))))));')

# plugging everything in to the lichen-only tree

Lichen_tree <- ape::read.tree(text='((Candelaria_concolor,(Trypethelium_eluteriae,(Lasallia_pensylvanica,((Pertusaria_dactylina, Stictis_radiata),(Porpidia_crustulata,((Calicium_adspersum,(Physcia_stellaris,(Physconia_muscigena, Anaptychia_crinalis))),((Peltigera_elizabethae, (Lobaria_pulmonaria, Lobaria_oregana)),((Xanthomendoza_fallax,(Xanthoria_elegans, Xanthoria_parietina)),(Micarea_prasina,((Letharia_columbiana,((Pseudevernia_consocians,(Cetrariella_delisei,(Evernia_mesomorpha, Evernia_perfragilis))),(Platismatia_tuckermanii,(Parmelia_sulcata,((Xanthoparmelia_cumberlandia, Xanthoparmelia_conspersa),(Punctelia_rudecta,(Parmotrema_crinitum,(Flavoparmelia_caperata, Flavoparmelia_baltimorensis)))))))), (Cladonia_grayi,(Cladonia_cristatella, Cladonia_rangiferina)))))))))))));')

# species and genome size data -----

Full_species_info_data_for_R <- read_excel("~/Guelph/PhD/THESIS/Chapter 4 - Nuclear content data/Full_species_info_data_MASTER.xlsx",
col_types = c("text", "text", "text",
"text", "text", "text", "text", "numeric",
"text", "text", "text"))

# need to subset lichens from full genome size dataset for this

LichenDataOnly <-Full_species_info_data_for_R[Full_species_info_data_for_R$Symbiotic_state == 'lichen',]
subset(Full_species_info_data_for_R, Symbiotic_state == 'lichen')

# average by species

LichenDataOnly %>% group_by(Species, Genus) %>% summarize(Genome_size= mean(Genome_size)) -> LGS_range

# group genus and species together to make binomial names

LGS_range$Full_names <- paste(LGS_range$Genus, LGS_range$Species, sep="_")

# putting the lichen tree and data together ----

#quick test to make sure the data is reading properly
ggplot(LGS_range, aes(x=Genome_size, y=Full_names)) + geom_point()

# I think this is necessary? Was throwing errors before
LGS_range$Genome_size <- as.numeric(LGS_range$Genome_size)
LGS_range$Full_names <- as.vector(LGS_range$Full_names)

# have to subset the data again, it doesn't like the extra columns

LGS_reduced <- subset(LGS_range, select= -c(Genus, Species))

#___matching the tree and data together

lichen_tree_gs <- LGS_reduced %>% select(Full_names, Genome_size) %>% deframe() %>% match.phylo.data(Lichen_tree, .) # this made it as a list, need to convert to a dataframe

LTgs.df <- as.data.frame(lichen_tree_gs$data) %>% rownames_to_column("Species") %>% rename('Genome_size' = "lichen_tree_gs$data")

# now plotting

LichenOnlyTree <- ggtree(lichen_tree_gs$phy,
              branch.length = 0.5) +
  geom_tiplab(mapping=NULL) +
  geom_fruit(data=LTgs.df,
             geom=geom_bar,
             mapping=aes(x=Genome_size, y=Species),
             orientation="y",
             stat="identity",
             offset = 0.4,
             axis.params = list(axis="x", vjust=1.5, text.size=2)) +
  xlab(expression(paste("Genome size (pg)"))) +
  theme(axis.title.x = element_text(hjust=0.96))

ggsave(LichenOnlyTree, filename="LichenOnlyTree2.png", width=10, height=5.5, units=c("in"))

#___Comparing structural forms ----

LichenDataOnly  %>% group_by(Full_name) %>% mutate(Avg_gs= mean(Genome_size)) -> LichenData_Ready

LF_Structure <- LichenData_Ready %>% 
  filter(!(LF_structure =="NA")) %>% 
ggplot(data=., aes(x=LF_structure, y=Avg_gs, fill=LF_structure)) +
  geom_boxplot() +
  xlab(expression(paste("Lichen structural form"))) +
  ylab(expression(paste("Genome size (pg)")))+
  theme_Publication() +
  scale_fill_brewer(palette="Accent") +
  theme(legend.position="none")

ggsave(LF_Structure, filename="LF_structure2.png")

# is it significant? violates assumptions, but ...

aov(Avg_gs ~ LF_structure, data=LichenData_Ready)

# because the data are non normal, and is from 0 to +infinity, we need to choose a different error distribution to run the ANOVA on
# family = Gamma(link =log) looks to be best for positive continuous data

# normal ANOVA assumes gaussian - we need to change it to glm

LichenGLM <- glm(family=Gamma(link= 'inverse'), formula= Avg_gs ~ LF_structure, data=LichenData_Ready)

summary(LichenGLM) # missing crustose because that's what it's comparing against

# estimate = parameter in the equation; t-value = test statistic to see if comparison is different than zero

summary(aov(LichenGLM))
# output:
#             Df   Sum Sq   Mean Sq F value Pr(>F)
#LF_structure  3 0.000259 8.641e-05   0.666  0.579
#Residuals    31 0.004020 1.297e-04               
#2 observations deleted due to missingness
# So, this indicates no effect of LF structure - likely due to sample size and variance in this case

plot(LichenGLM)

# residuals vs fitted plot - shows "cone" that indicates a lack of constant variance (which is one of the assumptions, so yea we violated it)

# Joining the big tree and with data ----

# average genome size per species into a new column
# using mutate function instead of summarize so I don't lose the other columns

Full_species_info_data_for_R  %>% group_by(Full_name) %>% mutate(Avg_gs= mean(Genome_size)) -> AllFungiData

# base plot

ggplot(AllFungiData, aes (x=Avg_gs, y=Order)) + geom_boxplot()

# need a vector to order the data so it matches the tree

AllFungiData$Order <- factor(AllFungiData$Order, levels = c("Other_fungi","Basidiomycota","Taphrinales","Pneumocystidales","Schizosaccharomycetales","Saccharomycetales","Pezizales","Orbiliales","Rhytismatales","Helotiales","Erysiphales", "Xylariales","Microascales","Hypocreales","Glomerellales","Diaporthales","Magnaporthales","Ophiostomatales","Coniochaetales","Chaetosphaeriales","Sordariales","Boliniales","Candelariales","Capnodiales","Dothideales","Mytilinidiales","Hysteriales","Phaeotrichales","Pleosporales","Botryosphaeriales","Trypetheliales",   "Xylonales","Symbiotaphrinales","Mycocaliciales","Chaetothyriales","Ascosphaerales","Eurotiales","Onygenales","Umbilicariales","Pertusariales","Ostropales","Lecideales", "Caliciales","Peltigerales","Teloschistales","Lecanorales")) 
                             
# backup copy of old order; it changed due to fixing a tree error
# levels = c("Other_fungi","Basidiomycota","Taphrinales","Pneumocystidales","Schizosaccharomycetales","Saccharomycetales","Pezizales","Orbiliales","Rhytismatales","Helotiales","Erysiphales", "Xylariales","Microascales","Hypocreales","Glomerellales","Diaporthales","Magnaporthales","Ophiostomatales","Coniochaetales","Chaetosphaeriales","Sordariales","Boliniales","Candelariales","Symbiotaphrinales","Capnodiales","Dothideales","Mytilinidiales","Hysteriales","Phaeotrichales","Pleosporales","Botryosphaeriales","Trypetheliales",   "Xylonales","Mycocaliciales","Chaetothyriales","Ascosphaerales","Eurotiales","Onygenales","Umbilicariales","Pertusariales","Ostropales","Lecideales", "Caliciales","Peltigerales","Teloschistales","Lecanorales")) 


#___gs across orders -----

ggplot(AllFungiData, aes (x=log(Avg_gs), y=Order, fill=Class)) + 
  geom_boxplot() + xlab(expression(paste("log Genome size (pg)"))) + 
  theme_Publication() + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position = "none", axis.text.y=element_text(size=8))

#___symbiotic state counts across orders ----

# don't need outgroups for this, since they have no data and their NAs show up weird
AllFungiData2 <- AllFungiData[-c(1365,1366),]

# plot - colour by Symbiotic state
  
ggplot(AllFungiData2, aes(x=Symbiotic_state, y=Order, colour=Symbiotic_state)) +
  geom_count() +xlab(expression(paste("Symbiotic state")))+
  theme_Publication() + 
  theme(axis.text.x=element_text(angle=90, size=8), axis.text.y=element_text(size=8), axis.title.y=element_blank(),axis.title.x=element_text(size=10), legend.position = "right",legend.direction="vertical")+ 
  guides(colour= FALSE) + 
  scale_colour_colorblind()

# same plot, but with colour by Class

ggplot(AllFungiData2, aes(x=Symbiotic_state, y=Order, colour=Class)) +
  geom_count() +xlab(expression(paste("Symbiotic state")))+
  theme_Publication() + 
  theme(axis.text.x=element_text(angle=90, size=8), axis.text.y=element_text(size=8), axis.title.y=element_blank(),axis.title.x=element_text(size=10), legend.position = "right",legend.direction="vertical")+ 
  guides(colour= FALSE) + 
  scale_colour_colorblind()

# need to export these with transparent backgrounds

#colourblind palette for more than 8 categories

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#333888", "#AA4499", 
                                      "#44AA99", "#999933","#0069FF", "#882255", "#661100", "#6699CC", "#888888")
                                      scales::show_col(safe_colorblind_palette)

#___exporting final plots -----                                      
#tree 
AscoTree_nonames <-ggtree(FullAscoTree, mapping=NULL) + 
                                        theme(plot.background = element_rect(fill="transparent",color=NA), panel.background = element_rect(fill="transparent"))

ggsave(AscoTree_nonames, filename ="FullAScoTree_transp2.png", bg="transparent")

#boxplot
AllF_boxplot <- ggplot(AllFungiData, aes (x=log(Avg_gs), y=Order, fill=Class)) + 
  geom_boxplot() + 
  xlab(expression(paste("log Genome size (pg)"))) + 
  theme_Publication() + 
  theme(axis.title.y = element_blank()) + 
  theme(legend.position = "none", axis.text.y=element_text(size=8)) + 
  scale_fill_manual(values=safe_colorblind_palette) + 
  theme(plot.background = element_rect(fill="transparent",color=NA), panel.background = element_rect(fill="transparent"))

ggsave(AllF_boxplot, filename ="AllFgs_transp1.png", bg="transparent")

# count of symbiotic states
AllF_count <- ggplot(AllFungiData2, aes(x=Symbiotic_state, y=Order, colour=Class)) +
  geom_count() +
  xlab(expression(paste("Symbiotic state")))+
  theme_Publication() + 
  theme(axis.text.x=element_text(angle=90, size=8), axis.text.y=element_text(size=8), axis.title.y=element_blank(),axis.title.x=element_text(size=10), legend.position = "right",legend.direction="vertical")+ guides(colour= FALSE) + 
  scale_colour_manual(values=safe_colorblind_palette) + 
  theme(plot.background = element_rect(fill="transparent",color=NA), panel.background = element_rect(fill="transparent"))

ggsave(AllF_count, filename ="AllFcount_transp1.png", bg="transparent")

# genome size by symbiotic state

gs_symstate <- AllFungiData %>% 
  filter(!(Class %in% c("Basidiomycota","Other"))) %>% 
ggplot(data=., aes(x=log(Avg_gs), y=Symbiotic_state, fill=Symbiotic_state)) +
  geom_boxplot() +
  xlab(expression(paste("log Genome size (pg)"))) +
  ylab(expression(paste("Symbiotic state"))) +
  theme_Publication() +
  theme(legend.position = "none") +
  scale_fill_colorblind()

ggsave(gs_symstate, filename="gs_symstate1.png", width = 6, height=4, units=c("in"))

# if we want no NAs
AllFungiData %>% 
  filter(!(Symbiotic_state =="NA")) %>% 
  ggplot(data=., aes(x=log(Avg_gs), y=Symbiotic_state)) +
  geom_boxplot()


# Comparing lichenized and non-lichenized lineages -----
# need to prepare data frame first



# first, check distribution of data to see if normal - Shapiro-Wilk test

ggqqplot(AllFungiData$Genome_size)

AllFungiData %>% 
  filter(!(Order =="Pezizales")) %>% 
  ggplot(data=., mapping=aes(sample=log(Avg_gs), fill=Order)) +
  stat_qq_band(alpha=0.15) +
  stat_qq_line(alpha=0.15) +
  stat_qq_point(size=0.5) + 
  facet_wrap(~ Order)

# unpaired t-tests
# note: here both assumptions are violated (not normally distributed, variances are likely not equal)

#___Lecanoromycetes vs Eurotiomycetes -----

# F test
AllFungiData %>% 
  filter(Class %in% c("Lecanoromycetes", "Eurotiomycetes")) %>% 
  var.test(Avg_gs ~ Class, data=.)

# results: F test to compare two variances
#data:  Avg_gs by Class
#F = 0.18646, num df = 214, denom df = 32, p-value = 3.413e-14
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
#  0.1038174 0.3014222
#  ratio of variances 
#0.1864602 

# unpaired two-tailed Welch's T-test

LFF_comparison <- AllFungiData %>% 
  filter(Class %in% c("Lecanoromycetes", "Eurotiomycetes")) %>% 
  t.test(Avg_gs ~ Class, data=., var.equal = (FALSE))

# Eurotiomycetes mean genome size is significantly larger Lecanoromycetes

#___Trypetheliales vs Botryosphaeriales -----

# F test
AllFungiData %>% 
  filter(Order %in% c("Trypetheliales", "Botryosphaeriales")) %>% 
  var.test(Avg_gs ~ Order, data=.)
# not enough data - Trypetheliales has only 1 data point

# no test for me :(

# Candelariales vs Dothideomycetes

# only 1 data point for Candelariales - no tests :(

#___boxplots to visually compare all three pairs ----

# Comp 1

AllFungiData$Class <- factor(AllFungiData$Class, levels = c("Eurotiomycetes","Lecanoromycetes","Dothideomycetes","Candelariomycetes"))

class.colors <- c("Lecanoromycetes" = "#88CCEE", "Eurotiomycetes" = "#CC6677", "Candelariomycetes" = "#88CCEE", "Dothideomycetes"="#CC6677")

Comp1 <- AllFungiData %>% 
  filter(Class %in% c("Lecanoromycetes", "Eurotiomycetes","Candelariomycetes", "Dothideomycetes")) %>% 
  ggplot(data=., aes(x=Class, y=log(Avg_gs), fill=Class)) + 
  geom_boxplot() + 
  labs(x="Lineage", y="log Genome size (pg)") +
  scale_fill_manual(values=class.colors)+
  theme_Publication() +
  theme(legend.position="none", axis.title.x=element_blank(),axis.text.x=element_text(size=8))

ggsave(Comp1, filename ="Comp1.2.png", width=8, height=6, units=c("in"))

# Comp 2

order.colours <- c("Trypetheliales"="#88CCEE", "Botryosphaeriales"="#CC6677")

Comp2 <- AllFungiData %>% 
  filter(Order %in% c("Trypetheliales", "Botryosphaeriales")) %>% 
  ggplot(data=., aes(x=Order, y=log(Avg_gs), fill=Order)) + 
  geom_boxplot() +
  labs(x="Lineage", y="log Genome size (pg)") +
  scale_fill_manual(values=order.colours) +
  theme_Publication() +
  theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_text(size=8)) +
  coord_cartesian(ylim=c(0,0.22)) 

ggsave(Comp2, filename ="Comp2.2.png", width=4.5, height=6, units=c("in"))






