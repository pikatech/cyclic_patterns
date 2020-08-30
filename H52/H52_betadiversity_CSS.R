rm(list = ls())
setwd("/home/yo39qed/time-series analysis/output")

#saveRDS()
library(phyloseq)
library(reshape2)
library(readr)
library(RADanalysis)
require(gridExtra)
library(ggalluvial)
library(gganimate)
library(gifski)
library(devtools)
library(pairwiseAdonis)
library(vegan)
library(ggplot2)
library(ggrepel)
library(plotrix)

set_colors_Year <- get_palette("lancet", 7)


set_colors_recharge <- c("Recharge"="orange", "Discharge"="grey")
set_shapes_recharge <- c("Recharge"=15, "Discharge"=1)

set_colors_Period <-c("P13"= "brown","P15"="red", "P16"="deeppink", "P17"="orange",
                      "P18"="purple", "P19"="blue", "P14"="darkgreen")

set_colors_PeriodRecharge <- c("P13Recharge" ="brown","P15Recharge" = "red", "P15Discharge" ="red", 
                               "P16Recharge" ="deeppink", "P16Discharge"="deeppink", "P17Recharge"="orange",    
                               "P17Discharge"="orange", "P18Recharge" ="purple",  "P18Discharge"= "purple",
                               "P19Recharge"= "blue", 
                               "P19Discharge"= "blue", "P14Recharge"="darkgreen", "P14Discharge"="darkgreen")


color2 <- c("Alphaproteobacteria" = "purple4", "Gammaproteobacteria"="mediumpurple2", "Deltaproteobacteria" = "mediumorchid2", "Proteobacteria_unclassified"="magenta1",
            "Actinobacteria"="brown",
            "Bacteria_unclassified"="gray90",
            "Bacteroidia"="orange1", "Ignavibacteria" = "yellow",
            "Dehalococcoidia"="lightpink",
            "Elusimicrobia" = "burlywood",
            "Nitrospira"="blue", "Thermodesulfovibrionia"="royalblue", 
            "Saccharimonadia" = "green4", "Parcubacteria"="green2", "Gracilibacteria"="darkolivegreen4", "ABY1"="lightgreen",
            "Omnitrophicaeota_cl"="darkorange1", 
            "Brocadiae"="firebrick2",
            "Others"="gray",
            "Chlamydiae" = "black"
) 

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


# Load data
load("/home/yo39qed/time-series analysis/output/H52/data_norm52.phyloseq")
#################################################################################

#################################################################################
###           PERMANOVA using adonis function in Vegan package              #####
#################################################################################
# recharge and Period
## remove no-water-level data
# data_norm52_select <- subset_samples(data_norm52, Period !="Period0")
data_norm52_select <- data_norm52
data_norm52_select <- filter_taxa(data_norm52_select , function(x) sum(x) > 0, TRUE)

set.seed(123)
data_norm52_bray <- phyloseq::distance(data_norm52_select, method = "bray") 
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(data_norm52_select))

# Adonis test (Number of permutations: 999)
set.seed(123)
ht_RechargePeriod <- adonis(data_norm52_bray ~ Recharge + Period + Recharge*Period, data = sampledf) #  ***
# This output tells us that our adonis test is significant 
# so we can reject the null hypothesis that our three countries have the same centroid.
ht_RechargePeriod
ht_RechargePeriod$aov.tab$"Pr(>F)" # use the R2 values for VennDiagram
write.table(data.frame(ht_RechargePeriod$aov.tab), "H52/bray/permanova_ht_RechargePeriod_H52.csv", sep=";", dec=",", col.names=NA)

### Pairwise comparison ##
set.seed(123)
adonis_pair_Period <- pairwise.adonis(vegan_otu(data_norm52_select), sampledf$Period)
adonis_pair_Period
write.table(data.frame(adonis_pair_Period), "H52/bray/permanova_Period_pairwise_H52.csv", sep=";", dec=",", col.names=NA)

####
set.seed(123)
adonis_pair_PeriodRecharge <- pairwise.adonis(vegan_otu(data_norm52_select), sampledf$PeriodRecharge)
adonis_pair_PeriodRecharge
write.table(data.frame(adonis_pair_PeriodRecharge), "H52/bray/permanova_PeriodRecharge_pairwise_H52.csv", sep=";", dec=",", col.names=NA)

# Homogeneity of dispersion test (Number of permutations: 999)
library(vegan)
dis_Period<- betadisper(data_norm52_bray, sampledf$Period)
set.seed(123)
permutest(dis_Period) # 0.003**

# Homogeneity of dispersion test (Number of permutations: 999)
library(vegan)
dis_PeriodRecharge<- betadisper(data_norm52_bray, sampledf$PeriodRecharge)
set.seed(123)
permutest(dis_PeriodRecharge) # 0.003**


##################################################
# Unconstrained ordination: PCoA using bray      #
##################################################

set.seed(123)
ordu = ordinate(data_norm52, "PCoA", "bray")
scree <- plot_scree(ordu, "Scree plot for bacterial community, 0.2 filter, bray")
scree$data
tbd <- cbind(ordu$vectors[, 1:2], sample_data(data_norm52))

#### define error bars: not scientific, because bray-curtis is non-metric #####
ag1 <- aggregate(Axis.1 ~ Period + Recharge, data = tbd, 
          FUN = function(x) c(mean = mean(x), se = std.error(x)))


ag2 <- aggregate(Axis.2 ~ Period + Recharge, data = tbd, 
                 FUN = function(x) c(mean = mean(x), se = std.error(x)))


ag <- data.frame(merge(ag1, ag2))
ag$Axis.1.mean <- ag$Axis.1[,1]
ag$Axis.1.se <- ag$Axis.1[,2]
ag$Axis.2.mean <- ag$Axis.2[,1]
ag$Axis.2.se <- ag$Axis.2[,2]

tp_final <- ggscatter(ag, x = "Axis.1.mean", y = "Axis.2.mean",
                 color = "Period", 
                 shape = "Recharge",
                 palette = set_colors_Period,
                 #  ellipse = TRUE, slipse.type="norm",
                 size=6,alpha=0.9#,
                # ellipse.border.remove = TRUE,
                 # ellipse.alpha = 0.2,
                
)+theme_bw()+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=0.5) +
  scale_shape_manual(values = set_shapes_recharge)+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  theme(legend.title=element_text(size=14, face="bold"))+
  geom_errorbar(aes(ymin = Axis.2.mean-Axis.2.se,ymax = Axis.2.mean+Axis.2.se, color=factor(Period)), width=0.01) +
  geom_errorbarh(aes(xmin = Axis.1.mean-Axis.1.se,xmax = Axis.1.mean+Axis.1.se,  color=factor(Period)), height = 0.01) +
  
  # geom_text_repel(aes(label = substr(rownames(tbd), 8, 11), color=PeriodRecharge), size=3) +
  xlab(paste0("PCoA1 (", round(scree$data[1,2]*100,1), "%)"))+
  ylab(paste("PCoA2 (", round(scree$data[2,2]*100,1), "%)", sep=""))+
  geom_line(aes(group = Period
                , colour = Period
                ), lty = 1, size=1.25, 
             color="black", alpha=0.6,
            arrow = arrow(length = unit(0.04, "npc")))+
  annotate("label", label = paste0("adonis(bray ~ Recharge + Period + Recharge*Period)\n", "p < ", ht_RechargePeriod$aov.tab$"Pr(>F)"[1], ", R-squared = ", 
                                   round((ht_RechargePeriod$aov.tab$R2[1]+ht_RechargePeriod$aov.tab$R2[2]+ht_RechargePeriod$aov.tab$R2[3]),2)), x = -0.2, y = -0.15, size = 4.5, 
           colour = "black",hjust = 0)+
  ggtitle("H52: PCoA_Period_recharge")

tp_final #5*6

ggsave("H52/bray/pcoa_H52_PeriodRecharge_mean_linked.pdf", width = 7, height = 5)


