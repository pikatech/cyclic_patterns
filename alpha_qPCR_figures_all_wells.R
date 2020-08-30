setwd("/home/yo39qed/time-series analysis/output")
#saveRDS()
library(phyloseq)
require(gridExtra)
library(ggplot2)
library(ggpubr)

set_colors_Year <- get_palette("lancet", 7)

set_colors_Period <-c("P13"= "brown","P15"="red", "P16"="deeppink", "P17"="orange",
                      "P18"="purple", "P19"="blue", "P14"="darkgreen")

set_colors_PeriodRecharge <- c("P13Recharge" ="brown","P15Recharge" = "red", "P15Discharge" ="red", 
                               "P16Recharge" ="deeppink", "P16Discharge"="deeppink", "P17Recharge"="orange",    
                               "P17Discharge"="orange", "P18Recharge" ="purple",  "P18Discharge"= "purple",
                               "P19Recharge"= "blue", 
                               "P19Discharge"= "blue", "P14Recharge"="darkgreen", "P14Discharge"="darkgreen")


set_colors_recharge <- c("Recharge"="orange", "Discharge"="grey")
set_shapes_recharge <- c("Recharge"=15, "Discharge"=1)


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
load("/home/yo39qed/time-series analysis/output/H41/data_norm41.phyloseq")
load("/home/yo39qed/time-series analysis/output/H41/data_norm43.phyloseq")
load("/home/yo39qed/time-series analysis/output/H41/data_norm52.phyloseq")

#######   alpha diversity  analysis ##########

summary_H41 <- data.frame(sample_data(data_norm41))
summary_H41$Recharge <- factor(summary_H41$Recharge, levels = c("Discharge", "Recharge"))

summary_H43<- data.frame(sample_data(data_norm43))
summary_H43$Recharge <- factor(summary_H43$Recharge, levels = c("Discharge", "Recharge"))

summary_H52 <- data.frame(sample_data(data_norm52))
summary_H52$Recharge <- factor(summary_H52$Recharge, levels = c("Discharge", "Recharge"))

summary <- rbind(summary_H41, summary_H43, summary_H52)
write.table(summary, "summary_univariate_Wells.csv", sep=";", dec=",", col.names = NA)


#################################### Shannon #################################
# 1) Recharge
res <- aggregate(Shannon ~ Recharge + Well, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res
write.table(res, "mean_Shannon_recharge_Wells.csv", sep=";", dec=",", col.names = NA)

p_recharge <- ggboxplot(summary, x="Recharge", y = "Shannon", fill = "Recharge",
                        palette = set_colors_recharge,
                        add = "mean" , facet.by = "Well"
                        #shape = "Aquifer",
)+ 
  # stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  #ylim(0,10) +
  theme_bw()+
  stat_compare_means(label.x = 1, label.y = 9, size=5) +
  # facet_wrap(.~Period,scales = "fixed", nrow = 1)+
  theme(strip.text = element_text(size=14))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) 

p_recharge
ggexport(p_recharge, filename = 'shannon_Recharge_wells.pdf', width = 8, height = 4)

##############################  Recharge + period  ########################
res <- aggregate(Shannon ~ Recharge + Period + Well, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res
write.table(res, "mean_Shannon_recharge_period_Wells.csv", sep=";", dec=",", col.names = NA)

res <- aggregate(Shannon ~ Period + Well, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res
write.table(res, "mean_Shannon_period_Wells.csv", sep=";", dec=",", col.names = NA)

p_final <- ggboxplot(summary, x="Recharge", y = "Shannon", fill = "Recharge",
                     palette = set_colors_recharge,
                     add = "mean" 
                     #shape = "Aquifer",
)+ 
#  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  #ylim(0,10) +
  theme_bw()+
  stat_compare_means(label.x = 1, label.y = 4, size=4) +
  facet_grid(Well~Period,scales = "fixed")+
  theme(strip.text = element_text(size=14))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  +
 # geom_hline(yintercept=mean(summary$Shannon, na.rm=TRUE), linetype="dashed",  color = "black", size=1) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) 

p_final
ggexport(p_final, filename = 'shannon_PeriodRecharge_wells.pdf', width = 14, height = 7)

##################### Effect of Period on each well ###############################
# 1) H41
well <- "H41"
summary1 <- subset(summary, Well == well)
# To perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(summary1$Shannon)
# Use non-parametric method
kruskal.test(Shannon ~ Period, data = summary1) 
# anova_shannon <- aov(Shannon ~ Period + Recharge + Period*Recharge, data = summary1)
# summary(anova_shannon)
res3 <- aggregate(Shannon ~ Period, data=summary1, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res3
write.table(res3, paste0(well, "/Shannon_Period_", well, ".csv"), sep=";", dec=",", col.names = NA)

# Dunn test for multiple comparisons 
#Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations
#install.packages("FSA")
library(FSA)
PT<-dunnTest(Shannon~Period,summary1, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.unadj, # or PT$P.unadj
        threshold  = 0.05)

# 2) H43
well <- "H43"
summary1 <- subset(summary, Well == well)
# To perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(summary1$Shannon)
# Use non-parametric method
kruskal.test(Shannon ~ Period, data = summary1) 
# anova_shannon <- aov(Shannon ~ Period + Recharge + Period*Recharge, data = summary1)
# summary(anova_shannon)
res3 <- aggregate(Shannon ~ Period, data=summary1, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res3
write.table(res3, paste0(well, "/Shannon_Period_", well, ".csv"), sep=";", dec=",", col.names = NA)

# Dunn test for multiple comparisons 
#Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations
#install.packages("FSA")
library(FSA)
PT<-dunnTest(Shannon~Period,summary1, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.unadj, # or PT$P.unadj
        threshold  = 0.05)

# 3) H52
well <- "H52"
summary1 <- subset(summary, Well == well)
# To perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(summary1$Shannon)
# Use non-parametric method
kruskal.test(Shannon ~ Period, data = summary1) 
# anova_shannon <- aov(Shannon ~ Period + Recharge + Period*Recharge, data = summary1)
# summary(anova_shannon)
res3 <- aggregate(Shannon ~ Period, data=summary1, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res3
write.table(res3, paste0(well, "/Shannon_Period_", well, ".csv"), sep=";", dec=",", col.names = NA)

# Dunn test for multiple comparisons 
#Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations
#install.packages("FSA")
library(FSA)
PT<-dunnTest(Shannon~Period,summary1, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.unadj, # or PT$P.unadj
        threshold  = 0.05)


##################### Effect of Well ###############################
summary1 <- summary
# To perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(summary1$Shannon)
# Use non-parametric method
kruskal.test(Shannon ~ Well, data = summary1) 

res3 <- aggregate(Shannon ~ Well, data=summary1, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res3
write.table(res3, "Shannon_wells.csv", sep=";", dec=",", col.names = NA)

# Dunn test for multiple comparisons 
#Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations
#install.packages("FSA")
library(FSA)
PT<-dunnTest(Shannon~Well,summary1, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.unadj, # or PT$P.unadj
        threshold  = 0.05)


##################################################################################################################
#                                                     1) qPCR
##################################################################################################################

##############################  Recharge   ########################
res <- aggregate(qPCR ~ Recharge + Well, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res
write.table(res, "mean_qPCR_recharge_wells.csv", sep=";", dec=",", col.names = NA)

res <- aggregate(qPCR ~ Recharge + Period + Well, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res
write.table(res, "mean_qPCR_recharge_period_Wells.csv", sep=";", dec=",", col.names = NA)

res <- aggregate(qPCR ~ Period + Well, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res
write.table(res, "mean_qPCR_period_Wells.csv", sep=";", dec=",", col.names = NA)

p_recharge1 <- ggboxplot(summary, x="Recharge", y = "qPCR", fill = "Recharge",
                        palette = set_colors_recharge,
                        add = "mean" , facet.by = "Well"
                        #shape = "Aquifer",
)+ 
  # stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  #ylim(0,1200000000) +
  theme_bw()+
  stat_compare_means(label.x = 1, label.y = 9.3, size=5) +
  # facet_wrap(.~Period,scales = "fixed", nrow = 1)+
  theme(strip.text = element_text(size=14))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                      labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) 

p_recharge1
ggexport(p_recharge1, filename = 'qPCR_Recharge_wells.pdf', width = 8, height = 4)

##############################  Recharge + period  ########################

p_final <- ggboxplot(summary, x="Recharge", y = "qPCR", fill = "Recharge",
                     palette = set_colors_recharge,
                     add = "mean" 
                     #shape = "Aquifer",
)+ 
 # stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  #ylim(0,10) +
  theme_bw()+
  stat_compare_means(label.x = 1, label.y = 9.3, size=4) +
  facet_grid(Well~Period,scales = "fixed")+
  theme(strip.text = element_text(size=14))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) 

p_final
ggexport(p_final, filename = 'qPCR_PeriodRecharge_wells.pdf', width = 14, height = 7)

#dev.off()

### Pairwise comparison
library(FSA)
PT<-dunnTest(qPCR~PeriodRecharge,summary, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.adj, # or PT$P.unadj
        threshold  = 0.05)



################## Effect of Period ####
well <- "H52"
summary1 <- subset(summary, Well == well)
# To perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(summary1$qPCR)
# Use non-parametric method
kruskal.test(qPCR ~ Period, data = summary1) 
# anova_qPCR <- aov(qPCR ~ Period + Recharge + Period*Recharge, data = summary1)
# summary(anova_qPCR)
res3 <- aggregate(qPCR ~ Period, data=summary1, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res3
write.table(res3, paste0(well, "/qPCR_Period_", well, ".csv"), sep=";", dec=",", col.names = NA)

# Dunn test for multiple comparisons 
#Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations
#install.packages("FSA")
library(FSA)
PT<-dunnTest(qPCR~Period,summary1, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.unadj, # or PT$P.unadj
        threshold  = 0.05)

#############################################################
## Combine two figures
grid.arrange(p_recharge, p_recharge1, nrow=2)
dev.copy(pdf,"Shannon_qPCR_Period.pdf", width = 9, height =5)
dev.off()


### Effect of well
################## Effect of Period ####
summary1 <- summary
# To perform the Shapiro-Wilk test of normality for one variable (univariate):
shapiro.test(summary1$qPCR)
# Use non-parametric method
kruskal.test(qPCR ~ Well, data = summary1) 
# anova_qPCR <- aov(qPCR ~ Period + Recharge + Period*Recharge, data = summary1)
# summary(anova_qPCR)
res3 <- aggregate(qPCR ~ Well, data=summary1, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res3
write.table(res3, "qPCR_wells.csv", sep=";", dec=",", col.names = NA)

# Dunn test for multiple comparisons 
#Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations
#install.packages("FSA")
library(FSA)
PT<-dunnTest(qPCR~Well,summary1, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.unadj, # or PT$P.unadj
        threshold  = 0.05)

### Recharge
res2 <- aggregate(qPCR ~ PeriodRecharge, data=summary, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res2
write.table(res2, "H41/mean_qPCR_recharge_Period.csv", sep=";", dec=",", col.names = NA)
