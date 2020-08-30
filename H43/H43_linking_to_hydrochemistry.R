rm(list = ls())
setwd("/home/yo39qed/time-series analysis/output")
#saveRDS()
library(phyloseq)
library(reshape2)
library(readr)
require(gridExtra)
library(ggplot2)
library(magrittr)
library(digest)
library(corrplot)

set_colors_well <- c(
 "H13"= "green1", "H14" ="green4" , # Well H1 
  # "#A2A980", "wheat", # Well H2
"H31"=  "cyan1","H32"= "cyan4", # well H3
"H43"=  "#A2D9F7", "H42"= "lightgoldenrod1","H43" = "lightgoldenrod4",# well H4
"H51"=  "#FBC58C","H52"= "mediumpurple3","H53"= "plum" # Well H5
) 
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

set_colors <- c("Alphaproteobacteria" = "purple4", "Gammaproteobacteria"="mediumpurple2", "Deltaproteobacteria" = "mediumorchid2", "Proteobacteria_unclassified"="magenta1",
                "Actinobacteria"="brown","Thermoleophilia"="brown3","Acidimicrobiia"="brown4","Actinobacteria_unclassified"="brown2",
                "Bacteria_unclassified"="gray90",
                "Bacteroidia"="orange1", "Ignavibacteria" = "yellow",
                "Dehalococcoidia"="lightpink", "KD4-96"="hotpink", "Anaerolineae" = "pink1",  "JG30-KF-CM66"="pink2", # Chloroflexi
                "Elusimicrobia" = "burlywood","Elusimicrobia_unclassified" ="yellow4",
                "Nitrospira"="blue", "Thermodesulfovibrionia"="royalblue",  "HDB-SIOI1093"="blue3", # Nitrospirae
                "Saccharimonadia" = "green4", "Parcubacteria"="green2", "Gracilibacteria"="darkolivegreen4", "ABY1"="lightgreen", "Patescibacteria_unclassified" ="greenyellow", "Berkelbacteria"="#5CA595",
                "Omnitrophicaeota_cl"="darkorange2", "Omnitrophia"="chocolate",
                "Brocadiae"="firebrick2","OM190"="brown1", "Phycisphaerae"="orangered3",
                "Others"="gray",
                "Gemmatimonadetes" = "black",
                "Chlamydiae"="bisque",
                "Kiritimatiellae"="chocolate2",
                "Lineage_llc"="blue2",  "Verrucomicrobiae"="red2","Deinococci" ="yellow2",
                "MVP-15" ="darkturquoise",
                "Blastocatellia_(Subgroup_4)"="cyan2","Holophagae"="cornflowerblue", "Subgroup_6"="cyan1", # acidobacteria
                "Planctomycetacia" = "deepskyblue4",
                "vadinHA49"="mediumpurple4",
                "Nitrospinia"="violetred4",
                "Clostridia"="gray", # Firmicutes
                "NC10"="yellow2", # Rokubacteria 
                "Margulisbacteria_cl"="yellow4",
                "Sericytochromatia"="lightcyan", 
                "Babeliae"= "grey70",
                "Cyanobacteria_unclassified" ="#00FFBF"
                
)

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# Load data
load("/home/yo39qed/time-series analysis/output/H43/data_norm43.phyloseq")


###############################################################################
###                correlation between environemntal parameters             ###
###############################################################################

dat_all <- data.frame(sample_data(data_norm43))
u <- dat_all # for all env samples

## Remove variables that contain too many missing or 0 values: e.g. NO2.1, FE, MN, FE.1, MN.1
u1 <- u[ , -which(names(u) %in% c("NO2_1","FE", "MN", "FE.1", "MN.1", "P", 
                                  #"DO", "WL_mamsl", # For mcluster 4
                                  #"FE2","NH4", # FOr mcluster 1
                                  "PO4"))]

# Remove the unnecessary variables
u2 <- u1[ , -which(names(u1) %in% c("EC25","air_temperature__minimum", "air_temperature__maximum", "sunshine_duration",
                                   "soil_temperature__minimum","relative_humidity__minimum","relative_humidity__mean","relative_humidity__maximum",
                                   "wind_speed__mean","gust_speed__maximum",
                                   "Chao1","se.chao1", "ACE","se.ACE", "Simpson", "InvSimpson",
                                   "Fisher","WL_mbgl", "air_temperature__mean", "precipitation__sum"))]


# select the necessary variables
df <-subset(u2, select=EC:qPCR)

variable.names(df)
df <- df[,-which(names(df)%in% c("Year", "Month" , "Recharge", "Period", "PeriodRecharge", "PD", "Observed", "YearMonth"))]

#### Change the variable names
variable.names(df)
colnames(df)[which(names(df) == "PH")] <- "pH"
colnames(df)[which(names(df) == "ORP_NHE")] <- "ORP"
colnames(df)[which(names(df) == "FE2")] <- "Iron2"
colnames(df)[which(names(df) == "NO2")] <- "Nitrite"
colnames(df)[which(names(df) == "NH4")] <- "Ammonium"
colnames(df)[which(names(df) == "NO3")] <- "Nitrate"
colnames(df)[which(names(df) == "PO4")] <- "Orthophosphate"
colnames(df)[which(names(df) == "SO4")] <- "Sulphate"
colnames(df)[which(names(df) == "Cl")] <- "Chloride"
colnames(df)[which(names(df) == "CA")] <- "Calcium"
colnames(df)[which(names(df) == "K")] <- "Potassium"
colnames(df)[which(names(df) == "MG")] <- "Magnesium"
colnames(df)[which(names(df) == "NA.")] <- "Sodium"
colnames(df)[which(names(df) == "P")] <- "Phosphorus"
colnames(df)[which(names(df) == "S")] <- "Sulfur"
colnames(df)[which(names(df) == "WL_mamsl")] <- "Water_level"
#colnames(df)[which(names(df) == "precipitation__sum")] <- "precipitation"
colnames(df)[which(names(df) == "Temp")] <- "Water_temperature"


#### correlation plot (Pearson correlation)

p.mat <- cor.mtest(df)$p
corr_mat=round(cor(df, use = "pairwise.complete.obs"),2)

write.table(p.mat, "H43/env_correlation_pvalue_H43.csv", sep=";", dec=",", col.names = NA)
write.table(corr_mat, "H43/env_correlation_pearson_H43.csv", sep=";", dec=",", col.names = NA)

corrplot(corr_mat, type = "upper", order = "hclust", mar = c(0,0,0,0), method = "square", outline = T, 
         addgrid.col = "darkgray", addrect = 8, rect.col = "black", rect.lwd = 5, tl.col = "black", 
         tl.cex = 0.75, cl.cex = 0.75,
         p.mat = p.mat, sig.level = 0.05, insig = "blank")

dev.copy(pdf,"H43/env_pearson_H43.pdf", width = 6, height = 6)
dev.off()

#############################################################################
##                       iv) Constrained Ordinations                       ##
#############################################################################
# constrained ordination: to see how environmental variables are associated with these changes in community composition. 
# 1) We constrain the ordination axes to linear combinations of environmental variables. 
# 2) We then plot the environmental scores onto the ordination

########################### CAP ##########################################
# see: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

# Remove water temp and water level (because large data missing for H42 and H43)
# remove Water_temperature, since it didn't affect the community compostion but got many missing values
# remove shannon, sodium

u3 <- cbind(u2[,which(names(u2)%in% c("Samplename",  "Date",  "Sample.ID",  "Sample", "Sam_cam", "Well",
                                      "Date_Ref" ,"Year", "Month","YearMonth", "Recharge", "Period",        
                                      "PeriodRecharge" ))],df[,-which(names(df)%in% c("Water_temperature" ,  "Water_level", "Shannon" ))])


# write.table(u3, "H43/env_cca_H43.csv", sep=";", dec=",", col.names=NA)

df1 <- na.omit(u3) # 74 reduced to 46 samples (with water temp), 53 (without water temp)

my_data_prop_sqrt_not_na <- data_norm43

sample_data(my_data_prop_sqrt_not_na)<-df1
# Don't do this. It filtered out too many samples

my_data_prop_sqrt_not_na <- prune_taxa(taxa_sums(my_data_prop_sqrt_not_na) > 0, my_data_prop_sqrt_not_na)

bray_not_na <- phyloseq::distance(physeq = my_data_prop_sqrt_not_na, method = "bray")

############### 1) Environmental variables #########################
##### CAP for quantitative data #####
# Run DCA (Detrended correspondence analysis) to determine whether to use RDA or CCA
dca <- ordinate(physeq = my_data_prop_sqrt_not_na, 
                method = "DCA")
dca # or summary(dca) to check the "Axis lengths" of the DCA1 (Important) # 5.5844 or 4.0760 (after removing NA)
# If this value < 3, it is better to use RDA
# If this value > 4, it is better to use CCA
# If this value is between 3 and 4, either use CCA or RDA
# Here we use RDA for all wells

####################################################################################

##  Envfit + CCA/RDA: much more powerful when more complex factors are tested  ##

####################################################################################

library(vegan)
otutable <- vegan_otu(my_data_prop_sqrt_not_na) # import otu table
variable.names(df1)
dat <- data.frame(scale(sample_data(my_data_prop_sqrt_not_na)[,-c(1:12)])) # remove categorical data and Shannon

################# variance inflation factors ################### 
# Linear dependencies between constraints can be investigated via the variance inflation factor or VIF
# VIF is a measure of how much the variance of $\hat{\beta}_j$ is inflated by presence of other covariates
# Lots of rules of thumb
# VIF >= 20 indicates strong collinearity in constraints
# VIF >= 10 potnetially of concern & should be looked at
# they will be completely removed from the estimation, and no biplot scores or centroids are calculated for these aliased constraints. 

# vif.cca(ord)
ord_all <- capscale(bray_not_na~., dat)
temp <- vif.cca(ord_all)
temp
select_para <- names(temp[temp < 10])
select_para
# test all environmental factors

##### 2. RDA #####
dat1 <- dat[,select_para]
sample_data(my_data_prop_sqrt_not_na)<- dat1

bray_not_na <- phyloseq::distance(physeq = my_data_prop_sqrt_not_na, method = "bray")

cap_ord <- ordinate(
  physeq = my_data_prop_sqrt_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~.)

# Significance test
# anova(cap_ord, perm=999)
set.seed(123)
anova(cap_ord, by="margin", perm=999) # marginal effects of the terms 
drop <- drop1(cap_ord, test="perm")

#### significant factor ####
select_para1 <- row.names(subset(data.frame(drop), Pr..F. < 0.05))

###############################################################################################################
############################ prepare final model ##########################################################
# remove NA values (less env parameter, more samples are left)
u4 <- u3[,select_para1]
df2 <- na.omit(u4) # 74
df3 <- data.frame(scale(df2)) # remove categorical data and Shannon

my_data_prop_sqrt_not_na1 <- data_norm43

sample_data(my_data_prop_sqrt_not_na1)<-df3
# Don't do this. It filtered out too many samples

my_data_prop_sqrt_not_na1 <- prune_taxa(taxa_sums(my_data_prop_sqrt_not_na1) > 0, my_data_prop_sqrt_not_na1)

bray_not_na1 <- phyloseq::distance(physeq = my_data_prop_sqrt_not_na1, method = "bray")


cap_ord <- ordinate(
  physeq = my_data_prop_sqrt_not_na1, 
  method = "CAP",
  distance = bray_not_na1,
  formula = ~ .)

# Significance test
ht <- anova(cap_ord, perm=999)
set.seed(123)
ht_env <- anova(cap_ord, by="margin", perm=999) # marginal effects of the terms 
anova(cap_ord, by="terms") # sequential
drop1(cap_ord, test="perm")
ordistep(cap_ord, perm.max = 999)
summary(cap_ord)
RsquareAdj(cap_ord) # ## RsquareAdj gives the same result as component [a] of varpart

################### Biplot #####################

# CAP plot
sample_data(my_data_prop_sqrt_not_na1) <- u # add categorical data for plotting

cap_plot_env <- plot_ordination(
  physeq = my_data_prop_sqrt_not_na1, 
  ordination = cap_ord, 
# color="mCluster",
 color="Period", 
 shape="Recharge",
  axes = c(1,2))+ 
  theme_bw() +
 # geom_point(aes(colour = mCluster), alpha = 0.8, size = 4) +
  geom_point(aes(colour = Period), alpha = 0.8, size = 4, stroke = 1) +
 # scale_color_manual(values = set_colors_mcluster) +
  scale_color_manual(values = set_colors_Period) +
  scale_shape_manual(values = set_shapes_recharge) +
 #stat_ellipse(aes(fill=Recharge),type = "norm") +
  scale_fill_manual(values = set_colors_recharge)+
  annotate("label", label = paste0("P < ", ht$"Pr(>F)"[1],"\n" , "adj. R-squared = ", 
                                   round(RsquareAdj(cap_ord)$adj.r.squared*100,1), "%"), 
           x = 0.65, y = -1.5, size = 6, fill="transparent",
           colour = "black",hjust = 0)+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14))
  
cap_plot_env
cap_plot_env$layers<-cap_plot_env$layer [-1]

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1*1.4, 
                 yend = CAP2*1.4, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.7 * CAP1, 
                 y = 1.7 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
p_H43 <- cap_plot_env + 
  geom_segment(
    mapping = arrow_map, 
    size = 0.75, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(mapping = label_map,  size = 6, color="black",  data = arrowdf,  show.legend = FALSE)+
  ggtitle("H43: Linking to environmental variables") +
 # xlim(-1.5, NA)+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=0.5) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=0.5) +
 # scale_y_continuous(breaks = c(-1.5, -1, -0.5,0,  0.5, 1, 1.5))+
  theme(legend.key = element_blank(),  #removes the box around each legend item
        legend.position = "right")+ #legend at the bottom
  theme(axis.title=element_text(size=18), 
        axis.text=element_text(size=18), legend.text=element_text(size=18)) +
  theme(legend.title=element_text(size=18, face="bold"))
p_H43 #8*7

ggsave("H43/cap_linking_env_H43.pdf", width = 9, height = 7)


