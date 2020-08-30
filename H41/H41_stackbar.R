rm(list = ls())

setwd("/home/yo39qed/time-series analysis/output")

#saveRDS()
library(phyloseq)
require(gridExtra)
library(ggplot2)
library(magrittr)
library(lubridate)
library(ggalluvial)
library(dplyr)

#####################
set_colors <- c("Alphaproteobacteria" = "purple4", "Gammaproteobacteria"="mediumpurple2", "Deltaproteobacteria" = "mediumorchid2", "Proteobacteria_unclassified"="magenta1",
                "Actinobacteria"="brown","Thermoleophilia"="brown3","Acidimicrobiia"="brown4","Actinobacteria_unclassified"="brown2",
                "Bacteria_unclassified"="gray90",
                "Bacteroidia"="orange1", "Ignavibacteria" = "yellow",
                "Dehalococcoidia"="lightpink", "KD4-96"="hotpink", "Anaerolineae" = "pink1",  "JG30-KF-CM66"="pink2", # Chloroflexi
                "Elusimicrobia" = "burlywood","Elusimicrobia_unclassified" ="yellow4",
                "Nitrospira"="blue", "Thermodesulfovibrionia"="royalblue",  "HDB-SIOI1093"="blue3", # Nitrospirae
                "Saccharimonadia" = "green4", "Parcubacteria"="green2", "Gracilibacteria"="darkolivegreen4", "ABY1"="lightgreen", "Patescibacteria_unclassified" ="greenyellow", "Berkelbacteria"="#5CA595",
                "Omnitrophicaeota_cl"="darkorange2", 
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
                "NC10"="orange1", # Rokubacteria 
                "Margulisbacteria_cl"="yellow4",
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
load("/home/yo39qed/time-series analysis/output/H41/data_norm41.phyloseq")

###########################################################################
#########                    2. Stack bar                        ##########
###########################################################################

my_data_prop <- transform_sample_counts(data_norm41, function(x) 100*x/sum(x))
#my_data_prop <- transform_sample_counts(my_data41, function(x) 100*x/sum(x))
#my_data_prop <- transform_sample_counts(prunedSet, function(x) 100*x/sum(x))


RANK="Class"# relative abundance
my_data_prop_rank <- tax_glom(my_data_prop, taxrank=RANK)

sample_data(my_data_prop_rank)$Filter <- "Filter2"
my_data_prop_abun = filter_taxa(my_data_prop_rank, function(x) mean(x) > 1, TRUE)

# select class as rank
## Calculate the abundance sum of other rare taxa, We want to calculate a new phyla called "Others"
others <- tax_glom(my_data_prop_abun, "Kingdom" )
otu_table(others) <- 100 - otu_table(others)
tax_table(others)@.Data[,2:6] <- "Others"# Define all taxanomic levels as "Others"
taxa_names(others) <- "Others" # Define the taxa name as "Others"

# Combine the abundant taxa with the sum of rare taxa
OTU1 <- otu_table(t(rbind(otu_table(others), otu_table(my_data_prop_abun) )), taxa_are_rows = FALSE)
TAX1 <- tax_table(rbind(tax_table(others), tax_table(my_data_prop_abun) )) 
METADATA1 <- sample_data(my_data_prop_abun)

final <- phyloseq(OTU1, TAX1, METADATA1)

# save(final, file ="H41/final_H41.phyloseq")

##### melt ####
my_data_rank_melt <- psmelt (final)
#nsamples(my_data_rank)
# order the rank by abundance
my_data_rank_melt <- transform(my_data_rank_melt, Phylum=reorder(Class, order(Abundance, decreasing=TRUE)))
write.table(my_data_rank_melt, "H41/Class_over1percent_per_well.csv", sep=";", dec=",", col.names=NA)
#me <- aggregate(Abundance ~ Well+Phylum, my_data_rank_melt, function(x) c(mean = mean(x), sd = sd(x)))

############################## 1) simple figure #################################
############# Used ######## sampling campaign ########################

my_data_rank_melt$Class<- as.factor(my_data_rank_melt$Class)
my_data_rank_melt$Date_Ref <- format(as.Date(my_data_rank_melt$Date, origin="1899-12-30"), "%m/%d/%Y")
my_data_rank_melt$Date_Ref <- as.Date(my_data_rank_melt$Date_Ref, "%m/%d/%Y")
my_data_rank_melt$Year <- format(as.Date(my_data_rank_melt$Date_Ref, format="%d/%m/%Y"),"%Y")
my_data_rank_melt$Month <- format(as.Date(my_data_rank_melt$Date_Ref, format="%d/%m/%Y"),"%m")
my_data_rank_melt$Year <- as.factor(my_data_rank_melt$Year)
my_data_rank_melt$Month <- as.factor(my_data_rank_melt$Month)

my_data_rank_melt_select <- my_data_rank_melt
my_data_rank_melt_select [is.na(my_data_rank_melt_select )] <- 0 
# This is important, empty cells will cause failure of Alluvial plots

# Bar
my_data_rank_melt_select$Class<- factor(my_data_rank_melt_select$Class, levels = 
                                          c("Others", "Bacteria_unclassified",
                                            "Parcubacteria", "ABY1",  "Saccharimonadia" , "Gracilibacteria", 
                                            "Alphaproteobacteria" , "Gammaproteobacteria", "Deltaproteobacteria" , 
                                            "Nitrospira", "Thermodesulfovibrionia", 
                                            "Brocadiae","Phycisphaerae",
                                            "Actinobacteria", 
                                            "Bacteroidia", 
                                            "KD4-96","Dehalococcoidia",
                                            "Elusimicrobia",
                                            "Omnitrophicaeota_cl", 
                                            "Gemmatimonadetes" 
                                            ))



##############################################################################################

## 2) stacked area (Date vs abundance)

#################################################################################################


a_H41 <- ggplot(my_data_rank_melt_select,
       aes(x = Date_Ref,  y=Abundance,
           fill = Class, label = Class)) +
  scale_fill_manual(values=set_colors) +
  geom_area(alpha=0.6, size=0.5, color="white") +
  # facet_wrap(Well~Year)+
  theme_bw()+
  # geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "bottom") +
  ylab("Relative abundance (%)") +
  xlab("Date") +
  scale_x_date(limits = as.Date(c('18/2/2013', '29/5/2019'), format="%d/%m/%Y"),
               date_breaks = "2 months", date_minor_breaks = "1 month", 
               labels = date_format("%b %Y") )+
  guides(fill=guide_legend(nrow=5))+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=12),legend.title=element_text(size=14)) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0.5))+
  ggtitle("Bacterial community changes over time")

a_H41
ggsave("H41/stacked_area_class_H41_mothur-unsorted_norm.pdf", width = 9.45, height = 7)



##############################################################################################
#### calculate mean and sd of relative abundance

# Assume `ps` is a phyloseq object that has been glommed to the Family level
# and converted to proportions.
tb_OTU <- my_data_prop %>%
  # tax_glom(taxrank="Genus") %>%
  psmelt() %>%
  group_by(OTU, Genus, Family, Order, Class, Phylum) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
write.table(tb_OTU, "H41/mean_SD_relativE_abundance_OTU.csv", sep=";", dec=",", col.names=NA)

tb_genus <- my_data_prop %>%
  tax_glom(taxrank="Genus") %>%
  psmelt() %>%
  group_by(Genus, Family, Order, Class, Phylum) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
write.table(tb_genus, "H41/mean_SD_relativE_abundance_genus.csv", sep=";", dec=",", col.names=NA)

tb_family <- my_data_prop %>%
  tax_glom(taxrank="Family") %>%
  psmelt() %>%
  group_by(Family, Order, Class, Phylum) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
write.table(tb_family, "H41/mean_SD_relativE_abundance_family.csv", sep=";", dec=",", col.names=NA)

tb_order <- my_data_prop %>%
  tax_glom(taxrank="Order") %>%
  psmelt() %>%
  group_by(Order, Class, Phylum) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
write.table(tb_order, "H41/mean_SD_relativE_abundance_order.csv", sep=";", dec=",", col.names=NA)

tb_class <- my_data_prop %>%
  tax_glom(taxrank="Class") %>%
  psmelt() %>%
  group_by(Class, Phylum) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
write.table(tb_class, "H41/mean_SD_relativE_abundance_class.csv", sep=";", dec=",", col.names=NA)


tb_phylum <- my_data_prop %>%
  tax_glom(taxrank="Phylum") %>%
  psmelt() %>%
  group_by(Phylum) %>%
  summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
write.table(tb_phylum, "H41/mean_SD_relative_abundance_phylum.csv", sep=";", dec=",", col.names=NA)

