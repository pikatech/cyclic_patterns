###########################################################################
#########                1. Data preparation                     ##########
###########################################################################
setwd("/home/yo39qed/time-series analysis")

### Import another metadata file which only contains the metadata of H52 well
load("my_data_dna.phyloseq") # we have already produced this file when we prepared for H41 data
metadata <- read.csv("input/metadata/output/metadata_H52_redo.csv", sep=";", dec = ",", row.names=1, header=TRUE)

my_data52 <- my_data_dna
sample_data(my_data52) <- metadata

my_data52 <- filter_taxa(my_data52, function(x) sum(x) > 0, TRUE)
levels(sample_data(my_data52)$Sam_cam)

sample_data(my_data52)$Year <- format(as.Date(sample_data(my_data52)$Date_Ref, format="%d/%m/%Y"),"%Y")
sample_data(my_data52)$Month <- format(as.Date(sample_data(my_data52)$Date_Ref, format="%d/%m/%Y"),"%m")
sample_data(my_data52)$Year <- as.factor(sample_data(my_data52)$Year)
sample_data(my_data52)$Month <- as.factor(sample_data(my_data52)$Month)

####################### Data normalization #####################
# check: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# Estimating alpha diversity of microbial communities is problematic no matter what you do. 
# My best stab at it is to subsample the libraries with replacement to estimate the species abundance of the real population while standardizing sampling effort.
sort(sample_sums(my_data52)) 
min_lib <- min(sample_sums(my_data52))
min_lib # 2502 for all samples, if we wanna compare

####################################################
# MetagenomeSeq: CSS normalization 
# "Remove the OTUs, of which the sum across all samples is greater than 0.005% of all OTUs based on Bokulich et al. (2013)."

minTotRelAbun = 1e-5
x = taxa_sums(my_data52)
keepTaxa = taxa_names(my_data52)[which((x / sum(x)) > minTotRelAbun)]
prunedSet = prune_taxa(keepTaxa, my_data52)
sort(sample_sums(prunedSet))

test <- phyloseq_to_metagenomeSeq(prunedSet)
p <- cumNormStatFast(test)

## Default value being used.
alphadata <- cumNorm(test, p = p)

# Export normalised and not-log matrix for use 
otu_norm <- MRcounts(alphadata, norm = TRUE, log = FALSE)

data_norm52 <- prunedSet
otu_table(data_norm52) <- otu_table(otu_norm, taxa_are_rows = T)

# CCS somehow fixed the variation in sampling depth. 
# CSS will sometimes decrease the fold-difference in sampling depth but not always. 
max(sample_sums(my_data52))/min(sample_sums(my_data52)) #  70.43006
max(sample_sums(data_norm52))/min(sample_sums(data_norm52)) # 8.869802


############################
# Estimate alpha diversity #
############################

setwd("/home/yo39qed/time-series analysis")
alpha <- estimate_richness(my_data52)
metadata <- sample_data(my_data52)

library(plyr)
alpha$Samplename <- row.names(alpha)
metadata$Samplename <- row.names(metadata)
alphadiv <- merge(data.frame(metadata),data.frame(alpha),
                  by = "Samplename", all = TRUE)
row.names(alphadiv) <- alphadiv$Samplename
write.table(alphadiv, "output/H52/alphadiv52_not_normalized.csv", sep=";", dec=",", col.names=NA) # The final metadata file with sequence =3000 /sample for depth for alphadiv estimation

sample_data(my_data52) <- alphadiv

##### save all DNA raw H52 data 
save(my_data52, file ="my_data52.phyloseq")

sample_data(data_norm52) <- alphadiv

####################################
########  phylogenetic tree ######## 
####################################

# get ref sequences
ref <- refseq(data_norm52)
library(Biostrings)
writeXStringSet(ref, "output/H52/rep_seq_CSS.fasta")

######## Calculate the phylogenetic tree in mothur ######

#dist.seqs(fasta=rep_seq_CSS.fasta,output=lt, processors=4)
#clearcut(phylip=current)
#quit()

#### integrate the tree file into the phyloseq file
tree <- read.tree("output/H52/rep_seq_CSS.phylip.tre") 
# plot(basinbactree) # this step will take very long time due to too many taxa
phytree <- phy_tree(tree)

data_norm52 <- merge_phyloseq(data_norm52, phytree)
save(data_norm52, file ="output/H52/data_norm52.phyloseq")

#####################################################################
##### Phylogenetic diversity (Faiths PD-index) ######################
library(picante)
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}

pd <- estimate_pd(data_norm52)
sample_data(data_norm52)$PD <- pd$PD

save(data_norm52, file ="output/H52/data_norm52.phyloseq")

#################################################################
### ##### Now define the recharge and Discharge period ###### ###
#################################################################

# Based on water level
metadata <- sample_data(data_norm52)

metadata$YearMonth<-paste0(metadata$Year,metadata$Month)
metadata$Recharge <- as.factor(metadata$YearMonth)
levels(metadata$Recharge)

metadata$Recharge <- revalue(metadata$Recharge, c( "201305"="Recharge","201306"="Recharge", 
                                                  
                                                  "201307"="Discharge", "201308"="Discharge","201309" ="Discharge", 
                                                  "201310"="Discharge", "201311"="Discharge", "201312"="Discharge",
                                                  "201401"="Discharge", "201402"="Discharge","201403"="Discharge",
                                                  "201404"="Discharge","201405"="Discharge",
                                                  "201406"="Discharge", "201407"="Discharge","201408"="Discharge",
                                                  "201409"="Discharge","201410"="Discharge","201411"="Discharge",
                                                  "201412"="Discharge",

                                                  "201501"="Recharge", "201502"="Recharge","201503"="Recharge",
                                                  "201504"="Recharge", "201505"="Recharge",
                                                  
                                                  "201506"="Discharge", "201507"="Discharge", "201508"="Discharge",
                                                  "201509"="Discharge","201510"="Discharge","201511"="Discharge",
                                                  
                                                  "201512"="Recharge","201601"="Recharge", "201602"="Recharge",
                                                  "201603"="Recharge","201604"="Recharge","201605"="Recharge",
                                                  
                                                  "201606"="Discharge","201607"="Discharge", "201608"="Discharge",
                                                  "201609"="Discharge","201610"="Discharge","201611"="Discharge",
                                                  "201612"="Discharge","201701"="Discharge","201702"="Discharge",
                                                  
                                                  "201703"="Recharge","201704"="Recharge","201705"="Recharge",
                                                  
                                                  "201706"="Discharge", "201707"="Discharge", "201708"="Discharge", 
                                                  "201709"="Discharge","201710"="Discharge",
                                                  
                                                  "201711"="Recharge",
                                                  "201712"="Recharge","201801"="Recharge", "201802"="Recharge",
                                                  "201803"="Recharge","201804"="Recharge", "201805"="Recharge",
                                                  
                                                  "201806"="Discharge",
                                                  "201807"="Discharge", "201808"="Discharge","201809"="Discharge",
                                                  "201810"="Discharge","201811"="Discharge","201812"="Discharge",
                                                  
                                                  "201901"="Recharge",
                                                  "201902"="Recharge","201903"="Recharge","201904"="Recharge",
                                                  "201905"="Recharge"
                                                  
))

metadata$Period <- as.factor(metadata$YearMonth) 
metadata$Period <- revalue(metadata$Period, c("201305"="P13","201306"="P13", 
                                              
                                              "201307"="P15", "201308"="P15","201309" ="P15", 
                                              "201310"="P15", "201311"="P15", "201312"="P15",
                                              "201401"="P15", "201402"="P15","201403"="P15",
                                              "201404"="P15","201405"="P15",
                                              "201406"="P15", "201407"="P15","201408"="P15",
                                              "201409"="P15","201410"="P15","201411"="P15",
                                              "201412"="P15",
                                              
                                              "201501"="P15", "201502"="P15","201503"="P15",
                                              "201504"="P15", "201505"="P15",
                                              
                                              "201506"="P16", "201507"="P16", "201508"="P16",
                                              "201509"="P16","201510"="P16","201511"="P16",
                                              
                                              "201512"="P16","201601"="P16", "201602"="P16",
                                              "201603"="P16","201604"="P16","201605"="P16",
                                              
                                              "201606"="P17","201607"="P17", "201608"="P17",
                                              "201609"="P17","201610"="P17","201611"="P17",
                                              "201612"="P17","201701"="P17","201702"="P17",
                                              
                                              "201703"="P17","201704"="P17","201705"="P17",
                                              
                                              "201706"="P18", "201707"="P18", "201708"="P18", 
                                              "201709"="P18","201710"="P18",
                                              
                                              "201711"="P18",
                                              "201712"="P18","201801"="P18", "201802"="P18",
                                              "201803"="P18","201804"="P18", "201805"="P18",
                                              
                                              "201806"="P19",
                                              "201807"="P19", "201808"="P19","201809"="P19",
                                              "201810"="P19","201811"="P19","201812"="P19",
                                              
                                              "201901"="P19",
                                              "201902"="P19","201903"="P19","201904"="P19",
                                              "201905"="P19"
))

metadata$PeriodRecharge<-paste0(metadata$Period,metadata$Recharge) 

# Change PNK089: to Period4Discharge
reP15 <- subset(metadata, Sam_cam!="PNK089") # 30/05/2017, it is not recharge at all
reP16 <- subset(metadata, Sam_cam=="PNK089")
reP16[,c("Recharge", "Period", "PeriodRecharge")] <- c("Discharge",  "P18", "P19Discharge")
metadata <- rbind(reP15, reP16)

# Manual correction -- > Change PNK102: to Period5Discharge
reP15 <- subset(metadata, Sam_cam!="PNK102") # 30/05/2018, it is not recharge at all
reP16 <- subset(metadata, Sam_cam=="PNK102")
reP16[,c("Recharge", "Period", "PeriodRecharge")] <- c("Discharge",  "P19", "P19Discharge")
metadata <- rbind(reP15, reP16)

sample_data(data_norm52)<-metadata

###### Add qPCR data ####
# Add qpcr data
setwd("/home/yo39qed/time-series analysis/input")
qpcr <- read.csv("qPCR.csv", sep=";", dec = ",",  header=TRUE)

setwd("/home/yo39qed/time-series analysis/output")
qpcr_H52 <- qpcr[,c("Samp_Cam", "H52")]
colnames(qpcr_H52) <-  c("Sam_cam", "qPCR")
summary1 <- data.frame(sample_data(data_norm52))

summary <- merge(data.frame(summary1), data.frame(qpcr_H52),
                 by = "Sam_cam", all=T) 


summary<-summary[!is.na(summary$Samplename),] # remove empty rows with NA values (samples with missing microbial data)

row.names(summary)<-summary$Samplename

sample_data(data_norm52)<-summary

save(data_norm52, file ="H52/data_norm52.phyloseq")

