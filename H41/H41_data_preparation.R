###########################################################################
#########                1. Data preparation                     ##########
###########################################################################

setwd("/home/yo39qed/time-series analysis")

library(phyloseq)
library(Biostrings)
library(reshape2)
library(ape) # to read tree/fasta file
library(readr)
library(vegan)
library(plyr)
library(picante)
library(metagenomeSeq);packageVersion("metagenomeSeq")

# Import data
my_data <- import_biom("input/Hainich6_final.pcr.filter.agc.0.03.biom", refseqfilename="input/Hainich6.rep.txt") #331 samples
colnames(tax_table(my_data))
colnames(tax_table(my_data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

### to check the reads
lib <- data.frame(sample_sums(my_data))
write.table(lib, "output/library size_all.csv", sep=";", col.names=NA)

##### save all data 
save(my_data, file ="aquadiva_all.phyloseq")

# To import metadata, this metadata contains no RNA, no H21 and H23 samples, 0.2 micrometer filter ####
metadata2 <- read.csv("input/general_metadata_updated_run22.csv", sep=";", dec = ",", row.names=1, header=TRUE)
# Add another variable as "Cluster"
attach(metadata2) 
metadata2$mCluster[Well == "H13" |Well == "H14" |Well == "H21" ] <- "mCluster 1"
metadata2$mCluster[Well == "H31" |Well == "H32"] <- "mCluster 2"
metadata2$mCluster[Well == "H51"] <- "mCluster 3"
metadata2$mCluster[Well == "H42" |Well == "H43" ] <- "mCluster 4"
metadata2$mCluster[Well == "H52" |Well == "H53" ] <- "mCluster 5"
metadata2$mCluster[Well == "H41"] <- "mCluster 6"
detach(metadata2)

my_data_dna <- my_data
sample_data(my_data_dna) <- metadata2

##### save all DNA data 
save(my_data_dna, file ="aquadiva_DNA.phyloseq")

sample_names(my_data_dna) <- sample_data(my_data_dna)$Sample
sample_data(my_data_dna)$Sample[duplicated(sample_data(my_data_dna)$Sample)] # check if there are duplicate cells

# library size of each sample
lib <- data.frame(sample_sums(my_data_dna))
write.table(lib, "output/library size_dna.csv", sep=";", col.names=NA)

### Import another metadata file which only contains the metadata of H41 well
metadata <- read.csv("input/metadata/output/metadata_H41_redo.csv", sep=";", dec = ",", row.names=1, header=TRUE)
my_data41 <- my_data_dna
sample_data(my_data41) <- metadata

# Remove samples that had no counts
my_data41 <- filter_taxa(my_data41, function(x) sum(x) > 0, TRUE)

# check the year and month of each sample
sample_data(my_data41)$Year <- format(as.Date(sample_data(my_data41)$Date_Ref, format="%d/%m/%Y"),"%Y")
sample_data(my_data41)$Month <- format(as.Date(sample_data(my_data41)$Date_Ref, format="%d/%m/%Y"),"%m")
sample_data(my_data41)$Year <- as.factor(sample_data(my_data41)$Year)
sample_data(my_data41)$Month <- as.factor(sample_data(my_data41)$Month)


####################### Data normalization #####################
# check: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# Estimating alpha diversity of microbial communities is problematic no matter what you do. 
# My best stab at it is to subsample the libraries with replacement to estimate the species abundance of the real population while standardizing sampling effort.
sort(sample_sums(my_data41)) 
min_lib <- min(sample_sums(my_data41))
min_lib # 3032 for all samples, if we wanna compare

####################################################
# MetagenomeSeq: CSS normalization 
# "Remove the OTUs, of which the sum across all samples is greater than 0.005% of all OTUs based on Bokulich et al. (2013)."

minTotRelAbun = 1e-5
x = taxa_sums(my_data41)
keepTaxa = taxa_names(my_data41)[which((x / sum(x)) > minTotRelAbun)]
prunedSet = prune_taxa(keepTaxa, my_data41)
sort(sample_sums(prunedSet))

test <- phyloseq_to_metagenomeSeq(prunedSet)
p <- cumNormStatFast(test)

## Default value being used.
alphadata <- cumNorm(test, p = p)

# Export normalised and not-log matrix for use 
otu_norm <- MRcounts(alphadata, norm = TRUE, log = FALSE)

data_norm41 <- prunedSet
otu_table(data_norm41) <- otu_table(otu_norm, taxa_are_rows = T)

# CCS somehow fixed the variation in sampling depth. 
# CSS will sometimes decrease the fold-difference in sampling depth but not always. 
max(sample_sums(my_data41))/min(sample_sums(my_data41)) #  41.74472 # raw
max(sample_sums(data_norm41))/min(sample_sums(data_norm41)) # 6.552239 # CSS normalized


############################
# Estimate alpha diversity #
############################

setwd("/home/yo39qed/time-series analysis")
alpha <- estimate_richness(my_data41)
metadata <- sample_data(my_data41)

alpha$Samplename <- row.names(alpha)
metadata$Samplename <- row.names(metadata)
alphadiv <- merge(data.frame(metadata),data.frame(alpha),
                  by = "Samplename", all = TRUE)
row.names(alphadiv) <- alphadiv$Samplename

write.table(alphadiv, "output/H41/alphadiv41_not_normalized.csv", sep=";", dec=",", col.names=NA) # The final metadata file with sequence =3000 /sample for depth for alphadiv estimation

sample_data(my_data41) <- alphadiv
sample_data(data_norm41) <- alphadiv

##### save the raw data for H41 well
save(my_data41, file ="my_data41.phyloseq")

####################################
########  phylogenetic tree ######## 
####################################

# get ref sequences
ref <- refseq(data_norm41)
library(Biostrings)
writeXStringSet(ref, "output/H41/rep_seq_CSS.fasta")

######## Calculate the phylogenetic tree in mothur program ######

# dist.seqs(fasta=rep_seq_CSS.fasta,output=lt, processors=4)
# clearcut(phylip=current)
# quit()

#### integrate the tree file into the phyloseq file
tree <- read.tree("output/H41/rep_seq_CSS.phylip.tre") 
# plot(basinbactree) # this step will take very long time due to too many taxa
phytree <- phy_tree(tree)

data_norm41 <- merge_phyloseq(data_norm41, phytree)
save(data_norm41, file ="output/H41/data_norm41.phyloseq")

#####################################################################
##### Phylogenetic diversity (Faiths PD-index) ######################

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

pd <- estimate_pd(data_norm41)
sample_data(data_norm41)$PD <- pd$PD

save(data_norm41, file ="output/H41/data_norm41.phyloseq")

#################################################################
### ##### Now define the recharge and Discharge period ###### ###
#################################################################

# Based on water level
metadata <- sample_data(data_norm41)

metadata$YearMonth<-paste0(metadata$Year,metadata$Month)
metadata$Recharge <- as.factor(metadata$YearMonth)

metadata$Recharge <- revalue(metadata$Recharge, c("201302"="Recharge", "201303"="Recharge", "201305"="Recharge",
                                                  
                                                  "201307"="Discharge", "201308"="Discharge","201309" ="Discharge", 
                                                  "201310"="Discharge", "201311"="Discharge", "201312"="Discharge",
                                                  "201401"="Discharge", "201402"="Discharge","201403"="Discharge","201404"="Discharge",
                                                  "201406"="Discharge", "201407"="Discharge","201408"="Discharge",
                                                  "201409"="Discharge","201410"="Discharge","201411"="Discharge",
                                                  
                                                  "201501"="Recharge","201502"="Recharge","201503"="Recharge","201504"="Recharge", "201505"="Recharge",
                                                  
                                                  "201506"="Discharge", "201507"="Discharge", "201508"="Discharge",
                                                  "201509"="Discharge","201510"="Discharge","201511"="Discharge",
                                                  
                                                  "201512"="Recharge", "201602"="Recharge","201603"="Recharge","201604"="Recharge","201605"="Recharge",
                                                  
                                                  "201606"="Discharge","201607"="Discharge", "201608"="Discharge",
                                                  "201609"="Discharge","201610"="Discharge","201611"="Discharge",
                                                  "201612"="Discharge","201701"="Discharge",
                                                  
                                                  "201703"="Recharge","201704"="Recharge","201705"="Recharge",
                                                  
                                                  "201706"="Discharge", "201707"="Discharge", "201708"="Discharge", "201710"="Discharge",
                                                  
                                                  "201712"="Recharge","201801"="Recharge", "201802"="Recharge",
                                                  "201803"="Recharge","201804"="Recharge", "201805"="Recharge",
                                                  
                                                  "201807"="Discharge", "201808"="Discharge","201809"="Discharge",
                                                  "201810"="Discharge","201811"="Discharge","201812"="Discharge",
                                                  "201901"="Discharge",
                                                  
                                                  "201902"="Recharge","201903"="Recharge","201904"="Recharge",
                                                  "201905"="Recharge"
                                                  
))
metadata$Period <- as.factor(metadata$YearMonth) 
metadata$Period <- revalue(metadata$Period, c("201302"="P13", "201303"="P13", "201305"="P13",
                                              
                                              "201307"="P15", "201308"="P15","201309" ="P15", 
                                              "201310"="P15", "201311"="P15", "201312"="P15",
                                              "201401"="P15", "201402"="P15","201403"="P15","201404"="P15",
                                              "201406"="P15", "201407"="P15","201408"="P15",
                                              "201409"="P15","201410"="P15","201411"="P15",
                                              
                                              "201501"="P15","201502"="P15","201503"="P15","201504"="P15", "201505"="P15",
                                              
                                              "201506"="P16", "201507"="P16", "201508"="P16",
                                              "201509"="P16","201510"="P16","201511"="P16",
                                              
                                              "201512"="P16", "201602"="P16","201603"="P16","201604"="P16","201605"="P16",
                                              
                                              "201606"="P17","201607"="P17", "201608"="P17",
                                              "201609"="P17","201610"="P17","201611"="P17",
                                              "201612"="P17","201701"="P17",
                                              
                                              "201703"="P17","201704"="P17","201705"="P17",
                                              
                                              "201706"="P18", "201707"="P18", "201708"="P18", "201710"="P18",
                                              
                                              "201712"="P18","201801"="P18", "201802"="P18",
                                              "201803"="P18","201804"="P18", "201805"="P19",
                                              
                                              "201807"="P19", "201808"="P19","201809"="P19",
                                              "201810"="P19","201811"="P19","201812"="P19",
                                              "201901"="P19",
                                              
                                              "201902"="P19","201903"="P19","201904"="P19",
                                              "201905"="P19"
))

metadata$PeriodRecharge<-paste0(metadata$Period,metadata$Recharge) 

#View(data.frame(metadata))
# Manual correction -- > Change PNK102: to Period5Discharge
reP15 <- subset(metadata, Sam_cam!="PNK102") # May 30, 2018, it is not recharge at all
reP16 <- subset(metadata, Sam_cam=="PNK102")
reP16[,c("Recharge", "Period", "PeriodRecharge")] <- c("Discharge",  "P19", "P19Discharge")
metadata <- rbind(reP15, reP16)

sample_data(data_norm41)<-metadata

###### Add qPCR data ####
# Add qpcr data
setwd("/home/yo39qed/time-series analysis/input")
qpcr <- read.csv("qPCR.csv", sep=";", dec = ",",  header=TRUE)

setwd("/home/yo39qed/time-series analysis/output")
qpcr_H41 <- qpcr[,c("Samp_Cam", "H41")]
colnames(qpcr_H41) <-  c("Sam_cam", "qPCR")
summary1 <- data.frame(sample_data(data_norm41))

summary <- merge(data.frame(summary1), data.frame(qpcr_H41),
                 by = "Sam_cam", all=T) 

summary<-summary[!is.na(summary$Samplename),] # remove empty rows with NA values (samples with missing microbial data)
row.names(summary)<-summary$Samplename
sample_data(data_norm41)<-summary

##### save the normalized data for H41 well
# Note: we used the normalized counts for all statistical analysis
save(data_norm41, file ="H41/data_norm41.phyloseq")

