rm(list = ls())

################################################
### FEAST ### Note: not working for proportional data
################################################


#saveRDS()
library(phyloseq)
require(gridExtra)
library(vegan)
library(ggplot2)
library(plyr)
library(FSA)

vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

set_colors_Period <-c("P13"= "brown","P15"="red", "P16"="deeppink", "P17"="orange",
                      "P18"="purple", "P19"="blue", "P14"="darkgreen")

####################################################################################
load("/home/yo39qed/time-series analysis/output/H41/data_norm41.phyloseq")
setwd("/home/yo39qed/time-series analysis/FEAST")

# Preparing for metadata and OTU data
# "The sink is recharge samples per period, the source are discharge samples per period
period_option <- c("P15","P16","P17", "P18", "P19")

############################################################################################################
## For loop to calculate the recharge effect (unknown source that contribute to the microbiomes) per phase
############################################################################################################
for (i in 1:length(period_option)){
selected_period <- period_option[i]
data_norm41_select <- subset_samples(data_norm41, Period==selected_period)

metadata <- data.frame(Env=sample_data(data_norm41_select)$Sam_cam, SourceSink1=sample_data(data_norm41_select)$Recharge, SourceSink=sample_data(data_norm41_select)$Recharge) # T3

metadata$SourceSink <- "Source"
metadata$SourceSink[metadata$SourceSink1=="Recharge"] <- "Sink"
row.names(metadata) <- sample_names(data_norm41_select) # T3
metadata <- metadata[, c("Env", "SourceSink")]
otus <- vegan_otu(data_norm41_select)# T3

###############################################
## Now run FEAST_example_Multiple_sink.R
getwd()
print("Change directory path")
dir_path = paste("/home/yo39qed/time-series analysis/FEAST/FEAST-master/")
setwd(paste0(dir_path, "FEAST_src"))
source("src.R")

#Set the arguments of your data
EM_iterations = 1000 #default value
##if you use different sources for each sink, different_sources_flag = 1, otherwise = 0
different_sources_flag = 0

# Extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# Double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}


if(different_sources_flag == 0){
  
  metadata$id[metadata$SourceSink == 'Source'] = NA
  metadata$id[metadata$SourceSink == 'Sink'] = c(1:length(which(metadata$SourceSink == 'Sink')))
}


envs <- metadata$Env
Ids <- na.omit(unique(metadata$id))
Proportions_est <- list()


for(it in 1:length(Ids)){
  
  
  # Extract the source environments and source/sink indices
  if(different_sources_flag == 1){
    
    train.ix <- which(metadata$SourceSink=='Source' & metadata$id == Ids[it])
    test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])

  }
  
  else{
    
    train.ix <- which(metadata$SourceSink=='Source')
    test.ix <- which(metadata$SourceSink=='Sink' & metadata$id == Ids[it])
  }
  
  num_sources <- length(train.ix)
  COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user
  
  # Define sources and sinks
  
  sources <- as.matrix(rarefy(otus[train.ix,], COVERAGE))
  sinks <- as.matrix(rarefy(t(as.matrix(otus[test.ix,])), COVERAGE))
  
  
  print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
  print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
  print(paste("The sink is:", envs[test.ix]))
  
  # Estimate source proportions for each sink
  set.seed(123)
  FEAST_output<-FEAST(source=sources, sinks = t(sinks), env = envs[train.ix], em_itr = EM_iterations, COVERAGE = COVERAGE)
  Proportions_est[[it]] <- FEAST_output$data_prop[,1]
  
  
  names(Proportions_est[[it]]) <- c(as.character(envs[train.ix]), "unknown")
  
  if(length(Proportions_est[[it]]) < num_sources +1){
    
    tmp = Proportions_est[[it]]
    Proportions_est[[it]][num_sources] = NA
    Proportions_est[[it]][num_sources+1] = tmp[num_sources]
  }
  
  print("Source mixing proportions")
  print(Proportions_est[[it]])
  

}

print(Proportions_est)


## Unlist the list data
# unlist(Proportions_est)
prop_merged <- data.frame(Source = names(unlist(Proportions_est)), Proportions_est=unlist(Proportions_est))
prop_merged$Period <-selected_period
setwd("/home/yo39qed/time-series analysis/FEAST")
write.table(prop_merged, paste0("H41/feast_recharge_effect_", selected_period, ".csv"), sep=";", dec=",", col.names = NA)

}

#########################################################################################################
#########################################################################################################
#### Read and combine tables
setwd("/home/yo39qed/time-series analysis/FEAST")
t1 <- read.csv("H41/feast_recharge_effect_P15.csv", sep=";", dec = ",",  header=TRUE)
t2 <- read.csv("H41/feast_recharge_effect_P16.csv", sep=";", dec = ",",  header=TRUE)
t3 <- read.csv("H41/feast_recharge_effect_P17.csv", sep=";", dec = ",",  header=TRUE)
t4 <- read.csv("H41/feast_recharge_effect_P18.csv", sep=";", dec = ",",  header=TRUE)
t5 <- read.csv("H41/feast_recharge_effect_P19.csv", sep=";", dec = ",",  header=TRUE)

t <- do.call("rbind", list(t1, t2, t3, t4, t5))
t$Proportions_est<-t$Proportions_est*100 # change to 100%
prop_merged_unknown <-subset(t, Source=="unknown")

################################## Mean and SD ######################################
res <- aggregate(Proportions_est ~ Period, data=prop_merged_unknown, function(x) c(mean = mean(x), sd = sd(x), se=sd(x)/sqrt(length(x)))) 
res
# write.table(res, "H41/feast_unknown_source_recharge_H41.csv", sep=";", dec=",", col.names = NA)
mean(prop_merged_unknown$Proportions_est)
sd(prop_merged_unknown$Proportions_est)
# Dunn test for multiple comparisons 
#Zar (2010) states that the Dunn test is appropriate for groups 
# with unequal numbers of observations
#install.packages("FSA")

PT<-dunnTest(Proportions_est~Period,prop_merged_unknown, method="bh")
PT<-PT$res

library(rcompanion) #To have R convert this table to a compact letter display for us
cldList(comparison = PT$Comparison,
        p.value    = PT$P.unadj, # or PT$P.unadj
        threshold  = 0.05)

########################################################################################################
prop_merged_unknown$Well <- "H41"
# Add one empty row for P14

test <- prop_merged_unknown[1,]
test$Period <- "P14"
test$Proportions_est <- NA
prop_merged_unknown <-rbind(test, prop_merged_unknown)

prop_merged_unknown$Period <- factor(prop_merged_unknown$Period, levels = c("P14", "P15", "P16", "P17", "P18", "P19"))

p_recharge_feast_H41 <- ggboxplot(prop_merged_unknown, x="Period", y = "Proportions_est", fill = "Period",
                        palette = set_colors_Period,
                        add = "mean" 
                        #shape = "Aquifer",
)+ 
 # stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  ylim(0,80) +
  theme_bw()+
  facet_grid(Well~., scales='free')+
#  ylim(0,80)+
  ylab("Unknown source contribution (%)")+
  stat_compare_means(label.x = 1.5, label.y = 3, size=6.5) +
  # facet_wrap(.~Period,scales = "fixed", nrow = 1)+
  theme(strip.text = element_text(size=14))  +
  geom_hline(yintercept=mean(prop_merged_unknown$Proportions_est, na.rm=TRUE), linetype="dashed", 
             color = "black", size=1) +
  ggtitle("Unknown source contributing to H41 recharge communities")+
  annotate("label", x = 0.8, y = mean(prop_merged_unknown$Proportions_est, na.rm=TRUE), size=6.5,
           label = round(mean(prop_merged_unknown$Proportions_est, na.rm=TRUE),1), color="black")+
  theme(strip.text = element_text(size=20))  +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20), legend.text=element_text(size=20),
        legend.title=element_text(size=20, face="bold")) 

p_recharge_feast_H41

setwd("/home/yo39qed/time-series analysis/output")
ggsave('H41/feast_Recharge_H41.pdf', width = 8, height = 4)
ggsave('H41/feast_Recharge_H41.jpg', width = 8, height = 4)

################################### Combine all wells ###################################
prop_merged_unknow_H41 <- prop_merged_unknown
prop_merged_unknow_H41$Well <- "H41"


