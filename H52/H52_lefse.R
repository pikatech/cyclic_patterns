rm(list = ls())
setwd("/home/yo39qed/time-series analysis/output")

#saveRDS()
library(phyloseq)
library(reshape2)
library(readr)
require(gridExtra)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(data.table)
# library(remotes)
# remotes::install_github("ying14/yingtools2")
# library(yingtools2) # not working, so I have to use the original functions
set_colors_Year <- get_palette("lancet", 7)

set_colors_Period <-c("P0"= "brown","P1"="red", "P2"="deeppink", "P3"="orange",
                      "P4"="purple", "P5"="blue")

set_colors_recharge <- c("Recharge"="orange", "Discharge"="grey")
set_shapes_recharge <- c("Recharge"=15, "Discharge"=1)

set_colors_PeriodRecharge <- c("P0Recharge" ="brown","P1Recharge" = "red", "P1Discharge" ="red", 
                               "P2Recharge" ="deeppink", "P2Discharge"="deeppink", "P3Recharge"="orange",    
                               "P3Discharge"="orange", "P4Recharge" ="purple",  "P4Discharge"= "purple",
                               "P5Recharge"= "blue", 
                               #"P6Discharge"="blue", 
                               "P5Discharge"= "blue")


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


set_colors_phylum <- c("Acidobacteria" ="cyan1",    "Actinobacteria" = "brown",   "Bacteroidetes"  ="orange1",    "BHI80_139" = "lightgreen",        
                       "Chlamydiae" ="bisque",       "Chloroflexi" ="lightpink",       "Cyanobacteria"  ="#00FFBF",    
                       "Elusimicrobia"   = "burlywood",    "Entotheonellaeota" = "pink2",  "Firmicutes" ="gray",        "Hydrogenedentes"="lightgreen",    "Kiritimatiellaeota"="chocolate2", 
                       "Nitrospinae" = "mediumpurple4",       "Nitrospirae" ="blue",     
                       "Omnitrophicaeota" ="darkorange2",   "Patescibacteria" ="green",   "Planctomycetes"="firebrick2",     "Proteobacteria"  = "purple4",   "Rokubacteria" ="orange1",    
                       "Spirochaetes" ="yellow",      "Verrucomicrobia"="red2"
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

################################################################################################################
#############################  lefse ###########################################################################

get.samp <- function(phy,stats=FALSE,measures=c("Observed","InvSimpson","Shannon")) {
  requireNamespace(c("phyloseq","tibble"),quietly=TRUE)
  
  if (is.null(sample_data(phy,FALSE))) {
    #if no sample_data, return single data frame with sample column
    sdata <- data.frame(sample=sample_names(phy),stringsAsFactors=FALSE)
  } else {
    if ("sample" %in% phyloseq::sample_variables(phy)) {stop("YTError: phyloseq sample_data already contains the reserved variable name \"sample\"")}
    sdata <- sample_data(phy) %>% data.frame(stringsAsFactors=FALSE) %>% tibble::rownames_to_column("sample")
  }
  if (stats) {
    dup.names <- intersect(c("nseqs",measures),names(sdata))
    if (length(dup.names)>0) {
      sdata <- sdata[,setdiff(names(sdata),dup.names)]
      warning("YTWarning: Following variables are duplicated. Deleting old values from phyloseq: ",paste(dup.names,collapse=", "))
    }
    sdata$nseqs <- phyloseq::sample_sums(phy)
    sdata <- cbind(sdata,estimate_richness(phy,measures=measures))
  }
  return(sdata)
}

get.otu.melt = function(phy,filter.zero=TRUE,sample_data=TRUE) {
  requireNamespace(c("phyloseq","data.table"),quietly=TRUE)
  # supports "naked" otu_table as `phy` input.
  otutab = as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    otutab <- t(otutab)
  }
  otudt = data.table(otutab, keep.rownames = TRUE)
  data.table::setnames(otudt, "rn", "otu")
  # Enforce character otu key
  # note that .datatable.aware = TRUE needs to be set for this to work well.
  otudt[, otuchar:=as.character(otu)]
  otudt[, otu := NULL]
  data.table::setnames(otudt, "otuchar", "otu")
  # Melt count table
  mdt = data.table::melt.data.table(otudt, id.vars = "otu", variable.name = "sample",value.name = "numseqs")
  if (filter.zero) {
    # Remove zeroes, NAs
    mdt <- mdt[numseqs > 0][!is.na(numseqs)]
  } else {
    mdt <- mdt[!is.na(numseqs)]
  }
  # Calculate relative abundance
  mdt[, pctseqs := numseqs / sum(numseqs), by = sample]
  if(!is.null(tax_table(phy, errorIfNULL=FALSE))) {
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(phy, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "otu")
    # Enforce character otu key
    taxdt[, otuchar := as.character(otu)]
    taxdt[, otu := NULL]
    setnames(taxdt, "otuchar", "otu")
    # Join with tax table
    setkey(taxdt, "otu")
    setkey(mdt, "otu")
    mdt <- taxdt[mdt]
  }
  if (sample_data & !is.null(sample_data(phy, errorIfNULL = FALSE))) {
    # If there is a sample_data, join with it.
    sampledt = data.table(as(sample_data(phy, errorIfNULL = TRUE), "data.frame"),keep.rownames=TRUE)
    setnames(sampledt, "rn", "sample")
    # Enforce character sample key
    sampledt[, samplechar := as.character(sample)]
    sampledt[, sample := NULL]
    setnames(sampledt, "samplechar", "sample")
    # Join with tax table
    setkey(sampledt, "sample")
    setkey(mdt, "sample")
    mdt <- sampledt[mdt]
  }
  return(tbl_df(mdt))
  # return(mdt)
}

lefse2 <- function (phy, class, subclass = NA, subject = NA, anova.alpha = 0.05, 
                    wilcoxon.alpha = 0.05, lda.cutoff = 2, wilcoxon.within.subclass = FALSE, 
                    one.against.one = FALSE, mult.test.correction = 0, make.lefse.plots = FALSE, 
                    by_otus = FALSE, levels = rank_names(phy)) { 
  keepvars <- c(class, subclass, subject, "sample")
  keepvars <- unique(keepvars[!is.na(keepvars)])
  samp <- get.samp(phy)[, keepvars]
  if (by_otus) {
    otu <- get.otu.melt(phy, sample_data = FALSE)
    otu.levels <- otu %>% 
      mutate(taxon = otu) %>% 
      group_by(sample,taxon) %>% 
      summarize(pctseqs = sum(pctseqs)) %>% 
      mutate(taxon = gsub(" ", "_", taxon))
  }else{
    otu <- get.otu.melt(phy, sample_data = FALSE)
    otu.list <- lapply(1:length(levels), function(i) {
      lvls <- levels[1:i]
      lvl <- levels[i]
      otu.level <- otu
      otu.level$taxon <- do.call(paste, c(lapply(lvls, function(l) otu[[l]]), sep = "|"))
      otu.level$rank <- lvl
      otu.level2 <- otu.level %>% 
        group_by(sample, taxon,  rank) %>% 
        summarize(pctseqs = sum(pctseqs)) %>% 
        ungroup()
      return(otu.level2)
    })
    otu.levels <- bind_rows(otu.list) %>% 
      mutate(taxon = gsub(" ","_", taxon))
  }
  otu.tbl <- otu.levels %>% 
    dcast(sample ~ taxon, value.var = "pctseqs",fill = 0) %>% 
    left_join(samp, by = "sample") %>% 
    select_(.dots = c(keepvars,lazyeval::interp(~everything())))
  if (is.na(subject) | subject != "sample") {
    otu.tbl <- otu.tbl %>% select(-sample)
  }
  tbl <- otu.tbl %>% t()
  write.table(tbl, "lefse.txt", quote = FALSE, sep = "\t", 
              col.names = FALSE)
  opt.class <- paste("-c", which(keepvars %in% class))
  opt.subclass <- ifelse(is.na(subclass), "", paste("-s", which(keepvars %in% subclass)))
  opt.subject <- ifelse(is.na(subject), "", paste("-u", which(keepvars %in% subject)))
  format.command <- paste("python format_input.py lefse.txt lefse.in",  opt.class, opt.subclass, opt.subject, "-o 1000000")
  system(format.command)
  lefse.command <- paste("python run_lefse.py lefse.in lefse.res", 
                         "-a", anova.alpha, "-w", wilcoxon.alpha, "-l", lda.cutoff, 
                         "-e", as.numeric(wilcoxon.within.subclass), "-y", as.numeric(one.against.one), 
                         "-s", mult.test.correction)
  system(lefse.command)
  print("Wrote lefse.res")
  lefse.out <- read.table("lefse.res", header = FALSE, sep = "\t") %>% 
    rename(taxon = V1, log.max.pct = V2, direction = V3,lda = V4, p.value = V5)
  if (make.lefse.plots) {
    system("python plot_res.py lefse.res lefse_lda.png")
    print("Wrote lefse_lda.png")
    system("python plot_cladogram.py lefse.res lefse_clado.pdf --format pdf")
    print("Wrote lefse_clado.pdf")
  }
  return(lefse.out)
}


############################## Lefse to test the effect of recharge on each taxonomic group ##############

data_norm52_prop <- data_norm52 %>%
  subset_samples(Period !="P13")%>%
 # transform_sample_counts( function(x) x/sum(x)*100)%>%
  filter_taxa( function(x) sum(x) > 0, TRUE)

# lefse_data52 <- filter_taxa(data_norm52_prop, function(x) mean(x) > 1e-3, TRUE)
lefse_data52 <- data_norm52_prop
#move to directory where you downloaded the lefse zip, in my case it was in downloads:
setwd("/home/yo39qed/lefse/code")

#run the new function on the phyloseq object, specifying class and subclass variable names that must be in your sample_data slot of phyloseq object
lefse.tbl <- lefse2(lefse_data52, class="Recharge" ,subclass = "Period"
)
# with subclass (Period): P13 is removed
# without subclass (Period): P13 is kept

############################################################################################################################
#### Now do it in galaxy ###
# I chose the result with LDA > 3, edit the .lefse_internal_res file (remove the insig. taxa in Notepad ++)
# 1) Search->mark..-> input "charge"-> tick "bookmark line" -> Mark All 
# 2) Search->Bookmark -> Remove Unmarked lines
# upload the edited .lefse_internal_res file into galaxy, specify the format type as "lefse_internal_res"

# 2, 5, 4, 5 # LDAover3 ----- H41
# 2, 4, 4, 4 # LDAover2, not working well, too many features - H41 ----- till order level
# default setting: 2, 6, 4, 6 # LDAover2 ------ till genus level


# Note: change colors of the groups in .svg file (notepad+++, replace color code)!!!
# orange: #FFA500, grey: #BEBEBE
# green: #008000, red: #ff0000
# If necessary, compress the cladogram figures (.svg) online, otherwise inkscape works extremely slowly: https://www.svgminify.com/
# Edit the cladogram figures (e.g. legends and font sizes) in inkscape
############################################################################################################################

LDA_H52 <- read.delim("~/time-series analysis/input/lefse/Galaxy92-[B)_LDA_Effect_Size_(LEfSe)_on_data_91]_H52_LDAover2.lefse_internal_res", header=F)

names(LDA_H52) <- c("Feature", "log_of_the_highest_class_average", "class", "LDA_effect_size", "p-value")

#make a lefse style barplot.. here i filtered for 2 consistencies, big LDA values:

LDA_H52 %>%
  filter(!is.na(LDA_effect_size),
         LDA_effect_size>3,
         class %in% c("Recharge","Discharge")) %>%
  mutate(LDA_effect_size=ifelse(class=="Recharge",-LDA_effect_size,LDA_effect_size)) %>%
  ggplot() +
  geom_bar(aes(x=reorder(Feature,LDA_effect_size),y=LDA_effect_size,fill=class),color="black",stat="identity") +
  theme_pubr()+
  scale_fill_manual(values=set_colors_recharge) +
  geom_hline(yintercept=2, linetype="dashed", 
             color = "black", size=0.5) +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "black", size=0.5) +
  geom_hline(yintercept=-2, linetype="dashed", 
             color = "black", size=0.5) +
  geom_hline(yintercept=-1, linetype="dashed", 
             color = "black", size=0.5) +
  geom_hline(yintercept=3, linetype="dashed", 
             color = "black", size=0.5) +
  geom_hline(yintercept=-3, linetype="dashed", 
             color = "black", size=0.5) +
  geom_hline(yintercept=-4, linetype="dashed", 
             color = "black", size=0.5) +
  geom_hline(yintercept=4, linetype="dashed", 
             color = "black", size=0.5) +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), 
        legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold"))+
  scale_y_continuous( breaks = c(-4,-3,-2,-1,0,1,2,3,4)) + 
  ylab("LDA effect size")+
  coord_flip()

#### Extract the significant LDA values
LDA_H52 %>%
  filter(!is.na(LDA_effect_size),
         LDA_effect_size>2,
         class %in% c("Recharge"
                    # ,
                    # "Discharge"
         )) -> filter_3 # with absolute LDA score > 2
dim(filter_3) # Recharge:33, Discharge:40

############### Or we can separate at different taxonomic levels and plot ##### 
LDA_H52 %>%
  filter(!is.na(LDA_effect_size),
         LDA_effect_size>2,
         class %in% c("Recharge","Discharge")) -> tet


library(stringr)
tet$Kingdom <- str_split_fixed(tet$Feature, "\\.", 6)[ ,1]
tet$Phylum <- str_split_fixed(tet$Feature, "\\.", 6)[ ,2]
tet$Class <- str_split_fixed(tet$Feature, "\\.", 6)[ ,3]
tet$Order <- str_split_fixed(tet$Feature, "\\.", 6)[ ,4]
tet$Family <- str_split_fixed(tet$Feature, "\\.", 6)[ ,5]
tet$Genus <- str_split_fixed(tet$Feature, "\\.", 6)[ ,6]

#write.table(tet, "H52/lefse_result_for_bar_plot.csv", sep=";", dec=",", col.names=NA)

tet[tet==""] <- NA #Recoding Values to Missing

tet_phylum  <- tet[is.na(tet$Class),]
tet_class  <- tet[is.na(tet$Order),]
tet_order  <- tet[is.na(tet$Family),]
tet_family  <- tet[is.na(tet$Genus),]
tet_genus  <- tet[!is.na(tet$Genus),]

### Yeah!! Now we can plot the results at different levels!!
# 1) at genus level

tet_genus %>%
  mutate(LDA_effect_size=ifelse(class=="Discharge",-LDA_effect_size,LDA_effect_size)) -> tet2 # reversed LDA values in 2 directions

# Change the names a bit, because some genus share the same name but belong to different families
attach(tet2)
#detach(tet2)
lb <- str_split_fixed(tet2[order(LDA_effect_size),]$Feature , "\\.", 6)[ ,6] # [, 6] is genus level, [, 5]: family, etc 

#########################################################################################
tet2 %>%
  ggplot()   +
  theme_pubr()+
  scale_fill_manual(values=set_colors) +
  geom_bar(aes(x=reorder(Feature,LDA_effect_size),y=LDA_effect_size,fill=Class),color="black",stat="identity")+ 
  geom_hline(yintercept=2, linetype="dashed", 
             color = "orange", size=1) +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "orange", size=1) +
  geom_hline(yintercept=-2, linetype="dashed", 
             color = "grey", size=1) +
  geom_hline(yintercept=-1, linetype="dashed", 
             color = "grey", size=1) +
  geom_hline(yintercept=3, linetype="dashed", 
             color = "orange", size=1) +
  geom_hline(yintercept=-3, linetype="dashed", 
             color = "grey", size=1) +
  geom_hline(yintercept=-4, linetype="dashed", 
             color = "grey", size=1) +
  geom_hline(yintercept=4, linetype="dashed", 
             color = "orange", size=1) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=2) +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20), 
        legend.text=element_text(size=20),legend.title=element_text(size=20, face="bold"))+
  scale_x_discrete(labels=lb) +
  scale_y_continuous( breaks = c(-4,-3,-2,-1,0,1,2,3,4)) + 
  theme(legend.position = "right") +  
  #scale_x_discrete(expand = c(0, 2))+
  xlab("Genus")+
  ylab("LDA effect size")+
  # ggtitle("LEfSe")+
  # theme(plot.margin=unit(c(8,1,1,1), "lines"))+ # top, right, bottom,left
  coord_flip() + annotate("label", y = -2, x = 14, size=10,
                          label ="Recharge-depressed", fill="grey", color="white")+
  annotate("label", y = 2, x = 9, size=10,
           label ="Recharge-favored", fill="orange", color="white") -> p2

p2

setwd("/home/yo39qed/time-series analysis/output")
ggsave("H52/lefse_with_subclass_genus_H52.pdf", width = 18, height = 12)
ggsave("H52/lefse_with_subclass_genus_H52.png", width = 18, height = 12)

