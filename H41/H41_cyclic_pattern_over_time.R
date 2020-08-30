rm(list = ls())

setwd("/home/yo39qed/time-series analysis/output")

#saveRDS()
library(phyloseq)
require(gridExtra)
library(ggplot2)
library(ade4)
library(mgcv)


# Load data

load("/home/yo39qed/time-series analysis/output/H41/data_norm41.phyloseq")

#########################################################################################################################################
##################################### Linking to date ########################################################################################

dist_H41 <- phyloseq::distance(data_norm41, "bray")
m_H41 <- as.matrix(dist_H41)
m_H41[,which(m_H41 == max(dist_H41), arr.ind = T)]
sub_otu <- otu_table(data_norm41) # to get the otu table 
sub_map <- sample_data(data_norm41) # to get the metadata table
sub_map$day <- sub_map$Date - (min(sub_map$Date)-1) # 41323 is the date (numeric) of PNK33, so the sampling date of PNK33 is Day 1

Day <- sub_map[,"day"] 

variable.names(sub_map)
mean(dist_H41)  
sd(dist_H41) 
median(dist_H41) 


#### Plot #####

mantel_r41 <- mantel.rtest(dist_H41, dist(Day), nrepet = 9999)
mantel_r41

df41_1 <- data.frame(cbind(dist_H41, dist(Day)))
names(df41_1) <- c("Bray", "Day")
df41_1$Well <- "H41"

pt41_bray <- ggplot(df41_1, aes(x=Day, y=Bray)) + 
  geom_point(color="#A2D9F7")+
  geom_smooth()+
  theme_classic()+
 # stat_cor(label.y = 0.95, size=5 ) +         # Add correlation coefficient
  ylab("Ecological distance (Bray-Curtis)") + 
  xlab("Day changes")+
  #geom_hline(yintercept = mean(df41_1$Bray), linetype =2, lwd=1, color = "orange")+
  annotate("text", y = 0.28, x = 1200, size=5,fill="transparent",
           label = paste0("min_dist = ", round(min(df41_1$Bray),2),"; max_dist = ", round(max(df41_1$Bray),2),"; median_dist = ", round(median(df41_1$Bray),2)), color="blue")+
  #  stat_regline_equation(label.y = 1.0, size=5)+
 # facet_grid(Well~., scales="free")+
  theme(strip.text.x = element_text(size = 16))+
  theme(strip.background =element_rect(fill="#A2D9F7"))+
  annotate("text", y = 0.34, x = 1550, size=5,
           label ="Mantel test: simulated p < 0.001", fill="transparent", color="black")+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14))+
  scale_x_continuous(breaks=seq(0, max(df41_1$Day), 365))+
  geom_vline(xintercept = seq(0, max(df41_1$Day), 365), linetype="dashed", 
             color = "grey", size=0.5)+
  theme(axis.title=element_blank(),
        axis.text.x=element_blank()#,
        # axis.ticks=element_blank() 
  )

pt41_bray

ggsave("H41/bray/correlation_water_day_bray_H41.pdf", width = 7, height = 5.1)
ggsave("H41/bray/correlation_water_day_bray_H41.jpg", width = 7, height = 5.1)


## Remove all points, leaving only the curve
t41_bray_no_point <- ggplot(df41_1, aes(x=Day, y=Bray)) + 
  #  geom_point(color="#A2D9F7")+
  geom_smooth()+
  theme_classic()+
  # stat_cor(label.y = 0.95, size=5 ) +         # Add correlation coefficient
  # ylab("Ecological distance") + 
  # xlab("Day changes")+
  #geom_hline(yintercept = mean(df41_1$Bray), linetype =2, lwd=1, color = "orange")+
  #  annotate("label", y = 0.28, x = 1200, size=5,fill="transparent",
  #       label = paste0("min_dist = ", round(min(df41_1$Bray),2),"; max_dist = ", round(max(df41_1$Bray),2),"; median_dist = ", round(median(df41_1$Bray),2)), color="blue")+
  #  stat_regline_equation(label.y = 1.0, size=5)+
  facet_grid(~Well)+
  theme(strip.text = element_text(size = 16))+
  theme(strip.background =element_rect(fill="#A2D9F7"))+
  # annotate("label", y = 0.34, x = 1550, size=5,
  #          label ="Mantel test: simulated P < 0.001", fill="transparent", color="black")+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14))+
  theme(axis.title=element_blank()#,
        # axis.text=element_blank(),
        # axis.ticks=element_blank() 
  )+
  scale_x_continuous(breaks=seq(0, max(df41_1$Day), 365))+
  geom_vline(xintercept = seq(0, max(df41_1$Day), 365), linetype="dashed", 
               color = "grey", size=0.5)

t41_bray_no_point



t41_bray_no_point1 <- ggplot(df41_1, aes(x=Day, y=Bray)) + 
  #  geom_point(color="#A2D9F7")+
  geom_smooth(size=2)+
  theme_classic()+
  # stat_cor(label.y = 0.95, size=5 ) +         # Add correlation coefficient
  # ylab("Ecological distance") + 
  # xlab("Day changes")+
  #geom_hline(yintercept = mean(df41_1$Bray), linetype =2, lwd=1, color = "orange")+
  #  annotate("label", y = 0.28, x = 1200, size=5,fill="transparent",
  #       label = paste0("min_dist = ", round(min(df41_1$Bray),2),"; max_dist = ", round(max(df41_1$Bray),2),"; median_dist = ", round(median(df41_1$Bray),2)), color="blue")+
  #  stat_regline_equation(label.y = 1.0, size=5)+
  #facet_grid(~Well)+
  theme(strip.text = element_text(size = 16))+
  theme(strip.background =element_rect(fill="#A2D9F7"))+
  # annotate("label", y = 0.34, x = 1550, size=5,
  #          label ="Mantel test: simulated P < 0.001", fill="transparent", color="black")+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14))+
  theme(axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank() 
  )+
  scale_x_continuous(breaks=seq(0, max(df41_1$Day), 365))+
  geom_vline(xintercept = seq(0, max(df41_1$Day), 365), linetype="dashed", 
             color = "grey", size=0.5)

t41_bray_no_point1

ggsave("H41/bray/correlation_water_day_bray_no_point_H41_1.pdf", width = 4, height = 4)
ggsave("H41/bray/correlation_water_day_bray_no_point_H41_1.jpg", width = 4, height = 4)

### fitting gam
# https://people.maths.bris.ac.uk/~sw15190/mgcv/check-select.pdf

fit <- gam(Bray ~ s(Day, bs = "cs"), data = df41_1) # family: Gassian

summary(fit) # needed

fit$coefficients 
# We can use the GCV scores a bit like AIC, smaller values indicated better fitting models.
# The deviance explained is a bit like R2 for models where sums of squares 
# doesn't make much sense as a measure of discrepancy between the observations and the fitted values.
gam.check(fit)
plot(fit)

plot(Bray ~ Day, data =df41_1, col = "red",  pch = 16)
curve(predict(fit, newdata = data.frame(Day = x)), add = TRUE, 
      from = min(df41_1$Day), to = max(df41_1$Day), n = 1e3, lwd = 2)

## Residual (Tukey Anscombe) plot:
plot(residuals(fit) ~ fitted(fit))
abline(h = 0, col = "gray")

