rm(list = ls())

setwd("/home/yo39qed/time-series analysis/output")

#saveRDS()
library(phyloseq)
require(gridExtra)
library(ggplot2)
library(ade4)
library(mgcv)


# Load data

load("/home/yo39qed/time-series analysis/output/H52/data_norm52.phyloseq")

#########################################################################################################################################
##################################### Linking to date ########################################################################################

dist_H52 <- phyloseq::distance(data_norm52, "bray")
m_H52 <- as.matrix(dist_H52)
m_H52[,which(m_H52 == max(dist_H52), arr.ind = T)]
sub_otu <- otu_table(data_norm52) # to get the otu table 
sub_map <- sample_data(data_norm52) # to get the metadata table
sub_map$day <- sub_map$Date - (min(sub_map$Date)-1) # 52323 is the date (numeric) of PNK33, so the sampling date of PNK33 is Day 1

Day <- sub_map[,"day"] 

variable.names(sub_map)
mean(dist_H52)  
sd(dist_H52) 
median(dist_H52) 


#### Plot #####

mantel_r52 <- mantel.rtest(dist_H52, dist(Day), nrepet = 9999)
mantel_r52

df52_1 <- data.frame(cbind(dist_H52, dist(Day)))
names(df52_1) <- c("Bray", "Day")
df52_1$Well <- "H52"

pt52_bray <- ggplot(df52_1, aes(x=Day, y=Bray)) + 
  geom_point(color="#A2D9F7")+
  geom_smooth()+
  theme_classic()+
 # stat_cor(label.y = 0.95, size=5 ) +         # Add correlation coefficient
  ylab("Ecological distance (Bray-Curtis)") + 
  xlab("Day changes")+
  #geom_hline(yintercept = mean(df52_1$Bray), linetype =2, lwd=1, color = "orange")+
  annotate("text", y = 0.28, x = 1200, size=5,fill="transparent",
           label = paste0("min_dist = ", round(min(df52_1$Bray),2),"; max_dist = ", round(max(df52_1$Bray),2),"; median_dist = ", round(median(df52_1$Bray),2)), color="blue")+
  #  stat_regline_equation(label.y = 1.0, size=5)+
 # facet_grid(Well~., scales="free")+
  theme(strip.text.x = element_text(size = 16))+
  theme(strip.background =element_rect(fill="#A2D9F7"))+
  annotate("text", y = 0.34, x = 1550, size=5,
           label ="Mantel test: simulated p < 0.001", fill="transparent", color="black")+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14))+
  scale_x_continuous(breaks=seq(0, max(df52_1$Day), 365))+
  geom_vline(xintercept = seq(0, max(df52_1$Day), 365), linetype="dashed", 
             color = "grey", size=0.5)+
  theme(axis.title=element_blank(),
        axis.text.x=element_blank()#,
        # axis.ticks=element_blank() 
  )

pt52_bray

ggsave("H52/bray/correlation_water_day_bray_H52.pdf", width = 7, height = 5.1)
ggsave("H52/bray/correlation_water_day_bray_H52.jpg", width = 7, height = 5.1)


## Remove all points, leaving only the curve
t52_bray_no_point <- ggplot(df52_1, aes(x=Day, y=Bray)) + 
  #  geom_point(color="#A2D9F7")+
  geom_smooth()+
  theme_classic()+
  # stat_cor(label.y = 0.95, size=5 ) +         # Add correlation coefficient
  # ylab("Ecological distance") + 
  # xlab("Day changes")+
  #geom_hline(yintercept = mean(df52_1$Bray), linetype =2, lwd=1, color = "orange")+
  #  annotate("label", y = 0.28, x = 1200, size=5,fill="transparent",
  #       label = paste0("min_dist = ", round(min(df52_1$Bray),2),"; max_dist = ", round(max(df52_1$Bray),2),"; median_dist = ", round(median(df52_1$Bray),2)), color="blue")+
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
  scale_x_continuous(breaks=seq(0, max(df41_1$Day), 365),  labels = seq(0, 6))+
  geom_vline(xintercept = seq(0, max(df52_1$Day), 365), linetype="dashed", 
               color = "grey", size=0.5)

t52_bray_no_point



t52_bray_no_point1 <- ggplot(df52_1, aes(x=Day, y=Bray)) + 
  #  geom_point(color="#A2D9F7")+
  geom_smooth(size=2)+
  theme_classic()+
  # stat_cor(label.y = 0.95, size=5 ) +         # Add correlation coefficient
  # ylab("Ecological distance") + 
  # xlab("Day changes")+
  #geom_hline(yintercept = mean(df52_1$Bray), linetype =2, lwd=1, color = "orange")+
  #  annotate("label", y = 0.28, x = 1200, size=5,fill="transparent",
  #       label = paste0("min_dist = ", round(min(df52_1$Bray),2),"; max_dist = ", round(max(df52_1$Bray),2),"; median_dist = ", round(median(df52_1$Bray),2)), color="blue")+
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
  scale_x_continuous(breaks=seq(0, max(df41_1$Day), 365),  labels = seq(0, 6))+
  geom_vline(xintercept = seq(0, max(df52_1$Day), 365), linetype="dashed", 
             color = "grey", size=0.5)

t52_bray_no_point1

ggsave("H52/bray/correlation_water_day_bray_no_point_H52_1.pdf", width = 4, height = 4)
ggsave("H52/bray/correlation_water_day_bray_no_point_H52_1.jpg", width = 4, height = 4)

### fitting gam
# https://people.maths.bris.ac.uk/~sw15190/mgcv/check-select.pdf

fit <- gam(Bray ~ s(Day, bs = "cs"), data = df52_1) # family: Gassian
summary(fit) # needed

fit$coefficients 
fit$sp
# We can use the GCV scores a bit like AIC, smaller values indicated better fitting models.
# The deviance explained is a bit like R2 for models where sums of squares 
# doesn't make much sense as a measure of discrepancy between the observations and the fitted values.
gam.check(fit)
# Basis dimension (k): the maximum allowable degrees of freedom for smooth terms 
par(mfrow = c(2,2))
gam.check(fit)
plot(fit)

par(mfrow = c(1,1))
plot(Bray ~ Day, data =df52_1, col = "red",  pch = 16)
curve(predict(fit, newdata = data.frame(Day = x)), add = TRUE, 
      from = min(df52_1$Day), to = max(df52_1$Day), n = 1e3, lwd = 2)

## Residual (Tukey Anscombe) plot:
plot(residuals(fit) ~ fitted(fit))
abline(h = 0, col = "gray")

