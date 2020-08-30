rm(list = ls())

setwd("/home/yo39qed/time-series analysis/input/metadata/")
#set_colors_Year <-c("red", "brown", "orange", "pink3", "purple", "blue", "black")
set_colors_Year <- get_palette("lancet", 7)

require(zoo)
library(magrittr)
library(gganimate)


############## making metadata table ###########################################
# import all weather, hydrochemistry, temperature and water level data

hydro <- read.csv("hydrochemistry_updated.csv", sep=";", dec = ",",  header=TRUE)
weather <- read.table("weather_parameters_Weberstedt.csv", sep=";", dec = ",", header=TRUE)
tempr <- read.csv("temperature_updated.csv", sep=";", dec = ",", row.names=1, header=TRUE)
waterlevel <- read.csv("water_level_updated.csv", sep=";", dec = ",",  header=TRUE)
# import all sampling date
sampdates <- read.csv("sampling_dates.csv", sep=";", dec = ",",  header=TRUE)
## For specific well
specific_well <- subset(hydro, Well == "H52") # take out all specific well

df_weather <-subset(weather, Date %in% sampdates$Date )  # match the date 

spwell <- subset(specific_well, Date %in% sampdates$Date )  # match the date 

df_tempr<-subset(tempr[,c("Date_Ref", "Date", "H52_Temp")], Date %in% sampdates$Date )  # match the date 

df_waterlevel <-subset(waterlevel[,c("Date_Ref", "Date","H52_depth_to_water_m_bgl", "H52_water_level_m_amsl")], Date %in% sampdates$Date )  # match the date 

ts1 <- merge(data.frame(df_weather, row.names=NULL), data.frame(df_tempr, row.names=NULL),
             by = "Date", all = TRUE)
ts2 <- merge(data.frame(spwell, row.names=NULL), data.frame(df_waterlevel, row.names=NULL), by = "Date", all = TRUE)

ts <- merge(data.frame(ts2, row.names=NULL), data.frame(ts1, row.names=NULL), by = "Date", all = TRUE)

# check if the date after merge of different datasets were consistent
ts[,c("Date_Ref.x.x", "Date_Ref.x.y","Date_Ref.y.x", "Date_Ref.y.y")]

# remove the unnecessary columns
ts <- subset(ts, select=-Date_Ref.x.y)
ts <- subset(ts, select=-Date_Ref.y.y)
ts <- subset(ts, select=-Date_Ref.y.x)

# rename some columns
colnames(ts)[colnames(ts)=="Date_Ref.x.x"] <- "Date_Ref"
colnames(ts)[colnames(ts)=="H52_Temp"] <- "Temp"
colnames(ts)[colnames(ts)=="H52_depth_to_water_m_bgl"] <- "WL_mbgl"
colnames(ts)[colnames(ts)=="H52_water_level_m_amsl"] <- "WL_mamsl"

metadata_H52 <- ts
write.table(metadata_H52, "output/metadata_H52.csv", sep=";", dec=",", col.names=NA)

############# Then we can do necessary editing in excel
# There are missing values ("Sample" and "Sam_Cam", because we don't have hydrochemistry data in there)
# data.frame(sample_sums(my_data))
## Note: there were two spwell PNK69 in the .biom file (H52_PNK69 and H52PNK69_02)
# we need to decide to use the one with bigger library size: H52PNK69_02
### H52_PNK69: 5663 and H52PNK69_02: 11018

# Note: H52PNK66 (H-21) doesn't have hydrochemical data

############################## Plot figures ##################################
##### water temp, water level and precipitation ###
#####################################################
setwd("/home/yo39qed/time-series analysis")

weather$Date_Ref <- format(as.Date(weather$Date, origin="1899-12-30"), "%m/%d/%Y")
weather$Date_Ref <- as.Date(weather$Date_Ref, "%m/%d/%Y")

tempr$Date_Ref <- format(as.Date(tempr$Date, origin="1899-12-30"), "%m/%d/%Y")
tempr$Date_Ref <- as.Date(tempr$Date_Ref, "%m/%d/%Y")

waterlevel$Date_Ref <- format(as.Date(waterlevel$Date, origin="1899-12-30"), "%m/%d/%Y")
waterlevel$Date_Ref <- as.Date(waterlevel$Date_Ref, "%m/%d/%Y")

min=41323 #min(sample_data(data_norm41)$Date)
max=43614 #max(sample_data(data_norm41)$Date)
min
max

# sample_data(data_norm52)$Date_Ref1
weather_H52 <- subset(weather, Date >= min & Date <= max,
                      select=c(Date, Date_Ref, precipitation__sum)) 

tempr_H52 <- subset(tempr, Date >= min & Date <= max,
                    select=c(Date, Date_Ref, H52_Temp)) 


waterlevel_H52 <- subset(waterlevel, Date >= min & Date <= max,
                         select=c(Date, Date_Ref, H52_water_level_m_amsl))

write.table(waterlevel_H52, "output/H52/water_level_H52.csv", sep=";", dec=",", col.names=NA)
######################################################################################################
##########################################      Plot    ##############################################
######################################################################################################

###  Do not do start ## 
# because the precipitation is the same for all wells
# 1) Precipiation

weather_H52[is.na(weather_H52)] <- 0
write.table(weather_H52, "output/H52/weather_H52.csv", sep=";", dec=",", col.names=NA)
# I calculated the accumuated precipitation per year
weather_H52 <- read.csv("output/H43/weather_H43_accumulated.csv", sep=";", dec = ",", row.names=1, header=TRUE)
weather_H52$Date_Ref <- format(as.Date(weather_H52$Date, origin="1899-12-30"), "%m/%d/%Y")
weather_H52$Date_Ref <- as.Date(weather_H52$Date_Ref, "%m/%d/%Y")
weather_H52$Year <- format(as.Date(weather_H52$Date_Ref, format="%d/%m/%Y"),"%Y")
weather_H52$Month <- format(as.Date(weather_H52$Date_Ref, format="%d/%m/%Y"),"%m")
weather_H52$Year <- as.factor(weather_H52$Year)
weather_H52$Month <- as.factor(weather_H52$Month)
yearly_precipitation<-
  aggregate( weather_H52$precipitation__sum, 
             by=list(Category=weather_H52$Year), FUN=sum)

rain_H52 <- ggplot(weather_H52, aes(x=Date_Ref))+ 
# geom_line(aes(y = Yearly_accumulated,color=Year),size=1, show.legend=FALSE)+
 # geom_bar(aes(y = Yearly_accumulated), color="grey", stat="identity", width=0.5, show.legend=FALSE) + 
  geom_area(aes(y=Yearly_accumulated), fill="grey", position = 'stack', show.legend=FALSE)+
  geom_line(aes(y = precipitation__sum*10,color=Year),size=0.5)+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0.5))+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) +
  xlab("Date")+
  ylab("Accumulated precipitation [mm]")+
 scale_y_continuous(sec.axis = sec_axis(~./10, name = "Daily precipitation [mm]"))+ # second axix
  scale_color_manual(values=set_colors_Year)+
  scale_fill_manual(values=set_colors_Year)+
 # xlim(as.Date(c('18/2/2013', '29/5/2019'), format="%d/%m/%Y") ) +
  scale_x_date(limits = as.Date(c('18/2/2013', '29/5/2019'), format="%d/%m/%Y"),
               date_breaks = "2 months", date_minor_breaks = "1 month", 
               labels = date_format("%b %Y") )+
  annotate("label", x = as.Date(c("2013-07-01", "2014-07-01", "2015-07-01", 
                                 "2016-07-01", "2017-07-01", "2018-07-01", 
                                 "2019-04-01")), y = 880, size=5,
           label = round(yearly_precipitation$x,0), color=set_colors_Year)+
  theme( axis.line.y.left = element_line(color = "grey21"), 
         axis.text.y.left = element_text(color = "grey21"),
         axis.title.y.left = element_text(color = "grey21"),
           axis.ticks.y.left = element_line(color = "grey21"))+
  theme(legend.position="none")
rain_H52  

ggsave("output/H52/rain_time_H52.pdf", width = 10, height = 4)



#####################################################################################################
######################              water level and temperature              ########################
#####################################################################################################
n <- 22.7

min(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE)
max(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE)
max(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE)-min(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE)
min(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE)/min(tempr_H52$H52_Temp, na.rm=TRUE)
max(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE)/max(tempr_H52$H52_Temp, na.rm=TRUE)


wl_t_H52 <- ggplot(waterlevel_H52
                   [!is.na(waterlevel_H52$H52_water_level_m_amsl),]
                   , aes(x=Date_Ref, y=H52_water_level_m_amsl))+ 
  geom_line(size=1, color="blue")+
  geom_line(data=tempr_H52
           # [!is.na(tempr_H52$H52_Temp),]
            , aes(x=Date_Ref, y=H52_Temp*n),
            size=1, color="red")+
  theme_pubr()+
  geom_hline(yintercept=mean(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE), linetype="dashed", 
               color = "blue", size=1) +
  geom_hline(yintercept=mean(tempr_H52$H52_Temp, na.rm=TRUE)*n, linetype="dashed", 
             color = "red", size=1) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0.5))+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) +
  xlab("Date")+
  ylab("Water level [MAMSL]")+
  scale_y_continuous(limits = c(215, 238),
                     #breaks = c(220, 230, 240, 250), 
                     sec.axis = sec_axis(~./n, 
                                         breaks= c(9.4,9.5, 9.6,9.7, 9.8),
                                         # breaks= pretty_breaks(),
                                         name = "Water temperature [°C]"))+ # second axix
  scale_x_date(limits = as.Date(c('18/2/2013', '29/5/2019'), format="%d/%m/%Y"),
               date_breaks = "2 months", date_minor_breaks = "1 month", 
               labels = date_format("%b %Y") )+
  annotate("label", x = as.Date("2013-03-01"), y = mean(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE), size=5,
           label = round(mean(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE),0), color="blue")+
  annotate("label", x = as.Date("2013-03-01"), y = mean(tempr_H52$H52_Temp, na.rm=TRUE)*n, size=5,
           label = round(mean(tempr_H52$H52_Temp, na.rm=TRUE),1), color="red")+
  theme( axis.line.y.right = element_line(color = "red"), 
         axis.text.y.right = element_text(color = "red"),
         axis.title.y.right = element_text(color = "red"),
         axis.ticks.y.right = element_line(color = "red"))+
theme( axis.line.y.left = element_line(color = "blue"), 
       axis.text.y.left = element_text(color = "blue"),
       axis.title.y.left = element_text(color = "blue"),
       axis.ticks.y.left = element_line(color = "blue"))

wl_t_H52 


ggsave("output/H52/waterlevel_temperature_H52.pdf", width = 10, height = 4)



#### Water level only  ############
wl_H52 <- ggplot(waterlevel_H52
                 [!is.na(waterlevel_H52$H52_water_level_m_amsl),]
                 , aes(x=Date_Ref, y=H52_water_level_m_amsl))+ 
  geom_line(size=1, color="black")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0.5))+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) +
  xlab("Date")+
  geom_hline(yintercept=mean(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE), linetype="dashed", 
             color = "black", size=1) +
  annotate("label", x = as.Date("2013-03-01"), y = mean(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE), size=5,
           label = round(mean(waterlevel_H52$H52_water_level_m_amsl, na.rm=TRUE),0), color="black")+
  ylab("Water level [MAMSL]")+
  scale_x_date(limits = as.Date(c('18/2/2013', '29/5/2019'), format="%d/%m/%Y"),
               date_breaks = "2 months", date_minor_breaks = "1 month", 
               labels = date_format("%b %Y") )

wl_H52

ggsave("output/H52/waterlevel_time_H52.pdf", width = 9.45, height = 4)

##################################################################################################################
# Plot water temperature

t_H52 <- ggplot(data=tempr_H52
                #[!is.na(tempr_H52$H52_Temp),]
                , aes(x=Date_Ref, y=H52_Temp)
)+
  geom_line(size=1, color="black" )+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0.5))+
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14),legend.title=element_text(size=14, face="bold")) +
  xlab("Date")+
  geom_hline(yintercept=round(mean(tempr_H52$H52_Temp, na.rm=TRUE),1), linetype="dashed", 
             color = "black", size=1) +
  annotate("label", x = as.Date("2013-03-01"), y = mean(tempr_H52$H52_Temp, na.rm=TRUE), size=5,
           label = round(mean(tempr_H52$H52_Temp, na.rm=TRUE),1), color="black")+
  ylab("Water temperature [°C]")+
  scale_x_date(limits = as.Date(c('18/2/2013', '29/5/2019'), format="%d/%m/%Y"),
               date_breaks = "2 months", date_minor_breaks = "1 month", 
               labels = date_format("%b %Y") )

t_H52 

ggsave("output/H52/temp_time_H52.pdf", width = 9.45, height = 4)
