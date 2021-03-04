library(tidyverse)
library(lubridate)
library(adehabitatLT)
library(bayesmove)
library(move)

options(scipen = 999)

#STEP 1. SEGMENT BY SEASON USING NET-SQUARED DISPLACEMENT####

#Load in raw data & filter out pts without enough satellite
dat <- read.csv("PinPoint 2217 2019-07-10 11-21-47.csv") %>% 
  dplyr::filter(Status=="Valid") %>% 
  mutate(doy=yday(ymd(RTC.date)))

#Wrangle for as.traj function. Note filtered to just one point per day
locs <- dat %>% 
  dplyr::select(RTC.date, RTC.time, Latitude, Longitude, Altitude)  %>% 
  mutate(date=as.POSIXct(ymd(RTC.date)),
         time=hms(RTC.time),
         ID=2217,
         burst=1) %>% 
  dplyr::filter(hour(time)==15)

#Use as.traj function to get NSD
traj <- as.ltraj(xy=locs[,c("Longitude", "Latitude")],
                 id=locs$ID,
                 date=locs$date,
                 burst=locs$burst,
                 proj4string = CRS("+proj=longlat +datum=WGS84"))
#Look at it
traj
head(traj[[1]])

#Join back to dataframe and segment
trajdf <- as.data.frame(traj[[1]]) %>% 
  mutate(date=as_date(date),
         doy=yday(date),
         index=row_number(),
         season=case_when(R2n<1 ~ "breed",
                          R2n>8355 ~ "winter",
                          R2n>1 & R2n<8355 & year(date)==2018 ~ "fallmig",
                          R2n>1 & R2n<8355 & year(date)==2019 ~ "springmig"))

#Check segmentation
table(trajdf$season)

#Visualize
ggplot(trajdf) +
  geom_line(aes(x=index, y=R2n)) +
  geom_point(aes(x=index, y=R2n, colour=season))
#Some funnyness in winter classifications

#Segment per location for wintering grounds instead
trajdf <- as.data.frame(traj[[1]]) %>% 
  mutate(date=as_date(date),
         doy=yday(date),
         index=row_number(),
         season=case_when(R2n<1 ~ "breed",
                          Longitude > -53.33 ~ "winter",
                          R2n>1 & Longitude < -53.33 & year(date)==2018 ~ "fallmig",
                          R2n>1 & Longitude < -53.3 & year(date)==2019 ~ "springmig"))

#Join back to full datset
datseg <- 

#Write out
write.csv(datseg, "PinPoint2217_NSDsegmented.csv", row.names = FALSE)

#STEP 2. SEGMENT STOPOVER FROM MIGRATION FOR DAILY POINTS####

datmig <- datseg %>% 
  filter(hour==15,
         season %in% c("fallmig", "springmig"))

#Visualize elevation data to determine whether it has the potential to inform segmentation
ggplot(datseg) +
  geom_point(aes(x=Index, y=R2n, colour=log(Altitude))) +
  scale_colour_viridis_c()

ggplot(datseg) +
  geom_violin(aes(x=season, y=Altitude)) +
  geom_jitter(aes(x=season, y=Altitude, colour=hour)) +
  facet_wrap(~season, scales="free") +
  scale_colour_viridis_c()

#Yes, especially for winter. Explore bayesmove instead of HMMs
#Try observation level segmentation instead of using foiegras and then segmenting

 #2a. Fall migration----

#2b. Spring migration----

#STEP 3. SEGMENT FORAGING FROM ROOSTING FROM MIGRATION FOR BURST POINTS####
#Do this separately for each season?