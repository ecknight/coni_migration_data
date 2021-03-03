library(tidyverse)
library(lubridate)
library(adehabitatLT)

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

#Segment as per NSD
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
  geom_point(aes(x=index, y=R2n))

#Join back to full datset
datseg <- dat %>% 
  left_join(trajdf %>% 
              dplyr::select(doy, R2n, season)) %>% 
  mutate(season=ifelse(is.na(season), "breed", season))

#Write out
write.csv(datseg, "/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/PinPoint2217_NSDsegmented.csv", row.names = FALSE)

#STEP 2. SEGMENT STOPOVER FROM MIGRATION FOR DAILY POINTS####


#STEP 3. SEGMENT FORAGING FROM ROOSTING FROM MIGRATION FOR BURST POINTS####
#Do this separately for each season?