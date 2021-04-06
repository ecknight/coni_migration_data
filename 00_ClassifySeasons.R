#title: Behavioural segmentation of Common Nighthawk migration data
#author: Elly C. Knight
#created: March 4, 2021

#remotes::install_github("joshcullen/bayesmove")

#PREAMBLE####

library(tidyverse) #for data wrangling
library(lubridate) #for date manipulation
library(adehabitatLT) #for NSD calculation
library(bayesmove) #for behavioural segmentation
library(sf) #for gis
library(furrr) #for segmentation in bayesmove
library(purrr) #for segmentation in bayesmove
library(moveHMM) #for behavioural segmentation
library(suncalc) #for time relative to sunset calculations
library(moveNT) #for network analysis
library(meanShiftR) #for density cluster analysis

options(scipen = 999)


#Get background map data
whemi <- map_data("world", region=c("Canada", 
                                    "USA", 
                                    "Mexico",
                                    "Guatemala", 
                                    "Belize", 
                                    "El Salvador",
                                    "Honduras", 
                                    "Nicaragua", 
                                    "Costa Rica",
                                    "Panama", 
                                    "Jamaica", 
                                    "Cuba", 
                                    "The Bahamas",
                                    "Haiti", 
                                    "Dominican Republic", 
                                    "Antigua and Barbuda",
                                    "Dominica", 
                                    "Saint Lucia", 
                                    "Saint Vincent and the Grenadines", 
                                    "Barbados",
                                    "Grenada",
                                    "Trinidad and Tobago",
                                    "Colombia",
                                    "Venezuela",
                                    "Guyana",
                                    "Suriname",
                                    "Ecuador",
                                    "Peru",
                                    "Brazil",
                                    "Bolivia",
                                    "Paraguay",
                                    "Chile",
                                    "Argentina",
                                    "Uruguay")) %>% 
  dplyr::filter(!group%in%c(258:264))

#STEP 1. WRANGLING####

#1a. Segment by season----

#Load in raw data & filter out pts without enough satellite
dat <- read.csv("PinPoint 2217 2019-07-10 11-21-47.csv") %>% 
  dplyr::filter(Status=="Valid")

#Wrangle for as.traj function. Note filtered to just one point per day
locs <- dat %>% 
  dplyr::select(RTC.date, RTC.time, Latitude, Longitude, Altitude)  %>% 
  mutate(date=as.POSIXct(ymd_hms(paste0(RTC.date, RTC.time))),
         time=hms(RTC.time),
         id=2217)

#Use as.traj function to get NSD
traj <- as.ltraj(xy=locs[,c("Longitude", "Latitude")],
                 id=locs$id,
                 date=locs$date,
                 proj4string = CRS("+proj=longlat +datum=WGS84"))
#Look at it
traj
head(traj[[1]])

#Segment per location for wintering grounds instead
datseg <- dat %>% 
  left_join(as.data.frame(traj[[1]]) %>%
              mutate(doy=yday(date)) %>% 
              dplyr::rename(Longitude=x, Latitude=y)) %>% 
  mutate(season=case_when(R2n<1 ~ "breed",
                          Longitude > -53.33 ~ "winter",
                          R2n>1 & Longitude < -53.33 & year(date)==2018 ~ "fallmig",
                          R2n>1 & Longitude < -53.3 & year(date)==2019 ~ "springmig"),
         date=as.POSIXct(ymd_hms(paste0(RTC.date, RTC.time))),
         time=hms(RTC.time),
         hour=hour(time))

#Check segmentation
table(datseg$season)
summary(datseg)

#Visualize
ggplot(datseg) +
  geom_line(aes(x=Index, y=R2n)) +
  geom_point(aes(x=Index, y=R2n, colour=season))

ggplot(subset(datseg, season %in% c("springmig", "fallmig"))) +
  geom_path(aes(x=Longitude, y=Latitude)) +
  geom_point(aes(x=Longitude, y=Latitude, colour=doy), size=4) +
  scale_colour_viridis_c() +
  facet_wrap(~season)

#1b. Calculate time relative to sunset----
datsun <- getSunlightTimes(data=datseg %>% 
                             rename(lat=Latitude, lon=Longitude) %>% 
                             mutate(date=as.Date(date))) %>% 
  dplyr::select(sunrise, sunset) %>% 
  cbind(datseg) %>% 
  mutate(sethrs = as.numeric((sunset - date)/3600),
         sethrs = ifelse(sethrs > 5, sethrs-24, sethrs),
         risehrs = as.numeric((sunrise - date)/3600),
         risehrs = ifelse(risehrs > 12, risehrs-24, risehrs),)

#Write out
write.csv(datsun, "PinPoint2217_NSDsegmented.csv", row.names = FALSE)