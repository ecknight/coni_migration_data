#title: Behavioural segmentation of Common Nighthawk migration data
#author: Elly C. Knight
#created: April 6, 2021

#CONI are a fly-and-forage migrator, meaning that they do not make formal stopovers, but move every day of migration. This continuous motion combined with the somewhat tortuous migration path means that classification of migration and stopover behaviours via traditional methods (e.g., moveHMM, Bayesmove) is not effective. I've instead using a simple clustering analysis of geographic location to identify locations where migration is slower and this individual is obviously spending more time. I've retained only cluster wtih at least 3 points (i.e., >= 72 hours), as per Warnock et al. 2010 J Av. Biol 41:621-626.

#PREAMBLE####

library(tidyverse) #for data wrangling
library(lubridate) #for date manipulation
library(sf) #for gis
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

datsun <- read.csv("00_PinPoint2217_SeasonClassified.csv")

#STEP 2. SEGMENT STOPOVER FROM MIGRATION - MEANSHIFT####

#2a. Fall----
#2ai. Wrangle----
datfall <- datsun %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datsun) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="fallmig", hour==15)

#2aii. Mean shift classification----
mat.fall1 <- matrix(datfall$X)
mat.fall2 <- matrix(datfall$Y)
mat.fall <- cbind(mat.fall1, mat.fall2)

shiftfall <- meanShift(mat.fall,
                       algorithm="KDTREE",
                       bandwidth=c(1,1))

clustfall <- datfall %>% 
  mutate(cluster = shiftfall[[1]]) %>% 
  group_by(cluster) %>% 
  mutate(count=n()) %>% 
  ungroup() %>% 
  mutate(stopover = ifelse(count > 2, 1, 0))

table(clustfall$cluster)

ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(data=clustfall, aes(x=Longitude, y=Latitude)) +
  geom_point(data=clustfall, aes(x=Longitude, y=Latitude, fill=factor(stopover)), pch=21, size=3, alpha=0.7) +
  labs(x = "", y = "") +
  xlim(c(-170, -30)) +
  theme_bw() +
  scale_fill_manual(values=c("orange", "blue"))

ggsave("MeanShift_Map_Fall.jpeg", width=15, height=22, unit="in")

#2b. Spring----
#2bi. Wrangle----
datspring <- datsun %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datsun) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="springmig", hour==15)

#2bii. Mean shift classification----
mat.spring1 <- matrix(datspring$X)
mat.spring2 <- matrix(datspring$Y)
mat.spring <- cbind(mat.spring1, mat.spring2)

shiftspring <- meanShift(mat.spring,
                         algorithm="KDTREE",
                         bandwidth=c(1,1))

clustspring <- datspring %>% 
  mutate(cluster = shiftspring[[1]]) %>% 
  group_by(cluster) %>% 
  mutate(count=n()) %>% 
  ungroup() %>% 
  mutate(stopover = ifelse(count > 2, 1, 0))

table(clustspring$cluster)

ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(data=clustspring, aes(x=Longitude, y=Latitude)) +
  geom_point(data=clustspring, aes(x=Longitude, y=Latitude, fill=factor(stopover)), pch=21, size=3, alpha=0.7) +
  labs(x = "", y = "") +
  xlim(c(-170, -30)) +
  theme_bw() +
  scale_fill_manual(values=c("orange", "blue"))

ggsave("MeanShift_Map_Spring.jpeg", width=15, height=22, unit="in")

#2c. Put together----
datstop <- rbind(clustfall, clustspring) %>% 
  dplyr::select(Latitude, Longitude, stopover) %>% 
  full_join(datsun) %>% 
  mutate(stopover=ifelse(is.na(stopover), 0, stopover))

write.csv(datstop, "01_PinPoint2217_StopoverClassified.csv", row.names = FALSE)
