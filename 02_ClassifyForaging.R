#title: Behavioural segmentation of Common Nighthawk migration data
#author: Elly C. Knight
#created: April 6, 2021

#This script is only for days when tags took one point per hour over night (i.e, every 10 days)

#PREAMBLE####

library(tidyverse) #for data wrangling
library(lubridate) #for date manipulation
library(sf) #for gis

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

datstop <- read.csv("01_PinPoint2217_StopoverClassified.csv")

#STEP 3. SEGMENT FORAGING FROM ROOSTING FROM MIGRATION FOR BURST POINTS####

#Visualize elevation data to determine whether it has the potential to inform segmentation
ggplot(datstop) +
  geom_violin(aes(x=season, y=Altitude)) +
  geom_jitter(aes(x=season, y=Altitude, colour=hour)) +
  facet_wrap(~season, scales="free") +
  scale_colour_viridis_c()

#Yes, especially for winter foraging vs roosting

#I think need to do each day separately because are 10 days apart
doys <- datstop %>% 
  dplyr::filter(hour==4) %>% 
  dplyr::select(doy, season) %>% 
  unique()

#3a. Stationary seasons----
doys.stat <- doys %>% 
  dplyr::filter(season %in% c("winter", "breed"))

datburst <- datstop %>% 
  dplyr::filter(hour!=15,
                doy %in% doys.stat$doy)

for(i in 1:nrow(doys.stat)){
  
  doy.i <- doys.stat$doy[i]
  
  #2ai. Wrangling----
  #Filter to fall migration points, prepare for bayesmove data prep
  dat.i <- datburst %>% 
    st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
    st_transform(crs=3857) %>% 
    st_coordinates() %>% 
    cbind(datburst) %>% 
    mutate(id=2217) %>% 
    dplyr::filter(doy==doy.i) %>% 
    dplyr::select(id, date, X, Y, Altitude, sethrs) %>% 
    mutate(X=X/1000,
           Y=Y/1000)
  
  prep.i <- prepData(dat.i, type="UTM", coordNames = c("X", "Y"))
  
  plot(prep.i, compact=T)
  
  #2aii. Model----
  mu0 <- c(0.5,2) #(state1,state2) # mean step lengths for state 1 (roosting) and state 2 (foraging)
  sigma0 <- c(1,300) #(state1,state2)
  stepPar0 <- c(mu0,sigma0)
  angleMean0 <- c(pi,0) # angle mean
  kappa0 <- c(1,1) # angle concentration
  anglePar0 <- c(angleMean0,kappa0)
  formula <- ~Altitude + sethrs
  
  m <- fitHMM(data=prep.i,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0, formula=formula)
  m
  CI(m)
  plot(m, plotCI=TRUE)
  
  #2aiii. Classify----
  states <- viterbi(m)
  sp <- stateProbs(m)
  
  plotStates(m,animals="Animal1")
  
  
}