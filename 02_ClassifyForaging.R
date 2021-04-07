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

#Bayesmove by point----
for(i in 1:length(doys.stat)){
  
  doy.i <- doys.stat$doy[i]
  
  #Filter to fall migration points, prepare for bayesmove data prep
  dat.i <- datburst %>% 
    st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
    st_transform(crs=3857) %>% 
    st_coordinates() %>% 
    cbind(datburst) %>% 
    mutate(id=2217) %>% 
    dplyr::filter(doy==doy.i) %>% 
    dplyr::select(id, date, X, Y, Altitude, sethrs) 
  
  #Calulate step lengths, turning angles, and time steps
  steps.i <- prep_data(dat = dat.i, coord.names = c("X","Y"), id = "id")
  
  #Visualize
  # View distributions of data streams
  hist(steps.i$step, breaks=20)
  hist(steps.i$angle, breaks=20)
  hist(steps.i$Altitude, breaks=20)
  hist(steps.i$sethrs, breaks=20)
  
  # Create list from data frame
  round.i <- round_track_time(dat = steps.i, id = "id", int = 3600, tol = 120, time.zone = "UTC", units="secs") %>% 
    mutate(dt = ifelse(row_number()==nrow(steps.i), 3600, dt))
  list.i<- df_to_list(dat = round.i, ind = "id")
  filter.i <- filter_time(dat=list.i, int = 3600) %>% 
    bind_rows
  
  # Define bin number and limits for turning angles
  angle.bin.lims<- seq(from=-pi, to=pi, by=pi/4)  #8 bins
  
  # Define bin number and limits for step lengths
  dist.bin.lims<- quantile(filter.i$step, c(0,0.10,0.25,0.5,0.75,0.9,1), na.rm=T)
  
  # Define bin number and limits for altitude
  alt.bin.lims<- quantile(filter.i$Altitude, c(0,0.10,0.25,0.5,0.75,0.9,1), na.rm=T)
  
#  discrete.i<- discrete_move_var(filter.i, lims = list(dist.bin.lims, angle.bin.lims, alt.bin.lims), varIn = c("step", "angle", "Altitude"), varOut = c("SL","TA", "AL"))
  
  discrete.i<- discrete_move_var(filter.i, lims = list(dist.bin.lims, alt.bin.lims), varIn = c("step",  "Altitude"), varOut = c("SL", "AL"))
  
  # Only retain id and discretized step length (SL), turning angle (TA), altitude
#  input.i <- subset(discrete.i, select = c(id, SL, TA, AL))
  input.i <- subset(discrete.i, select = c(id, SL, AL))

  #2aii. Segmentation----
  
  inputlist.i <-  df_to_list(dat = input.i, ind = "id")
  
  set.seed(1)
  
  # Define hyperparameter for prior distribution
  alpha<- 1
  # Set number of iterations for the Gibbs sampler
  ngibbs<- 10000
  # Set the number of bins used to discretize each data stream
#  nbins<- c(7,8,7)
  nbins<- c(7,7)
  
  future::plan(multisession)  #run all MCMC chains in parallel
  #refer to future::plan() for more details
  dat.res.i <- segment_behavior(data = inputlist.i, ngibbs = ngibbs, nbins = nbins,
                              alpha = alpha)
  
  # Trace-plots for the number of breakpoints per ID
  traceplot(data = dat.res.i$nbrks, ngibbs = ngibbs, type = "nbrks")
  
  # Trace-plots for the log marginal likelihood (LML) per ID
  traceplot(data = dat.res.i$LML, ngibbs = ngibbs, type = "LML")
  
  # Determine MAP for selecting breakpoints
  MAP.est.i<- get_MAP(dat = dat.res.i$LML, nburn = 5000)
  MAP.est.i
  
  brkpts.i<- get_breakpts(dat = dat.res.i$brkpts, MAP.est = MAP.est.i)
  
  # How many breakpoints estimated per ID?
  apply(brkpts.i[,-1], 1, function(x) length(purrr::discard(x, is.na)))
  
  # Plot breakpoints over the data
  discretelist.i<- df_to_list(dat = discrete.i, ind = "id")
  
#  plot_breakpoints(data = discretelist.i, as_date = FALSE, var_names = c("step","angle", "Altitude"), var_labels = c("Step Length (km)", "Turning Angle (rad)", "Altitude (m)"), brkpts = brkpts.i)
  plot_breakpoints(data = discretelist.i, as_date = FALSE, var_names = c("step","Altitude"), var_labels = c("Step Length (km)", "Altitude (m)"), brkpts = brkpts.i)
  
  # Assign track segments to all observations by ID
  segfall<- assign_tseg(dat = discretelistfall, brkpts = brkpts)
  
  head(tracks.seg)
  
  ggsave(paste0("Figures/BurstClassification/BurstClassification_BayesMove_",i,".jpeg"), width=8, height=6)
  
  
}

#MoveHMM----

doys <- datstop %>% 
  dplyr::filter(hour==4) %>% 
  dplyr::select(doy, season) %>% 
  unique()

datburst <- datstop %>% 
  dplyr::filter(hour!=15,
                season=="winter",
                doy %in% doys$doy)

doys.list <- unique(datburst$doy)

for(i in 1:length(doys.list)){
  
  doy.i <- doys.list[i]
  
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
  mu0 <- c(1,10) #(state1,state2) # mean step lengths for state 1 (stopover) and state 2 (migration)
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