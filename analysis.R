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

#STEP 2. SEGMENT STOPOVER FROM MIGRATION - BAYESMOVE - POINT LEVEL####

#https://joshcullen.github.io/bayesmove/articles/Cluster-observations.html

#2a. Fall migration----

#2ai. Wrangling----
#Filter to fall migration points, prepare for bayesmove data prep
datfall <- datsun %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datsun) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="fallmig", hour==15) %>% 
  dplyr::select(id, date, X, Y, Altitude) 

#Calculate step lengths, turning angles, and time steps
stepsfall <- prep_data(dat = datfall, coord.names = c("X","Y"), id = "id")

#Visualize
# View distributions of data streams
#Step length
q <- quantile(stepsfall$step, seq(0,1,0.1), na.rm=T) %>% 
  data.frame() 
q$quantile <- paste0("percent", rownames(q))
colnames(q) <- c("value", "quantile")

#Look at jenks breaks
j <- getJenksBreaks(stepsfall$step, 1) #not useful

ggplot(stepsfall) +
  geom_histogram(aes(x=step)) +
  geom_vline(data=q, aes(xintercept=value, colour=quantile)) +
  geom_vline(aes(xintercept=j[1]), linetype="dashed")

ggplot(stepsfall) +
  geom_histogram(aes(x=log(step))) +
  geom_vline(data=q, aes(xintercept=log(value), colour=quantile)) +
  geom_vline(aes(xintercept=log(j[1])), linetype="dashed")

#Step angle
ggplot(stepsfall) +
  geom_histogram(aes(x=angle))

# Create list from data frame
roundfall <- round_track_time(dat = stepsfall, id = "id", int = 86400, tol = 600, time.zone = "UTC", units="secs")
listfall<- df_to_list(dat = roundfall, ind = "id")
filterfall <- filter_time(dat=listfall, int = 86400) %>% 
  bind_rows

# Define bin number and limits for turning angles
angle.bin.lims<- seq(from=-pi, to=pi, by=pi/4)  #8 bins

# Define bin number and limits for step lengths
dist.bin.lims<- quantile(filterfall$step, c(0,0.1,0.25,0.50,0.75,1), na.rm=T) #6 bins

# Define bin number and limits for altitude
alt.bin.lims<- quantile(filterfall$Altitude, c(0,0.25,0.50,0.75,0.90,1), na.rm=T)

discretefall<- discrete_move_var(filterfall, lims = list(dist.bin.lims, angle.bin.lims), varIn = c("step","angle"), varOut = c("SL","TA"))
#discretefall<- discrete_move_var(filterfall, lims = list(dist.bin.lims), varIn = c("step"), varOut = c("SL"))
#discretefall<- discrete_move_var(filterfall, lims = list(dist.bin.lims, alt.bin.lims), varIn = c("step", "Altitude"), varOut = c("SL", "AL"))

# Only retain id and discretized step length (SL), turning angle (TA), altitude
inputfall <- subset(discretefall, select = c(id, SL, TA))
#inputfall <- subset(discretefall, select = c(id, SL))
#inputfall <- subset(discretefall, select = c(id, SL, AL))

#2aii. Classification----

set.seed(1)

# Define model params
alpha=0.1  #prior
ngibbs=5000  #number of Gibbs sampler iterations
nburn=ngibbs/2  #number of burn-in iterations
nmaxclust=2  #number of maximum possible states (clusters) present

# Run model
dat.res<- bayesmove::cluster_obs(dat=inputfall, alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust,
                      nburn=nburn)

# Inspect traceplot of log-likelihood
plot(dat.res$loglikel, type = "l")

# Determine the MAP estimate of the posterior
MAP.iter<- get_MAP_internal(dat = dat.res$loglikel, nburn = nburn)

# Determine the likely number of behavioral states (i.e., what are the fewest states that account for >= 90% of all observations?)
theta<- dat.res$theta[MAP.iter,]
names(theta)<- 1:length(theta)
theta<- sort(theta, decreasing = TRUE)
theta %>% cumsum()

# Store cluster order for plotting and behavioral state extraction
ord<- as.numeric(names(theta))

# Extract bin estimates for each possible state from the `phi` matrix of the model results
behav.res<- get_behav_hist(dat = dat.res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle"), ord = ord,
                           MAP.iter = MAP.iter)

# Plot state-dependent distributions 
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(3), "grey35","grey35",
                               "grey35","grey35"), guide = FALSE) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ var, scales = "free_x")

# Extract MAP estimate of behavioral states
z<- factor(dat.res$z[[MAP.iter]])
levels(z)<- ord

# Relabel factor levels
statesfall <- discretefall %>% 
  st_as_sf(coords=c("x", "y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  cbind(discretefall) %>% 
  mutate(state = case_when(as.character(z) == 1 ~ "Migration",
                           as.character(z) == 2 ~ "Stopover"))

#Plot output
ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(data = statesfall, aes(x=X, y=Y), color="gray60", size=0.25) +
  geom_point(data = statesfall, aes(X, Y, fill=state), size=3, pch=21, alpha=0.7) +
  labs(x = "", y = "") +
  xlim(c(-170, -30)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14)))

ggsave("FallMigration_DailyClassification_BayesMove.jpeg", width=8, height=6)

#STEP 2. SEGMENT STOPOVER FROM MIGRATION - BAYESMOVE - WITH BREAKPOINTS####

#https://joshcullen.github.io/bayesmove/articles/Segment-the-tracks.html
#https://joshcullen.github.io/bayesmove/articles/Prepare-data-for-analysis.html
#https://joshcullen.github.io/bayesmove/articles/Cluster-track-segments.html

#2a. Fall migration----

#Run step 2ai above first

#2aii. Segmentation----

inputlistfall <-  df_to_list(dat = inputfall, ind = "id")

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1
# Set number of iterations for the Gibbs sampler
ngibbs<- 10000
# Set the number of bins used to discretize each data stream
nbins<- c(6,8)
#nbins<- c(6)

future::plan(multisession)  #run all MCMC chains in parallel
#refer to future::plan() for more details
dat.res <- segment_behavior(data = inputlistfall, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")

# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")

# Determine MAP for selecting breakpoints
MAP.est<- get_MAP(dat = dat.res$LML, nburn = 5000)
MAP.est

brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))

# Plot breakpoints over the data
discretelistfall<- df_to_list(dat = discretefall, ind = "id")

plot_breakpoints(data = discretelistfall, as_date = FALSE, var_names = c("step","angle"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)"), brkpts = brkpts)

#This is not looking super useful. Seems like step length might be the only informative variable. Should consider removing turning angle.

# Assign track segments to all observations by ID
segfall<- assign_tseg(dat = discretelistfall, brkpts = brkpts)

head(tracks.seg)


#STEP 2. SEGMENT STOPOVER FROM MIGRATION - MOVEHMM####

#https://cran.r-project.org/web/packages/moveHMM/vignettes/moveHMM-guide.pdf

#2a. Fall migration----
#2ai. Wrangling----
#Filter to fall migration points, prepare for bayesmove data prep
datfall <- datsun %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datsun) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="fallmig", hour==15) %>% 
  dplyr::select(id, date, X, Y, Altitude) %>% 
  mutate(X=X/1000,
         Y=Y/1000)

prepfall <- prepData(datfall, type="UTM", coordNames = c("X", "Y"))

plot(prepfall, compact=T)

#2aii. Model----
mu0 <- c(10,100) #(state1,state2) # mean step lengths for state 1 (stopover) and state 2 (migration)
sigma0 <- c(1,300) #(state1,state2)
stepPar0 <- c(mu0,sigma0)
angleMean0 <- c(pi,0) # angle mean
kappa0 <- c(1,1) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

m <- fitHMM(data=prepfall,nbStates=2,stepPar0=stepPar0,
            anglePar0=anglePar0)
m
CI(m)
plot(m, plotCI=TRUE)

#2aiii. Classify----
states <- viterbi(m)
sp <- stateProbs(m)

plotStates(m,animals="Animal1")

#STEP 2. SEGMENT STOPOVER FROM MIGRATION - MOVENT####
library(devtools)
install_github("BastilleRousseau/moveNT")

#2ai. Wrangle----
datfall <- datsun %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datsun) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="fallmig", hour==15) %>% 
  dplyr::select(id, X, Y, date)

#make an ltraj
traj <- as.ltraj(xy=datfall[,c("X", "Y")],
                 id=datfall$id,
                 date=datfall$date,
                 proj4string = CRS("+init=epsg:3857"))
adj <- traj2adj(traj, res=100000)
dim(adj[[1]])
plot(adj[[2]])
plot(adj[[3]])

stk <- adj2stack(adj, grph=F)
plot(stk)
graphmet(stk)
cv(val(stk, 4))

clust2<-clustnet(stk, id=3, nclust=2, grph=F)
summary(clust2[[1]])
plot(clust2[[2]])

#Nope

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
  mutate(stopover = ifelse(count > 1, 1, 0))

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

ggplot(clustfall) +
  geom_line(aes(x=doy, y=lag(dist))) +
  geom_point(aes(x=doy, y=lag(dist), colour=factor(stopover)))

#2b. Spring----

#2ai. Wrangle----
datspring <- datsun %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datsun) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="springmig", hour==15)

#2aii. Mean shift classification----
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
  mutate(stopover = ifelse(count > 1, 1, 0))

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

ggplot(clustspring) +
  geom_line(aes(x=doy, y=lag(dist))) +
  geom_point(aes(x=doy, y=lag(dist), colour=factor(stopover)))


#STEP 3. SEGMENT FORAGING FROM ROOSTING FROM MIGRATION FOR BURST POINTS####

#Visualize elevation data to determine whether it has the potential to inform segmentation
ggplot(datsun) +
  geom_point(aes(x=Index, y=R2n, colour=log(Altitude))) +
  scale_colour_viridis_c()

ggplot(datsun) +
  geom_violin(aes(x=season, y=Altitude)) +
  geom_jitter(aes(x=season, y=Altitude, colour=hour)) +
  facet_wrap(~season, scales="free") +
  scale_colour_viridis_c()

#Yes, especially for winter foraging vs roosting

#I think need to do each day separately because are 10 days apart

doys <- datsun %>% 
  dplyr::filter(hour==4) %>% 
  dplyr::select(doy, season) %>% 
  unique()

datburst <- datsun %>% 
  dplyr::filter(hour!=15,
                season=="winter",
                doy %in% doys$doy)

doys.list <- unique(datburst$doy)

#Bayesmove by point----
for(i in 1:lenght(doys.list)){
  
  doy.i <- doys.list[i]
  
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
  dist.bin.lims<- quantile(filter.i$step, c(0,0.05,0.10,0.5,0.75,1), na.rm=T)  #5 bins
  
  # Define bin number and limits for altitude
  alt.bin.lims<- quantile(filter.i$Altitude, c(0,0.05,0.10,0.5,0.75,1), na.rm=T)  #5 bins
  
  set.bin.lims <- quantile(filter.i$sethrs, c(0,0.1,0.25,0.50,0.75,0.9,1), na.rm=T)
  
#  discrete.i<- discrete_move_var(filter.i, lims = list(dist.bin.lims, angle.bin.lims, alt.bin.lims, set.bin.lims), varIn = c("step", "angle", "Altitude", "sethrs"), varOut = c("SL","TA", "AL", "HR"))
  discrete.i<- discrete_move_var(filter.i, lims = list(dist.bin.lims, alt.bin.lims, set.bin.lims), varIn = c("step", "Altitude", "sethrs"), varOut = c("SL", "AL", "HR"))
  
  # Only retain id and discretized step length (SL), turning angle (TA), altitude
#  input.i <- subset(discrete.i, select = c(id, SL, TA, AL, HR))
  input.i <- subset(discrete.i, select = c(id, SL, AL, HR))
  
  #2aii. Classification----
  
  set.seed(1)
  
  # Define model params
  alpha=0.1  #prior
  ngibbs=5000  #number of Gibbs sampler iterations
  nburn=ngibbs/2  #number of burn-in iterations
  nmaxclust=2  #number of maximum possible states (clusters) present
  
  # Run model
  dat.res<- bayesmove::cluster_obs(dat=input.i, alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust,
                                   nburn=nburn)
  
  # Inspect traceplot of log-likelihood
  plot(dat.res$loglikel, type = "l")
  
  # Determine the MAP estimate of the posterior
  MAP.iter<- get_MAP_internal(dat = dat.res$loglikel, nburn = nburn)
  
  # Determine the likely number of behavioral states (i.e., what are the fewest states that account for >= 90% of all observations?)
  theta<- dat.res$theta[MAP.iter,]
  names(theta)<- 1:length(theta)
  theta<- sort(theta, decreasing = TRUE)
  theta %>% cumsum()
  
  # Store cluster order for plotting and behavioral state extraction
  ord<- as.numeric(names(theta))
  
  # Extract bin estimates for each possible state from the `phi` matrix of the model results
  behav.res<- get_behav_hist(dat = dat.res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                             var.names = c("Step Length", "Altitude", "Time since\nsunset"), ord = ord,
                             MAP.iter = MAP.iter)
  
  # Plot state-dependent distributions 
  ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
    geom_bar(stat = 'identity') +
    labs(x = "\nBin", y = "Proportion\n") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x.bottom = element_text(size = 12),
          strip.text = element_text(size = 14),
          strip.text.x = element_text(face = "bold")) +
    scale_fill_manual(values = c(viridis::viridis(3), "grey35","grey35",
                                 "grey35","grey35"), guide = FALSE) +
    scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
    scale_x_continuous(breaks = 1:8) +
    facet_grid(behav ~ var, scales = "free_x")
  
  # Extract MAP estimate of behavioral states
  z<- factor(dat.res$z[[MAP.iter]])
  levels(z)<- ord
  
  # Relabel factor levels
  states.i <- discrete.i %>% 
    st_as_sf(coords=c("x", "y"), crs=3857) %>% 
    st_transform(crs=4326) %>% 
    st_coordinates() %>% 
    cbind(discrete.i) %>% 
    mutate(state = case_when(as.character(z) == 1 ~ "Migrate",
                             as.character(z) == 2 ~ "Forage",
                             as.character(z) == 3 ~ "Roost"))
  
  
  #Plot output
  ggplot() +
    geom_path(data = states.i, aes(x=X, y=Y), color="gray60", size=0.25) +
    geom_point(data = states.i, aes(X, Y, fill=state), size=3, pch=21, alpha=0.7) +
    labs(x = "", y = "") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14)))
  
  ggsave(paste0("Figures/BurstClassification/BurstClassification_BayesMove_",i,".jpeg"), width=8, height=6)
  
  
}

#MoveHMM----

doys <- datsun %>% 
  dplyr::filter(hour==4) %>% 
  dplyr::select(doy, season) %>% 
  unique()

datburst <- datsun %>% 
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