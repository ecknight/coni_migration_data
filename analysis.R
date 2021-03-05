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

#STEP 1. SEGMENT BY SEASON USING NET-SQUARED DISPLACEMENT####

#Load in raw data & filter out pts without enough satellite
dat <- read.csv("PinPoint 2217 2019-07-10 11-21-47.csv") %>% 
  dplyr::filter(Status=="Valid")

#Wrangle for as.traj function. Note filtered to just one point per day
locs <- dat %>% 
  dplyr::select(RTC.date, RTC.time, Latitude, Longitude, Altitude)  %>% 
  mutate(date=as.POSIXct(ymd(RTC.date)),
         time=hms(RTC.time),
         id=2217,
         burst=1) %>% 
  dplyr::filter(hour(time)==15)

#Use as.traj function to get NSD
traj <- as.ltraj(xy=locs[,c("Longitude", "Latitude")],
                 id=locs$id,
                 date=locs$date,
                 burst=locs$burst,
                 proj4string = CRS("+proj=longlat +datum=WGS84"))
#Look at it
traj
head(traj[[1]])

#Segment per location for wintering grounds instead
datseg <- dat %>% 
  left_join(as.data.frame(traj[[1]]) %>%
              mutate(date=as_date(date),
                     doy=yday(date),
                     index=row_number()) %>% 
              dplyr::rename(Longitude=x, Latitude=y)) %>% 
  mutate(season=case_when(R2n<1 ~ "breed",
                          Longitude > -53.33 ~ "winter",
                          R2n>1 & Longitude < -53.33 & year(date)==2018 ~ "fallmig",
                          R2n>1 & Longitude < -53.3 & year(date)==2019 ~ "springmig"),
         season=ifelse(is.na(season), "breed", season),
         R2n=ifelse(is.na(R2n), 0, R2n),
         time=hms(RTC.time),
         hour=hour(time))

#Check segmentation
table(datseg$season)

#Visualize
ggplot(datseg) +
  geom_line(aes(x=index, y=R2n)) +
  geom_point(aes(x=index, y=R2n, colour=season))

#Write out
write.csv(datseg, "PinPoint2217_NSDsegmented.csv", row.names = FALSE)

#STEP 2. SEGMENT STOPOVER FROM MIGRATION - BAYESMOVE - POINT LEVEL####

#https://joshcullen.github.io/bayesmove/articles/Cluster-observations.html

#2a. Fall migration----

#2ai. Wrangling----
#Filter to fall migration points, prepare for bayesmove data prep
datfall <- datseg %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datseg) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="fallmig", hour==15) %>% 
  dplyr::select(id, date, X, Y, Altitude) 

#Calulate step lengths, turning angles, and time steps
stepsfall <- prep_data(dat = datfall, coord.names = c("X","Y"), id = "id")

#Visualize
# View distributions of data streams
hist(stepsfall$step)
hist(stepsfall$angle)

# Create list from data frame
listfall<- df_to_list(dat = stepsfall, ind = "id")
filterfall <- filter_time(dat=listfall, int = 86400) %>% 
  bind_rows

# Define bin number and limits for turning angles
angle.bin.lims<- seq(from=-pi, to=pi, by=pi/4)  #8 bins

# Define bin number and limits for step lengths
dist.bin.lims<- quantile(filterfall$step, c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins

# Define bin number and limits for altitude
alt.bin.lims<- quantile(filterfall$Altitude, c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins

discretefall<- discrete_move_var(filterfall, lims = list(dist.bin.lims, angle.bin.lims), varIn = c("step","angle"), varOut = c("SL","TA"))


# Only retain id and discretized step length (SL), turning angle (TA)
inputfall <- subset(discretefall, select = c(id, SL, TA))

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
theta %>% cumsum()  #first 3 states likely (represent 97.6% of all assigned states)

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
#mehhhh not sold on this output
#should play with bins for step length

#2b. Spring migration----

#2bi. Wrangling----
#Filter to spring migration points, prepare for bayesmove data prep
datspring <- datseg %>% 
  st_as_sf(crs=4326, coords=c("Longitude", "Latitude")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(datseg) %>% 
  mutate(id=2217) %>% 
  dplyr::filter(season=="springmig", hour==15) %>% 
  dplyr::select(id, date, X, Y, Altitude) 

#Calulate step lengths, turning angles, and time steps
stepsspring <- prep_data(dat = datspring, coord.names = c("X","Y"), id = "id")

#Visualize
# View distributions of data streams
hist(stepsspring$step)
hist(stepsspring$angle)

# Create list from data frame
listspring<- df_to_list(dat = stepsspring, ind = "id")
filterspring <- filter_time(dat=listspring, int = 86400) %>% 
  bind_rows

# Define bin number and limits for turning angles
angle.bin.lims<- seq(from=-pi, to=pi, by=pi/4)  #8 bins

# Define bin number and limits for step lengths
dist.bin.lims<- quantile(filterspring$step, c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins

# Define bin number and limits for altitude
alt.bin.lims<- quantile(filterspring$Altitude, c(0,0.25,0.50,0.75,0.90,1), na.rm=T)  #5 bins

discretespring<- discrete_move_var(filterspring, lims = list(dist.bin.lims, angle.bin.lims), varIn = c("step","angle"), varOut = c("SL","TA"))

# Only retain id and discretized step length (SL), turning angle (TA)
inputspring <- subset(discretespring, select = c(id, SL, TA))

#2bii. Classification----

set.seed(1)

# Define model params
alpha=0.1  #prior
ngibbs=5000  #number of Gibbs sampler iterations
nburn=ngibbs/2  #number of burn-in iterations
nmaxclust=2  #number of maximum possible states (clusters) present

# Run model
dat.res<- bayesmove::cluster_obs(dat=inputspring, alpha=alpha, ngibbs=ngibbs, nmaxclust=nmaxclust,
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
statesspring <- discretespring %>% 
  st_as_sf(coords=c("x", "y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  cbind(discretespring) %>% 
  mutate(state = case_when(as.character(z) == 1 ~ "Migration",
                           as.character(z) == 2 ~ "Stopover"))

#Plot output
ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(data = statesspring, aes(x=X, y=Y), color="gray60", size=0.25) +
  geom_point(data = statesspring, aes(X, Y, fill=state), size=1.5, pch=21, alpha=0.7) +
  labs(x = "", y = "") +
  xlim(c(-170, -30)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14)))

#STEP 2. SEGMENT STOPOVER FROM MIGRATION - BAYESMOVE - WITH BREAKPOINTS####

#https://joshcullen.github.io/bayesmove/articles/Segment-the-tracks.html
#https://joshcullen.github.io/bayesmove/articles/Prepare-data-for-analysis.html
#https://joshcullen.github.io/bayesmove/articles/Cluster-track-segments.html

#2a. Fall migration----

#Run step 2bi above first

#2aii. Segmentation----

inputlistfall <-  df_to_list(dat = inputfall, ind = "id")

set.seed(1)

# Define hyperparameter for prior distribution
alpha<- 1
# Set number of iterations for the Gibbs sampler
ngibbs<- 10000
# Set the number of bins used to discretize each data stream
nbins<- c(5,8)

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


#2aiii. Classification----


#STEP 3. SEGMENT FORAGING FROM ROOSTING FROM MIGRATION FOR BURST POINTS####

#Do this separately for each season? Yes, because foraging behaviour could differ between seasons

#Visualize elevation data to determine whether it has the potential to inform segmentation
ggplot(datseg) +
  geom_point(aes(x=Index, y=R2n, colour=log(Altitude))) +
  scale_colour_viridis_c()

ggplot(datseg) +
  geom_violin(aes(x=season, y=Altitude)) +
  geom_jitter(aes(x=season, y=Altitude, colour=hour)) +
  facet_wrap(~season, scales="free") +
  scale_colour_viridis_c()

#Yes, especially for winter foraging vs roosting