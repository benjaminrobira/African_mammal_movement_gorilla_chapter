# ~~~~~~~~~~~~
# WHAT: R code for book chapter "Do seasonal frugivory and cognition shape foraging movements in wild western gorillas?"
# WHO: Benjamin Robira, benjamin.robira@normalesup.org (copyright Benjamin Robira)
# WHEN: November 2022
# ~~~~~~~~~~~~

# What does this script do?
# This script computes all analyses based on R software used in this book chapter. It, however, lacks all data information as those cannot be freely available, due to the vulnerability of the study species.
# It is intended for having access to codes only.

# Found an issue? 
# Please contact me: benjamin.robira@normalesup.org

rm(list = ls())

# Movement data -----------------------------------------------------------

load("Renvironment/Intergroup.RData")

data <- merged
data$Date_group <- paste(data$Ymd, data$Group, sep = "_")

#Preparing data to compute ID/RD
library(lubridate)
library(dplyr)
dataForIDRD_df <- data
dataForIDRD_df$Season <- NA
dataForIDRD_df$Season[month(dataForIDRD_df$Datetime) %in% c(6,7,8,9)] <- "High-frugivory"
dataForIDRD_df$Season[month(dataForIDRD_df$Datetime) %in% c(1,2,3,4,11,12)] <- "Low-frugivory"

for(group in 1:3){
  for(season in 1:3){
    write.table(dataForIDRD_df %>% 
                  filter(
                    Group == unique(Group)[group],
                    Season == unique(Season)[season]
                  ) %>% 
                  dplyr::select(Datetime, Longitude, Latitude,Group,Season), paste("Data/dataForIDRDContinuous_txygs_", unique(dataForIDRD_df$Group)[group], unique(dataForIDRD_df$Season)[season], ".txt", sep = ""), sep =",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
}

library(Rcpp)
Rcpp::sourceCpp("Functions/discretisePathInSpace2.cpp")

library(svMisc)
rediscretisedData <- c()
idToLoop <- unique(data$Date_group)
library(lubridate)
for (i in 1:length(idToLoop)){
  test_transitory <- discretisePathInSpace2(
    longitudesVector = as.vector(data$Longitude[data$Date_group == idToLoop[i]]),
    latitudesVector = as.vector(data$Latitude[data$Date_group == idToLoop[i]]),
    datesVector = decimal_date(data$Datetime[data$Date_group == idToLoop[i]]),
    distanceSampling = 30,
    timeIntervalToConsiderPoints = 30*60,
    displayProgress = FALSE
  )
  
  test_transitory <- as.data.frame(test_transitory)
  colnames(test_transitory) <- c("Longitude", "Latitude", "Datetime")
  test_transitory$Group <- data$Group[data$Date_group == idToLoop[i]][1]
  
  rediscretisedData <- rbind(rediscretisedData,test_transitory)
  progress(i/length(idToLoop)*100)
}

# rediscretisedData <- as.data.frame(rediscretisedData)
# colnames(rediscretisedData) <- c("Longitude", "Latitude", "Datetime", "Group")

rediscretisedData$Datetime <- date_decimal(as.numeric(rediscretisedData$Datetime))
rediscretisedData$Date_group <- paste(date(rediscretisedData$Datetime), rediscretisedData$Group, sep = "_")
head(rediscretisedData)

# Compute weekly UD ------------------------------------------------------

rediscretisedData_forCalculatingArea <- rediscretisedData %>% 
  mutate(
    monthYear = paste(month(Datetime), year(Datetime)),
    monthYearGroup = paste(month(Datetime), year(Datetime), Group),
    weekMonthYearGroup = paste(week(Datetime), month(Datetime), year(Datetime), Group)
    
  ) %>% 
  group_by(date(Datetime), Group) %>% 
  mutate(
    diffTimePrevious = difftime(Datetime, lag(Datetime), units = "hours")
  )
rediscretisedData_forCalculatingArea$Season <- NA
rediscretisedData_forCalculatingArea$Season[month(rediscretisedData_forCalculatingArea$Datetime) %in% c(6,7,8,9)] <- "High-frugivory"
rediscretisedData_forCalculatingArea$Season[month(rediscretisedData_forCalculatingArea$Datetime) %in% c(1,2,3,4,11,12)] <- "Low-frugivory"

library(adehabitatHR)
library(sf)

#Computing area for each week as well as cumulative observation time
toRunThrough <- unique(rediscretisedData_forCalculatingArea$weekMonthYearGroup)
area_v <- rep(NA, times = length(toRunThrough))
area_movement_v <- rep(NA, times = length(toRunThrough))
hoursObs_v <- rep(NA, times = length(toRunThrough))
for(i in 1:length(toRunThrough)){
  print(i/length(toRunThrough)*100)
  dataAsDf <- rediscretisedData_forCalculatingArea %>% 
    filter(
      weekMonthYearGroup == toRunThrough[i]
    )
  if(length(unique(date(dataAsDf$Datetime))) == 7){#Select complete weeks of monitoring
    hoursObs_v[i] = sum(dataAsDf$diffTimePrevious, na.rm = TRUE)
    dataAsDf.sp <- SpatialPoints(coords = dataAsDf[, c(1,2)], 
                                 proj4string = CRS(paste("+proj=utm","+zone=33","+ellps=WGS84", "+datum=WGS84", "+units=m", "+towgs84:0,0,0", sep=" "))
    )
    
    kernel_data <- adehabitatHR::kernelUD(dataAsDf.sp)
    kernel_volume <- getvolumeUD(kernel_data)
    kernel_contour <- getverticeshr(kernel_volume, 95)
    area_v[i] <- st_area(st_as_sf(kernel_contour)) 
    
    
    library(lubridate)
    dataAsDf <- rediscretisedData_forCalculatingArea %>% 
      filter(
        weekMonthYearGroup == toRunThrough[i]
      ) %>% 
      as.data.frame()
    dataAsDf <- as.ltraj(
      xy = dataAsDf[, c(1,2)],
      date = dataAsDf$Datetime,
      id = dataAsDf$Group,
      burst = dataAsDf$Group,
      typeII = TRUE,
      proj4string = CRS(
        paste(
          "+proj=utm",
          "+zone=33",
          "+ellps=WGS84",
          "+datum=WGS84",
          "+units=m",
          "+towgs84:0,0,0",
          sep = " "
        )
      )
    )
    
    coeffD <- BRB.D(dataAsDf, Tmax = 60*60, Lmin = 15)
    kernel_data <- adehabitatHR::BRB(dataAsDf, D = coeffD, Tmax = 60*60, Lmin = 10, hmin = 100)
    kernel_volume <- getvolumeUD(kernel_data)
    kernel_contour <- getverticeshr(kernel_volume, 95)
    
    area_movement_v[i] <- st_area(st_as_sf(kernel_contour)) 
    
  }
}

cor.test(area_movement_v, area_v)

#Testing the difference

library(tidyr)
toBind <- toRunThrough %>% 
  as.data.frame() %>% 
  rename(
    id = "."
  ) %>% 
  separate(
    id, into = c("week", "month", "year", "group"), sep = " ", remove = FALSE
  ) %>%
  mutate(
    season = NA,
    season = ifelse(month %in% c(6,7,8,9), "High-frugivory", season),
    season = ifelse(month %in% c(1,2,3,4,11,12), "Low-frugivory", season)
  )

dfAnalysisUD <- cbind(toBind, area_movement_v, hoursObs_v) %>% 
  as.data.frame() %>% 
  filter(!is.na(season)) %>% 
  mutate(   
    area.log = log(area_movement_v),
    hoursObs.log = log(hoursObs_v)
  ) %>% 
  filter(complete.cases(.)) %>% 
  as.data.frame()

library(lme4)

nrow(dfAnalysisUD)
table(dfAnalysisUD$group)
table(dfAnalysisUD$season)
hist(dfAnalysisUD$hoursObs_v)

dfAnalysisUD$area_movement_v <- dfAnalysisUD$area_movement_v / 100 / 100 #in ha

table(dfAnalysisUD$season)

modelTestUD <- lmer(formula = area_movement_v ~ season + (1|group), offset = hoursObs_v, data = dfAnalysisUD)
diagnostics.plot.dharma(modelTestUD)

summary(modelTestUD)$coefficients
confint(modelTestUD)

# Extract DPL Bai hokou our study ----------------------------------------------

# DPL data
library(dplyr)
DPL.df <- distance.df %>% 
  mutate(
    Loc = ifelse(Group == RC2 | Group == RC1, "RC", "CAR") 
  ) %>% 
  filter(Durationday > 6 & !is.na(DD) & Loc == "CAR")

head(DPL.df)
tail(DPL.df)

## Correct the DPL data with "cleaner" data

library(dplyr)
library(tidyr)
summaryDPL_ourstudy <- DPL.df %>% 
  group_by(
    Loc
  ) %>%
  summarise(
    DDMeanFrug = mean(DD[Season == "FRUGIVORY" & !is.na(Season)]),
    DDMeanFol = mean(DD[Season == "FOLIVORY" & !is.na(Season)]),
    DDMean = mean(DD),
    DDMinFrug = min(DD[Season == "FRUGIVORY" & !is.na(Season)]),
    DDMinFol = min(DD[Season == "FOLIVORY" & !is.na(Season)]),
    DDMin = min(DD),
    DDMaxFrug = max(DD[Season == "FRUGIVORY" & !is.na(Season)]),
    DDMaxFol = max(DD[Season == "FOLIVORY" & !is.na(Season)]),
    DDMax = max(DD),
    sd = sd(DD),
    NobsFrug = length(DD[Season == "FRUGIVORY" & !is.na(Season)]),
    NobsFol = length(DD[Season == "FOLIVORY" & !is.na(Season)]),
    Nobs = length(DD)
  )

toBind1 <- summaryDPL_ourstudy %>% dplyr::select(Loc, DDMinFrug, DDMaxFrug, DDMeanFrug, sd, NobsFrug)
colnames(toBind1) <- c("Loc", "Min", "Max", "Mean", "sd", "Nobs")
toBind2 <- summaryDPL_ourstudy %>% dplyr::select(Loc, DDMinFol, DDMaxFol, DDMeanFol, sd, NobsFol)
colnames(toBind2) <- c("Loc", "Min", "Max", "Mean", "sd", "Nobs")
summaryDPL_ourstudy_season <- rbind(toBind1, toBind2)
summaryDPL_ourstudy_season$colour = c("salmon", 
                                      #"salmon", 
                                      #"darkgreen", 
                                      "darkgreen")
summaryDPL_ourstudy_season$BackgroundColour <- summaryDPL_ourstudy_season$colour

testCAR <- t.test(x = DPL.df$DD[DPL.df$Season == "FRUGIVORY" & DPL.df$Loc == "CAR"], 
                  y = DPL.df$DD[DPL.df$Season == "FOLIVORY" & DPL.df$Loc == "CAR"])

library(lme4)
DPL.df$DD.log <- log(DPL.df$DD)
testCARlin <- lmer(DD.log ~ Season + (1 | Group), 
                   data = DPL.df[!is.na(DPL.df$Season) & DPL.df$Loc == "CAR",])
#Difference in m
exp(summary(testCARlin)$coefficients[1,1] + summary(testCARlin)$coefficients[2,1]) - exp(summary(testCARlin)$coefficients[1,1])
#CI difference
exp(confint(testCARlin)[3,1] + confint(testCARlin)[4,1]) - exp(confint(testCARlin)[3,1])
exp(confint(testCARlin)[3,2] + confint(testCARlin)[4,2]) - exp(confint(testCARlin)[3,2])

diagnostics.plot.dharma(testCARlin)

# Plotting the DPL at various field sites ---------------------------------

## Enter the data
comparativeTableDPL <- data.frame(
  Author = c("Doran-Sheehy et al.",
             "Etiendem & Tagg",
             "Cipolletta",
             "Seiler & Robbins",
             "Wangue et al.",
             "Tutin",
             "Goldsmith"),
  Year = c("2004",
           "2013",
           "2004",
           "2020",
           "2015",
           "1996",
           "1999"),
  Species = c("Gorilla g. gorilla",
              "Gorilla g. diehli",
              "Gorilla g. gorilla",
              "Gorilla g. gorilla",
              "Gorilla g. gorilla",
              "Gorilla g. gorilla",
              "Gorilla g. gorilla"),
  Fieldsite = c("Mondika",
                "Mawambi",
                "Bai Hokou",
                "Loango",
                "Dipikar Island",
                "Lopé",
                "Bai Hokou"),
  Country = c("Republic of Congo",
              "Cameroon",
              "Central African Republic",
              "Gabon",
              "Cameroon",
              "Gabon",
              "Central African Republic"),
  Nobs = c("334",
           "86",
           "649",
           "239",
           "502",
           "80",
           "95"),
  Mean = c("2014",
           "1642",
           "1717",
           "2562",
           "1600",
           "1105",
           "2600"),
  Min = c("400",
          "451",
          "250",
          "354",
          "180",
          "220",
          "300"),
  Max = c("4860",
          "6400",
          "4300",
          "7575",
          "5300",
          "2790",
          "5000"),
  Seasonaleffect = c("yes",
                     "no",
                     "yes",
                     "yes",
                     NA,
                     "yes",
                     "yes")
)

## Process data

library(dplyr)
comparativeTableDPL <- comparativeTableDPL %>% 
  mutate(
    Nobs = as.numeric(Nobs),
    Mean = as.numeric(Mean),
    Min = as.numeric(Min),
    Max = as.numeric(Max),
    AuthorYear = paste(Author, Year, sep=", "),
    Colour = ifelse(!is.na(Seasonaleffect) & Seasonaleffect == "yes", "black", "darkgrey"),
    BackgroundColour = ifelse(is.na(Seasonaleffect), "white", Colour)
    
  ) %>% 
  arrange(
    Species,
    Country,
    Fieldsite,
    Year
  ) %>% 
  mutate(
    y = rev(1:length(Author)) + 2
  )

## Add our own data 
toAdd1 <- unlist(as.vector(summaryDPL_ourstudy[1, c(12,3,6,9)+1]))
names(toAdd1) <- NULL

comparativeTableDPL <- rbind(
  comparativeTableDPL,
  c("Robira et al", "this study", "Gorilla g. gorilla", "Bai Hokou", "Central African Republic", floor(toAdd1), "yes", "Robira et al., this study", "black", "black", 2)
)
comparativeTableDPL <- comparativeTableDPL %>% 
  mutate(
    Nobs = as.numeric(Nobs),
    Mean = as.numeric(Mean),
    Min = as.numeric(Min),
    Max = as.numeric(Max),
    y = as.numeric(y)
  )
library(ggplot2)
DPLPlot <- ggplot(
  data = comparativeTableDPL,
  aes(
    x = Mean,
    y = y,
    fill = BackgroundColour,
    color = Colour
  )
) +
  # Range
  geom_errorbar(aes(xmin  = Min,
                    xmax  = Max,
                    color = Colour),
                width = 0.15,
                size  = 0.7) +
  # Mean
  geom_point(aes(size = log(Nobs)*2,
                 fill = BackgroundColour,
                 color = Colour),
             pch = 21) +
  # Seasonal effect
  geom_point(
    data = summaryDPL_ourstudy_season,
    aes(x = Mean,
        y = rev(as.numeric(as.factor(Loc))) + 1 + c(0.2, #Given I have withdrawn mondika data, I have added the +1
                                                    -0.2),
        size = log(Nobs)
    ),
    fill = c("salmon",
             "darkgreen"), 
    colour = c("salmon",
               "darkgreen"), 
    pch = 21,
    show.legend = FALSE
  ) +
  geom_errorbar(data = summaryDPL_ourstudy_season, 
                aes(xmin  = Mean - qt(0.05/2, Nobs)*sd/sqrt(Nobs - 1),
                    xmax  = Mean + qt(0.05/2, Nobs)*sd/sqrt(Nobs - 1),
                    y = rev(as.numeric(as.factor(Loc))) + 1 + c(0.2, 
                                                                -0.2), fill = NA),
                color = c("salmon",
                          "darkgreen"), 
                width = 0.15,
                size  = 0.7) +
  # Min, Max, Mean text
  geom_text(
    aes(
      x = Min,
      y = y,
      label = Min
    ),
    size = 3,
    nudge_y = 0.25
  ) +
  geom_text(
    aes(
      x = Max,
      y = y,
      label = Max
    ),
    size = 3,
    nudge_y = 0.25
  ) +
  geom_text(
    aes(
      x = Mean,
      y = y,
      label = Mean
    ),
    size = 4,
    nudge_y = 0.3
  ) +
  scale_colour_manual(values = c("black", "darkgrey")) +
  scale_fill_manual(values = c("black", "darkgrey", "white"),
                    labels = c("Yes", "No", "Not available"),
                    name = "Seasonal effect") +
  guides(colour = "none",
         size = "none") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  
  scale_y_continuous(breaks=1:nrow(comparativeTableDPL) + 1, labels = rev(paste(# Given i have withdrawn Mondika data, I have added + !
    comparativeTableDPL$AuthorYear,
    paste(comparativeTableDPL$Fieldsite, comparativeTableDPL$Country, sep = ", "),
    sep = "\n"
  ))) +
  xlim(0, round(max(comparativeTableDPL$Max/1000, na.rm = TRUE))*1000) +
  ylab("") +
  xlab("Daily path length (m)")

DPLPlot

# Analysing the movement pattern as a function of season ------------------

## Preparing data as tracks ------------------------------------------------

DPL.df$dateGroup <- paste(DPL.df$Date, DPL.df$Group, sep = "_")
rediscretisedData <- rediscretisedData %>% 
  filter(
    Date_group %in% DPL.df$dateGroup#Keep days for which we have sufficient monitoring
  )


library(amt)
crs = CRS(paste("+proj=utm","+zone=33","+ellps=WGS84", "+datum=WGS84", "+units=m", "+towgs84:0,0,0", sep = " ")))

## Extract movement metrics

### Group data
trk_grouped <- trk %>% nest(data = -"id")

#Rediscretize data at constant time step and calculate movement var (speed, heading angle and step length)
trk_rediscretized_bygroup <- trk_grouped %>% 
  dplyr::mutate(steps = map(data, function(x) {
    if(nrow(x) > 1){#must include 2 locs
      x %>% amt::steps()} else{
      }
  }
  ))
trk_rediscretized_bygroup_unnest <- trk_rediscretized_bygroup[
  sapply(trk_rediscretized_bygroup$steps, function(x)length(x)) > 0,
]

### un-nest the data
trk_rediscretized_bygroup_unnest <- trk_rediscretized_bygroup %>% unnest(cols = steps) 

### Match season
trk_rediscretized_bygroup_unnest$Season <- NA
trk_rediscretized_bygroup_unnest$Season[month(trk_rediscretized_bygroup_unnest$t1_) %in% c(6,7,8,9)] <- "High-frugivory"
trk_rediscretized_bygroup_unnest$Season[month(trk_rediscretized_bygroup_unnest$t1_) %in% c(1,2,3,4,11,12)] <- "Low-frugivory"

## Save data per group (for calculating ID/RD) ---------------

dataForIDRD_df <- trk_rediscretized_bygroup_unnest %>% 
  separate(id, into = c("Date", "Group"), sep = "_") %>% 
  dplyr::select(
    t1_, x1_, y1_, Group, Season
  )

for(group in 1:3){
  for(season in 1:3){
    write.table(dataForIDRD_df %>% 
                  filter(
                    Group == unique(Group)[group],
                    Season == unique(Season)[season]
                  ), paste("Data/DataforIDRD_txygs_", unique(dataForIDRD_df$Group)[group], unique(dataForIDRD_df$Season)[season], ".txt", sep = ""), sep =",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
}

## Analysing sinuosity of movements ----------------------------------------

library(dplyr)
trk_rediscretized_bygroup_unnest_boxplot <- trk_rediscretized_bygroup_unnest %>% 
  mutate(
    costa_ = cos(ta_),
    costa_ = ifelse(costa_ == 0, NA, costa_) #Remove if due to linear interpolation
  ) %>% 
  separate(id, into = c("Date", "Group"), sep = "_") %>% 
  group_by(Date, Group, Season) %>% 
  summarise(
    meancosta_ = mean(costa_, na.rm = TRUE),
    Nobs = length(costa_[!is.na(costa_)]),
    Displacement = sqrt((x2_[length(x2_)] - x1_[1])**2 + (y2_[length(y2_)] - y1_[1])**2)
  ) %>% 
  ungroup() %>% 
  filter(
    !is.na(Season) &
      Nobs > 10
  )

min(trk_rediscretized_bygroup_unnest_boxplot$meancosta_, na.rm = TRUE)
max(trk_rediscretized_bygroup_unnest_boxplot$meancosta_, na.rm = TRUE)

trk_rediscretized_bygroup_unnest_boxplot[which(trk_rediscretized_bygroup_unnest_boxplot$meancosta_ == min(trk_rediscretized_bygroup_unnest_boxplot$meancosta_, na.rm = TRUE)),]

## generalised linear model

linear_models <- trk_rediscretized_bygroup_unnest_boxplot %>% 
  mutate(
    Nobssqrt = sqrt(Nobs)
  ) %>% 
  filter(complete.cases(.)) %>% 
  group_by(Group) %>% 
  do(anova(lm(meancosta_ ~ Season + Nobssqrt, data = .), lm(meancosta_ ~ Nobssqrt, data = .), test = "F")) %>% 
  filter(!is.na(Df)) 

#Check assumptions
trk_rediscretized_bygroup_unnest_boxplot %>% 
  mutate(
    Nobssqrt = sqrt(Nobs)
  ) %>% 
  filter(complete.cases(.)) %>% 
  group_by(Group) %>% 
  do(diagnostics.plot.dharma(lm(meancosta_ ~ Season + Nobssqrt, data = .)))

source("Functions/toolbox.R")
labelsTest <- linear_models %>% 
  mutate(
    PvalueText = pvalueToText(`Pr(>F)`)
  ) 
# Plot boxplot
n_fun <- function(x){
  return(data.frame(y = min(x) - 0.05,
                    label = length(x)))
}

Linearityplot <- ggplot(data = trk_rediscretized_bygroup_unnest_boxplot,
                        aes(x = Group, y = meancosta_, fill = Season)) +
  stat_boxplot(geom ='errorbar', width = 0.6, size = 0.25) +
  geom_boxplot(width = 0.6) +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group = Season),
               hjust = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = "mean", aes(group = Season), pch = 21, fill = "white", size = 1, position = position_dodge(0.6)) +
  annotate("text", x = 1:3, y = 1.05, label = labelsTest$PvalueText, size = 5) +
  scale_fill_manual(values = c("salmon", "darkgreen")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  scale_x_discrete(labels = c("A", "B", "C")) +
  xlab("Gorilla group") +
  # ylim(-0.5,1.1) +
  scale_y_continuous(minor_breaks = seq(0, 60, 0.05)) +
  ylab("Daily mean cosine\nof turning angles")

Linearityplot

# Analysis of displacement ------------------------------------------------

## Displacement as a function of season ------------------------------------

n_fun <- function(x){
  return(data.frame(y = -100,
                    label = length(x)))
}

linear_models <- trk_rediscretized_bygroup_unnest_boxplot %>% 
  mutate(
    Nobssqrt = sqrt(Nobs)
  ) %>% 
  filter(complete.cases(.)) %>% 
  group_by(Group) %>% 
  do(anova(lm(Displacement ~ Season + Nobssqrt, data = .), lm(Displacement ~ Nobssqrt, data = .), test = "F")) %>% 
  filter(!is.na(Df)) 

#Check assumptions
trk_rediscretized_bygroup_unnest_boxplot %>% 
  mutate(
    Nobssqrt = sqrt(Nobs)
  ) %>% 
  filter(complete.cases(.)) %>% 
  group_by(Group) %>% 
  do(diagnostics.plot.dharma(lm(Displacement ~ Season + Nobssqrt, data = .)))

labelsTest <- linear_models %>% 
  mutate(
    PvalueText = pvalueToText(`Pr(>F)`)
  ) 

Displacementplot <- ggplot(data = trk_rediscretized_bygroup_unnest_boxplot,
                           aes(x = Group, y = Displacement, fill = Season)) +
  stat_boxplot(geom ='errorbar', width = 0.6, size = 0.25) +
  geom_boxplot(width = 0.6) +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group = Season), 
               hjust = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = "mean", aes(group = Season), pch = 21, fill = "white", size = 1, position = position_dodge(0.6)) +
  annotate("text", x = 1:3, y = 4000, label = labelsTest$PvalueText, size = 5) +
  scale_fill_manual(values = c("salmon", "darkgreen")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  scale_x_discrete(labels = c("A", "B", "C")) +
  xlab("Gorilla group") +
  #ylim(-100,4000) +
  scale_y_continuous(minor_breaks = seq(-200,4200, 200)) +
  ylab("Displacement (m)")

Displacementplot

## Displacement as a function of daily path length (DD) --------------------

### Fitting Displacement ~ DD: LOESS
DPL.df$Date_group <- paste(DPL.df$Date, DPL.df$Group, sep = "_")
trk_rediscretized_bygroup_unnest_boxplot$Date_group <- paste(trk_rediscretized_bygroup_unnest_boxplot$Date, trk_rediscretized_bygroup_unnest_boxplot$Group, sep = "_")

dataForDisplacement <- left_join(
  trk_rediscretized_bygroup_unnest_boxplot,
  DPL.df %>%  dplyr::select(DD, Date_group),
  by = c("Date_group" = "Date_group")
)
head(dataForDisplacement) 

#### Frugivory
loessDisplacementFrug <- loess(Displacement ~ DD, dataForDisplacement %>%  filter(Season == "High-frugivory"), span = 0.50)
predictedLoessFrug <- predict(loessDisplacementFrug,  data = dataForDisplacement, se = TRUE)

loessDisplacementFrug <- cbind(
  loessDisplacementFrug$x,
  loessDisplacementFrug$y,
  loessDisplacementFrug$fitted,
  predictedLoessFrug$se.fit
) %>%  as.data.frame()
colnames(loessDisplacementFrug) <- c("x", "y", "fitted", "se")
loessDisplacementFrug <- loessDisplacementFrug %>%  arrange(x)

loessDisplacementFrugCI <- cbind(
  c(loessDisplacementFrug$x, rev(loessDisplacementFrug$x)),
  c(loessDisplacementFrug$fitted, rev(loessDisplacementFrug$fitted)) + c(loessDisplacementFrug$se, - rev(loessDisplacementFrug$se))
) %>% as.data.frame()
colnames(loessDisplacementFrugCI) <- c("x", "y")

#### Folivory

loessDisplacementFol <- loess(Displacement ~ DD, dataForDisplacement %>%  filter(Season == "Low-frugivory"), span = 0.50)
predictedLoessFol <- predict(loessDisplacementFol,  data = dataForDisplacement, se = TRUE)

loessDisplacementFol <- cbind(
  loessDisplacementFol$x,
  loessDisplacementFol$y,
  loessDisplacementFol$fitted,
  predictedLoessFol$se.fit
) %>%  as.data.frame()
colnames(loessDisplacementFol) <- c("x", "y", "fitted", "se")
loessDisplacementFol <- loessDisplacementFol %>%  arrange(x)

loessDisplacementFolCI <- cbind(
  c(loessDisplacementFol$x, rev(loessDisplacementFol$x)),
  c(loessDisplacementFol$fitted, rev(loessDisplacementFol$fitted)) + c(loessDisplacementFol$se, - rev(loessDisplacementFol$se))
) %>% as.data.frame()
colnames(loessDisplacementFolCI) <- c("x", "y")

### Plotting

displacementAndDPLplot <- ggplot(dataForDisplacement, aes(x = Displacement, y = DD)) +
  
  #Frugivorous season
  geom_polygon(loessDisplacementFrugCI, mapping = aes(x = x,
                                                      y = y), fill = "salmon", alpha = 0.25) + # Se
  geom_point(loessDisplacementFrug, mapping = aes(x = x, y = y), pch = 21, col = "salmon", fill = adjustcolor("salmon", alpha.f = 0.25)) + #Raw data
  geom_line(loessDisplacementFrug, mapping = aes(x = x, y = fitted), col="salmon") + # Model
  
  #Folivorous season
  geom_polygon(loessDisplacementFolCI, mapping = aes(x = x,
                                                     y = y), fill = "darkgreen", alpha = 0.25) + # Se
  geom_point(loessDisplacementFol, mapping = aes(x = x, y = y), pch = 21, col = "darkgreen", fill = adjustcolor("darkgreen", alpha.f = 0.25)) + #Raw data
  geom_line(loessDisplacementFol, mapping = aes(x = x, y = fitted), col="darkgreen") + # Model
  
  #Diagonal line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  xlab("Daily path length (m)") +
  ylab("Displacement (m)") +
  scale_x_continuous(expand = c(0,0), 
                     minor_breaks = seq(-250, 9000, 250)) +
  scale_y_continuous(expand = c(0,0), 
                     minor_breaks = seq(-250, 9000, 250))

displacementAndDPLplot

# Travelling speed --------------------------------------------------------

distanceForMovingInM = 30

library(Rcpp)
Rcpp::sourceCpp("Functions/resampleConstantTime2.cpp")

library(svMisc)
rediscretisedDataTime <- c()
idToLoop <- unique(data$Date_group)
library(lubridate)
for (i in 1:length(idToLoop)){
  test_transitory <- resampleConstantTime2(
    longitudesVector = as.vector(data$Longitude[data$Date_group == idToLoop[i]]),
    latitudesVector = as.vector(data$Latitude[data$Date_group == idToLoop[i]]),
    datesVector = decimal_date(data$Datetime[data$Date_group == idToLoop[i]]),
    dateInterval = 1*60,
    dateThreshold = 30*60,
    displayProgress = FALSE
  )
  
  test_transitory <- as.data.frame(test_transitory)
  colnames(test_transitory) <- c("Longitude", "Latitude", "Datetime")
  test_transitory$Group <- data$Group[data$Date_group == idToLoop[i]][1]
  
  # check <- data[data$Date_group == idToLoop[i],]
  # plot(check$Longitude, check$Latitude, type = "l")
  # lines(test_transitory$Longitude, test_transitory$Latitude, col = "red")
  
  rediscretisedDataTime <- rbind(rediscretisedDataTime,test_transitory)
  progress(i/length(idToLoop)*100)
}

library(zoo)
library(lubridate)
dataForSpeed <-
  rediscretisedDataTime %>% 
  mutate(
    Datetime = date_decimal(Datetime),
    Ymd = date(Datetime)
  ) %>% 
  group_by(Group, Ymd) %>% 
  mutate(
    rowNB = 1:length(Longitude),
    
    #Moved if distance > 30m
    previousLongitude = dplyr::lag(Longitude),
    previousLatitude = dplyr::lag(Latitude),
    previousDatetime = dplyr::lag(Datetime),
    distanceToPrevious = sqrt((Longitude - previousLongitude)**2 + (Latitude - previousLatitude)**2),
    hasMoved = ifelse(distanceToPrevious <= distanceForMovingInM, "NO", "YES"),
    previousHasMoved = dplyr::lag(hasMoved),
    subsequentHasMoved = dplyr::lead(hasMoved),
    
    #Add start and end of move
    startMove = ifelse(is.na(hasMoved) | (hasMoved == "NO" & subsequentHasMoved == "YES"), "YES", "NO"),
    endMove = ifelse((rowNB == max(rowNB) | hasMoved == "YES" & subsequentHasMoved == "NO"), "YES", "NO"),
    
    #Correct first and end of move if necessary: double start/double end
    startMove = ifelse(length(which(startMove == "YES")) > 1 & rowNB == 1 & which(startMove == "YES")[min(2, length(which(startMove == "YES")) + 1)] < which(endMove == "YES")[1], "NO", startMove),
    endMove = ifelse(length(which(endMove == "YES")) > 1 & rowNB == max(rowNB) & which(startMove == "YES")[length(which(startMove == "YES"))] < which(endMove == "YES")[max(1, length(which(endMove == "YES")) - 1)], "NO", endMove),
    endMove = ifelse(is.na(endMove), "NO", endMove),#Correct first point
    startMove = ifelse(is.na(startMove), "NO", startMove)#Correct last point
  ) %>% 
  ungroup() %>% 
  filter(!(startMove == "YES" & endMove == "YES")) %>% 
  mutate(
    segmentID = ifelse(startMove == "YES" | endMove == "YES", "validSegment", "notValidSegment")
  ) %>% 
  group_by(segmentID) %>% 
  mutate(
    segmentID = 1:length(segmentID)
  ) %>% 
  ungroup() %>% 
  mutate(
    #Group variable by moving steps
    segmentID = ifelse(startMove == "YES" | endMove == "YES", segmentID, NA), #Keep only entrance exit
    row = 1:nrow(.),
    segmentID = ifelse(row == 1, 0, segmentID),
    segmentID = na.locf(segmentID), #Fill NA
    segmentID = ifelse(segmentID %in% unique(segmentID[startMove == "YES"]), segmentID, NA), #Remove points that are not movement
    segmentID = ifelse(is.na(segmentID) & !is.na(lag(segmentID)), lag(segmentID), segmentID)#Include last point of movement
  )

# check <- dataForSpeed[dataForSpeed$Group == CAR1 & dataForSpeed$Ymd == "2016-06-15",]
# View(check)

dataForSpeedSummarisedConstantTimeApproach <- dataForSpeed %>% 
  # Keep only movement
  filter(!is.na(segmentID)) %>% 
  mutate(
    travelledDistance = sqrt((Longitude - previousLongitude)**2 + (Latitude - previousLatitude)**2),
    timediff = difftime(Datetime, previousDatetime, units = "sec"),
    date = date(Datetime)
  ) %>% 
  group_by(segmentID, Group, date) %>% 
  summarise(
    travelledDistance = sum(travelledDistance),
    timediff = sum(as.numeric(timediff))
  ) %>% 
  ungroup() %>% 
  group_by(date, Group) %>% 
  summarise(
    speed = sum(travelledDistance, na.rm = TRUE)/sum(timediff, na.rm = TRUE)
  ) %>% 
  mutate(
    Season = NA,
    Season = ifelse(month(date) %in% c(6,7,8,9), "High-frugivory", Season),
    Season = ifelse(month(date) %in% c(11,12,1,2,3,4), "Low-frugivory", Season)
  )

dataForSpeedConstantStepApproach <- trk_rediscretized_bygroup_unnest %>% 
  mutate(
    speed = sl_/as.numeric(dt_)
  )
hist(log(dataForSpeedConstantStepApproach$speed))
exp(-2)

dataForSpeedConstantStepApproach <- trk_rediscretized_bygroup_unnest %>% 
  mutate(
    speed = sl_/as.numeric(dt_)
  ) %>% 
  filter(
    speed > distanceForMovingInM/60 #Move from more radius of canopy in 1 min
  ) %>% 
  group_by(id) %>% 
  summarise(
    speed = sum(sl_)/sum(as.numeric(dt_)),
    countPosition = length(sl_)
  ) %>% 
  filter(
    countPosition >= 10
  ) %>% 
  separate(col = id, into = c("date", "Group"), sep = "_", remove = FALSE) %>% 
  ungroup() %>% 
  mutate(
    Season = NA,
    Season = ifelse(month(date) %in% c(6,7,8,9), "High-frugivory", Season),
    Season = ifelse(month(date) %in% c(11,12,1,2,3,4), "Low-frugivory", Season)
  ) %>% 
  arrange(Group, date)

dataForSpeedSummarisedConstantTimeApproach <- dataForSpeedSummarisedConstantTimeApproach %>% 
  mutate(Dategroup = paste(date, Group, sep = "_"))
dataMergedSpeed <- left_join(dataForSpeedConstantStepApproach,
                             dataForSpeedSummarisedConstantTimeApproach,
                             by = c("id" = "Dategroup"))

#Check correlation speed between the two methods
plot(dataMergedSpeed$speed.x,
     dataMergedSpeed$speed.y)

dataForSpeedSummarised <- dataForSpeedConstantStepApproach #dataForSpeedSummarisedConstantTimeApproach

linear_models <- dataForSpeedSummarised %>% 
  mutate(speed_log = log(speed)) %>% 
  filter(!is.na(Season)) %>% 
  ungroup() %>% 
  filter(complete.cases(.)) %>% 
  group_by(Group) %>% 
  do(anova(lm(speed_log ~ Season, data = .), lm(speed_log ~ 1, data = .), test = "F")) %>% 
  filter(!is.na(Df)) 

#Check assumptions
dataForSpeedSummarised %>% 
  mutate(speed_log = log(speed)) %>% 
  filter(!is.na(Season)) %>% 
  ungroup() %>% 
  filter(complete.cases(.)) %>% 
  group_by(Group) %>% 
  do(diagnostics.plot.dharma(lm(speed_log ~ Season, data = .)))

source("Functions/toolbox.R")
labelsTest <- linear_models %>% 
  mutate(
    PvalueText = pvalueToText(`Pr(>F)`)
  ) 
# Plot boxplot
n_fun <- function(x){
  return(data.frame(y = min(x) - 0.025,
                    label = length(x)))
}

Speedplot <- ggplot(data = dataForSpeedSummarised %>% 
                      filter(!is.na(Season)),
                    aes(x = Group, y = speed, fill = Season)) +
  stat_boxplot(geom = 'errorbar', width = 0.6, size = 0.25) +
  geom_boxplot(width = 0.6) +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group = Season),
               hjust = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = "mean", aes(group = Season), pch = 21, fill = "white", size = 1, position = position_dodge(0.6)) +
  annotate("text", x = 1:3, y = 1.35, label = labelsTest$PvalueText, size = 5) +
  scale_fill_manual(values = c("salmon", "darkgreen")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  scale_x_discrete(labels = c("A", "B", "C")) +
  xlab("Gorilla group") +
  #ylim(0.4,1.5) +
  scale_y_continuous(minor_breaks = seq(0.4,1.5, 0.05)) +
  ylab("Travelling speed (m/s)")

Speedplot


Speedplot_constantTime <- ggplot(data = dataForSpeedSummarisedConstantTimeApproach %>% 
                                   filter(!is.na(Season)),
                                 aes(x = Group, y = speed, fill = Season)) +
  stat_boxplot(geom ='errorbar', width = 0.6, size = 0.25) +
  geom_boxplot(width = 0.6) +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group = Season),
               hjust = 0.5, position = position_dodge(0.6)) +
  stat_summary(fun = "mean", aes(group = Season), pch = 21, fill = "white", size = 1, position = position_dodge(0.6)) +
  annotate("text", x = 1:3, y = 1.35, label = labelsTest$PvalueText, size = 5) +
  scale_fill_manual(values = c("salmon", "darkgreen")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  scale_x_discrete(labels = c("A", "B", "C")) +
  xlab("Gorilla group") +
  #ylim(0.4,1.5) +
  scale_y_continuous(minor_breaks = seq(0.4,1.5, 0.05)) +
  ylab("Travelling speed (m/s)")

Speedplot_constantTime

library(ggpubr)
ggarrange(Speedplot,
          Speedplot_constantTime,
          common.legend = TRUE)

# Analysing recursion to feeding locations --------------------------------

### Extract the feeding locations

library(readr)
data_behaviour <- read_delim("Data/Data_behaviour.txt", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE, guess_max=10000)
head(data_behaviour)

visitedTrees_df <- data_behaviour %>% # Select only fruit, young leaves qnd trees
  filter(
    Food_type == "fruit" |
      Food_type == "unripe_fruit" |
      Food_type == "young_leaves" |
      Food_type == "leaves"
  ) %>% 
  filter(!is.na(Longitude))

visitedTrees_df$Food_type2 <- visitedTrees_df$Food_type
visitedTrees_df$Food_type2[visitedTrees_df$Food_type == "unripe_fruit"] <- "fruit"

### Barycentring feeding locations associated to same food that are close

#Cluster trees
visitedTrees_df$Rownb <- 1:nrow(visitedTrees_df)
visitedTrees_df$Assigned_tree <- NA
treeLoc.df <- matrix(NA, ncol=3, nrow=nrow(visitedTrees_df))

clusterTree = 30 #in m
library(svMisc)
for(i in 1:nrow(visitedTrees_df)){
  if(is.na(visitedTrees_df$Assigned_tree[i])){
    visitedTrees_dftoUtilise <- visitedTrees_df[is.na(visitedTrees_df$Assigned_tree)&
                                                  visitedTrees_df$Food_type2 == visitedTrees_df$Food_type2[i],]
    
    visitedTrees_dftoUtilise <- visitedTrees_dftoUtilise[
      (visitedTrees_dftoUtilise$Latitude - visitedTrees_df$Latitude[i])**2 +
        (visitedTrees_dftoUtilise$Longitude - visitedTrees_df$Longitude[i])**2 <= clusterTree*clusterTree,
    ]
    
    treeLoc.df[i,] <- c(visitedTrees_df$Food_type2[i], 
                        mean(visitedTrees_dftoUtilise$Longitude, na.rm=TRUE),
                        mean(visitedTrees_dftoUtilise$Latitude, na.rm=TRUE))
    visitedTrees_df$Assigned_tree[visitedTrees_dftoUtilise$Rownb] <- "YES"
  }
  progress(i/nrow(visitedTrees_df)*100)
}

treeLoc.df <- as.data.frame(treeLoc.df)
colnames(treeLoc.df) <- c("Food", "Longitude", "Latitude")
treeLoc.df <- treeLoc.df[complete.cases(treeLoc.df),]
nrow(treeLoc.df)

treeLoc.df$Longitude <- as.numeric(treeLoc.df$Longitude)
treeLoc.df$Latitude <- as.numeric(treeLoc.df$Latitude)

### Extract recursion to feeding location
library(recurse)
recursionsTransitory.df <- getRecursionsAtLocations(x = data %>% dplyr::select(Longitude, Latitude, Datetime, Group) %>% filter(Group %in% c(CAR3, CAR1, CAR2)) %>% as.data.frame(), 
                                                    locations = treeLoc.df %>% dplyr::select(Longitude, Latitude) %>% as.data.frame(), 
                                                    radius = 30, threshold = 0.5, timeunits = "days")

### Distribution of visit number for each food item
obsDays_df <- merged %>%
  group_by(Group) %>% 
  summarise(
    Nobscompletedays = length(unique(Ymd))
  )

dataForPlotDensity <-
  recursionsTransitory.df$revisitStats %>%
  group_by(id, x, y) %>% 
  summarise(
    revisits = length(x),
    Group = unique(id)[1],
    Tree = paste(floor(unique(x))[1], floor(unique(y))[1], sep = "_")
  ) %>% 
  #Add food type
  left_join(
    .,
    treeLoc.df %>% mutate(Tree = paste(floor(Longitude), floor(Latitude), sep = "_")) %>% dplyr::select(Food, Tree),
    by = c("Tree" = "Tree")
  ) %>% 
  #Add monitoring duration
  left_join(
    .,
    obsDays_df,
    by = c("Group" = "Group")
  ) %>% 
  filter(!is.na(Food)) %>% 
  # Weight visit number by sampling effort
  mutate(
    revisits = revisits/Nobscompletedays
  )

round(summary(dataForPlotDensity$revisits)*365.25, digits = 1)

densityDistributionVisit_total <- ggplot(dataForPlotDensity, 
                                         aes(x = revisits, fill = Food, col = Food)) +
  geom_density(alpha = 0.4) +
  geom_density(alpha = 0) +
  guides(colour = "none") +
  scale_fill_manual(name="Food type", values = c("salmon", "darkgreen", "darkolivegreen3"), labels = c("Fruit", "(Mature) leaves", "Young leaves")) +
  scale_colour_manual(name="Food type", values = c("salmon", "darkgreen", "darkolivegreen3"), labels = c("Fruit", "(Mature) leaves", "Young leaves")) +
  theme_bw() +
  theme(legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 17,face = "bold"),
        panel.grid.minor=element_line(colour = "grey93"),
        panel.grid.major=element_line(colour = "grey93"),
        # Facet wrap layout
        strip.text = element_text(size = 13, face="bold"),
        strip.background = element_blank()) +
  labs(y = "Proportion of sample", x = "Frequency of revisits at\nthe feeding site (per obs. days)") +
  scale_x_continuous(expand = c(0, 0), minor_breaks = seq(0, 100, 0.005)) +
  scale_y_continuous(expand = c(0, 0), minor_breaks = seq(0, 100, 2.5))


densityDistributionVisit_total

### Calculate indices of recursion as in Jodie Martin et al. on frequently visited trees, within season

min(recursionsTransitory.df$revisitStats$entranceTime)
max(recursionsTransitory.df$revisitStats$entranceTime)

#Should have: 17 months (or 16 if substraction of one), 7 of frug, 7 of fol

dfVisitsSummary <- recursionsTransitory.df$revisitStats %>% 
  filter(entranceTime > dmy_hms("01-07-2016 00:00:00"))  %>% #Remove first month to avoid biasing calculations
  mutate(
    x_y = paste(x, y, sep = "_"),
    Month = month(entranceTime),
    Month_date = paste(month(entranceTime), year(entranceTime)),
    Nstudyperiod = length(unique(Month_date)),
    NmonthsFrug = length(unique(Month_date[Month %in% c(6,7,8,9)])), # I suppose I have at least one obs for each month, thus this should count correctly (otherwise some month might be missing)
    NmonthsFol = length(unique(Month_date[Month %in% c(1,2,3,4,11,12)])) #Same, eems correct from verification
  ) %>% 
  group_by(
    x_y,
    id
  ) %>% 
  summarise(
    #Number of visits
    Nvisits = length(entranceTime),
    NvisitsFrug = sum(ifelse(Month %in% c(6,7,8,9), 1, 0)),
    NvisitsFol = sum(ifelse(Month %in% c(1,2,3,4,11,12), 1, 0)),
    
    #Frequency of visits (per month)
    FV = Nvisits/Nstudyperiod,
    FVFrug = NvisitsFrug/NmonthsFrug,
    FVFol = NvisitsFol/NmonthsFol,
    
    #Median time visits (in month)
    
    Nintervals = length(timeSinceLastVisit),
    NIntervalsFrug = length(timeSinceLastVisit[Month %in% c(6,7,8,9) & lag(Month) %in% c(6,7,8,9)]),
    NIntervalsFol = length(timeSinceLastVisit[Month %in% c(1,2,3,4,11,12) & lag(Month) %in% c(1,2,3,4,11,12)]),
    medianTimeInterval = median(timeSinceLastVisit, na.rm = TRUE)/365*12,
    medianTimeIntervalFrug = median(timeSinceLastVisit[Month %in% c(6,7,8,9) & lag(Month) %in% c(6,7,8,9)], na.rm = TRUE)/mean(c(30,31,31,30)),
    medianTimeIntervalFol = median(timeSinceLastVisit[Month %in% c(1,2,3,4,11,12) & lag(Month) %in% c(1,2,3,4,11,12)], na.rm = TRUE)/mean(c(31,28,31,30,30,31)),
    
    #Recursion index
    RV = FV * medianTimeInterval/log(2),
    RVFrug = FVFrug * medianTimeIntervalFrug/log(2),
    RVFol = FVFol * medianTimeIntervalFol/log(2)
  ) %>% 
  unique() %>% 
  left_join(
    ., visitedTrees_df %>% mutate(long_lat = paste(Longitude, Latitude, sep = "_")) %>% dplyr::select(Food_type2, long_lat), by = c("x_y" = "long_lat")
  ) %>% 
  mutate(
    sufficientSampling = "no",
    sufficientSampling = ifelse(
      Food_type2 == "leaves" &
        NIntervalsFol >= 4,
      "yes", 
      sufficientSampling
    ),
    sufficientSampling = ifelse(
      Food_type2 == "fruit" &
        NIntervalsFrug >= 4,
      "yes", 
      sufficientSampling
    ),
    sufficientSampling = ifelse(
      Food_type2 == "young_leaves" &
        NIntervalsFol >= 4,
      "yes", 
      sufficientSampling
    )
  )  %>% 
  unique() %>% 
  filter(
    !is.na(Food_type2) &
      sufficientSampling == "yes"
  ) # Filter if not food of interest, or insufficcient sampling

head(dfVisitsSummary) 
dfVisitsSummary$RV_seasonAdapted <- dfVisitsSummary$RVFol
dfVisitsSummary$RV_seasonAdapted[dfVisitsSummary$Food_type2 == "Fruit"] <- dfVisitsSummary$RVFrug

# plotrecursivityIndex <- ggplot(data = dfVisitsSummary,
#        aes(x = id, y = RV, fill = Food_type2)) +
#   stat_boxplot(geom ='errorbar', width = 0.6, size = 0.25) +
#   geom_boxplot(width = 0.6) +
#   stat_summary(fun.data = n_fun, geom = "text",
#                aes(group = Food_type2),
#                hjust = 0.5, position = position_dodge(0.6)) +
#   # annotate("text", x = 1:3, y = 1.05, label = labelsTest$PvalueText, size = 5) +
#   scale_fill_manual(name = "Food type", values=c("salmon", "darkgreen", "darkolivegreen3"), labels = c("Fruit", "(Mature) leaves", "Young leaves")) +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", size = 16),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12),
#         panel.grid.minor = element_line(colour = "grey93"),
#         panel.grid.major = element_line(colour = "grey93")) +
#   scale_x_discrete(labels = c("A", "B", "C")) +
#   xlab("Gorilla group") +
#   ylim(-0.1,1.1) +
#   ylab("Recursivity index")
# plotrecursivityIndex

n_fun <- function(x){
  return(data.frame(y = min(x) - 0.1,
                    label = length(x)))
}

linear_models <- dfVisitsSummary  %>% 
  ungroup() %>% 
  #group_by(id) %>% 
  do(anova(lm(RV_seasonAdapted ~ Food_type2, data = .), lm(RV_seasonAdapted ~ 1, data = .), test = "F")) %>% 
  filter(!is.na(Df)) 

labelsTest <- linear_models %>% 
  mutate(
    PvalueText = pvalueToText(`Pr(>F)`)
  ) 

# dataSegment <-
#   data.frame(
#     x = 1:3 - 0.3, 
#     xend = 1:3 + 0.3, 
#     y = rep(1.275, times = 3), 
#     yend = rep(1.275, times = 3)
#   )
# dataSegment$Food_type2 <- "leaves" #fake it because otherwise error

dataSegment2 <-
  data.frame(
    x = 1 - 0.3, 
    xend = 3 + 0.3, 
    y = rep(1.7, times = 1), 
    yend = rep(1.7, times = 1)
  )
dataSegment2$Food_type2 <- "leaves" #fake it because otherwise error


plotrecursivityIndex_seasonAdapted <- ggplot(data = dfVisitsSummary,
                                             aes(x = Food_type2,#id, 
                                                 y = RV_seasonAdapted, 
                                                 fill = Food_type2)) +
  stat_boxplot(geom = 'errorbar', width = 0.6, size = 0.25) +
  geom_boxplot(width = 0.6) +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group = Food_type2),
               hjust = 0.5, position = position_dodge(0.6)) +
  # annotate("text", x = 1:3, y = 1.3, label = labelsTest$PvalueText, size = 5) +
  # geom_segment(data = dataSegment, aes(x = x, xend = xend, y = y, yend = yend), 
  # color = "black", fill = "black") +
  annotate("text", x = 2, y = 1.8, label = labelsTest$PvalueText, size = 5) +
  geom_segment(data = dataSegment2, aes(x = x, xend = xend, y = y, yend = yend),
               color = "black", fill = "black") +
  scale_fill_manual(name = "Food type", values=c("salmon", "darkgreen", "darkolivegreen3"), labels = c("Fruit", "(Mature) leaves  ", "Young leaves")) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed")+
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  # scale_x_discrete(labels = c("A", "B", "C")) +
  # xlab("Gorilla group") +
  scale_x_discrete(labels = c("", "", "")) +
  xlab("Food type") +
  # ylim(-0.1,1.3) +
  ylab("Revisit regularity index") +
  scale_y_continuous(minor_breaks = seq(-10, 60, 0.1))
plotrecursivityIndex_seasonAdapted 

recursionPlot <- ggarrange(
  densityDistributionVisit_total,
  plotrecursivityIndex_seasonAdapted,
  labels = c("A)", "B)"),
  common.legend = TRUE,
  legend = "right",
  legend.grob = get_legend(plotrecursivityIndex_seasonAdapted)
)
recursionPlot 

# Analysis of efficiency --------------------------------------------------

numberSimulations = 100
radiusVision = 30

trk_rediscretized_bygroup_unnest <- trk_rediscretized_bygroup_unnest %>% 
  separate(
    col = id, 
    into = c("date", "group"),
    sep = "_",
    remove = FALSE
  )

library(sf)
mcpGroups <- lapply(sort(unique(trk_rediscretized_bygroup_unnest$group)), function(x){
  
  library(adehabitatHR)
  dataAsDf <- trk_rediscretized_bygroup_unnest %>% filter(group == x) %>% dplyr::select(x1_, y1_)
  dataAsDf.sp <- SpatialPoints(coords = dataAsDf)
  
  mcp_data <- adehabitatHR::mcp(dataAsDf.sp, percent = 95, unin = "m", unout = "m")
  
  # Intersect road network with MCP
  mcp_data.sf <- st_as_sf(mcp_data)
  return(mcp_data.sf)
}
)
# options(error=recover)
uniqueValues <- unique(trk_rediscretized_bygroup_unnest$id)

listComparisonEfficiency <- lapply(unique(trk_rediscretized_bygroup_unnest$id), function(x){
  print(which(uniqueValues == x))
  dataToWork <- trk_rediscretized_bygroup_unnest %>% 
    filter(id == x) %>% 
    dplyr::select(id, x1_, y1_, t1_, group) %>% 
    st_as_sf(., coords = c("x1_", "y1_")) %>% 
    mutate(
      month = month(t1_),
      season = NA,
      season = ifelse(month %in% c(6,7,8,9), "fruit", season),
      season = ifelse(month %in% c(11,12,1,2,3,4), "leaves", season),
      geometry_lagged = lag(geometry),
      line =
        st_sfc(purrr::map2(
          .x = geometry,
          .y = geometry_lagged,
          .f = ~{if(!st_is_empty(.y)){st_union(c(.x, .y)) %>% st_cast("LINESTRING")}}
        ))
    )
  
  if(is.na(dataToWork$season[1]) | nrow(dataToWork) < 10){#Not good season or small data = error in tracking
    
  }else{
    #DPL
    DPL <- sum(st_length(dataToWork$line))
    
    #Number of trees encounter
    treeLoc.sf <- treeLoc.df %>% 
      mutate(
        Food = ifelse(Food == "young_leaves", "leaves", Food)
      ) %>% 
      st_as_sf(coords = c("Longitude", "Latitude"))
    
    NTreesEncountered = nrow(st_intersection(
      treeLoc.sf,
      st_buffer(
        dataToWork,
        dist = radiusVision
      )
    )
    )
    
    #Number of trees theoretically expected
    randomAngle = runif(numberSimulations, 0, 360) * pi / 180
    dataBallistic = data.frame(
      startPointx = rep(trk_rediscretized_bygroup_unnest$x1_[1], times = numberSimulations),
      startPointy = rep(trk_rediscretized_bygroup_unnest$y1_[1], times = numberSimulations),
      endPointx = trk_rediscretized_bygroup_unnest$x1_[1] + DPL * cos(randomAngle),
      endPointy = trk_rediscretized_bygroup_unnest$y1_[1] + DPL * sin(randomAngle)
    ) %>% 
      mutate(id = paste0(endPointx, endPointy))
    
    #Remove if points not in mcp
    dataBallistic_sf <- st_as_sf(dataBallistic, coords = c("endPointx", "endPointy"))
    mcpToConsider <- mcpGroups[[as.numeric(as.factor(unique(dataToWork$group)))[1]]]
    intersectionDf <- st_intersection(dataBallistic_sf, mcpToConsider)
    count = 0
    continueProcessing = "YES"
    while(nrow(intersectionDf) < 100){
      randomAngle = runif(100 - nrow(intersectionDf), 0, 360) * pi / 180
      toAdd = data.frame(
        startPointx = rep(trk_rediscretized_bygroup_unnest$x1_[1], times = 100 - nrow(intersectionDf)),
        startPointy = rep(trk_rediscretized_bygroup_unnest$y1_[1], times = 100 - nrow(intersectionDf)),
        endPointx = trk_rediscretized_bygroup_unnest$x1_[1] + DPL * cos(randomAngle),
        endPointy = trk_rediscretized_bygroup_unnest$y1_[1] + DPL * sin(randomAngle)
      ) %>% 
        mutate(id = paste0(endPointx, endPointy))
      
      dataBallistic <- rbind(dataBallistic %>% filter(id %in% intersectionDf$id), toAdd)
      dataBallistic_sf <- st_as_sf(dataBallistic, coords = c("endPointx", "endPointy"))
      intersectionDf <- st_intersection(dataBallistic_sf, mcpToConsider)
      count = count + 1
      if(count > 30){
        continueProcessing = "NO"
        break
      }
      # print(count)
    }
    
    if(continueProcessing == "YES"){#When the DPL is too long, no such ballistic walk can be computed, so we discard these cases
      NTree_simulations_v <- sapply(1:numberSimulations,
                                    function(i){
                                      dataBallistic.sf1 <- st_as_sf(dataBallistic[i,], coords = c("startPointx", "startPointy"))
                                      dataBallistic.sf2 <- st_as_sf(dataBallistic[i,], coords = c("endPointx", "endPointy"))
                                      dataBallistic.sf <- st_union(dataBallistic.sf1, dataBallistic.sf2) %>% st_cast("LINESTRING") %>%  st_buffer(., dist = radiusVision)
                                      NTree <- nrow(
                                        st_intersection(
                                          treeLoc.sf,
                                          dataBallistic.sf
                                        )
                                      )
                                      return(ifelse(!is.null(NTree), NTree, 0))
                                    }
      )
      NTree_simulations <- mean(NTree_simulations_v, na.rm=TRUE)
      return(c(
        ifelse(dataToWork$season[1] == "fruit", "FRUGIVORY", "FOLIVORY"),
        DPL,
        ifelse(is.null(NTreesEncountered), 0, NTreesEncountered),
        ifelse(is.null(NTree_simulations), 0, NTree_simulations)
      ))
    }
  }
})

comparisonEfficiency_df <- do.call("rbind", listComparisonEfficiency) %>%  as.data.frame()
colnames(comparisonEfficiency_df) <- c("Season", "DPL", "NTreesVisited_observed", "NTreesVisited_simulated")
comparisonEfficiency_df$DPL <- as.numeric(comparisonEfficiency_df$DPL)
comparisonEfficiency_df$NTreesVisited_observed <- as.numeric(comparisonEfficiency_df$NTreesVisited_observed)
comparisonEfficiency_df$NTreesVisited_simulated <- as.numeric(comparisonEfficiency_df$NTreesVisited_simulated)

#save df
saveRDS(comparisonEfficiency_df, "Renvironment/comparisonEfficiency_df.rds")

n_fun <- function(x){
  return(data.frame(y = min(x) - 5,
                    label = length(x)/2))
}

labelsTest <- lapply(c("FRUGIVORY", "FOLIVORY"), function(x){
  data <-  comparisonEfficiency_df %>% ungroup() %>% filter(Season == x)
  ttest <- t.test(as.numeric(data$NTreesVisited_observed), as.numeric(data$NTreesVisited_simulated), paired = TRUE)
  return(pvalueToText(ttest$p.value))
})

comparisonEfficiencyTidy_df <-
  comparisonEfficiency_df %>% 
  pivot_longer(
    cols = c("NTreesVisited_observed", "NTreesVisited_simulated"),
    names_to = "Type",
    values_to = c("NvisitedTrees")
  ) %>% 
  mutate(
    Season_type = paste0(Season, Type),
    treeByKm = NvisitedTrees * 1000 / DPL,
    alpha = ifelse(Type == "NTreesVisited_simulated", 0.75, 1)
  )

plotEfficiency <- ggplot(data = comparisonEfficiencyTidy_df,
                         aes(x = Season, 
                             y = treeByKm, 
                             fill = Season,
                             alpha = Type)) +
  stat_boxplot(geom ='errorbar', width = 0.6, size = 0.25) +
  geom_boxplot(width = 0.6) +
  stat_summary(fun.data = n_fun, geom = "text",
               aes(group = Season),
               hjust = 0.5, position = position_dodge(0)) +
  stat_summary(fun = "mean", aes(group = Season_type), pch = 21, fill = "white", size = 1,
               position = position_dodge(0.6), alpha = 1) +
  annotate("text", x = 1:2, y = ceiling(max(comparisonEfficiencyTidy_df$treeByKm)) + 5, 
           label = labelsTest, size = 5) +
  scale_fill_manual(values = c("salmon", "darkgreen")) +
  scale_alpha_manual(values = c(1, 0.25)) +
  guides(fill = "none", alpha = "none") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93")) +
  scale_x_discrete(labels = c("Obs.  Sim.\nHigh-frugivory", "Obs.  Sim.\nLow-frugivory")) +
  xlab("Season") +
  # ylim(-5,ceiling(max(comparisonEfficiencyTidy_df$treeByKm)) + 5) +
  ylab("Number of visited trees\nper travelled km")  +
  scale_y_continuous(minor_breaks = seq(-100, 100, 2))
plotEfficiency

library(ggpubr)
LinearitySpeedEfficiencyDisplacementPlot <- ggarrange(
  ggarrange(Linearityplot,
            Speedplot,
            labels = c("A)", "B)"),
            common.legend = TRUE,
            ncol = 2,
            nrow =1
  ),
  ggarrange(displacementAndDPLplot, labels = c("C)")),
  ggarrange(
    Displacementplot,
    plotEfficiency,
    labels = c("D)", "E)"),
    legend = "none",
    ncol = 2,
    nrow =1
  ),
  nrow = 3,
  heights = c(0.3, 0.4, 0.3)
)
LinearitySpeedEfficiencyDisplacementPlot


# Save plots --------------------------------------------------------------

ggsave("Graphics/DPL.png", DPLPlot, dpi = 800, width = 12, height = 7)
ggsave("Graphics/Movement.png", LinearitySpeedEfficiencyDisplacementPlot, dpi = 800, width = 12, height = 16)
ggsave("Graphics/Recursion.png", recursionPlot, dpi = 800, width = 12, height = 7)

# Samples -----------------------------------------------------------------

# Complete obs
merged %>%
  group_by(Group) %>% 
  summarise(
    Nobscompletedays = length(unique(Ymd))
  )

# Those with sufficient monitoring during the day
DPL.df  %>% 
  filter(Loc == "CAR") %>% 
  group_by(Group) %>% 
  summarise(
    Nobscompletedays = length(Date)
  )

# Behavioural data
data_behaviour %>%
  group_by(Group) %>% 
  summarise(
    Nobscompletedays = length(unique(Date))
  )

#Provide sampling size: individual duration to Shelly on the 18/11/2022
data_behaviour <- read_delim("Data/Data_behaviour.txt", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE, guess_max=10000)
head(data_behaviour)
saveObsTime <- data_behaviour %>%
  group_by(Focal, Group, Date) %>% 
  summarise(
    obsTime = difftime(Time[length(Time)], Time[1], units = "min")
  ) %>% 
  group_by(Focal, Group) %>% 
  summarise(
    obsTimeHours = as.numeric(sum(obsTime))/60
  )

write.table(saveObsTime, "Data/Observation_duration_focal.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Plotting ID/RD --------------------------------------

library(raster)
library(sf)
library(stars)

groupToPlot <- c(CAR1, CAR3, CAR2)
raster_l <- vector("list", length = length(groupToPlot))

raster_to_df <- function(x){
  test_spdf <- as(x, "SpatialPixelsDataFrame")
  test_df <- as.data.frame(test_spdf)
  colnames(test_df) <- c("value", "x", "y")
  test_df$value[test_df$value == -100] <- NA
  return(test_df)
}

for(i in 1:length(groupToPlot)){
  toPlotID_frugivory <- raster(paste("Data/CID_", groupToPlot[i], "High-frugivory.asc", sep = ""))
  proj4string(toPlotID_frugivory) <- CRS(paste("+proj=utm","+zone=33","+ellps=WGS84", "+datum=WGS84", "+units=m", "+towgs84:0,0,0", sep=" "))
  
  toPlotID_folivory <- raster(paste("Data/CID_", groupToPlot[i], "Low-frugivory.asc", sep = ""))
  proj4string(toPlotID_folivory) <- CRS(paste("+proj=utm","+zone=33","+ellps=WGS84", "+datum=WGS84", "+units=m", "+towgs84:0,0,0", sep=" "))
  
  toPlotRD_frugivory <- raster(paste("Data/CRD_", groupToPlot[i], "High-frugivory.asc", sep = ""))
  proj4string(toPlotRD_frugivory) <- CRS(paste("+proj=utm","+zone=33","+ellps=WGS84", "+datum=WGS84", "+units=m", "+towgs84:0,0,0", sep=" "))
  
  toPlotRD_folivory <- raster(paste("Data/CRD_", groupToPlot[i], "Low-frugivory.asc", sep = ""))
  proj4string(toPlotRD_folivory) <- CRS(paste("+proj=utm","+zone=33","+ellps=WGS84", "+datum=WGS84", "+units=m", "+towgs84:0,0,0", sep=" "))
  
  maxX <- max(extent(toPlotID_frugivory)@xmax,
              extent(toPlotID_folivory)@xmax,
              extent(toPlotRD_frugivory)@xmax,
              extent(toPlotRD_folivory)@xmax)
  
  minX <- min(extent(toPlotID_frugivory)@xmin,
              extent(toPlotID_folivory)@xmin,
              extent(toPlotRD_frugivory)@xmin,
              extent(toPlotRD_folivory)@xmin)
  
  maxY <- max(extent(toPlotID_frugivory)@ymax,
              extent(toPlotID_folivory)@ymax,
              extent(toPlotRD_frugivory)@ymax,
              extent(toPlotRD_folivory)@ymax)
  
  minY <- min(extent(toPlotID_frugivory)@ymin,
              extent(toPlotID_folivory)@ymin,
              extent(toPlotRD_frugivory)@ymin,
              extent(toPlotRD_folivory)@ymin)
  
  #Put raster to same resolution and size
  toPlotID_folivory <- extend(toPlotID_folivory, extent(c(minX, maxX, minY, maxY)), value = NA)
  toPlotID_frugivory <- extend(toPlotID_frugivory, extent(c(minX, maxX, minY, maxY)), value = NA)
  toPlotRD_frugivory <- extend(toPlotRD_frugivory, extent(c(minX, maxX, minY, maxY)), value = NA)
  toPlotRD_folivory <- extend(toPlotRD_folivory, extent(c(minX, maxX, minY, maxY)), value = NA)
  
  #Change NA in -100 that will later be backtransformed to NA
  toPlotID_frugivory[is.na(toPlotID_frugivory)] <- -100
  toPlotID_folivory[is.na(toPlotID_folivory)] <- -100
  toPlotRD_frugivory[is.na(toPlotRD_frugivory)] <- -100
  toPlotRD_folivory[is.na(toPlotRD_folivory)] <- -100
  
  # toPlot_sf <- toPlot %>%  st_as_stars() %>%  st_as_sf
  raster_l[[i]] <-
    lapply(list(toPlotID_frugivory, toPlotID_folivory, toPlotRD_frugivory, toPlotRD_folivory), raster_to_df)
}

##Plot tools
library(ggplot2)
library(ragg)#To have better display
library(ggpattern)
library(ggrepel)
library(grid)
##Spatial
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

#Function to plot from raster
plotFolivory <- function(test_df, titley = ""){
  ggplot() +  
    geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8, na.rm = TRUE) + 
    scale_fill_gradient(low = "gold", high = "darkgreen", na.value = NA) +
    coord_equal() +
    theme_bw() +
    theme(axis.title = element_text(face = "bold", size = 16),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"),
          panel.grid.minor = element_line(colour = "grey93"),
          panel.grid.major = element_line(colour = "grey93")) +
    annotation_scale(location = "bl",
                     width_hint = 0.25,
                     height = unit(0.2, "cm")) +
    annotation_north_arrow(
      height = unit(1, "cm"),
      width = unit(1, "cm"),
      location = "tr",
      which_north = "true",
      pad_x = unit(0.05, "in"),
      pad_y = unit(0.1, "in"),
      style = north_arrow_fancy_orienteering
    ) +
    xlim(c(min(test_df[,2], na.rm = TRUE), max(test_df[,2], na.rm = TRUE))) +
    ylim(c(min(test_df[,3], na.rm = TRUE), max(test_df[,3], na.rm = TRUE))) +
    xlab("") +
    ylab(titley) +
    guides(fill = "none")
}

plotFrugivory <- function(test_df, titley = ""){
  ggplot() +  
    geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8, na.rm = TRUE) + 
    scale_fill_gradient(low = "gold", high = "salmon", na.value = NA) +
    coord_equal() +
    theme_bw() +
    theme(axis.title = element_text(face = "bold", size = 16),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"),
          panel.grid.minor = element_line(colour = "grey93"),
          panel.grid.major = element_line(colour = "grey93")) +
    annotation_scale(location = "bl",
                     width_hint = 0.25,
                     height = unit(0.2, "cm")) +
    annotation_north_arrow(
      location = "tr",
      height = unit(1, "cm"),
      width = unit(1, "cm"),
      which_north = "true",
      pad_x = unit(0.05, "in"),
      pad_y = unit(0.1, "in"),
      style = north_arrow_fancy_orienteering
    ) +
    xlim(c(min(test_df[,2], na.rm = TRUE), max(test_df[,2], na.rm = TRUE))) +
    ylim(c(min(test_df[,3], na.rm = TRUE), max(test_df[,3], na.rm = TRUE))) +
    xlab("") +
    ylab(titley) +
    guides(fill = "none")
}

# Create plot for each group and season
counter = 0
plot_l <- lapply(raster_l, function(x){
  counter <<- counter + 1
  return(list(
    plotFrugivory(x[[1]], titley = paste0("Group ", counter)),
    plotFolivory(x[[2]]),
    plotFrugivory(x[[3]]),
    plotFolivory(x[[4]])
  )
  )
})

#Rorganise plot to gather group and species
library(ggpubr)
arrange1 <- ggarrange(plot_l[[1]][[1]], plot_l[[1]][[2]], 
                      labels = c("Low-frugivory", "High-frugivory"),#"Low-frug.", "High-frug.", "Low-frug.", "High-frug."), 
                      label.x = c(0.5, 0.5), label.y = c(0.975, 0.975), hjust = 0.25) 
arrange1
arrange2 <- ggarrange(plot_l[[1]][[3]], plot_l[[1]][[4]], 
                      labels = c("Low-frugivory", "High-frugivory"),#"Low-frug.", "High-frug.", "Low-frug.", "High-frug."), 
                      label.x = c(0.5, 0.5), label.y = c(0.975, 0.975), hjust = 0.25) 
arrange2
arrange3  <- ggarrange(arrange1, arrange2, 
                       labels = c("ID", "RD"),#"Low-frug.", "High-frug.", "Low-frug.", "High-frug."), 
                       label.x = c(0.5, 0.5), label.y = c(1, 1)) 
arrange3

arrange4 <- ggarrange(plot_l[[2]][[1]], plot_l[[2]][[2]], plot_l[[2]][[3]], plot_l[[2]][[4]], nrow = 1) 
arrange5 <- ggarrange(plot_l[[3]][[1]], plot_l[[3]][[2]], plot_l[[3]][[3]], plot_l[[3]][[4]], nrow = 1) 

plotIDRD <- ggarrange(arrange3, arrange4, arrange5, nrow = 3, heights = c(0.35,0.35, 0.3))

ggsave("Graphics/IDRD.png", plotIDRD, dpi = 800, width = 12, height = 12)

save.image("Renvironment/analysisBook.RData")