######             CCFPS TECHNICAL REPORT 1 ANALYSIS            #######
library(FSA)
library(nlstools)
library(lubridate)
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(purrr)
library(tidyverse)
library(rlist)
library(Hmisc)
library(ggplot2)
library(survival)
library(ggfortify)
library(gam)


# Functions ##############################################################
#function to turn a date into a week of the study
date_to_studyweek <- function (date) {
  studyweek <- floor(as.numeric(date - as_date("2013-05-03"))/7)+1
  return(studyweek)
}

#function to find the start date of a study week
studyweek_startdate <- function (Week) {
  startdate <- as_date("2013-05-03")+weeks(Week-1)
  return(startdate)
}

#Function to find the end date of a study week
studyweek_enddate <- function (Week) {
  enddate <- as_date("2013-05-09")+weeks(Week-1)
  return(enddate)
}

# Create helper function to split the interval of the CRE into a given unit size (ie weeks)
# Remember that the original weekly bin was only for detections, we are now
# assuming that a fish is remaining in a location and need to code each week
# it is at-large as such
split_interval_weeks <- function(start, end, unit) {
  breaks <- seq(floor_date(start, "day"), ceiling_date(end, "day"), by = unit)
  timeline <- c(start, breaks[breaks > start & breaks < end], end)
  tibble(.start = head(timeline, -1), .end = tail(timeline, -1))
}

# Format Data #############################################################
# This script picks up wfile.path the BKM_QAQC script ends

## Read in the detections and hydrophone data =============================
all_bkm <- readRDS(file.path("1. Data","Inputs","all_bkm.rds"))

hydrophones <- read.csv(file.path("1. Data","Inputs","HydrophoneSites.csv"), header = T, 
                        stringsAsFactors = T)
hydrophones$SiteCode <- factor(hydrophones$SiteCode, 
                               levels = c("IC3","IC2","IC1","RGD1","RGU1",
                                          "ORS1","WC1","ORN1","WC2","ORN2",
                                          "GL1","CVP1","WC3","ORS3","CLRS"))
Hydrophones <- hydrophones[,c(1,5)]
Hydrophones <- Hydrophones[!is.na(Hydrophones$SiteCode),]

## Format the data and turn into a weekly detection history ===============
# Extract the start of each bookmark into a single table
all_bkm_week <- all_bkm %>%
  mutate(Week = date_to_studyweek(as_date(CorrLocal_DateTime)),
         WeekEnd = date_to_studyweek(as_date(CorrLocal_DateTime)))

# Create a vector of all study weeks, all tags, and all hydrophones
Weeks <- seq(1,257, by = 1)
Tags <- all_bkm_week %>%
  distinct(TagID)

# Cross all the vectors and name the columns correctly
DateHydTagMatrix<- crossing(Weeks, hydrophones, Tags)
DateHydTagMatrix <- DateHydTagMatrix %>%
  select(Weeks, SiteCode, HydrophoneID, TagID)
names(DateHydTagMatrix) <- c("Week","SiteCode","HyfonNo", "TagID")

### Hydrophone Checks -----------------------------------------------------
#Determine whether the hydrophone heard its assigned beacon tag that week
beacons <- read.csv(file.path("1. Data","Inputs","BeaconTags.csv"), header = T, stringsAsFactors = T)
beacons$species <- "Beacon"

HydBeaconMatrix <- filter(DateHydTagMatrix, 
                          DateHydTagMatrix$TagID %in% beacons$TagID)
HydBeaconMatrix <- HydBeaconMatrix %>%
  left_join(.,beacons, by = "TagID")
HydBeaconMatrix <- HydBeaconMatrix %>%
  select(Week, SiteCode, HyfonNo, TagID)
HydBeaconBin <- HydBeaconMatrix %>%
  left_join(.,select(all_bkm_week, TagID, HyfonNo, Week, MarkType), 
            by = c("Week","TagID","HyfonNo")
            ) %>%
  mutate(Detected = !is.na(MarkType)) %>%
  distinct(Week,SiteCode,HyfonNo,TagID,Detected) %>%
  mutate(AssignedTag = beacons[match(SiteCode,beacons$Location),2]) %>%
  mutate(BeaconDetected = ifelse(AssignedTag == TagID & Detected == T,T,F)) %>%
  filter(BeaconDetected == TRUE) %>%
  distinct(Week,SiteCode,BeaconDetected) %>%
  arrange(SiteCode, Week)

### Build the weekly detection histories ----------------------------------
# Create a week-binned tag record
all_bkm_week <- all_bkm_week %>%
  mutate(SiteCode = hydrophones[match(HyfonNo, hydrophones$HydrophoneID),1])

TagBKM_Bin <- DateHydTagMatrix %>%
  left_join(.,
            select(all_bkm_week, TagID, HyfonNo, Week, 
                   CorrLocal_DateTime, MarkType), 
            by = c("TagID", "HyfonNo", "Week")
            ) %>%
  mutate(Detected = !is.na(MarkType))%>%
  filter(Detected = TRUE) %>%
  filter(!is.na(CorrLocal_DateTime))

TagBKM_Bin <- TagBKM_Bin %>%
  group_by(Week, SiteCode, TagID) %>%
  summarise(first_detect = min(CorrLocal_DateTime,na.rm = TRUE),
            last_detect = max(CorrLocal_DateTime, na.rm = TRUE)) %>%
  arrange(TagID, Week, last_detect)

# Create a column for whether the beacon tag was heard that week
TagBKM_Bin <- TagBKM_Bin %>%
  left_join(.,select(HydBeaconBin, Week, SiteCode, BeaconDetected), 
            by = c("Week","SiteCode"))

# Create a column for whether any tag was heard on that week
# First create a matrix for Date, SiteCode, 
AnyDetected <- TagBKM_Bin %>%
  select(Week, SiteCode, TagID) %>%
  distinct(Week, SiteCode) %>%
  mutate(Detected = TRUE)

# Then join that matrix to the Full Date-Binned Matrix
TagBKM_Bin <- TagBKM_Bin %>%
  left_join(.,select(AnyDetected, Week, SiteCode, AnyDetected=Detected), 
            by = c("Week","SiteCode"))

# Both of the above joins creat "NA" for Non-detection, so we code NA as FALSE
TagBKM_Bin$BeaconDetected[is.na(TagBKM_Bin$BeaconDetected)] <- FALSE
TagBKM_Bin$AnyDetected[is.na(TagBKM_Bin$AnyDetected)] <- FALSE

### Add in the Release Data from the CSV
release_dates <- read.csv(file.path("1. Data","Inputs","Release_Dates.csv"), 
                          stringsAsFactors = F, header = T)
release_dat <- read.csv(file.path("1. Data","Inputs","ReleaseData.csv"), 
                        stringsAsFactors = F, header = T)
release <- release_dates %>%
  left_join(.,release_dat, by = "TagID")
release$SiteCode <- "Release"
release$Week <- date_to_studyweek(parse_date(release$Date, format = "%m/%d/%Y"))
release$Detected <- TRUE
release$last_detect = lubridate::mdy(release$Date)

TagBKM_Bin <- TagBKM_Bin %>%
  left_join(.,select(release, TagID, Species), by = "TagID")
TagBKM_Bin <- bind_rows(TagBKM_Bin, 
                        select( release, Species,TagID, SiteCode, last_detect,
                                Length_cm, WT_lbs, Week))

TagBKM_Bin$SiteCode <- ordered(TagBKM_Bin$SiteCode, 
                               levels = c("IC3","IC2","IC1","Release", "RGD1",
                                          "RGU1","ORS1","WC1","ORN1","WC2",
                                          "ORN2","GL1","CVP1","WC3","ORS3",
                                          "CLRS"))

TagBKM_Bin <- TagBKM_Bin %>%
  group_by(TagID) %>%
  arrange(TagID, Week, last_detect) %>%
  mutate(Species = ifelse(TagID > 10000, "Beacon",Species)) %>%
  select(Species, TagID, Week, last_detect, SiteCode, 
         BeaconDetected, AnyDetected)

# Clean up
rm(list = c("release_dat","release_dates","HydBeaconMatrix",
            "HydBeaconBin", "DateHydTagMatrix"))

### Summarise Tag Release Data --------------------------------------------
release_summary <- release %>%
  group_by(Species, Week) %>%
  tally()

releases <- ggplot()+
  geom_bar(data = release_summary,aes(x = studyweek_startdate(Week), y = n), 
           stat = "identity", width = 5)+
  scale_x_date(breaks = "3 months", expand = c(0,0))+
  facet_grid(Species~., scales = "free_y")+
  theme(
    axis.text = element_text(angle = -45, vjust = 0.5, hjust = 0)
  )+
  labs(x = "Date", y = "Total Number of Released Individuals")

pdf(file = file.path("1. Data","Figures","ReleaseSummary.pdf"), width = 15, height = 5)
releases
dev.off()


# Check a tag
ggplot()+
  geom_point(data = TagBKM_Bin[TagBKM_Bin$TagID %in% c(5004.24),], 
             aes(x = Week, y = SiteCode, group = TagID))+
  geom_step(data = TagBKM_Bin[TagBKM_Bin$TagID %in% c(5004.24),], 
            aes(x = Week, y = SiteCode, group = TagID))+
  scale_y_discrete()+
  theme(panel.spacing = unit(0, "lines"))+
  facet_grid(TagID ~.)



# Begin Analysis ##########################################################
## Shed Tag Review ========================================================
HydVisit <- TagBKM_Bin %>%
  mutate(Location = case_when(
    SiteCode %in% c("IC1","IC2","IC3") ~ "Intake", 
    SiteCode %in% c("RGD1","RGU1") ~ "Radial Gates",
    SiteCode %in% c("ORS1","WC1") ~ "Intersection",
    SiteCode %in% c("WC1","WC2","WC3") ~ "West Canal",
    SiteCode %in% c("ORN1","ORN2") ~ "Old River North",
    SiteCode %in% c("Release") ~ "Release",
    SiteCode %in% c("ORS3") ~ "Old River South",
    SiteCode %in% c("CVP1") ~ "CVP Intake",
    SiteCode %in% c("GL1") ~ "Grant Line Canal"))

Site_CRE <- HydVisit %>% 
  group_by(Species,TagID)%>%
  arrange(TagID,Week)%>%
  filter(SiteCode != "Release") %>%
  mutate(
    SiteChange = coalesce(Location != lag(Location), FALSE) %>% cumsum()
  ) %>% 
  group_by(Species,TagID, SiteChange, Location) %>% 
  summarise(
    StartWeek = min(Week),
    EndWeek = max(Week)
  ) %>% 
  ungroup()%>%
  arrange(TagID, StartWeek)

Site_CRE <- Site_CRE %>%
  mutate(duration = EndWeek - StartWeek + 1)

Site_CRE_Quantiles <- Site_CRE %>%
  filter(Species != "Beacon") %>%
  select(Species, duration) %>%
  group_by(Species) %>%
  mutate( qnt_0   = quantile(duration, probs= 0),
          qnt_10   = quantile(duration, probs= 0.1),
          qnt_25  = quantile(duration, probs= 0.25),
          qnt_50  = quantile(duration, probs= 0.5),
          qnt_75  = quantile(duration, probs= 0.75),
          qnt_90 = quantile(duration, probs= 0.9),
          qnt_95 = quantile(duration, probs= 0.95),
          qnt_99 = quantile(duration, probs= 0.99),
          qnt_100 = quantile(duration, probs= 1),
          mean = mean(duration),
          sd = sd(duration)
  ) %>%
  distinct(qnt_0, qnt_10, qnt_25, qnt_50, qnt_75, qnt_90, 
           qnt_95, qnt_99, qnt_100 ,mean ,sd)

ShedTags_plot <- ggplot()+
  geom_point(data = Site_CRE[Site_CRE$Species != "Beacon",], 
             aes(x = Species, y = duration), 
             alpha = 0.25, shape = 1, color = "blue")+
  geom_boxplot(data = Site_CRE_Quantiles, 
               aes(x = Species, ymin = qnt_10, lower = qnt_25, 
                   middle = qnt_50, upper = qnt_75, ymax=qnt_95), 
               stat = "identity", lwd = 0.5)+
  geom_point(data = Site_CRE %>%
               filter(Species != "Beacon") %>%
               left_join(Site_CRE_Quantiles) %>%
               filter(duration > qnt_95), 
             aes(x = Species, y = duration), 
             alpha = 0.25, shape = 4, color = "red")+
  scale_y_log10(name = "Uninterrupted Detection Periods (weeks)")+
  ggtitle("Analysis of Uninterrupted Detection Periods at a Single Receiver")+
  theme_classic()
  
pdf(file = file.path("1. Data", "Figures","ShedTagReview.pdf"),
    width = 4, height = 4.5)
ShedTags_plot
dev.off()

ShedTags <- Site_CRE %>%
  left_join(.,Site_CRE_Quantiles,by = "Species")%>%
  mutate(Shed = ifelse(duration > qnt_95,TRUE,FALSE)) %>%
  filter(Shed == T) %>%
  group_by(TagID)%>%
  summarise(ShedWeek = min(StartWeek)+1)

## Estimating TagLife =====================================================
# Tag life based on maximum Tag Life for specific tag type
TagTypes <- read.csv(file.path("1. Data","Inputs","TagTypes.csv"), header = T, 
                     stringsAsFactors = F)
TagLife <- TagBKM_Bin %>% 
  left_join(., TagTypes, by = "TagID") %>%
  mutate(Date = studyweek_startdate(Week)) %>%
  group_by(TagID, TagType) %>%
  summarise(TagLife = as.numeric(max(Date)-min(Date))) %>%
  group_by(TagType) %>%
  summarise(TagLife = max(TagLife)) #Which summary should we use? 

TagFail <- release %>%
  select(TagID, Species, Date) %>%
  left_join(.,TagTypes, "TagID") %>%
  mutate(Date = parse_date(Date, format = "%m/%d/%Y")) %>%
  mutate(TagType = ifelse(is.na(TagType),"LG",TagType)) %>%
  left_join(.,TagLife, "TagType") %>%
  mutate(faildate = Date + days(TagLife)) %>%
  mutate(Week = date_to_studyweek(faildate), Detected = T, 
         SiteCode = "Tag Failure") %>%
  select(Week, SiteCode, TagID, Detected, Species)

TagBKM_Bin <- bind_rows(TagBKM_Bin,TagFail)
TagBKM_Bin$SiteCode <- ordered(TagBKM_Bin$SiteCode, 
                               levels = c("Release", "IC3","IC2","IC1","RGD1",
                                          "RGU1","ORS1","WC1","ORN1","WC2",
                                          "ORN2","GL1","CVP1","WC3","ORS3",
                                          "CLRS","Tag Failure"))
TagBKM_Bin <- TagBKM_Bin %>%
  left_join(TagTypes, by = "TagID")

# Check some tags 
ggplot(TagBKM_Bin[TagBKM_Bin$TagID %in% c(5956.24) & 
                    TagBKM_Bin$SiteCode != "Tag Failure",])+
  geom_point(aes(x = Week, y = SiteCode, group = TagID))+
  geom_step(aes(x = Week, y = SiteCode, group = TagID))+
  scale_y_discrete()+
  theme(panel.spacing = unit(0, "lines"))

## Redo weekly detection histories including tag failures =================
# Determine how many and which sites a tag visited on a single week
# Redo Hydvisit to include Tag Failure
HydVisit <- TagBKM_Bin %>%
  mutate(Location = case_when(
    SiteCode %in% c("IC1","IC2","IC3") ~ "Intake", 
    SiteCode %in% c("RGD1","RGU1") ~ "Radial Gates",
    SiteCode %in% c("ORS1","WC1") ~ "Intersection",
    SiteCode %in% c("WC1","WC2","WC3") ~ "West Canal",
    SiteCode %in% c("ORN1","ORN2") ~ "Old River North",
    SiteCode %in% c("Release") ~ "Release",
    SiteCode %in% c("ORS3") ~ "Old River South",
    SiteCode %in% c("CVP1") ~ "CVP Intake",
    SiteCode %in% c("GL1") ~ "Grant Line Canal",
    SiteCode %in% c("Tag Failure") ~ "Tag Failure"))

WeeklySiteVisit <- TagBKM_Bin %>%
  filter(Species != "Beacon") %>%
  group_by(Species, Week, TagID) %>%
  summarise(SiteVisits = list(unique(as.character(SiteCode))),
            StartSite = as.character(SiteCode)[1],
            EndSite = as.character(SiteCode)[length(as.character(SiteCode))]) %>%
  mutate(NumSites = lengths(SiteVisits)) %>%
  ungroup()

inside <- c("Release","IC1", "IC2", "IC3", "RGD1")
outside <- c("ORS1", "WC1", "WC2", "WC3", "ORN1", 
             "ORN2", "ORS3", "GL1", "CVP1", "CLRS","RGU1")

TotalSiteVisitSummary <- WeeklySiteVisit %>%
  group_by(Species,TagID) %>%
  mutate(TAL = max(Week)-min(Week)+1) %>%
  filter(!(grepl("Tag Failure",SiteVisits))) %>%
  ungroup() %>%
  group_by(Species, TagID, TAL) %>%
  summarise(
    inside_sites = sum(
      unlist(
        lapply(SiteVisits,function(x){unlist(x, 
                                  recursive = TRUE) %in% inside}))),
    outside_sites = sum(
      unlist(
        lapply(SiteVisits,function(x){unlist(x, 
                                             recursive = TRUE) %in% outside}))),
    numSites = sum(NumSites)-1) %>%
  ungroup()

no_detect_tags <- TotalSiteVisitSummary %>%
  filter(numSites == 0) %>%
  pull(TagID)

### Code Location as either Inside, Outside, or Transiting ----------------

WeeklySiteVisit$Location <- lapply(WeeklySiteVisit$SiteVisits, 
                             function(x) case_when(
                               # Weeks wfile.path ALL detections are "Inside"
                               all(x %in% inside) ~ "INSIDE",
                               # Weeks wfile.path ALL detections are "Outside"
                               all(x %in% outside) ~ "OUTSIDE",
                               # Weeks with detections at the Intake Canal 
                               # AND beyond the radial gates are "FULL TRANSIT"
                               any(x %in% c("Release", "IC1", "IC2", "IC3")) &
                                     any(x %in% c("ORS1", "WC1", "WC2", "WC3", 
                                              "ORN1", "ORN2", "ORS3", "GL1", 
                                              "CVP1", "CLRS")) ~ "FULL TRANSIT",
                               # Any week wfile.path ALL detections are "Tag Failure"
                               all(x == "Tag Failure") ~ "Tag Failure",
                               # Anything else is a partial transit (i.e. may
                               # or may not have actually crossed the gates)
                               TRUE ~ "PARTIAL TRANSIT"))

# Determine if a partial transit can be coded as inside or outside based on its
# final detection that week
WeeklySiteVisit <- WeeklySiteVisit %>%
  mutate(
    Location = ifelse(
      #As long as the last detection that week was not "RGD1" or "RGU1"...
      Location == "PARTIAL TRANSIT" & 
        !(EndSite %in% c("RGD1","RGU1")), 
      case_when(
        # If it was outside, code location to outside
        EndSite %in% outside ~ "OUTSIDE",
        # If it was inside, code location to inside
        EndSite %in% inside ~ "INSIDE",
        # Backup in case I made a mistake, if both start and end are radial
        # gate detects, it's unresolved
        StartSite == "RGU1" & EndSite == "RGD1" ~ "UNRESOLVED TRANSIT",
        StartSite == "RGD1" & EndSite == "RGU1" ~ "UNRESOLVED TRANSIT"),
      Location))

# For the last Location, since it is tag failure
# Code location to be equal to the last known location
WeeklySiteVisit <- WeeklySiteVisit %>%
  group_by(TagID)%>%
  arrange(TagID, Week) %>%
  mutate(Location = 
           ifelse(EndSite == "Tag Failure",
                  case_when(
                    (lag(Location) == "FULL TRANSIT" &
                       lag(EndSite) %in% outside) ~ "OUTSIDE",
                    (lag(Location) == "FULL TRANSIT" &
                       lag(EndSite) %in% inside) ~ "INSIDE",
                    TRUE ~ "Tag Failure"),  
                  Location)) %>%
  mutate(Location = ifelse(EndSite == "Tag Failure" & 
                             lag(Location) == "PARTIAL TRANSIT",
                           "UNRESOLVED",Location),
         Location = ifelse(Location == "Tag Failure",lag(Location), Location))

#Remove Data from periods after the tag is believed to be shed
WeeklySiteVisit <- WeeklySiteVisit %>%
  left_join(., ShedTags, by = "TagID") %>%
  #week 257 is the end of the study
  mutate(ShedWeek = ifelse(is.na(ShedWeek),257,ShedWeek)) %>% 
  filter(Week < ShedWeek)

## Group Residency periods into Continuous Events =========================
CCF_Residency_CRE <- WeeklySiteVisit %>% 
  filter(!(TagID %in% no_detect_tags)) %>% #filter out undetected tags
  ungroup() %>%
  group_by(Species,TagID)%>%
  arrange(TagID, Week) %>%
  mutate(Location = as.character(Location)) %>%
  mutate(
    LocationChange = coalesce(Location != lag(Location) | 
                                # Do not combine transits!
                                grepl("TRANSIT", Location),
                              FALSE) %>% 
      cumsum()
  ) %>% 
  group_by(Species,TagID, LocationChange, Location) %>% 
  summarise(
    StartWeek = min(Week),
    EndWeek = max(Week), 
  ) %>% 
  ungroup() %>%
  arrange(TagID, StartWeek)

### Extend the Residency periods for fish with gaps -----------------------
# Gap is assumed to spent at the same location as the next site following
# Movement
CCF_Residency_CRE <- CCF_Residency_CRE %>%
  group_by(TagID) %>%
  arrange(TagID, StartWeek) %>%
  mutate(gap = lead(StartWeek)-EndWeek-1) %>%
  #mutate(gap = ifelse(is.na(gap),0,gap))
  mutate(
    StartWeek = case_when(
      #If tfile.path isn't a transit, just add time to the next location
      !grepl("TRANSIT", Location) &
        is.finite(lag(gap)) & lag(gap) > 0 ~ StartWeek - lag(gap), 
      #If tfile.path is a transit at both locations, do not extend
      grepl("TRANSIT", lag(Location)) & 
        grepl("TRANSIT", Location) ~ StartWeek,
      #If tfile.path is a transit at one location, extend the non-transit location
      grepl("TRANSIT", lag(Location)) &
        !grepl("TRANSIT",Location) &
        is.finite(lag(gap)) & lag(gap) > 0 ~ StartWeek - lag(gap),
      TRUE ~ StartWeek),
    EndWeek = case_when(
      grepl("TRANSIT",lead(Location)) & 
        !grepl("TRANSIT",Location) &
        is.finite(gap) & gap > 0 ~ EndWeek + gap,
      TRUE ~ EndWeek)
    )

# Estimate Age based on length ############################################
## Retrieve Length at capture from release data ===========================
CCF_Residency_CRE <- CCF_Residency_CRE %>%
  ungroup() %>%
  left_join(.,release[,c(1,7)], by = "TagID") %>%
  ungroup()

# Develop Age-length key based on a methods from 
# fishR Vignette (D.Ogle 2013), based on methods from 
# Iserman and Knight (2015)

## Striped Bass Aging =====================================================
# Building the age-length key table ####
set.seed(42069)

AgeData <- read.csv(file.path("1. Data","Inputs","AgedFish.csv"), header = T, stringsAsFactors = F)

stb.age <- AgeData[AgeData$Species == "Striped Bass" & AgeData$Length_cm >20 & AgeData$Age < 8,] %>%
  select(Length_cm, Age) %>%
  arrange(Age)

x <- stb.age$Age
y <- stb.age$Length_cm

alk.lm <- lm(y ~ I(log(x)))

predicted.intervals <- predict(alk.lm,data.frame(x=stb.age$Age),interval='confidence', level=0.995)

#plot(stb.age$Age, stb.age$Length_cm)
#lines(stb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
#lines(stb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1)
#lines(stb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1)

x_new <- rep(seq(2,7, by = 1), 100000)

alk.lm$fitted.values <- predict(alk.lm, data.frame(x = x_new))


y_new <- simulate(alk.lm)[,1]

#plot(x_new, y_new, col = 'red')
#points(stb.age$Age, stb.age$Length_cm, pch = 16, )
#lines(stb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
#lines(stb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1, lty = 2)
#lines(stb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1, lty = 2)


alk_new <- tibble(Age = x_new, Length_cm = y_new)

# Find the minimum size for use in lencat function
FSA::Summarize(~Length_cm, data = alk_new, digits = 1)

# use lencat function to define length categories for each fish starting at min 24, and increasing by 2 
stb.age1 <- FSA::lencat(~Length_cm, data = alk_new, startcat=0, w=1)

stb.raw <- with(stb.age1, table(LCat,Age))
stb.key <- prop.table(stb.raw, margin = 1)
round(stb.key,2)

### Assign Ages to individuals at capture using age-length key ------------
stb.len <- CCF_Residency_CRE %>%
  filter(Species == "Striped Bass") %>%
  distinct(TagID, Length_cm) %>%
  filter(!is.na(Length_cm))

stb.len1 <- FSA::alkIndivAge(stb.key, ~Length_cm, data = stb.len, type = "SR", seed = 42069)
FSA::alkPlot(stb.key)

stb.key.tib <- as_tibble(stb.key)

stb.key.expanded <- stb.key.tib %>%
  mutate(LCat = as.numeric(LCat),Age = as.numeric(Age)) %>%
  mutate(n = n*5) %>%
  map_df(., rep, .$n)

## Largemouth Bass Aging ==================================================
lmb.age <- AgeData[AgeData$Species == "Largemouth Bass",] %>%
  select(Length_cm, Age) %>%
  arrange(Age)

x <- lmb.age$Age
y <- lmb.age$Length_cm

alk.lm <- lm(y ~ I(log(x)))

predicted.intervals <- predict(alk.lm,data.frame(x=lmb.age$Age),interval='confidence', level=0.995)

#plot(lmb.age$Age, lmb.age$Length_cm)
#lines(lmb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
#lines(lmb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1)
#lines(lmb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1)

x_new <- rep(seq(2,4, by = 1), 100000)

alk.lm$fitted.values <- predict(alk.lm, data.frame(x = x_new))

y_new <- simulate(alk.lm)[,1]

#plot(x_new, y_new, col = 'red')
#points(lmb.age$Age, lmb.age$Length_cm, pch = 16, )
#lines(lmb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
#lines(lmb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1, lty = 2)
#lines(lmb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1, lty = 2)


alk_new <- tibble(Age = x_new, Length_cm = y_new)

# Find the minimum size for use in lencat function
FSA::Summarize(~Length_cm, data = alk_new, digits = 1)

# use lencat function to define length categories for each fish starting at min 24, and increasing by 2 
lmb.age1 <- FSA::lencat(~Length_cm,data=alk_new,startcat=18,w=1)

lmb.raw <- with(lmb.age1, table(LCat,Age))
lmb.key <- prop.table(lmb.raw, margin = 1)
round(lmb.key,2)

## Assign Ages to individuals at capture using age-length key =============
lmb.len <- CCF_Residency_CRE %>%
  filter(Species == "Largemouth Bass") %>%
  distinct(TagID, Length_cm) %>%
  filter(!is.na(Length_cm))

lmb.len1 <- FSA::alkIndivAge(lmb.key, ~Length_cm, data = lmb.len, type = "SR", seed = 42069)
FSA::alkPlot(lmb.key)

AssignedAge <- bind_rows(stb.len1, lmb.len1)

AssignedAge <- AssignedAge %>%
  rename("EstCaptureAge" = age)

CCF_Residency_CRE <- CCF_Residency_CRE %>%
  ungroup() %>%
  left_join(., AssignedAge, by = c("TagID"))

CCF_Residency_CRE <- CCF_Residency_CRE %>%
  group_by(TagID) %>%
  summarise(CaptureDate = min(studyweek_startdate(StartWeek))) %>%
  left_join(CCF_Residency_CRE,., by = "TagID")%>%
  mutate(CurrentDate = studyweek_startdate(StartWeek)) %>%
  mutate(TimePassed = CurrentDate - CaptureDate) %>%
  mutate(EstimatedAge =  floor(TimePassed/dyears(1))+EstCaptureAge) %>%
  mutate(CaptureWeek = date_to_studyweek(CaptureDate)) %>%
  select(Species, TagID, LocationChange, Location, StartWeek, EndWeek, Length_cm.x, EstCaptureAge, CaptureWeek, EstimatedAge)

# Cross each tag with each week (universal dataframe)
tag_x_weeks <- crossing(Tags,Week = Weeks)

## Determine the Age of an Individual at Week n and Weeks at Age a.
weekly_tag_age <- CCF_Residency_CRE %>%
  group_by(TagID) %>%
  summarise(CaptureWeek = min(StartWeek),
            CaptureAge = min(EstCaptureAge),
            EndWeek = max(EndWeek)) %>%
  left_join(tag_x_weeks) %>%
  mutate(time_from_capture = Week - CaptureWeek) %>%
  filter(time_from_capture >= 0 | Week > EndWeek) %>%
  mutate(EstimatedAge = CaptureAge + floor(time_from_capture/52)) %>%
  mutate(EstimatedAge = ifelse(is.na(EstimatedAge), -1, EstimatedAge),
         AgeBin = case_when(
           between(EstimatedAge,0,2) ~ "1-2",
           between(EstimatedAge,3,5) ~ "3-5",
           EstimatedAge >= 6 ~ "6+"
         )) %>%
  dplyr::select(TagID,Week,AgeBin)

weekly_tag_age <- weekly_tag_age %>%
  group_by(TagID,AgeBin) %>%
  summarise(StartWeek = min(Week),
            EndWeek = max(Week)) %>%
  mutate(time_at_age = EndWeek - StartWeek) %>%
  select(TagID, AgeBin, time_at_age) %>%
  right_join(weekly_tag_age)

CCF_Residency_CRE <- CCF_Residency_CRE %>%
  mutate(EstimatedAge = ifelse(is.na(EstimatedAge), -1, EstimatedAge),
         AgeBin = case_when(
           between(EstimatedAge,0,2) ~ "1-2",
           between(EstimatedAge,3,5) ~ "3-5",
           EstimatedAge >= 6 ~ "6+"
         )) %>%
  arrange(TagID, StartWeek)

# Examine Residency ########################################################
# Residency defined as being consecutively detected
# inside or outside of CCF without detection at the opposite location

## Determine the Duration of Continuous Events =============================
CCF_Residency_CRE <- CCF_Residency_CRE %>% 
  mutate(Duration_Weeks = as.integer(EndWeek - StartWeek + 1)) %>%
  filter(Species != "Beacon" | Location != "Tag Failure") %>%
  mutate(Location = factor(Location, levels = c("INSIDE","ESTIMATED INSIDE",
                                                "FULL TRANSIT", "PARTIAL TRANSIT",
                                                "ESTIMATED OUTSIDE","OUTSIDE",
                                                "UNRESOLVED")))

cre_plot <- ggplot(CCF_Residency_CRE[CCF_Residency_CRE$EstimatedAge > 0,] %>%
                     filter(Species == "Striped Bass")
                   )+
  geom_boxplot(aes(x = Location, y = Duration_Weeks,color = Location, group = Location),
               fill = "grey95")+
  facet_grid(Species~AgeBin)+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0)
  )

cre_plot

### Create a Summary of Residency -----------------------------------------
CRE_summary <- CCF_Residency_CRE %>%
  group_by(Species,TagID,Location) %>%
  mutate(MeanResidency_Weeks = mean(Duration_Weeks),
         TotalResidency_Weeks = sum(Duration_Weeks),
         events = n()) %>%
  select(Species, TagID,Location,MeanResidency_Weeks,TotalResidency_Weeks, events) %>%
  add_tally() %>%
  distinct(Species, TagID, Location,MeanResidency_Weeks,TotalResidency_Weeks, events)

CRE_Age_summary <- CCF_Residency_CRE %>%
  group_by(Species,TagID,Location,AgeBin) %>%
  mutate(MeanResidency_Weeks = mean(Duration_Weeks),
         TotalResidency_Weeks = sum(Duration_Weeks)) %>%
  select(Species, TagID,AgeBin,Location,MeanResidency_Weeks,TotalResidency_Weeks) %>%
  add_tally() %>%
  distinct(Species, TagID, AgeBin,Location,MeanResidency_Weeks,TotalResidency_Weeks)

# Examine Movement ########################################################
CRE_Movement <- CCF_Residency_CRE %>%
  ungroup() %>%
  arrange(TagID,StartWeek)%>%
  group_by(TagID) %>%
  mutate(move_direction = case_when(
    # If the Current Week is some form of transit, but prior is outside, and 
    # after is inside, then this is an entry
    Location %in% c("FULL TRANSIT","PARTIAL TRANSIT") & 
      lag(Location) == "OUTSIDE" & 
      lead(Location) == "INSIDE" ~ "ENTRY",
    # Opposite of above for exit
    Location %in% c("FULL TRANSIT","PARTIAL TRANSIT") & 
      lag(Location) == "INSIDE" & 
      lead(Location) == "OUTSIDE" ~ "EXIT",
    # If tfile.path is a Full Transit and both before and after are the same,
    # then assume a complete entry and exit or exit and entry
    Location %in% c("FULL TRANSIT") & 
      lag(Location) == lead(Location) ~ "FULL TRANSIT",
    # If its a Partial Transit and the both before and after are the same, 
    # then do not assume an actual transit, unresolved. 
    Location %in% c("PARTIAL TRANSIT") & 
      lag(Location) == lead(Location) ~ "UNRESOLVED TRANSIT",
    # If none of the above are true and tfile.path is a Full Transit, and the next 
    # location is outside, "EXIT"
    Location == "FULL TRANSIT" & 
      lead(Location) == "OUTSIDE" ~ "EXIT",
    # If none of the above are true and tfile.path is a Full Transit, and the next 
    # location is inside, "ENTRY"
    Location == "FULL TRANSIT" & 
      lead(Location) == "INSIDE" ~ "ENTRY",
    # If the location is not a full transit and the current location is
    # different from the last location ...
    !(Location %in% c("FULL TRANSIT")) & 
      Location != lag(Location) ~ case_when(
        # If the last location was "Outside" and now "Inside", Entry
        lag(Location) == "OUTSIDE" & Location == "INSIDE" ~ "ENTRY",
        # If the last location was "Inside" and now "Outside", Exit
        lag(Location) == "INSIDE" & Location == "OUTSIDE" ~ "EXIT",
        # If the next location is the same as the current, No Transit
        lead(Location) == Location ~ "NOT TRANSIT"),
    #Otherwise it's unresolved
    TRUE ~ "UNRESOLVED TRANSIT"))

# Edits for special cases
# Location 0 is "New Tag" if tfile.path is not movement
CRE_Movement$move_direction <- ifelse(
  CRE_Movement$move_direction == "UNRESOLVED TRANSIT" & 
    CRE_Movement$LocationChange == 0,
  "New Tag", 
  CRE_Movement$move_direction)
# If tfile.path is an NA, or an Unresolved transit and current tag is different from 
# the next tag, this is an Unresolved transit movement
CRE_Movement$move_direction <- ifelse(
  (is.na(CRE_Movement$move_direction)|
     CRE_Movement$move_direction == "UNRESOLVED TRANSIT") & 
    CRE_Movement$TagID != lead(CRE_Movement$TagID) & 
    CRE_Movement$Location %in% c("FULL TRANSIT", "PARTIAL TRANSIT"),
  "UNRESOLVED TRANSIT", 
  CRE_Movement$move_direction)

# finally all the rest must be fish remaining in one location
CRE_Movement$move_direction <- ifelse(
  is.na(CRE_Movement$move_direction),
  "NOT TRANSIT", 
  CRE_Movement$move_direction)

# Develop a table of Migrants #############################################
#Exits and Full Transits (at least 1 exit)
emigrants <- CRE_Movement %>%
  filter(Species %in% c("Striped Bass","Largemouth Bass")) %>%
  ungroup() %>%
  filter(move_direction %in% c("EXIT", "FULL TRANSIT")) %>%
  group_by(Species, TagID, AgeBin) %>%
  tally() %>%
  rename("Exits" = n)
#Entries and Full Transits (at least 1 entry)
immigrants <- CRE_Movement %>%
  filter(Species %in% c("Striped Bass","Largemouth Bass")) %>%
  ungroup() %>%
  filter(move_direction %in% c("ENTRY", "FULL TRANSIT")) %>%
  group_by(TagID, AgeBin) %>%
  tally() %>%
  rename("Entries" = n)

migrants <-emigrants %>%
  left_join(immigrants)

#fill 0s 
migrants$Exits <- replace_na(migrants$Exits,0)
migrants$Entries <- replace_na(migrants$Entries,0)
migrants <- migrants %>%
  pivot_longer(cols = c(Entries,Exits), names_to = "Movement", values_to = "Count")

ggplot(migrants %>% filter(Species == "Striped Bass")) +
  #geom_jitter(aes(x = AgeBin, y = Count, color = Movement))+
  geom_boxplot(aes(x = AgeBin, y = Count, color = Movement, 
                   group = interaction(AgeBin, Movement)))+
  facet_grid(Species ~.)+
  theme_classic()

# Splitting Each CRE by Time Periods ######################################
# Split all CREs by weeks to increase resolution
tbl_split_weeks <- CRE_Movement %>%
  mutate(StartDate = studyweek_startdate(StartWeek), 
         EndDate = studyweek_enddate(EndWeek)) %>%
  mutate(periods = map2(StartDate, EndDate,
                        split_interval_weeks, unit = "7 days")) %>%
  unnest(periods) %>%
  #Cut the last week at March 31, 2018
  mutate(.start = as.Date(ifelse(.start > as.Date("2018-03-31"),
                                 as.Date("2018-03-31"), 
                                 as.Date(.start)), 
                          origin = origin),
         .end = as.Date(ifelse(.end > as.Date("2018-03-31"),
                               as.Date("2018-03-31"), 
                               as.Date(.end)), 
                        origin = origin))

# CONNECTIVITY METRICS ###################################################
# Metrics created:
#  >Transit Frequency
#  >Weeks until first transit

## Transit Frequency =====================================================
# Transit frequency is the number of transits a fish makes
# Divided by their total time at large

### Individual Transit Frequency -----------------------------------------
# Count the total number of resolved transits a fish makes
transits <- CRE_Movement %>%
  filter(Species != "Catfish") %>%
  ungroup() %>%
  filter(move_direction %in% c("ENTRY","EXIT", "FULL TRANSIT") & Species != "Beacon") %>%
  group_by(Species, TagID) %>%
  tally %>%
  rename("transits" = n)

emmigrant_tags <- transits %>%
  ungroup() %>%
  distinct(TagID) %>%
  mutate(emmigrant = T)

# Determine the total number of weeks the fish was "at large"
t_large <- CRE_Movement %>%
  filter(Species != "Catfish") %>%
  ungroup() %>%
  group_by(TagID) %>%
  summarise(t_large = max(EndWeek)-min(StartWeek)+1)

# Calculate the transit frequency
trans_freq <- transits %>%
  left_join(.,t_large, by = "TagID") %>%
  mutate(transit_freq = transits/t_large) %>%
  mutate(emmigrant = T)

# For non-emmigrating fish, set transits to 0 and transit frequency to 0
trans_freq_nonemmigrant <- CRE_Movement %>%
  filter(Species != "Catfish") %>%
  ungroup() %>%
  group_by(Species, TagID) %>%
  left_join(.,emmigrant_tags,by = "TagID")%>%
  filter(is.na(emmigrant), Species != "Beacon") %>%
  summarise(t_large = max(EndWeek)-min(StartWeek)+1) %>%
  mutate(transit_freq = 0, transits = 0) %>%
  mutate(emmigrant = F)

trans_freq <- bind_rows(trans_freq,trans_freq_nonemmigrant)

trans_freq <- trans_freq %>%
  left_join(., AssignedAge, by = "TagID") %>%
  ungroup() %>%
  mutate(EstCaptureAge = ifelse(is.na(EstCaptureAge),-1,EstCaptureAge)) %>%
  mutate(AgeBin = case_when(
    between(EstCaptureAge,0,2) ~ "1-2",
    between(EstCaptureAge,3,5) ~ "3-5",
    EstCaptureAge >= 6 ~ "6+"
  ))%>%
  filter(!is.na(Length_cm)) %>%
  mutate(transit_freq = ifelse(transit_freq == 0,0.00000000001,transit_freq))

# Transit Frequency Plots
freq_boxplot_tot <- ggplot(trans_freq)+
  geom_boxplot(aes(x = Species, y = transit_freq, group = Species))+
  labs(title = "All Predator Transit Frequencies 2013-2018", y = "Transit Frequency (Transits / Week)")+
  theme_classic()

freq_boxplot_migrants <- ggplot(data = trans_freq[trans_freq$emmigrant == T,])+
  geom_boxplot(aes(AgeBin,y = transit_freq, group = AgeBin))+
  labs(title = "Emigrant Predator Transit Frequencies by Age 2013-2018", y = "Transit Frequency (Transits / Week)")+
  facet_grid(Species~.)+
  theme_classic()

freq_boxplot_age <- ggplot(trans_freq %>% filter(Species == "Striped Bass"))+
  geom_boxplot(aes(x = AgeBin, y = transit_freq, group = AgeBin))+
  labs(title = "Transit Frequencies by Age 2013-2018", y = "Transit Frequency (Transits / Week)")+
  facet_grid(Species ~.)+
  theme_classic()

freq_boxplot_tot

freq_boxplot_migrants

freq_boxplot_age

## Emmigrant Analysis =====================================================
Exits <- CRE_Movement %>%
  group_by(Species,TagID, move_direction) %>%
  filter(move_direction %in% c("EXIT","FULL TRANSIT")) %>%
  ungroup() %>% group_by(Species, TagID, move_direction) %>%
  summarise(Exit_1 = dplyr::first(StartWeek),
            Exit_2 = dplyr::nth(StartWeek,2),
            Exit_3 = dplyr::nth(StartWeek,3),
            Exit_4 = dplyr::nth(StartWeek,4),
            Exit_5 = dplyr::nth(StartWeek,5),
            Exit_6 = dplyr::nth(StartWeek,6),
            Exit_7 = dplyr::nth(StartWeek,7),
            Exit_8 = dplyr::nth(StartWeek,8)) %>%
  ungroup() %>%
  pivot_longer(cols = Exit_1:Exit_8, 
               names_to = c("Movement", "ID"),
               names_sep = "_",
               values_to = "Week") %>%
  arrange(TagID, Week) %>%
  filter(!is.na(Week)) %>%
  group_by(TagID) %>%
  mutate(ID = row_number())

Entrys <- CRE_Movement %>%
  group_by(Species,TagID, move_direction) %>%
  filter(move_direction %in% c("ENTRY","FULL TRANSIT")) %>%
  ungroup() %>% group_by(Species, TagID, move_direction) %>%
  summarise(Entry_1 = dplyr::first(StartWeek),
            Entry_2 = dplyr::nth(StartWeek,2),
            Entry_3 = dplyr::nth(StartWeek,3),
            Entry_4 = dplyr::nth(StartWeek,4),
            Entry_5 = dplyr::nth(StartWeek,5),
            Entry_6 = dplyr::nth(StartWeek,6),
            Entry_7 = dplyr::nth(StartWeek,7),
            Entry_8 = dplyr::nth(StartWeek,8)) %>%
  ungroup() %>%
  pivot_longer(cols = Entry_1:Entry_8, 
               names_to = c("Movement", "ID"),
               names_sep = "_",
               values_to = "Week") %>%
  arrange(TagID, Week) %>%
  filter(!is.na(Week)) %>%
  group_by(TagID) %>%
  mutate(ID = row_number())

NewTags <- CRE_Movement %>%
  group_by(Species,TagID, move_direction) %>%
  filter(move_direction %in% c("New Tag")) %>%
  ungroup() %>% group_by(Species, TagID, move_direction) %>%
  summarise(NewTag_1 = dplyr::first(StartWeek)) %>%
  ungroup() %>%
  pivot_longer(cols = NewTag_1, 
               names_to = c("Movement", "ID"),
               names_sep = "_",
               values_to = "Week") %>%
  arrange(TagID, Week) %>%
  filter(!is.na(Week)) %>%
  group_by(TagID) %>%
  mutate(ID = row_number())

emmigration_table <- NewTags %>%
  bind_rows(Exits) %>%
  bind_rows(Entrys) %>%
  arrange(TagID,Week,ID) %>%
  mutate(jdate = yday(studyweek_startdate(Week)),
         month = month(studyweek_startdate(Week)),
         wks_btwn = Week-lag(Week))

### Time to First Exit -----------------------------------------------------
first_exit <- emmigration_table %>%
  filter(Movement %in% c("NewTag","Exit") & ID == 1) %>%
  filter(!is.na(wks_btwn)) %>%
  select(Species, TagID, ttFirstExit = wks_btwn)

### Average Time between movements? ----------------------------------------
time_btwn_movements <- emmigration_table %>%
  left_join(weekly_tag_age) %>%
  filter(Movement != "New Tag") %>%
  filter(!is.na(wks_btwn)) %>%
  ungroup() %>%
  group_by(Species, TagID, AgeBin) %>%
  summarise(mean_time_btwn_mvmts = mean(wks_btwn),
            max_time_btwn_mvmts = max(wks_btwn))

### Weeks from release to first exit ---------------------------------------

#Find all fish that ever changed location at least once and create a new column "emmigrant" and set the value = 1
CRE_emmigrants <- CCF_Residency_CRE %>%
  group_by(TagID)%>%
  filter(LocationChange > 0) %>%
  distinct(TagID)%>%
  mutate(emmigrant = 1)

# For all fish, find their first CRE (Location 0) and pull the Tag ID, Species, and how long the CRE lasts
cre_surv <- CCF_Residency_CRE  %>%
  group_by(TagID)%>%
  filter(LocationChange == 0) %>%
  select(TagID, Species, AgeBin, Duration_Weeks) %>%
  left_join(.,CRE_emmigrants, by = "TagID")

# For any Fish which hasn't been coded as an emmigrant above
# set emmigrant status to 0
cre_surv$emmigrant[is.na(cre_surv$emmigrant)]<-0

#Create a new column time which transforms duration_sec to days
cre_surv$time <- as.numeric(cre_surv$Duration_Weeks)

#### Fit a Kaplan Meier Survivor Model to the data -------------------------

# the column emmigrant is being used to censor data
# if a fish was not "witnessed" emigrating it is censored
km <- with(cre_surv[cre_surv$Species != "Beacon",], survival::Surv(time,emmigrant))
head(km,1000)


## The First Survival Curve is for All fish (censored)
km_fit <- survfit(Surv(time, emmigrant)~1, data = cre_surv[cre_surv$Species != "Beacon",])
summary(km_fit, times = c(seq(1, 1550, 100)))

autoplot(km_fit) 

## Next try to separate by Species
km_fit_sp <- survival::survfit(survival::Surv(time, emmigrant)~Species, data = cre_surv[cre_surv$Species != "Beacon",])
autoplot(km_fit_sp)

pdf(file.path("1. Data", "Figures", "TimeToEmmigration_All.pdf"), width = 10, height =7)
autoplot(km_fit_sp)
dev.off()

## Fit our Age Bins for Striped Bass
cre_surv_age <- cre_surv %>%
  filter(Species == "Striped Bass")

km_fit_age <- survival::survfit(survival::Surv(time, emmigrant)~AgeBin, data = cre_surv_age)
autoplot(km_fit_age)+theme_classic()

pdf(file.path("1. Data", "Figures", "TimeToEmmigration_StripedBass_x_Age.pdf"), 
    width = 12, height =5)
autoplot(km_fit_age)+theme_classic()+
  labs(title = "Time to First Exit", y = "Percent of Population Non-Emmigrated",
       x = "Time (weeks)")
dev.off()

# Prepare the data for Cluster analysis ###################################

## Total Receivers detected On / Time at Large ============================
DetectedRcvrs <- TotalSiteVisitSummary %>%
  select(Species, TagID, InsideSiteVisits = inside_sites, 
         OutsideSiteVisits = outside_sites, numSites, TAL) %>%
  mutate(lifetime_inside = InsideSiteVisits/TAL,
         lifetime_outside = OutsideSiteVisits/TAL,
         lifetime_all = numSites/TAL)

## Receivers/Wk detected On ================================================
WkDetectedRcvrs <- WeeklySiteVisit %>%
  filter(!(grepl("Tag Failure", SiteVisits))) %>%
  ungroup() %>%
  left_join(weekly_tag_age, by = c("TagID","Week")) %>%
  group_by(TagID, AgeBin) %>%
  summarise(mean_WkRcvrs = mean(NumSites, na.rm = TRUE),
            max_WkRcvrs = max(NumSites, na.rm = TRUE))

## Estimated Capture Age ===================================================
CaptureAge <- AssignedAge %>%
  select(TagID, EstCaptureAge) %>%
  mutate(CaptureAgeBin = case_when(
    between(EstCaptureAge,0,2) ~ "1-2",
    between(EstCaptureAge,3,5) ~ "3-5",
    EstCaptureAge >= 6 ~ "6+"
  )) %>%
  select(TagID, CaptureAgeBin)

## Weeks Undetected ========================================================
Wks_Undetected <- TagBKM_Bin %>%
  filter(SiteCode != "Tag Failure") %>%
  ungroup() %>%
  left_join(weekly_tag_age, by = c("TagID","Week")) %>%
  group_by(Species, TagID, AgeBin) %>%
  distinct(Species, TagID, AgeBin, Week) %>%
  group_by(Species, TagID, AgeBin) %>%
  mutate(time_btwn = Week - lag(Week)) %>%
  filter(!is.na(time_btwn)) %>%
  summarise(max_time_btwn_detects = max(time_btwn, na.rm = TRUE),
            avg_time_btwn_detects = mean(time_btwn,na.rm = TRUE))

## Weeks btwn Movements ====================================================
time_btwn_movements

## Quantiles of Distance per Week ==========================================
### Assign minimum distance travelled in one week --------------------------

#Read in Linear Distance
distances <- read.csv(file.path("1. Data", "Inputs", "HydDistance.csv"))

#Filter Detections
TagBKM_Bin %>%
  filter(Species != "Beacon", !(SiteCode %in% c("Release",NA,"Tag Failure"))) %>%
  left_join(weekly_tag_age, by = c("TagID","Week")) %>%
  group_by(Species, Week, TagID, AgeBin) %>%
  summarise(SiteVisits = unique(list(as.character(SiteCode)))) %>%
  right_join(TagBKM_Bin %>% filter(!(SiteCode %in% c("Release",NA,"Tag Failure")))) %>%
  left_join(distances, by = c("SiteCode" = "SiteA"), relationship = "many-to-many") %>%
  mutate(NumSites = lengths(SiteVisits)) %>%
  group_by(Species, Week, TagID, AgeBin) %>%
  filter(SiteB %in% unlist(SiteVisits)) %>%
  group_by(Species, Week, TagID, AgeBin, NumSites) %>%
  summarise(weekly_max_distance = max(Distance),
            weekly_mean_distance = mean(Distance)) -> WeeklyDistance

### Summarize Distance -----------------------------------------------------
WeeklyDistance %>%
  ungroup() %>%
  group_by(Species, TagID, AgeBin) %>%
  summarise(mean_distance = mean(weekly_mean_distance),
            q25_distance = quantile(weekly_max_distance, 0.25),
            q50_distance = quantile(weekly_max_distance, 0.50),
            q75_distance = quantile(weekly_max_distance, 0.75),
            q99_distance = quantile(weekly_max_distance, 0.99)) -> Distances

## Emigration Status =======================================================
emmigrant_tags

## Estimated Age at First Exit =============================================
FirstExitAge <- AssignedAge %>%
  left_join(first_exit) %>%
  filter(!is.na(ttFirstExit)) %>%
  mutate(FirstExitAge = as.integer(EstCaptureAge+(ttFirstExit/52))) %>%
  mutate(FirstExitAgeBin = case_when(
    between(FirstExitAge,0,2) ~ "1-2",
    between(FirstExitAge,3,5) ~ "3-5",
    FirstExitAge >= 6 ~ "6+"
  )) %>%
  select(TagID, FirstExitAgeBin)

## Avg Week of Year of Exit ================================================
avg_exit_wk <- emmigration_table %>%
  filter(Movement == "Exit") %>%
  select(TagID, Week, jdate, month) %>%
  mutate(woy = lubridate::week(studyweek_startdate(Week))) %>%
  group_by(TagID) %>%
  summarise(mean_exit_woy = mean(woy))

## Avg Week of Year of Entry ===============================================
avg_entry_wk <- emmigration_table %>%
  filter(Movement == "Entry") %>%
  select(TagID, Week, jdate, month) %>%
  mutate(woy = lubridate::week(studyweek_startdate(Week))) %>%
  group_by(TagID) %>%
  summarise(mean_entry_woy = mean(woy))
  
## Number of Exits and Entries =============================================
TotalMovement <- emmigration_table %>%
  left_join(weekly_tag_age, by = c("TagID","Week")) %>%
  filter(Movement %in% c("Exit", "Entry")) %>%
  group_by(TagID, AgeBin, Movement) %>%
  count() %>%
  pivot_wider(names_from = "Movement", values_from = "n", values_fill = 0) %>%
  rowwise() %>%
  mutate(All_Moves = rowSums(across(Entry:Exit))) %>%
  filter(!is.na(AgeBin))

MovementFrequency <- TotalMovement %>%
  left_join(
    distinct(
      select(weekly_tag_age, TagID, AgeBin, time_at_age)
      )) %>%
  mutate(Exit_freq = Exit/time_at_age,
         Entry_freq = Entry/time_at_age,
         All_Moves_freq = All_Moves/time_at_age) %>%
  select(TagID, AgeBin, Exit_freq:All_Moves_freq)

## Average Residence Time ==================================================
avg_residence <- CRE_Age_summary %>%
  filter(Location %in% c("INSIDE","OUTSIDE")) %>%
  select(TagID, Location, MeanResidency_Weeks) %>%
  pivot_wider(names_from ="Location",
              names_prefix = "AvgRes_",
              values_from = "MeanResidency_Weeks",
              values_fill = 0)

total_residence <- CRE_Age_summary %>%
  filter(Location %in% c("INSIDE","OUTSIDE")) %>%
  select(TagID, Location, TotalResidency_Weeks) %>%
  pivot_wider(names_from ="Location",
              names_prefix = "TotRes_",
              values_from = "TotalResidency_Weeks",
              values_fill = 0)

## Join the Data
dat <- DetectedRcvrs %>%
  left_join(WkDetectedRcvrs) %>%
  left_join(CaptureAge,) %>%
  left_join(Wks_Undetected) %>%
  left_join(time_btwn_movements) %>%
  left_join(Distances) %>%
  left_join(emmigrant_tags) %>%
  left_join(FirstExitAge) %>%
  left_join(avg_exit_wk) %>%
  left_join(avg_entry_wk) %>%
  left_join(TotalMovement) %>%
  left_join(MovementFrequency) %>%
  left_join(avg_residence) %>%
  left_join(total_residence) %>%
  filter(Species %in% c("Striped Bass",
                        "Largemouth Bass"
  )) %>%
  filter(!(TagID %in% no_detect_tags)) %>% # Remove tags with no detections
  replace_na(list(Species = NA,
                  TagID = NA,
                  lifetime_inside = 0,
                  lifetime_outside = 0,
                  InsideSiteVisits = 0,
                  OutsideSiteVisits = 0,
                  numSites = 0,
                  mean_WkRcvrs = 0,
                  max_WkRcvrs = 0,
                  AgeBin = "0", 
                  max_time_btwn_detects = 0,
                  avg_time_btwn_detects = 0,
                  mean_time_btwn_mvmts = 0,
                  max_time_btwn_mvmts = 0,
                  mean_distance = 0,
                  q25_distance = 0, 
                  q50_distance = 0, 
                  q75_distance = 0, 
                  q99_distance = 0, 
                  emmigrant = FALSE, 
                  FirstExitAgeBin = "0", 
                  mean_exit_woy = 0, 
                  mean_entry_woy = 0,
                  Entry = 0, 
                  Exit = 0, 
                  All_Moves = 0,
                  Exit_freq = 0,
                  Entry_freq = 0,
                  All_Moves_freq = 0,
                  AvgRes_INSIDE = 0, 
                  AvgRes_OUTSIDE = 0, 
                  TotRes_INSIDE = 0, 
                  TotRes_OUTSIDE = 0)) %>%
  mutate(TagAge = paste(TagID, AgeBin, sep = "-"),
         Species = ifelse(Species == "Striped Bass",1,2),
         AgeBin = as.numeric(factor(AgeBin, levels = c("0","1-2","3-5","6+"))),
         FirstExitAgeBin = as.numeric(factor(FirstExitAgeBin, 
                                             levels = c("0","1-2","3-5","6+"))),
         CaptureAgeBin = as.numeric(factor(CaptureAgeBin, levels = c("0","1-2","3-5","6+"))),
  ) %>%
  filter(!is.na(CaptureAgeBin)) %>%
  column_to_rownames(var = "TagAge") %>%
  select(-c(TagID)) %>%
  filter(Species == 1)

## Conduct Cluster Analysis ===============================================
#Following http://stratigrafia.org/8370/lecturenotes/clusterAnalysis.html

library(vegan)
set.seed(123)

### Formatting ------------------------------------------------------------
datT1 <- decostand(dat, method = "total")
datT2 <- decostand(datT1, method = "max")
datBray <- vegdist(datT2)

### Create the Dendogram --------------------------------------------------
bray_agnes <- cluster::agnes(datBray, method = "ward")

plot(bray_agnes, which.plots = 2, main="Movement Groups", cex=0.1)
abline(h=4, col="red", lwd = 2)
abline(h=5, col="yellow2", lwd = 2)
abline(h=8, col="green", lwd = 2)

### Choosing the right clustering algorithm --------------------------------
# https://www.datanovia.com/en/lessons/choosing-the-best-clustering-algorithms/

# Compute clValid
library(clValid)

clmethods <- c("kmeans", "diana", "fanny", "pam", "clara", "agnes")

valid <- clValid::clValid(datT2, nClust= 3:5, clMethods = clmethods,
                             validation = c("internal","stability"), verbose = TRUE)
measures <- as_tibble(clValid::measures(valid))
measures$Measure = c("APN","AD","ADM","FOM","Connectivity","Dunn","Silhouette")

measures <- measures %>%
  pivot_longer(cols = first(names(measures)):nth(names(measures), n = -2), 
               names_to = c("k","method"), 
               names_sep = "\\.", values_to = "score") %>%
  pivot_wider(id_cols = c(method,k), names_from = Measure, values_from = score)

# Check all possible sorting of the ranking measures to determine which
# method and k appear in each rank and the amount of representation

# Dunn and Silhouette measurements are meant to be maximized, so multiplying by
# -1 reverts the scale for comparing ranks
measures_scaled <- measures %>%
  mutate(Silhouette = -1*Silhouette,
         Dunn = -1*Dunn)

ggplot(measures_scaled %>%
         mutate(k = as.character(k)) %>%
         select(-c(method,k)) %>%
         scale() %>%
         bind_cols(measures %>% select(method, k)) %>%
         pivot_longer(cols = APN:Silhouette, 
                      names_to = "test", 
                      values_to = "values"))+
  geom_bar(aes(x = interaction(method,k), y = values), stat = "identity")+
  facet_wrap(vars(test))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

names <- c("Connectivity", "Silhouette", "Dunn", "APN", "AD", "ADM", "FOM")

cross <- crossing("a" = names, "b" = names,
         "c" = names, "d" = names,
         "e" = names, "f" = names,
         "g" = names) %>%
  filter(a != b & a != c & a != d & a != e & a != f & a != g &
           b != c & b != d & b != e & b != f & b != g &
           c != d & c != e & c != f & c != g &
           d != e & d != f & d != g &
           e != f & e != g &
           f != g)

#This loop runs a while
gof <- data.frame("T1" = rep(NA,times = nrow(measures)))
for( i in 1:nrow(cross)){
  v <- measures_scaled %>%
    arrange(across(all_of(as.character(cross[i,])))) %>%
    unite("method_k", method, k, sep = "_") %>%
    pull(method_k)


  gof <- gof %>%
    bind_cols(.,v)
}
saveRDS(gof,"gof.rds")
gof <- readRDS("gof.rds")

gof_pvt <- gof %>%
  as_tibble() %>%
  select(-c(T1)) %>%
  rowid_to_column(var = "rank") %>%
  pivot_longer(cols = c(2:5041), names_to = "iteration", 
               values_to = "method_k")

gof_summary <- gof_pvt %>%
  group_by(rank, method_k) %>%
  count() %>%
  ungroup() %>%
  arrange(rank, desc(n))

gof_summary
# kmeans_2, pam_5, kmeans_3,

### Review clusters --------------------------------------------------------
# kmeans with k set to 2 from above
kmeans <- factoextra::eclust(datT2, "kmeans", k = 2, 
                          nboot = 5000,
                          hc_metric = "manhattan", 
                          hc_method = "ward.D2", 
                          graph = FALSE)

# Visualize cluster 
factoextra::fviz_cluster(kmeans, labelsize = 0,
                      as.ggplot = TRUE)+
  theme_classic()+
  labs(title = "Divisive Analysis, k = 2")

#Hierarchical K-means with k=2 (similar to diana k = 2 above, different result)
hkmeans_2 <- factoextra::hkmeans(datT2,
                                 k = 2,
                                 hc.metric = "manhattan",
                                 iter.max = 5000)
factoextra::fviz_cluster(hkmeans_2,
                         as.ggplot = TRUE, stand = TRUE,
                         labelsize = 0)+
  theme_classic()+
  labs(title = "Hierarchical K-means, k = 5")

# Clara with k set to 3 from above
clara_3 <- factoextra::eclust(datT2, "clara", k = 3, 
                              nboot = 5000,
                              hc_metric = "manhattan", 
                              hc_method = "ward.D2", 
                              graph = FALSE)

# Visualize cluster
factoextra::fviz_cluster(clara_3,
                         as.ggplot = TRUE, labelsize = 0)+
  theme_classic()+
  labs(title = "Clustering Large Applications, k = 3")

#Kmeans with k set to 3 from above
kmeans <- factoextra::eclust(datT2, "kmeans", k = 3, 
                             nboot = 5000,
                             hc_metric = "manhattan", 
                             hc_method = "ward.D2", 
                             graph = FALSE)
# Visualize cluster
factoextra::fviz_cluster(kmeans,
                         as.ggplot = TRUE, labelsize = 0)+
  theme_classic()+
  labs(title = "K-means, k = 3")

#Hierarchical K-means with k=3 (similar to K-means k = 3, above)
hkmeans_3 <- factoextra::hkmeans(datT2,
                                 k = 3,
                                 hc.metric = "manhattan",
                                 iter.max = 5000)

factoextra::fviz_cluster(hkmeans_3,
                         as.ggplot = TRUE, stand = TRUE,
                         labelsize = 0)+
  theme_classic()+
  labs(title = "Hierarchical K-means, k = 3")

#Hierarchical K-means with k=3 (similar to K-means k = 3, above)
hkmeans_4 <- factoextra::hkmeans(datT2,
                                 k = 4,
                                 hc.metric = "manhattan",
                                 iter.max = 5000)

factoextra::fviz_cluster(hkmeans_4,
                         as.ggplot = TRUE, stand = TRUE,
                         labelsize = 0)+
  theme_classic()+
  labs(title = "Hierarchical K-means, k = 4")
# PAM with k set to 5 from above
pam_4 <- factoextra::eclust(datT2, "pam", k = 4, 
                          nboot = 5000,
                          hc_metric = "manhattan", 
                          hc_method = "ward.D2", 
                          graph = FALSE)
factoextra::fviz_cluster(pam_4,
                         as.ggplot = TRUE, stand = TRUE,
                         labelsize = 0)+
  theme_classic()+
  labs(title = "PAM, k = 4")
# PAM with k set to 5 from above
pam_5 <- factoextra::eclust(datT2, "pam", k = 5, 
                            nboot = 5000,
                            hc_metric = "manhattan", 
                            hc_method = "ward.D2", 
                            graph = FALSE)

# Visualize cluster
factoextra::fviz_cluster(pam_5,
                      as.ggplot = TRUE, labelsize = 0)+
  theme_classic()+
  labs(title = "PAM, k = 5")

# Clara with k set to 4 from above
clara_4 <- factoextra::eclust(datT2, "clara", k = 4, 
                              nboot = 5000,
                              hc_metric = "manhattan", 
                              hc_method = "ward.D2", 
                              graph = FALSE)

# Visualize cluster
factoextra::fviz_cluster(clara_4,
                         as.ggplot = TRUE, labelsize = 0)+
  theme_classic()+
  labs(title = "Clustering Large Applications, k = 5")

# PAM with k set to 5 from above
clara_5 <- factoextra::eclust(datT2, "clara", k = 5, 
                            nboot = 5000,
                            hc_metric = "manhattan", 
                            hc_method = "ward.D2", 
                            graph = FALSE)

# Visualize cluster
factoextra::fviz_cluster(clara_5,
                         as.ggplot = TRUE, labelsize = 0)+
  theme_classic()+
  labs(title = "Clustering Large Applications, k = 5")

#Hierarchical K-means with k=5 (similar to PAM k = 5, above)
hkmeans_5 <- factoextra::hkmeans(datT2,
                                 k = 5,
                                 hc.metric = "manhattan",
                                 iter.max = 100000)
factoextra::fviz_cluster(hkmeans_5,
                         as.ggplot = TRUE,
                         labelsize = 0)+
  theme_classic()+
  labs(title = "Hierarchical K-means, k = 5")

## Select Appropriate Cluster ================================================
# 3 or 5 clusters has relatively good support (long "legs", similarity is
# determined by leg height, short legs = more similar, similar height bifurcation
# mean similar points of separation)
  
# Test Cluster AU (Takes FOREVER) Gives an idea of statistical confidence
# Biological behavior data, understood to have a lower statistical confidence
# AU = confidence in assignment (max 100, min 0)
# library(pvclust)
# set.seed(123)
# res.pv <- pvclust::pvclust(datT2 %>% t(), method.hclust = "ward.D2",
#                            method.dist = "manhattan", nboot = 1000,
#                            parallel = TRUE, r= seq(0.5,1.5, by = .1))
# plot(res.pv, hang = -1, cex = 0.5)


# #Decent support for 3 or 5 bifurcations
cluster_assign <- data.frame("cluster" = factor(hkmeans_5$cluster,
                                                levels = c(1,2,3,4,5))) %>%
  rownames_to_column("TagID")

clustered_dat <- dat %>%
  rownames_to_column("TagID") %>%
  left_join(cluster_assign)

### Examine the variables by cluster --------------------------------------
species_plot <- ggplot(clustered_dat %>%
                         group_by(Species) %>%
                         mutate(total = n()) %>%
                         group_by(Species, cluster) %>%
                         summarise(prop = n()/total) %>%
                         distinct()
                         )+
  geom_bar(aes(fill = factor(Species),
               x = factor(Species), y = prop),
           position = "dodge",
           stat = "identity")+
  facet_grid(cluster ~ .)+
  coord_flip()+
  scale_fill_discrete(name = "Species",
                      breaks = c(1,2),
                      labels = c("Striped Bass", "Largemouth Bass"))+
  scale_y_continuous(name = "Proportion of Species Total",
                     labels = scales::percent)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position = "top")
species_plot

age_plot <- ggplot(clustered_dat %>%
                     group_by(Species, AgeBin) %>%
                     mutate(total = n()) %>%
                     group_by(Species, AgeBin, cluster) %>%
                     summarise(prop = n()/total) %>%
                     distinct()
)+
  geom_bar(aes(fill = factor(AgeBin),
               x = factor(AgeBin), y = prop),
           position = "dodge",
           stat = "identity")+
  facet_grid(cluster ~ .)+
  coord_flip()+
  scale_fill_discrete(name = "Age",
                      breaks = c(2,3,4),
                      labels = c("1-2","3-5","6+"))+
  scale_y_continuous(name = "Proportion of Species Total",
                     labels = scales::percent)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position = "top")
age_plot

exit_age_plot <- ggplot(clustered_dat %>%
                     group_by(Species, FirstExitAgeBin) %>%
                     mutate(total = n()) %>%
                     group_by(Species, FirstExitAgeBin, cluster) %>%
                     summarise(prop = n()/total) %>%
                     distinct()
)+
  geom_bar(aes(fill = factor(FirstExitAgeBin),
               x = factor(FirstExitAgeBin), y = prop),
           position = "dodge",
           stat = "identity")+
  facet_grid(cluster ~ .)+
  coord_flip()+
  scale_fill_discrete(name = "Age at First Exit",
                      breaks = c(1,2,3,4),
                      labels = c("NA","1-2","3-5","6+"))+
  scale_y_continuous(name = "Proportion of Species Total",
                     labels = scales::percent)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position = "top")
exit_age_plot

sites_plot1 <- ggplot(clustered_dat %>%
                       select(cluster,TagID,"Inside" = InsideSiteVisits,
                              "Outside" = OutsideSiteVisits,
                              "All" = numSites) %>%
                       pivot_longer(cols = Inside:All,
                                    names_to = "Site",
                                    values_to = "Total"))+
  geom_boxplot(aes(y = Total, x = Site, fill = Site))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  labs(y = "Total Number of Site Visits")
sites_plot1

sites_plot2 <- ggplot(clustered_dat %>%
                        select(cluster,TagID,"Inside" = lifetime_inside,
                               "Outside" = lifetime_outside,
                               "All" = lifetime_all) %>%
                        pivot_longer(cols = Inside:All,
                                     names_to = "Site",
                                     values_to = "Total"))+
  geom_boxplot(aes(y = Total, x = Site, fill = Site))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  labs(y = "Average Number of Sites Visits per Week at Large")
sites_plot2

wk_rcvrs_plot <- ggplot(clustered_dat %>%
                          select(cluster,TagID,mean_WkRcvrs,max_WkRcvrs) %>%
                          pivot_longer(cols = mean_WkRcvrs:max_WkRcvrs,
                                       names_to = "stat",
                                       values_to = "value"))+
  geom_violin(aes(y = value, x = stat, fill = stat))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  scale_x_discrete(breaks = c("mean_WkRcvrs","max_WkRcvrs"),
                   labels = c("Mean", "Max"))+
  labs(x = "Statistic", y = "Number of Receivers per Week")
wk_rcvrs_plot

time_detects_plot <- ggplot(clustered_dat%>%
                              select(cluster,TagID,avg_time_btwn_detects,
                                     max_time_btwn_detects) %>%
                              pivot_longer(cols = 
                                             avg_time_btwn_detects:max_time_btwn_detects,
                                           names_to = "stat",
                                           values_to = "value"))+
  geom_boxplot(aes(y = value, x = stat, fill = stat))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  labs(x = "Statistic", y = "Time between Detections")+
  scale_y_log10()+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
time_detects_plot

time_move_plot <- ggplot(clustered_dat %>%
                          select(cluster,TagID,mean_time_btwn_mvmts,
                                 max_time_btwn_mvmts) %>%
                          pivot_longer(cols = mean_time_btwn_mvmts:
                                       max_time_btwn_mvmts,
                                       names_to = "stat",
                                       values_to = "value"))+
  geom_boxplot(aes(y = value, x = stat, fill = stat))+
  facet_grid(cluster ~ ., scales = 'free_y')+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  scale_x_discrete(breaks = c("mean_time_btwn_mvmts",
                              "max_time_btwn_mvmts"),
                   labels = c("Mean", "Max"))+
  labs(x = "Statistic", y = "Weeks between Transits")
time_move_plot

distance_plot <- ggplot(clustered_dat %>% 
         select(cluster, q25_distance:q99_distance) %>%
         pivot_longer(cols = q25_distance:q99_distance, 
                      names_to = "metric",
                      values_to = "distance"),
         aes(x = metric, y = distance, group = cluster))+
  stat_summary(fun.data=mean_se,geom="linerange",color="black")+
  stat_summary(fun.y=mean,geom="line")+
  facet_grid(cluster ~ .)+
  scale_x_discrete(name = "Quantile",
                   breaks = c("q25_distance","q50_distance","q75_distance",
                              "q99_distance"),
                   labels = c("25%","50%","70%","99%"))+
  scale_y_continuous(name = "Mean Weekly Distance Travelled")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.background = element_rect(color = "grey50"))
distance_plot

emmigrant_plot <- ggplot(clustered_dat %>%
                           group_by(cluster) %>%
                           mutate(total = n()) %>%
                           group_by(emmigrant, cluster) %>%
                           summarise(prop = n()) %>%
                           distinct())+
  geom_bar(aes(x = emmigrant, y = prop, fill = emmigrant), stat = "identity")+
  scale_fill_discrete(name = "",
                      breaks = c(TRUE,FALSE),
                      labels = c("Emigrant","Non-Emigrant"))+
  facet_grid(cluster ~ .)+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        legend.position = "top")+
  labs(x = "", y = "Number of Individuals")
emmigrant_plot

mean_exit_woy <- ggplot(clustered_dat%>%
                          separate(TagID,into = c("TagID","AgeBin"), sep = "-",
                                   extra = "merge") %>%
                          filter(TagID %in% c(emmigrant_tags$TagID)))+
  geom_violin(aes(y = lubridate::ymd("2017-12-31")+(mean_exit_woy*7), x = 1), 
              fill = "grey50",
              trim = FALSE, bw = 7)+
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip()+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b-%d")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(title = "Mean Week of Exit")
mean_exit_woy

mean_entry_woy <- ggplot(clustered_dat%>%
                           separate(TagID,into = c("TagID","AgeBin"), sep = "-",
                                    extra = "merge") %>%
                           filter(mean_entry_woy > 0))+
  geom_violin(aes(y = lubridate::ymd("2017-12-31")+(mean_entry_woy*7), x = 1), fill = "grey50",
              scale = "width",
              trim = FALSE, bw = 7)+
  facet_grid(cluster ~ .,drop = FALSE)+
  coord_flip()+
  theme_classic()+
  scale_y_date(name = "Date", date_labels = "%b-%d")+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(title = "Mean Week of Entry")
mean_entry_woy 

movements <- ggplot(clustered_dat %>%
                       select(cluster,Exit,Entry, 
                              All_Moves) %>%
                       pivot_longer(cols = Exit:All_Moves,
                                    names_to = "movement",
                                    values_to = "count"))+
  geom_boxplot(aes(y = count,
                  x = 1,
               fill = movement),
           )+
  facet_grid(cluster ~ .)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
        )+
  scale_y_continuous(name = "Number of Transits",
                     breaks = seq(0,10,2))+
  scale_x_continuous(name = "", limits= c(0.3,1.6))+
  scale_fill_discrete(name = "")+
  coord_flip()
movements

move_freq <- ggplot(clustered_dat %>%
                      select(cluster,Exit_freq,Entry_freq, 
                             All_Moves_freq) %>%
                      pivot_longer(cols = Exit_freq:All_Moves_freq,
                                   names_to = "movement",
                                   values_to = "count"))+
  geom_boxplot(aes(y = count,
                   x = 1,
                   fill = movement),
  )+
  facet_grid(cluster ~ .)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "top"
  )+
  scale_y_continuous(name = "Transit Frequency (per Week)")+
  scale_x_continuous(name = "")+
  coord_flip()
move_freq

residency <- ggplot(clustered_dat %>%
                    select(cluster, AvgRes_INSIDE, AvgRes_OUTSIDE,
                           TotRes_INSIDE, TotRes_OUTSIDE) %>%
                    pivot_longer(cols = AvgRes_INSIDE:TotRes_OUTSIDE,
                                 names_to = c("stat","location"),
                                 values_to = "value",
                                 names_sep = "_"))+
  geom_boxplot(aes(x = location, y = value, fill = stat), 
              bw = 0.5, trim = FALSE)+
  facet_grid(cluster ~ .)+
  scale_x_discrete(name = "",
                   breaks = c("INSIDE","OUTSIDE"),
                   labels = c("Inside","Outside"))+
  scale_fill_discrete(name = "Statistic",
                      breaks = c("AvgRes","TotRes"),
                      labels = c("Average","Total"))+
  scale_y_continuous(name = "Weeks of Residency")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "top")+
  coord_flip()
residency

pdf("All_Plots.pdf", width = 15, height = 8)
factoextra::fviz_cluster(hkmeans_5,
                         as.ggplot = TRUE, labelsize = 0)+
  theme_classic()+
  labs(title = "Hierarchical k-means, k = 5")

gridExtra::grid.arrange(species_plot, age_plot, emmigrant_plot, nrow = 1)
gridExtra::grid.arrange(sites_plot1 + theme(legend.position = "top"), 
                        sites_plot2 + theme(legend.position = "top"), 
                        residency, nrow = 1)
gridExtra::grid.arrange(time_detects_plot, time_move_plot, distance_plot, nrow = 1)
gridExtra::grid.arrange(movements, mean_exit_woy, mean_entry_woy,
                        nrow = 1)
dev.off()

for(i in unique(clustered_dat$TagID)){
  tmp <- TagBKM_Bin %>%
    filter(TagID == i) %>%
    filter(SiteCode != "Tag Failure") %>%
    mutate(TagID = as.character(TagID)) %>%
    left_join(cluster_assign)
  
  if(nrow(tmp)>1){
    ggplot()+
      geom_point(data = tmp, 
                 aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
      geom_step(data = tmp, 
                aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
      scale_x_date(name = "Date",
                   breaks = scales::breaks_pretty(n = 8),
                   date_labels = "%b-%d")+
      scale_y_discrete()+
      theme_classic()+
      theme(panel.spacing = unit(0, "lines"),
            axis.text.x = element_text(angle = -45, hjust = 0))+
      facet_grid(TagID ~.)
    ggsave(file.path("1. Data","Figures","DetectionPlots",paste(tmp$TagID[1],tmp$cluster[1],"ex.pdf",sep = "-")),
                device = "pdf", width = 5, height = 2.5)
  } else {
    print("Tag not in cluster")
  }
}
