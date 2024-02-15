######      CCFPS TECHNICAL REPORT: Build Detection Histories           #######

source(file.path("2. Code","0_Setup.R"))

# Format Data #############################################################
# This script picks up where original BKM_QAQC script ended (not included here)

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

## Read in the Release Data and Format
release <- read.csv(file.path("1. Data","Inputs","ReleaseData.csv"), 
                    stringsAsFactors = F, header = T)
release$SiteCode <- "Release"
release$Week <- date_to_studyweek(parse_date(release$Date, format = "%m/%d/%Y"))
release$Detected <- TRUE
release$last_detect = lubridate::mdy(release$Date)

## Format the data and turn into a weekly detection history ===============
# Extract the start of each bookmark into a single table
all_bkm_week <- all_bkm %>%
  mutate(Week = date_to_studyweek(as_date(CorrLocal_DateTime)),
         WeekEnd = date_to_studyweek(as_date(CorrLocal_DateTime)))

# Create a vector of all study weeks, all tags, and all hydrophones
Weeks <- seq(1,257, by = 1)
ReleasedTags <- release %>%
  filter(Species == "Striped Bass") %>%
  distinct(TagID) %>% arrange() %>% pull()
DetectedTags <- all_bkm_week %>%
  distinct(TagID) %>% arrange() %>% pull()

# Cross all the vectors and name the columns correctly
DateHydTagMatrix<- crossing(Weeks, hydrophones, "TagID" = DetectedTags)
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


rm(all_bkm_week, hydrophones, Hydrophones, HydBeaconBin, HydBeaconMatrix, 
   DateHydTagMatrix)
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

### Add in the Release Data
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
         BeaconDetected, AnyDetected) %>%
  filter(Species == "Striped Bass")

### Summarise Tag Release Data --------------------------------------------
release_summary <- release %>%
  filter(Species == "Striped Bass") %>%
  group_by(Species, Week) %>%
  tally()

releases <- ggplot()+
  geom_bar(data = release_summary,aes(x = studyweek_startdate(Week), y = n), 
           stat = "identity", width = 5)+
  scale_x_date(breaks = "3 months", expand = c(0,0))+
  facet_grid(Species~., scales = "free_y")+
  theme_classic()+
  theme(
    axis.text = element_text(angle = -45, vjust = 0.5, hjust = 0)
  )+
  labs(x = "Date", y = "Total Number of Released Individuals")+
  theme_classic()

png(file = file.path("1. Data","Figures","ReleaseSummary.png"), width = 15, height = 5,
    units = "in", res = 300)
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
    SiteCode %in% c("GL1") ~ "Grant Line Canal",
    SiteCode %in% c("CLRS") ~ "Curtis Landing")) %>%
  filter(SiteCode != "Tag Failure")

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
  select(Species, Location, duration) %>%
  group_by(Species, Location) %>%
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

(ShedTags_plot <- ggplot()+
    geom_boxplot(data = Site_CRE_Quantiles, 
                 aes(x = Location, ymin = qnt_10, lower = qnt_25, 
                     middle = qnt_50, upper = qnt_75, ymax=qnt_95), 
                 stat = "identity", lwd = 0.5)+
    geom_jitter(data = Site_CRE[Site_CRE$Species != "Beacon",] %>%
                  left_join(Site_CRE_Quantiles), 
                aes(x = Location, y = duration, shape = duration > qnt_99,
                    color = duration > qnt_99), 
                alpha = 0.4, width = 0.25, height = 0.01)+
    scale_color_manual(name = "Estimated Fate",
                       breaks = c(TRUE, FALSE),
                       values = c("red","blue"),
                       labels = c("Shed","Not Shed"))+
    scale_shape_manual(name = "Estimated Fate",
                       breaks = c(TRUE, FALSE),
                       values = c(4,1),
                       labels = c("Shed","Not Shed"))+
    scale_y_log10(name = "Uninterrupted Detection Periods (weeks)")+
    ggtitle("Analysis of Uninterrupted Detection Periods \nat a Single Receiver")+
    theme_classic())
  
png(file = file.path("1. Data", "Figures","ShedTagReview.png"),
    width = 5, height = 4.5, units ="in", res = 300)
ShedTags_plot
dev.off()

ShedTags <- Site_CRE %>%
  left_join(.,Site_CRE_Quantiles)%>%
  mutate(Shed = ifelse(duration > qnt_99,TRUE,FALSE)) %>%
  filter(Shed == T) %>%
  group_by(TagID)%>%
  summarise(ShedWeek = min(StartWeek)+1)

TagBKM_Bin %>%
  left_join(ShedTags) %>%
  mutate(ShedWeek = ifelse(is.na(ShedWeek),257,ShedWeek)) %>%
  filter(Week < ShedWeek) -> TagBKM_Bin

## Estimating TagLife =====================================================
# Tag life based on maximum Tag Life for specific tag type
TagTypes <- read.csv(file.path("1. Data","Inputs","TagTypes.csv"), header = T, 
                     stringsAsFactors = F)
TagLife <- TagBKM_Bin %>% 
  left_join(., TagTypes, by = "TagID") %>%
  mutate(Date = studyweek_startdate(Week)) %>%
  group_by(TagID, TagType) %>%
  summarise(TagLife_obs = as.numeric(max(Date)-min(Date))) %>%
  group_by(TagType) %>%
  summarise(med_TagLife_obs = median(TagLife_obs),
            mean_TagLife_obs = mean(TagLife_obs),
            max_TagLife_obs = max(TagLife_obs)) %>% #Which summary should we use? 
  mutate(man_TagLife = case_when(
    TagType == "LG" ~ as.integer(median(c(220,400))), # median of range from 220 to 400 days
    TagType == "LY" ~ as.integer(median(c(912,1460))), # median of range from 2.5 to 4 years
    TagType == "LZ" ~ as.integer(median(c(1460,1825))),# median of range from 4 to 5 years
    is.na(TagType) ~ as.integer(median(c(912,1460,1825))) # All NA tags were from years where LY and LZ tags used
  ))

TagLife_obs <- TagBKM_Bin %>%
  group_by(TagID) %>%
  mutate(Date = studyweek_startdate(Week)) %>%
  summarise(TagLife_obs = as.numeric(max(Date)-min(Date)))
  

TagFail <- release %>%
  filter(Species == "Striped Bass") %>%
  select(TagID, Species, Date) %>%
  left_join(.,TagTypes, "TagID") %>%
  mutate(Date = parse_date(Date, format = "%m/%d/%Y")) %>%
  mutate(TagType = ifelse(is.na(TagType),"LG",TagType)) %>%
  left_join(.,TagLife_obs, "TagID") %>%
  left_join(.,TagLife, "TagType") %>%
  mutate(sel_TagLife = case_when(
    TagLife_obs > man_TagLife ~ TagLife_obs, # Observed tag life exceeds manufacturer estimate
    TagLife_obs <= man_TagLife ~ man_TagLife, # Observed tag life is less than manufacturer estimate
  )) %>%
  mutate(faildate = Date + days(sel_TagLife)) %>%
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
png(file = file.path("1. Data", "Figures","ExampleTagDetection.png"),
    width = 10, height = 4.5, units ="in", res = 300)
(ggplot(TagBKM_Bin[TagBKM_Bin$TagID %in% c(5004.24) & 
                    TagBKM_Bin$SiteCode != "Tag Failure",])+
  geom_point(aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
  geom_step(aes(x = studyweek_startdate(Week), y = SiteCode, group = TagID))+
  geom_hline(aes(yintercept = 5.5), color = "blue", linetype = "dashed", 
             linewidth = 1)+
  scale_y_discrete()+
  scale_x_date(breaks = "1 month",
               date_labels = "%b %d")+
  theme(panel.spacing = unit(0, "lines"))+
  labs(x = "Date", y = "Site")+
  theme_classic())
dev.off()

## Redo weekly detection histories including tag failures =================
# Determine how many and which sites a tag visited on a single week
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
  filter(Species == "Striped Bass") %>%
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
  filter(!(grepl("Tag Failure",SiteVisits)) & !(grepl("Release",SiteVisits))) %>%
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

no_detect_tags <- ReleasedTags[!(ReleasedTags %in% DetectedTags)]
saveRDS(no_detect_tags,file.path("1. Data","Outputs","no_detect_tags.rds"))

# Remove Fish with no detections
WeeklySiteVisit <- WeeklySiteVisit %>%
filter(!(TagID %in% no_detect_tags)) 

### Code Location as either Inside, Outside, or Transiting ----------------

WeeklySiteVisit$Location <- lapply(WeeklySiteVisit$SiteVisits, 
                             function(x) case_when(
                               # Weeks where all detections are "Inside"
                               all(x %in% inside) ~ "INSIDE",
                               # Weeks where all detections are "Outside"
                               all(x %in% outside) ~ "OUTSIDE",
                               # Weeks with detections at the Intake Canal 
                               # AND beyond the radial gates are "FULL TRANSIT"
                               any(x %in% c("Release", "IC1", "IC2", "IC3")) &
                                 any(x %in% c("ORS1", "WC1", "WC2", "WC3", 
                                              "ORN1", "ORN2", "ORS3", "GL1", 
                                              "CVP1", "CLRS")) ~ "FULL TRANSIT",
                               # Weeks where the last detections are "Inside"
                               last(unlist(x)) %in% c("Release", "IC1",
                                                      "IC2", "IC3") ~ "INSIDE",
                               # Weeks where the last detections are "Outside"
                               last(unlist(x)) %in% c("ORS1", "WC1", "WC2", 
                                                      "WC3", "ORN1", "ORN2", 
                                                      "ORS3", "GL1", "CVP1", 
                                                      "CLRS") ~ "OUTSIDE",
                               # Any week where ALL detections are "Tag Failure"
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
saveRDS(WeeklySiteVisit, file.path("1. Data","Outputs","WeeklySiteVisit.rds"))
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
  filter(!(TagID %in% no_detect_tags)) %>%
  group_by(TagID) %>%
  arrange(TagID, StartWeek) %>%
  mutate(gap = lead(StartWeek)-EndWeek-1) %>%
  #mutate(gap = ifelse(is.na(gap),0,gap))
  mutate(
    StartWeek = case_when(
      #If there isn't a transit, just add time to the next location
      !grepl("TRANSIT", Location) &
        is.finite(lag(gap)) & lag(gap) > 0 ~ StartWeek - lag(gap), 
      #If there is a transit at both locations, do not extend
      grepl("TRANSIT", lag(Location)) & 
        grepl("TRANSIT", Location) ~ StartWeek,
      #If there is a transit at one location, extend the non-transit location
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

# Add in the Age Analysis data ############################################
# Get the Data from the Age Analysis
AssignedAge <- readRDS(file.path("1. Data","Outputs","Assigned_Ages.rds"))

CCF_Residency_CRE <- CCF_Residency_CRE %>%
  left_join(select(AssignedAge, TagID, EstCaptureAge), by = "TagID")

# Cross each tag with each week (universal dataframe)
tag_x_weeks <- crossing("TagID" = DetectedTags,Week = Weeks)

## Determine the Age of an Individual at Week n and Weeks at Age a.
weekly_tag_age <- CCF_Residency_CRE %>%
  filter(Species == "Striped Bass") %>%
  group_by(TagID) %>%
  summarise(CaptureWeek = min(StartWeek),
            CaptureAge = min(EstCaptureAge),
            EndWeek = max(EndWeek)) %>%
  left_join(tag_x_weeks) %>%
  mutate(time_from_capture = Week - CaptureWeek) %>%
  filter(time_from_capture >= 0 & Week <= EndWeek) %>%
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
  left_join(weekly_tag_age, by = c("TagID","StartWeek" = "Week")) %>%
  arrange(TagID, StartWeek)

saveRDS(CCF_Residency_CRE, file.path("1. Data","Outputs","CCF_Residency_CRE.rds"))
saveRDS(TagBKM_Bin, file.path("1. Data", "Outputs", "TagBKM_Bin.rds"))
saveRDS(weekly_tag_age, file.path("1. Data","Outputs","weekly_tag_age.rds"))
saveRDS(inside, file.path('1. Data',"Outputs","inside_sites.rds"))
saveRDS(outside, file.path('1. Data',"Outputs","outside_sites.rds"))

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
CCF_Residency_CRE %>% 
  filter(EstCaptureAge > 0 & !grepl("TRANSIT",Location)) %>% 
  group_by(AgeBin, Location) %>%
  summarise(n = n(),
            pos = -20)-> res_tally
  
cre_plot <- ggplot(CCF_Residency_CRE[CCF_Residency_CRE$EstCaptureAge > 0,] %>%
                     filter(Species == "Striped Bass" & !grepl("TRANSIT",Location))
                   )+
  geom_boxplot(aes(x = Location,
                   y = Duration_Weeks,color = Location, group = Location),
               fill = "grey95", size = 1)+
  geom_text(data = res_tally, aes(x = Location, y = pos, group = Location,
                label = paste("n = ", n, sep = "")))+
  facet_grid(Species~AgeBin)+
  theme_classic()+
  theme(
    #axis.text.x = element_text(angle = -45, hjust = 0)
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
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
    # If current week is a Full Transit and both before and after are the same,
    # then assume a complete entry and exit or exit and entry
    Location %in% c("FULL TRANSIT") & 
      lag(Location) == lead(Location) ~ "FULL TRANSIT",
    # If its a Partial Transit and the both before and after are the same, 
    # then do not assume an actual transit, unresolved. 
    Location %in% c("PARTIAL TRANSIT") & 
      lag(Location) == lead(Location) ~ "UNRESOLVED TRANSIT",
    # If none of the above are true and there is a Full Transit, and the next 
    # location is outside, "EXIT"
    Location == "FULL TRANSIT" & 
      lead(Location) == "OUTSIDE" ~ "EXIT",
    # If none of the above are true and there is a Full Transit, and the next 
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
# Location 0 is "New Tag" if there is not movement
CRE_Movement$move_direction <- ifelse(
  CRE_Movement$move_direction == "UNRESOLVED TRANSIT" & 
    CRE_Movement$LocationChange == 0,
  "New Tag", 
  CRE_Movement$move_direction)
# If there is an NA, or an Unresolved transit and current tag is different from 
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
                   group = interaction(AgeBin, Movement)),
               size = 1, fill = "grey95")+
  facet_grid(Species ~.)+
  theme_classic()
ggsave(file.path("1. Data", "Figures","Movements_by_Age.png"), device = "png",
       width = 6, height = 4)

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
ggsave(file.path("1. Data", "Figures","TransitFrequency_by_Age.png"), device = "png",
       width = 6, height = 4)

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

autoplot(km_fit)+theme_classic()

## Next try to separate by Species
km_fit_sp <- survival::survfit(survival::Surv(time, emmigrant)~Species, data = cre_surv[cre_surv$Species != "Beacon",])
autoplot(km_fit_sp)+theme_classic()

png(file.path("1. Data", "Figures", "TimeToEmmigration_All.png"), width = 10, 
    height =7, units = "in", res = 300)
autoplot(km_fit_sp)
dev.off()

## Fit our Age Bins for Striped Bass
cre_surv_age <- cre_surv %>%
  filter(Species == "Striped Bass")

km_fit_age <- survival::survfit(survival::Surv(time, emmigrant)~AgeBin, data = cre_surv_age)
autoplot(km_fit_age)+theme_classic()

png(file.path("1. Data", "Figures", "TimeToEmmigration_StripedBass_x_Age.png"), 
    width = 12, height =5,units = "in", res = 300)
autoplot(km_fit_age)+theme_classic()+
  labs(title = "Time to First Exit", y = "Percent of Population Non-Emmigrated",
       x = "Time (weeks)")
dev.off()

saveRDS(tbl_split_weeks,file.path("1. Data","Outputs","WeeklyTable.rds"))
saveRDS(CRE_Movement,file.path("1. Data","Outputs","Movement_CREs.rds"))
saveRDS(emmigrant_tags, file.path("1. Data","Outputs","emmigrants.rds"))
saveRDS(emmigration_table, file.path("1. Data","Outputs","emmigrant_table.rds"))
