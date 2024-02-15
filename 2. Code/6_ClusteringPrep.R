# Prepare the data for Cluster analysis ###################################
source(file.path("2. Code","0_Setup.R"))

WeeklySiteVisit <- readRDS(file.path("1. Data","Outputs","WeeklySiteVisit.rds"))
tbl_split_weeks <- readRDS(file.path("1. Data","Outputs","WeeklyTable.rds"))
CRE_Movement <- readRDS(file.path("1. Data","Outputs","Movement_CREs.rds"))
CCF_Residency_CRE <- readRDS(file.path("1. Data","Outputs","CCF_Residency_CRE.rds"))
TagBKM_Bin <- readRDS(file.path("1. Data", "Outputs", "TagBKM_Bin.rds"))
stb_mig <- read.csv(file.path("1. Data", "Inputs", "striped_bass_migration.csv"))
weekly_tag_age <- readRDS(file.path("1. Data","Outputs","weekly_tag_age.rds"))
AssignedAge <- readRDS(file.path("1. Data","Outputs","Assigned_Ages.rds"))
inside <- readRDS(file.path('1. Data',"Outputs","inside_sites.rds"))
outside <- readRDS(file.path('1. Data',"Outputs","outside_sites.rds"))
emmigrant_tags  = readRDS(file.path("1. Data","Outputs","emmigrants.rds"))
emmigration_table = readRDS(file.path("1. Data","Outputs","emmigrant_table.rds"))
no_detect_tags = readRDS(file.path("1. Data","Outputs","no_detect_tags.rds"))

# Migration Season #########################################################
# Generalised from data from Calhoun 1952, Radtke 1966, Sabal 2018, and Goertler 2021
seasons <- data.frame(
  season = c("Immigration","Spawn","Emmigration","Residence","Immigration"),
  start =  c(ymd("2013-01-01"),ymd("2013-03-15"),ymd("2013-05-15"),ymd("2013-08-01"), ymd("2013-11-01")),
  end = c(ymd("2013-03-14"),ymd("2013-05-14"),ymd("2013-07-31"),ymd("2013-10-31"),ymd("2013-12-31"))
  ) %>%
  mutate(season_duration = date_to_studyweek(end)-date_to_studyweek(start),
         start = yday(start),
         end = yday(end))

#
tbl_split_weeks %>%
  select(TagID, LocationChange, Location, AgeBin, time_at_age, move_direction, .start, .end) %>%
  mutate(doy = yday(.start)) %>%
  left_join(seasons, by = join_by(between(x$doy, y$start, y$end))) -> season_week_split

WeeklySiteVisit %>%
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  left_join(seasons, by = join_by(between(x$doy, y$start, y$end))) -> season_site_visits

TagBKM_Bin %>%
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  left_join(seasons, by = join_by(between(x$doy, y$start, y$end))) -> season_bkm_bin

## Total Number of Receivers Detected per year per season
season_site_visits %>%
group_by(TagID) %>%
  ungroup() %>%
  left_join(weekly_tag_age) %>%
  filter(!(grepl("Tag Failure",SiteVisits) & NumSites == 1)) %>%
  filter(!(grepl("Release",SiteVisits) & NumSites == 1)) %>%
  mutate(NumSites = case_when(
    grepl("Tag Failure", SiteVisits) ~ NumSites-1,
    grepl("Release", SiteVisits) ~ NumSites-1,
    TRUE ~ NumSites)) %>%
  mutate(year = year(studyweek_startdate(Week))) %>%
  group_by(TagID) %>%
  mutate(first_year = min(year),
         first_season = first(season)) %>%
  group_by(year, TagID) %>%
  mutate(season_duration = ifelse(year == first_year & season == first_season,
                                  ceiling((end-doy)/7)+1,
                                  end-start)) %>% 
  group_by(year, TagID, AgeBin, season) %>%
  summarise(
    season_duration = max(season_duration),
    inside_sites = sum(
      unlist(
        lapply(SiteVisits,function(x){unlist(x, 
                                             recursive = TRUE) %in% inside}))),
    outside_sites = sum(
      unlist(
        lapply(SiteVisits,function(x){unlist(x, 
                                             recursive = TRUE) %in% outside}))),
    all_sites = inside_sites+outside_sites) %>%
  group_by(TagID, AgeBin, season) %>% 
  summarise(season_duration = sum(season_duration),
            inside_sites_total = sum(inside_sites),
            outside_sites_total = sum(outside_sites),
            all_sites_total = sum(all_sites)) %>%
  mutate(inside_site_freq = inside_sites_total/season_duration,
         outside_site_freq = outside_sites_total/season_duration,
         all_site_freq = all_sites_total/season_duration) %>%
  ungroup() %>% select(-season_duration) -> seasonal_site_visit_sum

## Receivers/Wk detected On ================================================
WkDetectedRcvrs <- season_site_visits %>%
  ungroup() %>%
  left_join(weekly_tag_age) %>%
  filter(!(grepl("Tag Failure",SiteVisits) & NumSites == 1)) %>%
  filter(!(grepl("Release",SiteVisits) & NumSites == 1)) %>%
  mutate(NumSites = case_when(
    grepl("Tag Failure", SiteVisits) ~ NumSites-1,
    grepl("Release", SiteVisits) ~ NumSites-1,
    TRUE ~ NumSites)) %>%
  group_by(TagID, AgeBin, season) %>%
  summarise(weekly_receivers_mean = round(mean(NumSites, na.rm = TRUE),0),
            weekly_receivers_max = max(NumSites, na.rm = TRUE))

## Estimated Capture Age ===================================================
CaptureAge <- AssignedAge %>%
  select(TagID, EstCaptureAge) %>%
  mutate(CaptureAgeBin = case_when(
    between(EstCaptureAge,0,2) ~ "1-2",
    between(EstCaptureAge,3,5) ~ "3-5",
    EstCaptureAge >= 6 ~ "6+"
  )) %>%
  select(TagID, CaptureAgeBin)

# CaptureLength <- release %>%
#   select(TagID, Length_cm) %>%
#   mutate(length_class = (floor(Length_cm/5)*5)) %>%
#   select(TagID, length_class)

## Weeks Undetected ========================================================
Wks_Undetected <- season_bkm_bin %>%
  filter(SiteCode != "Tag Failure") %>%
  ungroup() %>%
  left_join(weekly_tag_age, by = c("TagID","Week")) %>%
  group_by(TagID, AgeBin, season, Week) %>%
  distinct(TagID, AgeBin, season, Week) %>%
  group_by(TagID, AgeBin) %>%
  mutate(time_btwn = Week - lag(Week)) %>% 
  filter(!is.na(time_btwn)) %>%
  group_by(TagID,AgeBin, season) %>%
  summarise(time_btwn_detects_max = max(time_btwn, na.rm = TRUE),
            time_btwn_detects_mean = round(mean(time_btwn,na.rm = TRUE),0))

## Weeks btwn Movements ====================================================
season_week_split %>%
  select(TagID, AgeBin, season, LocationChange, move_direction, .start) %>%
  filter(!(move_direction %in% c("New Tag","NOT TRANSIT"))) %>%
  group_by(TagID, AgeBin, season, LocationChange, move_direction) %>%
  summarise(Week = min(date_to_studyweek(.start))) %>%
  group_by(TagID, AgeBin) %>%
  arrange(TagID, Week) %>%
  mutate(time_btwn_movements = Week - lag(Week)) %>%
  group_by(TagID, AgeBin, season) %>% 
  summarise(time_btwn_mvmts_max = max(time_btwn_movements, na.rm = TRUE),
            time_btwn_mvmts_mean = round(mean(time_btwn_movements,na.rm = TRUE),0)) %>%
  mutate(across(time_btwn_mvmts_max:time_btwn_mvmts_mean, 
                ~ifelse(is.na(.x)|!is.finite(.x),100,.x))) -> time_btwn_movements

## Quantiles of Distance per Week ==========================================
### Assign minimum distance travelled in one week --------------------------

#Read in Linear Distance
distances <- read.csv(file.path("1. Data", "Inputs", "HydDistance.csv")) 

#Filter Detections
season_site_visits %>%
  left_join(weekly_tag_age) %>%
  select(TagID, AgeBin, season, Week, "SiteA" = SiteVisits) %>%
  unnest(cols = SiteA) %>%
  filter(!(grepl("Tag Failure",SiteA)) & !(grepl("Release",SiteA))) %>%
  left_join(season_site_visits %>% 
              select(TagID, season, Week, "SiteB" = SiteVisits) %>%
              unnest(cols = SiteB) %>% 
              filter(!(SiteB %in% c("Release",NA,"Tag Failure")))) %>%
  left_join(distances, relationship = "many-to-many") %>%
  # If one of the visits is CLRS but the visit before was inside (especially IC)
  # Likely fish was transported to CLRS not just detected. Only one fish appears 
  # To have swum out the RGS, passed GL1 and then been detected on CLRS
  mutate(Distance = case_when(
    SiteA == "CLRS" & SiteB %in% inside ~ 0, 
    TRUE ~ Distance
  )) %>%
  group_by(season, Week, TagID, AgeBin) %>%
  summarise(weekly_max_distance = max(Distance),
            weekly_mean_distance = mean(Distance)) -> WeeklyDistance

### Summarize Distance -----------------------------------------------------
WeeklyDistance %>%
  ungroup() %>%
  group_by(TagID, AgeBin, season) %>%
  summarise(mean_distance = mean(weekly_mean_distance, na.rm = TRUE),
            q05_distance = quantile(weekly_max_distance, 0.01, na.rm = TRUE),
            q25_distance = quantile(weekly_max_distance, 0.25, na.rm = TRUE),
            q50_distance = quantile(weekly_max_distance, 0.50, na.rm = TRUE),
            q75_distance = quantile(weekly_max_distance, 0.75, na.rm = TRUE),
            q99_distance = quantile(weekly_max_distance, 0.99, na.rm = TRUE)) -> Distances

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
season_week_split %>%
  select(TagID, AgeBin, season, season_duration, LocationChange, move_direction,
         start, end,.start) %>%
  mutate(year = year(.start)) %>%
  group_by(TagID) %>%
  mutate(first_year = min(year),
         first_date = min(.start)) %>%
  group_by(year, TagID) %>%
  mutate(season_duration = ifelse(year == first_year & .start == first_date,
                                  ceiling((end-yday(.start))/7),
                                  end-start)) %>%
  filter(!(move_direction %in% c("New Tag","NOT TRANSIT","UNRESOLVED TRANSIT", 
                                 "PARTIAL TRANSIT"))) %>%
  group_by(year, TagID, AgeBin, season) %>%
  mutate(season_duration = max(season_duration)) %>%
  distinct(TagID, AgeBin, season, season_duration, move_direction) %>%
  mutate(count = ifelse(move_direction == "FULL TRANSIT", 2,1))%>%
  pivot_wider(names_from = "move_direction",values_from = "count", values_fill = 0) %>%
  mutate(entries = ifelse(`FULL TRANSIT` == 2, ENTRY + 1, ENTRY),
         exits = ifelse(`FULL TRANSIT` == 2, EXIT + 1, EXIT),
         total_mvmts = entries+exits) %>%
  select(-c(ENTRY,EXIT,`FULL TRANSIT`)) %>%
  ungroup() %>%
  group_by(TagID, AgeBin, season) %>%
  summarise(entry_total = sum(entries),
            exit_total = sum(exits),
            all_total = sum(total_mvmts),
            season_duration = sum(season_duration)) %>%
  mutate(entry_freq = entry_total/season_duration,
            exit_freq = exit_total/season_duration,
            all_freq = all_total/season_duration) -> movement_metrics

## Residence Time ==================================================
season_week_split %>%
  select(TagID, AgeBin, season, season_duration, LocationChange, Location, end,.start) %>%
  mutate(year = year(.start),
         Week = date_to_studyweek(.start)) %>%
  group_by(year, TagID, AgeBin, season) %>%
  filter(Location %in% c("INSIDE","OUTSIDE")) %>%
  group_by(year, TagID, AgeBin, season, LocationChange, Location) %>%
  mutate(residence = n()) %>%
  ungroup() %>%
  distinct(TagID, AgeBin, season, Location, residence) %>%
  group_by(TagID, AgeBin, season, Location) %>%
  summarise(residence_total = sum(residence),
            residence_mean = mean(residence),
            residence_max = max(residence)) %>%
  pivot_wider(names_from = Location,
              values_from = residence_total:residence_max,
              names_sep = "_",
              values_fill = 0)-> residence

## Gate Operations and Water Year Type =====================================
gate_ops <- readRDS(file.path("1. Data","Outputs","radial_gate_summary.rds"))

season_week_split %>%
  mutate(week = date_to_studyweek(.start)) %>% 
  ungroup() %>%
  left_join(gate_ops, by = "week") %>% 
  select(TagID, AgeBin, season, week, wyt, perc_open) %>%
  filter(!is.na(wyt)) %>% # Few couple weeks are missing data
  ungroup() %>% 
  group_by(TagID, AgeBin, season) %>%
  summarise(wyt = getmode(wyt),
            gate_perc_open = mean(perc_open, na.rm = TRUE)) -> gate_ops

## Join the Data
dat <- seasonal_site_visit_sum %>%
  left_join(WkDetectedRcvrs) %>%
  left_join(CaptureAge) %>%
  left_join(Wks_Undetected) %>%
  left_join(time_btwn_movements) %>%
  left_join(Distances) %>%
  left_join(emmigrant_tags) %>%
  left_join(avg_exit_wk) %>%
  left_join(avg_entry_wk) %>%
  left_join(movement_metrics) %>%
  left_join(residence) %>%
  left_join(gate_ops, by = c("TagID","AgeBin", "season")) %>%
  filter(!(TagID %in% no_detect_tags)) %>% # Remove tags with no detections
  replace_na(list(TagID = NA,
                  CaptureAge = NA,
                  inside_sites_total = 0,
                  outside_sites_total = 0,
                  all_sites_total = 0,
                  inside_site_freq = 0,
                  outside_site_freq = 0,
                  all_site_freq = 0,
                  weekly_receivers_mean = 0,
                  weekly_receivers_max = 0,
                  time_btwn_detects_mean = 100,
                  time_btwn_detects_max = 100,
                  time_btwn_mvmts_mean = 100,
                  time_btwn_mvmts_max = 100,
                  mean_distance = 0,
                  q05_distance = 0,
                  q25_distance = 0, 
                  q50_distance = 0, 
                  q75_distance = 0, 
                  q99_distance = 0,
                  emmigrant = FALSE,
                  mean_exit_woy = 0, 
                  mean_entry_woy = 0,
                  entry_total = 0, 
                  exit_total = 0, 
                  all_total = 0,
                  entry_freq = 0,
                  exit_freq = 0,
                  all_freq = 0,
                  residence_total_INSIDE = 0,
                  residence_total_OUTSIDE = 0,
                  residence_mean_INSIDE = 0,
                  residence_mean_OUTSIDE = 0,
                  residence_max_INSIDE = 0,
                  residence_max_OUTSIDE = 0,
                  gate_perc_open = 0,
                  wyt = "Critical")) %>%
  mutate(ID = paste(TagID, AgeBin, season, sep = "-"),
         #Species = ifelse(Species == "Striped Bass",1,2),
         AgeBin = as.numeric(factor(AgeBin, levels = c("0","1-2","3-5","6+"))),
         CaptureAgeBin = as.numeric(factor(CaptureAgeBin, levels = c("0","1-2","3-5","6+"))),
         wyt = as.numeric(wyt),
         season = as.numeric(factor(season, levels = c("Immigration","Spawn","Emmigration","Residence")))
  ) %>%
  distinct() %>%
  filter(!is.na(CaptureAgeBin)) %>%
  column_to_rownames(var = "ID") %>%
  select(-TagID, -season_duration)

saveRDS(dat, file.path("1. Data", "Outputs" , "Cluster_Analysis_Data.rds"))

