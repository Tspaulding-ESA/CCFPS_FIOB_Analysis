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
  season = c("Immigration_a","Spawn","Emmigration","Residence","Immigration_b"),
  start =  c(ymd("2013-01-01"),ymd("2013-03-15"),ymd("2013-05-15"),ymd("2013-08-01"), ymd("2013-11-01")),
  end = c(ymd("2013-03-14"),ymd("2013-05-14"),ymd("2013-07-31"),ymd("2013-10-31"),ymd("2013-12-31"))
  ) %>%
  mutate(season_duration = date_to_studyweek(end)-date_to_studyweek(start),
         start = yday(start),
         end = yday(end)) %>%
  mutate(season_sub = stringr::str_remove(season,"_.")) %>%
  group_by(season_sub) %>%
  mutate(season_duration = sum(season_duration, na.rm = TRUE))

#
tbl_split_weeks %>%
  select(TagID, LocationChange, Location, AgeBin, time_at_age, move_direction, .start, .end) %>%
  mutate(doy = yday(.start),
         Week = date_to_studyweek(.start))%>%
  left_join(seasons, by = join_by(between(x$doy, y$start, y$end))) %>%
  mutate(season = case_when(
    season %in% c("Immigration_a","Immigration_b") ~ "Immigration",
    TRUE ~ season)
    )-> season_week_split

WeeklySiteVisit %>%
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  left_join(seasons, by = join_by(between(x$doy, y$start, y$end))) %>%
  mutate(season = case_when(
    season %in% c("Immigration_a","Immigration_b") ~ "Immigration",
    TRUE ~ season)) -> season_site_visits

TagBKM_Bin %>%
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  left_join(seasons, by = join_by(between(x$doy, y$start, y$end)))%>%
  mutate(season = case_when(
    season %in% c("Immigration_a","Immigration_b") ~ "Immigration",
    TRUE ~ season)) -> season_bkm_bin

## Total Number of Receivers Detected per year per season
season_site_visits %>%
  full_join(dplyr::select(season_week_split,TagID,Week, season, season_duration)) %>%
  ungroup() %>%
  left_join(weekly_tag_age) %>%
  filter(!(grepl("Tag Failure",SiteVisits) & NumSites == 1)) %>%
  filter(!(grepl("Release",SiteVisits) & NumSites == 1)) %>%
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  mutate(year = as.numeric(year(studyweek_startdate(Week)))) %>%
  mutate(year = ifelse(season == "Immigration" & doy > 300,
                       year + 1, year)) %>%
  group_by(TagID) %>%
  mutate(first_year = min(year, na.rm = TRUE),
         first_season = first(season)) %>% 
  mutate(season_duration = ifelse(year == first_year & season == first_season,
                                  floor((end-doy)/7)+1,
                                  season_duration)) %>%
  group_by(TagID, AgeBin, season,year) %>%
  summarise(
    season_duration = case_when(
      season == first_season & year == first_year ~ max(season_duration, na.rm = TRUE),
      TRUE ~ max(season_duration, na.rm = TRUE)),
    inside_sites_total = 
      sum(unlist(lapply(SiteVisits,
                        function(x){unlist(x, recursive = TRUE) %in% inside})
                 )==TRUE),
    outside_sites_total =  
      sum(unlist(lapply(SiteVisits,
                        function(x){unlist(x, recursive = TRUE) %in% outside})
                 )==TRUE),
    all_sites_total = inside_sites_total + outside_sites_total) %>% 
  group_by(TagID, AgeBin, year, season) %>% 
  mutate(inside_site_freq = inside_sites_total/season_duration,
         outside_site_freq = outside_sites_total/season_duration,
         all_site_freq = all_sites_total/season_duration) %>%
  ungroup() %>% select(-season_duration) %>%
  replace_na(list(inside_site_freq = 0,
                  outside_site_freq = 0,
                  all_site_freq = 0))-> seasonal_site_visit_sum

## Receivers/Wk detected On ================================================
WkDetectedRcvrs <- season_site_visits %>%
  full_join(dplyr::select(season_week_split,TagID,Week, season)) %>%
  ungroup() %>%
  left_join(weekly_tag_age) %>%
  filter(!(grepl("Tag Failure",SiteVisits) & NumSites == 1)) %>%
  filter(!(grepl("Release",SiteVisits) & NumSites == 1)) %>%
  mutate(NumSites = case_when(
    grepl("Tag Failure", SiteVisits) ~ NumSites-1,
    grepl("Release", SiteVisits) ~ NumSites-1,
    TRUE ~ NumSites))  %>%
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  mutate(year = as.numeric(year(studyweek_startdate(Week)))) %>%
  mutate(year = ifelse(season == "Immigration" & doy > 300,
                       year + 1, year)) %>%
  group_by(TagID, AgeBin, year, season) %>%
  summarise(weekly_receivers_mean = round(mean(NumSites, na.rm = TRUE),0),
            weekly_receivers_max = max(NumSites, na.rm = TRUE)) %>%
  mutate(weekly_receivers_mean = ifelse(is.na(weekly_receivers_mean),0,weekly_receivers_mean),
         weekly_receivers_max = ifelse(!is.finite(weekly_receivers_max),0,weekly_receivers_max))

## Weeks Undetected ========================================================
Wks_Undetected <- season_bkm_bin %>%
  filter(SiteCode != "Tag Failure") %>%
  ungroup() %>%
  full_join(select(season_week_split, -AgeBin)) %>%
  left_join(weekly_tag_age, by = c("TagID","Week")) %>%
  group_by(TagID, AgeBin, season, season_duration, Week) %>%
  summarise(SiteCode = max(SiteCode, na.rm = TRUE)) %>%
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  mutate(year = as.numeric(year(studyweek_startdate(Week)))) %>%
  mutate(year = ifelse(season == "Immigration"& doy > 300,
                       year + 1, year)) %>% 
  group_by(TagID, year, season, season_duration, AgeBin) %>%
  mutate(
    detection_change = coalesce(!is.na(SiteCode), TRUE) %>% 
      cumsum()
  ) %>% 
  mutate(SiteCode = !is.na(SiteCode)) %>%
  group_by(TagID,  AgeBin, year, season, season_duration, SiteCode, detection_change) %>%
  tally(name = "time_btwn") %>%
  mutate(time_btwn = case_when(
    time_btwn > season_duration ~ season_duration,
    time_btwn == 1 & SiteCode ~ 0,
    TRUE~time_btwn-1)) %>%
  ungroup() %>%
  group_by(TagID,AgeBin, year, season, season_duration) %>%
  summarise(time_btwn_detects_max = max(time_btwn, na.rm = TRUE),
            time_btwn_detects_total = sum(time_btwn,na.rm = TRUE))

## Weeks btwn Movements ====================================================
season_week_split %>%
  select(TagID, AgeBin, season, season_duration, LocationChange, move_direction, .start) %>%
  filter(!(move_direction %in% c("New Tag","NOT TRANSIT"))) %>%
  mutate(year = ifelse(season == "Immigration" & yday(.start) > 303, 
                       year(.start)+1,
                       year(.start))) %>%
  group_by(TagID, AgeBin, year, season, season_duration, LocationChange, move_direction) %>%
  summarise(Week = min(date_to_studyweek(.start))) %>%
  group_by(TagID, AgeBin, year, season, season_duration) %>%
  arrange(TagID, Week, year) %>%
  mutate(time_btwn_movements = Week - lag(Week)) %>%
  group_by(TagID, AgeBin, year, season, season_duration) %>% 
  summarise(time_btwn_mvmts_max = max(time_btwn_movements, na.rm = TRUE),
            time_btwn_mvmts_mean = round(mean(time_btwn_movements,na.rm = TRUE),0)) %>%
  mutate(across(time_btwn_mvmts_max:time_btwn_mvmts_mean, 
                ~ifelse(is.na(.x)|!is.finite(.x),season_duration,.x))) -> time_btwn_movements

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
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  mutate(year = as.numeric(year(studyweek_startdate(Week)))) %>%
  mutate(year = ifelse(season == "Immigration" & doy > 303,
                       year + 1, year)) %>%
  group_by(TagID, AgeBin, year, season) %>%
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
  mutate(woy = lubridate::week(studyweek_startdate(Week)),
         year = year(studyweek_startdate(Week)),
         doy = yday(studyweek_startdate(Week))) %>%
  left_join(seasons, by = join_by(between(x$Week, y$start, y$end))) %>%
  mutate(year = ifelse(season_sub == "Immigration" & doy > 300,
                       year + 1, year)) %>%
  group_by(TagID, year, "season" = season_sub) %>%
  summarise(mean_exit_woy = mean(woy))

## Avg Week of Year of Entry ===============================================
avg_entry_wk <- emmigration_table %>%
  filter(Movement == "Entry") %>%
  select(TagID, Week, jdate, month) %>%
  mutate(woy = lubridate::week(studyweek_startdate(Week)),
         year = year(studyweek_startdate(Week)),
         doy = yday(studyweek_startdate(Week))) %>%
  left_join(seasons, by = join_by(between(x$Week, y$start, y$end))) %>%
  mutate(year = ifelse(season_sub == "Immigration" & doy > 300,
                       year + 1, year)) %>%
  group_by(TagID, year, "season" = season_sub) %>%
  summarise(mean_entry_woy = mean(woy))

## Number of Exits and Entries =============================================
season_week_split %>%
  select(TagID, AgeBin, season, season_duration, LocationChange, move_direction,
         start, end,.start) %>%
  mutate(doy = yday(.start)) %>%
  mutate(year = as.numeric(year(.start))) %>%
  mutate(year = ifelse(season == "Immigration" & doy > 300,
                       year + 1, year)) %>%
  group_by(TagID) %>%
  mutate(first_year = min(year),
         first_season = season[which(.start == min(.start))]) %>%
  group_by(year, TagID) %>%
  mutate(season_duration = ifelse(year == first_year & season == first_season,
                                  ceiling((end-yday(.start))/7),
                                  season_duration)) %>%
  filter(!(move_direction %in% c("New Tag","NOT TRANSIT","UNRESOLVED TRANSIT"))) %>%
  group_by(year, TagID, AgeBin, season) %>%
  mutate(season_duration = max(season_duration, na.rm = TRUE)) %>% 
  distinct(TagID, AgeBin, season, season_duration, LocationChange, move_direction) %>%
  mutate(count = case_when(
    move_direction == "ENTRY AND EXIT x 2" ~ 4,
    move_direction == "ENTRY AND EXIT" ~ 2,
    TRUE ~ 1))%>%
  pivot_wider(names_from = "move_direction",values_from = "count", values_fill = 0) %>%
  mutate(
    entries = case_when(
      `ENTRY AND EXIT` == 2 ~ ENTRY + 1,
      TRUE ~ ENTRY),
    exits = case_when(
      `ENTRY AND EXIT` == 2 ~ EXIT + 1,
      TRUE ~ EXIT),
    total_mvmts = entries+exits) %>%
  select(-c(ENTRY,EXIT,`ENTRY AND EXIT`)) %>%
  ungroup() %>%
  group_by(TagID, AgeBin, year, season) %>%
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
  mutate(doy = yday(studyweek_startdate(Week))) %>%
  mutate(year = as.numeric(year(studyweek_startdate(Week)))) %>%
  mutate(year = ifelse(season == "Immigration" & doy > 303,
                       year + 1, year)) %>%
  group_by(year, TagID, AgeBin, season) %>%
  filter(!(Location %in% c("UNRESOLVED","FULL TRANSIT"))) %>%
  mutate(Location = case_when(
    grepl("NewTag", Location) ~ "INSIDE",
    grepl("INSIDE",Location) ~ "INSIDE",
    grepl("OUTSIDE",Location) ~ "OUTSIDE",
  )) %>%
  group_by(year, TagID, AgeBin, season, LocationChange, Location) %>%
  summarise(residence = n()) %>% 
  ungroup() %>%
  group_by(TagID, AgeBin, year, season, Location) %>%
  summarise(residence_total = sum(residence, na.rm = TRUE),
            residence_mean = mean(residence, na.rm = TRUE),
            residence_max = max(residence, na.rm = TRUE)) %>%
  pivot_wider(names_from = Location,
              values_from = residence_total:residence_max,
              names_sep = "_",
              values_fill = 0)-> residence

## Gate Operations and Water Year Type =====================================
gate_ops <- readRDS(file.path("1. Data","Outputs","radial_gate_summary.rds"))

season_week_split %>%
  mutate(year = year(.start),
         week = date_to_studyweek(.start)) %>% 
  mutate(doy = yday(.start)) %>%
  mutate(year = ifelse(season == "Immigration" & doy > 303,
                       year + 1, year)) %>%
  ungroup() %>%
  left_join(gate_ops, by = "week") %>% 
  select(TagID, AgeBin, year, season, week, wyt, perc_open) %>%
  filter(!is.na(wyt)) %>% # Few couple weeks are missing data
  ungroup() %>% 
  group_by(TagID, AgeBin, year, season) %>%
  summarise(wyt = getmode(wyt),
            gate_perc_open = mean(perc_open, na.rm = TRUE)) -> gate_ops

## Join the Data
dat <- seasonal_site_visit_sum %>%
  left_join(WkDetectedRcvrs) %>%
  left_join(Wks_Undetected) %>%
  left_join(time_btwn_movements) %>%
  left_join(Distances) %>%
  left_join(avg_exit_wk) %>%
  left_join(avg_entry_wk) %>%
  left_join(movement_metrics) %>%
  left_join(residence) %>%
  left_join(gate_ops, by = c("TagID","AgeBin", "season", "year")) %>%
  filter(!(TagID %in% no_detect_tags)) %>% # Remove tags with no detections
  replace_na(list(TagID = NA,
                  time_btwn_detects_total = 100,
                  time_btwn_detects_max = 100,
                  time_btwn_mvmts_mean = 100,
                  time_btwn_mvmts_max = 100,
                  mean_distance = 0,
                  q05_distance = 0,
                  q25_distance = 0, 
                  q50_distance = 0, 
                  q75_distance = 0, 
                  q99_distance = 0,
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
  mutate(ID = paste(TagID, AgeBin, season, year, sep = "-"),
         #Species = ifelse(Species == "Striped Bass",1,2),
         AgeBin = as.numeric(factor(AgeBin, levels = c("0","1-2","3-5","6+"))),
         wyt = as.numeric(wyt),
         season = as.numeric(factor(season, levels = c("Immigration","Spawn","Emmigration","Residence"))),
         time_btwn_detects_max = ifelse(time_btwn_detects_max == 100,
                                         season_duration,
                                         time_btwn_detects_max),
         time_btwn_detects_total = ifelse(time_btwn_detects_total == 100,
                                         season_duration,
                                         time_btwn_detects_total),
         time_btwn_mvmts_max = ifelse(time_btwn_mvmts_max == 100,
                                        season_duration,
                                        time_btwn_mvmts_max),
         time_btwn_mvmts_mean = ifelse(time_btwn_mvmts_mean == 100,
                                         season_duration,
                                         time_btwn_mvmts_mean)
         
  ) %>%
  distinct() %>%
  column_to_rownames(var = "ID") %>%
  select(-TagID, -season_duration, -year)

saveRDS(dat, file.path("1. Data", "Outputs" , "Cluster_Analysis_Data.rds"))

