############  FIOB Demographics Analysis Pull Data     ################

# Description
# Demographics from CCFPS, PRES, PFRS, EPFRRS

# Author:
# Taylor Spaulding
# tspaulding@esassoc.com

source(file.path("2. Code","0_Setup.R"))

# Bring in DELVE Data ######################################################
## CCFPS DELVE Project 5 ===================================================

# Fish Data
delver::get_dataset_file(project_id = 5,
                         dataset_id = 55,
                         version = 4,
                         filename = "ccfps_organism-collection_predators_20210507.csv"
) %>%
  fix_names() %>% 
  distinct(across(c(common_name:weight_g, cap_method))) %>%
  mutate(date = date(mdy_hm(sampling_datetime)),
         fork_length_mm = round(ifelse(length > 10, length*10, length*100)),
         weight_g = ifelse(weight <25, round(weight*453.592,2),
                           round(weight*45.3592,2)),
         sample_method = case_when(
           cap_method == "HL" ~ "hook-line",
           cap_method == "E" ~ "boat-electrofishing",
           cap_method == "GN" ~ "gill-net",
           TRUE ~ "hook-line")) %>%
  select(date, sample_method, common_name, fork_length_mm, weight_g) ->
  cfps_individual_final

ggplot(cfps_individual_final)+
  geom_point(aes(x = fork_length_mm, y = weight_g, color = common_name))

## PRES DELVE Project 8 ===================================================

delver::get_dataset_file(project_id = 8,
                         dataset_id = 120,
                         version = 2,
                         filename = "PRES_organism-sample_2024-01-16.csv") %>%
  fix_names() %>%
  mutate(date = date(sampling_datetime),
         sample_method = "boat-electrofishing",
         weight_g = ifelse(weight_g > 10000, weight_g/1000,
                           ifelse(weight_g > 4000, weight_g/4, weight_g))) %>%
  select(date, sample_method, common_name, fork_length_mm, total_length_mm, 
         weight_g) -> pres_individual_final

ggplot(pres_individual_final)+
  geom_point(aes(x = fork_length_mm, y = weight_g, color = common_name))

## PFRS DELVE Project 5 ===================================================
# This data not on DELVE

read.csv(file.path("1. Data","Inputs","pfrs_processing_data_20240118.csv")) %>%
  fix_names() %>% 
  mutate(date = mdy(date),
         common_name = case_when(
           species == "SB" ~ "striped-bass",
           species == "LMB" ~ "largemouth-bass",
           TRUE ~ species
         ),
         sample_method = case_when(
           gear.method == "Fyke" ~ "fyke-trap",
           gear.method == "Hoop" ~ "hoop-trap",
           gear.method == "Lampara" ~ "lampara",
           gear.method == "Kodiak" ~ "kodiak-trawl",
           gear.method == "Seine" ~ "beach-seine",
           TRUE ~ NA
         ),
         fork_length_mm = `fl..mm.`,
         weight_g = as.numeric(`ww..kg.`) * 1000) %>%
  select(date, sample_method, common_name, fork_length_mm, 
         weight_g) -> pfrs_individual_final

ggplot(pfrs_individual_final)+
  geom_point(aes(x = fork_length_mm, y = weight_g, color = common_name))


## EPFRRS DELVE Project 75 ===================================================

epfrrs_individual_fish <- delver::get_dataset_file(project_id = 75,
                                                   dataset_id = 222,
                                                   version = 1,
                                                   filename = "FishMeasurements.csv") %>%
  fix_names()

epfrrs_count_fish <- delver::get_dataset_file(project_id = 75,
                                              dataset_id = 223,
                                              version = 1,
                                              filename = "Processing_Counts.csv") %>%
  fix_names()

epfrrs_seine <- delver::get_dataset_file(project_id = 75,
                                         dataset_id = 224,
                                         version = 1,
                                         filename = "BeachSeine_Metadata.csv") %>%
  fix_names() %>%
  mutate(sample_method = "beach-seine")


epfrrs_efish <- delver::get_dataset_file(project_id = 75,
                                         dataset_id = 224,
                                         version = 1,
                                         filename = "Efishing_Transects_Metadata.csv") %>%
  fix_names()

epfrrs_hl <- delver::get_dataset_file(project_id = 75,
                                      dataset_id = 224,
                                      version = 1,
                                      filename = "HookLine_Metadata.csv") %>%
  fix_names()

epfrrs_hoop_retrieve <- delver::get_dataset_file(project_id = 75,
                                                 dataset_id = 224,
                                                 version = 1,
                                                 filename = "HoopTrap_Retrieval_Metadata.csv") %>%
  fix_names() %>%
  mutate(sample_method = "hoop-trap")

epfrrs_kodiak <- delver::get_dataset_file(project_id = 75,
                                          dataset_id = 224,
                                          version = 1,
                                          filename = "KodiakTrawl_Metadata.csv") %>%
  fix_names()
# Join The Data ###########################################################

## EPFRRS =================================================================
epfrrs_hoop_retrieve %>%
  select(processing_record_id, sample_method) -> epfrrs_hoop_red

epfrrs_kodiak %>%
  select(processing_record_id, sample_method) -> epfrrs_kodiak_red

epfrrs_seine %>%
  select(processing_record_id, sample_method) -> epfrrs_seine_red

epfrrs_hl %>%
  select(hl_record_id, sample_method) -> epfrrs_hl_red

epfrrs_efish %>%
  select(processing_record_id, sample_method) -> epfrrs_efish_red

epfrrs_individual_fish %>%
  mutate(date = date(sampling_datetime)) %>%
  select(date, common_name, fork_length_mm, total_length_mm, processing_record_id, 
         hl_record_id) -> epfrrs_ind_red

epfrrs_count_fish %>%
  select(date, common_name, "tallied" = obs_value, processing_record_id) -> epfrrs_count_red

epfrrs_ind_red %>%
  left_join(epfrrs_hl_red, by = "hl_record_id", multiple = "first") %>%
  filter(!is.na(sample_method))-> epfrrs_hl_fish
epfrrs_ind_red %>%
  left_join(epfrrs_kodiak_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_kodiak_fish
epfrrs_ind_red %>%
  left_join(epfrrs_hoop_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_hoop_fish
epfrrs_ind_red %>%
  left_join(epfrrs_seine_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_seine_fish
epfrrs_ind_red %>%
  left_join(epfrrs_efish_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_efish_fish


epfrrs_count_red %>%
  left_join(epfrrs_kodiak_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_kodiak_count
epfrrs_count_red %>%
  left_join(epfrrs_hoop_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_hoop_count
epfrrs_count_red %>%
  left_join(epfrrs_seine_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_seine_count
epfrrs_count_red %>%
  left_join(epfrrs_efish_red, by = "processing_record_id", multiple = "first") %>%
  filter(!is.na(sample_method)) -> epfrrs_efish_count

epfrrs_individual_final <- bind_rows(epfrrs_hl_fish, epfrrs_kodiak_fish, epfrrs_hoop_fish,
                                     epfrrs_seine_fish, epfrrs_efish_fish) %>%
  select(-processing_record_id, -hl_record_id)

epfrrs_count_final <- bind_rows(epfrrs_kodiak_count, epfrrs_hoop_count,
                                epfrrs_seine_count, epfrrs_efish_count) %>%
  select(-processing_record_id) %>%
  filter(common_name %in% c("striped-bass","striped_bass")) %>%
  group_by(date, sample_method, common_name) %>%
  summarise(tallied = sum(tallied, na.rm = TRUE))

# Cleanup
rm(epfrrs_seine_fish,epfrrs_hoop_fish,epfrrs_kodiak_fish,epfrrs_hl_fish,
   epfrrs_hl_red,epfrrs_seine_red,epfrrs_kodiak_red,epfrrs_hoop_red,
   epfrrs_count_fish, epfrrs_individual_fish, epfrrs_ind_red, epfrrs_count_red,
   epfrrs_kodiak_count, epfrrs_hoop_count, epfrrs_seine_count, epfrrs_efish_count,
   epfrrs_efish, epfrrs_efish_fish, epfrrs_efish_red, epfrrs_hl, epfrrs_hoop_retrieve,
   epfrrs_kodiak, epfrrs_seine)


## Combine All Data ================================================================
bind_rows(cfps_individual_final, pres_individual_final,
          pfrs_individual_final, epfrrs_individual_final
) %>%
  arrange(date)-> combined_individuals

saveRDS(combined_individuals, file.path("1. Data", "Outputs", "combined_sampling.rds"))
saveRDS(epfrrs_count_final, file.path("1. Data","Outputs","tallied_epfrrs.rds"))
