############  FIOB Demographics Analysis Demographics      ################

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
  filter(common_name == "striped-bass")

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

combined_individuals %>%
  mutate(common_name = case_when(
    common_name %in% c("SB","Striped Bass","striped_bass",
                       "striped-bass") ~ "striped_bass",
    TRUE ~ common_name
  )) %>%
  filter(common_name == "striped_bass",
         !is.na(fork_length_mm)) %>%
  mutate(measurement = "measured") -> measured_sb

combined_individuals %>%
  mutate(common_name = case_when(
    common_name %in% c("SB","Striped Bass","striped_bass",
                       "striped-bass") ~ "striped_bass",
    TRUE ~ common_name
  )) %>%
  filter(common_name == "striped_bass",
         is.na(fork_length_mm)) %>%
  group_by(date, sample_method, common_name) %>%
  tally(name = "tallied") %>% ungroup() -> tallied_sb

tallied_sb <- bind_rows(tallied_sb, epfrrs_count_final)

rm(epfrrs_individual_final, pres_individual_final, 
   pfrs_individual_final, cfps_individual_final)

ggplot(measured_sb)+
  geom_point(aes(x = fork_length_mm, y = weight_g))

# Begin Analysis ###########################################################

## Expanded Length Frequency
# A large number of fish were not assigned lengths, expand the tallied
# individuals by the length frequency observed in the measured individuals

measured_sb %>%
  ungroup() %>%
  group_by(date, sample_method, common_name) %>%
  mutate(measured = n()) %>%
  left_join(select(tallied_sb, date, sample_method, common_name, 
                   tallied)) %>%
  replace_na(list("tallied" = 0)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(total_count = measured + tallied,
         perc_measured = measured/total_count) %>%
  ungroup() %>%
  as_tibble() %>%
  select(date, sample_method, common_name, fork_length_mm, weight_g, 
         measured, tallied, total_count, perc_measured) %>%
  arrange(date) %>%
  mutate(doy = yday(date),
         year = year(date))-> combined_sb



expanded_lengths_list <- list()
for(y in unique(combined_sb$year)){
  for(d in unique(combined_sb$doy)){
    for(s in unique(combined_sb$sample_method)){
      tmp <- combined_sb %>%
        filter(year == y & doy == d & sample_method == s)
      lengths <- tmp %>% pull(fork_length_mm)
      tally <- tmp[1,] %>% pull(tallied)
      min_length = min(tmp$fork_length_mm)
      total <- tmp[1,] %>% pull(total_count)
      tryCatch(
        {
          expanded <-  FSA::expandLenFreq(x = lengths,
                                          w = 1,
                                          additional = tally,
                                          startcat = min_length,
                                          total = total,
                                          decimals = 0)
          expanded_df <- data.frame("year" = rep(y,times = length(expanded)),
                                    "doy" = rep(d,times = length(expanded)),
                                    "sample_method" = rep(s, times = length(expanded)),
                                    "common_name" = rep("striped_bass", times = length(expanded)),
                                    "fork_length_mm" = expanded,
                                    "measurement" = rep("simulated", times = length(expanded)))
        },
        error = function(e){
          expanded_df <- data.frame("year" = y,
                                    "doy" = d,
                                    "sample_method" = s,
                                    "common_name" = "striped_bass",
                                    "fork_length_mm" = NA_real_,
                                    "measurement" = "no simulation")
        }
      )
      expanded_lengths_list[[paste(y,d,s)]] <- expanded_df
    }
  }
}


expanded_sb <- bind_rows(measured_sb, expanded_lengths) %>%
  mutate(measurement = factor(measurement, levels = c("simulated","measured")))

# Find the fish which were not simulated
expanded_sb %>%
  group_by(date, measurement) %>%
  tally() %>%
  arrange(date) %>%
  pivot_wider(names_from = measurement, values_from = n, values_fill = 0) -> sim_compare

combined_sb %>%
  left_join(sim_compare) %>%
  filter(tallied != simulated) %>% View()

## Fork Length Distributions ===============================================
### Overall ----------------------------------------------------------------
ggplot(expanded_sb)+
  geom_density(aes(x = fork_length_mm, fill = measurement), 
               position = position_stack())+
  geom_rug(data=subset(expanded_sb, 
                       measurement == "measured"),
           aes(x = fork_length_mm, color = measurement)) +
  geom_rug(data=subset(expanded_sb,
                       measurement == "simulated"),
           aes(x = fork_length_mm, color = measurement), sides="t")+
  theme_classic()+
  facet_grid(measurement ~ .)+
  labs(x = "Forklength (mm)",
       y = "Catch Distribution")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

### By Gear Type -----------------------------------------------------------

expanded_sb %>%
  group_by(sample_method, measurement) %>%
  tally() %>%
  pivot_wider(names_from = measurement, values_from = `n`) %>%
  mutate(perc_measured = measured/(measured+simulated))
  

ggplot(expanded_sb %>%
         filter())+
  geom_histogram(aes(x = fork_length_mm, fill = measurement),
               position = position_stack(), binwidth = 5)+
  geom_rug(data=subset(expanded_sb, 
                       measurement == "measured"),
           aes(x = fork_length_mm, color = measurement)) +
  geom_rug(data=subset(expanded_sb,
                       measurement == "simulated"),
           aes(x = fork_length_mm, color = measurement), sides="t")+
  theme_classic()+
  labs(x = "Forklength (mm)",
       y = "Catch Distribution")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  facet_grid(sample_method ~ ., scales = "free_y")
    
### By Time of Year --------------------------------------------------------
ggplot(expanded_sb)+
  geom_density(aes(x = fork_length_mm, fill = measurement),
               position = position_stack())+
  geom_rug(data=subset(expanded_sb, 
                       measurement == "measured"),
           aes(x = fork_length_mm, color = measurement)) +
  geom_rug(data=subset(expanded_sb,
                       measurement == "simulated"),
           aes(x = fork_length_mm, color = measurement), sides="t")+
  theme_classic()+
  labs(x = "Forklength (mm)",
       y = "Catch Distribution")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  facet_grid(factor(case_when(
    week %in% c(1:8,48:53) ~ "Winter",
    week %in% c(9:21) ~ "Spring",
    week %in% c(22:34) ~ "Summer",
    week %in% c(35:47) ~ "Fall"),
    levels = c("Fall","Winter","Spring","Summer")) ~ ., scales = "free_y")

## Age Analysis from CFPS --------------------------------------------------
#source(file.path("2. Code","Age_Analysis.R"))
stb.key <- read_rds(file.path("1. Data","Outputs","SB_ALKey.rds"))



### Assign Ages to Striped Bass at capture using age-length key ------------
stb.len <- measured_sb %>%
  filter(common_name == "striped_bass") %>%
  mutate(fork_length_cm = fork_length_mm/10) %>%
  filter(fork_length_cm > 12)

stb.len1 <- FSA::alkIndivAge(key = stb.key, 
                             formula = ~fork_length_cm, 
                             data = stb.len, 
                             type = "SR", seed = 42069)
FSA::alkPlot(stb.key)

