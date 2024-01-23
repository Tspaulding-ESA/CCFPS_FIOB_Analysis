############  FIOB Demographics Analysis Pull Data     ################

# Description
# Demographics from CCFPS, PRES, PFRS, EPFRRS

# Author:
# Taylor Spaulding
# tspaulding@esassoc.com

source(file.path("2. Code","0_Setup.R"))

# Bring in the Combined Data from "Pull Demographic Data" script
combined_individuals <- readRDS(file.path("1. Data", "Outputs", 
                                          "combined_sampling.rds")) %>%
  mutate(
    doy = yday(date),
    year = year(date),
    common_name = case_when(
      common_name %in% c("SB","Striped Bass","striped_bass",
                         "striped-bass") ~ "striped_bass",
      TRUE ~ common_name
    ))

tallied_epfrrs <- readRDS(file.path("1. Data","Outputs","tallied_epfrrs.rds")) %>%
  mutate(
    doy = yday(date),
    year = year(date),
    common_name = case_when(
      common_name %in% c("SB","Striped Bass","striped_bass",
                         "striped-bass") ~ "striped_bass",
      TRUE ~ common_name
    ))

combined_individuals %>%
  filter(common_name == "striped_bass",
         !is.na(fork_length_mm)) %>%
  mutate(measurement = "measured")-> measured_sb

combined_individuals %>%
  filter(common_name == "striped_bass",
         is.na(fork_length_mm)) %>%
  group_by(year, doy, sample_method, common_name) %>%
  tally(name = "tallied") %>% ungroup() -> tallied_sb

tallied_sb <- bind_rows(tallied_sb, tallied_epfrrs) %>%
  group_by(year, doy, sample_method, common_name) %>%
  summarise(tallied = sum(tallied, na.rm = tallied))

ggplot(measured_sb)+
  geom_point(aes(x = fork_length_mm, y = weight_g))

# Begin Analysis ###########################################################

## Expanded Length Frequency ===============================================
# A large number of fish were not assigned lengths, expand the tallied
# individuals by the length frequency observed in the measured individuals

# Gerritsen and McGrath 2007 suggest 10x the number of size classes. 

# First create a dataframe that contains counts for the total number of measured
# and tallied individuals for each day for each method
measured_sb %>%
  ungroup() %>%
  group_by(year, doy, sample_method, common_name) %>%
  mutate(measured = n()) %>%
  left_join(select(tallied_sb, year, doy, sample_method, common_name, 
                   tallied)) %>%
  replace_na(list("tallied" = 0)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(total_count = measured + tallied,
         perc_measured = measured/total_count) %>%
  ungroup() %>%
  as_tibble() %>%
  select(year, doy, sample_method, common_name, fork_length_mm, weight_g, 
         measured, tallied, total_count, perc_measured) %>%
  arrange(year, doy)-> combined_sb

expanded_lengths_list <- list()
for(y in 1:length(unique(combined_sb$year))){
  for(d in 1:length(unique(combined_sb$doy))){
    for(s in 1:length(unique(combined_sb$sample_method))){
      tmp <- combined_sb %>%
        filter(year == unique(combined_sb$year)[y] & 
                 doy == unique(combined_sb$doy)[d] & 
                 sample_method == unique(combined_sb$sample_method)[s])
      if(nrow(tmp) > 0){
        year <- unique(combined_sb$year)[y]
        doy <- unique(combined_sb$doy)[d]
        method <- unique(combined_sb$sample_method)[s]
        lengths <- tmp %>% pull(fork_length_mm)
        tally <- tmp[1,] %>% pull(tallied)
        min_length = min(tmp$fork_length_mm)
        total <- tmp[1,] %>% pull(total_count)
        if(length(lengths) > 1 & 
           length(lengths) > tally*0.01 & 
           tally >= 1) {
          expanded <-  expandLenFreq(x = lengths,
                                      w = 5,
                                      additional = tally,
                                      startcat = min_length,
                                      total = total,
                                      decimals = 0)
        } else {
          expanded <- -9999
        }
        if(all(!is.null(expanded))){
          expanded_df <- data.frame("year" = year,
                                    "doy" = doy,
                                    "sample_method" = method,
                                    "common_name" = "striped_bass",
                                    "fork_length_mm" = expanded,
                                    "measurement" = "simulated")
        } else {
          expanded_df <- data.frame("year" = year,
                                    "doy" = doy,
                                    "sample_method" = method,
                                    "common_name" = "striped_bass",
                                    "fork_length_mm" = NA_real_,
                                    "measurement" = "simulated")
        }
        expanded_lengths_list[[paste(year,doy,method)]] <- expanded_df
        rm(year,doy, method, lengths, tally, min_length, total, expanded_df)
      }
    }
  }
}

expanded_lengths_list_pois <- list()
for(y in 1:length(unique(combined_sb$year))){
  for(d in 1:length(unique(combined_sb$doy))){
    for(s in 1:length(unique(combined_sb$sample_method))){
      tmp <- combined_sb %>%
        filter(year == unique(combined_sb$year)[y] & 
                 doy == unique(combined_sb$doy)[d] & 
                 sample_method == unique(combined_sb$sample_method)[s])
      if(nrow(tmp) > 0){
        year <- unique(combined_sb$year)[y]
        doy <- unique(combined_sb$doy)[d]
        method <- unique(combined_sb$sample_method)[s]
        lengths <- tmp %>% pull(fork_length_mm)
        tally <- tmp[1,] %>% pull(tallied)
        min_length = min(tmp$fork_length_mm)
        total <- tmp[1,] %>% pull(total_count)
        if(length(lengths) > 1 & 
           length(lengths) > tally*0.01 & 
           tally >= 1) {
          expanded <-  expandLenFreq2(x = lengths,
                                     w = 1,
                                     additional = tally,
                                     startcat = min_length,
                                     total = total,
                                     densfun = "Poisson",
                                     decimals = 0)
        } else {
          expanded <- -9999
        }
        if(all(!is.null(expanded))){
          expanded_df <- data.frame("year" = year,
                                    "doy" = doy,
                                    "sample_method" = method,
                                    "common_name" = "striped_bass",
                                    "fork_length_mm" = expanded,
                                    "measurement" = "simulated")
        } else {
          expanded_df <- data.frame("year" = year,
                                    "doy" = doy,
                                    "sample_method" = method,
                                    "common_name" = "striped_bass",
                                    "fork_length_mm" = NA_real_,
                                    "measurement" = "simulated")
        }
        expanded_lengths_list_pois[[paste(year,doy,method)]] <- expanded_df
        rm(year,doy, method, lengths, tally, min_length, total, expanded_df)
      }
    }
  }
}

expanded_lengths <- bind_rows(expanded_lengths_list) %>%
  filter(!is.na(fork_length_mm) & fork_length_mm != -9999) %>%
  arrange(year, doy)
expanded_lengths_pois <- bind_rows(expanded_lengths_list_pois) %>%
  filter(!is.na(fork_length_mm) & fork_length_mm != -9999) %>%
  arrange(year, doy)


expanded_sb <- bind_rows(measured_sb, expanded_lengths) %>%
  mutate(measurement = factor(measurement, levels = c("simulated","measured")),
         date = ymd(paste(year-1,12,31))+days(doy),
         season = factor(case_when( # Based on insolation (Solar Intensity)
           month(date) %in% c(11,12,1) ~ "Winter",
           month(date) %in% c(2,3,4) ~ "Spring",
           month(date) %in% c(5,6,7) ~ "Summer",
           month(date) %in% c(8,9,10) ~ "Fall"), levels = c("Winter", "Spring",
                                                            "Summer", "Fall")
         ))

expanded_sb_pois <- bind_rows(measured_sb, expanded_lengths_pois) %>%
  mutate(measurement = factor(measurement, levels = c("simulated","measured")),
         date = ymd(paste(year-1,12,31))+days(doy),
         season = factor(case_when( # Based on insolation (Solar Intensity)
           month(date) %in% c(11,12,1) ~ "Winter",
           month(date) %in% c(2,3,4) ~ "Spring",
           month(date) %in% c(5,6,7) ~ "Summer",
           month(date) %in% c(8,9,10) ~ "Fall"), levels = c("Winter", "Spring",
                                                            "Summer", "Fall")
         ))

# Find the samples which were not simulated
expanded_sb_pois %>%
  group_by(year, doy, sample_method, measurement) %>%
  tally() %>%
  arrange(year, doy, sample_method) %>%
  pivot_wider(names_from = measurement, values_from = n, values_fill = 0) -> sim_compare

# Check to make sure no fish came in above samples were from samples where fish 
# were tallied but not simulated
combined_sb %>%
  left_join(sim_compare) %>%
  select(year, doy, sample_method, measured, tallied, simulated) %>% 
  distinct() %>%
  mutate(perc_measured = measured/(measured+tallied),
         perc_tallied = tallied/(measured+tallied),
         perc_simulated = simulated/(measured+tallied)) %>%
  filter(tallied != simulated)

## Fork Length Distributions ===============================================
### Overall ----------------------------------------------------------------
ggplot(expanded_sb)+
  geom_histogram(aes(x = fork_length_mm, fill = measurement),
                 position = position_stack(), binwidth = 5)+
  theme_classic()+
  labs(x = "Forklength (mm)",
       y = "Catch Distribution")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  labs(title = "Distribution from multiplying existing LF Distribution")

ggplot(expanded_sb_pois)+
  geom_histogram(aes(x = fork_length_mm, fill = measurement),
                 position = position_stack(), binwidth = 5)+
  theme_classic()+
  labs(x = "Forklength (mm)",
       y = "Catch Distribution")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  labs(title = "Distribution from sampling the Poisson Fit of the existing LF Distribution")

### By Gear Type -----------------------------------------------------------

expanded_sb_pois %>%
  group_by(sample_method, measurement) %>%
  tally() %>%
  ungroup() %>%
  pivot_wider(names_from = measurement, values_from = `n`) %>%
  ungroup() %>%
  bind_rows(data.frame(sample_method = "TOTAL",
                       simulated = sum(.$simulated, na.rm = TRUE),
                       measured = sum(.$measured, na.rm = TRUE))) %>%
  mutate(perc_measured = measured/(measured+simulated))

ggplot(expanded_sb_pois %>%
         filter())+
  geom_histogram(aes(x = fork_length_mm, fill = measurement),
               position = position_stack(), binwidth = 5)+
  theme_classic()+
  labs(x = "Forklength (mm)",
       y = "Catch Distribution")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  facet_grid(sample_method ~ ., scales = "free_y")
    
### By Time of Year --------------------------------------------------------
ggplot(expanded_sb_pois)+
  geom_histogram(aes(x = fork_length_mm, fill = measurement),
                 position = position_stack(), binwidth = 5)+
  theme_classic()+
  labs(x = "Forklength (mm)",
       y = "Catch Distribution")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  facet_grid(season ~ ., scales = "free_y")

# Age Analysis #############################################################

## Get the Original Analyis from CFPS ======================================
#source(file.path("2. Code","Age_Analysis.R"))
stb.key <- read_rds(file.path("1. Data","Outputs","SB_ALKey.rds"))

## Assign Ages to all Striped Bass using age-length key ====================
stb.len <- expanded_sb %>%
  filter(common_name == "striped_bass") %>%
  mutate(fork_length_cm = fork_length_mm/10) %>%
  filter(fork_length_cm > 12)

expanded_sb <- FSA::alkIndivAge(key = stb.key, 
                             formula = ~fork_length_cm, 
                             data = stb.len, 
                             type = "SR", seed = 42069) %>%
  mutate(age_class = factor(case_when(
    age %in% c(2) ~ "0 - 2 yrs",
    age %in% c(3:5) ~ "3 - 5 yrs",
    age %in% c(6:7) ~ "6+ yrs"),
    levels = c("0 - 2 yrs","3 - 5 yrs","6+ yrs")),
    month = factor(month.abb[month(date)], levels = month.abb))
FSA::alkPlot(stb.key)

## Age Analysis ============================================================ 

### Age Distribution by Date -----------------------------------------------
ggplot(expanded_sb)+
  geom_bar(aes(x = month, group = age_class, fill = age_class), stat = "count",
           just = 0.5, position = position_fill(reverse = TRUE))+
  scale_fill_discrete(name = "Age Class")+
  scale_y_continuous(name = "Percent of Catch",
                     labels = scales::percent)+
  theme_classic()

### Age Distribution by Season ----------------------------------------------
ggplot(expanded_sb)+
  geom_bar(aes(x = season, group = age_class, fill = age_class), stat = "count",
           just = 0.5, position = position_fill(reverse = TRUE))+
  scale_fill_discrete(name = "Age Class")+
  scale_y_continuous(name = "Percent of Catch",
                     labels = scales::percent)+
  theme_classic()



