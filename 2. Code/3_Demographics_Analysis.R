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

# expanded_lengths_list <- list()
# for(y in 1:length(unique(combined_sb$year))){
#   for(d in 1:length(unique(combined_sb$doy))){
#     for(s in 1:length(unique(combined_sb$sample_method))){
#       tmp <- combined_sb %>%
#         filter(year == unique(combined_sb$year)[y] & 
#                  doy == unique(combined_sb$doy)[d] & 
#                  sample_method == unique(combined_sb$sample_method)[s])
#       if(nrow(tmp) > 0){
#         year <- unique(combined_sb$year)[y]
#         doy <- unique(combined_sb$doy)[d]
#         method <- unique(combined_sb$sample_method)[s]
#         lengths <- tmp %>% pull(fork_length_mm)
#         tally <- tmp[1,] %>% pull(tallied)
#         min_length = min(tmp$fork_length_mm)
#         total <- tmp[1,] %>% pull(total_count)
#         if(length(lengths) > 1 & 
#            length(lengths) > tally*0.01 & 
#            tally >= 1) {
#           expanded <-  expandLenFreq(x = lengths,
#                                       w = 5,
#                                       additional = tally,
#                                       startcat = min_length,
#                                       total = total,
#                                       decimals = 0)
#         } else {
#           expanded <- -9999
#         }
#         if(all(!is.null(expanded))){
#           expanded_df <- data.frame("year" = year,
#                                     "doy" = doy,
#                                     "sample_method" = method,
#                                     "common_name" = "striped_bass",
#                                     "fork_length_mm" = expanded,
#                                     "measurement" = "simulated")
#         } else {
#           expanded_df <- data.frame("year" = year,
#                                     "doy" = doy,
#                                     "sample_method" = method,
#                                     "common_name" = "striped_bass",
#                                     "fork_length_mm" = NA_real_,
#                                     "measurement" = "simulated")
#         }
#         expanded_lengths_list[[paste(year,doy,method)]] <- expanded_df
#         rm(year,doy, method, lengths, tally, min_length, total, expanded_df)
#       }
#     }
#   }
# }
# 
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

# expanded_lengths <- bind_rows(expanded_lengths_list) %>%
#   filter(!is.na(fork_length_mm) & fork_length_mm != -9999) %>%
#   arrange(year, doy)
expanded_lengths_pois <- bind_rows(expanded_lengths_list_pois) %>%
  filter(!is.na(fork_length_mm) & fork_length_mm != -9999) %>%
  arrange(year, doy)


# expanded_sb <- bind_rows(measured_sb, expanded_lengths) %>%
#   mutate(measurement = factor(measurement, levels = c("simulated","measured")),
#          date = ymd(paste(year-1,12,31))+days(doy),
#          season = factor(case_when( # Based on insolation (Solar Intensity)
#            month(date) %in% c(11,12,1) ~ "Winter",
#            month(date) %in% c(2,3,4) ~ "Spring",
#            month(date) %in% c(5,6,7) ~ "Summer",
#            month(date) %in% c(8,9,10) ~ "Fall"), levels = c("Winter", "Spring",
#                                                             "Summer", "Fall")
#          ))

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
# ggplot(expanded_sb)+
#   geom_histogram(aes(x = fork_length_mm, fill = measurement),
#                  position = position_stack(), binwidth = 5)+
#   theme_classic()+
#   labs(x = "Forklength (mm)",
#        y = "Catch Distribution")+
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )+
#   labs(title = "Distribution from multiplying existing LF Distribution")

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
stb.len <- expanded_sb_pois %>%
  filter(common_name == "striped_bass") %>%
  mutate(fork_length_cm = fork_length_mm/10) %>%
  filter(fork_length_cm > 12)

expanded_sb_pois <- FSA::alkIndivAge(key = stb.key, 
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
ggplot(expanded_sb_pois)+
  geom_bar(aes(x = month, group = age_class, fill = age_class), stat = "count",
           just = 0.5, position = position_fill(reverse = TRUE))+
  scale_fill_discrete(name = "Age Class")+
  scale_y_continuous(name = "Percent of Catch",
                     labels = scales::percent)+
  theme_classic()

### Age Distribution by Season ----------------------------------------------
ggplot(expanded_sb_pois)+
  geom_bar(aes(x = season, group = age_class, fill = age_class), stat = "count",
           just = 0.5, position = position_fill(reverse = TRUE))+
  scale_fill_discrete(name = "Age Class")+
  scale_y_continuous(name = "Percent of Catch",
                     labels = scales::percent)+
  theme_classic()

# Population Health #########################################################

## Relative Weight =========================================================
# Get output from the Bayes analysis of Fishbase data on Striped Bass
bayes_out <- read.csv(file.path("1. Data","Outputs","morone_saxatilis.csv"))
bayes_b <- bayes_out$mean_b
bayes_a <- bayes_out$true_geom_mean_a

# calculate the Relative Weight and Expected Weight for each fish
expanded_sb_pois |>
  mutate(fsa_name = "Striped Bass") |>
  mutate(logW = log10(weight_g),
         logL = log10(fork_length_mm),
         K = weight_g/(fork_length_mm^3)*100000,
         wr_fsa = wrAdd(weight_g, fork_length_mm, fsa_name)/100,
         bayes_ws = bayes_a*((fork_length_mm/10)^bayes_b),
         fsa_ws = weight_g/wr_fsa,
         wr_bayes = weight_g/bayes_ws) -> wr


### Length-Weight Regression Plot -------------------------------------------

# Regress the relationship
m = nls(weight_g ~ a*fork_length_mm^b, start = list(a =0.0000007,b=3), 
        data = wr);

# Create Equation Labels
eq1_label = quote(plain("Derived Equation: "))
eq1 = substitute( italic(y) == a  ~italic(x)^b, 
                           list(a = format(coef(m)[[1]], digits = 5), 
                                b = format(coef(m)[[2]], digits = 5)))
eq2_label = quote(plain("Standard Equation: "))
eq2 = substitute( italic(y) == a  ~italic(x)^b, 
                           list(a = format(bayes_a, digits = 5), 
                                b = format(bayes_b, digits = 5)))

eq1 <- substitute(eq1_label * eq1, list(eq1_label = eq1_label,eq1 = eq1))
eq2 <- substitute(eq2_label * eq2, list(eq2_label = eq2_label, eq2 = eq2))
eq <- c(as.character(as.expression(eq1)), as.character(as.expression(eq2)))

eq_label = data.frame(x = c(quantile(wr$fork_length_mm, 0.50),
                            quantile(wr$fork_length_mm, 0.50)),
                      y = c(quantile(wr$weight_g, 0.9999, na.rm = TRUE),
                            quantile(wr$weight_g, 0.9999, na.rm = TRUE)),
                      label = as.character(eq))
# Plot 
condition_plot <- ggplot()+
  geom_point(data = measured_sb, aes(x = fork_length_mm, y = weight_g))+
  # geom_point(data = wr, aes(x = fork_length_mm, y = bayes_ws), 
  #            inherit.aes = FALSE, color = "red")+
  stat_smooth(data = measured_sb, aes(x = fork_length_mm, y = weight_g),
              method = 'nls', formula = 'y~a*x^b', 
              method.args=list(start = list(a = 1,b=1), 
                               control=nls.control(maxiter=16000)),se=FALSE,
              color = "blue")+
  geom_line(data = wr %>% distinct(fork_length_mm, bayes_ws) %>% arrange(), 
            aes(x = fork_length_mm, y = bayes_ws),
            color = "red", linewidth = 1)+
  scale_x_continuous(expand=c(0,0), limits = c(150,1100))+
  scale_y_continuous(expand=c(0,0), limits = c(0,15000))+
  labs(x = "Length (mm)", y= "Weight (g)")+
  theme_classic()+
  coord_cartesian(clip = 'off')+
  geom_text(data = eq_label,
            aes(x, y, label = label),
            color = c("blue","red"),
            parse = TRUE, inherit.aes = FALSE,
            size = 5, nudge_y = c(0,1000))
condition_plot

### Relative Weight Point Plot ---------------------------------------------
ggplot(wr)+
  # geom_point(aes(x=fork_length_mm, y= wr_fsa),shape = 21, 
  #            color = "red", alpha = 0.1)+
  geom_point(aes(x=fork_length_mm, y= wr_bayes),shape = 21, 
             color = "blue", alpha = 0.1)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  lims(y = quantile(wr$wr_fsa,c(.01,.99), na.rm = TRUE),
       x = quantile(wr$fork_length_mm, c(0,1), na.rm = TRUE))+
  theme_classic()

### Relative Weight Distribution by Size Class (bin = 10mm) ----------------
ggplot(wr %>%
         mutate(size_class = floor(fork_length_mm/10)))+
  geom_boxplot(aes(x = size_class, y = wr_bayes, 
                  group = factor(size_class), 
                  fill = factor(size_class)),
               outlier.alpha = 0,
               alpha = 0.25)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "tomato",
             size = 1.5)+
  geom_smooth(data = wr %>%
              mutate(size_class = floor(fork_length_mm/10)) %>%
                group_by(size_class) %>%
                dplyr::summarise(mean_wr = mean(wr_bayes, na.rm = TRUE)),
              aes(x = size_class, y = mean_wr), 
              method = "gam", 
              formula = y ~ s(x, bs = "tp", k = 21), 
              se = FALSE,
              linetype = "dashed",
              size = 1.5,
              color = "blue")+
  coord_cartesian(ylim = c(0,2))+
  theme_classic()

## Proportional Stock Density

