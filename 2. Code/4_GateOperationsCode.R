#### GATE OPS FROM DWR FILE HANDLER ####
require(lubridate)
require(dplyr)

# Read in the Gate Open Data
dat_gates <- read.csv(file.path("1. Data","Inputs", "Gate_CLOSE_OPEN.csv"), 
                      stringsAsFactors = F, header = T) |>
  fix_names() %>%
  rename("closed" = certainclose,
         "open" = certainopen) %>%
  mutate(across(everything(),~lubridate::mdy_hms(., tz = tz_loc))) %>%
  mutate(date = as_date(closed),
         week = date_to_studyweek(date),
         day_begin = with_tz(ymd_hms(paste(date, "00:00:01"),
                                         tz = tz_loc), tzone = tz_loc),
         day_end = with_tz(ymd_hms(paste(date, "23:59:59"),
                                       tz = tz_loc),
                           tzone = tz_loc)) %>%
  mutate(next_closed = as_datetime(lead(closed), tz = tz_loc),
         last_opened = as_datetime(lag(open), tz = tz_loc)) %>%
  mutate(next_closed = as_datetime(ifelse(next_closed > day_end, 
                              day_end, 
                              next_closed),
                              tz = tz_loc),
         last_opened = as_datetime(ifelse(last_opened < day_begin, 
                              as_datetime(day_begin, tz = tz_loc),
                              as_datetime(last_opened, tz = tz_loc)),
                              tz = tz_loc)
  ) %>%
  mutate(cross_date = lead(date) != date)

# Calculate the various amounts of time the gates were open in a given day
dat_gates %>% 
  mutate(open_1 = ifelse(date(closed) == date(open),
                         ifelse(last_opened > day_begin & !is.na(last_opened),
                          difftime(closed, last_opened, units = "hours"),
                          difftime(closed, day_begin, units = "hours")),
                          NA),
         open_2 = ifelse(date(closed) != date(open),
                         ifelse(last_opened > day_begin,
                                difftime(closed, last_opened, units = "hours"),
                                NA),
                         NA),
         open_3 = ifelse(date(open) != lead(date),
                         difftime(day_end, open, units = "hours"),
                         NA)) %>%
  mutate(across(open_1:open_3, .fns = function(x) ifelse(x < 0, 0, x))) %>%
  mutate(total_open = rowSums(pick(open_1:open_3), na.rm = T)) %>%
  mutate(wy = ifelse(month(date) >  10, year(date) + 1, year(date)),
         wyt = esaRmisc::get_water_year_type(wy, valley = "SJR")) %>%
  group_by(wy, wyt, week, date) %>%
  summarise(total_open = sum(total_open),
            perc_open = total_open/24) -> gate_sum

saveRDS(gate_sum, file.path("1. Data","Outputs","radial_gate_summary.rds"))

gate_sum %>%
  ungroup() %>%
  group_by(wyt, wy) %>%
  summarise(perc_open = mean(perc_open))