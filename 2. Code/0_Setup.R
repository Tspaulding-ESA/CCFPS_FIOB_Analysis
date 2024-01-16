options(dplyr.summarise.inform = FALSE)
library(lubridate)
library(tidyverse)
library(delver)
library(sf)

# delver::set_user_token([insert user token here])

source(file.path("2.Code", "functions.R"))
tz_loc = "America/Los_Angeles"
