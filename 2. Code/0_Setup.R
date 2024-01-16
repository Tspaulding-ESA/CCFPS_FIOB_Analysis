options(dplyr.summarise.inform = FALSE)
library(lubridate)
library(tidyverse)
library(delver)
library(sf)

token = read_file(file.path("1. Data","Inputs","DelveToken.txt"))

delver::set_user_token(token)

source(file.path("2. Code", "functions.R"))
tz_loc = "America/Los_Angeles"
