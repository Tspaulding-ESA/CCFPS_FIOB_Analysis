
fix_names <- function(data){
  new_names = gsub("-|\\s+", "_", tolower(names(data)))
  setNames(data, new_names)
}

from_canonical <- function(x){
  stringr::str_to_title(gsub("-", " ", x))
}

create_links <- function(data, parent, child){
  # assumes that parent-child will either be one-to-one or one-to-many
  # i.e., parent column will contain delimited lists
  out = list()
  for (i in seq_along(data[[parent]])){
    child_split = trimws(strsplit(data[[child]][i], ",")[[1]])
    out[[i]] = data.frame(parent = data[[parent]][i],
                          child = child_split)
  }
  bind_rows(out)
}

# Define a function to fix the bbox to be in EPSG:26910
ggmap_bbox <- function(map, crs) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  # Coonvert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_trans <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), crs))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_trans["ymin"]
  attr(map, "bb")$ll.lon <- bbox_trans["xmin"]
  attr(map, "bb")$ur.lat <- bbox_trans["ymax"]
  attr(map, "bb")$ur.lon <- bbox_trans["xmax"]
  map
}

#function to turn a date into a week of the study
date_to_studyweek <- function (date) {
  studyweek <- floor(as.numeric(date - as_date("2013-05-03"))/7)+1
  return(studyweek)
}

#function to find the start date of a study week
studyweek_startdate <- function (Week) {
  startdate <- as_date("2013-05-03")+weeks(Week-1)
  return(startdate)
}

#Function to find the end date of a study week
studyweek_enddate <- function (Week) {
  enddate <- as_date("2013-05-09")+weeks(Week-1)
  return(enddate)
}

# Create helper function to split the interval of the CRE into a given unit size (ie weeks)
# Remember that the original weekly bin was only for detections, we are now
# assuming that a fish is remaining in a location and need to code each week
# it is at-large as such
split_interval_weeks <- function(start, end, unit) {
  breaks <- seq(floor_date(start, "day"), ceiling_date(end, "day"), by = unit)
  timeline <- c(start, breaks[breaks > start & breaks < end], end)
  tibble(.start = head(timeline, -1), .end = tail(timeline, -1))
}

