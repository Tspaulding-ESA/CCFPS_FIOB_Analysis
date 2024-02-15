getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode2 <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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

#Recreate the FSA expand Len Freq function but with sampling from a poisson distribution
expandLenFreq2 <- function(x, w, additional, startcat = NULL, total = additional + 
                             length(x), decimals = decs$wdec, 
                           densfun = "Poisson",  show.summary = TRUE)
{
  if (!is.vector(x)) 
    STOP("'x' must be a vector.")
  if (!is.numeric(x)) 
    STOP("'x' must be numeric.")
  if (w <= 0) 
    STOP("'w' must be positive")
  if (!is.null(startcat)) 
    if (startcat <= 0) 
      STOP("'startcat' must be positive")
  if (total < length(x)) 
    STOP("Total number to expand to must be greater than number in 'x'.")
  if (is.null(startcat)) 
  startcat <- floor(min(x, na.rm = TRUE)/w) * w
  decs <- FSA:::iCheckStartcatW(startcat, w, x)
  num <- total - length(x)
  lcat <- lencat(x, w = w, startcat = startcat)
  lenfreq <- prop.table(table(lcat, dnn = NULL))
  lambda = MASS::fitdistr(lengths, densfun = densfun)$estimate
  cats <- seq(startcat,max(x),by = w)
  new.lenfreq <- numeric()
  for(i in 1:length(cats)){
    prob = ppois(cats[i], lambda = lambda) - ppois(cats[i-1], lambda = lambda)
    new.lenfreq[i] <- ifelse(!is_empty(prob),prob,0)
  }
  names(new.lenfreq) <- cats
  reps = floor(num*new.lenfreq)
  nrand.lens <- rep(cats, reps)
  rand.lens <- sample(cats, num - sum(reps), replace = TRUE, 
                      prob = new.lenfreq)
  new.lens <- c(nrand.lens, rand.lens)
  maxval <- w - 1/(10^decimals)
  if (maxval > 0) {
    new.lens <- new.lens + stats::runif(length(new.lens), 
                                        min = 0, max = maxval)
  }
  new.lens <- round(new.lens, decimals)
  if (show.summary) {
    cat("Length Frequency Expansion using:\n", "Measured length frequency of", 
        length(x), "individuals:\n")
    print(round(lenfreq, 4))
    cat("\nPoisson allocations of", length(nrand.lens), 
        "individuals by length category\n")
    tmp <- reps
    print(tmp)
    cat("\nRandom allocations of", length(rand.lens), "individuals from resampling\n", 
        "With final length frequency table of:\n")
    final.lens <- lencat(new.lens, w = w, startcat = startcat)
    print(table(final.lens, dnn = NULL))
  }
  invisible(sort(new.lens))
}

# Pairwise PERMANOVA
pairwise.adonis2 <- function(resp, fact, p.method = "none", nperm = 999) {
  require(vegan)
  resp <- as.matrix(resp)
  fact <- factor(fact)
  fun.p <- function(i, j) {
    fact2 <- droplevels(fact[as.numeric(fact) %in% c(i, j)])
    index <- which(fact %in% levels(fact2))
    resp2 <- as.dist(resp[index, index])
    result <- adonis2(resp2 ~ fact2, permutations = nperm)
    result$`Pr(>F)`[1]
  }
  multcomp <- pairwise.table(fun.p, levels(fact), p.adjust.method = p.method)
  return(list(fact = levels(fact), p.value = multcomp, p.adjust.method = p.method))
}
