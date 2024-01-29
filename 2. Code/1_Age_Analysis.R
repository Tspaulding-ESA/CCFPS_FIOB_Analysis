###################### Estimate Age based on length #######################
source(file.path("2. Code","0_Setup.R"))
## Retrieve Length at capture from release data ===========================
release_data <- read.csv(file.path("1. Data","Inputs","ReleaseData.csv"))

# Develop Age-length key based on a methods from 
# fishR Vignette (D.Ogle 2013), based on methods from 
# Iserman and Knight (2015)

## Striped Bass Aging =====================================================
# Building the age-length key table ####
set.seed(42069)

AgeData <- read.csv(file.path("1. Data","Inputs","AgedFish.csv"), header = T, stringsAsFactors = F)

stb.age <- AgeData[AgeData$Species == "Striped Bass" & AgeData$Length_cm >20 & AgeData$Age < 8,] %>%
  select(Length_cm, Age) %>%
  arrange(Age)

x <- stb.age$Age
y <- stb.age$Length_cm

alk.lm <- lm(y ~ I(log(x)))

predicted.intervals <- predict(alk.lm,data.frame(x=stb.age$Age),interval='confidence', level=0.995)

plot(stb.age$Age, stb.age$Length_cm)
lines(stb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
lines(stb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1)
lines(stb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1)

x_new <- rep(seq(2,7, by = 1), 100000)

alk.lm$fitted.values <- predict(alk.lm, data.frame(x = x_new))

y_new <- simulate(alk.lm)[,1]

#plot(x_new, y_new, col = 'red')
#points(stb.age$Age, stb.age$Length_cm, pch = 16, )
#lines(stb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
#lines(stb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1, lty = 2)
#lines(stb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1, lty = 2)

alk_new <- tibble(Age = x_new, Length_cm = y_new)

# Find the minimum size for use in lencat function
FSA::Summarize(~Length_cm, data = alk_new, digits = 1)

# use lencat function to define length categories for each fish starting at min 24, and increasing by 2 
stb.age1 <- FSA::lencat(~Length_cm, data = alk_new, startcat=0, w=1)

stb.raw <- with(stb.age1, table(LCat,Age))
stb.key <- prop.table(stb.raw, margin = 1)
round(stb.key,2)
saveRDS(stb.key, file.path("1. Data","Outputs","SB_ALKey.rds"))

## Largemouth Bass Aging ==================================================
lmb.age <- AgeData[AgeData$Species == "Largemouth Bass",] %>%
  select(Length_cm, Age) %>%
  arrange(Age)

x <- lmb.age$Age
y <- lmb.age$Length_cm

alk.lm <- lm(y ~ I(log(x)))

predicted.intervals <- predict(alk.lm,data.frame(x=lmb.age$Age),interval='confidence', level=0.995)

#plot(lmb.age$Age, lmb.age$Length_cm)
#lines(lmb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
#lines(lmb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1)
#lines(lmb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1)

x_new <- rep(seq(2,4, by = 1), 100000)

alk.lm$fitted.values <- predict(alk.lm, data.frame(x = x_new))

y_new <- simulate(alk.lm)[,1]

#plot(x_new, y_new, col = 'red')
#points(lmb.age$Age, lmb.age$Length_cm, pch = 16, )
#lines(lmb.age$Age, predicted.intervals[,1], col = 'red', lwd = 3)
#lines(lmb.age$Age, predicted.intervals[,2], col = 'black', lwd = 1, lty = 2)
#lines(lmb.age$Age, predicted.intervals[,3], col = 'black', lwd = 1, lty = 2)


alk_new <- tibble(Age = x_new, Length_cm = y_new)

# Find the minimum size for use in lencat function
FSA::Summarize(~Length_cm, data = alk_new, digits = 1)

# use lencat function to define length categories for each fish starting at min 24, and increasing by 2 
lmb.age1 <- FSA::lencat(~Length_cm,data=alk_new,startcat=18,w=1)

lmb.raw <- with(lmb.age1, table(LCat,Age))
lmb.key <- prop.table(lmb.raw, margin = 1)
round(lmb.key,2)

saveRDS(lmb.key, file.path("1. Data","Outputs","LMB_ALKey.rds"))

# Cleanup #################################################################
rm(alk_new, alk.lm, AgeData, lmb.age, lmb.age1,
   predicted.intervals, release_data, stb.age, stb.age1,
   lmb.raw, lmb.key, stb.raw, stb.key, x, x_new, y, y_new)
