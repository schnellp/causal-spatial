library(tidyverse)
library(maps)
library(readxl)
library(spam)


STATES <- str_to_lower(c(state.name, "District of Columbia"))

##############################
### Circulatory death data ###
##############################


data.raw <- read.csv("data-mortality.csv")

states.vector <- str_split_fixed((as.character(data.raw$County)), ", ", n = 2)[, 2]
state.abb.dc <- c(state.abb, "DC")
state.name.dc <- c(str_to_lower(state.name), "district of columbia")

data.circ <- data.raw %>%
  mutate(State = state.name.dc[match(states.vector, state.abb.dc)]) %>%
  filter(!(State %in% c("alaska", "hawaii")),
         !is.na(State),
         County != "Clifton Forge city, VA") %>%
  mutate(
    Deaths = as.numeric(if_else(Deaths == "Suppressed", as.character(NA), as.character(Deaths))),
    Population = as.numeric(as.character(Population)),
    county_fips = str_pad(County.Code, 5, pad = "0")
  ) %>%
  select(-c(County.Code, Crude.Rate, Notes))

overall.crude.rate <- as.numeric(as.character(data.raw %>%
                                                filter(Notes == "Total") %>%
                                                pull(Crude.Rate))) / 100000

data.circ <- data.circ %>%
  mutate(Expected = overall.crude.rate * Population)


#################
### food data ###
#################

data.food <- as_tibble(read.csv("data-access.csv")) %>%
  select(FIPS, PCT_HHNV1MI) %>%
  mutate(FIPS = str_pad(FIPS, 5, pad = "0"))

data <- left_join(data.circ, data.food, by = c("county_fips" = "FIPS"))

##################
### other data ###
##################

dta <- read.csv("data-census.csv")

dta <- dta %>%
  select(region, State, county, county_fips,
         smokerate,
         TotPop, PctUrban, PctWhite, PctBlack, PctHisp,
         PctHighSchool, MedianHHInc, PctPoor, PctFemale, PctOccupied,
         PctMovedIn5, MedianHValue, PopPerSQM) %>%
  mutate(county_fips = str_pad(county_fips, 5, pad = "0"),
         State = as.character(State))

data <- left_join(data, dta, by = c("county_fips", "State"))

##################
### merge data ###
##################

data <- data %>%
  mutate(
    relative_risk = Deaths / Expected
  )

data <- data[rowSums(is.na(data[, -c(2, ncol(data))])) == 0, ] %>%
  select(FullName = County, FIPS = county_fips,
         Region = region, State, County = county,
         Population, Deaths, Expected, RelativeRisk = relative_risk,
         PCT_HHNV1MI,
         SmokeRate = smokerate,
         TotPop, PopPerSQM,
         PctUrban, PctWhite, PctBlack, PctHisp, PctFemale,
         PctHighSchool, PctPoor, PctOccupied, PctMovedIn5,
         MedianHHInc, MedianHValue
         )

  
###########################
### spatial information ###
###########################

adj <- read.csv("data-adjacency.csv")

N <- nrow(data)
CAR.W <- matrix(0, nrow = N, ncol = N)

for (i in 1 : N) {
  fips.nbs <- adj$fipsneighbor[adj$fipscounty ==
                                 as.numeric(data$FIPS[i])]
  CAR.W[i, as.numeric(data$FIPS) %in% fips.nbs] <- 1
}
diag(CAR.W) <- 0
CAR.D <- diag(pmax(rowSums(CAR.W), 1))

NB <- list()
for (i in 1 : nrow(CAR.W)) {
  NB[[i]] <- which(CAR.W[i, ] == 1)
}

CAR.W <- as.spam(CAR.W)
CAR.D <- as.spam(CAR.D)

write.csv(data, file = "data-combined.csv")
