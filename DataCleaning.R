setwd("C:/Users/Lara/Desktop/ADNIMERGE")

## for arrange function
#install.packages("dplyr")
## for corrplot
#install.packages("corrplot")
#install.packages("Hmisc")
#install.packages('ADNIMERGE_0.0.1.tar.gz', repos = NULL, type = 'source')

library(ADNIMERGE) 
library(dplyr)
library(magrittr)
library(corrplot)
library(tableone)
library(xtable)
library(ggplot2)
library(gtools)
library(mstate)
library(colorspace) 

data(adnimerge)
data(datadic)
data(treatdis)

countpatients<-unique(adnimerge$RID)
# Data cleaning ####

#Before moving on to the descriptive statistics and data analysis, some data cleaning will be performed. This is mostly aimed at reducing the number of rows with missing diagnoses.

## Save original dataframe 
original_adnimerge <- adnimerge
Dxconfr<-original_adnimerge[,c("DX","DX.bl","RID")]

## Remove patients without baseline measurement
adnimerge <- adnimerge[!(is.na(adnimerge$DX) & adnimerge$VISCODE == "bl"),]

countpatients1<-unique(adnimerge$RID)
## Remove patients with only 1 measurement (censoring time = baseline)
#rid is Participant roster ID
once <- table(adnimerge$RID) == 1
once <- names(once[once == T])

adnimerge <- adnimerge[!(adnimerge$RID %in% once),]

rm(once)

countpatients2<-unique(adnimerge$RID)


#### Identify deaths####
##grepl: searches for matches in characters
##tolower:convert uppercase in lowercase
death1 <- grepl('death', tolower(treatdis$WDREASON))
death2 <- grepl('death', tolower(treatdis$WDPARTCOM))
death3 <- grepl('die', tolower(treatdis$WDPARTCOM))
death4 <- grepl('expire', tolower(treatdis$WDREASON))
death5 <- grepl('expire', tolower(treatdis$WDPARTCOM))
treatdis$death <- death1 | death2 | death3 | death4 | death5

deaths <- treatdis[treatdis$death == T,]
deaths <- deaths[!(deaths$RID %in% c(222, 393, 5129, 1368, 2110, # manual check 
                                     2205, 4530, 4575, 4614, # showed these were 
                                     4652, 4888, 4947, 1260, 1646, 1694)), ] # not actual deaths

#describe(deaths)
#describe(treatdis)
#The unique() function in R is used to eliminate or delete the duplicate values or the 
#rows present in the vector, data frame, or matrix as well.
deaths_bool <- deaths[deaths$RID %in% unique(adnimerge$RID),] %>% 
  select(RID, EXAMDATE, death) %>%
  mutate(RID = as.numeric(RID),
         EXAMDATE = as.Date(EXAMDATE)) %>%
  arrange(RID, desc(EXAMDATE)) %>%
  filter(!duplicated(RID))

# add death data to all data 
baseline.dates <- adnimerge[adnimerge$VISCODE == "bl",] %>%
  mutate(RID = as.numeric(RID),
         EXAMDATE = as.Date(EXAMDATE)) %>%
  select(RID, EXAMDATE) %>%
  rename(bl.date = EXAMDATE)
deaths <- left_join(deaths_bool, baseline.dates, by = "RID") %>%
  mutate(Month.bl = as.numeric(difftime(EXAMDATE, bl.date, unit = 'days'))/(365.25/12),
         death = "Death",
         VISCODE = NA) %>%
  rename(DX = death) %>%
  select(RID, EXAMDATE, DX, VISCODE, Month.bl)

rm(baseline.dates, death1, death2, death3, death4, death5, deaths_bool)



### Convert data to long format ####
#To start with, data should be converted to long format with only the first and last date of each diagnosis available. 

# correct long format containing only diagnoses and dates 
long <- adnimerge %>%
  select(RID, EXAMDATE, DX, VISCODE, EXAMDATE.bl, Month.bl) %>%
  mutate(Month.bl = round(Month.bl, 2),
         DX = as.character(DX))

# assume: if after dementia missing state, then stays in dementia 
for(i in 1:(nrow(long)-1)){
  if(long[i, "RID"] == long[i+1, "RID"]){
    if(!is.na(long[i, "DX"]) &
       is.na(long[i+1, "DX"]) &
       long[i, "DX"] == "Dementia")
      long[i+1, "DX"] = "Dementia"
  }
}
#arrange() orders the rows of a data frame by the values of selected columns
long <- long %>%
  arrange(RID, DX, EXAMDATE)

original_long<-long
# remove in between measurements of the same diagnosis 
long$remove <- NA

for(i in 2:(nrow(long)-1)){
  if(long[i, "RID"] == long[i+1, "RID"] &
     !is.na(long[i, "DX"]) &
     !is.na(long[i+1, "DX"]) &
     long[i, "DX"] == long[i+1, "DX"]){
    
    if(long[i-1, "RID"] == long[i, "RID"]  &
       !is.na(long[i-1, "DX"]) &
       long[i-1, "DX"] == long[i, "DX"]){
      long[i, "remove"] = T
    }
    
  }
} 

long <- long[is.na(long$remove),] %>%
  arrange(RID, EXAMDATE)

# keep only the first date of diagnosis 
first_date <- long %>%
  filter(!is.na(DX)) %>%
  mutate(RID_DX = paste0(RID, "_", DX)) %>%
  arrange(RID, EXAMDATE) %>%
  filter(!duplicated(RID_DX)) %>%
  select(RID, DX, EXAMDATE) %>%
  rename(first = EXAMDATE)
first_date <- reshape(first_date, idvar = "RID", timevar = "DX", direction = "wide")

# keep only the last date of diagnosis
last_date <- long %>%
  filter(!is.na(DX)) %>%
  mutate(RID_DX = paste0(RID, "_", DX)) %>%
  arrange(RID, desc(EXAMDATE)) %>%
  filter(!duplicated(RID_DX)) %>%
  select(RID, DX, EXAMDATE) %>%
  rename(last = EXAMDATE)
last_date <- reshape(last_date, idvar = "RID", timevar = "DX", direction = "wide")

long <- left_join(long, first_date)
long <- left_join(long, last_date)
rm(first_date, last_date)

# remove patients that regress (MCI -> CN, dem -> MCI, dem -> CN) stages 
#The R function setdiff indicates which elements of a vector or data frame X are 
#not existent in a vector or data frame Y
#rid to remove are change into NA then a comparison is made between previous rid and present 
#rid and the different ones are removed
pt.before <- unique(long$RID)
long <- long %>%
  mutate(RID = ifelse((!is.na(first.MCI) & !is.na(last.CN) & first.MCI < last.CN) |
                        (!is.na(first.Dementia) & !is.na(last.CN) & first.Dementia < last.CN) |
                        (!is.na(first.Dementia) & !is.na(last.MCI) & first.Dementia < last.MCI), 
                      NA, RID))
pt.after <- unique(long$RID)[!is.na(unique(long$RID))]
length(setdiff(pt.before, pt.after))
rm(pt.before, pt.after)

long <- long[!is.na(long$RID),]
longcout <- unique(long$RID)
# remove NA values if they are in an interval between two of the same diagnoses
for(i in 1:(nrow(long))){
  if(is.na(long[i, "DX"])){
    if(!is.na(long[i, "first.CN"]) & 
       long[i, "first.CN"] <= long[i, "EXAMDATE"] &
       long[i, "EXAMDATE"] <= long[i, "last.CN"]){
      long[i, "RID"] = NA
    }else if(!is.na(long[i, "first.MCI"]) & 
             long[i, "first.MCI"] <= long[i, "EXAMDATE"] &
             long[i, "EXAMDATE"] <= long[i, "last.MCI"]){
      long[i, "RID"] = NA
    }else if(!is.na(long[i, "first.Dementia"]) & 
             long[i, "first.Dementia"] <= long[i, "EXAMDATE"] &
             long[i, "EXAMDATE"] <= long[i, "last.Dementia"]){
      long[i, "RID"] = NA
    }
  }
} 

long <- long[!is.na(long$RID), ] %>%
  arrange(RID, EXAMDATE) %>%
  select(RID, EXAMDATE, DX, VISCODE, Month.bl)

long_withna<-long

# keep only first and last NA value
for(i in 2:(nrow(long)-1)){
  if(long[i, "RID"] == long[i+1, "RID"]){
    if(is.na(long[i-1, "DX"]) & 
       is.na(long[i, "DX"]) &
       is.na(long[i+1, "DX"]))
      long[i, "RID"] = NA
  }
}
long <- long[!is.na(long$RID), ]

#alternative code:
#long2<-long_withna
#long2$drops=NA
#for(i in 2:(nrow(long2)-1)){
#  if(long2[i, "RID"] == long2[i+1, "RID"]){
#    if(long2[i-1,"RID"]==long2[i, "RID"]){
#      if(is.na(long2[i-1, "DX"]) & 
#         is.na(long2[i, "DX"]) &
#         is.na(long2[i+1, "DX"]))
#        long2[i, "drops"] = "elim"
#    }
#  }
#}
#long2 <- long2[is.na(long2$drops), ]
#long2$drops<-NULL
#COMPARISON
#long.before <- unique(long$RID)
#long2.after <- unique(long2$RID)
#length(setdiff(long.before, long2.after))
#rm(long.before, long2.after)

countpatients4<-unique(long$RID)
# add deaths 
deaths <- deaths[deaths$RID %in% long$RID,]
long <- rbind(long, deaths) %>%
  arrange(RID, EXAMDATE)


### Looking at transitions ####
#We are specifically interested in the transitions between states. In order to take a quick 
#look at these transitions from a long format dataframe, we will write a function that produces 
#a new dataframe with as output the patient ID, consecutive diagnoses and type of transitions.

# function for crude look at transitions

transitions <- function(data){
  # keep only one line per diagnosis 
  for(i in 1:(nrow(data)-1)){
    if(data[i, "RID"] == data[i+1, "RID"]){
      if((!is.na(data[i, "DX"]) &
          !is.na(data[i+1, "DX"]) &
          data[i, "DX"] == data[i+1, "DX"]) | 
         (is.na(data[i, "DX"]) &
          is.na(data[i+1, "DX"]))){
        data[i, "RID"] = NA
      }
    }
  }
  data <- data[!is.na(data$RID), ] %>% 
    arrange(RID, Month.bl) %>%
    select(RID, DX, Month.bl)
  
  # add person indexing
  data$nr <- 1
  for(i in 1:(nrow(data)-1)){
    if(data[i, "RID"] == data[i+1, "RID"]){
      data[i+1, "nr"] = data[i, "nr"] + 1
    }
  }
  
  # name transitions 
  data$transitiontype <- as.character(NA)
  for(i in 1:(nrow(data)-1)){
    if(data[i, "RID"] == data[i+1, "RID"]){
      data[i, "transitiontype"] = paste0(data[i, "DX"], "_to_", data[i+1, "DX"])
    }
  }
  
  # add information about starting and final states
  # starting state
  ## unknown baseline diagnosis
  data$transitiontype <- ifelse(data$nr == 1 & is.na(data$transitiontype) 
                                & is.na(data$DX),
                                "NoTransitions_BaselineUnknown", data$transitiontype)
  ## known baseline diagnosis
  data$transitiontype <- ifelse(data$nr == 1 & is.na(data$transitiontype) 
                                & !is.na(data$DX),
                                paste0("NoTransitions_Baseline_", data$DX), data$transitiontype)
  
  # final state
  data$transitiontype <- ifelse(data$nr != 1 & is.na(data$transitiontype) 
                                & is.na(data$DX),
                                "FinalStateUnknown", data$transitiontype)
  
  data$transitiontype <- ifelse(data$nr != 1 & is.na(data$transitiontype) 
                                & !is.na(data$DX),
                                paste0("Finalstate_", data$DX), data$transitiontype)
  
  data$nr <- NULL
  return(data)
}


## Transitions
trans <- as.data.frame(table(transitions(long)$transitiontype))


#### Adding states ####
#Since we are dealing with many NA values, we are introducing two extra states: 
#  - MCI-: a state that comes between CN and MCI
# - MCI+: a state that comes between MCI and dementia


# adding MCI- and MCI+
AddingStates <- left_join(long, transitions(long))
AddingStates$DX <- ifelse(is.na(AddingStates$DX) & AddingStates$transitiontype == "NA_to_MCI", "MCI-", AddingStates$DX)
AddingStates$DX <- ifelse(is.na(AddingStates$DX) & AddingStates$RID %in% c(520, 1190, 6001, 6327), "MCI-", AddingStates$DX)
#manual added? others: if there's Na and then MCI

for(i in 1:(nrow(AddingStates)-2)){
  if(AddingStates[i, "RID"] == AddingStates[i+1, "RID"] &
     !is.na(AddingStates[i, "transitiontype"]) &
     AddingStates[i, "transitiontype"] == "MCI_to_NA"){
    AddingStates[i+1, "DX"] = "MCI+"
    
    if(is.na(AddingStates[i+2, "DX"])) AddingStates[i+2, "DX"] = "MCI+"
  }
}

long2 <- AddingStates
long2$transitiontype <- NULL
rm(AddingStates)


## Censor patients that only have CN observations 

#Since for patients with only CN observations, we cannot reasonably assume they 
#will experience cognitive decline, we will be censoring these patients on their last known CN observation. 

CN_to_NA_pat <- transitions(long2)
CN_to_NA_pat <- CN_to_NA_pat$RID[CN_to_NA_pat$transitiontype == "CN_to_NA"]
#ridorig<-subset(original_adnimerge, select = c(RID, DX))
CN_to_NA <- long2[long2$RID %in% CN_to_NA_pat,]
CN_to_NA <- CN_to_NA[!is.na(CN_to_NA$DX),]

long3 <- long2[!(long2$RID %in% CN_to_NA_pat),]
long3 <- rbind(long3, CN_to_NA) %>%
  arrange(RID, EXAMDATE)

rm(CN_to_NA_pat, CN_to_NA)


# remove patients that only have one observation (baseline) left 
once <- table(long3$RID) == 1
once <- names(once[once == T])

long3 <- long3[!(long3$RID %in% once),]

rm(once)

countpatients5<-unique(long3$RID)
# numbers per transition
trans3 <- as.data.frame(table(transitions(long3)$transitiontype))

## Take out MCI- state

#Since patients have a probability of $1$ to leave the MCI- state to MCI, as this is how the MCI- state was defined, 
#we create a new dataframe where we remove the MCI- states. 

long4 <- long3[long3$DX != "MCI-",]
#!!!!!!long4 has all NA values in the last row
long4<-head(long4, -1)

# numbers per transition
trans4 <- as.data.frame(table(transitions(long4)$transitiontype))

rm(long, long2, trans)

lon3couter <- unique(long3$RID)