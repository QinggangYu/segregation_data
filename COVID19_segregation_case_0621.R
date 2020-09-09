#load packages
library(lme4)
library(lmerTest)
library(plyr)
library(dplyr)
library(ggplot2)
library(knitr)
library(forcats)
library(scales)
library(ggcorrplot)
library(readxl)
library(gtools)
library(emmeans)
library(MuMIn)


#change this
setwd('/Users/qinggang/Desktop/BCG_data/Irene_data/')

#read in files
cases <- read.csv('time_series_covid19_confirmed_US_0621.csv')
cases_exl <- read_excel('time_series_covid19_confirmed_US_0621.xlsx')
countycitydata <- read.csv('Data_county_city_updated.csv')
lockdown <- read_excel('DateofLockdown_States_county.xlsx')
reopen <- read_excel('Dateofreopen_States_county.xlsx')

#Remove columns before Combined_Key
casesbycounty <- cases[11:ncol(cases)]

#Set baseline -- for each county, day 1 = first day of 20 confirmed infections, to a max of 30 days
cases_baseline <- data.frame(matrix(ncol = 30, nrow = nrow(casesbycounty)))
cases_baseline <- cbind(casesbycounty[,1], cases_baseline)

newrows <- data.frame(matrix(ncol = 30, nrow = 0))

for (i in 1:nrow(casesbycounty)){
  row <- casesbycounty[i,2:ncol(casesbycounty)]
  nonzero <- row[row>=20]
  #use the line below if use 100 cases as start point
  #nonzero <- row[row>=100] 
  len <- length(nonzero)
  if (len >= 30) { #this is for countries with more than 30 timepoints - we want just the first 30 points
    newrow <- nonzero[1:30]
  } else if (len >= 15) { #this picks only the counties with at least 15 non-zero days (15 days of data > 20 cases)
    newrow <- c(nonzero, rep(NA, 30-len))
  } else {
    newrow <- rep(NA, 30)
  }
  newrows <- rbind(newrows, newrow)
}

cases_baseline[,2:31] <- newrows

#Fix NAs
countycitydata[countycitydata=="#N/A"] <- NA
#Retains counties under the 100 largest metro areas (should results in 577 counties)
countycitydata <- subset(countycitydata, !is.na(countycitydata$City))

#Join dataframes
cases_baseline <- plyr::rename(cases_baseline, replace = c("casesbycounty[, 1]" = "County_State"))
cases_baseline <- left_join(countycitydata, cases_baseline, by = "County_State")

#Remove Counties with insufficient amount of data (NA for x1, or < 15 days of data)
cases_baseline <- subset(cases_baseline, !is.na(cases_baseline$X1))

######################################
### calculate lockdown/reopen data ###
######################################

#prepare data
casesbycounty <- cases_exl[11:ncol(cases_exl)]
casesbycounty <- plyr::rename(casesbycounty, replace = c("Combined_Key" = "County_State"))
matched_data <- left_join(countycitydata, casesbycounty, by = "County_State")
matched_data_s <- matched_data[, c(3:5, 29:ncol(matched_data))]

#hold data from loop
location <- data.frame(matrix(ncol = 3, nrow = 0))
time_list <- data.frame(matrix(ncol = 4, nrow = 0))
count_list <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in 1:nrow(matched_data_s)){
  state <- matched_data_s[i, 1]
  #grab lockdown date of the state the county is in
  lockdown_date <- lockdown[lockdown$State == state, "Date"]
  lockdown_date <- as.Date(lockdown_date$Date)
  reopen_date <- reopen[reopen$State == state, "Date"]
  reopen_date <- as.Date(reopen_date$Date)
  #grab the actual case data
  row <- matched_data_s[i,4:ncol(matched_data_s)]
  nonzero <- row[row>=20]
  #use the line below if use 100 cases as start point
  #nonzero <- row[row>=100]
  num_day <- length(nonzero)
  #Only care about counties with more than 15 days of data
  if (num_day >= 15){
    #First day is straighforward, just count the day backward from the last day
    first_day_ind <- length(row) - num_day + 1
    first_day <- names(row)[first_day_ind]
    first_day <- as.Date(as.numeric(first_day), origin = "1899-12-30")
    if (num_day <= 30){ #if less than 30 days, the last day would be the last day of the dataset
      last_day <- names(row)[length(row)]
      last_day <- as.Date(as.numeric(last_day), origin = "1899-12-30")
    } else { #if more than 30 days, the last day is the 30th day after the first day
      last_day_ind <- first_day_ind + 29
      last_day <- names(row)[last_day_ind]
      last_day <- as.Date(as.numeric(last_day), origin = "1899-12-30")
      num_day <- 30 #reset the number of days from > 30 to 30
    }
    diff_last_lock <- as.numeric(difftime(lockdown_date, last_day, units = c("days")))
    diff_lock_first <- as.numeric(difftime(first_day, lockdown_date, units = c("days")))
    if (diff_last_lock >= 0){ #if lockdown happens after the last day of data (condition 1)
      count <- 0  #number of days after lockdown would be zero
      coding <- -1  #code as -1
    } else if (diff_lock_first >= 0){  #if lockdown happens before the first day of data (condition 2)
      count <- num_day   #number of days after lockdown would be the same as number of days of data
      if (diff_lock_first < 14){  
        coding <- -1  # if lockdown happens less than 14 days before the first day, code as -1
      } else {
        coding <- 1
      }
    } else { # if lockdown happens between the first day and the last day (condition 3)
      count <- as.numeric(difftime(last_day, lockdown_date, units = c("days")))
      count <- count + 1   #count the difference between the lockdown date and the last day
      coding <- -1
    }
    diff_last_open <- as.numeric(difftime(reopen_date, last_day, units = c("days")))
    diff_open_first <- as.numeric(difftime(first_day, reopen_date, units = c("days")))
    if (diff_last_open >= 0){ #if reopen happens after the last day of data (condition 1)
      count2 <- 0  #number of days after lockdown would be zero
    } else if (diff_open_first >= 0){  #if reopen happens before the first day of data (condition 2)
      count2 <- num_day   #number of days after lockdown would be the same as number of days of data
    } else { # if reopen happens between the first day and the last day (condition 3)
      count2 <- as.numeric(difftime(last_day, reopen_date, units = c("days")))
      count2 <- count2 + 1   #count the difference between the lockdown date and the last day
    }
    time_var <- c(lockdown_date, reopen_date, first_day, last_day)
    location <- rbind(location, matched_data_s[i,1:3])
    time_list <- rbind(time_list, time_var)
    count_var <- c(count, count2, coding)
    count_list <- rbind(count_list, count_var)
  }
}
df <- cbind.data.frame(location, time_list, count_list)
names(df) <- c("State", "County", "City", "lockdown_date", "reopen_date", "first_day", "last_day", "count", "count2", "coding")
# Since R resets the date as numeric, switch it back to the actual day
df$first_day <- as.Date(df$first_day, origin = "1970-1-1")
df$last_day <- as.Date(df$last_day, origin = "1970-1-1")
df$lockdown_date <- as.Date(df$lockdown_date, origin = "1970-1-1")
df$reopen_date <- as.Date(df$reopen_date, origin = "1970-1-1")

#Add the lockdown data to the dataframe  
#test <- cases_baseline$County_State == df$County (double check if the county lines up between the two df)
cases_baseline$ld_count <- df$count
cases_baseline$ld_code <- df$coding
cases_baseline$rop_count <- df$count2

#Reshape to long format - 1 day per row
cases_baseline_long <- reshape(data = cases_baseline, varying = list(29:58),
                               direction = "long")
#rename
cases_baseline_long <- plyr::rename(cases_baseline_long, replace = c("X1" = "cases"))
cases_baseline_long <- plyr::rename(cases_baseline_long, replace = c("time" = "day"))

#Natural log of cases (to account for exponential growth)
cases_baseline_long$lncases <- log(cases_baseline_long$cases)

#Center day (so we can look at main effects - at the mean day)
cases_baseline_long$c_day <- as.vector(scale(cases_baseline_long$day, center = T, scale = F))

#Convert Segregation data and covariates to numeric
cases_baseline_long$Black_2005.9 <- as.numeric(as.character(cases_baseline_long$Black_2005.9))
cases_baseline_long$Black_2000 <- as.numeric(as.character(cases_baseline_long$Black_2000))
cases_baseline_long$Black_Change <- as.numeric(as.character(cases_baseline_long$Black_Change))
cases_baseline_long$Black_Share <- as.numeric(as.character(cases_baseline_long$Black_Share))

cases_baseline_long$Hispanic_2005.9 <- as.numeric(as.character(cases_baseline_long$Hispanic_2005.9))
cases_baseline_long$Hispanic_2000 <- as.numeric(as.character(cases_baseline_long$Hispanic_2000))
cases_baseline_long$Hispanic_Change <- as.numeric(as.character(cases_baseline_long$Hispanic_Change))
cases_baseline_long$Hispanic_Share <- as.numeric(as.character(cases_baseline_long$Hispanic_Share))

cases_baseline_long$Asian_2005.9 <- as.numeric(as.character(cases_baseline_long$Asian_2005.9))
cases_baseline_long$Asian_2000 <- as.numeric(as.character(cases_baseline_long$Asian_2000))
cases_baseline_long$Asian_Change <- as.numeric(as.character(cases_baseline_long$Asian_Change))
cases_baseline_long$Asian_Share <- as.numeric(as.character(cases_baseline_long$Asian_Share))

cases_baseline_long$GINI <- as.numeric(as.character(cases_baseline_long$GINI))
cases_baseline_long$Population <- as.numeric(as.character(cases_baseline_long$Population))
cases_baseline_long$PopDensity <- as.numeric(as.character(cases_baseline_long$PopDensity))
cases_baseline_long$MedianAge <- as.numeric(as.character(cases_baseline_long$MedianAge))
cases_baseline_long$Income_MedianHousehold <- as.numeric(as.character(cases_baseline_long$Income_MedianHousehold))
cases_baseline_long$Poverty_Percent <- as.numeric(as.character(cases_baseline_long$Poverty_Percent))
cases_baseline_long$Over65_Percent <- as.numeric(as.character(cases_baseline_long$Over65_Percent))
cases_baseline_long$Perc_White <- as.numeric(as.character(cases_baseline_long$Percent_white_county))
cases_baseline_long$Perc_nonwhite <- 100 - cases_baseline_long$Perc_White

#Convert Segregation data to z scores

cases_baseline_long$z_Segregation_Black <- scale(cases_baseline_long$Black_2005.9, center = T, scale = T)
cases_baseline_long$z_Segregation_Hispanic <- scale(cases_baseline_long$Hispanic_2005.9, center = T, scale = T)
cases_baseline_long$z_Segregation_Asian <- scale(cases_baseline_long$Asian_2005.9, center = T, scale = T)

cases_baseline_long$z_SegregationChange_Black <- scale(cases_baseline_long$Black_Change, center = T, scale = T)
cases_baseline_long$z_SegregationChange_Hispanic <- scale(cases_baseline_long$Hispanic_Change, center = T, scale = T)
cases_baseline_long$z_SegregationChange_Asian <- scale(cases_baseline_long$Asian_Change, center = T, scale = T)

cases_baseline_long$z_BlackShare <- scale(cases_baseline_long$Black_Share, center = T, scale = T)
cases_baseline_long$z_HispanicShare <- scale(cases_baseline_long$Hispanic_Share, center = T, scale = T)
cases_baseline_long$z_AsianShare <- scale(cases_baseline_long$Asian_Share, center = T, scale = T)


#Rescale covariates
cases_baseline_long$Income_MedianHousehold_thous <- cases_baseline_long$Income_MedianHousehold/1000

cases_baseline_long$Pop_thous <- cases_baseline_long$Population/1000
cases_baseline_long$lnPop_thous <- log(cases_baseline_long$Pop_thous)

#convert covariates to z scores
cases_baseline_long$z_Income_MedianHousehold_thous <- scale(cases_baseline_long$Income_MedianHousehold_thous,
                                                            center = T, scale = T)
cases_baseline_long$z_Pop_thous <- scale(cases_baseline_long$Pop_thous, center = T, scale = T)
cases_baseline_long$z_lnPop_thous <- scale(cases_baseline_long$lnPop_thous, center = T, scale = T)

cases_baseline_long$z_GINI <- scale(cases_baseline_long$GINI, center = T, scale = T)
cases_baseline_long$z_PopDensity <- scale(cases_baseline_long$PopDensity, center = T, scale = T)
cases_baseline_long$z_MedianAge <- scale(cases_baseline_long$MedianAge, center = T, scale = T)
cases_baseline_long$z_PercentPoverty <- scale(cases_baseline_long$Poverty_Percent, center = T, scale = T)
cases_baseline_long$z_PercentOver65 <- scale(cases_baseline_long$Over65_Percent, center = T, scale = T)
cases_baseline_long$z_Perc_nonwhite <- scale(cases_baseline_long$Perc_nonwhite, center = T, scale = T)
cases_baseline_long$z_ld_count <- scale(cases_baseline_long$ld_count, center = T, scale = T)
cases_baseline_long$z_rop_count <- scale(cases_baseline_long$rop_count, center = T, scale = T)

cases_baseline_long$lncases_percap <- log(cases_baseline_long$cases/cases_baseline_long$Pop_thous)
#cases_baseline_long$lncases <- cases_baseline_long$lncases/cases_baseline_long$Pop_thous

#Black segregation (if want to control for number of days after lockdown, add in z_ld_count + c_day:z_ld_count to the formula)
cases_mod_B <- lmer(lncases_percap ~ c_day + z_Segregation_Black +
                      c_day:z_Segregation_Black + z_Perc_nonwhite + c_day:z_Perc_nonwhite  +
                      z_GINI + c_day:z_GINI + z_lnPop_thous + c_day:z_lnPop_thous + z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous +
                      z_PopDensity + c_day:z_PopDensity + z_PercentOver65 + c_day:z_PercentOver65 + 
                      z_Segregation_Black:z_GINI + c_day:z_Segregation_Black:z_GINI +
                      (c_day||City) + (c_day||City:County_State) , data=cases_baseline_long,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(cases_mod_B)
#if not converging, grab the estimate
ss1 <- getME(cases_mod_B, c("theta"))
# and try again
c_mod_B_updated <- update(cases_mod_B, start=ss1, 
                          control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(c_mod_B_updated)

#if it still does not work, get a slightly perturbed estimate
#see details here: https://rdrr.io/cran/lme4/man/convergence.html
ss2 <- runif(length(ss1),ss1/1.1,ss1*1.1)
c_mod2B_updated <- update(cases_mod_B, start=ss2, 
                          control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))

modB_trd <- emtrends(c_mod_B_updated, pairwise ~ z_GINI*z_Segregation_Black, 
                     at = list(z_GINI = c(-1, 1), z_Segregation_Black = c(-1, 1)), var = "c_day")
summary(modB_trd)

#Hispanic segregation
cases_mod_H <- lmer(lncases ~ c_day + z_Segregation_Hispanic + 
                      c_day:z_Segregation_Hispanic + z_Perc_nonwhite + c_day:z_Perc_nonwhite + 
                      z_GINI + c_day:z_GINI + z_lnPop_thous + c_day:z_lnPop_thous + z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous +
                      z_PopDensity + c_day:z_PopDensity + z_PercentOver65 + c_day:z_PercentOver65 +
                      z_Segregation_Hispanic:z_GINI + c_day:z_Segregation_Hispanic:z_GINI +
                      (c_day||City) + (c_day||City:County_State), data=cases_baseline_long,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(cases_mod_H)

#if not converging, grab the estimate
ss1 <- getME(cases_mod_H, c("theta"))
# and try again
c_mod_H_updated <- update(cases_mod_H, start=ss1, 
                          control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(c_mod_H_updated)

ss2 <- runif(length(ss1),ss1/1.01,ss1*1.01)
c_mod2H_updated <- update(cases_mod_H, start=ss2, 
                          control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(c_mod2H_updated)

modH_trd <- emtrends(c_mod_H_updated, pairwise ~ z_GINI*z_Segregation_Hispanic, 
                     at = list(z_GINI = c(-1, 1), z_Segregation_Hispanic = c(-1, 1)), var = "c_day")
summary(modH_trd)


#Asian segregation
cases_mod_A <- lmer(lncases ~ c_day + z_Segregation_Asian + c_day:z_Segregation_Asian + 
                      z_Perc_nonwhite +c_day:z_Perc_nonwhite + 
                      z_ld_count + c_day:z_ld_count + z_rop_count + c_day:z_rop_count +
                      z_GINI + c_day:z_GINI + z_lnPop_thous + c_day:z_lnPop_thous + z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous + z_PopDensity + c_day:z_PopDensity +
                      z_PercentOver65 + c_day:z_PercentOver65 + 
                      z_Segregation_Asian:z_GINI + c_day:z_Segregation_Asian:z_GINI +
                      (c_day||City) + (c_day||City:County_State), data=cases_baseline_long,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(cases_mod_A)
ss1 <- getME(cases_mod_A, c("theta"))
# and try again
c_mod_A_updated <- update(cases_mod_A, start=ss1, 
                          control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
#if it still does not work, get a slightly perturbed estimate
#see details here: https://rdrr.io/cran/lme4/man/convergence.html
ss2 <- runif(length(ss1),ss1/1.1,ss1*1.1)
c_mod2A_updated <- update(cases_mod_A, start=ss2, 
                          control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(c_mod2A_updated)

modA_trd <- emtrends(c_mod2A_updated, pairwise ~ z_GINI*z_Segregation_Asian, 
                     at = list(z_GINI = c(-1, 1), z_Segregation_Asian = c(-1, 1)), var = "c_day")
summary(modA_trd)

#Base model
cases_Basemod <- lmer(lncases ~ c_day + z_Perc_nonwhite + c_day:z_Perc_nonwhite +
                      z_lnPop_thous + c_day:z_lnPop_thous + z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous +
                      z_PopDensity + c_day:z_PopDensity + z_PercentOver65 + c_day:z_PercentOver65 + 
                      (c_day||County_State) , data=cases_baseline_long,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))

cases_Basemod_nopop <- lmer(lncases ~ c_day + z_Perc_nonwhite + c_day:z_Perc_nonwhite + 
                      z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous +
                      z_PopDensity + c_day:z_PopDensity + z_PercentOver65 + c_day:z_PercentOver65 +
                      (c_day||City) + (c_day||City:County_State), data=cases_baseline_long,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
ss1 <- getME(cases_Basemod_nopop, c("theta"))
# and try again
c_updated <- update(cases_Basemod_nopop, start=ss1, 
                          control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))


r.squaredGLMM(c_mod_B_updated)
r.squaredGLMM(cases_Basemod)
r.squaredGLMM(c_mod_H_updated)
r.squaredGLMM(c_updated)

### SANITY CHECK (first 15 days)

cases_baseline_long15 <- subset(cases_baseline_long, cases_baseline_long$day <= 15)
cases_baseline_long15$c_day <- as.vector(scale(cases_baseline_long15$day, center = T, scale = F))

cases_mod_B <- lmer(lncases ~ c_day + z_Segregation_Black +
                      c_day:z_Segregation_Black + z_Perc_nonwhite + c_day:z_Perc_nonwhite +
                      z_GINI + c_day:z_GINI + z_lnPop_thous + c_day:z_lnPop_thous + z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous +
                      z_PopDensity + c_day:z_PopDensity + z_PercentOver65 + c_day:z_PercentOver65 + 
                      z_Segregation_Black:z_GINI + c_day:z_Segregation_Black:z_GINI +
                      (c_day||City) + (c_day||City:County_State) , data=cases_baseline_long15,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(cases_mod_B)

cases_mod_H <- lmer(lncases ~ c_day + z_Segregation_Hispanic + 
                      c_day:z_Segregation_Hispanic + z_Perc_nonwhite + c_day:z_Perc_nonwhite +
                      z_GINI + c_day:z_GINI + z_lnPop_thous + c_day:z_lnPop_thous + z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous +
                      z_PopDensity + c_day:z_PopDensity + z_PercentOver65 + c_day:z_PercentOver65 +
                      z_Segregation_Hispanic:z_GINI + c_day:z_Segregation_Hispanic:z_GINI +
                      (c_day||City) + (c_day||City:County_State), data=cases_baseline_long15,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(cases_mod_H)

cases_mod_A <- lmer(lncases ~ c_day +  z_Segregation_Asian + c_day:z_Segregation_Asian +
                      z_Perc_nonwhite + c_day:z_Perc_nonwhite +
                      z_GINI + c_day:z_GINI + z_lnPop_thous + c_day:z_lnPop_thous + z_Income_MedianHousehold_thous +
                      c_day:z_Income_MedianHousehold_thous + z_PopDensity + c_day:z_PopDensity +
                      z_PercentOver65 + c_day:z_PercentOver65 + 
                      z_Segregation_Asian:z_GINI + c_day:z_Segregation_Asian:z_GINI +
                      (c_day||City) + (c_day||City:County_State), data=cases_baseline_long15,
                    control = lmerControl(optCtrl = list(maxeval = 100000000, xtol_abs=1e-8, ftol_abs=1e-8)))
summary(cases_mod_A)