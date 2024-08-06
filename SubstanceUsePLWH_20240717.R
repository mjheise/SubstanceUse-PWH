#################################################
#                                               #  
#        CNICS: Substance Use Analyses          #
#                                               #
#                  MJ Heise                     #
#                June 14, 2024                  #
#                                               #
#################################################

# CODE DESCRIPTION: This code analyzes data from the Centers for AIDS Research 
# Network of Integrated Clinical Systems (CNICS) cohort of people living with HIV. 
# Data are from 2 years before and after COVID-19 onset (March 2020). The focal 
# analysis examines predictors of substance use. 
#
# 1. DATA CLEANING: 
# Organize and create variables for analysis.
#
# 2. DESCRIPTIVE STATISTICS:
# Summary of descriptive statistics for sample (N = 7126).
#
# 3. MISSING DATA:
# Examine whether variables predicted missingness in the post-shelter-in-place 
# period and visualize missingness across the sample.
#
# 4. PREDICTORS OF MODERATE/HIGH SUD POST-SIP:
# Examine predictors of whether a participant reported moderate/high SUD risk at
# post-SIP clinic visits in a generalized mixed effects model.
#
# 5. PREDICTORS OF WHO INCREASED IN SUD RISK:
# Examine predictors of who increased from none/mild to moderate/severe SUD risk
# in a generalized linear model. 
#
# 6. CHANGES IN SPECIFIC SUBSTANCES:
# Examine changes in specific substances pre/post shelter-in-place.

# Libraries
library(psych) # v.2.3.9, summary statistics
library(plotrix) # v.3.8-2, std.error function
library(tidyverse) #v.2.0.0, data management
library(lme4) # v.1.1-34, linear mixed effects models
library(lmerTest) # v.3.1-3, p-values from linear mixed effects models
library(emmeans) # v.1.8.9, moderation marginal means
library(haven) # v.2.5.3, read stata data
library(officer) # V.0.3.15, powerpoint ggplot
library(rvg) # v.0.2.5, powerpoint ggplot
library(finalfit) # v.1.0.7, missing data
library(ggpubr) # v.0.6.0, p-values to ggplot
library(naniar) # v.1.0.0, missing data
library(gtsummary) # v. 1.7.2, formatted regression tables

# Functions
# Output odds ratios and 95% confidence intervals from model output
# -Input: Fit object (e.g., fit1) from a regression model.
# -Notes: Additional estimates will appear at the end of model output. Low 95% CI
# corresponds to the lower 95% confidence interval, Up 95% CI corresponds to the 
# upper 95% confidence interval. 
summary_or <- function(m) {
  s <- summary(m)
  
  lowCI <- (exp(s$coefficients[,1] - 1.96*s$coefficients[,2]))
  highCI <- (exp(s$coefficients[,1] + 1.96*s$coefficients[,2]))
  or <- exp(s$coefficients[,1])
  
  s$coefficients <- cbind(s$coefficients, or)
  s$coefficients <- cbind(s$coefficients, lowCI)
  s$coefficients <- cbind(s$coefficients, highCI)
  
  colnames(s$coefficients)[5] <- "Odds Ratio"
  colnames(s$coefficients)[6] <- "Low 95% CI"
  colnames(s$coefficients)[7] <- "Up 95% CI"
  return(s)
}

# Save ggplot object in powerpoint slide
# -Input: The ggplot object that you want to save in a powerpoint.
# -Optional inputs: Specified width and height of the outputted graph.
#  If no arguments are specified, the graph will encompass the entire
#  powerpoint slide.
# -Notes: After running the function, a window will open that 
#  allows you to select the powerpoint file. The graph will
#  save on a new slide at the end of the powerpoint.
create_pptx <- function(plt = last_plot(), path = file.choose(), width = 0, height = 0){
  if(!file.exists(path)) {
    out <- read_pptx()
  } else {
    out <- read_pptx(path)
  }
  
  if (width != 0 & height != 0) {
    out %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = dml(ggobj = plt), location = ph_location(left = 0, top = 0,
                                                               width = width, height = height)) %>%
      print(target = path)
  } else {
    out %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = dml(ggobj = plt), location = ph_location_fullsize()) %>%
      print(target = path)
    
  }
  
}


#### 1. DATA CLEANING ####
# Read in CNICS PRO data
cnics <- read.csv('CNICS_SubstanceUseDat_20231024.csv',
                na = 'NA')

# Rename subject ID variable
cnics %>% rename(subNo = study_id) -> cnics

# Read in demographic data
stata <- read_dta('VLposthousing2noodup.dta')

# Select variables to merge into final dataset
stata %>%
  rename(subNo = studyid) %>%
  select(subNo, transgender, hispanic) -> stata

# Merge data
dat <- merge(cnics, stata, by = 'subNo', all.x = T)

# Remove raw CNICS PRO and demographic data
rm(cnics, stata)

# Convert date into months since lockdown (March 2018) and tidy variables
dat %>%
  mutate(weeks = as.numeric(round(difftime(dat$surveydate,"2018-03-19", units = "weeks"), 2) %>%
                          str_replace(" weeks", "")),
         months = weeks*0.23,
         monthsRound = round(months, 0),
         race = case_when(grepl(x = Race, 'Black') ~ 'Black',
                                      grepl(x = Race, 'White') ~ 'White',
                                      grepl(x = Race, 'American Indian') ~ 'Other',
                                      grepl(x = Race, 'Other') ~ 'Other',
                                      grepl(x = Race, 'Asian/Pacific Islander') ~ 'Other',
                                      grepl(x = Race, 'Multiracial') ~ 'Other',
                                      grepl(x = Race, 'Asian') ~ 'Other',
                                      grepl(x = Race, 'Pacific Islander') ~ 'Other'),
         gender = case_when(transgender == 'FTM Transgender' ~ 'Transgender man',
                            transgender == 'MTF Transgender' ~ 'Transgender woman',
                            transgender != 'FTM Transgender' & transgender != 'MTF Transgender' & BirthSex == 'Male' ~ 'Cis man',
                            transgender != 'FTM Transgender' & transgender != 'MTF Transgender' & BirthSex == 'Female' ~ 'Cis woman',
                            .default = NA),
         gender_categorical = case_when(gender == 'Transgender man' ~ NA,
                                        .default = gender),         
         sudTreatment = case_when(DependenceTreat_derivated == 'Yes' ~ 1,
                                  DependenceTreat_derivated == 'No' ~ 0, 
                                  .default = NA),
         SSI_binary = case_when(AnySevMod == 'Any Severe | Moderate Risk' ~ 1,
                                AnySevMod == 'Non-Severe & Non-Moderate Risk' ~ 0,
                                .default = NA),
         PHQ9_categorical = case_when(PHQ9_Score < 10 ~ 'Mild',
                                      dat$PHQ9_Score > 9 ~'Mod/Severe', 
                                      .default = NA),
         geographicRegion = case_when(site == 'JH' ~ 'Northeast',
                                      site == 'UNC' ~ 'South',
                                      site == 'UW' ~ 'West',
                                      site == 'UAB' ~ 'South',
                                      site == 'UCSD' ~ 'West',
                                      site == 'FENWAY' ~ 'Northeast',
                                      site == 'CWRU' ~ 'South',
                                      site == 'UCSF' ~ 'West')) -> dat


# Create variables for pre-COVID-19 substance use and depression
dat %>%
  filter(months<24 & !is.na(SSI_binary)) %>%  
  group_by(subNo) %>%
  mutate(SSI_preSIP_numObservations = sum(!is.na(SSI_binary)),
         preSIP_numVisits = sum(!is.na(surveydate)),
         PHQ9_preSIPmean = mean(PHQ9_Score, na.rm = T),
         preCovidPHQ9_binaryAvg = case_when(PHQ9_preSIPmean < 10 ~ 0, 
                                            PHQ9_preSIPmean > 0 ~ 1,
                                            .default = NA), 
         blq_binary = case_when(blq == 'UNdetectable' ~ 0,
                                blq == 'DEtectable' ~ 1,
                                .default = NA),
         preCovidVSsum = sum(blq_binary, na.rm = T),
         viralSupp_preSIP_everUnsuppressed = case_when(preCovidVSsum > 0 ~ 'not suppressed', 
                                                       preCovidVSsum == 0 ~ 'suppressed', 
                                                       .default = NA),
         preCovidPHQnum = case_when(PHQ9_categorical == 'Mild' ~ 0,
                                    PHQ9_categorical == 'Mod/Severe' ~ 1, 
                                    .default = NA),
         preCovidPHQsum = sum(preCovidPHQnum),
         PHQ9_preSIP_everSevere = case_when(preCovidPHQsum == 0 ~ 0,
                                            preCovidPHQsum > 0 ~ 1, 
                                            .default = NA),
         preCovidSSIsum = sum(SSI_binary),
         SSI_preSIP_everSevere = case_when(preCovidSSIsum == 0 ~ 0,
                                           preCovidSSIsum > 0 ~ 1,
                                           .default = NA), 
         SUDtreatmentEver_preSIPSum = sum(sudTreatment, na.rm = T), 
         SUDtreatmentEver_preSIP = case_when(SUDtreatmentEver_preSIPSum == 0 ~ 0,
                                             SUDtreatmentEver_preSIPSum > 0 ~ 1, 
                                             .default = NA)) %>%
  select(subNo, preSIP_numVisits, SUDtreatmentEver_preSIP, SSI_preSIP_everSevere,
         PHQ9_preSIP_everSevere, viralSupp_preSIP_everUnsuppressed, SSI_preSIP_numObservations) %>%
  slice(1) -> preC

dat <- merge(dat, preC, by = 'subNo', all = T)

# Remove pre-SIP dataset
rm(preC)

# Create variables for post-COVID-19 substance use and depression
dat %>%
  filter(months>24 & !is.na(SSI_binary)) %>%  
  group_by(subNo) %>%
  mutate(postCovidUse = sum(SSI_binary),
         SSI_postSIP_numObservations = sum(!is.na(SSI_binary)),
         postSIP_numVisits = sum(!is.na(surveydate)),
         SSI_postSIP_everSevere = case_when(postCovidUse > 0 ~ 1,
                                            postCovidUse == 0 ~ 0, 
                                            .default = NA),
         postCovidRS = case_when(blq == 'UNdetectable' ~ 0,
                                 blq == 'DEtectable' ~ 1, 
                                 .default = NA),
         postCovidVSsum = sum(postCovidRS, na.rm = T),
         viralSupp_postSIP_everUnsuppressed = case_when(postCovidVSsum > 0 ~ 'not suppressed', 
                                                        postCovidVSsum == 0 ~ 'suppressed', 
                                                        .default = NA),
         postCovidPHQnum = case_when(PHQ9_categorical == 'Mild' ~ 0, 
                                     PHQ9_categorical == 'Mod/Severe' ~ 1, 
                                     .default = NA),
         postCovidPHQsum = sum(postCovidPHQnum),
         PHQ9_postSIP_everSevere = case_when(postCovidPHQsum == 0 ~ 0,
                                             postCovidPHQsum > 0 ~ 1, 
                                             .default = NA),
         SUDtreatmentEver_postSIPSum = sum(sudTreatment, na.rm = T), 
         SUDtreatmentEver_postSIP = case_when(SUDtreatmentEver_postSIPSum == 0 ~ 0,
                                              SUDtreatmentEver_postSIPSum > 0 ~ 1, 
                                              .default = NA),
         postSIP_data = case_when(!is.na(SSI_binary) ~ 1,
                                  .default = NA)) %>%
  select(subNo, postSIP_numVisits, SSI_postSIP_everSevere, SUDtreatmentEver_postSIP,
         postSIP_data, PHQ9_postSIP_everSevere, viralSupp_postSIP_everUnsuppressed, SSI_postSIP_numObservations) %>%
  slice(1) -> postC

dat <- merge(dat, postC, by = 'subNo', all = T)


# Create variables based on pre- & post-SIP data for whether participants increased in substance use or not,
# and number of pre- & post-SIP visits
dat %>%
  mutate(SSI_postSIPincrease = case_when(SSI_preSIP_everSevere == 0 & SSI_postSIP_everSevere == 0 ~ 0,
                                        SSI_preSIP_everSevere == 1 & SSI_postSIP_everSevere == 0 ~ 0,
                                        SSI_preSIP_everSevere == 0 & SSI_postSIP_everSevere == 1 ~ 1,
                                        SSI_preSIP_everSevere == 1 & SSI_postSIP_everSevere == 1 ~ 0, 
                                        .default = NA),
         SSI_preSIP_numObservations = case_when(is.na(SSI_preSIP_numObservations) ~ 0, 
                                                .default = SSI_preSIP_numObservations),
         SSI_postSIP_numObservations = case_when(is.na(SSI_postSIP_numObservations) ~ 0, 
                                                 .default = SSI_postSIP_numObservations),
         preSIP_numVisits = case_when(is.na(preSIP_numVisits) ~ 0, 
                                      .default = preSIP_numVisits),
         postSIP_numVisits = case_when(is.na(postSIP_numVisits) ~ 0, 
                                       .default = postSIP_numVisits)) -> dat

# Create a dataframe for examining participant-level data (datUnique)
dat %>%
  group_by(subNo) %>%
  slice(1) -> datUnique

# Create a dataframe of only post-SIP data
datPostC <- dat %>%
  filter(months>24)


#### 2. DESCRIPTIVE STATISTICS ####
# Set total sample size for datUnique
total = length(datUnique$subNo)

# Age
describe(datUnique$Age)

# Ethnicity
datUnique %>%
  group_by(hispanic) %>%
  summarise(n = n(),
            percent = n()/total*100)

# Race
datUnique %>%
  group_by(race) %>%
  summarise(n = n(),
            percent = n()/total*100)

# Gender
datUnique %>%
  group_by(gender) %>%
  summarise(n = n(),
            percent = n()/total*100)

# Geographic region
datUnique %>%
  group_by(geographicRegion) %>%
  summarise(n = n(),
            percent = n()/total*100)

# Substance Use increase
datUnique %>%
  group_by(SSI_postSIPincrease) %>%
  summarise(n = n(),
            percent = n()/total*100)


#### 3. MISSING DATA ####
fit1 <- lm(postSIP_data ~ PHQ9_preSIP_everSevere + Age + geographicRegion + gender_categorical + 
             race + viralSupp_preSIP_everUnsuppressed, data = datUnique)
summary(fit1)

# Number of substance use PROs
describe(dat$SSI_preSIP_numObservations)
describe(dat$SSI_postSIP_numObservations)

# Number of pre- and post-SIP clinic visits
describe(dat$preSIP_numVisits)
describe(dat$postSIP_numVisits)

# Enrollment cascade
# Participant with a pre-SIP and post-SIP visit
dat %>%
  filter(preSIP_numVisits > 0 & postSIP_numVisits > 0) %>%
  group_by(subNo) %>%
  slice(1) -> temp

# Participants with only a post-SIP visit
dat %>%
  filter(postSIP_numVisits > 0 & preSIP_numVisits == 0) %>%
  group_by(subNo) %>%
  slice(1) -> temp

# Participants with only a pre-SIP visit
dat %>%
  filter(preSIP_numVisits > 0 & postSIP_numVisits == 0) %>%
  group_by(subNo) %>%
  slice(1) -> temp

# Visualize missing data across variables
vis_miss(dat[c('months', 'Age', 'gender', 'race', 'geographicRegion', 
               'situation', 'blq', 'SSI_binary', 'PHQ9_Score')])


#### 4. PREDICTORS OF MODERATE/HIGH SUD POST-SIP ####
# Rescale variables for model convergence
datPostC %>%
  mutate(monthsPostCovid_scale = scale(months-24),
         monthsPostCovid = (months-24),
         age_10yrs = (Age/10),
         gender_refM = relevel(factor(gender_categorical), ref = 'Cis man'),
         blq_refUn = relevel(factor(blq), ref = 'UNdetectable'),
         race_refW = relevel(factor(race), ref = 'White'),
         race_refO = relevel(factor(race), ref = 'Other'),
         PHQ9_categorical_refSev = relevel(factor(PHQ9_categorical), ref = 'Mod/Severe'),
         gender_refTW = relevel(factor(gender_categorical), ref = 'Transgender woman'),
         geographicRegion_refW = relevel(factor(geographicRegion), ref = 'West')) -> datPostC

# Generalized mixed effects model, moderate/high SUD risk predicted by an interaction
# between time and depression
fit1 <- glmer(SSI_binary ~ monthsPostCovid_scale*PHQ9_categorical_refSev + age_10yrs + gender_refM + race_refW + 
                geographicRegion_refW + blq_refUn +
                (1|subNo),
              family = binomial,
              data = datPostC,
              glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5)))
summary_or(fit1)

# Format table with odds ratios
fit1 %>%
  tbl_regression(exp = T) 

# Plot interaction from fit1
# Set levels of depression to facet wrap
PHQ_modSev <- 'Mod/Severe'
PHQ_mild <- 'Mild'

# Save estimates from the model for scaled months
mylist <- list(monthsPostCovid_scale=seq(-2,2,by=.35), PHQ9_categorical=c(PHQ_mild,PHQ_modSev))
plotInteract <- emmip(fit1, PHQ9_categorical_refSev~monthsPostCovid_scale,at=mylist, CIs=TRUE, pbkrtest.limit = 4691) + 
  theme_minimal() + 
  ylab('SSI (severe/moderate)')

# Unscale months for plot
plotInteract[["data"]][["monthsPostCovid_scale"]] = c(0.06, 0.06, 
                                                      2.6, 2.6,
                                                      5.1, 5.1,
                                                      7.67, 7.67, 
                                                      10.2, 10.2, 
                                                      12.74, 12.74, 
                                                      15.28, 15.28, 
                                                      17.82, 17.82, 
                                                      20.35, 20.35, 
                                                      22.89, 22.89, 
                                                      25.42, 25.42,
                                                      27.97, 27.97)

plotDat <- plotInteract[['data']]

# Convert coefficients to probability for plot
plotDat$p = exp(plotDat$yvar)/(1+exp(plotDat$yvar))
plotDat$p_up95 = exp(plotDat$yvar + plotDat$SE)/(1+exp(plotDat$yvar + plotDat$SE))
plotDat$p_low95 = exp(plotDat$yvar - plotDat$SE)/(1+exp(plotDat$yvar - plotDat$SE))

# Create dotplot of model estimates of the probability of SUD risk by PHQ9 severity
plot1 <- ggplot(plotDat, aes(x = monthsPostCovid_scale, y = p, color = PHQ9_categorical_refSev)) + 
  geom_errorbar(aes(ymin=p_low95, ymax=p_up95), width=.4, linewidth = 1) + 
  geom_point(size = 3) +
  scale_color_manual(values=c("#7fcdbb", "#2c7fb8")) + 
  xlim(-1, 25) + ylim(.25, .75) + 
  ylab('Probability of moderate/high SUD risk') + theme_minimal()

# Save to powerpoint
create_pptx(plot1)


# Generalized mixed effects model, moderate/high SUD risk predicted by a 3-way interaction
# between time, depression (PHQ-9) and gender
fit2 <- glmer(SSI_binary ~ monthsPostCovid_scale*PHQ9_categorical*gender_refTW + 
                age_10yrs + geographicRegion + blq + race_refW +
                (1|subNo),
              family = binomial(link = "logit"),
              data = datPostC,
              glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5)))
summary(fit2)

# Format table with odds ratios
fit2 %>%
  tbl_regression(exp = T) 

# Fit an additive model (rather than an interaction) for model comparison with 
# and without the 3-way interaction
fit2b <- glmer(SSI_binary ~ monthsPostCovid_scale + PHQ9_categorical + gender_refTW + 
                age_10yrs + geographicRegion + blq + race +
                (1|subNo),
              family = binomial(link = "logit"),
              data = datPostC,
              glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5)))

anova(fit2, fit2b)


# Plot interaction from fit2
# Set levels of depression to facet wrap
PHQ_modSev <- 'Mod/Severe'
PHQ_mild <- 'Mild'

# Save estimates from the model for scaled months
mylist <- list(monthsPostCovid_scale=seq(-2,2,by=.35), PHQ9_categorical=c(PHQ_mild,PHQ_modSev))
plotInteract <- emmip(fit2, PHQ9_categorical~monthsPostCovid_scale*gender_refM,at=mylist, CIs=TRUE, pbkrtest.limit = 4691) + 
  theme_minimal() + 
  ylab('SSI (severe/moderate)')

# Unscale months for plot
plotInteract[["data"]][["monthsPostCovid_scale"]] = c(0.06, 0.06,
                                                      2.6, 2.6,
                                                      5.13, 5.13,
                                                      7.67, 7.67,
                                                      10.21, 10.21,
                                                      12.74, 12.74,
                                                      15.28, 15.28,
                                                      17.82, 17.82,
                                                      20.35, 20.35,
                                                      22.89, 22.89,
                                                      25.42, 25.42,
                                                      27.96, 27.96,
                                                      0.06, 0.06,
                                                      2.6, 2.6,
                                                      5.13, 5.13,
                                                      7.67, 7.67,
                                                      10.21, 10.21,
                                                      12.74, 12.74,
                                                      15.28, 15.28,
                                                      17.82, 17.82,
                                                      20.35, 20.35,
                                                      22.89, 22.89,
                                                      25.42, 25.42,
                                                      27.96, 27.96,
                                                      0.06, 0.06,
                                                      2.6, 2.6,
                                                      5.13, 5.13,
                                                      7.67, 7.67,
                                                      10.21, 10.21,
                                                      12.74, 12.74,
                                                      15.28, 15.28,
                                                      17.82, 17.82,
                                                      20.35, 20.35,
                                                      22.89, 22.89,
                                                      25.42, 25.42,
                                                      27.96, 27.96)

plotDat <- plotInteract[['data']]

plotDat$p = exp(plotDat$yvar)/(1+exp(plotDat$yvar))
plotDat$p_up95 = exp(plotDat$yvar + plotDat$SE)/(1+exp(plotDat$yvar + plotDat$SE))
plotDat$p_low95 = exp(plotDat$yvar - plotDat$SE)/(1+exp(plotDat$yvar - plotDat$SE))

# Plot generalized mixed effects model, interaction with gender and depression over time
plot2 <- ggplot(plotDat, aes(x = monthsPostCovid_scale, y = p, color = PHQ9_categorical)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin=p_low95, ymax=p_up95), width=.4, linewidth = 1) + 
  scale_color_manual(values=c("#7fcdbb", "#2c7fb8")) + 
  xlim(0, 25) + 
  ylim(0, 1) +
  facet_wrap(vars(gender_refM)) + ylab('Probability of moderate/high SUD risk') + theme_minimal()

create_pptx(plot2)


# Generalized mixed effects model, moderate/high SUD risk predicted by a 3-way
# interaction between time, depression, and race
fit2b <- glmer(SSI_binary ~ monthsPostCovid_scale*PHQ9_categorical*race_refO + 
                age_10yrs + geographicRegion + blq + gender_categorical + 
                (1|subNo),
              family = binomial(link = "logit"),
              data = datPostC,
              glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e5)))
summary(fit2b)

# Format table with odds ratios
fit2b %>%
  tbl_regression(exp = T) 


#### 5. PREDICTORS OF WHO INCREASED IN SUD RISK ####
datUnique %>%
  mutate(age_10yrs = (Age/10),
         race_factor = factor(race),
         race_refW = relevel(race_factor, ref = 'White'),
         geographicRegion_factor = factor(geographicRegion),
         geographicRegion_refW = relevel(geographicRegion_factor, ref = 'West'),
         geographicRegion_refN = relevel(geographicRegion_factor, ref = 'Northeast')) -> datUnique

# Number of participants who increased in SUD risk (from none/mild to moderate/severe SUD risk)
datUnique %>%
  filter(preSIP_numVisits > 0 & postSIP_numVisits > 0) %>%
  group_by(SSI_postSIPincrease) %>%
  summarise(n = n())

# Generalized linear model in which an increase in SUD risk was predicted 
# by demographic and structural variables
fit3 <- glm(SSI_postSIPincrease ~ age_10yrs + gender + race_refW + geographicRegion_refW +
              PHQ9_preSIP_everSevere + viralSupp_preSIP_everUnsuppressed, data = datUnique)

summary(fit3)

# Format table with odds ratios
fit3 %>%
  tbl_regression(exp = T) 

# Examine whether there is an effect of geographic region comparing Northeast and South
fit3a <- glm(SSI_postSIPincrease ~ age_10yrs + gender_refM + race_refW + geographicRegion_refN +
              PHQ9_preSIP_everSevere + viralSupp_preSIP_everUnsuppressed, data = datUnique)

summary(fit3a)


# Plot mean SUD by geographic region
datUnique %>%
  group_by(geographicRegion) %>%
  summarise(mean = mean(SSI_postSIPincrease, na.rm = T),
            se = std.error(SSI_postSIPincrease, na.rm = T)) -> plotDat

# Plot model-extracted mean SUD by geographic region
margmeans <- as.data.frame(emmeans(fit3, pairwise~geographicRegion_refW, mode = 'satterthwaite', adjust = 'none')[["emmeans"]])

margmeans %>%
  mutate(OR = exp(emmean),
         Low95Conf = exp(emmean - 1.96*SE),
         Up95Conf = exp(emmean + 1.96*SE)) -> margmeans

ggplot(plotDat, aes(x=geographicRegion, y=mean)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  ylim(0, .15) + 
  theme_minimal()


#### 6. CHANGES IN SUBSTANCES PRE-/POST-SIP ####
# Create a new dataframe to examine increases/changes in substance use pre/post-COVID
# by subsetting participants who completed the substance use measure (SSI)
# (datIncrease)
dat %>%
  subset(!is.na(SSI_binary)) -> datIncrease

# Participants in datIncrease all completed the substance use screener (SSI)
# so replace NA (did not see question because did not report any substance use)
# with 0
datIncrease %>% 
  mutate_at(c('ast_pot2_Bn', 'ast_coc2_Bn', 'ast_pst2_Bn','ast_meth2_Bn', 'ast_inh2_Bn',
              'ast_sed2_Bn', 'ast_hal2_Bn', 'ast_opi2_Bn', 'ast_opip2_Bn', 'ast_her2_Bn',
              'ast_fen2_Bn', 'ast_narc2_Bn'), ~replace_na(.,0)) -> datIncrease

# Create composites for stimulants and opiates
datIncrease %>%
  mutate(ast_stimulants_mean = rowMeans(subset(datIncrease, select = c('ast_coc2_Bn', 'ast_pst2_Bn', 'ast_meth2_Bn')), na.rm=T),
         ast_opi_mean = rowMeans(subset(datIncrease, select = c('ast_opi2_Bn', 'ast_opip2_Bn', 'ast_her2_Bn', 'ast_fen2_Bn', 'ast_narc2_Bn')), na.rm=T)) -> datIncrease

# Create summary variables for the proportion of visits in which the participant
# reported using the substance
datIncrease %>%
  group_by(subNo, Lockdown) %>%
  summarize(sub_cannM = mean(ast_pot2_Bn, na.rm = T),
            sub_cocM = mean(ast_coc2_Bn, na.rm = T),
            sub_pstM = mean(ast_pst2_Bn, na.rm = T),
            sub_methM = mean(ast_meth2_Bn, na.rm = T),
            sub_opisM = mean(ast_opi2_Bn, na.rm = T),
            sub_opiPxM = mean(ast_opip2_Bn, na.rm = T),
            sub_herM = mean(ast_her2_Bn, na.rm = T),
            sub_fenM = mean(ast_fen2_Bn, na.rm = T),
            sub_narcM = mean(ast_narc2_Bn, na.rm = T),
            sub_inhM = mean(ast_inh2_Bn, na.rm = T),
            sub_sedM = mean(ast_sed2_Bn, na.rm = T)) %>%
  group_by(Lockdown) %>%
  summarize(cann_M = mean(sub_cannM, na.rm = T),
            coc_M = mean(sub_cocM, na.rm = T),
            pst_M = mean(sub_pstM, na.rm = T),
            meth_M = mean(sub_methM, na.rm = T),
            opis_M = mean(sub_opisM, na.rm = T),
            opiPx_M = mean(sub_opiPxM, na.rm = T),
            her_M = mean(sub_herM, na.rm = T),
            fen_M = mean(sub_fenM, na.rm = T),
            narc_M = mean(sub_narcM, na.rm = T),
            inh_M = mean(sub_inhM, na.rm = T),
            sed_M = mean(sub_sedM, na.rm = T),
            cann_SE = std.error(sub_cannM, na.rm = T),
            coc_SE = std.error(sub_cocM, na.rm = T),
            pst_SE = std.error(sub_pstM, na.rm = T),
            meth_SE = std.error(sub_methM, na.rm = T),
            opis_SE = std.error(sub_opisM, na.rm = T),
            opiPx_SE = std.error(sub_opiPxM, na.rm = T),
            her_SE = std.error(sub_herM, na.rm = T),
            fen_SE = std.error(sub_fenM, na.rm = T),
            narc_SE = std.error(sub_narcM, na.rm = T),
            inh_SE = std.error(sub_inhM, na.rm = T),
            sed_SE = std.error(sub_sedM, na.rm = T)) %>%
  pivot_longer(cols = 'cann_M':'sed_SE',
               names_to = c('substance','test'), 
               names_sep = '_',
               values_to = c('value')) %>%
  pivot_wider(id_cols = c('Lockdown', 'substance'), 
              names_from = 'test',
              values_from = 'value') %>%
  filter(substance != 'opis') -> plotTable 

# Set factor levels for plotting
plotTable %>%
  mutate(substance = factor(plotTable$substance, levels = c('cann', 'coc', 'pst', 'meth', 
                                                                     'fen', 'her', 'narc', 'opiPx', 
                                                                     'sed', 'inh')),
         Lockdown = factor(plotTable$Lockdown, levels = c('Before', 'After'))) -> plotTable

# Plot
plot3 <- ggplot(data=plotTable, aes(x=substance, y=M, fill = Lockdown)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=M-SE, ymax=M+SE), width=.2,
                position=position_dodge(.9)) + 
  labs(y = 'Substance use pre/post COVID onset') +
  scale_fill_manual(values=c("#fec44f", "#f03b20")) + 
  theme_minimal()

create_pptx(plot3)

# Pivot data to wide for paired t-tests for each substance
datIncrease %>%
  group_by(subNo, Lockdown) %>%
  summarize(sub_cannM = mean(ast_pot2_Bn, na.rm = T),
            sub_stimM = mean(ast_stimulants_mean, na.rm = T),
            sub_cocM = mean(ast_coc2_Bn, na.rm = T),
            sub_pstM = mean(ast_pst2_Bn, na.rm = T),
            sub_methM = mean(ast_meth2_Bn, na.rm = T),
            sub_opiM = mean(ast_opi_mean, na.rm = T),
            sub_opiPxM = mean(ast_opip2_Bn, na.rm = T),
            sub_herM = mean(ast_her2_Bn, na.rm = T),
            sub_fenM = mean(ast_fen2_Bn, na.rm = T),
            sub_narcM = mean(ast_narc2_Bn, na.rm = T),
            sub_inhM = mean(ast_inh2_Bn, na.rm = T),
            sub_sedM = mean(ast_sed2_Bn, na.rm = T)) %>%
  pivot_wider(id_cols = 'subNo',
              names_from = 'Lockdown', 
              values_from = 'sub_cannM':'sub_sedM') -> datT


# Conduct paired t-tests the change in each substance before and after COVID
# List of substances to loop through
substances <- c("cannM", "cocM", "fenM", "herM", "inhM", "methM", "narcM", "pstM", "sedM", "stimM")

# Initialize an empty dataframe to store the results
results_df <- data.frame(
  Substance = character(),
  Mean.Before = numeric(),
  Mean.After = numeric(),
  Statistic = numeric(),
  P.Value = numeric(),
  Confidence.Lower = numeric(),
  Confidence.Upper = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each substance and perform the paired t-test
for (substance in substances) {
  before <- datT[[paste0("sub_", substance, "_Before")]]
  after <- datT[[paste0("sub_", substance, "_After")]]
  test_result <- t.test(before, after, paired = TRUE)
  
  # Store the results in the dataframe
  results_df <- rbind(results_df, data.frame(
    Substance = substance,
    Mean.Before = mean(before, na.rm = T),
    Mean.After = mean(after, na.rm = T),
    Statistic = test_result$statistic,
    P.Value = test_result$p.value,
    Confidence.Lower = test_result$conf.int[1],
    Confidence.Upper = test_result$conf.int[2],
    stringsAsFactors = FALSE
  ))
}

