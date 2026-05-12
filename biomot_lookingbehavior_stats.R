# Clean the environment
rm(list = ls())

# Load required packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggeffects)
library(parameters)

# Read csv
bm_looking_data <- read_csv("//fs.univie.ac.at/homedirs/reiseb98/Documents/BioMot/final_analysis/BM_looking.csv") #enter pathway to BM_looking.csv file

# Adjust variable coding
bm_looking_data <- bm_looking_data %>% mutate(across(where(is.character), as.factor))
bm_looking_data$age_group <- as.factor(bm_looking_data$age_group)
bm_looking_data$ID <- as.factor(bm_looking_data$ID)


# Create df with per trial type averages
bm_looking_avg <- bm_looking_data %>%
  select(-trial_duration, -age_group)%>%
  group_by(ID,trial_type) %>%
  summarise(mean_looking_time = mean(looking_time),
            mean_looking_perc = mean(looking_perc))
# Add condition logic
bm_looking_avg <- bm_looking_avg %>%
  mutate(condition = case_when(
    trial_type == 'S 21' ~ 'upright',
    trial_type == 'S 23' ~ 'upright',
    trial_type == 'S 22' ~ 'inverted',
    trial_type == 'S 24' ~ 'inverted',
  ),
  walking_direction = case_when(
    trial_type == 'S 21' ~ 'right',
    trial_type == 'S 23' ~ 'left',
    trial_type == 'S 22' ~ 'right',
    trial_type == 'S 24' ~ 'left',
  ) ) 




# descriptives for looking per age group
bm_looking_data %>% 
  group_by(age_group) %>%
  summarise(mean_looking = mean(looking_time),
            sd_looking = sd(looking_time),
            md_looking = median(looking_time))

# descriptives for looking per age group and per condition 
bm_looking_data %>% 
  group_by(age_group, condition) %>%
  summarise(mean_looking = mean(looking_time),
            sd_looking = sd(looking_time),
            md_looking = median(looking_time))



# distribution of looking time
hist(bm_looking_data$looking_time) # data left skewed 

# log-transform looking times
bm_looking_data <- bm_looking_data %>% filter(looking_time > 0 )
bm_looking_data$looking_time_transf <- log(bm_looking_data$looking_time) 



# linear mixed effects model
# contrast (sum) coding for categorical variables 
bm_looking_data$condition <- factor(bm_looking_data$condition,levels = c("upright", "inverted"))
bm_looking_data$walking_direction <- factor(bm_looking_data$walking_direction,levels = c("right", "left"))
contrasts(bm_looking_data$condition) = contr.sum(2)
contrasts(bm_looking_data$walking_direction) = contr.sum(2)
contrasts(bm_looking_data$age_group) = contr.sum(2)


model <- lmer(looking_time_transf ~ condition*walking_direction + age_group + (1 |ID), data = bm_looking_data)
summary(model)

# model assumptions
plot(model, which = 1) # homoscedasticity
plot(resid(model))
qqnorm(resid(model)) # normal distribution of residuals 
qqline(resid(model))
car::vif(model) # multicollinearity


#for standardized parameters:
model_parameters(model, standardize = "refit")
