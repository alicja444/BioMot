# Clean the environment
rm(list = ls())

# Load required packages
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(MuMIn)
library(dplyr)
library(ggplot2)
library(ggdist)
library(ggridges)
library(ggsignif)
library(tidyverse)
library(broom)
library(mice)
library(car)
library(psych)
library(irr)
library(MASS)
library(boot)
library(lm.beta)
library(parameters)
library(ggeffects)
library(performance)
library(effects)
library(sjPlot)
library(interactions)
library(dplyr)
library(knitr)


DB <- read_excel("//fs.univie.ac.at/homedirs/reiseb98/Documents/BioMot/final_analysis/BM_SNRs.xlsx") #enter pathway to BM_SNRs.xlsx file
DB_touch <- read_excel("//fs.univie.ac.at/homedirs/reiseb98/Documents/BioMot/final_analysis/BM_touch.xlsx") #enter pathway to BM_touch.xlsx file


DB_touch_1row <- DB_touch %>%
  group_by(participant_ID) %>%
  summarise(
    age_group = first(na.omit(age_group)),
    STQ           = first(na.omit(STQ)),
    PICTS         = first(na.omit(PICTS)),
    PICTS_carry   = first(na.omit(PICTS_carry)),
    .groups = "drop"
  )

DB <- dplyr::left_join(DB, DB_touch_1row, by = "participant_ID")


DB$participant_ID <- as.factor(DB$participant_ID)
DB$age <- as.factor(DB$age)
DB$orientation <- as.factor(DB$orientation)
DB$direction <- as.factor(DB$direction)
DB$harmonic <- as.factor(DB$harmonic)
DB$hemisphere <- as.factor(DB$hemisphere)
DB$no_epochs <- as.numeric(DB$no_epochs)
DB$SNR <- as.numeric(DB$SNR)
DB$STQ <- as.numeric(DB$STQ)
DB$PICTS <- as.numeric(DB$PICTS)
DB$PICTS_carry <- as.numeric(DB$PICTS_carry)

DB$STQr <- (0 + 68) - DB$STQ #reverse-scoring the STQ score (min + max) - score
DB_touch_1row$STQr <- (0 + 68) - DB_touch_1row$STQ #reverse-scoring the STQ score (min + max) - score


DB$SNR_centered <- DB$SNR - mean(DB$SNR, na.rm = TRUE)
DB$STQ_centered <- DB$STQ - mean(DB$STQ, na.rm = TRUE)
DB$STQr_centered <- DB$STQr - mean(DB$STQr, na.rm = TRUE)
DB$PICTS_centered <- DB$PICTS - mean(DB$PICTS, na.rm = TRUE)
DB$PICTS_carry_centered <- DB$PICTS_carry - mean(DB$PICTS_carry, na.rm = TRUE)
DB$no_epochs_centered <- DB$no_epochs - mean(DB$no_epochs, na.rm = TRUE)



#correlations between touch measures

sub <- DB_touch_1row[, 4:6, drop = FALSE]
ct <- corr.test(sub, use = "pairwise", method = "pearson")
ct$r        # correlations
ct$p        # adjusted p-values (Benjamini–Hochberg here)
ct$n        # pairwise N
ct$ci       # confidence intervals

# with significance levels
r <- round(ct$r, 2)
p <- ct$p
rstar <- matrix(paste0(format(r, nsmall = 2),
                       ifelse(p < 0.001, "***",
                              ifelse(p < 0.01,  "**",
                                     ifelse(p < 0.05,  "*",
                                            ifelse(p < 0.1,   "·", ""))))),
                nrow = nrow(r), dimnames = dimnames(r))
rstar[upper.tri(rstar, diag = TRUE)] <- ""  # keep lower triangle
kable(rstar, caption = "Psych::corr.test with BH-adjusted p-values")



# descriptive statistics (Table 1)

per_ID_averages <- DB %>%
  filter(!(harmonic == '7_2Hz'))%>%
  group_by(participant_ID) %>%
  mutate(SNR_summed = SNR + lead(SNR)) %>%  # add next row first
  filter(row_number() %% 2 == 1) %>%  # then keep only odd rows
  ungroup() %>%
  group_by(participant_ID, orientation) %>%
  summarise(avg_SNR = mean(SNR_summed, na.rm = TRUE),
            age = first(age)
            )%>%
  ungroup() 

SNR_averages <- per_ID_averages %>%
  group_by(age, orientation) %>%
  summarise(n = n_distinct(participant_ID),
            mean_SNR = mean(avg_SNR, na.rm = TRUE),
            sd_SNR   = sd(avg_SNR, na.rm = TRUE),
            min_SNR  = min(avg_SNR, na.rm = TRUE),
            max_SNR  = max(avg_SNR, na.rm = TRUE),
            .groups = "drop")

questionnaire_averages <- DB %>%
  group_by(age) %>%
  drop_na(PICTS, PICTS_carry, STQr)%>%
  summarise(n = n_distinct(participant_ID),
            across(c(PICTS, PICTS_carry, STQr), 
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x, na.rm = TRUE),
                        min  = \(x) min(x, na.rm = TRUE),
                        max  = \(x) max(x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")


#test significance of harmonics

#3 months, 2.4 Hz

data3M_2_4 <- DB %>%
  filter(
    age == "3",
    harmonic == "2_4Hz",
    orientation %in% c("up", "inverted")
  ) %>%
  group_by(participant_ID, orientation) %>%
  summarise(SNR = mean(SNR, na.rm = TRUE), .groups = "drop")

t_up <- data3M_2_4 %>%
  filter(orientation == "up") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# For 'inverted' condition
t_inverted <- data3M_2_4 %>%
  filter(orientation == "inverted") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# Print results
t_up
t_inverted

#3 months, 4.8 Hz

data3M_4_8 <- DB %>%
  filter(
    age == "3",
    harmonic == "4_8Hz",
    orientation %in% c("up", "inverted")
  ) %>%
  group_by(participant_ID, orientation) %>%
  summarise(SNR = mean(SNR, na.rm = TRUE), .groups = "drop")

# For 'up' condition
t_up2 <- data3M_4_8 %>%
  filter(orientation == "up") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# For 'inverted' condition
t_inverted2 <- data3M_4_8 %>%
  filter(orientation == "inverted") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# Print results
t_up2
t_inverted2

#3 months, 7.2 Hz

data3M_7_2 <- DB %>%
  filter(
    age == "3",
    harmonic == "7_2Hz",
    orientation %in% c("up", "inverted")
  ) %>%
  group_by(participant_ID, orientation) %>%
  summarise(SNR = mean(SNR, na.rm = TRUE), .groups = "drop")


# For 'up' condition
t_up3 <- data3M_7_2 %>%
  filter(orientation == "up") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# For 'inverted' condition
t_inverted3 <- data3M_7_2 %>%
  filter(orientation == "inverted") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# Print results
t_up3
t_inverted3

#6 months, 2.4 Hz

data6M_2_4 <- DB %>%
  filter(
    age == "6",
    harmonic == "2_4Hz",
    orientation %in% c("up", "inverted")
  ) %>%
  group_by(participant_ID, orientation) %>%
  summarise(SNR = mean(SNR, na.rm = TRUE), .groups = "drop")


# For 'up' condition
t_up4 <- data6M_2_4 %>%
  filter(orientation == "up") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# For 'inverted' condition
t_inverted4 <- data6M_2_4 %>%
  filter(orientation == "inverted") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# Print results
t_up4
t_inverted4

#6 months, 4.8 Hz

data6M_4_8 <- DB %>%
  filter(
    age == "6",
    harmonic == "4_8Hz",
    orientation %in% c("up", "inverted")
  ) %>%
  group_by(participant_ID, orientation) %>%
  summarise(SNR = mean(SNR, na.rm = TRUE), .groups = "drop")

# For 'up' condition
t_up5 <- data6M_4_8 %>%
  filter(orientation == "up") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# For 'inverted' condition
t_inverted5 <- data6M_4_8 %>%
  filter(orientation == "inverted") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# Print results
t_up5
t_inverted5

data6M_7_2 <- DB %>%
  filter(
    age == "6",
    harmonic == "7_2Hz",
    orientation %in% c("up", "inverted")
  ) %>%
  group_by(participant_ID, orientation) %>%
  summarise(SNR = mean(SNR, na.rm = TRUE), .groups = "drop")

# For 'up' condition
t_up6 <- data6M_7_2 %>%
  filter(orientation == "up") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# For 'inverted' condition
t_inverted6 <- data6M_7_2 %>%
  filter(orientation == "inverted") %>%
  pull(SNR) %>%
  t.test(mu = 0)

# Print results
t_up6
t_inverted6




DB <- subset(DB, harmonic != '7_2Hz')

#check if number of segments similar between conditions:

model_test <- lmer(no_epochs ~orientation +  (1|participant_ID), data = DB)
summary(model_test) #not significant: same number of epochs for both conditions
model_parameters(model_test, standardize = "refit")

#preregistered analyses
#1 H1: Neural entrainment to biological motion is observable at 3 months and at 6 months, but stronger at 6 months.
model1 <- lmer(SNR ~ orientation*age + (1|participant_ID), data = DB)
summary(model1)
vif(model1) #no multicollinearity problem
#for standardized parameters:
model_parameters(model1, standardize = "refit")


#No evidence. No main effect of orientation or age, no interaction.


#2 H2: Levels of infant carrying and other infant-directed touching behaviors reported by the primary caregiver(s) will be positively correlated with neural entrainment to biological motion at 3 and at 6 months.

model2 <- lmer(SNR ~ orientation*PICTS + orientation*PICTS_carry + orientation*STQr + (1|participant_ID), data = DB)
summary(model2)
vif(model2) #problem with multicollinearity: use centered variables

model2c <- lmer(SNR_centered ~ orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model2c)
model_parameters(model2c, standardize = "refit")

#check interactions with age?
#model2ca <- lmer(SNR_centered ~ orientation*PICTS_centered*age + orientation*PICTS_carry_centered*age + orientation*STQr_centered*age + (1|participant_ID), data = DB, REML = FALSE)
#summary(model2ca)


vif(model2c)

p2 <- ggpredict(model2, terms = "PICTS_carry")  # marginalizes over orientation
plot(p2) + labs(title = "Main (marginal) effect of PICTS_carry")

p3 <- ggpredict(model2, terms = "STQr")  # marginalizes over orientation
plot(p3) + labs(title = "Main (marginal) effect of STQ")


p4 <- ggpredict(model2, terms = c("PICTS_carry", "orientation"), length = 100)
plot(p4) +
  labs(title = " ", 
       x = "Infant Carrying", y = expression(paste("Amplitude [", mu, "V", "]"))) +
  scale_fill_manual(values = c("#D41159", "#1A85FF")) +
  scale_color_manual(values = c("#D41159", "#1A85FF"), breaks = c("inverted", "up"),
                     labels = c("inverted", "upright")) +
  theme(
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

  
p5 <- ggpredict(model2, terms = c("STQr", "orientation"), length = 100)
plot(p5) +
  labs(title = " ", 
       x = "Social Touch Questionnaire - reversed", y = expression(paste("Amplitude [", mu, "V", "]"))) +
  scale_fill_manual(values = c("#D41159", "#1A85FF")) +
  scale_color_manual(values = c("#D41159", "#1A85FF"), breaks = c("inverted", "up"),
                     labels = c("inverted", "upright")) +
  theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20), 
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20))


#exploratory analysis follows

#building the best "exploratory" model
#model2c is now the base model M0 to which we will add additional fixed effects of interest

#add three-way interaction orientation*PICTS_carry_centered*age
model1e <- lmer(SNR_centered ~ orientation*PICTS_centered + orientation*PICTS_carry_centered*age + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model1e)
vif(model1e) #no problem?

anova(model2c,model1e) #model2c better  

#add three-way interaction orientation*STQr_centered*age
model2e <- lmer(SNR_centered ~ orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered*age + (1|participant_ID), data = DB, REML = FALSE)
summary(model2e)
vif(model2e) #no problem

anova(model2c,model2e) #model2c better

#add number of available segments
model3e <- lmer(SNR_centered ~ no_epochs_centered + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model3e)
vif(model3e) #no problem

anova(model2c,model3e) #model2c better

#add direction to the model
model4e <- lmer(SNR_centered ~ direction + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model4e)
vif(model4e)

anova(model2c,model4e) #model4e better! becomes new bases to add effects to

#add direction*orientation interaction
model5e <- lmer(SNR_centered ~ direction*orientation + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model5e)
vif(model5e)

anova(model4e,model5e) #model4e better

#add interaction between direction and age
model6e <- lmer(SNR_centered ~ direction*age + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model6e)
vif(model6e)

anova(model4e,model6e) #model4e better

#add main effect of hemisphere
model7e <- lmer(SNR_centered ~ direction + hemisphere + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model7e)
vif(model7e)

anova(model4e,model7e) #model4e better

#add interaction of orientation with hemisphere
model8e <- lmer(SNR_centered ~ direction + orientation*hemisphere + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model8e)
vif(model8e)

anova(model4e,model8e) #model4e better

#add interaction of direction with hemisphere
model9e <- lmer(SNR_centered ~ direction*hemisphere + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model9e)
vif(model9e)

anova(model4e,model9e) #model9e better! BEST MODEL

#add three-way interaction direction*hemisphere*age
model10e <- lmer(SNR_centered ~ direction*hemisphere*age + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model10e)
vif(model10e)

anova(model9e,model10e) #model9e better! BEST MODEL

#add three-way interaction direction*hemisphere*orientation
model11e <- lmer(SNR_centered ~ direction*hemisphere*orientation + orientation*PICTS_centered + orientation*PICTS_carry_centered + orientation*STQr_centered + (1|participant_ID), data = DB, REML = FALSE)
summary(model11e)
vif(model11e)

anova(model9e,model11e) #model9e better! BEST MODEL


#for standardized parameters:
model_parameters(model9e, standardize = "refit")

p9e <- ggpredict(model9e, terms = c("hemisphere", "direction"), length = 100)
plot(p9e) +
  labs(title = " ", 
       x = "Hemisphere", y = expression(paste("Amplitude [", mu, "V", "]"))) +
  scale_fill_manual(values = c("#1AFF1A", "#4B0092")) +
  scale_color_manual(values = c("#1AFF1A", "#4B0092")) +
  theme(
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20))

# Effect of hemisphere within each level of direction
pairs(emmeans(model9e, ~ hemisphere | direction), adjust = "holm")

# Effect of direction within each level of hemisphere
pairs(emmeans(model9e, ~ direction | hemisphere), adjust = "holm")



p9e <- ggpredict(model9e, terms = c("hemisphere", "direction"), length = 100)
plot(p9e) +
  labs(title = "Interaction: direction × hemisphere",
       x = "hemisphere", y = "Predicted SNR")

p9e2 <- ggpredict(model9e, terms = c("PICTS_carry_centered", "orientation"), length = 100)
plot(p9e2) +
  labs(title = "Interaction: direction × hemisphere",
       x = "PICTS_carry_centered", y = "Predicted SNR")

p9e3 <- ggpredict(model9e, terms = c("STQr_centered", "orientation"), length = 100)
plot(p9e3) +
  labs(title = "Interaction: direction × hemisphere",
       x = "STQr_centered", y = "Predicted SNR")



tr_carry <- emtrends(model9e, ~ orientation, var = "PICTS_carry_centered")
pairs(tr_carry, adjust = "tukey")

tr_STQ <- emtrends(model9e, ~ orientation, var = "STQr_centered")
pairs(tr_STQ, adjust = "tukey")

emm <- emmeans(model9e, ~ direction*hemisphere)
emm



