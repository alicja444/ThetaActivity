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
library(gghalves)


DB <- read_excel("/Users/alicja/Documents/manuscripts/Frontal theta modulation/stats_consultation/dataset.xlsx")
DB_long <- read_excel("/Users/alicja/Documents/manuscripts/Frontal theta modulation/stats_consultation/dataset_long.xlsx")


#Variables preparation

DB$ID <- as.factor(DB$ID)
DB$modulation_indexT1theta <- as.numeric(DB$modulation_indexT1theta)
DB$modulation_r1T1theta <- as.numeric(DB$modulation_r1T1theta)
DB$modulation_r2T1theta <- as.numeric(DB$modulation_r2T1theta)
DB$number_of_segmentsT1 <- as.numeric(DB$number_of_segmentsT1)
DB$absolute_powerT1theta <- as.numeric(DB$absolute_powerT1theta)
DB$relative_powerT1theta <- as.numeric(DB$relative_powerT1theta)
DB$modulation_indexT1alpha <- as.numeric(DB$modulation_indexT1alpha)
DB$absolute_powerT1alpha <- as.numeric(DB$absolute_powerT1alpha)
DB$relative_powerT1alpha <- as.numeric(DB$relative_powerT1alpha)
DB$modulation_indexT2theta <- as.numeric(DB$modulation_indexT2theta)
DB$modulation_r1T2theta <- as.numeric(DB$modulation_r1T2theta)
DB$modulation_r2T2theta <- as.numeric(DB$modulation_r2T2theta)
DB$number_of_segmentsT2 <- as.numeric(DB$number_of_segmentsT2)
DB$absolute_powerT2theta <- as.numeric(DB$absolute_powerT2theta)
DB$relative_powerT2theta <- as.numeric(DB$relative_powerT2theta)
DB$modulation_indexT2alpha <- as.numeric(DB$modulation_indexT2alpha)
DB$absolute_powerT2alpha <- as.numeric(DB$absolute_powerT2alpha)
DB$relative_powerT2alpha <- as.numeric(DB$relative_powerT2alpha)
DB$Bayley <- as.numeric(DB$Bayley)
DB$Visual_search <- as.numeric(DB$Visual_search)
DB$Elfra_prod <- as.numeric(DB$Elfra_prod)
DB$Elfra_syntax <- as.numeric(DB$Elfra_syntax)
DB$Elfra_morph <- as.numeric(DB$Elfra_morph)
DB$Age_at_T3 <- as.numeric(DB$Age_at_T3)


#Variables preparation LONG

DB_long$ID <- as.factor(DB_long$ID)
DB_long$timepoint <- as.factor(DB_long$timepoint)
DB_long$theta_modulation <- as.numeric(DB_long$theta_modulation)
DB_long$theta_modr1 <- as.numeric(DB_long$theta_modr1)
DB_long$number_of_segments <- as.numeric(DB_long$number_of_segments)
DB_long$absolute_theta <- as.numeric(DB_long$absolute_theta)
DB_long$relative_theta <- as.numeric(DB_long$relative_theta)
DB_long$alpha_modulation <- as.numeric(DB_long$alpha_modulation)
DB_long$absolute_alpha <- as.numeric(DB_long$absolute_alpha)
DB_long$relative_alpha <- as.numeric(DB_long$relative_alpha)
DB_long$Bayley <- as.numeric(DB_long$Bayley)
DB_long$Visual_search <- as.numeric(DB_long$Visual_search)
DB_long$Elfra_prod <- as.numeric(DB_long$Elfra_prod)
DB_long$Elfra_syntax <- as.numeric(DB_long$Elfra_syntax)
DB_long$Elfra_morph <- as.numeric(DB_long$Elfra_morph)
DB_long$Age_at_T3 <- as.numeric(DB_long$Age_at_T3)

#Can we observe overall increase in theta power over the course of video viewing?

#6 months:
t.test(DB$modulation_indexT1theta, mu = 0)
#12 months:
t.test(DB$modulation_indexT2theta, mu = 0)


#Examining the reliability of theta modulation at 6 and 12 months

icc_dataT1 <- DB[, c("modulation_r1T1theta", "modulation_r2T1theta")]
icc_dataT2 <- DB[, c("modulation_r1T2theta", "modulation_r2T2theta")]

icc_resultT1 <- icc(icc_dataT1, model = "twoway", type = "agreement", unit = "single")
print(icc_resultT1)

icc_resultT2 <- icc(icc_dataT2, model = "twoway", type = "agreement", unit = "single")
print(icc_resultT2)

#Change in absolute theta power from 6 to 12 months

t.test(DB$absolute_powerT2theta,DB$absolute_powerT1theta, paired = TRUE, alternative = "two.sided")

#Change in relative theta from 6 to 12 months
t.test(DB$relative_powerT2theta,DB$relative_powerT1theta, paired = TRUE, alternative = "two.sided")

#Change in theta modulation index from 6 to 12 months
t.test(DB$modulation_indexT2theta, DB$modulation_indexT1theta, paired = TRUE, alternative = "two.sided")

#correlation between 6 month and 12 month absolute theta power
cor.test(DB$absolute_powerT2theta,DB$absolute_powerT1theta, use = "complete.obs", method = "pearson")

#correlation between 6 month and 12 month relative theta power
cor.test(DB$relative_powerT2theta,DB$relative_powerT1theta, use = "complete.obs", method = "pearson")

#correlation between 6 month and 12 month theta modulation index
cor.test(DB$modulation_indexT1theta,DB$modulation_indexT2theta, use = "complete.obs", method = "pearson")

# plots for 6 to 12 month changes
# Plot 1: Absolute theta change from 6 months to 12 months

DB_clean <- DB_long %>%
  filter(!is.na(absolute_theta))
labels_x <- c("6 months", "12 months")

plot1 <- ggplot(DB_clean, aes(x = timepoint, y = absolute_theta, fill = timepoint)) +
  geom_half_violin(side = "l", alpha = 0.6, trim = FALSE, color = NA) +  # Half violin on left
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +  # Boxplot
  geom_jitter(aes(color = timepoint), width = 0.1, size = 2, alpha = 0.6) +  # Raw data
  geom_line(aes(group = ID), color = "black", alpha = 0.3) +  # Line for paired data
  theme_minimal() +
  labs(title = " ", x = " ", y = expression(paste("Absolute theta power [", mu, "V"^2, "]"))) + 
  scale_fill_manual(values = c("plum", "olivedrab3")) +
  scale_color_manual(values = c("plum", "olivedrab3")) +  # For jitter points
  scale_x_discrete(labels = labels_x) +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))

plot1

#ggsave("absolutetheta_rain.png", plot = plot1, width = 6, height = 6, dpi = 300)

#Plot 2: Theta modulation indext at 6 and 12 months

DB_clean2 <- DB_long %>%
  filter(!is.na(theta_modulation))
labels_x <- c("6 months", "12 months")

plot2 <- ggplot(DB_clean2, aes(x = timepoint, y = theta_modulation, fill = timepoint)) +
  geom_half_violin(side = "l", alpha = 0.6, trim = FALSE, color = NA) +  # Half violin on left
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +  # Boxplot
  geom_jitter(aes(color = timepoint), width = 0.1, size = 2, alpha = 0.6) +  # Raw data
  geom_line(aes(group = ID), color = "black", alpha = 0.3) +  # Line for paired data
  theme_minimal() +
  labs(title = " ", x = " ", y = "Theta modulation index") + 
  scale_fill_manual(values = c("plum", "olivedrab3")) +
  scale_color_manual(values = c("plum", "olivedrab3")) +  # For jitter points
  scale_x_discrete(labels = labels_x) +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))

plot2
#ggsave("theta_modulation_rain.png", plot = plot2, width = 6, height = 6, dpi = 300)


#Plot 3: relative theta
DB_clean3 <- DB_long %>%
  filter(!is.na(relative_theta))
labels_x <- c("6 months", "12 months")

plot3 <- ggplot(DB_clean3, aes(x = timepoint, y = relative_theta, fill = timepoint)) +
  geom_half_violin(side = "l", alpha = 0.6, trim = FALSE, color = NA) +  # Half violin on left
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +  # Boxplot
  geom_jitter(aes(color = timepoint), width = 0.1, size = 2, alpha = 0.6) +  # Raw data
  geom_line(aes(group = ID), color = "black", alpha = 0.3) +  # Line for paired data
  theme_minimal() +
  labs(title = " ", x = " ", y = "Relative theta power") + 
  scale_fill_manual(values = c("plum", "olivedrab3")) +
  scale_color_manual(values = c("plum", "olivedrab3")) +  # For jitter points
  scale_x_discrete(labels = labels_x) +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))

plot3
#ggsave("relative_theta_rain.png", plot = plot3, width = 6, height = 6, dpi = 300)

#correlation between absolute theta at 6 and 12 months

plot4 <- ggplot(DB, aes(x = absolute_powerT1theta, y = absolute_powerT2theta)) +
  geom_point() +  # Plot the data points
  geom_smooth(method = "lm", col = "grey10") +  # Add the regression line
  labs(title = " ",
       x = expression(paste("Absolute theta power at 6 months [", mu, "V"^2, "]")),
       y = expression(paste("Absolute theta power at 12 months [", mu, "V"^2, "]"))) +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))

plot4
#ggsave("absolute_th_correlationnew.png", plot = plot4, width = 6, height = 6, dpi = 300)

# plot5:  correlation between theta modulation at 6 and 12 months
plot5 <- ggplot(DB, aes(x = modulation_indexT1theta, y = modulation_indexT2theta)) +
  geom_point() +  # Plot the data points
  geom_smooth(method = "lm", col = "grey10") +  # Add the regression line
  labs(title = " ",
       x = "Theta modulation index at 6 months",
       y = "Theta modulation index at 12 months") +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))

plot5
#ggsave("modulation_correlation.png", plot = plot5, width = 6, height = 6, dpi = 300)

# plot6: correlation between relative theta  at 6 and 12 months
plot6 <- ggplot(DB, aes(x = relative_powerT1theta, y = relative_powerT2theta)) +
  geom_point() +  # Plot the data points
  geom_smooth(method = "lm", col = "grey10") +  # Add the regression line
  labs(title = " ",
       x = "Relative theta power at 6 months",
       y = "Relative theta power at 12 months") +
  theme_minimal() + 
  xlim(0.4,0.62) +
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20))

plot6
#ggsave("relative_correlationnew.png", plot = plot6, width = 6, height = 6, dpi = 300)


#add a composite language score

# calculate z-scores for prod, morph & syntax
DB$z_Elfraprod <- as.numeric(scale(DB$Elfra_prod))
DB$z_Elframorph <- as.numeric(scale(DB$Elfra_morph))
DB$z_Elfrasyntax <- as.numeric(scale(DB$Elfra_syntax))
DB$Language <- rowMeans(DB[, c("z_Elfraprod", "z_Elframorph", "z_Elfrasyntax")], na.rm = TRUE)

#add composite abs theta, theta modul & rel theta scores

DB$z_absthetaT1 <- as.numeric(scale(DB$absolute_powerT1theta))
DB$z_absthetaT2 <- as.numeric(scale(DB$absolute_powerT2theta))
DB$absthetaz <- rowMeans(DB[, c("z_absthetaT1", "z_absthetaT2")], na.rm = TRUE)

DB$z_thetamodulT1 <- as.numeric(scale(DB$modulation_indexT1theta))
DB$z_thetamodulT2 <- as.numeric(scale(DB$modulation_indexT2theta))
DB$thetamodulz <- rowMeans(DB[, c("z_thetamodulT1", "z_thetamodulT2")], na.rm = TRUE)

DB$z_relthetaT1 <- as.numeric(scale(DB$relative_powerT1theta))
DB$z_relthetaT2 <- as.numeric(scale(DB$relative_powerT2theta))
DB$relthetaz <- rowMeans(DB[, c("z_relthetaT1", "z_relthetaT2")], na.rm = TRUE)

#add composite abs alphaa, alpha modul & rel alpha scores
DB$z_absalphaT1 <- as.numeric(scale(DB$absolute_powerT1alpha))
DB$z_absalphaT2 <- as.numeric(scale(DB$absolute_powerT2alpha))
DB$absalphaz <- rowMeans(DB[, c("z_absalphaT1", "z_absalphaT2")], na.rm = TRUE)

DB$z_alphamodulT1 <- as.numeric(scale(DB$modulation_indexT1alpha))
DB$z_alphamodulT2 <- as.numeric(scale(DB$modulation_indexT2alpha))
DB$alphamodulz <- rowMeans(DB[, c("z_alphamodulT1", "z_alphamodulT2")], na.rm = TRUE)

DB$z_relalphaT1 <- as.numeric(scale(DB$relative_powerT1alpha))
DB$z_relalphaT2 <- as.numeric(scale(DB$relative_powerT2alpha))
DB$relalphaz <- rowMeans(DB[, c("z_relalphaT1", "z_relalphaT2")], na.rm = TRUE)

#edit datasets for differing sample sizes (with the composite scores):

DB_subset_Bayley <- DB[, c("ID", "absthetaz", "thetamodulz", "relthetaz", 
                           "Bayley", "Age_at_T3")]
DB_clean_Bayley <- na.omit(DB_subset_Bayley)

DB_subset_VSearch <- DB[, c("ID", "absthetaz", "thetamodulz", "relthetaz", 
                            "Age_at_T3","Visual_search")]
DB_clean_VSearch <- na.omit(DB_subset_VSearch)

DB_subset_Language <- DB[, c("ID", "absthetaz", "thetamodulz", "relthetaz",  
                             "Age_at_T3","Language")]
DB_clean_Language <- na.omit(DB_subset_Language)

#now for alpha:

#Can we observe overall increase in alpha power over the course of video viewing?

#6 months:
t.test(DB$modulation_indexT1alpha, mu = 0)
#12 months:
t.test(DB$modulation_indexT2alpha, mu = 0)

#Change in absolute alpha power from 6 to 12 months

t.test(DB$absolute_powerT2alpha,DB$absolute_powerT1alpha, paired = TRUE, alternative = "two.sided")

#Change in relative alpha from 6 to 12 months
t.test(DB$relative_powerT2alpha,DB$relative_powerT1alpha, paired = TRUE, alternative = "two.sided")

#Change in alpha modulation index from 6 to 12 months
t.test(DB$modulation_indexT2alpha, DB$modulation_indexT1alpha, paired = TRUE, alternative = "two.sided")

#correlation between 6 month and 12 month absolute alpha power
cor.test(DB$absolute_powerT2alpha,DB$absolute_powerT1alpha, use = "complete.obs", method = "pearson")

#correlation between 6 month and 12 month relative theta power
cor.test(DB$relative_powerT2alpha,DB$relative_powerT1alpha, use = "complete.obs", method = "pearson")

#correlation between 6 month and 12 month theta modulation index
cor.test(DB$modulation_indexT1alpha,DB$modulation_indexT2alpha, use = "complete.obs", method = "pearson")

#predicting cognitive outcomes with alpha activity

DB_subset_BayleyA <- DB[, c("ID", "absalphaz", "alphamodulz", "relalphaz", 
                           "Bayley", "Age_at_T3")]
DB_clean_BayleyA <- na.omit(DB_subset_BayleyA)

DB_subset_VSearchA <- DB[, c("ID", "absalphaz", "alphamodulz", "relalphaz", 
                            "Age_at_T3","Visual_search")]
DB_clean_VSearchA <- na.omit(DB_subset_VSearchA)

DB_subset_LanguageA <- DB[, c("ID", "absalphaz", "alphamodulz", "relalphaz",  
                             "Age_at_T3","Language")]
DB_clean_LanguageA <- na.omit(DB_subset_LanguageA)



#predicting Bayley with composite scores:
full_model_Bayley <- lm(Bayley ~ absthetaz + thetamodulz + relthetaz + 
                          Age_at_T3, data = DB_clean_Bayley)

best_model_Bayley <- stepAIC(full_model_Bayley, direction = "both")

summary(best_model_Bayley)

model_std <- standardize_parameters(best_model_Bayley)
print(model_std)


# visualize the relationship between composite absolute theta & Bayley

# Fit the full model
model <- lm(Bayley ~ absthetaz + Age_at_T3, data = DB_clean_Bayley)

# Get residuals from partial regression
resid_y <- residuals(lm(Bayley ~ Age_at_T3, data = DB_clean_Bayley))
resid_x <- residuals(lm(absthetaz ~ Age_at_T3, data = DB_clean_Bayley))

# Create a dataframe for plotting
partial_data <- data.frame(resid_x, resid_y)

# Plot partial regression (added variable plot)
plotxx <- ggplot(partial_data, aes(x = resid_x, y = resid_y)) +
  geom_point() +
  geom_smooth(method = "lm", col = "grey10") +
  labs(title = " ",
       x = "Residuals of the Absolute theta composite score",
       y = "Residuals of Bayley-III Cognitive Index at 24 months") +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16))
plotxx
#ggsave("abstheta_bayleynew.png", plot = plotxx, width = 6, height = 6, dpi = 300)


#now for visual search
#predicting Visual search with composite scores:
full_model_VSearch <- lm(Visual_search ~ absthetaz + thetamodulz + relthetaz + 
                          Age_at_T3, data = DB_clean_VSearch)

best_model_VSearch <- stepAIC(full_model_VSearch, direction = "both")

summary(best_model_VSearch) #null model

#now for language

full_model_Language <- lm(Language ~ absthetaz + thetamodulz + relthetaz + 
                           Age_at_T3, data = DB_clean_Language)

best_model_Language <- stepAIC(full_model_Language, direction = "both")

summary(best_model_Language) 

model_std11 <- standardize_parameters(best_model_Language)
print(model_std11)


#now with T1 & T2 variables separately
#Bayley

DB_subset_Bayley2 <- DB[, c("ID", "modulation_indexT1theta", "modulation_indexT2theta", "absolute_powerT1theta", 
                           "absolute_powerT2theta",  
                           "Bayley", "relative_powerT1theta", "relative_powerT2theta", "Age_at_T3")]
DB_clean_Bayley2 <- na.omit(DB_subset_Bayley2)


full_model_Bayley2 <- lm(Bayley ~ modulation_indexT1theta + modulation_indexT2theta + absolute_powerT1theta +
                           absolute_powerT2theta +
                           relative_powerT1theta + relative_powerT2theta + 
                           Age_at_T3, data = DB_clean_Bayley2)

best_model_Bayley2 <- stepAIC(full_model_Bayley2, direction = "both")

summary(best_model_Bayley2)

model_std2 <- standardize_parameters(best_model_Bayley2)
print(model_std2)

#Visual search
DB_subset_VSearch2 <- DB[, c("ID", "modulation_indexT1theta", "modulation_indexT2theta", "absolute_powerT1theta", "absolute_powerT2theta",
                            "Age_at_T3", "relative_powerT1theta", "relative_powerT2theta","Visual_search")]
DB_clean_VSearch2 <- na.omit(DB_subset_VSearch2)

full_model_VSearch2 <- lm(Visual_search ~ modulation_indexT1theta + modulation_indexT2theta + absolute_powerT1theta +
                           absolute_powerT2theta + 
                           relative_powerT1theta + relative_powerT2theta + 
                           Age_at_T3, data = DB_clean_VSearch2)

best_model_VSearch2 <- stepAIC(full_model_VSearch2, direction = "both")

summary(best_model_VSearch2)

model_std3 <- standardize_parameters(best_model_VSearch2)
print(model_std3)

#Language

DB_subset_Language2 <- DB[, c("ID", "modulation_indexT1theta", "modulation_indexT2theta", "absolute_powerT1theta", "absolute_powerT2theta",
                             "Age_at_T3", "relative_powerT1theta", "relative_powerT2theta","Language")]
DB_clean_Language2 <- na.omit(DB_subset_Language2)

full_model_Language2 <- lm(Language ~ modulation_indexT1theta + modulation_indexT2theta + absolute_powerT1theta +
                            absolute_powerT2theta + 
                            relative_powerT1theta + relative_powerT2theta + 
                            Age_at_T3, data = DB_clean_Language2)

best_model_Language2 <- stepAIC(full_model_Language2, direction = "both")

summary(best_model_Language2)

model_std4 <- standardize_parameters(best_model_Language2)
print(model_std4)

#redoing analyses with larger datasets (including only T1/T2 data for the respective models)
#Bayley

DB_subset_Bayley4 <- DB[, c("ID", "modulation_indexT1theta", "absolute_powerT1theta", 
                            "Bayley",  "Age_at_T3")]
DB_clean_Bayley4 <- na.omit(DB_subset_Bayley4)

chosen_model_Bayley <- lm(Bayley ~ modulation_indexT1theta + absolute_powerT1theta +
                           Age_at_T3, data = DB_clean_Bayley4)

summary(chosen_model_Bayley)

#Visual search
DB_subset_VSearch4 <- DB[, c("ID", "modulation_indexT1theta", 
                            "Visual_search")]
DB_clean_VSearch4 <- na.omit(DB_subset_VSearch4)

chosen_model_VSearch <- lm(Visual_search ~ modulation_indexT1theta, data = DB_clean_VSearch4)
summary(chosen_model_VSearch)

#Language

DB_subset_Language4 <- DB[, c("ID",  "absolute_powerT2theta", 
                            "Language",  "Age_at_T3")]
DB_clean_Language4 <- na.omit(DB_subset_Language4)

chosen_model_Language <- lm(Language ~ absolute_powerT2theta + 
                            Age_at_T3, data = DB_clean_Language4)

summary(chosen_model_Language)

#for the alpha frequency band
#predicting Bayley with composite scores
full_model_BayleyA <- lm(Bayley ~ absalphaz + alphamodulz + relalphaz + 
                           Age_at_T3, data = DB_clean_BayleyA)

best_model_BayleyA <- stepAIC(full_model_BayleyA, direction = "both")

summary(best_model_BayleyA)

#visual search
full_model_VSearchA <- lm(Visual_search ~ absalphaz + alphamodulz + relalphaz + 
                           Age_at_T3, data = DB_clean_VSearchA)

best_model_VSearchA <- stepAIC(full_model_VSearchA, direction = "both")

summary(best_model_VSearchA) #null model

#Language:
full_model_LanguageA <- lm(Language ~ absalphaz + alphamodulz + relalphaz + 
                            Age_at_T3, data = DB_clean_LanguageA)

best_model_LanguageA <- stepAIC(full_model_LanguageA, direction = "both")

summary(best_model_LanguageA) 


#analyses: change from T1 to T2 as a predictor (Supplementary material)

DB$change_abstheta <- DB$absolute_powerT2theta - DB$absolute_powerT1theta #increase
DB$change_thetamodul <- DB$modulation_indexT1theta - DB$modulation_indexT2theta #decrease
DB$change_reltheta <- DB$relative_powerT1theta - DB$relative_powerT2theta #decrease

#regressing Change on initial score
model_change <- lm(change_abstheta ~ absolute_powerT1theta, data = DB, na.action = na.exclude)
model_change2 <- lm(change_thetamodul ~ modulation_indexT1theta, data = DB, na.action = na.exclude)
model_change3 <- lm(change_reltheta ~ relative_powerT1theta, data = DB, na.action = na.exclude)

# Extract residuals (adjusted change scores)
DB$adjchange_abstheta <- residuals(model_change)
DB$adjchange_modultheta <- residuals(model_change2)
DB$adjchange_reltheta <- residuals(model_change3)

DB_subset_Bayley3 <- DB[, c("ID", 
                           "Bayley",  
                           "adjchange_abstheta", "adjchange_modultheta", "adjchange_reltheta",
                           "Age_at_T3")]
DB_clean_Bayley3 <- na.omit(DB_subset_Bayley3)


full_model_Bayley3 <- lm(Bayley ~ 
                          adjchange_abstheta + adjchange_modultheta + adjchange_reltheta +
                          Age_at_T3, data = DB_clean_Bayley3)

best_model_Bayley3 <- stepAIC(full_model_Bayley3, direction = "both")

summary(best_model_Bayley3)

model_std5 <- standardize_parameters(best_model_Bayley3)
print(model_std5)

DB_subset_VSearch3 <- DB[, c("ID", 
                            "Visual_search",  
                            "adjchange_abstheta", "adjchange_modultheta", "adjchange_reltheta",
                            "Age_at_T3")]
DB_clean_VSearch3 <- na.omit(DB_subset_VSearch3)

full_model_VSearch3 <- lm(Visual_search ~ 
                           adjchange_abstheta + adjchange_modultheta + adjchange_reltheta +
                           Age_at_T3, data = DB_clean_VSearch3)

best_model_VSearch3 <- stepAIC(full_model_VSearch3, direction = "both")

summary(best_model_VSearch3)

model_std55 <- standardize_parameters(best_model_VSearch3)
print(model_std55)

DB_subset_Language3 <- DB[, c("ID", 
                             "Language",  
                             "adjchange_abstheta", "adjchange_modultheta", "adjchange_reltheta",
                             "Age_at_T3")]
DB_clean_Language3 <- na.omit(DB_subset_Language3)

full_model_Language3 <- lm(Language ~ 
                            adjchange_abstheta + adjchange_modultheta + adjchange_reltheta +
                            Age_at_T3, data = DB_clean_Language3)

best_model_Language3 <- stepAIC(full_model_Language3, direction = "both")

summary(best_model_Language3)
vif(best_model_Language3)

model_std66 <- standardize_parameters(best_model_Language3)
print(model_std66)
