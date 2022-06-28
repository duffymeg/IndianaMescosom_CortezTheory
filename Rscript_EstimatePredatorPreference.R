# AUTHOR:   Michael Cortez mcortez@fsu.edu
# DATE:     March 2022
# PURPOSE:  Analyze data from predator preference experiments

# INPUTS:   Data from predation trials
#           Column 1 is host treatment:
#               control = no exposure to spores
#               exposed = exposed to spores but not visibly infected
#               infected = exposed to spores and visibly infected
#           Column 8 is day.block 
#               Blocks used because hosts are not visibly infected until day 10 or later
#               Block1 is measurements on days 1-7
#               Block 2 is measurements on days 10-13
#           Column 9 is percentage of dead hosts (out of 10)

# METHOD:   Step 1: Import data
#           Step 2: Two-way ANOVA
#           Step 3: Multiple comparison tests using Tukey's HSD 



library(dplyr)

#########################################################################
### STEP 1: Import data
Preference_data = read.csv("PredatorPreferenceData/predationtrialresultsclean_mhc.csv")

#########################################################################
### STEP 2: Two-way ANOVA

# Statistical model is percent.dead ~ day.block + treatment + day.block:treatment
# Fit the statistical model using lm
Preference_model <- lm(percent.dead ~ day.block + treatment + day.block : treatment, data = Preference_data)
# Create anova table
anova(Preference_model)

#########################################################################
### STEP 3: Multiple comparison tests using Tukey's HSD 

# Convert object from lm into aov object
Preference_aov <- aov(Preference_model)

# Perform Tukey HSD test
# Only "treatment" is used because day.block and interaction were not significant 
TukeyHSD(Preference_aov, which = 'treatment')

#########################################################################
### STEP 4: Compute increase in predation 

means = Preference_data %>% filter(treatment=='control' | treatment=='exposed' | treatment=='infected') %>% group_by(treatment) %>%
  summarise(mean=mean(percent.dead))

# Difference in attack rates between control and infected individuals is 
# the negative natural log of their ratio
attack_rate_difference = -log(means$mean[1]/means$mean[3])
attack_rate_difference



