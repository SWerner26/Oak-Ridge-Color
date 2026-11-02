# Load necessary libraries
library(dplyr)
library(lme4)
library(lmerTest) # Computes p-values for lmer models
library(ggplot2)

# 1. Load and merge the datasets
design_df <- read.csv("design.csv")
features_df <- read.csv("features.csv")

merged_df <- inner_join(design_df, features_df, 
                        by = c("sample_key", "experiment_id", "plant_id", "round_id"))

# 2. Calculate EXR, survival, and scale the time variable
model_df <- merged_df %>%
  mutate(
    EXR = 1.4 * proportional_red_mean - proportional_green_mean,
    # Standardize round_id to fix the LMM scaling warning
    round_id_scaled = as.numeric(scale(round_id))
  ) %>%
  group_by(plant_id) %>%
  mutate(
    # Check if all area_px values are NA to prevent the -Inf warning
    max_area_px = if(all(is.na(area_px))) NA_real_ else max(area_px, na.rm = TRUE),
    survived = as.factor(ifelse(max_area_px > 20000, 1, 0)),
    treatment = as.factor(treatment) 
  ) %>%
  ungroup()

# 3. Fit the Linear Mixed-Effects Model with the scaled variable
exr_model <- lmer(EXR ~ treatment * survived * round_id_scaled + (1 | plant_id), 
                  data = model_df,
                  na.action = na.exclude)

# 4. Output the model summaries
summary(exr_model)
anova(exr_model)

# 5. (Optional) Visualize the model's fitted values to verify the trajectories
model_df$fitted_EXR <- predict(exr_model)

ggplot(model_df, aes(x = round_id, y = fitted_EXR, color = treatment, linetype = survived)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    title = "LMM Fitted EXR Trajectories by Nickel Treatment and Survival",
    x = "Imaging Round (Time)",
    y = "Fitted Excess Red (EXR)",
    color = "Nickel (ppm)",
    linetype = "Survived"
  )

########################EXG#############################

# 1. Load and merge the datasets
design_df <- read.csv("design.csv")
features_df <- read.csv("features.csv")

merged_df <- inner_join(design_df, features_df, 
                        by = c("sample_key", "experiment_id", "plant_id", "round_id"))

# 2. Calculate EXR, survival, and scale the time variable
model_df <- merged_df %>%
  mutate(
    EXG = 2 * proportional_green_mean - proportional_red_mean - proportional_blue_mean,
    # Standardize round_id to fix the LMM scaling warning
    round_id_scaled = as.numeric(scale(round_id))
  ) %>%
  group_by(plant_id) %>%
  mutate(
    # Check if all area_px values are NA to prevent the -Inf warning
    max_area_px = if(all(is.na(area_px))) NA_real_ else max(area_px, na.rm = TRUE),
    survived = as.factor(ifelse(max_area_px > 20000, 1, 0)),
    treatment = as.factor(treatment) 
  ) %>%
  ungroup()

# 3. Fit the Linear Mixed-Effects Model with the scaled variable
exg_model <- lmer(EXG ~ treatment * survived * round_id_scaled + (1 | plant_id), 
                  data = model_df,
                  na.action = na.exclude)

# 4. Output the model summaries
summary(exg_model)
anova(exg_model)

# 5. (Optional) Visualize the model's fitted values to verify the trajectories
model_df$fitted_EXG <- predict(exg_model)

ggplot(model_df, aes(x = round_id, y = fitted_EXG, color = treatment, linetype = survived)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    title = "LMM Fitted EXG Trajectories by Nickel Treatment and Survival",
    x = "Imaging Round (Time)",
    y = "Fitted Excess Green (EXG)",
    color = "Nickel (ppm)",
    linetype = "Survived"
  )




# 1. CHECK THE ORIGINAL RESIDUALS
# Plot 1: Residuals vs Fitted (Checks for homoscedasticity / equal variance)
plot(exr_model, main = "Original Model: Residuals vs Fitted")

# Plot 2: Q-Q Plot (Checks for normality of residuals)
# If the dots curve wildly off the red line at the top/bottom, the outliers are skewing the model
qqnorm(resid(exr_model), main = "Original Model: Q-Q Plot")
qqline(resid(exr_model), col = "red")


model_df <- merged_df %>%
  mutate(
    EXR = 1.4 * proportional_red_mean - proportional_green_mean,
    round_id_scaled = as.numeric(scale(round_id))
  ) %>%
  group_by(plant_id) %>%
  mutate(
    max_area_px = if(all(is.na(area_px))) NA_real_ else max(area_px, na.rm = TRUE),
    survived = as.factor(ifelse(max_area_px > 20000, 1, 0)),
    treatment = as.factor(treatment) 
  ) %>%
  ungroup()

# 2. Calculate the shift value to ensure no negative numbers
min_exr <- min(model_df$EXR, na.rm = TRUE)
shift_value <- abs(min(0, min_exr)) + 0.01 

# 3. Apply the Log Transformation
model_df <- model_df %>%
  mutate(
    EXR_log = log(EXR + shift_value)
  )

# 4. Fit the Log-Transformed LMM
exr_model_log <- lmer(EXR_log ~ treatment * survived * round_id_scaled + (1 | plant_id), 
                      data = model_df,
                      na.action = na.exclude)

# 5. Check the Corrected Residuals and New Results
qqnorm(resid(exr_model_log), main = "Log-Transformed Model: Q-Q Plot")
qqline(resid(exr_model_log), col = "red")

anova(exr_model_log)