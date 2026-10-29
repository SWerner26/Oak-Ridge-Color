library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(tidyr)
library(lme4)
library(lmerTest)


##########################LOAD DATA AND MERGE DATASETS BY SAMPLE_KEY###################################
#Load both datasets
df_features <- read_csv("features.csv")
df_design <- read_csv("design.csv")

#Subset only takes sample_key, genotype and treatment from the design file
df_design_subset <- df_design %>%
  select(sample_key, genotype, treatment) %>%
  distinct()

#Merges feature and design_subset files sorted by sample_key
df_merged <- df_features %>% 
  inner_join(df_design_subset, by = "sample_key")

#########################PLANT COLOR CHARACTERISTICS OVER TIME#########################################
#Create a df for calculating and storing daily average color values and area
df_trends <- df_merged %>%
  mutate(day = as.Date(datetime)) %>% #Takes datetime column and changes it to just the day
  drop_na(red_mean, green_mean, blue_mean, area_px) %>% #Drops rows where red, green or blue are missing
  group_by(day) %>% #Groups the df by the day
  summarize(
    avg_red = mean(red_mean), #Calculates the mean red pixel intensity for the day, removing missing values
    avg_green = mean(green_mean), #Does the same as above for green
    avg_blue = mean(blue_mean), #Does the same as above for blue
    avg_ngrdi = mean(ngrdi_mean), #Does the same as above for ngrdi
    avg_area = mean(area_px), #Calculates the mean area for the day, removing missing values
    .groups = "drop"
  )

#Reshapes df so ggplot can graph multiple lines
df_color_long <- df_trends %>%
  select(day, avg_red, avg_green, avg_blue) %>% 
  pivot_longer(
    cols = c(avg_red, avg_green, avg_blue),
    names_to = "channel", #Names the column for each color "channel"
    values_to = "intensity" #Names the column with the pixel intensity values as "intensity"
  )

#Plots red, green and blue means 
p1 <- ggplot(df_color_long, aes(x = day, y = intensity, color = channel)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("avg_red" = "#E41A1C", "avg_green" = "#4DAF4A", "avg_blue" = "#377EB8"),
    labels = c("avg_red" = "Red", "avg_green" = "Green", "avg_blue" = "Blue")
  ) +
  theme_minimal() +
  labs(
    title = "Average Plant RGB Color Intensity Over Time",
    x = NULL, # Remove x-axis label here to avoid redundancy with the bottom plot
    y = "Mean Pixel Intensity",
    color = "Color Channel"
  )+
  theme(axis.text.x = element_blank()) # Hide top x-axis text to clean up the stack

#Plots the NGRDI separate since it is on a different scale
ggplot(df_trends, aes(x = day, y = avg_ngrdi)) +
  geom_line(color = "#984EA3", linewidth = 1) + 
  geom_point(color = "#984EA3", size = 2) +
  theme_minimal() +
  labs(
    title = "Average Normalized Green-Red Difference Index (NGRDI) Over Time",
    x = "Date",
    y = "Mean NGRDI"
  )

#Create the Area Plot (p2)
p2 <- ggplot(df_trends, aes(x = day, y = avg_area)) +
  geom_line(color = "#333333", linewidth = 1) +
  geom_point(color = "#333333", size = 2) +
  theme_minimal() +
  labs(
    subtitle = "Average Plant Area",
    x = "Date",
    y = "Mean Area (Pixels)"
  )

#Stack the plots using patchwork, The '/' operator stacks p1 on top of p2
combined_plot <- p1 / p2

#Display the combined plot
print(combined_plot)

#Display summary statistics for each average value
df_trends %>% 
  select(avg_area, avg_red, avg_green, avg_blue, avg_ngrdi, day) %>% 
  summary()

#####################PLANT COLOR ASSOCIATION WITH GROWTH###############################################
# Calculate Correlation Coefficients
# We isolate the area and color columns, drop missing values, and generate a correlation matrix
cor_matrix <- df_merged %>%
  drop_na(area_px, red_mean, green_mean, blue_mean, ngrdi_mean) %>%
  select(area_px, red_mean, green_mean, blue_mean, ngrdi_mean) %>%
  cor()

# Print the specific correlations between color metrics and plant area
print("Correlation with Plant Area (area_px):")
print(cor_matrix[, "area_px"])

#Pivot the df to graph
df_scatter <- df_merged %>%
  drop_na(area_px, red_mean) %>% #Removes missing values
  select(area_px, red_mean, green_mean, blue_mean, ngrdi_mean) %>%
  pivot_longer(
    cols = c(red_mean, green_mean, blue_mean, ngrdi_mean),
    names_to = "color_metric", #What color the value belongs to
    values_to = "metric_value" #Color intensity value
  ) %>%
  # Clean up the metric names for nicer plot labels
  mutate(color_metric = factor(color_metric, 
                               levels = c("red_mean", "green_mean", "blue_mean", "ngrdi_mean"),
                               labels = c("Red Mean", "Green Mean", "Blue Mean", "NGRDI Mean")))

##This graph shows the correlation between plant color and growth. 
#       The trend line is based on correlation
#       The colored dots are a scatter plot of the individual plants
#Create Scatter Plots with Trend Lines
ggplot(df_scatter, aes(x = metric_value, y = area_px, color = color_metric)) +
  # Using alpha = 0.2 makes overlapping points semi-transparent so dense clusters are visible
  geom_point(alpha = 0.2, size = 1.5) + 
  # Add a standard linear regression (least squares) line to visualize the trend direction
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) + 
  # Create a separate panel for each color metric and allow the x-axis to scale independently
  facet_wrap(~ color_metric, scales = "free_x") +
  # Map specific colors to match the metrics (purple for NGRDI)
  scale_color_manual(values = c("Red Mean" = "#E41A1C", 
                                "Green Mean" = "#4DAF4A", 
                                "Blue Mean" = "#377EB8", 
                                "NGRDI Mean" = "#984EA3")) +
  theme_minimal() +
  labs(
    title = "Association Between Plant Color Metrics and Growth",
    subtitle = "Plant Area (Pixels) vs. Color Averages",
    x = "Color Metric Value",
    y = "Plant Area (Pixels)"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold")
  )

######################MEAN COLOR AND AREA BY TREATMENT AND GENOTYPE###############################
#Determines how many different genotypes are in the subset file
unique_genotypes <- unique(df_design_subset$genotype)

#Loops through each genotype and calculates mean ngrdi and green then graphs them
for (sub in unique_genotypes) {
  df_summary_long <- df_merged %>%
    filter(!is.na(area_px), genotype == sub) %>% #Takes out rows with missing data, and only looks at rows with the genotype currently being used for the loop
    mutate(date_parsed = as.Date(date, format = "%m/%d/%Y")) %>% #Changes the date to an object ggplot will recognize
    group_by(date_parsed, treatment) %>% #Groups the data by date and treatment
    summarize(
      `Red Mean` = mean(red_mean, na.rm = TRUE), #Takes the mean red_mean for each day/treatment, removing missing values
      `Green Mean` = mean(green_mean, na.rm = TRUE),
      `Blue Mean` = mean(blue_mean, na.rm = TRUE),
      `Mean Area` = mean(area_px, na.rm = TRUE),
      .groups = "drop" #Ungroups them after math is done so no grouping errors occur later
    ) %>% 
    pivot_longer(cols = c(`Red Mean`, `Green Mean`, `Blue Mean`, `Mean Area`),
                 names_to = "metric", values_to = "value")
  
  #Puts similar genotypes next to each other for comparison 
  combined_plot <- ggplot(df_summary_long, aes(x = date_parsed, y = value, color = as.factor(treatment))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_y") + # This replaces the need for patchwork!
    theme_minimal() +
    labs(
      title = "Color Means & Area Over Time by Nickel Exposure Level",
      subtitle = paste("Genotype:", sub),
      x = "Date",
      y = "Value",
      color = "Nickel Level"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(combined_plot)
}

###########################SUMMARY STATISTICS#################################
summaryStats <- df_trends[c("avg_red", "avg_green", "avg_blue", "avg_area")]
summary(summaryStats)

range(df_trends$day)

unique_plant_id <- unique(df_design$plant_id)
length(unique_plant_id)
#######################DEATH HAVE A CERTAIN COLOR?#########################
##Shows distribution of area from smallest to largest 
# 1. Calculate max area and rank the plants from smallest to largest
plant_max_area <- df_merged %>%
  group_by(plant_id) %>%
  summarize(max_area = max(area_px, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_area) %>%
  # Add a sequential rank column (1 to 360) based on their sorted size
  mutate(rank = row_number())

# 2. Plot the ranked distribution
ggplot(plant_max_area, aes(x = rank, y = max_area)) +
  geom_point(color = "#377EB8", alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "Distribution of Maximum Plant Areas",
    subtitle = "Visualizing the growth gap across the cohort",
    x = "Plant Rank (1 = Smallest, 360 = Largest)",
    y = "Maximum Area Reached (Pixels)"
  )

# 1. Identify plants that plateaued/died early
# Group by plant_id to find the maximum area each individual plant achieved.
dead_plants_summary <- df_merged %>%
  group_by(plant_id) %>%
  summarize(max_area = max(area_px, na.rm = TRUE), .groups = "drop")

# Set the ceiling threshold to isolate the bottom quartile of growth
threshold <- 20000 
dead_plant_ids <- dead_plants_summary %>%
  filter(max_area < threshold) %>%
  pull(plant_id)

# 2. Filter the main dataset for only these 93 stunted plants
df_dead_plants <- df_merged %>%
  filter(plant_id %in% dead_plant_ids)

# 3. Calculate daily color averages for just this isolated cohort
df_dead_trends <- df_dead_plants %>%
  mutate(day = as.Date(datetime)) %>%
  drop_na(red_mean, green_mean, blue_mean) %>%
  group_by(day) %>%
  summarize(
    avg_red = mean(red_mean),
    avg_green = mean(green_mean),
    avg_blue = mean(blue_mean),
    avg_ngrdi = mean(ngrdi_mean),
    .groups = "drop"
  )

# 4. Prepare data for the RGB plot
df_dead_long <- df_dead_trends %>%
  select(day, avg_red, avg_green, avg_blue) %>%
  pivot_longer(
    cols = c(avg_red, avg_green, avg_blue),
    names_to = "channel",
    values_to = "intensity"
  )

# 5. Visualize the color trends for the plateaued plants
pDead <- ggplot(df_dead_long, aes(x = day, y = intensity, color = channel)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("avg_red" = "#E41A1C", "avg_green" = "#4DAF4A", "avg_blue" = "#377EB8"),
    labels = c("avg_red" = "Red", "avg_green" = "Green", "avg_blue" = "Blue")
  ) +
  theme_minimal() +
  labs(
    title = "RGB Color Trends for Plants that Plateaued Early",
    subtitle = paste("Cohort size:", length(dead_plant_ids), "plants"),
    x = "Date",
    y = "Mean Pixel Intensity",
    color = "Color Channel"
  )

# 1. Isolate just the unique plants and their treatments
unique_dead_plants <- df_dead_plants %>%
  distinct(plant_id, treatment)

# 2. Count the actual number of dead plants per treatment
deathTreatment <- table(unique_dead_plants$treatment)
print(deathTreatment)

combined_plot <- p1 / pDead
print(combined_plot)

##Does the trend line for each of the colors change if dead plants are removed? 
# 1. Identify plants that survived (max area >= 20,000 pixels)
alive_plant_ids <- df_merged %>%
  group_by(plant_id) %>%
  summarize(max_area = max(area_px, na.rm = TRUE), .groups = "drop") %>%
  filter(max_area >= 20000) %>%
  pull(plant_id)

# 2. Filter the main dataset to keep only the alive plants
df_alive_plants <- df_merged %>%
  filter(plant_id %in% alive_plant_ids)

# 3. Calculate daily color averages for just this surviving cohort
df_alive_trends <- df_alive_plants %>%
  mutate(day = as.Date(datetime)) %>%
  drop_na(red_mean, green_mean, blue_mean) %>%
  group_by(day) %>%
  summarize(
    avg_red = mean(red_mean),
    avg_green = mean(green_mean),
    avg_blue = mean(blue_mean),
    .groups = "drop"
  )

# 4. Pivot the data to a long format for plotting
df_alive_long <- df_alive_trends %>%
  select(day, avg_red, avg_green, avg_blue) %>%
  pivot_longer(
    cols = c(avg_red, avg_green, avg_blue),
    names_to = "channel",
    values_to = "intensity"
  )

# 5. Visualize the color trends for the alive plants
pAlive <- ggplot(df_alive_long, aes(x = day, y = intensity, color = channel)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("avg_red" = "#E41A1C", "avg_green" = "#4DAF4A", "avg_blue" = "#377EB8"),
    labels = c("avg_red" = "Red", "avg_green" = "Green", "avg_blue" = "Blue")
  ) +
  theme_minimal() +
  labs(
    title = "Average Plant RGB Color Intensity Over Time",
    subtitle = paste("Alive Cohort Only (", length(alive_plant_ids), " plants)", sep=""),
    x = "Date",
    y = "Mean Pixel Intensity",
    color = "Color Channel"
  )

combined_plot <- pAlive / pDead
print(combined_plot)

combined_plot <- p1 / pAlive
print(combined_plot)

# 1. Determine the status of each plant based on the 20,000 pixel threshold
plant_status <- df_merged %>%
  group_by(plant_id) %>%
  summarize(max_area = max(area_px, na.rm = TRUE), .groups = "drop") %>%
  mutate(status = ifelse(max_area >= 20000, "Alive", "Dead"))

# 2. Attach this status back to the main dataset
df_classified <- df_merged %>%
  inner_join(plant_status %>% select(plant_id, status), by = "plant_id")

# 3. Generate summary statistics for ALL plants
print("--- Summary Statistics: WHOLE DATASET ---")
df_classified %>% 
  select(red_mean, green_mean, blue_mean) %>% 
  summary() %>% 
  print()

# 4. Generate summary statistics for ALIVE plants only
print("--- Summary Statistics: ALIVE PLANTS ---")
df_classified %>% 
  filter(status == "Alive") %>% 
  select(red_mean, green_mean, blue_mean) %>% 
  summary() %>% 
  print()

# 5. Generate summary statistics for DEAD plants only
print("--- Summary Statistics: DEAD PLANTS ---")
df_classified %>% 
  filter(status == "Dead") %>% 
  select(red_mean, green_mean, blue_mean) %>% 
  summary() %>% 
  print()

# 1. Prepare the three daily trend datasets with a 'cohort' label
df_all_trends <- df_classified %>%
  mutate(day = as.Date(datetime)) %>%
  filter(!is.na(red_mean) & !is.na(green_mean) & !is.na(blue_mean)) %>%
  group_by(day) %>%
  summarize(
    avg_red = mean(red_mean, na.rm = TRUE),
    avg_green = mean(green_mean, na.rm = TRUE),
    avg_blue = mean(blue_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cohort = "All Plants")

df_alive_trends <- df_classified %>%
  filter(status == "Alive") %>%
  mutate(day = as.Date(datetime)) %>%
  filter(!is.na(red_mean) & !is.na(green_mean) & !is.na(blue_mean)) %>%
  group_by(day) %>%
  summarize(
    avg_red = mean(red_mean, na.rm = TRUE),
    avg_green = mean(green_mean, na.rm = TRUE),
    avg_blue = mean(blue_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cohort = "Alive")

df_dead_trends <- df_classified %>%
  filter(status == "Dead") %>%
  mutate(day = as.Date(datetime)) %>%
  filter(!is.na(red_mean) & !is.na(green_mean) & !is.na(blue_mean)) %>%
  group_by(day) %>%
  summarize(
    avg_red = mean(red_mean, na.rm = TRUE),
    avg_green = mean(green_mean, na.rm = TRUE),
    avg_blue = mean(blue_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cohort = "Dead")

# 2. Combine and reshape the data for Plot 1 (Alive vs. Dead)
df_compare_ad <- bind_rows(df_alive_trends, df_dead_trends) %>%
  pivot_longer(
    cols = c(avg_red, avg_green, avg_blue),
    names_to = "channel",
    values_to = "intensity"
  )

# 3. Create Plot 1: Alive vs. Dead
plot_alive_vs_dead <- ggplot(df_compare_ad, aes(x = day, y = intensity, color = channel, linetype = cohort)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("avg_red" = "#E41A1C", "avg_green" = "#4DAF4A", "avg_blue" = "#377EB8"),
    labels = c("avg_red" = "Red", "avg_green" = "Green", "avg_blue" = "Blue")
  ) +
  scale_linetype_manual(values = c("Alive" = "solid", "Dead" = "dotted")) +
  theme_minimal() +
  labs(
    title = "RGB Trends: Alive vs. Dead Cohorts",
    x = "Date",
    y = "Mean Pixel Intensity",
    color = "Color Channel",
    linetype = "Plant Status"
  )

# 4. Combine and reshape the data for Plot 2 (Alive vs. All Plants)
df_compare_all <- bind_rows(df_alive_trends, df_all_trends) %>%
  pivot_longer(
    cols = c(avg_red, avg_green, avg_blue),
    names_to = "channel",
    values_to = "intensity"
  )

# 5. Create Plot 2: Alive vs. All Plants
plot_alive_vs_all <- ggplot(df_compare_all, aes(x = day, y = intensity, color = channel, linetype = cohort)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("avg_red" = "#E41A1C", "avg_green" = "#4DAF4A", "avg_blue" = "#377EB8"),
    labels = c("avg_red" = "Red", "avg_green" = "Green", "avg_blue" = "Blue")
  ) +
  # Using a dotted line for the "All Plants" to distinguish it from the dashed "Dead" lines
  scale_linetype_manual(values = c("Alive" = "solid", "All Plants" = "dotted")) +
  theme_minimal() +
  labs(
    title = "RGB Trends: Alive vs. All Plants",
    x = "Date",
    y = "Mean Pixel Intensity",
    color = "Color Channel",
    linetype = "Cohort"
  )

# Display the plots (you can also use patchwork to stack them: plot_alive_vs_dead / plot_alive_vs_all)
print(plot_alive_vs_dead)
print(plot_alive_vs_all)

#########################EXG, EXR, VARI#################################
# Add EXG, EXR, and VARI to your merged dataframe
df_indices <- df_merged %>%
  mutate(day = as.Date(datetime)) %>%
  # Filter out rows with missing RGB data to prevent division-by-zero errors
  drop_na(red_mean, green_mean, blue_mean) %>%
  mutate(
    # 1. Calculate the sum of the RGB channels
    rgb_sum = red_mean + green_mean + blue_mean,
    
    # 2. Normalize the RGB values (creating coordinates r, g, b)
    r = red_mean / rgb_sum,
    g = green_mean / rgb_sum,
    b = blue_mean / rgb_sum,
    
    # 3. Calculate Excess Green (EXG) and Excess Red (EXR) using the normalized coordinates
    EXG = 2 * g - r - b,
    EXR = 1.4 * r - g,
    
    # 4. Calculate VARI (typically calculated using the raw/unnormalized channel means)
    VARI = (green_mean - red_mean) / (green_mean + red_mean - blue_mean)
  )

# Preview the new columns
head(df_indices %>% select(plant_id, day, EXG, EXR, VARI))

# 1. Calculate the daily averages for the three new indices
df_index_trends <- df_indices %>%
  group_by(day, treatment) %>%
  summarize(
    avg_EXG = mean(EXG, na.rm = TRUE),
    avg_EXR = mean(EXR, na.rm = TRUE),
    avg_VARI = mean(VARI, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Reshape the data from wide to long format for ggplot
df_index_long <- df_index_trends %>%
  pivot_longer(
    cols = c(avg_EXG, avg_EXR, avg_VARI),
    names_to = "index_name",
    values_to = "index_value"
  ) %>%
  # Rename the factors for cleaner facet labels in the plot
  mutate(index_name = factor(index_name, 
                             levels = c("avg_EXG", "avg_EXR", "avg_VARI"),
                             labels = c("Excess Green (EXG)", "Excess Red (EXR)", "VARI")))

# 3. Plot the indices using facet_wrap to accommodate their different scales
ggplot(df_index_long, aes(x = day, y = index_value, color = as.factor(treatment))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # Create a separate, vertically stacked panel for each index and let the y-axis adjust automatically
  facet_wrap(~ index_name, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    title = "Average Vegetation Indices Over Time",
    x = "Date",
    y = "Mean Index Value"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold")
  )

#########################EXG, EXR, VARI by Status #################################
# Add EXG, EXR, and VARI to your classified dataframe
df_indices <- df_classified %>%
  mutate(day = as.Date(datetime)) %>%
  # Filter out rows with missing RGB data to prevent division-by-zero errors
  drop_na(red_mean, green_mean, blue_mean) %>%
  mutate(
    # 1. Calculate the sum of the RGB channels
    rgb_sum = red_mean + green_mean + blue_mean,
    
    # 2. Normalize the RGB values (creating coordinates r, g, b)
    r = red_mean / rgb_sum,
    g = green_mean / rgb_sum,
    b = blue_mean / rgb_sum,
    
    # 3. Calculate Excess Green (EXG) and Excess Red (EXR) using the normalized coordinates
    EXG = 2 * g - r - b,
    EXR = 1.4 * r - g,
    
    # 4. Calculate VARI (typically calculated using the raw/unnormalized channel means)
    VARI = (green_mean - red_mean) / (green_mean + red_mean - blue_mean)
  )

# Preview the new columns
head(df_indices %>% select(plant_id, day, status, EXG, EXR, VARI))

# 1. Calculate the daily averages for the three new indices
df_index_trends <- df_indices %>%
  # CHANGED: Group by status instead of treatment
  group_by(day, status) %>%
  summarize(
    avg_EXG = mean(EXG, na.rm = TRUE),
    avg_EXR = mean(EXR, na.rm = TRUE),
    avg_VARI = mean(VARI, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Reshape the data from wide to long format for ggplot
df_index_long <- df_index_trends %>%
  pivot_longer(
    cols = c(avg_EXG, avg_EXR, avg_VARI),
    names_to = "index_name",
    values_to = "index_value"
  ) %>%
  # Rename the factors for cleaner facet labels in the plot
  mutate(index_name = factor(index_name, 
                             levels = c("avg_EXG", "avg_EXR", "avg_VARI"),
                             labels = c("Excess Green (EXG)", "Excess Red (EXR)", "VARI")))

# 3. Plot the indices using facet_wrap to accommodate their different scales
# CHANGED: Map color to 'status'
ggplot(df_index_long, aes(x = day, y = index_value, color = status)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # Create a separate, vertically stacked panel for each index and let the y-axis adjust automatically
  facet_wrap(~ index_name, scales = "free_y", ncol = 1) +
  # ADDED: Manually assign green to Alive and red to Dead for clarity
  scale_color_manual(values = c("Alive" = "#4DAF4A", "Dead" = "#E41A1C")) +
  theme_minimal() +
  labs(
    title = "Average Vegetation Indices Over Time by Plant Status",
    x = "Date",
    y = "Mean Index Value",
    color = "Plant Status" # Updated legend title
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold")
  )

####################Modeling###################################

#Logistic Regression (Predicting Survival via Color)
# 1. Extract the full suite of color indices for the first 5 days

early_color_predictors <- df_indices %>%
  filter(as.numeric(day - min(day, na.rm = TRUE)) <= 5) %>%
  group_by(plant_id, status) %>%
  summarize(
    early_EXG = mean(EXG, na.rm = TRUE),
    early_EXR = mean(EXR, na.rm = TRUE),
    early_VARI = mean(VARI, na.rm = TRUE),
    early_NGRDI = mean(ngrdi_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(is_alive = ifelse(status == "Alive", 1, 0))

# 2. Fit the model using ONLY color variables
color_survival_glm <- glm(is_alive ~ early_EXG + early_EXR + early_VARI + early_NGRDI, 
                          data = early_color_predictors, 
                          family = "binomial")

print("--- Logistic Regression: Predicting Survival from Early Color ---")
summary(color_survival_glm)

# 1. Test VARI purely on its own
vari_only_glm <- glm(is_alive ~ early_VARI, 
                     data = early_color_predictors, 
                     family = "binomial")

print("--- Logistic Regression: Survival Predicted by Early VARI ---")
summary(vari_only_glm)


# 2. Test NGRDI purely on its own
ngrdi_only_glm <- glm(is_alive ~ early_NGRDI, 
                      data = early_color_predictors, 
                      family = "binomial")

print("--- Logistic Regression: Survival Predicted by Early NGRDI ---")
summary(ngrdi_only_glm)

# 2. Test EXG purely on its own
exg_only_glm <- glm(is_alive ~ early_EXG, 
                      data = early_color_predictors, 
                      family = "binomial")

print("--- Logistic Regression: Survival Predicted by Early EXG ---")
summary(exg_only_glm)

# 2. Test EXR purely on its own
exr_only_glm <- glm(is_alive ~ early_EXR, 
                    data = early_color_predictors, 
                    family = "binomial")

print("--- Logistic Regression: Survival Predicted by Early EXR ---")
summary(exr_only_glm)


##Linear Mixed-Effects Model (Modeling Color Over Time)
# 1. Prepare the dataset with a numeric day and factorized status
df_color_lmm <- df_indices %>%
  mutate(
    day_numeric = as.numeric(day - min(day, na.rm = TRUE)),
    status = as.factor(status) # Treat Alive/Dead as categorical groups
  )

# 2. Fit the LMM predicting Excess Green (EXG) over time based on survival status
# (1 | plant_id) still accounts for natural plant-to-plant baseline color differences
exg_lmm <- lmer(EXG ~ day_numeric * status + (1 | plant_id), 
                data = df_color_lmm)

print("--- Linear Mixed Model: EXG (Greenness) Trajectory ---")
summary(exg_lmm)
anova(exg_lmm)


##Linear Mixed-Effects Model (Modeling Color Over Time) for NGRDI
# 1. Prepare the dataset with a numeric day and factorized status
df_color_lmm <- df_indices %>%
  mutate(
    day_numeric = as.numeric(day - min(day, na.rm = TRUE)),
    status = as.factor(status) # Treat Alive/Dead as categorical groups
  )

# 2. Fit the LMM predicting Excess Green (EXG) over time based on survival status
# (1 | plant_id) still accounts for natural plant-to-plant baseline color differences
NGRDI_lmm <- lmer(ngrdi_mean ~ day_numeric * status + (1 | plant_id), 
                data = df_color_lmm)

print("--- Linear Mixed Model: NGRDI Trajectory ---")
summary(NGRDI_lmm)
anova(NGRDI_lmm)