# Load required libraries
library(sp)
library(imputeTS)
library(ggplot2)
library(ggmap)
library(data.table)
library(cowplot)
library(corrplot)
library(keras)
library(tensorflow)
library(zoo)
library(data.table)

# Function to process each file and generate plots
processFile <- function(file, input_dir, output_dir, plot_dir) {
  # Read in tracking data and get names of variables in the data frame
  input_file <- file.path(input_dir, file)
  Tracking <- ReadDLCDataFromCSV(file = input_file, fps = 30)
  
  # Extract the input file name without the extension
  input_file_name <- sub(".csv$", "", basename(input_file))
  
  # Replace NAs in the x and y columns of the nose data frame with the last known values
  Tracking$data$nose$x <- zoo::na.locf(Tracking$data$nose$x)
  Tracking$data$nose$y <- zoo::na.locf(Tracking$data$nose$y)
  
  # Calibrate the tracking data, add zones
  Tracking <- CalibrateTrackingData(Tracking, method = "area", in.metric = 49*49, points = c("tl","tr","br","bl"))
  Tracking <- AddOFTZones(Tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.2, points = c("tl","tr","br","bl"))
  Tracking$px.to.cm
  Tracking <- OFTAnalysis(Tracking, points = "bodycentre", movement_cutoff = 1, integration_period = 5)
  
  # Calculate frequency of zone visits for objL and nose
  FreqRear <- cumsum(GetDistances(Tracking, "spine1", "bodycentre") <= 1) > 0 &
    cumsum(GetDistances(Tracking, "bodycentre", "spine2") <= 1) > 0 &
    c(0, diff(GetDistances(Tracking, "spine1", "bodycentre") <= 1)) == 1 &
    c(0, diff(GetDistances(Tracking, "bodycentre", "spine2") <= 1)) == 1
  FrequencyRear <- sum(FreqRear)
  
  # Create a data frame with the time_Left, time_Right, time_nose_Left, and time_nose_Right values
  df <- data.frame(
    file = input_file_name,
    center_time = Tracking$Report$bodycentre.center.total.time,
    center_distance = Tracking$Report$bodycentre.center.raw.distance,
    periphery_distance = Tracking$Report$bodycentre.periphery.raw.distance,
    center_transitions = Tracking$Report$bodycentre.center.transitions,
    raw_speed = Tracking$Report$bodycentre.raw.speed,
    raw_speed_center = Tracking$Report$bodycentre.center.raw.speed,
    raw_speed_periphery = Tracking$Report$bodycentre.periphery.raw.speed,
    stationary = Tracking$Report$bodycentre.time.stationary,
    center_stationary = Tracking$Report$bodycentre.center.time.stationary,
    periphery_stationary = Tracking$Report$bodycentre.periphery.time.stationary,
    FrequencyRear = FrequencyRear
  )
  
  # Construct the output file name
  output_file <- paste0(input_file_name, "_output.csv")
  
  # Write the data frame to a csv file with the constructed file name
  write.csv(df, file.path(output_dir, output_file), row.names = FALSE)
  
  # Plot the density paths for bodycentre, nose, and tailbase
  plots <- PlotDensityPaths(Tracking, points = c("bodycentre"))
  
  # Set the file name for the plot image
  plot_file <- file.path(plot_dir, paste0(input_file_name, "_plot.png"))
  
  # Save the plot image with the file name
  ggsave(plot_file, plot = plots$bodycentre, width = 7, height = 6)
  
  # Plot the zone visits for bodycentre
  zone_visits_plot <- PlotZoneVisits(Tracking, point = c("bodycentre"), zones = c("center", "periphery"))
  
  # Set the file name for the zone visits plot image
  zone_visits_file <- file.path(plot_dir, paste0(input_file_name, "_zone_visits.png"))
  
  # Save the zone visits plot image with the file name
  ggsave(zone_visits_file, plot = zone_visits_plot, width = 7, height = 2)
  
  # Return the data frame and plots
  return(list(df = df, plots = plots, zone_visits_plot = zone_visits_plot))
}

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Set input and output directories
input_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B4/OFT/SLEAP/formatted"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B4/OFT/SLEAP/output"

# Create subdirectories within the output directory
plot_dir <- file.path(output_dir, "DensityPlot")
dir.create(plot_dir, showWarnings = FALSE)

# Get a list of CSV files in the input directory
file_list <- list.files(path = input_dir, pattern = "*.csv")

# Process each file and generate plots using lapply for efficiency
results_list <- lapply(file_list, processFile, input_dir = input_dir, output_dir = output_dir, plot_dir = plot_dir)

# Extract data frames and plots from the results list
df_list <- lapply(results_list, function(result) result$df)
plots_list <- lapply(results_list, function(result) result$plots)
zone_visits_list <- lapply(results_list, function(result) result$zone_visits_plot)

# Combine all data frames into a single data frame
df_combined <- do.call(rbind, df_list)

# Write the combined data frame to a csv file
write.csv(df_combined, file.path(output_dir, "combined_output.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")

