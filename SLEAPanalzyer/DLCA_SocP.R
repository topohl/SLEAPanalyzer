# install tensor flow and prepare for libraries and add new environment r-reticulate
install.packages("tensorflow")
library(reticulate)
path_to_python <- install_python()
virtualenv_create("r-reticulate", python = path_to_python)
library(tensorflow)
install_tensorflow(envname = "r-reticulate")
install.packages("keras")
library(keras)
install_keras(envname = "r-reticulate")
library(tensorflow)

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

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Set input and output directories
input_dir <- "C:/Users/topohl/Documents/GitHub/sleap/docs/notebooks"
output_dir <- "C:/Users/topohl/Documents/GitHub/sleap/docs/notebooks/output"

# Get a list of CSV files in the input directory
file_list <- list.files(path = input_dir, pattern = "*.csv")

# create new list df_list
df_list <- list()

# Loop through each file in the input directory
for (file in file_list) {
  
  # Read in tracking data and get names of variables in the data frame
  input_file <- file.path(input_dir, file)
  Tracking <- ReadDLCDataFromCSV(file = input_file, fps = 30)

  # Extract the input file name without the extension
  input_file_name <- sub(".csv$", "", basename(input_file))

  # Replace NAs in the x and y columns of the nose data frame with the last known values
  Tracking$data$nose$x <- zoo::na.locf(Tracking$data$nose$x)
  Tracking$data$nose$y <- zoo::na.locf(Tracking$data$nose$y)
  
  # Calibrate the tracking data, add zones, and plot the zones
  Tracking <- CalibrateTrackingData(Tracking, method = "area",in.metric = 44*24, points = c("tl","tr","br","bl"))
  Tracking <- AddOFTZones(Tracking, scale_center = 0.5,scale_periphery  = 0.8 ,scale_corners = 0.4, points = c("tl","tr","br","bl"))
  PlotZones(Tracking)
  PlotPointData(Tracking, points = c("nose"))

  # Calculate time_Left and time_Right direct
  total_time <- 10 * 60
  IsInZoneLeft <- GetDistances(Tracking, "socl", "nose") < 6
  IsInZoneRight <- GetDistances(Tracking, "socr", "nose") < 6
  ContactLeft <- sum(IsInZoneLeft) * total_time / length(IsInZoneLeft)
  ContactRight <- sum(IsInZoneRight) * total_time / length(IsInZoneRight)
  
  # Calculate time with nose between 6-10 cm of either socl or socr
  IsInProxLeft <- GetDistances(Tracking, "socl", "nose") >= 6 & GetDistances(Tracking, "socl", "nose") <= 10
  IsInProxRight <- GetDistances(Tracking, "socr", "nose") >= 6 & GetDistances(Tracking, "socr", "nose") <= 10
  ProxLeft <- sum(IsInProxLeft) * total_time / length(IsInProxLeft)
  ProxRight <- sum(IsInProxRight) * total_time / length(IsInProxRight)

  # Create a data frame with the time_Left, time_Right, time_nose_Left, and time_nose_Right values
  df <- data.frame(ContactLeft = ContactLeft, ContactRight = ContactRight, 
                 ProxLeft = ProxLeft, ProxRight = ProxRight)

  # Construct the output file name
  output_file <- paste0(input_file_name, "_output.csv")

  # Write the data frame to a csv file with the constructed file name
  write.csv(df, output_file, row.names = FALSE)

  # Add the current data frame to the list
  df_list[[length(df_list) + 1]] <- df

  # Print progress message
  cat(sprintf("Processed file %s\n", file))
}

# Combine all data frames into a single data frame
df_combined <- do.call(rbind, df_list)

# Write the combined data frame to a csv file
write.csv(df_combined, file.path(output_dir, "combined_output.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")




