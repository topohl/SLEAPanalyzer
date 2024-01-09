
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
input_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B3/NOR/SLEAP/formatted"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B3/NOR/SLEAP/testoutput"
# Create a new folder for saving the plots
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# Get a list of CSV files in the input directory
file_list <- list.files(path = input_dir, pattern = "*.csv")

# create a new list df_list
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
  Tracking$data$bodycentre$x <- zoo::na.locf(Tracking$data$bodycentre$x)
  Tracking$data$bodycentre$y <- zoo::na.locf(Tracking$data$bodycentre$y)
  
  # Calibrate the tracking data, add zones, and plot the zones
  Tracking <- CalibrateTrackingData(Tracking, method = "area", in.metric = 49 * 49, points = c("tl", "tr", "br", "bl"))
  Tracking <- AddOFTZones(Tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.4, points = c("tl", "tr", "br", "bl"))
  PlotZones(Tracking)
  PlotPointData(Tracking, points = c("nose"))
  
  Tracking <- OFTAnalysis(Tracking, points = "bodycentre", movement_cutoff = 5, integration_period = 5)
  
  # Create a new data frame to store vector, distance, angle, and circle values
  df_tracking <- data.frame(
    vectorNoseBodycentre = cbind((Tracking$data$nose$x - Tracking$data$bodycentre$x), (Tracking$data$nose$y - Tracking$data$bodycentre$y)),
    vectorNoseObjL = cbind((Tracking$data$nose$x - Tracking$data$objL$x), (Tracking$data$nose$y - Tracking$data$objL$y)),
    vectorNoseObjR = cbind((Tracking$data$nose$x - Tracking$data$objR$x), (Tracking$data$nose$y - Tracking$data$objR$y)),
    angleL = acos(sum((Tracking$data$nose$x - Tracking$data$bodycentre$x) * (Tracking$data$objL$x - Tracking$data$bodycentre$x) + (Tracking$data$nose$y - Tracking$data$bodycentre$y) * (Tracking$data$objL$y - Tracking$data$bodycentre$y)) / (sqrt(sum((Tracking$data$nose$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$nose$y - Tracking$data$bodycentre$y)^2)) * sqrt(sum((Tracking$data$objL$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$objL$y - Tracking$data$bodycentre$y)^2)))),
    angleR = acos(sum((Tracking$data$nose$x - Tracking$data$bodycentre$x) * (Tracking$data$objR$x - Tracking$data$bodycentre$x) + (Tracking$data$nose$y - Tracking$data$bodycentre$y) * (Tracking$data$objR$y - Tracking$data$bodycentre$y)) / (sqrt(sum((Tracking$data$nose$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$nose$y - Tracking$data$bodycentre$y)^2)) * sqrt(sum((Tracking$data$objR$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$objR$y - Tracking$data$bodycentre$y)^2)))),
    circleL = sqrt((Tracking$data$objL$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$objL$y - Tracking$data$bodycentre$y)^2),
    circleR = sqrt((Tracking$data$objR$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$objR$y - Tracking$data$bodycentre$y)^2)
  )
  
  # Create new variables in the data frame to indicate if the circle is within 4cm
  df_tracking$circleL <- ifelse(df_tracking$circleL <= 4, 1, 0)
  df_tracking$circleR <- ifelse(df_tracking$circleR <= 4, 1, 0)
  
  # Initialize vectors to store angles
  anglesDegreesObjL <- numeric(length = nrow(df_tracking))
  anglesDegreesObjR <- numeric(length = nrow(df_tracking))
  
  # Loop through each row for vectorNoseObjL
for (i in 1:nrow(df_tracking)) {
  # Extract vectors a and c from vectorNoseObjL and vector matrices
  a <- c(df_tracking$vectorNoseObjL.1[i], df_tracking$vectorNoseObjL.2[i])
  c <- c(df_tracking$vectorNoseBodycentre.1[i], df_tracking$vectorNoseBodycentre.2[i])

  # Calculate dot product
  dotProductAC <- sum(a * c)

  # Calculate magnitudes
  magnitudeA <- sqrt(sum(a^2))
  magnitudeC <- sqrt(sum(c^2))

  # Calculate the angle in radians
  angleRadiansAC <- acos(dotProductAC / (magnitudeA * magnitudeC))

  # Convert angle to degrees and store in the vector
  anglesDegreesObjL[i] <- angleRadiansAC * 180 / pi
}
  
  # Loop through each row for vectorNoseObjR
for (i in 1:nrow(df_tracking)) {
    # Extract vectors a and c from vectorNoseObjR and vector matrices
    a <- c(df_tracking$vectorNoseObjR.1[i], df_tracking$vectorNoseObjR.2[i])
    c <- c(df_tracking$vectorNoseBodycentre.1[i], df_tracking$vectorNoseBodycentre.2[i])

    # Calculate dot product
    dotProductAC <- sum(a * c)

    # Calculate magnitudes
    magnitudeA <- sqrt(sum(a^2))
    magnitudeC <- sqrt(sum(c^2))

    # Calculate the angle in radians
    angleRadiansAC <- acos(dotProductAC / (magnitudeA * magnitudeC))

    # Convert angle to degrees and store in the vector
    anglesDegreesObjR[i] <- angleRadiansAC * 180 / pi
}
   
  # Calculate time_Left and time_Right only if angle L or R is between 90 and 270 degrees
  total_time <- 10 * 60
  IsInZoneLeft <- GetDistances(Tracking, "objL", "nose") <= 4 & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(anglesDegreesObjL) >= 90 & abs(anglesDegreesObjL) <= 270
  IsInZoneRight <- GetDistances(Tracking, "objR", "nose") <= 4 & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(anglesDegreesObjR) >= 90 & abs(anglesDegreesObjR) <= 270
  ContactLeft <- sum(IsInZoneLeft) * total_time / length(IsInZoneLeft)
  ContactRight <- sum(IsInZoneRight) * total_time / length(IsInZoneRight)
  
  # Calculate time with nose between 6-10 cm of either objL or objR
  IsInProxLeft <- GetDistances(Tracking, "objL", "nose") > 4 & GetDistances(Tracking, "objL", "nose") <= 8
  IsInProxRight <- GetDistances(Tracking, "objR", "nose") > 4 & GetDistances(Tracking, "objR", "nose") <= 8
  ProxLeft <- sum(IsInProxLeft) * total_time / length(IsInProxLeft)
  ProxRight <- sum(IsInProxRight) * total_time / length(IsInProxRight)
  
  # Calculate frequency of zone visits for objL and nose
  IsInZoneVisitL <- cumsum(GetDistances(Tracking, "objL", "nose") <= 4) > 0 & c(0, diff(GetDistances(Tracking, "objL", "nose") <= 4)) == 1
  FrequencyL <- sum(IsInZoneVisitL)
  IsInZoneVisitR <- cumsum(GetDistances(Tracking, "objR", "nose") <= 4) > 0 & c(0, diff(GetDistances(Tracking, "objR", "nose") <= 4)) == 1
  FrequencyR <- sum(IsInZoneVisitR)
  
  # Calculate frequency of rearings
  FreqRear <- cumsum(GetDistances(Tracking, "spine1", "bodycentre") <= 1) > 0 &
    cumsum(GetDistances(Tracking, "bodycentre", "spine2") <= 1) > 0 &
    c(0, diff(GetDistances(Tracking, "spine1", "bodycentre") <= 1)) == 1 &
    c(0, diff(GetDistances(Tracking, "bodycentre", "spine2") <= 1)) == 1
  FrequencyRear <- sum(FreqRear)
  
  # Create a data frame with the time_Left, time_Right, time_nose_Left, and time_nose_Right values
  df <- data.frame(
    file = input_file_name,
    ContactLeft = ContactLeft,
    ContactRight = ContactRight,
    ProxLeft = ProxLeft,
    ProxRight = ProxRight,
    FrequencyL = FrequencyL,
    FrequencyR = FrequencyR,
    distance = Tracking$Report$bodycentre.raw.distance,
    stationary = Tracking$Report$bodycentre.time.stationary,
    speed_moving = Tracking$Report$bodycentre.speed.moving,
    speed_raw = Tracking$Report$bodycentre.raw.speed,
    FrquencyRear = FrequencyRear
  )
  
  # Construct the output file name
  output_file <- paste0(input_file_name, "_output.csv")
  
  # Write the data frame to a csv file with the constructed file name
  write.csv(df, output_file, row.names = FALSE)
  
  # Add the current data frame to the list
  df_list[[length(df_list) + 1]] <- df
  
  # Plot the density paths for bodycentre
  plots <- PlotDensityPaths(Tracking, points = c("bodycentre"))
  
  # Set the file name for the plot image
  plot_file <- file.path(plot_dir, paste0(input_file_name, "_DensityPath.png"))
  
    # Save the plot image with the file name
    ggsave(plot_file, plot = plots$bodycentre, width = 7, height = 6)
    
    # Print progress message
    cat(sprintf("Processed file %s\n", file))
}

# Combine all data frames into a single data frame
df_combined <- do.call(rbind, df_list)

# Write the combined data frame to a csv file
write.csv(df_combined, file.path(output_dir, "combined_output.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")