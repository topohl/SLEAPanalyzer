#' @title Open Field Test (OFT) Data Processing Script
#' @description This script processes Open Field Test (OFT) data from multiple experimental batches.
#' It reads tracking data from CSV files, calibrates and analyzes the data, assigns zones, and generates plots.
#' The results are saved in structured output directories.
#' 
#' @author Tobias Pohl
#' 
#' @date 2025-01-30
#' 
#' @details
#' The script performs the following steps:
#' 
#' 1. Installs and loads required R packages.
#' 2. Sets the working directory and sources custom analysis functions.
#' 3. Defines experimental batches to process.
#' 4. Reads tracking data, calibrates measurements, assigns zones, and performs Open Field Test analysis.
#' 5. Extracts behavioral metrics and saves results in CSV format.
#' 6. Generates and saves density plots of movement patterns.

# Load required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  stringr, sp, imputeTS, ggplot2, ggmap, data.table, cowplot, corrplot, keras, tensorflow, zoo, dplyr, reticulate
)

# Set working directory (modify path as needed)
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")

# Load external functions required for analysis
source('R/DLCAnalyzer_Functions_final.R')

# Define batches to process
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")

# Loop through each batch for processing
for (batch in batches) {

  # Define input and output directories for the current batch
  inputDir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior", batch, "OFT/SLEAP/formatted")
  outputDir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior", batch, "OFT/SLEAP/output")
  
  # Create output directory if it does not exist
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  
  # Load animal ID and Code mappings
  animalIDCode <- read.table("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt", header = TRUE)
  
  # Create directory for saving plots
  plotDir <- file.path(outputDir, "plots")
  dir.create(plotDir, showWarnings = FALSE)
  
  # Retrieve list of tracking data files
  fileList <- list.files(path = inputDir, pattern = "*.csv")
  
  # Initialize list to store processed data frames
  dfList <- list()
  
  # Process each tracking data file in the batch
  for (file in fileList) {
    # Read in tracking data from CSV file
    inputFile <- file.path(inputDir, file)
    tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)
    
    # Extract file name without extension
    inputFileName <- sub(".csv$", "", basename(inputFile))
    
    # Calibrate data and assign zones
    tracking <- CalibrateTrackingData(tracking, method = "area", in.metric = 49 * 49, points = c("tl", "tr", "br", "bl"))
    tracking <- AddOFTZones(tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.4, points = c("tl", "tr", "br", "bl"))
    
    # Generate and visualize tracking plots
    PlotZones(tracking)
    PlotPointData(tracking, points = c("nose"))
    
    # Perform behavioral analysis
    tracking <- OFTAnalysis(tracking, points = "bodycentre", movement_cutoff = 5, integration_period = 5)
    
    # Extract unique identifier code from filename
    code <- str_extract(inputFileName, "[A-Za-z0-9]{4}")
    
    # Compute the frequency of rear events
    frequencyRear <- sum(freqRear)
    
    # Construct data frame with key behavioral metrics
    df <- data.frame(
      file = inputFileName,
      ID = animalIDCode[animalIDCode$Code == code, "ID"],
      Batch = batch,
      Code = code,
      centerTime = tracking$Report$bodycentre.center.total.time,
      centerDistance = tracking$Report$bodycentre.center.raw.distance,
      peripheryDistance = tracking$Report$bodycentre.periphery.raw.distance,
      centerTransitions = tracking$Report$bodycentre.center.transitions,
      rawSpeed = tracking$Report$bodycentre.raw.speed,
      rawSpeedCenter = tracking$Report$bodycentre.center.raw.speed,
      rawSpeedPeriphery = tracking$Report$bodycentre.periphery.raw.speed,
      stationary = tracking$Report$bodycentre.time.stationary,
      centerStationary = tracking$Report$bodycentre.center.time.stationary,
      peripheryStationary = tracking$Report$bodycentre.periphery.time.stationary,
      frequencyRear = frequencyRear,
      rawDistance = tracking$Report$bodycentre.raw.distance
    )
    
    # Save processed data to CSV file
    outputFile <- file.path(outputDir, paste0(inputFileName, "_output.csv"))
    write.csv(df, outputFile, row.names = FALSE)
    
    # Append current data frame to the list
    dfList[[length(dfList) + 1]] <- df
    
    # Generate density path plots
    plots <- PlotDensityPaths(tracking, points = c("bodycentre"))
    
    # Set the file name for the plot image
    plotFile <- file.path(plotDir, paste0(inputFileName, "_DensityPath.png"))
    
    # Save plot as an image
    ggsave(plotFile, plot = plots$bodycentre, width = 7, height = 4)
    
    # Print processing status
    message(sprintf("Processed file %s", file))
  }
  
  # Combine all batch data into a single data frame
  dfCombined <- do.call(rbind, dfList)
  
  # Save the combined results
  write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)
  
  # Print completion message for the batch
  message(sprintf("Processing complete for batch %s", batch))
}