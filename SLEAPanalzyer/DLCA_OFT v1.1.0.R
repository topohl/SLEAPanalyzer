
# This script processes Open Field Test (OFT) data for multiple batches of experiments.
# It installs and loads required packages, sets the working directory, and sources necessary functions.
# The script reads tracking data from CSV files, calibrates the data, adds zones, performs analysis, 
# and generates plots. The results are saved to output directories.

# Install and load required packages
requiredPackages <- c("stringr", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "keras", "tensorflow", "zoo", "dplyr", "reticulate")

# Install missing packages and load all required packages
lapply(requiredPackages, function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    tryCatch({
      install.packages(package)
    }, error = function(e) {
      message(sprintf("Error installing package %s: %s", package, e$message))
    })
  }
  library(package, character.only = TRUE)
})

# Set working directory and load R script with custom functions
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Define batches to process
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")

# Loop through each batch
for (batch in batches) {
  # Set input and output directories for the current batch
  inputDir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior", batch, "OFT/SLEAP/formatted")
  outputDir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior", batch, "OFT/SLEAP/output")
  
  # Create output directory if it doesn't exist
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  
  # Read in animal ID and Code list from a text file
  animalIDCode <- read.table("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt", header = TRUE)
  
  # Create a new folder for saving the plots
  plotDir <- file.path(outputDir, "plots")
  dir.create(plotDir, showWarnings = FALSE)
  
  # Get a list of CSV files in the input directory
  fileList <- list.files(path = inputDir, pattern = "*.csv")
  
  # Initialize a list to store data frames
  dfList <- list()
  
  # Loop through each file in the input directory
  for (file in fileList) {
    # Read in tracking data from CSV file
    inputFile <- file.path(inputDir, file)
    tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)
    
    # Extract the input file name without the extension
    inputFileName <- sub(".csv$", "", basename(inputFile))
    
    # Calibrate the tracking data and add zones
    tracking <- CalibrateTrackingData(tracking, method = "area", in.metric = 49 * 49, points = c("tl", "tr", "br", "bl"))
    tracking <- AddOFTZones(tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.4, points = c("tl", "tr", "br", "bl"))
    
    # Plot zones and point data
    PlotZones(tracking)
    PlotPointData(tracking, points = c("nose"))
    
    # Perform OFT analysis on the tracking data
    tracking <- OFTAnalysis(tracking, points = "bodycentre", movement_cutoff = 5, integration_period = 5)
    
    # Extract code from filename
    code <- str_extract(inputFileName, "[A-Za-z0-9]{4}")
    
    # Extract the frequency of rear events
    frequencyRear <- sum(freqRear)
    
    # Create a data frame with the extracted data
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
    
    # Construct the output file name
    outputFile <- file.path(outputDir, paste0(inputFileName, "_output.csv"))
    
    # Write the data frame to a CSV file
    write.csv(df, outputFile, row.names = FALSE)
    
    # Add the current data frame to the list
    dfList[[length(dfList) + 1]] <- df
    
    # Plot the density paths for bodycentre
    plots <- PlotDensityPaths(tracking, points = c("bodycentre"))
    
    # Set the file name for the plot image
    plotFile <- file.path(plotDir, paste0(inputFileName, "_DensityPath.png"))
    
    # Save the plot image
    ggsave(plotFile, plot = plots$bodycentre, width = 7, height = 4)
    
    # Print progress message
    message(sprintf("Processed file %s", file))
  }
  
  # Combine all data frames into a single data frame
  dfCombined <- do.call(rbind, dfList)
  
  # Write the combined data frame to a CSV file
  write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)
  
  # Print "done" message for the current batch
  message(sprintf("Processing complete for batch %s", batch))
}
