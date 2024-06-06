# Install required packages
requiredPackages <- c("stringr", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "keras", "tensorflow", "zoo", "dplyr", "reticulate")

# Check if packages are installed, if not install and load them
for (package in requiredPackages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Define batches
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")

# Loop through each batch
for (batch in batches) {
  # Set input and output directories depending on the batch
  inputDir <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/OFT/SLEAP/formatted/")
  outputDir <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/OFT/SLEAP/output/")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  
  # Read in animal ID and Code list from .txt file
  animalIDCode <- read.table("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Planning/animalIDCode.txt", header = TRUE)
  
  # Create a new folder for saving the plots
  plotDir <- file.path(outputDir, "plots")
  dir.create(plotDir, showWarnings = FALSE)
  
  # Get a list of CSV files in the input directory
  fileList <- list.files(path = inputDir, pattern = "*.csv")
  
  # Create new list dfList
  dfList <- list()
  
  # Loop through each file in the input directory
  for (file in fileList) {
    # Read in tracking data and get names of variables in the data frame
    inputFile <- file.path(inputDir, file)
    tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)
    
    # Extract the input file name without the extension
    inputFileName <- sub(".csv$", "", basename(inputFile))
    
    # Calibrate the tracking data, add zones, and plot the zones
    tracking <- CalibrateTrackingData(tracking, method = "area", in.metric = 49 * 49, points = c("tl", "tr", "br", "bl"))
    tracking <- AddOFTZones(tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.4, points = c("tl", "tr", "br", "bl"))
    PlotZones(tracking)
    PlotPointData(tracking, points = c("nose"))
    
    tracking <- OFTAnalysis(tracking, points = "bodycentre", movement_cutoff = 5, integration_period = 5)
    
    # Extract code from filename
    code <- str_extract(inputFileName, "[A-Za-z0-9]{4}")
    
    # Calculate frequency of rearings
    freqRear <- cumsum(GetDistances(tracking, "spine1", "bodycentre") <= 1) > 0 &
      cumsum(GetDistances(tracking, "bodycentre", "spine2") <= 1) > 0 &
      c(0, diff(GetDistances(tracking, "spine1", "bodycentre") <= 1)) == 1 &
      c(0, diff(GetDistances(tracking, "bodycentre", "spine2") <= 1)) == 1
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
    outputFile <- paste0(inputFileName, "_output.csv")
    
    # Write the data frame to a csv file with the constructed file name
    write.csv(df, outputFile, row.names = FALSE)
    
    # Add the current data frame to the list
    dfList[[length(dfList) + 1]] <- df
    
    # Plot the density paths for bodycentre
    plots <- PlotDensityPaths(tracking, points = c("bodycentre"))
    
    # Set the file name for the plot image
    plotFile <- file.path(plotDir, paste0(inputFileName, "_DensityPath.png"))
    
    # Save the plot image with the file name
    ggsave(plotFile, plot = plots$bodycentre, width = 7, height = 4)
    
    # Print progress message
    cat(sprintf("Processed file %s\n", file))
  }
  
  # Combine all data frames into a single data frame
  dfCombined <- do.call(rbind, dfList)
  
  # Write the combined data frame to a csv file
  write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)
  
  # Print "done" message
  cat("Processing complete for batch", batch, "\n")
}
