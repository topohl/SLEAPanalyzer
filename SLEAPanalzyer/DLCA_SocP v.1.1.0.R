#' @title Batch Processing of Social Preference Test Data
#' @description
#' This script processes behavioral tracking data from a social preference test.
#' It loads required packages, iterates over different batches and test phases,
#' reads and processes tracking data, calculates time spent in different zones,
#' and generates output files with the extracted metrics and visualizations.
#'
#' @details
#' - Installs and loads required packages.
#' - Sets working directory and sources necessary functions.
#' - Iterates through experimental batches and social preference phases.
#' - Reads animal ID codes and novel location data.
#' - Processes tracking data using DLC output and extracts behavioral metrics.
#' - Calculates time in novel and familiar zones, proximity measures, and visit frequency.
#' - Saves processed data and visualizations.
#'
#' @note
#' - Ensure all required R packages are installed.
#' - Update the file paths to match the actual directory structure.
#' - This script assumes data is structured in specific formatted CSV and TXT files.
#'
#' @author Tobias Pohl
#' @date 2025-02-06
#' @version 1.1.0

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(stringr, sp, imputeTS, ggplot2, ggmap, data.table, cowplot, corrplot, keras, tensorflow, zoo, dplyr, reticulate)

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Define batch and socpPhase
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
socpPhases <- c("HAB", "S1", "S2")

# Loop through each batch
for (batch in batches) {
  # Loop through each socpPhase
  for (socpPhase in socpPhases) {
    # Set input and output directories depending on the batch and socpPhase
    inputDir <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/SocP/SLEAP/formatted/", socpPhase)
    outputDir <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/SocP/SLEAP/output/", socpPhase)
    novelLocDir <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/SocP")
    
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
    
    # Read novelLoc file
    novelLoc <- read.table(file.path(novelLocDir, paste0("novelLoc", socpPhase, ".txt")), header = TRUE, sep = "\t")
    
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
      tracking <- CalibrateTrackingData(tracking, method = "area", in.metric = 44 * 24, points = c("tl", "tr", "br", "bl"))
      tracking <- AddOFTZones(tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.4, points = c("tl", "tr", "br", "bl"))
      PlotZones(tracking)
      PlotPointData(tracking, points = c("nose"))
      
      tracking <- OFTAnalysis(tracking, points = "bodycentre", movement_cutoff = 5, integration_period = 5)
      
      # Extract code from filename
      code <- str_extract(inputFileName, "^[A-Za-z0-9]{4}")
      
      # Calculate time in contact with novel and familiar objects
      if (code %in% novelLoc$Code) {
        totalTime <- 10 * 60
        
        if (novelLoc[novelLoc$Code == code, "NovelLoc"] == "R") {
          isInZoneLeft <- GetDistances(tracking, "socl", "nose") <= 6 & GetDistances(tracking, "socl", "bodycentre") > 1
          isInZoneRight <- GetDistances(tracking, "socr", "nose") <= 6 & GetDistances(tracking, "socr", "bodycentre") > 1
          contactNov <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
          contactFam <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
        } else {
          isInZoneLeft <- GetDistances(tracking, "socl", "nose") <= 6 & GetDistances(tracking, "socl", "bodycentre") > 1
          isInZoneRight <- GetDistances(tracking, "socr", "nose") <= 6 & GetDistances(tracking, "socr", "bodycentre") > 1
          contactFam <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
          contactNov <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
        }
        
        contactLeft <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
        contactRight <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
      } else {
        contactLeft <- 0
        contactRight <- 0
      }
      
      # Check if the nose is in the contact zone of objR
      isContactsocr <- ifelse(isInZoneRight, "contact socr", "")
      # Check if the nose is in the contact zone of objL
      isContactsocl <- ifelse(isInZoneLeft, "contact socl", "")
      
      # Calculate time with nose between 6-10 cm of either objL or objR
      if (code %in% novelLoc$Code) {
        totalTime <- 10 * 60
        if (novelLoc[novelLoc$Code == code, "NovelLoc"] == "R") {
          isInProxLeft <- GetDistances(tracking, "socl", "nose") > 6 & GetDistances(tracking, "socl", "nose") <= 10
          isInProxRight <- GetDistances(tracking, "socr", "nose") > 6 & GetDistances(tracking, "socr", "nose") <= 10
          proxNov <- sum(isInProxLeft) * totalTime / length(isInProxLeft)
          proxFam <- sum(isInProxRight) * totalTime / length(isInProxRight)
        } else {
          isInProxLeft <- GetDistances(tracking, "socl", "nose") > 6 & GetDistances(tracking, "socl", "nose") <= 10
          isInProxRight <- GetDistances(tracking, "socr", "nose") > 6 & GetDistances(tracking, "socr", "nose") <= 10
          proxFam <- sum(isInProxLeft) * totalTime / length(isInProxLeft)
          proxNov <- sum(isInProxRight) * totalTime / length(isInProxRight)
        }
      } else {
        proxNov <- 0
        proxFam <- 0
      }
      
      # Calculate frequency of zone visits for objL and nose
      isInZoneVisitL <- cumsum(GetDistances(tracking, "socl", "nose") <= 6) > 0 & c(0, diff(GetDistances(tracking, "socl", "nose") <= 6)) == 1
      frequencyL <- sum(isInZoneVisitL)
      isInZoneVisitR <- cumsum(GetDistances(tracking, "socr", "nose") <= 6) > 0 & c(0, diff(GetDistances(tracking, "socr", "nose") <= 6)) == 1
      frequencyR <- sum(isInZoneVisitR)
      
      # Create a data frame with the timeLeft, timeRight, timeNoseLeft, and timeNoseRight values
      df <- data.frame(
        file = inputFileName,
        ID = animalIDCode[animalIDCode$Code == code, "ID"],
        Code = code,
        contactNovel = contactNov,
        contactFamiliar = contactFam,
        proxNovel = proxNov,
        proxFamiliar = proxFam,
        contactLeft = contactLeft,
        contactRight = contactRight,
        frequencyLeft = frequencyL,
        frequencyRight = frequencyR,
        distance = tracking$Report$bodycentre.raw.distance,
        stationary = tracking$Report$bodycentre.time.stationary,
        speedMoving = tracking$Report$bodycentre.speed.moving,
        speedRaw = tracking$Report$bodycentre.raw.speed
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
    cat("Processing complete for batch", batch, "and socpPhase", socpPhase, "\n")
  }
}
