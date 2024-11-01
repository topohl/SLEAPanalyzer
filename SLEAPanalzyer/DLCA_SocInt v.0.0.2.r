# Check if packages are installed, install them if needed, and load them
required_packages <- c("tensorflow", "reticulate", "keras", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "zoo", "stringr", "tidyverse")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Define length of video in minutes and translation to seconds
videoLength <- 2
totalTime <- videoLength * 60

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/SLEAPanalyzer/SLEAPanalzyer/")
source('DLCAnalyzer_Functions_final.R')

# Set input and output directories
inputDir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC"
outputDir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC/output"
# Create a new folder for saving the plots
plotDir <- file.path(outputDir, "plots")
dir.create(plotDir, showWarnings = FALSE)

# Get a list of CSV files in the input directory
fileList <- list.files(path = inputDir, pattern = "*.csv")

# create a new list dfList
dfList <- list()
# Function to calculate angle between two vectors
calcAngle <- function(vec1, vec2) {
  dotProd <- sum(vec1 * vec2)
  mag1 <- sqrt(sum(vec1^2))
  mag2 <- sqrt(sum(vec2^2))
  angleRad <- acos(dotProd / (mag1 * mag2))
  return(angleRad * 180 / pi)
}

# Loop through each file in the input directory
for (file in fileList) {
  # Read tracking data
  inputFile <- file.path(inputDir, file)
  Tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)

  # Initialize dfTracking
  dfTracking <- data.frame()  # Create an empty data frame to store results

  # extract length of Tracking$data$nose_1$x
  length <- length(Tracking$data$nose_1$x)
  seconds_length <- length / 30
  minutes_length <- seconds_length / 60

  # read in 017-1DLC_dlcrnetms5_SHOct29shuffle1_200000_el_filtered.csv from inpitDir
  Tracking <- ReadDLCDataFromCSV(file = "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC/017-1DLC_dlcrnetms5_SHOct29shuffle1_200000_el_filtered.csv", fps = 30)

  str(dfTracking)

  inputFileName <- sub(".csv$", "", basename(inputFile))

  bodyparts <- c("nose_1", "leftEar_1", "rightEar_1", "bodycentre_1", "leftSide_1", "rightSide_1", "tailBase_1", "tailEnd_1", "nose_2", "leftEar_2", "rightEar_2", "bodycentre_2", "leftSide_2", "rightSide_2", "tailBase_2", "tailEnd_2")
# assign bodyparts based on number behind separator "_"
bodyparts1 <- bodyparts[1:8]
bodyparts2 <- bodyparts[9:16]

  # Function to handle interpolation with specified behavior
interpolate_with_leading_trailing <- function(x) {
  # Handle leading NAs
  if (is.na(x[1])) {
    first_non_na <- which(!is.na(x))[1]
    if (!is.na(first_non_na)) {
      x[1:(first_non_na - 1)] <- x[first_non_na]
    }
  }
  
  # Handle trailing NAs
  if (is.na(x[length(x)])) {
    last_non_na <- tail(which(!is.na(x)), 1)
    if (!is.na(last_non_na)) {
      x[(last_non_na + 1):length(x)] <- x[last_non_na]
    }
  }
  
  # Interpolate internal NAs
  x <- na.approx(x, na.rm = FALSE)  # Linear interpolation between known values
  
  return(x)
}

# Interpolate missing values in x and y coordinates for all body parts
for (bodypart in bodyparts) {
  # Interpolate the x values
  Tracking$data[[bodypart]]$x <- interpolate_with_leading_trailing(Tracking$data[[bodypart]]$x)
  
  # Interpolate the y values
  Tracking$data[[bodypart]]$y <- interpolate_with_leading_trailing(Tracking$data[[bodypart]]$y)
}

  # Replace NAs in nose and bodycentre columns with last known values
  #Tracking$data$nose$x <- zoo::na.locf(Tracking$data$nose1$x)
  #Tracking$data$nose$y <- zoo::na.locf(Tracking$data$nose1$y)
  #Tracking$data$bodycentre$x <- zoo::na.locf(Tracking$data$bodycentre1$x)
  #Tracking$data$bodycentre$y <- zoo::na.locf(Tracking$data$bodycentre1$y)

  # Calibrate tracking data, add zones, and plot
  #Tracking <- CalibrateTrackingData(Tracking, method = "area", in.metric = 49 * 49, points = c("tl", "tr", "br", "bl"))
  #Tracking <- AddOFTZones(Tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.4, points = c("tl", "tr", "br", "bl"))
  #PlotZones(Tracking)
  #PlotPointData(Tracking, points = c("nose"))
  
  # Perform OFT analysis for both bodycentre_1 and bodycentre_2
  Tracking <- OFTAnalysis(Tracking, points = c("bodycentre_1", "bodycentre_2"), movement_cutoff = 5, integration_period = 5)
  
# Initialize a list to store distance values
distances <- list()

  # Create data frame for calculating angles
  dfTracking <- data.frame(
    vecNose1Bodycentre1 = cbind((Tracking$data$nose_1$x - Tracking$data$bodycentre_1$x), (Tracking$data$nose_1$y - Tracking$data$bodycentre_1$y)),
    vecNose2Bodycentre2 = cbind((Tracking$data$nose_2$x - Tracking$data$bodycentre_2$x), (Tracking$data$nose_2$y - Tracking$data$bodycentre_2$y))
  )

# Loop through each body part combination
for (i in 1:length(bodyparts1)) {
  for (j in 1:length(bodyparts2)) {
    bodypart1 <- bodyparts1[i]
    bodypart2 <- bodyparts2[j]
    
    # Extract the animal numbers
    animal1 <- substr(bodypart1, nchar(bodypart1) - 1, nchar(bodypart1))
    animal2 <- substr(bodypart2, nchar(bodypart2) - 1, nchar(bodypart2))
    
    # Check if body parts belong to different animals
    if (animal1 != animal2) {
      # Calculate distance between body parts
      # Calculate vectors between nose of animal 1 and body parts of animal 2
      for (bodypart in bodyparts2) {
        dfTracking[[paste0("vecNose1", bodypart)]] <- cbind((Tracking$data$nose_1$x - Tracking$data[[bodypart]]$x), (Tracking$data$nose_1$y - Tracking$data[[bodypart]]$y))
      }

      # Calculate vectors between nose of animal 2 and body parts of animal 1
      for (bodypart in bodyparts1) {
        dfTracking[[paste0("vecNose2", bodypart)]] <- cbind((Tracking$data$nose_2$x - Tracking$data[[bodypart]]$x), (Tracking$data$nose_2$y - Tracking$data[[bodypart]]$y))
      }
      
    }
  }
}

# calculate angles between vecNose1BodyCentre1 and bodyparts of animal 2, also calculate angles between vecNose2BodyCentre2 and bodyparts of animal 1
for (bodypart in bodyparts2) {
  dfTracking[[paste0("angleDegNose1", bodypart)]] <- numeric(length = nrow(dfTracking))
  for (i in 1:nrow(dfTracking)) {
    vecNose1Bodycentre1 <- c(dfTracking$vecNose1Bodycentre1.1[i], dfTracking$vecNose1Bodycentre1.2[i])
    vecNose1Bodypart2 <- c(dfTracking[[paste0("vecNose1", bodypart)]][i, 1], dfTracking[[paste0("vecNose1", bodypart)]][i, 2])
    
    dfTracking[[paste0("angleDegNose1", bodypart)]][i] <- calcAngle(vecNose1Bodycentre1, vecNose1Bodypart2)
  }
}

for (bodypart in bodyparts1) {
  dfTracking[[paste0("angleDegNose2", bodypart)]] <- numeric(length = nrow(dfTracking))
  for (i in 1:nrow(dfTracking)) {
    vecNose2Bodycentre2 <- c(dfTracking$vecNose2Bodycentre2.1[i], dfTracking$vecNose2Bodycentre2.2[i])
    vecNose2Bodypart1 <- c(dfTracking[[paste0("vecNose2", bodypart)]][i, 1], dfTracking[[paste0("vecNose2", bodypart)]][i, 2])
    
    dfTracking[[paste0("angleDegNose2", bodypart)]][i] <- calcAngle(vecNose2Bodycentre2, vecNose2Bodypart1)
  }
}


  # calculate interaction of animal 1 with animal 2 nose
    isInNoseToNoseContactanimal1 <- GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270
    contactNose1toNose2 <- sum(isInNoseToNoseContactanimal1) * seconds_length / length(isInNoseToNoseContactanimal1)
    
    isInNoseToBodycentreContactanimal1 <- GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270
    contactNose1toBodycentre2 <- sum(isInNoseToBodycentreContactanimal1) * seconds_length / length(isInNoseToBodycentreContactanimal1)

    isInNoseToTailBaseContactanimal1 <- GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270
    contactNose1toTailBase2 <- sum(isInNoseToTailBaseContactanimal1) * seconds_length / length(isInNoseToTailBaseContactanimal1)

    sideBySideContactAnimal1 <- GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80
    contactSideBySideAnimal1 <- sum(sideBySideContactAnimal1) * seconds_length / length(sideBySideContactAnimal1)

    sideBySideReverseContactAnimal1 <- GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80 
    contactSideBySideReverseAnimal1 <- sum(sideBySideReverseContactAnimal1) * seconds_length / length(sideBySideReverseContactAnimal1)


    # calculate interaction of animal 1 with animal 2 nose
    isInNoseToNoseContactanimal2 <- GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270
    contactNose1toNose1 <- sum(isInNoseToNoseContactanimal1) * seconds_length / length(isInNoseToNoseContactanimal1)
    
    isInNoseToBodycentreContactanimal2 <- GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270
    contactNose1toBodycentre2 <- sum(isInNoseToBodycentreContactanimal1) * seconds_length / length(isInNoseToBodycentreContactanimal1)

    isInNoseToTailBaseContactanimal2 <- GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270
    contactNose1toTailBase2 <- sum(isInNoseToTailBaseContactanimal2) * seconds_length / length(isInNoseToTailBaseContactanimal2)

    sideBySideContactAnimal2 <- GetDistances(Tracking, "nose_2", "nose_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "tailBase_1") <= 80
    contactSideBySideAnimal2 <- sum(sideBySideContactAnimal2) * seconds_length / length(sideBySideContactAnimal2)

    sideBySideReverseContactAnimal2 <- GetDistances(Tracking, "nose_2", "tailBase_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "nose_1") <= 80
    contactSideBySideReverseAnimal2 <- sum(sideBySideReverseContactAnimal2) * seconds_length / length(sideBySideReverseContactAnimal2)

  # calculate latency until first contact between nose_1 and any bodypart of animal 2 in seconds
    latencyNose1toNose2 <- min(which(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270)) * 1/30
    latencyNose1toBodycentre2 <- min(which(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270)) * 1/30
    latencyNose1toTailBase2 <- min(which(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270)) * 1/30
    latencySideBySideAnimal1 <- min(which(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80)) * 1/30
    latencySideBySideReverseAnimal1 <- min(which(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80)) * 1/30

  # calculate frequency of contact between nose_1 and any bodypart of animal 2
    isInNoseToNoseContactanimal1 <- cumsum(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270) == 1
    isInNoseToBodycentreContactanimal1 <- cumsum(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270) == 1
    isInNoseToTailBaseContactanimal1 <- cumsum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270) == 1
    sideBySideContactAnimal1 <- cumsum(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80) == 1
    sideBySideReverseContactAnimal1 <- cumsum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80) == 1
    frequencyNose1toNose2 <- sum(isInNoseToNoseContactanimal1)
    frequencyNose1toBodycentre2 <- sum(isInNoseToBodycentreContactanimal1)
    frequencyNose1toTailBase2 <- sum(isInNoseToTailBaseContactanimal1)
    frequencySideBySideAnimal1 <- sum(sideBySideContactAnimal1)
    frequencySideBySideReverseAnimal1 <- sum(sideBySideReverseContactAnimal1)

    # calculate proximity of animal 1 to animal 2
    isInProxAnimal1 <- GetDistances(Tracking, "nose_1", "nose_2") > 4 & GetDistances(Tracking, "nose_1", "nose_2") <= 8
    proxAnimal1 <- sum(isInProxAnimal1) * seconds_length / length(isInProxAnimal1)

  # calculate proximity with orientation towards other animal (angle)
    isInProxAnimal1Angle <- GetDistances(Tracking, "nose_1", "nose_2") > 4 & GetDistances(Tracking, "nose_1", "nose_2") <= 8 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270
    proxAnimal1Angle <- sum(isInProxAnimal1Angle) * seconds_length / length(isInProxAnimal1Angle)
  
  # check which animal initiated first contact between nose and bodypart
    isInNoseToNoseContactanimal1 <- cumsum(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270) == 1
    isInNoseToBodycentreContactanimal1 <- cumsum(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270) == 1
    isInNoseToTailBaseContactanimal1 <- cumsum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270) == 1
    sideBySideContactAnimal1 <- cumsum(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80) == 1
    sideBySideReverseContactAnimal1 <- cumsum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80) == 1
    isInNoseToNoseContactanimal2 <- cumsum(GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270) == 1
    isInNoseToBodycentreContactanimal2 <- cumsum(GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose2bodycentre_1) >= 90 & abs(dfTracking$angleDegNose2bodycentre_1) <= 270) == 1
    isInNoseToTailBaseContactanimal2 <- cumsum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose2tailBase_1) >= 90 & abs(dfTracking$angleDegNose2tailBase_1) <= 270) == 1
    sideBySideContactAnimal2 <- cumsum(GetDistances(Tracking, "nose_2", "nose_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "tailBase_1") <= 80) == 1
    sideBySideReverseContactAnimal2 <- cumsum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "nose_1") <= 80) == 1
    isInProxAnimal1 <- GetDistances(Tracking, "nose_1", "nose_2") > 4 & GetDistances(Tracking, "nose_1", "nose_2") <= 8
  
  # Calculate frequency of rearings
  #freqRear <- cumsum(GetDistances(Tracking, "spine1", "bodycentre1") <= 1) > 0 &
  #  cumsum(GetDistances(Tracking, "bodycentre", "spine2") <= 1) > 0 &
  #  c(0, diff(GetDistances(Tracking, "spine1", "bodycentre") <= 1)) == 1 &
  #  c(0, diff(GetDistances(Tracking, "bodycentre", "spine2") <= 1)) == 1
  #frequencyRear <- sum(freqRear)
  
  # Create data frame with results
  df <- data.frame(
    file = inputFileName,
    contactNose1toNose2 = contactNose1toNose2,
    contactNose1toBodycentre2 = contactNose1toBodycentre2,
    contactNose1toTailBase2 = contactNose1toTailBase2,
    contactSideBySideAnimal1 = contactSideBySideAnimal1,
    contactSideBySideReverseAnimal1 = contactSideBySideReverseAnimal1,
    latencyNose1toNose2 = latencyNose1toNose2,
    latencyNose1toBodycentre2 = latencyNose1toBodycentre2,
    latencyNose1toTailBase2 = latencyNose1toTailBase2,
    latencySideBySideAnimal1 = latencySideBySideAnimal1,
    latencySideBySideReverseAnimal1 = latencySideBySideReverseAnimal1,
    frequencyNose1toNose2 = frequencyNose1toNose2,
    frequencyNose1toBodycentre2 = frequencyNose1toBodycentre2,
    frequencyNose1toTailBase2 = frequencyNose1toTailBase2,
    frequencySideBySideAnimal1 = frequencySideBySideAnimal1,
    frequencySideBySideReverseAnimal1 = frequencySideBySideReverseAnimal1,
    proxAnimal1 = proxAnimal1,
    proxAnimal1Angle = proxAnimal1Angle,
    distance = Tracking$Report$bodycentre.raw.distance,
    stationary = Tracking$Report$bodycentre.time.stationary,
    speedMoving = Tracking$Report$bodycentre.speed.moving,
    speedRaw = Tracking$Report$bodycentre.raw.speed,
    frequencyRear = frequencyRear
  )
  
  # Write data frame to output file
  outputFile <- paste0(inputFileName, "_output.csv")
  write.csv(df, outputFile, row.names = FALSE)
  
  # Add data frame to list
  dfList[[length(dfList) + 1]] <- df
  
  # Save density paths plot
  plots <- PlotDensityPaths(Tracking, points = c("bodycentre"))
  plotFile <- file.path(plotDir, paste0(inputFileName, "_DensityPath.png"))
  ggsave(plotFile, plot = plots$bodycentre, width = 7, height = 6)
  
  # Print progress message
  cat(sprintf("Processed file %s\n", file))
}

# Combine all data frames into a single data frame
dfCombined <- do.call(rbind, dfList)

# Write combined data frame to csv file
write.csv(dfCombined, file.path(outputDir, "combined_output_test.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")
