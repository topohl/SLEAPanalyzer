# Check if packages are installed, install them if needed, and load them
required_packages <- c("tensorflow", "reticulate", "keras", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "zoo")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Set input and output directories
inputDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B3/NOR/SLEAP/formatted"
outputDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B3/NOR/SLEAP/testoutput"
# Create a new folder for saving the plots
plotDir <- file.path(outputDir, "plots")
dir.create(plotDir, showWarnings = FALSE)

# Get a list of CSV files in the input directory
fileList <- list.files(path = inputDir, pattern = "*.csv")

# create a new list dfList
dfList <- list()

calcAngle <- function(vecNose, vecBodyCenter) {
    # Calculate dot product
    dotProd <- sum(vecNose * vecBodyCenter)

    # Calculate magnitudes
    magA <- sqrt(sum(vecNose^2))
    magC <- sqrt(sum(vecBodyCenter^2))

    # Calculate the angle in radians
    angleRad <- acos(dotProd / (magA * magC))

    # Convert angle to degrees and return
    return(angleRad * 180 / pi)
}

# Loop through each file in the input directory
for (file in fileList) {
  
  # Read in tracking data and get names of variables in the data frame
  inputFile <- file.path(inputDir, file)
  Tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)
  
  # Extract the input file name without the extension
  inputFileName <- sub(".csv$", "", basename(inputFile))
  
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
  dfTracking <- data.frame(
    vecNoseBodycentre = cbind((Tracking$data$nose$x - Tracking$data$bodycentre$x), (Tracking$data$nose$y - Tracking$data$bodycentre$y)),
    vecNoseObjL = cbind((Tracking$data$nose$x - Tracking$data$objL$x), (Tracking$data$nose$y - Tracking$data$objL$y)),
    vecNoseObjR = cbind((Tracking$data$nose$x - Tracking$data$objR$x), (Tracking$data$nose$y - Tracking$data$objR$y)),
    circleL = sqrt((Tracking$data$objL$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$objL$y - Tracking$data$bodycentre$y)^2),
    circleR = sqrt((Tracking$data$objR$x - Tracking$data$bodycentre$x)^2 + (Tracking$data$objR$y - Tracking$data$bodycentre$y)^2)
  )
  
  # Create new variables in the data frame to indicate if the circle is within 4cm
  dfTracking$circleL <- ifelse(dfTracking$circleL <= 4, 1, 0)
  dfTracking$circleR <- ifelse(dfTracking$circleR <= 4, 1, 0)
  
  # Initialize vectors to store angles
  angleDegObjL <- numeric(length = nrow(dfTracking))
  angleDegObjR <- numeric(length = nrow(dfTracking))
  
# Loop through each row of dfTracking
for (i in 1:nrow(dfTracking)) {
    # Extract vectors a and c from vectorNoseObjL and vector matrices
    vecNoseObjL <- c(dfTracking$vecNoseObjL.1[i], dfTracking$vecNoseObjL.2[i])
    vecNoseObjR <- c(dfTracking$vecNoseObjR.1[i], dfTracking$vecNoseObjR.2[i])
    vecBodyCentre <- c(dfTracking$vecNoseBodycentre.1[i], dfTracking$vecNoseBodycentre.2[i])

    # Calculate angles and store in the vectors
    angleDegObjL[i] <- calcAngle(vecNoseObjL, vecBodyCentre)
    angleDegObjR[i] <- calcAngle(vecNoseObjR, vecBodyCentre)
}
   
  # Calculate time_Left and time_Right only if angle L or R is between 90 and 270 degrees
  totalTime <- 10 * 60
  isInZoneLeft <- GetDistances(Tracking, "objL", "nose") <= 4 & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(angleDegObjL) >= 90 & abs(angleDegObjL) <= 270
  isInZoneRight <- GetDistances(Tracking, "objR", "nose") <= 4 & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(angleDegObjR) >= 90 & abs(angleDegObjR) <= 270
  contactLeft <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
  contactRight <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
  
  # Calculate time with nose between 6-10 cm of either objL or objR
  isInProxLeft <- GetDistances(Tracking, "objL", "nose") > 4 & GetDistances(Tracking, "objL", "nose") <= 8
  isInProxRight <- GetDistances(Tracking, "objR", "nose") > 4 & GetDistances(Tracking, "objR", "nose") <= 8
  proxLeft <- sum(isInProxLeft) * totalTime / length(isInProxLeft)
  proxRight <- sum(isInProxRight) * totalTime / length(isInProxRight)
  
  # Calculate frequency of zone visits for objL and nose
  isInZoneVisitL <- cumsum(GetDistances(Tracking, "objL", "nose") <= 4) > 0 & c(0, diff(GetDistances(Tracking, "objL", "nose") <= 4)) == 1
  isInZoneVisitR <- cumsum(GetDistances(Tracking, "objR", "nose") <= 4) > 0 & c(0, diff(GetDistances(Tracking, "objR", "nose") <= 4)) == 1
  frequencyL <- sum(isInZoneVisitL)
  frequencyR <- sum(isInZoneVisitR)
  
  # Calculate frequency of rearings
  freqRear <- cumsum(GetDistances(Tracking, "spine1", "bodycentre") <= 1) > 0 &
    cumsum(GetDistances(Tracking, "bodycentre", "spine2") <= 1) > 0 &
    c(0, diff(GetDistances(Tracking, "spine1", "bodycentre") <= 1)) == 1 &
    c(0, diff(GetDistances(Tracking, "bodycentre", "spine2") <= 1)) == 1
  frequencyRear <- sum(freqRear)
  
  # Create a data frame with the time_Left, time_Right, time_nose_Left, and time_nose_Right values
  df <- data.frame(
    file = inputFileName,
    contactLeft = contactLeft,
    contactRight = contactRight,
    proxLeft = proxLeft,
    proxRight = proxRight,
    frequencyL = frequencyL,
    frequencyR = frequencyR,
    distance = Tracking$Report$bodycentre.raw.distance,
    stationary = Tracking$Report$bodycentre.time.stationary,
    speedMoving = Tracking$Report$bodycentre.speed.moving,
    speedRaw = Tracking$Report$bodycentre.raw.speed,
    frequencyRear = frequencyRear
  )
  
  # Construct the output file name
  outputFile <- paste0(inputFileName, "_output.csv")
  
  # Write the data frame to a csv file with the constructed file name
  write.csv(df, outputFile, row.names = FALSE)
  
  # Add the current data frame to the list
  dfList[[length(dfList) + 1]] <- df
  
  # Plot the density paths for bodycentre
  plots <- PlotDensityPaths(Tracking, points = c("bodycentre"))
  
  # Set the file name for the plot image
  plotFile <- file.path(plotDir, paste0(inputFileName, "_DensityPath.png"))
  
    # Save the plot image with the file name
    ggsave(plotFile, plot = plots$bodycentre, width = 7, height = 6)
    
    # Print progress message
    cat(sprintf("Processed file %s\n", file))
}

# Combine all data frames into a single data frame
dfCombined <- do.call(rbind, dfList)

# Write the combined data frame to a csv file
write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")
