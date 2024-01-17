# Check if packages are installed, install them if needed, and load them
required_packages <- c("tensorflow", "reticulate", "keras", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "zoo", "stringr", "tidyverse")

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
inputDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B4/NOR/SLEAP/formatted"
outputDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B4/NOR/SLEAP/output_new"
novelLocDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B4/NOR"
# Create a new folder for saving the plots
plotDir <- file.path(outputDir, "plots")
dir.create(plotDir, showWarnings = FALSE)

# Get a list of CSV files in the input directory
fileList <- list.files(path = inputDir, pattern = "*.csv")

# only use files with "NOV" in filename
# fileList <- fileList[str_detect(fileList, "NOV")]

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

# Read novelLoc file
novelLoc <- read.table(file.path(novelLocDir, "novelLoc.txt"), header = TRUE, sep = "\t")

# Loop through each file in the input directory
for (file in fileList) {
  # Read tracking data
  inputFile <- file.path(inputDir, file)
  Tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)
  inputFileName <- sub(".csv$", "", basename(inputFile))
  
  # Replace NAs in nose and bodycentre columns with last known values
  Tracking$data$nose$x <- zoo::na.locf(Tracking$data$nose$x)
  Tracking$data$nose$y <- zoo::na.locf(Tracking$data$nose$y)
  Tracking$data$bodycentre$x <- zoo::na.locf(Tracking$data$bodycentre$x)
  Tracking$data$bodycentre$y <- zoo::na.locf(Tracking$data$bodycentre$y)
  
  # Calibrate tracking data, add zones, and plot
  Tracking <- CalibrateTrackingData(Tracking, method = "area", in.metric = 49 * 49, points = c("tl", "tr", "br", "bl"))
  Tracking <- AddOFTZones(Tracking, scale_center = 0.5, scale_periphery = 0.8, scale_corners = 0.4, points = c("tl", "tr", "br", "bl"))
  PlotZones(Tracking)
  PlotPointData(Tracking, points = c("nose"))
  
  # Perform OFT analysis
  Tracking <- OFTAnalysis(Tracking, points = "bodycentre", movement_cutoff = 5, integration_period = 5)
  
  # Create data frame for calculating angles
  dfTracking <- data.frame(
    vecNoseBodycentre = cbind((Tracking$data$nose$x - Tracking$data$bodycentre$x), (Tracking$data$nose$y - Tracking$data$bodycentre$y)),
    vecNoseObjL = cbind((Tracking$data$nose$x - Tracking$data$objL$x), (Tracking$data$nose$y - Tracking$data$objL$y)),
    vecNoseObjR = cbind((Tracking$data$nose$x - Tracking$data$objR$x), (Tracking$data$nose$y - Tracking$data$objR$y))
  )
  
  # Initialize vectors to store angles
  angleDegObjL <- numeric(length = nrow(dfTracking))
  angleDegObjR <- numeric(length = nrow(dfTracking))
  
  # Calculate angles for each row in dfTracking
  for (i in 1:nrow(dfTracking)) {
    vecNoseObjL <- c(dfTracking$vecNoseObjL.1[i], dfTracking$vecNoseObjL.2[i])
    vecNoseObjR <- c(dfTracking$vecNoseObjR.1[i], dfTracking$vecNoseObjR.2[i])
    vecBodyCentre <- c(dfTracking$vecNoseBodycentre.1[i], dfTracking$vecNoseBodycentre.2[i])
    
    angleDegObjL[i] <- calcAngle(vecNoseObjL, vecBodyCentre)
    angleDegObjR[i] <- calcAngle(vecNoseObjR, vecBodyCentre)
  }
  
  # Define square coordinates around objL and objR
  squareWidth <- 9
  squareHeight <- 7
  objLSquare <- data.frame(
    x = c(Tracking$data$objL$x - squareWidth/2, Tracking$data$objL$x + squareWidth/2, Tracking$data$objL$x + squareWidth/2, Tracking$data$objL$x - squareWidth/2),
    y = c(Tracking$data$objL$y - squareHeight/2, Tracking$data$objL$y - squareHeight/2, Tracking$data$objL$y + squareHeight/2, Tracking$data$objL$y + squareHeight/2)
  )
  objRSquare <- data.frame(
    x = c(Tracking$data$objR$x - squareWidth/2, Tracking$data$objR$x + squareWidth/2, Tracking$data$objR$x + squareWidth/2, Tracking$data$objR$x - squareWidth/2),
    y = c(Tracking$data$objR$y - squareHeight/2, Tracking$data$objR$y - squareHeight/2, Tracking$data$objR$y + squareHeight/2, Tracking$data$objR$y + squareHeight/2)
  )
  
  # Extract code from filename
  code <- str_extract(inputFileName, "^[A-Za-z0-9]{4}")
  
  # Check if code exists in novelLoc file
  if (code %in% novelLoc$Code) {
    totalTime <- 10 * 60
    if (novelLoc[novelLoc$Code == code, "NovelLoc"] == "R") {
      isInZoneLeft <- point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objLSquare$x, objLSquare$y) & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(angleDegObjL) >= 70 & abs(angleDegObjL) <= 290
      isInZoneRight <- GetDistances(Tracking, "objR", "nose") <= 4 & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(angleDegObjR) >= 70 & abs(angleDegObjR) <= 290
          contactNov <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
          contactFam <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
    } else {
      isInZoneLeft <- GetDistances(Tracking, "objL", "nose") <= 4 & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(angleDegObjL) >= 70 & abs(angleDegObjL) <= 290
      isInZoneRight <- point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objRSquare$x, objRSquare$y) & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(angleDegObjR) >= 70 & abs(angleDegObjR) <= 290
          contactFam <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
          contactNov <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
    }
    contactLeft <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
    contactRight <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
  } else {
    contactLeft <- 0
    contactRight <- 0
  }
  
  # Calculate time with nose between 6-10 cm of objL or objR
  isInProxLeft <- GetDistances(Tracking, "objL", "nose") > 4 & GetDistances(Tracking, "objL", "nose") <= 8
  isInProxRight <- GetDistances(Tracking, "objR", "nose") > 4 & GetDistances(Tracking, "objR", "nose") <= 8
  proxLeft <- sum(isInProxLeft) * totalTime / length(isInProxLeft)
  proxRight <- sum(isInProxRight) * totalTime / length(isInProxRight)
  
  # Calculate time with nose between 6-10 cm of objL or objR with consideration of head angle
  isInProxLeftAngle <- GetDistances(Tracking, "objL", "nose") > 4 & GetDistances(Tracking, "objL", "nose") <= 8 & abs(angleDegObjL) >= 90 & abs(angleDegObjL) <= 270
  isInProxRightAngle <- GetDistances(Tracking, "objR", "nose") > 4 & GetDistances(Tracking, "objR", "nose") <= 8 & abs(angleDegObjL) >= 90 & abs(angleDegObjL) <= 270
  proxLeftAngle <- sum(isInProxLeftAngle) * totalTime / length(isInProxLeftAngle)
  proxRightAngle <- sum(isInProxRightAngle) * totalTime / length(isInProxRightAngle)
  
  # Calculate latency until first contact with objL or objR in seconds
  if (novelLoc[novelLoc$Code == code, "NovelLoc"] == "R") {
    latencyLeft <- min(which(point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objLSquare$x, objLSquare$y) & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(angleDegObjL) >= 90 & abs(angleDegObjL) <= 270)) * 1/30
    latencyRight <- min(which(GetDistances(Tracking, "objR", "nose") <= 4 & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(angleDegObjR) >= 70 & abs(angleDegObjR) <= 290)) * 1/30
  } else {
    latencyLeft <- min(which(GetDistances(Tracking, "objL", "nose") <= 4 & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(angleDegObjL) >= 70 & abs(angleDegObjL) <= 290)) * 1/30
    latencyRight <- min(which(point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objRSquare$x, objRSquare$y) & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(angleDegObjR) >= 90 & abs(angleDegObjR) <= 270)) * 1/30
  }
  
  # Check latencyLeft and latencyRight and compare both to note first contact, save latency into data frame
  if (is.na(latencyLeft) & is.na(latencyRight)) {
    latency <- NA
  } else if (is.na(latencyLeft)) {
    latency <- latencyRight
  } else if (is.na(latencyRight)) {
    latency <- latencyLeft
  } else {
    latency <- min(latencyLeft, latencyRight)
  }
  
  # Calculate frequency of zone visits for nose and objL or objR contact zones
  if (novelLoc[novelLoc$Code == code, "NovelLoc"] == "R") {
    isInZoneVisitL <- cumsum(point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objLSquare$x, objLSquare$y) & c(0, diff(point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objLSquare$x, objLSquare$y)))) == 1
    isInZoneVisitR <- cumsum(GetDistances(Tracking, "objR", "nose") <= 4 & c(0, diff(GetDistances(Tracking, "objR", "nose") <= 4))) == 1
  } else {
    isInZoneVisitL <- cumsum(GetDistances(Tracking, "objL", "nose") <= 4 & c(0, diff(GetDistances(Tracking, "objL", "nose") <= 4))) == 1
    isInZoneVisitR <- cumsum(point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objRSquare$x, objRSquare$y) & c(0, diff(point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objRSquare$x, objRSquare$y)))) == 1
  }
  frequencyL <- sum(isInZoneVisitL)
  frequencyR <- sum(isInZoneVisitR)
  
  # Calculate frequency of rearings
  freqRear <- cumsum(GetDistances(Tracking, "spine1", "bodycentre") <= 1) > 0 &
    cumsum(GetDistances(Tracking, "bodycentre", "spine2") <= 1) > 0 &
    c(0, diff(GetDistances(Tracking, "spine1", "bodycentre") <= 1)) == 1 &
    c(0, diff(GetDistances(Tracking, "bodycentre", "spine2") <= 1)) == 1
  frequencyRear <- sum(freqRear)
  
  # Create data frame with results
  df <- data.frame(
    file = inputFileName,
    contactLeft = contactLeft,
    contactRight = contactRight,
    contactNov = contactNov,
    contactFam = contactFam,
    proxLeft = proxLeft,
    proxRight = proxRight,
    proxLeftAngle = proxLeftAngle,
    proxRightAngle = proxRightAngle,
    latency = latency,
    latencyLeft = latencyLeft,
    latencyRight = latencyRight,
    frequencyL = frequencyL,
    frequencyR = frequencyR,
    distance = Tracking$Report$bodycentre.raw.distance,
    stationary = Tracking$Report$bodycentre.time.stationary,
    speedMoving = Tracking$Report$bodycentre.speed.moving,
    speedRaw = Tracking$Report$bodycentre.raw.speed,
    frequencyRear = frequencyRear,
    novelLoc = novelLoc[novelLoc$Code == code, "NovelLoc"]
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
write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")
