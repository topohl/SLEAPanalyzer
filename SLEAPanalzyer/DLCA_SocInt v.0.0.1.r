# Check if packages are installed, install them if needed, and load them
required_packages <- c("tensorflow", "reticulate", "keras", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "zoo", "stringr", "tidyverse")

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Define length of video in minutes and translation to seconds
videoLength <- 10
totalTime <- videoLength * 60

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Set input and output directories
inputDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B1/NOR/SLEAP/formatted"
outputDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B1/NOR/SLEAP/output_angle"
novelLocDir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B1/NOR"
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
    vecNose1Bodycentre1 = cbind((Tracking$data$nose1$x - Tracking$data$bodycentre1$x), (Tracking$data$nose1$y - Tracking$data$bodycentre1$y)),
    vecNose1Nose2 = cbind((Tracking$data$nose1$x - Tracking$data$nose2$x), (Tracking$data$nose1$y - Tracking$data$nose2$y)),
    vecNose2Bodycentre2 = cbind((Tracking$data$nose2$x - Tracking$data$bodycentre2$x), (Tracking$data$nose2$y - Tracking$data$bodycentre2$y)),
    vecNose2Nose1 = cbind((Tracking$data$nose2$x - Tracking$data$nose1$x), (Tracking$data$nose2$y - Tracking$data$nose1$y)),
  )
  
    # Calculate angles for each row in dfTracking
    for (i in 1:nrow(dfTracking)) {
        vecNose1Bodycentre1 <- c(dfTracking$vecNose1Bodycentre1.1[i], dfTracking$vecNose1Bodycentre1.2[i])
        vecNose1Nose2 <- c(dfTracking$vecNose1Nose2.1[i], dfTracking$vecNose1Nose2.2[i])
        vecNose2Bodycentre2 <- c(dfTracking$vecNose2Bodycentre2.1[i], dfTracking$vecNose2Bodycentre2.2[i])
        vecNose2Nose1 <- c(dfTracking$vecNose2Nose1.1[i], dfTracking$vecNose2Nose1.2[i])
        
        angleDegNose1 <- calcAngle(vecNose1Bodycentre1, vecNose1Nose2)
        angleDegNose2 <- calcAngle(vecNose2Bodycentre2, vecNose2Nose1)
        }

    # Initialize vectors to store angles
    angleDegNose1 <- numeric(length = nrow(dfTracking))
    angleDegNose2 <- numeric(length = nrow(dfTracking))

    # Loop through each row of dfTracking
    for (i in 1:nrow(dfTracking)) {
        vecNose1Bodycentre1 <- c(dfTracking$vecNose1Bodycentre1.1[i], dfTracking$vecNose1Bodycentre1.2[i])
        vecNose1Nose2 <- c(dfTracking$vecNose1Nose2.1[i], dfTracking$vecNose1Nose2.2[i])
        vecNose2Bodycentre2 <- c(dfTracking$vecNose2Bodycentre2.1[i], dfTracking$vecNose2Bodycentre2.2[i])
        vecNose2Nose1 <- c(dfTracking$vecNose2Nose1.1[i], dfTracking$vecNose2Nose1.2[i])

        angleDegNose1[i] <- calcAngle(vecNose1Bodycentre1, vecNose1Nose2)
        angleDegNose2[i] <- calcAngle(vecNose2Bodycentre2, vecNose2Nose1)
    }

  # calculate interaction of animal 1 with animal 2 nose
    isInNoseToNoseContactanimal1 <- GetDistances(Tracking, "nose1", "nose2") <= 2 & abs(angleDegNose1) >= 90 & abs(angleDegNose1) <= 270
    contactNose1 <- sum(isInNoseToNoseContact) * totalTime / length(isInNoseToNoseContact)
    isInNoseToNoseContactanimal2 <- GetDistances(Tracking, "nose2", "nose1") <= 2 & abs(angleDegNose2) >= 90 & abs(angleDegNose2) <= 270
    contactNose2 <- sum(isInNoseToNoseContact) * totalTime / length(isInNoseToNoseContact)

    # Read in bodyparts from names of Tracking$data names also, animal 1 has a 1 after the bodypart name, animal 2 a 2
# Calculate angle of vecNose1Bodycentre1 and every other body part
for (bodypart in bodyparts) {
    # Create new variables in the data frame to store the angle values
    dfTracking[[paste0("angleDegNose1", bodypart)]] <- numeric(length = nrow(dfTracking))
    dfTracking[[paste0("angleDegNose2", bodypart)]] <- numeric(length = nrow(dfTracking))
    
    # Loop through each row of dfTracking
    for (i in 1:nrow(dfTracking)) {
        # Extract vectors a and c from vectorNoseObjL and vector matrices
        vecNose1Bodycentre1 <- c(dfTracking$vecNose1Bodycentre1.1[i], dfTracking$vecNose1Bodycentre1.2[i])
        vecNose1Bodypart <- c(dfTracking[[paste0("vecNose1", bodypart)]][i, 1], dfTracking[[paste0("vecNose1", bodypart)]][i, 2])
        vecNose2Bodycentre2 <- c(dfTracking$vecNose2Bodycentre2.1[i], dfTracking$vecNose2Bodycentre2.2[i])
        vecNose2Bodypart <- c(dfTracking[[paste0("vecNose2", bodypart)]][i, 1], dfTracking[[paste0("vecNose2", bodypart)]][i, 2])
        
        # Calculate angles and store in the vectors
        dfTracking[[paste0("angleDegNose1", bodypart)]][i] <- calcAngle(vecNose1Bodycentre1, vecNose1Bodypart)
        dfTracking[[paste0("angleDegNose2", bodypart)]][i] <- calcAngle(vecNose2Bodycentre2, vecNose2Bodypart)
    }
}


    bodyparts <- names(Tracking$data)
    # exclude bodyparts with the names tl, tr, br, bl
    bodyparts <- bodyparts[!bodyparts %in% c("tl", "tr", "br", "bl")]

    # calculate angle of vecNose1Bodycentre1 and every other bodypart
    for (bodypart in bodyparts) {
        # Create new variables in the data frame to store the angle values
        assign(paste0("angleDegNose1", bodypart), numeric(length = nrow(dfTracking)))
        assign(paste0("angleDegNose2", bodypart), numeric(length = nrow(dfTracking)))
        # Loop through each row of dfTracking
        for (i in 1:nrow(dfTracking)) {
            # Extract vectors a and c from vectorNoseObjL and vector matrices
            vecNose1Bodycentre1 <- c(dfTracking$vecNose1Bodycentre1.1[i], dfTracking$vecNose1Bodycentre1.2[i])
            vecNose1Bodypart <- c(dfTracking$vecNose1Bodypart.1[i], dfTracking$vecNose1Bodypart.2[i])
            vecNose2Bodycentre2 <- c(dfTracking$vecNose2Bodycentre2.1[i], dfTracking$vecNose2Bodycentre2.2[i])
            vecNose2Bodypart <- c(dfTracking$vecNose2Bodypart.1[i], dfTracking$vecNose2Bodypart.2[i])

            # Calculate angles and store in the vectors
            assign(paste0("angleDegNose1", bodypart), calcAngle(vecNose1Bodycentre1, vecNose1Bodypart))
            assign(paste0("angleDegNose2", bodypart), calcAngle(vecNose2Bodycentre2, vecNose2Bodypart))
        }
    }

    # calculate distances of each bodypart of animal 1 to each bodypart of animal2 bodyparts hve the same names for both animals, only differentiated by either a 1 or a 2 after the bodypart eg nose1 and nose2
    for (bodypart in bodyparts) {
        # Create new variables in the data frame to store the distance values
        assign(paste0("distance", bodypart), numeric(length = nrow(dfTracking)))
        # Loop through each row of dfTracking
        for (i in 1:nrow(dfTracking)) {
            # Extract distances between bodyparts
            bodypart[i]

            # Calculate angles and store in the vectors
            assign(paste0("distance", bodypart), sqrt(sum((vecNose1Bodypart - vecNose2Bodypart)^2)))
        }
    }


  
  # calculate unisided interaction with noses between animal 1 and 2
    BP1toBP2Contact <- GetDistances(Tracking, "bp1", "bp2") <= 2 & abs(angleDegbp1) >= 90 & abs(angleDegbp1)
    onesidedContactBP1toBP2 <- sum(BP1toBP2Contact) * totalTime / length(BP1toBP2Contact)


    onesidedContactNose1 <- sum(isInNoseToNoseContact) * totalTime / length(isInNoseToNoseContact)
    isInNoseToNoseContactanimal2 <- GetDistances(Tracking, "nose2", "nose1") <= 2 & abs(angleDegNose2) >= 90 & abs(angleDegNose2) <= 270 & abs(angleDegNose1) <= 90 & abs(angleDegNose1) >= 270
    onesidedContactNose2 <- sum(isInNoseToNoseContact) * totalTime / length(isInNoseToNoseContact)


  # calculate unisided interaction with noses between animal 1 and 2
    isInNoseToNoseContactanimal1 <- GetDistances(Tracking, "nose1", "nose2") <= 2 & abs(angleDegNose1) >= 90 & abs(angleDegNose1) <= 270 & abs(angleDegNose2) <= 90 & abs(angleDegNose2) >= 270
    onesidedContactNose1 <- sum(isInNoseToNoseContact) * totalTime / length(isInNoseToNoseContact)
    isInNoseToNoseContactanimal2 <- GetDistances(Tracking, "nose2", "nose1") <= 2 & abs(angleDegNose2) >= 90 & abs(angleDegNose2) <= 270 & abs(angleDegNose1) <= 90 & abs(angleDegNose1) >= 270
    onesidedContactNose2 <- sum(isInNoseToNoseContact) * totalTime / length(isInNoseToNoseContact)

  # calculate reciprocal interaction between noses netween animal 1 and 2
    isInReciprocalNoseToNoseContactanimal1 <- GetDistances(Tracking, "nose1", "nose2") <= 2 & abs(angleDegNose1) >= 90 & abs(angleDegNose1) <= 270 & abs(angleDegNose2) >= 90 & abs(angleDegNose2) <= 270
    reciprocalContactNose1 <- sum(isInReciprocalNoseToNoseContact) * totalTime / length(isInReciprocalNoseToNoseContact)
    
  
  # Check if code exists in novelLoc file
  if (code %in% novelLoc$Code) {
    totalTime <- 10 * 60
    if (novelLoc[novelLoc$Code == code, "NovelLoc"] == "R") {
      isInZoneLeft <- point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objLSquare$x, objLSquare$y) & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(angleDegObjL) >= 70 & abs(angleDegObjL) <= 290
      isInZoneRight <- GetDistances(Tracking, "objR", "nose") <= 4 & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(angleDegObjR) >= 70 & abs(angleDegObjR) <= 290
    } else {
      isInZoneLeft <- GetDistances(Tracking, "objL", "nose") <= 4 & GetDistances(Tracking, "objL", "bodycentre") > 1 & abs(angleDegObjL) >= 70 & abs(angleDegObjL) <= 290
      isInZoneRight <- point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objRSquare$x, objRSquare$y) & GetDistances(Tracking, "objR", "bodycentre") > 1 & abs(angleDegObjR) >= 70 & abs(angleDegObjR) <= 290
    }
    contactLeft <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
    contactRight <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
  } else {
    contactLeft <- 0
    contactRight <- 0
  }

  # calculate interaction time of nose to nose contact between two mice
  isInNoseToNoseContact <- GetDistances(Tracking, "nose1", "nose2") <= 2 

  
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
    isInZoneVisitR <- cumsum((point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objRSquare$x, objRSquare$y) & c(0, diff(point.in.polygon(Tracking$data$nose$x, Tracking$data$nose$y, objRSquare$x, objRSquare$y)))) == 1
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
write.csv(dfCombined, file.path(outputDir, "combined_output_test.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")