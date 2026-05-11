# install and load packages
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  tensorflow, reticulate, keras, sp, imputeTS, ggplot2, ggmap,
  data.table, cowplot, corrplot, zoo, stringr, tidyverse
)

# Define length of video in minutes and translation to seconds
video_length <- 2
total_time <- videoLength * 60

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/SLEAPanalyzer/SLEAPanalzyer/")
source('DLCAnalyzer_Functions_final.R')

# Set input and output directories
inputDir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC"
outputDir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC/output"
plotDir <- file.path(outputDir, "plots")
dir.create(plotDir, showWarnings = FALSE)

# Get a list of CSV files in the input directory
fileList <- list.files(path = inputDir, pattern = "*.csv")

# Function to calculate angle between two vectors
calcAngle <- function(vec1, vec2) {
    dotProd <- sum(vec1 * vec2) 
    mag1 <- sqrt(sum(vec1^2))
    mag2 <- sqrt(sum(vec2^2))
    angleRad <- acos(dotProd / (mag1 * mag2))
    return(angleRad * 180 / pi)
}

# Function to handle interpolation with specified behavior
interpolate_with_leading_trailing <- function(x) { 
    if (is.na(x[1])) {
        first_non_na <- which(!is.na(x))[1]
        if (!is.na(first_non_na)) {
            x[1:(first_non_na - 1)] <- x[first_non_na]
        }
    }
    if (is.na(x[length(x)])) {
        last_non_na <- tail(which(!is.na(x)), 1)
        if (!is.na(last_non_na)) {
            x[(last_non_na + 1):length(x)] <- x[last_non_na]
        }
    }
    x <- na.approx(x, na.rm = FALSE)
    return(x)
}

# Process each file
dfList <- lapply(fileList, function(file) {
    inputFile <- file.path(inputDir, file)
    Tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)
    dfTracking <- data.frame()

    length <- length(Tracking$data$nose_1$x)
    seconds_length <- length / 30
    minutes_length <- seconds_length / 60

    bodyparts <- c("nose_1", "leftEar_1", "rightEar_1", "bodycentre_1", "leftSide_1", "rightSide_1", "tailBase_1", "tailEnd_1", "nose_2", "leftEar_2", "rightEar_2", "bodycentre_2", "leftSide_2", "rightSide_2", "tailBase_2", "tailEnd_2")
    bodyparts1 <- bodyparts[1:8]
    bodyparts2 <- bodyparts[9:16]

    # Interpolate missing values
    for (bodypart in bodyparts) {
        Tracking$data[[bodypart]]$x <- interpolate_with_leading_trailing(Tracking$data[[bodypart]]$x)
        Tracking$data[[bodypart]]$y <- interpolate_with_leading_trailing(Tracking$data[[bodypart]]$y)
    }

    # Perform OFT analysis
    Tracking <- OFTAnalysis(Tracking, points = c("bodycentre_1", "bodycentre_2"), movement_cutoff = 5, integration_period = 5)

    # Create data frame for calculating angles
    dfTracking <- data.frame(
        vecNose1Bodycentre1 = cbind((Tracking$data$nose_1$x - Tracking$data$bodycentre_1$x), (Tracking$data$nose_1$y - Tracking$data$bodycentre_1$y)),
        vecNose2Bodycentre2 = cbind((Tracking$data$nose_2$x - Tracking$data$bodycentre_2$x), (Tracking$data$nose_2$y - Tracking$data$bodycentre_2$y))
    )

    # Calculate vectors between nose and body parts
    for (bodypart in bodyparts2) {
        dfTracking[[paste0("vecNose1", bodypart)]] <- cbind((Tracking$data$nose_1$x - Tracking$data[[bodypart]]$x), (Tracking$data$nose_1$y - Tracking$data[[bodypart]]$y))
    }
    for (bodypart in bodyparts1) {
        dfTracking[[paste0("vecNose2", bodypart)]] <- cbind((Tracking$data$nose_2$x - Tracking$data[[bodypart]]$x), (Tracking$data$nose_2$y - Tracking$data[[bodypart]]$y))
    }

    # Calculate angles
    for (bodypart in bodyparts2) {
        dfTracking[[paste0("angleDegNose1", bodypart)]] <- apply(dfTracking, 1, function(row) {
            vecNose1Bodycentre1 <- c(row["vecNose1Bodycentre1.1"], row["vecNose1Bodycentre1.2"])
            vecNose1Bodypart2 <- c(row[[paste0("vecNose1", bodypart)]][1], row[[paste0("vecNose1", bodypart)]][2])
            calcAngle(vecNose1Bodycentre1, vecNose1Bodypart2)
        })
    }
    for (bodypart in bodyparts1) {
        dfTracking[[paste0("angleDegNose2", bodypart)]] <- apply(dfTracking, 1, function(row) {
            vecNose2Bodycentre2 <- c(row["vecNose2Bodycentre2.1"], row["vecNose2Bodycentre2.2"])
            vecNose2Bodypart1 <- c(row[[paste0("vecNose2", bodypart)]][1], row[[paste0("vecNose2", bodypart)]][2])
            calcAngle(vecNose2Bodycentre2, vecNose2Bodypart1)
        })
    }

    # Calculate interactions and metrics
    contactNose1toNose2 <- sum(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270) * seconds_length / length(Tracking$data$nose_1$x)
    contactNose1toBodycentre2 <- sum(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270) * seconds_length / length(Tracking$data$nose_1$x)
    contactNose1toTailBase2 <- sum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270) * seconds_length / length(Tracking$data$nose_1$x)
    contactSideBySideAnimal1 <- sum(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80) * seconds_length / length(Tracking$data$nose_1$x)
    contactSideBySideReverseAnimal1 <- sum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80) * seconds_length / length(Tracking$data$nose_1$x)

    contactNose2toNose1 <- sum(GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270) * seconds_length / length(Tracking$data$nose_2$x)
    contactNose2toBodycentre1 <- sum(GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose2bodycentre_1) >= 90 & abs(dfTracking$angleDegNose2bodycentre_1) <= 270) * seconds_length / length(Tracking$data$nose_2$x)
    contactNose2toTailBase1 <- sum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose2tailBase_1) >= 90 & abs(dfTracking$angleDegNose2tailBase_1) <= 270) * seconds_length / length(Tracking$data$nose_2$x)
    contactSideBySideAnimal2 <- sum(GetDistances(Tracking, "nose_2", "nose_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "tailBase_1") <= 80) * seconds_length / length(Tracking$data$nose_2$x)
    contactSideBySideReverseAnimal2 <- sum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "nose_1") <= 80) * seconds_length / length(Tracking$data$nose_2$x)

    latencyNose1toNose2 <- min(which(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270)) * 1/30
    latencyNose1toBodycentre2 <- min(which(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270)) * 1/30
    latencyNose1toTailBase2 <- min(which(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270)) * 1/30
    latencySideBySideAnimal1 <- min(which(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80)) * 1/30
    latencySideBySideReverseAnimal1 <- min(which(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80)) * 1/30

    latencyNose2toNose1 <- min(which(GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270)) * 1/30
    latencyNose2toBodycentre1 <- min(which(GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose2bodycentre_1) >= 90 & abs(dfTracking$angleDegNose2bodycentre_1) <= 270)) * 1/30
    latencyNose2toTailBase1 <- min(which(GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose2tailBase_1) >= 90 & abs(dfTracking$angleDegNose2tailBase_1) <= 270)) * 1/30
    latencySideBySideAnimal2 <- min(which(GetDistances(Tracking, "nose_2", "nose_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "tailBase_1") <= 80)) * 1/30
    latencySideBySideReverseAnimal2 <- min(which(GetDistances(Tracking, "nose_2", "tailBase_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "nose_1") <= 80)) * 1/30

    frequencyNose1toNose2 <- sum(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270)
    frequencyNose1toBodycentre2 <- sum(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270)
    frequencyNose1toTailBase2 <- sum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270)
    frequencySideBySideAnimal1 <- sum(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80)
    frequencySideBySideReverseAnimal1 <- sum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80)

    frequencyNose2toNose1 <- sum(GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270)
    frequencyNose2toBodycentre1 <- sum(GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose2bodycentre_1) >= 90 & abs(dfTracking$angleDegNose2bodycentre_1) <= 270)
    frequencyNose2toTailBase1 <- sum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose2tailBase_1) >= 90 & abs(dfTracking$angleDegNose2tailBase_1) <= 270)
    frequencySideBySideAnimal2 <- sum(GetDistances(Tracking, "nose_2", "nose_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "tailBase_1") <= 80)
    frequencySideBySideReverseAnimal2 <- sum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "nose_1") <= 80)

    proxAnimal1 <- sum(GetDistances(Tracking, "nose_1", "nose_2") > 4 & GetDistances(Tracking, "nose_1", "nose_2") <= 8) * seconds_length / length(Tracking$data$nose_1$x)
    proxAnimal1Angle <- sum(GetDistances(Tracking, "nose_1", "nose_2") > 4 & GetDistances(Tracking, "nose_1", "nose_2") <= 8 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270) * seconds_length / length(Tracking$data$nose_1$x)

    proxAnimal2 <- sum(GetDistances(Tracking, "nose_2", "nose_1") > 4 & GetDistances(Tracking, "nose_2", "nose_1") <= 8) * seconds_length / length(Tracking$data$nose_2$x)
    proxAnimal2Angle <- sum(GetDistances(Tracking, "nose_2", "nose_1") > 4 & GetDistances(Tracking, "nose_2", "nose_1") <= 8 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270) * seconds_length / length(Tracking$data$nose_2$x)

    # calculate following metrics based on nose to tailbase distance over frames (if one animal center is moving and the nose of the other animal is in the same position < 30 px maximum 1 second or 30 frames later)
    followingNose1toTailBase2 <- sum(sapply(1:(length(Tracking$data$nose_1$x) - 30), function(i) {
        if (any(GetDistances(Tracking, "nose_1", "tailBase_2")[max(1, i-30):i] <= 30)) {
            return(any(GetDistances(Tracking, "bodycentre_2", "bodycentre_2")[i:(i + 30)] > 5))
        } else {
            return(FALSE)
        }
    })) * seconds_length / length(Tracking$data$nose_1$x)

    followingNose2toTailBase1 <- sum(sapply(1:(length(Tracking$data$nose_2$x) - 30), function(i) {
        if (GetDistances(Tracking, "nose_2", "tailBase_1")[max(1, i-30):i] <= 30) {
            return(any(GetDistances(Tracking, "bodycentre_1", "bodycentre_1")[i:(i + 30)] > 5))
        } else {
            return(FALSE)
        }
    })) * seconds_length / length(Tracking$data$nose_2$x)

    df <- data.frame(
        file = sub(".csv$", "", basename(inputFile)),
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
        followingNose1toTailBase2 = followingNose1toTailBase2,
        followingNose2toTailBase1 = followingNose2toTailBase1,
        proxAnimal1 = proxAnimal1,
        proxAnimal1Angle = proxAnimal1Angle,
        contactNose2toNose1 = contactNose2toNose1,
        contactNose2toBodycentre1 = contactNose2toBodycentre1,
        contactNose2toTailBase1 = contactNose2toTailBase1,
        contactSideBySideAnimal2 = contactSideBySideAnimal2,
        contactSideBySideReverseAnimal2 = contactSideBySideReverseAnimal2,
        latencyNose2toNose1 = latencyNose2toNose1,
        latencyNose2toBodycentre1 = latencyNose2toBodycentre1,
        latencyNose2toTailBase1 = latencyNose2toTailBase1,
        latencySideBySideAnimal2 = latencySideBySideAnimal2,
        latencySideBySideReverseAnimal2 = latencySideBySideReverseAnimal2,
        frequencyNose2toNose1 = frequencyNose2toNose1,
        frequencyNose2toBodycentre1 = frequencyNose2toBodycentre1,
        frequencyNose2toTailBase1 = frequencyNose2toTailBase1,
        frequencySideBySideAnimal2 = frequencySideBySideAnimal2,
        frequencySideBySideReverseAnimal2 = frequencySideBySideReverseAnimal2,
        proxAnimal2 = proxAnimal2,
        proxAnimal2Angle = proxAnimal2Angle,
        distance = Tracking$Report$bodycentre.raw.distance,
        stationary = Tracking$Report$bodycentre.time.stationary,
        speedMoving = Tracking$Report$bodycentre.speed.moving,
        speedRaw = Tracking$Report$bodycentre.raw.speed
    )

    outputFile <- file.path(outputDir, paste0(sub(".csv$", "", basename(inputFile)), "_output.csv"))
    write.csv(df, outputFile, row.names = FALSE)

    plots <- PlotDensityPaths(Tracking, points = c("bodycentre"))
    plotFile <- file.path(plotDir, paste0(sub(".csv$", "", basename(inputFile)), "_DensityPath.png"))
    ggsave(plotFile, plot = plots$bodycentre, width = 7, height = 6)

    cat(sprintf("Processed file %s\n", file))
    return(df)
})

dfCombined <- do.call(rbind, dfList)
write.csv(dfCombined, file.path(outputDir, "combined_output_test.csv"), row.names = FALSE)
cat("done\n")
