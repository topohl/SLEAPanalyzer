# Load required packages
required_packages <- c("tensorflow", "reticulate", "keras", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "zoo", "stringr", "tidyverse")

install_and_load_packages <- function(packages) {
    for (package in packages) {
        if (!requireNamespace(package, quietly = TRUE)) {
            install.packages(package)
        }
        library(package, character.only = TRUE)
    }
}

install_and_load_packages(required_packages)

# Define constants
videoLength <- 2
totalTime <- videoLength * 60
inputDir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC"
outputDir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Rosalba/SocInteraction/DLC/output"
plotDir <- file.path(outputDir, "plots")
dir.create(plotDir, showWarnings = FALSE)

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/SLEAPanalyzer/SLEAPanalzyer/")
source('DLCAnalyzer_Functions_final.R')

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

# Function to process each file
process_file <- function(file) {
    inputFile <- file.path(inputDir, file)
    Tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)
    length <- length(Tracking$data$nose_1$x)
    seconds_length <- length / 30
    minutes_length <- seconds_length / 60

    bodyparts <- c("nose_1", "leftEar_1", "rightEar_1", "bodycentre_1", "leftSide_1", "rightSide_1", "tailBase_1", "tailEnd_1", "nose_2", "leftEar_2", "rightEar_2", "bodycentre_2", "leftSide_2", "rightSide_2", "tailBase_2", "tailEnd_2")

    for (bodypart in bodyparts1) {
        Tracking$data[[bodypart]]$x <- na.approx(Tracking$data[[bodypart]]$x, method = "linear", na.rm = FALSE)
        Tracking$data[[bodypart]]$y <- na.approx(Tracking$data[[bodypart]]$y, method = "linear", na.rm = FALSE)
    }

    for (bodypart in bodyparts2) {
        Tracking$data[[bodypart]]$x <- na.approx(Tracking$data[[bodypart]]$x, method = "linear", na.rm = FALSE)
        Tracking$data[[bodypart]]$y <- na.approx(Tracking$data[[bodypart]]$y, method = "linear", na.rm = FALSE)
    }

    Tracking <- OFTAnalysis(Tracking, points = "bodycentre", movement_cutoff = 5, integration_period = 5)

    dfTracking <- data.frame(
        vecNose1Bodycentre1 = cbind((Tracking$data$nose_1$x - Tracking$data$bodycentre_1$x), (Tracking$data$nose_1$y - Tracking$data$bodycentre_1$y)),
        vecNose2Bodycentre2 = cbind((Tracking$data$nose_2$x - Tracking$data$bodycentre_2$x), (Tracking$data$nose_2$y - Tracking$data$bodycentre_2$y))
    )

    for (bodypart in bodyparts2) {
        dfTracking[[paste0("vecNose1", bodypart)]] <- cbind((Tracking$data$nose_1$x - Tracking$data[[bodypart]]$x), (Tracking$data$nose_1$y - Tracking$data[[bodypart]]$y))
        dfTracking[[paste0("vecNose2", bodypart)]] <- cbind((Tracking$data$nose_2$x - Tracking$data[[bodypart]]$x), (Tracking$data$nose_2$y - Tracking$data[[bodypart]]$y))
    }

    for (bodypart in bodyparts2) {
        dfTracking[[paste0("angleDegNose1", bodypart)]] <- numeric(length = nrow(dfTracking))
        dfTracking[[paste0("angleDegNose2", bodypart)]] <- numeric(length = nrow(dfTracking))
        for (i in 1:nrow(dfTracking)) {
            vecNose1Bodycentre1 <- c(dfTracking$vecNose1Bodycentre1.1[i], dfTracking$vecNose1Bodycentre1.2[i])
            vecNose1Bodypart <- c(dfTracking[[paste0("vecNose1", bodypart)]][i, 1], dfTracking[[paste0("vecNose1", bodypart)]][i, 2])
            vecNose2Bodycentre2 <- c(dfTracking$vecNose2Bodycentre2.1[i], dfTracking$vecNose2Bodycentre2.2[i])
            vecNose2Bodypart <- c(dfTracking[[paste0("vecNose2", bodypart)]][i, 1], dfTracking[[paste0("vecNose2", bodypart)]][i, 2])

            dfTracking[[paste0("angleDegNose1", bodypart)]][i] <- calcAngle(vecNose1Bodycentre1, vecNose1Bodypart)
            dfTracking[[paste0("angleDegNose2", bodypart)]][i] <- calcAngle(vecNose2Bodycentre2, vecNose2Bodypart)
        }
    }

    # Calculate interactions and other metrics
    calculate_metrics <- function(Tracking, dfTracking, bodyparts1, bodyparts2, seconds_length) {
        metrics <- list()
        metrics$contactNose1toNose2 <- sum(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270) * seconds_length / length(dfTracking$angleDegNose1nose_2)
        metrics$contactNose1toBodycentre2 <- sum(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270) * seconds_length / length(dfTracking$angleDegNose1bodycentre_2)
        metrics$contactNose1toTailBase2 <- sum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270) * seconds_length / length(dfTracking$angleDegNose1tailBase_2)
        metrics$contactSideBySideAnimal1 <- sum(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80) * seconds_length / length(dfTracking$angleDegNose1nose_2)
        metrics$contactSideBySideReverseAnimal1 <- sum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80) * seconds_length / length(dfTracking$angleDegNose1tailBase_2)
        metrics$contactNose2toNose1 <- sum(GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270) * seconds_length / length(dfTracking$angleDegNose2nose_1)
        metrics$contactNose2toBodycentre1 <- sum(GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose2bodycentre_1) >= 90 & abs(dfTracking$angleDegNose2bodycentre_1) <= 270) * seconds_length / length(dfTracking$angleDegNose2bodycentre_1)
        metrics$contactNose2toTailBase1 <- sum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose2tailBase_1) >= 90 & abs(dfTracking$angleDegNose2tailBase_1) <= 270) * seconds_length / length(dfTracking$angleDegNose2tailBase_1)
        metrics$contactSideBySideAnimal2 <- sum(GetDistances(Tracking, "nose_2", "nose_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "tailBase_1") <= 80) * seconds_length / length(dfTracking$angleDegNose2nose_1)
        metrics$contactSideBySideReverseAnimal2 <- sum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "nose_1") <= 80) * seconds_length / length(dfTracking$angleDegNose2tailBase_1)
        metrics$latencyNose1toNose2 <- min(which(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270)) * 1/30
        metrics$latencyNose1toBodycentre2 <- min(which(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270)) * 1/30
        metrics$latencyNose1toTailBase2 <- min(which(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270)) * 1/30
        metrics$latencySideBySideAnimal1 <- min(which(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80)) * 1/30
        metrics$latencySideBySideReverseAnimal1 <- min(which(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80)) * 1/30
        metrics$latencyNose2toNose1 <- min(which(GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270)) * 1/30
        metrics$latencyNose2toBodycentre1 <- min(which(GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose2bodycentre_1) >= 90 & abs(dfTracking$angleDegNose2bodycentre_1) <= 270)) * 1/30
        metrics$latencyNose2toTailBase1 <- min(which(GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose2tailBase_1) >= 90 & abs(dfTracking$angleDegNose2tailBase_1) <= 270)) * 1/30
        metrics$latencySideBySideAnimal2 <- min(which(GetDistances(Tracking, "nose_2", "nose_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "tailBase_1") <= 80)) * 1/30
        metrics$latencySideBySideReverseAnimal2 <- min(which(GetDistances(Tracking, "nose_2", "tailBase_1") <= 80 & GetDistances(Tracking, "bodycentre_2", "bodycentre_1") <= 80 & GetDistances(Tracking, "tailBase_2", "nose_1") <= 80)) * 1/30
        metrics$frequencyNose1toNose2 <- sum(cumsum(GetDistances(Tracking, "nose_1", "nose_2") <= 30 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270) == 1)
        metrics$frequencyNose1toBodycentre2 <- sum(cumsum(GetDistances(Tracking, "nose_1", "bodycentre_2") <= 30 & abs(dfTracking$angleDegNose1bodycentre_2) >= 90 & abs(dfTracking$angleDegNose1bodycentre_2) <= 270) == 1)
        metrics$frequencyNose1toTailBase2 <- sum(cumsum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 30 & abs(dfTracking$angleDegNose1tailBase_2) >= 90 & abs(dfTracking$angleDegNose1tailBase_2) <= 270) == 1)
        metrics$frequencySideBySideAnimal1 <- sum(cumsum(GetDistances(Tracking, "nose_1", "nose_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "tailBase_2") <= 80) == 1)
        metrics$frequencySideBySideReverseAnimal1 <- sum(cumsum(GetDistances(Tracking, "nose_1", "tailBase_2") <= 80 & GetDistances(Tracking, "bodycentre_1", "bodycentre_2") <= 80 & GetDistances(Tracking, "tailBase_1", "nose_2") <= 80) == 1)
        metrics$frequencyNose2toNose1 <- sum(cumsum(GetDistances(Tracking, "nose_2", "nose_1") <= 30 & abs(dfTracking$angleDegNose2nose_1) >= 90 & abs(dfTracking$angleDegNose2nose_1) <= 270) == 1)
        metrics$frequencyNose2toBodycentre1 <- sum(cumsum(GetDistances(Tracking, "nose_2", "bodycentre_1") <= 30 & abs(dfTracking$angleDegNose2bodycentre_1) >= 90 & abs(dfTracking$angleDegNose2bodycentre_1) <= 270) == 1)
        metrics$frequencyNose2toTailBase1 <- sum(cumsum(GetDistances(Tracking, "nose_2", "tailBase_1") <= 30 & abs(dfTracking$angleDegNose2tailBase_1) >= 90 & abs(dfTracking$angleDegNose2tailBase_1) <= 270) == 1)
        metrics$proxAnimal1 <- sum(GetDistances(Tracking, "nose_1", "nose_2") > 4 & GetDistances(Tracking, "nose_1", "nose_2") <= 8) * seconds_length / length(dfTracking$angleDegNose1nose_2)
        metrics$proxAnimal1Angle <- sum(GetDistances(Tracking, "nose_1", "nose_2") > 4 & GetDistances(Tracking, "nose_1", "nose_2") <= 8 & abs(dfTracking$angleDegNose1nose_2) >= 90 & abs(dfTracking$angleDegNose1nose_2) <= 270) * seconds_length / length(dfTracking$angleDegNose1nose_2)
        metrics$distance <- Tracking$Report$bodycentre.raw.distance
        metrics$stationary <- Tracking$Report$bodycentre.time.stationary
        metrics$speedMoving <- Tracking$Report$bodycentre.speed.moving
        metrics$speedRaw <- Tracking$Report$bodycentre.raw.speed
        metrics$frequencyRear <- sum(cumsum(GetDistances(Tracking, "spine1", "bodycentre") <= 1) > 0 & cumsum(GetDistances(Tracking, "bodycentre", "spine2") <= 1) > 0 & c(0, diff(GetDistances(Tracking, "spine1", "bodycentre") <= 1)) == 1 & c(0, diff(GetDistances(Tracking, "bodycentre", "spine2") <= 1)) == 1)
        return(metrics)
    }

    metrics <- calculate_metrics(Tracking, dfTracking, bodyparts1, bodyparts2, seconds_length)

    df <- data.frame(
        file = sub(".csv$", "", basename(inputFile)),
        contactNose1toNose2 = metrics$contactNose1toNose2,
        contactNose1toBodycentre2 = metrics$contactNose1toBodycentre2,
        contactNose1toTailBase2 = metrics$contactNose1toTailBase2,
        contactSideBySideAnimal1 = metrics$contactSideBySideAnimal1,
        contactSideBySideReverseAnimal1 = metrics$contactSideBySideReverseAnimal1,
        latencyNose1toNose2 = metrics$latencyNose1toNose2,
        latencyNose1toBodycentre2 = metrics$latencyNose1toBodycentre2,
        latencyNose1toTailBase2 = metrics$latencyNose1toTailBase2,
        latencySideBySideAnimal1 = metrics$latencySideBySideAnimal1,
        latencySideBySideReverseAnimal1 = metrics$latencySideBySideReverseAnimal1,
        frequencyNose1toNose2 = metrics$frequencyNose1toNose2,
        frequencyNose1toBodycentre2 = metrics$frequencyNose1toBodycentre2,
        frequencyNose1toTailBase2 = metrics$frequencyNose1toTailBase2,
        frequencySideBySideAnimal1 = metrics$frequencySideBySideAnimal1,
        frequencySideBySideReverseAnimal1 = metrics$frequencySideBySideReverseAnimal1,
        proxAnimal1 = metrics$proxAnimal1,
        proxAnimal1Angle = metrics$proxAnimal1Angle,
        distance = metrics$distance,
        stationary = metrics$stationary,
        speedMoving = metrics$speedMoving,
        speedRaw = metrics$speedRaw,
        frequencyRear = metrics$frequencyRear
    )

    outputFile <- paste0(sub(".csv$", "", basename(inputFile)), "_output.csv")
    write.csv(df, file.path(outputDir, outputFile), row.names = FALSE)

    # Ensure data used for plotting is valid
    Tracking$data <- lapply(Tracking$data, function(df) {
        df <- df[complete.cases(df), ]
        df <- df[is.finite(df$x) & is.finite(df$y), ]
        return(df)
    })
    
    plots <- PlotDensityPaths(Tracking, points = c("bodycentre"))
    plotFile <- file.path(plotDir, paste0(sub(".csv$", "", basename(inputFile)), "_DensityPath.png"))
    ggsave(plotFile, plot = plots$bodycentre, width = 7, height = 6)

    return(df)
}

# Process all files and combine results
dfList <- lapply(fileList, process_file)
dfCombined <- do.call(rbind, dfList)
write.csv(dfCombined, file.path(outputDir, "combined_output_test.csv"), row.names = FALSE)

cat("done\n")
