# Check if packages are installed, install them if needed, and load them
requiredPackages <- c("tensorflow", "reticulate", "keras", "sp", "imputeTS", "ggplot2", "ggmap", "data.table", "cowplot", "corrplot", "zoo")

for (package in requiredPackages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Set input and output directories
inputDir <- "C:/Users/topohl/Documents/GitHub/sleap/docs/notebooks"
outputDir <- "C:/Users/topohl/Documents/GitHub/sleap/docs/notebooks/output"

# Get a list of CSV files in the input directory
fileList <- list.files(path = inputDir, pattern = "*.csv")

# create new list dfList
dfList <- list()

# Function to replace NAs in specified columns of a data frame with the last known values
replaceNAs <- function(tracking, columns) {
  for (col in columns) {
      tracking$data[[col]] <- zoo::na.locf(tracking$data[[col]])
  }
  return(tracking)
}

# Define the names of the columns to replace NAs
replaceNABodyparts <- c("nose", "bodycentre")

# Loop through each file in the input directory
for (file in fileList) {
  
  # Read in tracking data and get names of variables in the data frame
  inputFile <- file.path(inputDir, file)
  tracking <- ReadDLCDataFromCSV(file = inputFile, fps = 30)

  # Extract the input file name without the extension
  inputFileName <- sub(".csv$", "", basename(inputFile))

  # Call the replaceNAs function with the defined column names
  tracking <- replaceNAs(tracking, replaceNABodyparts)

  # Calibrate the tracking data, add zones, and plot the zones
  tracking <- CalibrateTrackingData(tracking, method = "area",in.metric = 44*24, points = c("tl","tr","br","bl"))
  tracking <- AddOFTZones(tracking, scale_center = 0.5,scale_periphery  = 0.8 ,scale_corners = 0.4, points = c("tl","tr","br","bl"))
  PlotZones(tracking)
  PlotPointData(tracking, points = c("nose"))

  # Calculate time_Left and time_Right direct
  totalTime <- 10 * 60
  isInZoneLeft <- GetDistances(tracking, "socl", "nose") < 6
  isInZoneRight <- GetDistances(tracking, "socr", "nose") < 6
  contactLeft <- sum(isInZoneLeft) * totalTime / length(isInZoneLeft)
  contactRight <- sum(isInZoneRight) * totalTime / length(isInZoneRight)
  
  # Calculate time with nose between 6-10 cm of either socl or socr
  isInProxLeft <- GetDistances(tracking, "socl", "nose") >= 6 & GetDistances(tracking, "socl", "nose") <= 10
  isInProxRight <- GetDistances(tracking, "socr", "nose") >= 6 & GetDistances(tracking, "socr", "nose") <= 10
  proxLeft <- sum(isInProxLeft) * totalTime / length(isInProxLeft)
  proxRight <- sum(isInProxRight) * totalTime / length(isInProxRight)

  # Create a data frame with the time_Left, time_Right, time_nose_Left, and time_nose_Right values
  df <- data.frame(file = inputFileName, contactLeft = contactLeft, contactRight = contactRight, 
                 proxLeft = proxLeft, proxRight = proxRight)

  # Construct the output file name
  outputFile <- paste0(inputFileName, "_output.csv")

  # Write the data frame to a csv file with the constructed file name
  write.csv(df, outputFile, row.names = FALSE)

  # Add the current data frame to the list
  dfList[[length(dfList) + 1]] <- df

  # Print progress message
  cat(sprintf("Processed file %s\n", file))
}

# Combine all data frames into a single data frame
dfCombined <- do.call(rbind, dfList)

# Write the combined data frame to a csv file
write.csv(dfCombined, file.path(outputDir, "combined_output.csv"), row.names = FALSE)

# Print "done" message
cat("done\n")
