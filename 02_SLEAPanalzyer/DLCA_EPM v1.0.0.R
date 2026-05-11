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

batches <- c("B1", "B2", "B3", "B4", "B5", "B6")

for (batch in batches) {
  inputFolder <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/EPM/SLEAP/formatted/")
  outputFolder <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/EPM/SLEAP/output/")
  overviewPlotFolder <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/", batch, "/EPM/SLEAP/OverviewPlot/")

  # Create the output folder if it doesn't exist
  if (!file.exists(outputFolder))
    dir.create(outputFolder)

  # Create the overview plot folder if it doesn't exist
  if (!file.exists(overviewPlotFolder))
    dir.create(overviewPlotFolder)

  files <- list.files(inputFolder)

  pipeline <- function(path){
    tracking <- ReadDLCDataFromCSV(file = path, fps = 30)
    tracking <- CalibrateTrackingData(tracking, method = "distance", in.metric = 60, points = c("tl", "br"))
    zoneInfo <- read.table("S:/Lab_Member/Tobi/Experiments/DLCAnalyzer/ArenaConfig/EPM_zoneinfo.csv", sep = ";", header = TRUE)
    tracking <- AddZones(tracking, zoneInfo)
    tracking <- EPMAnalysis(tracking, movementCutoff = 5, integrationPeriod = 5, points = "bodycentre", noseDips = TRUE)
    return(tracking)
  }

  trackingAll <- RunPipeline(files, inputFolder, FUN = pipeline)

  report <- MultiFileReport(trackingAll)

  # Save report as a .csv file
  reportFile <- paste0(outputFolder, "Report.csv")
  write.csv(report, file = reportFile, row.names = FALSE)
  
  # Notify completion
  cat("Batch", batch, "processing completed.\n")

  # Save overview plots as .tiff files
  library(ggplot2)
  for (file in files) {
    inputFile <- paste0(inputFolder, file)
    outputFile <- paste0(overviewPlotFolder, gsub(".csv", ".tiff", file))

    tracking <- pipeline(inputFile)
    plot <- OverviewPlot(tracking, "bodycentre")
    ggsave(filename = outputFile, plot = plot, device = "tiff", width = 8, height = 10)
    
    # Notify completion
    cat("Overview plot for", file, "generated.\n")
  }
}
