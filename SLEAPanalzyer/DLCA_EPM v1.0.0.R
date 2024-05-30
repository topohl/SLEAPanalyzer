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

inputFolder <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B5/EPM/SLEAP/formatted/"
outputFolder <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B5/EPM/SLEAP/output/"
overviewplot_folder <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/B5/EPM/SLEAP/OverviewPlot/"

# Create the output folder if it doesn't exist
if (!file.exists(outputFolder))
  dir.create(outputFolder)

# Create the overviewplot folder if it doesn't exist
if (!file.exists(overviewplot_folder))
  dir.create(overviewplot_folder)

files <- list.files(inputFolder)

pipeline <- function(path){
  Tracking <- ReadDLCDataFromCSV(file = path, fps = 30)
  Tracking <- CalibrateTrackingData(Tracking, method = "distance", in.metric = 60, points = c("tl", "br"))
  zoneinfo <- read.table("S:/Lab_Member/Tobi/Experiments/DLCAnalyzer/ArenaConfig/EPM_zoneinfo.csv", sep = ";", header = TRUE)
  Tracking <- AddZones(Tracking, zoneinfo)
  Tracking <- EPMAnalysis(Tracking, movement_cutoff = 5, integration_period = 5, points = "bodycentre", nosedips = TRUE)
  return(Tracking)
}

TrackingAll <- RunPipeline(files, inputFolder, FUN = pipeline)

Report <- MultiFileReport(TrackingAll)

# Save Report as a .csv file
report_file <- paste0(outputFolder, "Report.csv")
write.csv(Report, file = report_file, row.names = FALSE)

# Save OverviewPlots as .tiff files
library(ggplot2)
for (file in files) {
  input_file <- paste0(inputFolder, file)
  output_file <- paste0(overviewplot_folder, gsub(".csv", ".tiff", file))
  
  Tracking <- pipeline(input_file)
  plot <- OverviewPlot(Tracking, "bodycentre")
  ggsave(filename = output_file, plot = plot, device = "tiff", width = 8, height = 10)
}
