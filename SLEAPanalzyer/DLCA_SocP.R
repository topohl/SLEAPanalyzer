# -------------------------------------------------------------------
## install tensor flow and prepare for libraries and add new environment r-reticulate

install.packages("tensorflow")

library(reticulate)
path_to_python <- install_python()
virtualenv_create("r-reticulate", python = path_to_python)

library(tensorflow)
install_tensorflow(envname = "r-reticulate")

install.packages("keras")
library(keras)
install_keras(envname = "r-reticulate")

library(tensorflow)
tf$constant("Hello Tensorflow!")


# -------------------------------------------------------------------
## load dependencies

library(sp)         #tested with v1.3-2
library(imputeTS)   #tested with v3.0
library(ggplot2)    #tested with v3.1.0
library(ggmap)      #tested with v3.0.0
library(data.table) #tested with v1.12.8
library(cowplot)    #tested with v0.9.4
library(corrplot)   #tested with v0.84
library(keras)      #REQUIRES TENSORFLOW INSTALL. tested with v2.2.5.0


# -------------------------------------------------------------------
## for single file analysis of SocP time in zones data

# Set working directory and load R script
setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
source('R/DLCAnalyzer_Functions_final.R')

# Read in tracking data and get names of variables in the data frame
input_file <- "C:/Users/topohl/Documents/GitHub/sleap/docs/notebooks/SocP_F9L3_GX0788_S2_locs.csv"
Tracking <- ReadDLCDataFromCSV(file = input_file, fps = 30)

# Extract the input file name without the extension
input_file_name <- sub(".csv$", "", basename(input_file))

Tracking <- CalibrateTrackingData(Tracking, method = "area",in.metric = 44*24, points = c("tl","tr","br","bl"))

Tracking <- AddOFTZones(Tracking, scale_center = 0.5,scale_periphery  = 0.8 ,scale_corners = 0.4, points = c("tl","tr","br","bl"))
PlotZones(Tracking)

PlotPointData(Tracking, points = c("nose"))

# Calculate time_Left and time_Right
total_time <- 10 * 60
IsInZoneLeft <- GetDistances(Tracking, "socl", "nose") < 6
IsInZoneRight <- GetDistances(Tracking, "socr", "nose") < 6
IsInZoneLeft[is.na(IsInZoneLeft)] <- TRUE
IsInZoneRight[is.na(IsInZoneRight)] <- TRUE
time_Left <- sum(IsInZoneLeft) * total_time / length(IsInZoneLeft)
time_Right <- sum(IsInZoneRight) * total_time / length(IsInZoneRight)

# Create a data frame with the time_Left and time_Right values
df <- data.frame(time_Left = time_Left, time_Right = time_Right)

# Construct the output file name
output_file <- paste0(input_file_name, "_output.csv")

# Write the data frame to a csv file with the constructed file name
write.csv(df, output_file, row.names = FALSE)




