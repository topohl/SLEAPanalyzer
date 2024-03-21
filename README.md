# **SLEAP/DLC Tracking Data Analysis with Modified DLCAnalyzer**

This project provides a modified version of the DLCAnalyzer R package that supports the analysis of SLEAP/DLC tracking data. Additionally, this project adds functionality to DLCAnalyzer that is not provided out of the box.

## Overview

[DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer) is a great tool for analyzing and visualizing tracking data from various animal behavioral tests, written in R. However, it only supports analyzing data acquired using DeepLabCut. This project extends the functionality of DLCanalyzer to support analyzing tracking data acquired using [SLEAP](https://sleap.ai).

In addition to SLEAP support, this project adds custom analyses to DLCAnalyzer that are not provided out of the box. This allows researchers to perform custom analyses using the DLCAnalyzer framework.

### Features

* Modified version of DLCanalyzer that supports analyzing SLEAP/DLC tracking data.
* Custom analyses added to DLCanalyzer that are not provided out of the box. Made for batch-processing
  - Novel object recognition test (NOR)
    - Import experiment metadata
    - Define shape and size of objects
    - Automatic application of object shapes based on metadata
    - Define interaction angle
    - Calculate head angle between nose, body center, and objects
    - Calculation of contact time based on head angle towards object and shape of object
    - Time in proximity of objects
    - Time in proximity and oriented towards objects
    - Latency until first contact with objects
    - Distance
    - Speed
    - Automatic formatting and output of time spent with novel / familiar object
      
  - Social preference test (SocP)
    - Import experiment metadata
    - add interaction zones
    - add proximity zones
    - Time in contact with familiar / novel stimulus
    - Time in proximity with familiar / novel stimulus
    - Latency until first contact with stimulus
    - Distance
    - Speed
    - Automatic formatting and output of time spent with novel / familiar individual

  - Social interaction test (SocInt) - under construction
    - Define threshold for interaction
    - Define interaction angle
    - Bidirectional readout of interactions
    - Nose-nose, Nose-body, Nose-tail investigation
    - Side-by-side / Side-by-side reverse
    - Latencies and Frequencies of each interaction
  
* Pipeline to extract x and y coordinates and format SLEAP data using Python Jupyter Notebook.

## Installation

Installation processes for SLEAP and DLCanalyzer can be found here:
* [SLEAP](https://sleap.ai/installation.html)
* [DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer#getting-started)

To use the modified DLCAnalyzer package and the required packages, follow the installation guide by DLCAnalyzer.
For our purposes, you will also need to replace the original DLCAnalyzer package with the one provided here.

## Usage

After prediction of animal and geometric (arena) coordinates, extract the coordinates from respective .h5 files using the respective python notebook provided in SLEAPcoords.
Afterwards, merge and format the created coords files (animal and geom) using the DLCA_Dataform notebook.
The resulting file can be analyzed using the provided code within the SLEAPanalyzer directory which relies on basic functionality of DLCAnalyzer.
For this, update the respective analysis code with the path of your merged dataframes and adjust the output directories to where you would like to have the analysis file.
Following, change parameters to experimental settings or user-provided definitions of angles or distances.
