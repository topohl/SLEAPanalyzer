# **SLEAP/DLC Tracking Data Analysis with Modified DLCAnalyzer**

This project provides a modified version of the DLCAnalyzer R package that supports the analysis of SLEAP/DLC tracking data. Additionally, this project adds functionality to DLCAnalyzer that is not provided out of the box.

## Overview

[DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer) is a great tool for analyzing and visualizing tracking data from various animal behavioral tests, written in R. However, it only supports analyzing data acquired using DeepLabCut. This project extends the functionality of DLCanalyzer to support analyzing tracking data acquired using [SLEAP](https://sleap.ai).

In addition to SLEAP support, this project adds custom analyses to DLCAnalyzer that are not provided out of the box. This allows researchers to perform custom analyses using the DLCAnalyzer framework.

### Features

* Modified version of DLCanalyzer that supports analyzing SLEAP/DLC tracking data.
* Custom analyses added to DLCanalyzer that are not provided out of the box. Made for batch-processing
  - Novel object recognition test (NOR)
    - import experiment metadata
    - calculate head angle between nose, bodycentre, and objects
    - define interaction angle
    - define shape and size of objects
    - automatic application of object shapes based on metadata
    - calculation of contact time based on head angle towards object and shape of object
    - time in proximity of objects
    - time in proximity and oriented towards objects
    - latency until first contact with objects
    - distance
    - speed
    - ...
      
  - Social preference test (SocP)
    - add interaction zones
    - add proximity zones
    - calculate time of interaction

  - Social interaction test (SocInt)
  
* Pipeline to extract x and y coordinates and format SLEAP data using Python Jupyter Notebook.

## Installation

Installation processes for SLEAP and DLCanalyzer can be found here:
* [SLEAP](https://sleap.ai/installation.html)
* [DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer#getting-started)

To install the modified DLCanalyzer package and the required packages, run the following command in R:

## Usage

After prediction of animal and geometric (arena) coordinates, extract the coordinates from respective .h5 files using the respective python notebook provided in SLEAPcoords.
Afterwards, merge and format the created coords files (animal and geom) using the DLCA_Dataform notebook.
The resulting file can be analyzed using the provided code within the SLEAPanalyzer directory which relies on basic functionality of DLCAnalyzer.
