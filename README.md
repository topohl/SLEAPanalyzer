# **SLEAP/DLC Tracking Data Analysis with Modified DLCAnalyzer**

This project provides a modified version of the DLCAnalyzer R package that supports the analysis of SLEAP/DLC tracking data. Additionally, this project adds functionality to DLCAnalyzer that is not provided out of the box and extends the analysis to other behavioral tests, such as the NOR, Social PReference, Social Interactions between freely moving animals.

## Overview

[DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer) is a great tool for analyzing and visualizing tracking data from the Open Field Test, Elevated Plus Maze, and the Forced Swim Test, written in R. However, it only supports analyzing data acquired using DeepLabCut. This project extends the functionality of DLCanalyzer to support analyzing tracking data acquired using [SLEAP](https://sleap.ai).

In addition to SLEAP support, this project adds custom analyses of behavioral tests to DLCAnalyzer that are not provided out of the box. Those include latencies, frequencies, time spent in different zones when the individual is oriented towards specifc markers. This allows researchers to perform more elaborate and custom analyses using the DLCAnalyzer framework.

### Features

* **Modified DLCanalyzer:**
  - Supports analyzing SLEAP/DLC tracking data.
  - Includes custom analyses for batch-processing.

  #### Novel Object Recognition Test (NOR)
  - Import experiment metadata.
  - Define shape and size of objects.
  - Automatic application of object shapes based on metadata.
  - Define interaction angle.
  - Calculate head angle between nose, body center, and objects.
  - Calculation of contact time based on head angle towards object and shape of object.
  - Time in proximity of objects.
  - Time in proximity and oriented towards objects.
  - Latency until first contact with objects.
  - Distance.
  - Speed.
  - Automatic formatting and output of time spent with novel/familiar object.

  #### Social Preference Test (SocP)
  - Import experiment metadata.
  - Add interaction zones.
  - Add proximity zones.
  - Time in contact with familiar/novel stimulus.
  - Time in proximity with familiar/novel stimulus.
  - Latency until first contact with stimulus.
  - Distance.
  - Speed.
  - Automatic formatting and output of time spent with novel/familiar individual.

  #### Social Interaction Test (SocInt) - under construction
  - Define threshold for interaction.
  - Define interaction angle.
  - Bidirectional readout of interactions.
  - Nose-nose, nose-body, nose-tail investigation.
  - Side-by-side / side-by-side reverse.
  - Latencies and frequencies of each interaction.

* **Pipeline for SLEAP Data:**
  - Extracts x and y coordinates.
  - Formats SLEAP data using Python Jupyter Notebook.

## Installation

Installation processes for SLEAP and DLCanalyzer can be found here:
* [SLEAP](https://sleap.ai/installation.html)
* [DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer#getting-started)

To use the modified DLCAnalyzer package and the required packages, follow the installation guide by DLCAnalyzer.
For our purposes, you will also need to replace the original DLCAnalyzer package with the one provided here.

## Usage

1. **Prediction of Animal and Arena Coordinates**
   
   - Start by predicting animal and geometric (arena) coordinates.
   - Utilize the respective Python notebook provided in the `SLEAPcoords` directory to extract the coordinates from the corresponding .h5 files.

2. **Merge and Format Coordinates Files**
   
   - After extracting coordinates, proceed to merge and format the created coords files (animal and geom).
   - Use the `DLCA_Dataform` notebook for this purpose.

3. **Analysis of Merged Data**
   
   - Analyze the resulting merged file using the provided code within the `SLEAPanalyzer` directory.
   - This code relies on basic functionality of DLCAnalyzer.
   - Update the respective analysis code with the path of your merged dataframes and adjust the output directories to your preferred locations.

4. **Customization**
   
   - Modify parameters to match experimental settings or provide user-defined definitions of angles or distances as needed.
