# **SLEAP/DLC Tracking Data Analysis with Enhanced DLCAnalyzer**

This repository hosts an extended version of the **DLCAnalyzer** R package, providing support for analyzing tracking data from both SLEAP and DLC and introducing new analysis features for various behavioral tests. These tests include the **Novel Object Recognition (NOR) Test**, **Social Preference (SocP) Test**, and **Social Interaction (SocInt) Test**.

## Overview

[DLCAnalyzer](https://github.com/ETHZ-INS/DLCAnalyzer) is a powerful R package designed to analyze and visualize tracking data from behavioral tests like the Open Field Test, Elevated Plus Maze, and Forced Swim Test, specifically using data collected via DeepLabCut (DLC). This project builds on DLCAnalyzer to support tracking data from both [SLEAP](https://sleap.ai) and DLC, making it versatile for researchers who work with either platform.

In addition to SLEAP compatibility, this extended DLCAnalyzer version enables in-depth behavioral analyses. New metrics include latencies, frequencies, and time spent in specific zones relative to designated markers, expanding the toolkit for behavioral data analysis.

## Key Features

### Expanded DLCAnalyzer Functionality
- **Support for SLEAP/DLC Tracking Data**: Analyze and visualize data tracked by either SLEAP or DLC.
- **Batch Processing and Custom Analyses**: Simplify batch processing and enable custom analysis options for your experiments.

#### Novel Object Recognition Test (NOR)
- Import experimental metadata (e.g., object location, test duration).
- Define and assign object shapes (e.g., square, round) based on metadata.
- Analyze head angles relative to objects and calculate time spent in proximity.
- Compute object-oriented metrics, including:
  - **Latency** until first interaction with objects.
  - **Contact time** based on head angle relative to object.
  - **Proximity-based metrics**, including time in proximity and orientation toward objects.
  - **Distance and Speed** measures.
- Automatically format and export results for interactions with novel and familiar objects.

#### Social Preference Test (SocP)
- Import metadata for experiment details.
- Define **interaction** and **proximity zones** around familiar/novel stimuli.
- Calculate metrics such as:
  - **Time in contact** and **time in proximity** with each stimulus.
  - **Latency** to first contact with familiar or novel stimuli.
  - **Distance and Speed** measures.
- Automatic export of formatted results for interactions with familiar and novel individuals.

#### Social Interaction Test (SocInt) – *In Development*
- Define thresholds for interaction metrics.
- Analyze bidirectional interactions (e.g., **nose-to-nose**, **nose-to-body**, **nose-to-tail**).
- Calculate metrics for **side-by-side** and **side-by-side reverse** positioning.
- Track latencies and frequencies for each interaction type.

---

## Installation

1. **Install SLEAP and DLCAnalyzer**:
   - Follow the [SLEAP Installation Guide](https://sleap.ai/installation.html).
   - Follow the [DLCAnalyzer Installation Guide](https://github.com/ETHZ-INS/DLCAnalyzer#getting-started).

2. **Download and Set Up This Extended DLCAnalyzer**:
   - Clone or download this repository, which includes modified files needed to support SLEAP data and custom analyses.
   - Install any additional dependencies as outlined in the DLCAnalyzer guide.

---

## Usage Guide

Follow these steps to use the extended DLCAnalyzer for analyzing your tracking data:

### Step 1: Predict Animal and Arena Coordinates
   - Use the provided Python notebook in the `SLEAPcoords` directory to predict coordinates for animals and arena objects.
   - For SLEAP users:
      - Extract coordinates from `.h5` files for each experiment using the experiment-specific Jupyter notebook in `SLEAPcoords`.
   - For DLC users:
      - DLC provides a formatted .csv file, just merge the animal and geom files.

### Step 2: Merge and Format Coordinate Files
   - With coordinates extracted, use the `DLCA_Dataform` notebook to merge and format the coordinates for animals and arena components.
   - This notebook allows you to:
      - Combine data files from the animal and geom file.
      - Apply custom formatting options to standardize data across experiments.

### Step 3: Analyze the Merged Data
   - Load the merged data file into the `SLEAPanalyzer` code directory.
   - Set up the analysis script by specifying:
      - The path to your data files.
      - Output directories where results will be saved.
   - Run analyses according to your chosen behavioral test (e.g., NOR, SocP), using the modified DLCAnalyzer functions provided here.

### Step 4: Customize Analysis Parameters
   - Customize analysis parameters to fit your experimental settings. For example:
      - Define custom angles for orientation toward objects in NOR tests.
      - Set proximity distances for social preference interactions in SocP tests.
      - Adjust interaction thresholds for bidirectional metrics in SocInt tests.
   - Modify these parameters directly in the R scripts or configuration files to ensure they align with your study design.

### Step 5: Export and Interpret Results
   - Analysis results will be automatically saved in the output directory specified.
   - Results include formatted metrics like time in proximity, latency to interact, and movement speed.
   - Review and interpret these results based on the experimental goals for your behavioral study.

---

## Contributing

Contributions are welcome! If you’d like to help improve or extend this project, please feel free to submit a pull request or open an issue with ideas, bug reports, or suggestions.

## Support

For questions or feedback, please reach out via [issues](https://github.com/topohl/SLEAPanalyzer/issues) on GitHub or email us directly.

---

Feel free to reach out with questions or suggestions. This project aims to facilitate enhanced behavioral analysis and support the scientific community in generating high-quality insights from SLEAP/DLC tracking data.
