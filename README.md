# Estimating Biomass Using Imaging Spectroscopy
This repository provides a reproducible workflow for modeling vegetation traits using imaging spectroscopy data. The code can be executed using either:

The included randomly generated sample dataset, or User-provided imaging spectroscopy data formatted according to the specified structure.

## Features
Preprocessing of input spectral and trait data

Model training using Partial Least Squares Regression (PLSR)

Evaluation metrics and prediction visualization

Support for topographic corrections and phenological variability

Fully annotated and adaptable R code

## Getting Started
To run the workflow:

Clone this repository.

Open the R project or R script in RStudio.

Run the script using the provided example data or load your own dataset following the format of REF_DF.

## Requirements
R version = 4.2.2

Required packages: "pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra","spectratrait","reshape","agricolae",
                      "baseline","spectrolab","caret","tidyverse","readr","tidyr","ggrepel","directlabels","ggpubr",
                      "animation","pdftools","gghighlight","stringr","gtools","gsheet","GGally","rgdal","magick",
                      "pavo","mdatools","vip","pdp","magrittr","enpls","crayon", "highlight","seecolor"

## Data Format
Your dataset should include:

Metadata columns (e.g., plot ID, year, location, elevation)

Trait measurements (e.g., biomass)

Spectral bands named as X410, X420, ... up to X2400

## License
This project is licensed under the MIT License.
