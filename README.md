# Estimating Biomass Using Imaging Spectroscopy

This repository contains a reproducible workflow for estimating vegetation biomass in alpine grasslands using imaging spectroscopy data across multiple years.

The approach combines spectral analysis, trait data integration, and Partial Least Squares Regression (PLSR) modeling to explore relationships between spectral reflectance and plant trait, such as biomass.


## 🛠 Getting Started

This project can be run with:

Provided randomly generated sample data for demonstration

User-supplied imaging spectroscopy data following the expected input format

### Reproducing the Analysis

To reproduce the full workflow, refer to:

The PLSR\_Models.R script in this repository or

The Methods section of the paper
Estimating Biomass in Alpine Grasslands Using Multi-Year Imaging Spectroscopy Data


## 📁 Repository Structure


###### **Core Files**

* PLSR\_Models.R – Main modeling script
* Example dataset or template for user input
* README.md – Project overview


###### **Spatial Processing \& Visualization Notebooks**

These notebooks support the spatial preprocessing and visualization steps described in the manuscript.

* 020\_Slope\&Aspect\_v-shared.ipynb

Computes slope and aspect maps from a digital terrain model, generates hillshade layers for each flight year, and prepares terrain attributes for comparison with biomass estimates.

* 021\_Study\_Area\_Map\_v-shared.ipynb

*Produces the study area map showing park boundaries, APEX flight coverage, and field plot locations. Reproduces Figure 1 in the manuscript.*

* 022\_FL\_overlay\_v-shared.ipynb

*Overlays flight-line footprints onto the study area and visualizes acquisition geometry relevant for illumination and BRDF considerations.*

* 023\_Clipped\_FL\_Maps\_v-shared.ipynb

*Clips flight-line data to the analysis area and generates spatial layers used to compare reflectance and biomass across topographic corrections.*

* 024\_Masked\_Maps\_v-shared.ipynb
* 
*Applies terrain and vegetation masks to spatial layers for consistency across years and correction outputs. Produces masked maps used in Figures and Appendix materials.*


## 📦 Dependencies


###### **Make sure the following R packages are installed before running the code:**

pls, dplyr, reshape2, here, plotrix, ggplot2, gridExtra, spectratrait, reshape, agricolae,
baseline, spectrolab, caret, tidyverse, readr, tidyr, ggrepel, directlabels, ggpubr, animation,
pdftools, gghighlight, stringr, gtools, gsheet, GGally, rgdal, magick, pavo, mdatools, vip,
pdp, magrittr, enpls, crayon, highlight, seecolor

**Install all at once:**

packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra","spectratrait",
"reshape","agricolae","baseline","spectrolab","caret","tidyverse","readr",
"tidyr","ggrepel","directlabels","ggpubr","animation","pdftools","gghighlight",
"stringr","gtools","gsheet","GGally","rgdal","magick","pavo","mdatools","vip",
"pdp","magrittr","enpls","crayon","highlight","seecolor")

install.packages(setdiff(packages, installed.packages()\[,"Package"]))


###### **🐍 Python Environment (for spatial preprocessing notebooks)**

The Python Jupyter Notebooks in this repository (files 020\_\*.ipynb to 024\_\*.ipynb) require a dedicated geospatial environment.

A reproducible Conda environment file is provided:

* geopandas\_env.yml

**Create the environment:**

conda env create -f geopandas\_env.yml

conda activate geopandas\_env


## 📊 Input Data Format

Your input dataset should include:

Metadata (e.g., Year, Plotname, Elevation)

Trait values (e.g., Biom\_wet\_g)

Reflectance values across wavelengths (e.g., X410, X420, ..., X2400)

See Randomly\_generated\_sample\_data\_structure.csv for an example format.


## 📈 Outputs

* PLSR model summaries and selected components
* VIP scores (Variable Importance in Projection)
* Predicted vs. observed trait plots
* Trait-specific model diagnostics
* Spatial maps of terrain attributes, flight lines, and biomass distributions


## ✉️ Support

For questions or suggestions, feel free to open an issue or contact the repository maintainer.

