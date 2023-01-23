# The confluence of civil conflict and cholera in Yemen
Chad R. Wells*<sup>1</sup>, Abhishek Pandey*<sup>1</sup>, Meagan Fitzpatrick<sup>1,2</sup>, William Crystal <sup>1</sup>, Madhav Marathe <sup>3</sup>, Burton H. Singer<sup>4</sup>, Alison P. Galvani <sup>1</sup>

1. Center for Infectious Disease Modeling and Analysis (CIDMA), Yale School of Public Health, New Haven, CT 06520, USA
2. Center for Vaccine Development and Global Health, University of Maryland School of Medicine, Baltimore, MD 21201, USA
3. Network Systems Science and Advanced Computing Division, Biocomplexity Institute, University of Virginia, Virginia, USA
4. Emerging Pathogens Institute, University of Florida, P.O. Box 100009, Gainesville, FL 32610, USA

Copyright (C) <2023>, Chad R. Wells et. al. All rights reserved. Released under the GNU General Public License (GPL)

This repository contains codes and data used to analyze cholera incidence trends in the governorates and districts of Yemen

The model code is written in MATLAB and results are saved as MATLAB data files (extension .mat), with plots also being constructed in MATLAB. 

## OS System requirements
The codes developed here are tested on Windows operating system (Windows 10 Home: 64-bit). However as Matlab is available for most operating systems, codes should run on Mac OSX and Linux as well.

## Installation guide
### MATLAB
Installation instruction for MATLAB can be found at https://www.mathworks.com/help/install/install-products.html. Typical install time for MATLAB on a "normal" desktop is around 30-40 minutes. The current codes were developed and tested on MATLAB R2019b.

## Demo
Figure3 produces Figure3 in the main text showing the contribution of external factors for Houthi and Government controlled regions over the four epidemic waves

## Instructions for use
To generate the Figures and output of the calculations, select a script from Figures section to run in MATLAB and enter the name in the command line. All mat file are availble to generate figures and conduct the calculations. To run analysis on a different set of parameters, adjust the parameters in the script and enter the name of the script in the command line to run.

## Computational Code
### Main Figure Script
Figure1_Control: Temporal trends in stratified regions

Figure2_Control: Spatial heterogeneity among governorates

Figure3: the contribution of external factors for Houthi and Government controlled regions over the four epidemic waves

Figure4: the contribution of external factors for Houthi and Government controlled regions and disctricts of two selected governroates

### Supplementary Scripts
FigureS_Gov_VAL: Validation of incidence for governotes

FigureS_District_VAL: Validation of incidence for districts

FigureS_Contribution: Contribution of factors at the governorate level

Figure1_Control_COVID19: The temporal trend and spatial differences of model input data prior and during COVID pandemic

Table_Mediation_Analysis_Alt: Produces table for the mediation analysis 

### Analysis Scripts
WriteTableFitting: Summarazing model fit and model selection

WriteTableVal: Summarizes model validation 

TimeLineStat and NarrativeStat: Summarizes statistics in the narrative

ParameterTables: Produces the model estimated parameters

DataCorrFigure: Computes correlations

COVID19_Covariate_Effect: Provides values for supplementary table on the effect of covariates during pandemic period

Cholera_Incidence_Governorate_COVID19: Computes the incidence during pandemic

CalibrateModel: Calibrates the saturation coefficient of the model

Fit1, Fit2, Fit3, Fit4: Used to fit all the different models

Houthi_Map_Legend: Used to produce the houthi control legend for the maps

mediation_analysis_alt: Used to do the mediation analysis for diesel and conflict

mediation_analysis_alt_wheat: Used to do the mediation analysis for wheat and conflict

PrintSummaryTable: Produces LaTEX summary table of the governorate characteristics
