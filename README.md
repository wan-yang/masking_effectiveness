# Model code used in Yang & Shaman "Reconciling the efficacy and effectiveness of masking on epidemic outcomes"
author: Wan Yang\
date: March 2023\
note: this code is for research use only; may not be entirely optimized nor neatly organized. 

## Study Summary
During the COVID-19 pandemic, mask wearing in public settings has been a key control measure. However, the low effectiveness reported for masking has cast doubt on its validity. This study develops an agent-based model to interrogate influencing factors. Testing shows that transmission within-household where masks are rarely used can substantially lower effectiveness. Nonetheless, model results support the effectiveness of masking at both the individual and population levels, albeit at less-than-ideal levels. Overall, study findings indicate it is prudent for individuals to use masks during an epidemic, and for policy makers to recognize the less-than-ideal effectiveness of masking when devising interventions.

## data
This folder includes all relevant data needed to run the simulations.\
The POLYMOD contact data were downloaded from https://zenodo.org/record/3874557#.ZAiucy3MyF0\
The Google mobility data were downloaded from https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv\
The NYC vaccination data were downloaded from https://raw.githubusercontent.com/nychealth/covid-vaccine-data/main/people/trends-byage.csv

## scripts
This folder includes all model code used in the study. 
1. To compile and process relevant data, use the script "compile_*.R"
2. To run the model simulations, use the script "driver_ABMmasking.R":
- core models and functions are in the scripts "Fn_ABMmasking.R" and "Fns_util_AMBmasking.R"
- model settings for each analysis (main analysis, counterfactual assuming people mask at home, sensitivity on weights of nonphysical contact, and higher R0) are in the scripts "get_*settings.R"
3. To combine model runs, use the scripts "Fn_read_res*.R"
4. To compute the summary statistics and effectiveness estimates etc., use the scripts "Fn_cal_*.R"
5. To analyze results and plot the figures, use the script "ana_mask_runs_figs.R"
6. To generate the tables, use the script "output_tables.R"

### results
This folder includes the simulation results, as labeled. 

### figures
This folder includes all figures for the study.

### tables
This folder includes all tables for the study. 

### contact
Wan Yang: wy2202 at cumc.columbia.edu

### reference
Version 1: 
Yang W & Shaman J. 2023 Reconciling the efficacy and effectiveness of masking on epidemic outcomes. medRxiv 2023.05.10.23289803; doi: https://doi.org/10.1101/2023.05.10.23289803 

