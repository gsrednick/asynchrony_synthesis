# Code supporting manuscript by Srednick & Swearer in Nature Ecology & Evolution

**Title: Understanding diversity-synchrony-stability relationships in multitrophic communities.**

**Nature Ecology & Evolution 2024 - DOI: XXXX**

**Authors: Griffin Srednick & Stephen Swearer**

**Contents**

This repository contains two main scripts that analyze the separate components of this work: 
- Part_A_systematic_review.R
- Part_B_marine_synthesis.R

***Part_A_systematic_review.R*** analyzes data and plots results from the systematic review component. The anonymized data for this script are provided in "Data/curated_literature_data_anon.csv" which is directly linked to the script. 

***Part_B_marine_synthesis.R*** calculates synchrony and stability metrics across the 5 subtidal marine time series. This script intakes summarized data that are curated in "Data/Curated_data". There is a separate curation script for each dataset that filters sites for the minimum time period and classifies herbivores when necessary (as outlined in the supplement). This script also produces all plots associated with the marine synthesis component of the manuscript.

The raw data for marine synthesis component can be sourced from links at the top of each curation script within "Data/Curated_data". However this does not apply to data from the AIMS LTMP GBR dataset, which is only accessible through reasonable request to AIMS. As such, we are not permitted to share these data in raw or summarized form.


