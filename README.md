# Stock Assessment of Arctic Grayling in Beaver and Nome Creeks, 2023

This project plan outlines a two-event mark recapture experiment on
Arctic grayling _Thymallus arcticus_ in Beaver and Nome Creeks. Sampling
will be conducted in 36 river kilometers (rkm) of the mostly
road-accessible portion of Nome Creek and in the 44 rkm floatable
portion of Nome Creek to the confluence of Beaver and Wickersham Creeks
during July 2023. This follows a 2-year radiotelemetry study on Arctic
grayling in the Beaver Creek drainage, where information acquired on
life history such as migration timing, seasonal habitat preferences, and
distribution has informed the timing, index area, and future data
analyses for this mark-recapture experiment. The marking event will
occur in early July, where three 2-person crews will sample a total of 6
rkm each day (i.e., 2 rkm per crew) and deploy approximately 1,000 Floy
tags in Arctic grayling ≥250 mm FL. The recapture event will occur in
late July, where Arctic grayling will be captured throughout the index
area using the same methodologies employed during the first event and
fish will be examined for tags. The population estimates acquired in
2023 will be compared to the last assessment conducted in 2000.

## Operational Plan

The Operational Plan associated with this study can be found at

<https://www.adfg.alaska.gov/FedAidPDFs/ROP.SF.3F.2023.07.pdf>

## Repository Contents

### Data

This folder contains all raw data associated with fish marking and recapture, reformatted as .csv files: mark_float.csv, mark_hike.csv, recap_float.csv, and recap_hike.csv.  These data are read into R and additional data manipulation and error correction is performed in script.

Two additional files are included: site_boundaries.csv gives the starting and ending coordinates of all sampling reaches, and beaver_cr_rivernetwork_op.Rdata is an R workspace file that allows import of the rivernetwork associated with the sampling area, formatted for use with the riverdist package.

### R

**1_BCMR_data.R** loads all data and performs all data manipulation and reformatting and error correction.

**2_BCMR_assumptions.R** performs all statistical tests and checks to validate the mark-recapture assumptions.

**3_BCMR_estimates.R** produces output relevant to estimation of abundance and length composition.

### R_output

This folder contains a "cleaned" data file, **BCMR_alldata.csv**.

### BCMR_output.Rmd 

This Markdown file writes a Word document summarizing the results of assumption validation and estimation of abundance and length composition.
