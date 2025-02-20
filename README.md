# Eastern Indigo Snake Project
This project resulted in two manuscripts; an experimental study looking at the persistence and accumulation of environmental DNA (eDNA) and a field study comparing eDNA to camera traps and visual searches as monitoring methods for the eastern indigo snake. The datasets and scripts are described below. 

# Persistence of reptile DNA in a terrestrial substrate: a case study using the eastern indigo snake
Datasets and scripts used for Samuels, L. R. N., Chandler, H. C., Hoffman, M., Kronenberger, J. A., Elmore, M., Aldredge, R., Stegenga, B. S., Bogan, J. E., Davis, M. A., Hertz, S., Schwartz, M. K., & Wilcox, T. (2025). Persistence of Reptile DNA in a Terrestrial Substrate: A Case Study Using the Eastern Indigo Snake. Environmental DNA, 7(1), e70053. https://doi.org/10.1002/edn3.70053

# File guide

**Enclosure_qPCR.csv** - Dataset containing qPCR data for all 234 soil samples, required for Enclosure_study.R script.

**Enclosure_schedule.csv** - Dataset containing sampling schedule, quadrat and snake IDs, and time of snake presence in each quadrat, required for Enclosure_study.R script.

**Enclosure_study.R** - Script used for paper figures, data analysis, detection probability modeling, and bootstrapped 95% confidence intervals. 

**Enclosure_weather.csv** - Dataset containing all weather data compiled from HOBO lux and temperature loggers, and from the Leesburg International Airport NOAA station, required for Enclosure_study.R script.

# Comparison of three monitoring methods in situ for threatened species detection
Datasets and scripts used for ...

# File guide

**qPCR_mastersheet_field.csv** - Dataset containing qPCR data for all 306 eDNA soil samples, required for Field_study.R script.

**Field_study.R** - Script used to restructure data. Required for Positives.R.

**Field_weather.csv** - Weather data from FAWN tower and taken by field technicians. Required for Positives.R.

**camera_positives.csv** - eDNA samples confirmed positive by cameras. Required for Positives.R.

**Positive Sample Details_TOS.csv** - Details about location samples were taken at (burrows, fences, or ground).

**Positives.R** - Analysis of all known-positives, includes code for Figure 2 and Figure 3. Enclosure_study.R and Field_study.R are required for this script.

**EIS Survey Cost Estimates.xlsx** - Dataset with overview of all cost estimates, used for analysis in occupancy.R.

**eDNA_data**- Dataset of all eDNA detections from all surveys for all sites.

**VES_data** - Dataset of all visual encounter survey detections from all surveys for all sites.

**Camera_data** - Dataset of all camera detections from all surveys for all sites.

**presence_input** - Dataset used for R Presence occupancy modeling. Required for occupancy.R.

**rate_of_detection** - Dataset used for rate of detection calculations. Required for occupancy.R

**occupancy.R** - Script used for occupancy modeling and analysis. Includes code for Figure 4-7. 


# Contact information
Please reach out to us at the [National Genomic Center for Wildlife and Fish Conservation](https://www.fs.usda.gov/research/rmrs/projects/ngc) with any questions or comments. Data was analyzed by Leah Samuels at leah.samuels@usda.gov/leahrnsamuels@gmail.com and Taylor Wilcox at taylor.wilcox@usda.gov.
