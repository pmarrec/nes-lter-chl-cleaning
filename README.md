# nes-lter-chl-cleaning
You can find here all the raw fluorescence data from all nes-lter cruises and the matlab scripts to perform an automated cleaning of the chl-a data that are used for growth/grazing calculation. Each cruise has separated csv files and data for all cruises are processed simultaneously.

The raw fluorescence data for each cruise are stored in csv files in the "chl-grazing-experiment-raw" folder. The csv files contain the cast and niskin bottle numbers (in text format, with '' to avoid considering some multiple cast/niskins entries as dates), the measured Fo, Fa, blank, blank_acid data with treatment information (T0_TF, dilution, nutrient_treatment, light_level, replicate_bottle and replicate_chl). The files also include the vol_filtered, the vol_extracted, the calibration coeficient Fs and r and the start/end datetime of the incubations. Empty (NaN) columns for: date_time_utc_sampling, latitude, longitude, depth, temperature_sampling, nearest_sation, distance, incubation_tank, temperature_incubation_avg, temperature_incubation_std, are there for retrieving these additional data from the CTD files stored in the REST-API and from the csv files containing the temperature recorded in the incubation tanks during each cruise.

The **chl_grazing_experiment_retrieve_ctd_data** Matlab script is used to retrieve the date_time_utc_sampling, latitude, longitude, depth, temperature_sampling, nearest_sation and distance values ***(still need to find how to get the nearest_station and distance values)*** directly and automatically from the CTD files stored in the REST-API (except for AR66, data not available yet).\
*Inputs: CRUSIE-chl-grazing-experiments-raw.csv files*\
*Outputs: CRUISE-chla-grazing-experiments-ctd-raw.csv files stored in the chl-grazing-experiment-ctd-raw folder*

The **chl_grazing_experiment_retrieve_temp_inc_data.m** Matlab script is used to retrieve the mean/std temperature in the tanks during dilution experiments. For each cast/depth, get the mean/std temperature recorded by the HOBO data loggers in the different tanks. Because of the lack of records in terms of which bottles in which tanks (there are up to 3 CRUISE_temp_inc_X.csv files for each cruise), the mean/std in each tank are retrieved and the values with the lowest temperature difference will be kept and saved automatically.
_Input: CRUSIE-chl-grazing-experiments-ctd-raw.csv files and CRUISE_temp_inc_X.csv files_
_Outputs: CRUISE-chla-grazing-experiments-temp-inc-raw.csv files._

The **chl_grazing_experiment_chl_calc.m** Matlab script is used to calculate the Chl-a and Phaoe pigments concentration, after calculation of Fo-blank, Fa-blank and Fo/Fa ratios.
_Inputs: CRUSIE-chl-grazing-experiments-raw.csv files with Fo, Fa, blank values and calibration coefficients Fs and r_
_Outputs: CRUISE-chla-grazing-experiments-chl-calc.csv files stored in the chl-grazing-experiment-chl-calc folder_

The **chl_grazing_experiment_fofa_cleaning.m** Matlab script is used to perform the 1st step of data quality check (aka cleaning). This first step is based on 2 criteria:
1) Fo/Fa within 1-3 range,
2) Fo/Fa witihn +/- 2 StdDev confidence interval for a given type of filter
All the values that don't fit these criteria are flagged as questionable with a iode_quality_flag = 3
All the other values are flagged as good with a iode_quality_flag = 1
_Input: CRUSIE-chl-grazing-experiments-chl-calc.csv files with Chl-a_
_Outputs: CRUISE-chla-grazing-experiments-fo-fa-clean.csv files stored in the chl-grazing-experiment-fofa-clean folder_

The **chl_grazing_experiment_chl_conc_clean.m** Matlab script is used to perform the 2nd step of the cleaning based on the Chl-a concentrations. This second step is based on 2 criteria:
1) All negative Chl-a conc are flagged as questionnable/suspect with a iode_quality_flag = 3
2) For each station/depth/treatment/triplicate values:
Each triplicate value should stand in the +/- 2 x %CV of the mean values of the triplicate with a QC flag=1. %CV is considered as the mean %CV obtained on a given type of filter (GFF/10um) at T0 and at TF.
All the values that don't fit these criteria are flagged as questionable with a iode_quality_flag = 3
All the other values are flagged as good with a iode_quality_flag = 1
_Input: CRUSIE-chl-grazing-experiments-fofa-clean.csv files_
_Outputs: CRUISE-chla-grazing-experiments-chl-conc-clean.csv files._

The **chl_grazing_experiment_chl_conc_clean_special.m** Matlab script is to correct some T0 Chl-a concentration for few stations. In some cases (described below), the values of the T0 dil were missing, making the rates calculation impossible. In these case, the dilution values obatined from the other filters or from Flow-Cytometry (FCM) were used:
1) EN661 - L2: T0 dil >0&<200 (GFF) Chl-a conc values were way too high (>60%), while the dilution obatined from FCM and 10um filters gave the same dilution values = 23%. T0 dil >0&<200 Chl-a conc are then set as 23% of T0 wsw >0&<200 Chl-a conc.
2) EN668 - L6-D2: All T0 dil values (>0&<200, >10&<200 and >0&<10) were questionable with dilution ranging from 43% to 76%, while according to the FCM we had 20% dilution. T0 wsw >0&<200 values looks correct when compared to post-calibrated underway fluorescence. All the T0 wsw values are then considered are good and a T0 dil = 20% T0 wsw.
_Input: CRUSIE-chl-grazing-experiments-chl-conc-clean.csv files_
_Outputs: CRUISE-chla-grazing-experiments-clean.csv files._
