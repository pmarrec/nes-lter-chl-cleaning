# nes-lter-chl-cleaning
You can find here all the raw fluorescence data from all nes-lter cruises and the matlab scripts to perform an automated cleaning of the chl-a data that are used for growth/grazing calculation. Each cruise has separated csv files and data for all cruises are processed simultaneously

The raw fluorescence data for each cruise are stored in csv files in the "chl-grazing-experiment-raw" folder. The csv files contain the cast and niskin bottle numbers (in text format with '' to avoid considering some multiple cast/niskins entries as dates), the measured Fo, Fa, blank, blank_acid data with treatment information (T0_TF, dilution, nutrient_treatment, light_level, replicate_bottle and replicate_chl. The files also include the vol_filtered, the vol_extracted, the calibration coeficcient Fs and r and the start/end datetime of the incubations. date_time_utc_sampling, latitude, longitude, depth, nearest_sation and distance values will be retrieved later directly and automatically from the CTD files stored in the REST-API.

The chl_grazing_experiment_chl_calc.m Matlab script is used to calculate the Chl-a and Phaoe pigments concentration, after calculation of Fo-blank, Fa-blank and Fo/Fa ratios.
Inputs: CRUSIE-chl-grazing-experiments-raw.csv files with Fo, Fa, blank values and calibration coefficients Fs and r
Outputs: CRUISE-chla-grazing-experiments-chl-calc.csv files stored in the chl-grazing-experiment-chl-calc folder

The chl_grazing_experiment_fofa_cleaning.m Matlab script is used to perform the 1st step of data quality check (aka cleaning). This first step is based on 2 criteria:
1) Fo/Fa within 1-3 range,
2) Fo/Fa witihn +/- 2 StdDev confidence interval for a given type of filter
All the values that don't fit these criteria are flagged as questionable with a iode_quality_flag = 3
All the other values are flagged as good with a iode_quality_flag = 1
Input: CRUSIE-chl-grazing-experiments-chl-calc.csv files with Chl-a
Outputs: CRUISE-chla-grazing-experiments-clean.csv files stored in the chl-grazing-experiment-fofa-clean folder

The chla_grazing_experiment_chl_conc_clean.m Matlab script is used to perform the 2nd step of the cleaning based on the Chl-a concentrations. This second step is based on 2 criteria:
1) All negative Chl-a conc are flagged as questionnable/suspect with a iode_quality_flag = 3
2) For each station/depth/treatment/triplicate values:
Each triplicate value should stand in the +/- 2 x %CV of the mean values of the triplicate with a QC flag=1. %CV is considered as the mean %CV obtained on a given type of filter (GFF/10um) at T0 and at TF.
All the values that don't fit these criteria are flagged as questionable with a iode_quality_flag = 3
All the other values are flagged as good with a iode_quality_flag = 1
Input: CRUSIE-chl-grazing-experiments-fofa-clean.csv files 
Outputs: CRUISE-chla-grazing-experiments-chl-conc-clean.csv files.

