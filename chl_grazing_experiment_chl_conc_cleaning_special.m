%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for special Chl-a conc cleaning of the raw Chl-a
% obtained during dilution (grazing) experiments. 
%
% Special Chl-a conc cleaning:
% The goal of this process is to correct some T0 Chl-a concentration for
% few stations. In some cases (described below), the values of the T0 dil 
% were missing, making the rates calculation impossible. In these case, the
% dilution values obatined from the other filters or from Flow-Cytometry
% (FCM) were used.
% EN661 - L2: T0 dil >0&<200 (GFF) Chl-a conc values were way too high
% (>60%), while the dilution obatined from FCM and 10um filters gave the
% same dilution values = 23%. T0 dil >0&<200 Chl-a conc are then set as 23%
% of T0 wsw >0&<200 Chl-a conc.
% EN668 - L6-D2: All T0 dil values (>0&<200, >10&<200 and >0&<10) were
% quationable with dilution ranging from 43% to 76%, while according to the
% FCM we had 20% dilution. T0 wsw >0&<200 values looks correct when
% compared to post-calibrated underway fluorescence. All the T0 wsw values
% are then considered are good and a T0 dil = 20% T0 wsw.
%
% Input: CRUSIE-chl-grazing-experiments-chl-conc-clean.csv files 
%
% Outputs: CRUISE-chla-grazing-experiments-clean.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% June 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\Desktop\NES-LTER_Chla_Cleaning_Rates_Computation\';
%Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-chl-conc-clean\');
%Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-clean\');

%Find all the *.cnv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);%List all files of interest in the directory

% EN661 - L2 - cast 2:  T0 dil >0&<200 (GFF) Chl-a = 23% T0 wsw >0&<200 Chl-a conc
n1=10;
%load the .csv file of the corresponding cruise
tablename=strcat(rep1,list(n1).name);
T1=readtable(tablename);
b1=strcmp(T1.cast,'2''');
%Identify all values obtained with >0&<200 filters at T0 dil
c1=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'dil');
%Get the mean value of T0 wsw >0&<200
c2=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
chl_avg=mean(T1.chl(c2));
%Replace the T0 dil >0&<200 values by 23% of wsw
T1.chl(c1)=0.23*chl_avg;
%Assigned a iode_quality_flag of 1 for all these values
T1.iode_quality_flag(c1)=1;
%Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
%cruise
newname=strrep(list(n1).name,'chl-conc-clean','clean');%Replace raw by clean
newtablename=strcat(rep2,newname);%New tablename and path
writetable(T1,newtablename)

clear n1 b1 c1 c2 chl_avg

% EN668 - L6-D2 - cast 18, niskins 21-22-23-24:  
% T0 dil all filters Chl-a = 20% T0 wsw all filters Chl-a conc
n1=11;
%load the .csv file of the corresponding cruise
tablename=strcat(rep1,list(n1).name);
T1=readtable(tablename);
b1=strcmp(T1.cast,'18''') & strcmp(T1.niskin,'21-22-23-24''');
%Identify all values obtained with >0&<200 filters at T0 dil
c1=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'dil');
%Get the mean value of T0 wsw >0&<200
c2=b1 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
chl_avg=mean(T1.chl(c2));
%Replace the T0 dil >0&<200 values by 20% of wsw
T1.chl(c1)=0.20*chl_avg;
%Assigned a iode_quality_flag of 1 for all these values
T1.iode_quality_flag(c1)=1;
clear c1 c2 chl_avg
%Identify all values obtained with >10&<200 filters at T0 dil
c1=b1 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'dil');
%Get the mean value of T0 wsw >0&<200
c2=b1 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
chl_avg=mean(T1.chl(c2));
%Replace the T0 dil >0&<200 values by 20% of wsw
T1.chl(c1)=0.20*chl_avg;
%Assigned a iode_quality_flag of 1 for all these values
T1.iode_quality_flag(c1)=1;
clear c1 c2 chl_avg
%Identify all values obtained with >0&<10 filters at T0 dil
c1=b1 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'dil');
%Get the mean value of T0 wsw >0&<200
c2=b1 & strcmp(T1.filter_size,'>0&<10') & strcmp(T1.T0_TF,'T0') & strcmp(T1.dilution,'wsw') & T1.iode_quality_flag==1;
chl_avg=mean(T1.chl(c2));
%Replace the T0 dil >0&<200 values by 20% of wsw
T1.chl(c1)=0.20*chl_avg;
%Assigned a iode_quality_flag of 1 for all these values
T1.iode_quality_flag(c1)=1;
%Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
%cruise
newname=strrep(list(n1).name,'chl-conc-clean','clean');%Replace raw by clean
newtablename=strcat(rep2,newname);%New tablename and path
writetable(T1,newtablename)

clear n1 c1 c2 chl_avg

%Make some copies of all the other csv files in the new directory rep2,
%with a new name, without any change
for n1=1:8
    %load the .csv file of the corresponding cruise
    tablename=strcat(rep1,list(n1).name);
    T1=readtable(tablename);
    %Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    %cruise
    newname=strrep(list(n1).name,'chl-conc-clean','clean');%Replace raw by clean
    newtablename=strcat(rep2,newname);%New tablename and path
    writetable(T1,newtablename)
end
