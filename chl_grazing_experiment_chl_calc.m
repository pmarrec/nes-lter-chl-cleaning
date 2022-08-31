%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matlab script for Chl-a conc computation from the raw fluorescence Fo and
%Fa obtained during dilution (grazing) experiments.
%
%Chl-a computation:
%Chl-a = [Fs ([r/(r-1)] (Fo-Fa))] * vol_extracted / vol_filtered
%Phaeo = [Fs([r/(r-1)](r*Fa-Fo))] * vol_extracted / vol_filtered
%
%Input: CRUSIE-chl-grazing-experiments-raw.csv files with Fo, Fa,
%blank values and calibration coefficients Fs and r
%
%Outputs: CRUISE-chla-grazing-experiments-chl-calc.csv files.
%
%Written by Pierre Marrec
%
%pmarrec@uri.edu
%
%5/21/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\My Drive\NES-LTER_Chla_Cleaning_Rates_Computation\';
%Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-temp-inc-raw\');
%Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-chl-calc\');

%Find all the *.cnv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);%List all files of interest in the directory

for n1=1:numel(list)
    %load the .csv file of the corresponding cruise
    tablename=strcat(rep1,list(n1).name);
    T1=readtable(tablename);

    %Calculation of fo_blank
    T1.fo_blank=T1.fo-T1.blank;
    %Calculation of fa_blank
    T1.fa_blank=T1.fa-T1.blank_acid;
    %Calculation of fo_fa
    T1.fo_fa=T1.fo_blank./T1.fa_blank;
    %Calculation of chl
    T1.chl=T1.Fs.*(T1.r./(T1.r-1)).*(T1.fo_blank-T1.fa_blank)...
        .*T1.vol_extracted./T1.vol_filtered;
    %Calculation of phaeo
    T1.phaeo=T1.Fs.*(T1.r./(T1.r-1)).*(T1.r.*T1.fa_blank-T1.fo_blank)...
        .*T1.vol_extracted./T1.vol_filtered;

    %Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    %cruise
    newname=strrep(list(n1).name,'temp-inc-raw','chl-calc');%Replace raw by clean
    newtablename=strcat(rep2,newname);%New tablename and path
    writetable(T1,newtablename)

end