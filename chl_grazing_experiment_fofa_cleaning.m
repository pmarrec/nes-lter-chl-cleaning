%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for Chl-a conc cleaning of the raw Chl-a
% obtained during dilution (grazing) experiments based on Fo/Fa ratios.
%
% Fo/Fa cleaning:
% 2 criteria, for each cast/depth:
% 1) within 1-3 range,
% 2) witihn +/- 2 StdDev confidence interval for a given type of filter
% GFF after screening with 200um mesh = >0&<200um
% 10um after screening with 200um mesh = >10&<200um
% GFF without screening with 200um mesh = >0 (for EN627 L11-B)
% GFF with screening with 10um mesh = >0&<10um (for EN668)
% All the values that don't fit these criteria are flagged as questionable
% with a iode_quality_flag = 3
% All the other values are flagged as good with a iode_quality_flag = 1
%
% Input: CRUSIE-chl-grazing-experiments-chl-calc.csv files with Chl-a
%        calculated from the chl_grazing_experiment_chl_calc.m script
%
% Outputs: CRUISE-chla-grazing-experiments-fofa-clean.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 3/7/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\EDI_Growth_Grazing\DataPackage_GFF_and_10um\';
%Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-chl-calc\');
%Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-fofa-clean\');

%Find all the *.cnv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);%List all files of interest in the directory

for n1=1:numel(list)
    %load the .csv file of the corresponding cruise
    tablename=strcat(rep1,list(n1).name);
    T1=readtable(tablename);

    %Erase the extra " ' " in T1.cast and T1.niskin
    T1.cast=erase(T1.cast,"'");
    T1.niskin=erase(T1.niskin,"'");

    %identify each unique cast
    a1=unique(T1.cast);
    %find the rows corresponding to the corresponding cast
    for n2=1:length(a1)
        b1=strcmp(T1.cast,a1(n2));
        %for each cast, identify the unique sampling depth
        a2=unique(T1.niskin(b1));
        for n3=1:length(a2)
            b2=b1 & strcmp(T1.niskin,a2(n3));

            %for each sampling depth, indetify the T0/TF
            a3=unique(T1.T0_TF(b2));

            for n4=1:length(a3)
                b3=b2 & strcmp(T1.T0_TF,a3(n4));
                %for each sampling depth, indetify the unique filter sizes
                a4=unique(T1.filter_size(b3));

                for n5=1:length(a4)
                    %Get the index of all the values from a given size
                    %fraction
                    b4=b3 & strcmp(T1.filter_size,a4(n5));

                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    % Fo/Fa cleaning
                    % 2 criteria, for each cast/depth:
                    % 1) within 1-3 range,
                    % 2) witihn +/- 2 StdDev confidence interval for a given type of filter
                    % GFF after screening with 200um mesh = >0&<200um
                    % 10um after screening with 200um mesh = >10&<200um
                    %%%%%%%%%%%%%%%%%%%%%%%%%

                    %1st step, QC based on FoFa ratios
                    FoFa=T1.fo_fa(b4);
                    %Discard all values 1<FoFa<3
                    FoFa(FoFa>3)=[];
                    FoFa(FoFa<1)=[];
                    FoFa(isnan(FoFa))=[];
                    %Get the mean/stddev and the upper and lower limits of the
                    %confidence interval
                    FoFa_avg=mean(FoFa);
                    FoFa_std=std(FoFa);
                    FoFa_ulim=FoFa_avg+2*FoFa_std;
                    FoFa_llim=FoFa_avg-2*FoFa_std;

                    %Assigned a iode_quality flag (1=good, 3=questionable/suspect)
                    %to the data based on the nan values, the threshold (1<x<3)
                    % and the upper/lower limits defined (ulim and llim)
                    for n6=1:length(b4)
                        if b4(n6)==1
                            if (isnan(T1.fo_fa(n6))) || (T1.fo_fa(n6)>3) || (T1.fo_fa(n6)<1) || (T1.fo_fa(n6)>FoFa_ulim) || (T1.fo_fa(n6)<FoFa_llim)
                                T1.iode_quality_flag(n6)=3;
                            else
                                T1.iode_quality_flag(n6)=1;
                            end
                        end
                    end


                    clear FoFa FoFa_avg FoFa_std FoFa_ulim FoFa_llim


                end
            end
        end


    end

    %Make sure cast and niskin are in text format in the table
    T1.cast=strcat(T1.cast,"'");
    T1.niskin=strcat(T1.niskin,"'");

    %Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    %cruise
    newname=strrep(list(n1).name,'chl-calc','fofa-clean');%Replace raw by clean
    newtablename=strcat(rep2,newname);%New tablename and path

    writetable(T1,newtablename)

end


