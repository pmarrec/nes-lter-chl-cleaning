%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for Chl-a conc cleaning of the raw Chl-a
% obtained during dilution (grazing) experiments.
%
% Chl-a conc cleaning:
% All negative Chl-a conc are flagged as questionnable/suspect with a
% iode_quality_flag = 3
% Criteria, for each station/depth/treatment/triplicate values:
% Each triplicate value should stand in the +/- 2 x %CV of
% the mean values of the triplicate with a QC flag=1
% %CV is considered as the mean %CV obtained on a given type of
% filter (GFF/10um) at T0 and at TF.
% All the values that don't fit these criteria are flagged as questionable
% with a iode_quality_flag = 3
% All the other values are flagged as good with a iode_quality_flag = 1
%
% Input: CRUSIE-chl-grazing-experiments-fofa-clean.csv files 
%
% Outputs: CRUISE-chla-grazing-experiments-chl-conc-clean.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 3/7/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\EDI_Growth_Grazing\DataPackage_GFF_10um\';
%Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-fofa-clean\');
%Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-chl-conc-clean\');

%Find all the *.cnv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);%List all files of interest in the directory

for n1=1:numel(list)
    %load the .csv file of the corresponding cruise
    tablename=strcat(rep1,list(n1).name);
    T1=readtable(tablename);

    %Flag as questionable/suspect all <0 chl conc value
    for n=1:height(T1)
        if T1.chl(n)<0
           T1.iode_quality_flag(n)=3; 
        end
    end

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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %Chl-a conc cleaning
            % Criteria, for each station/depth/triplicate values:
            % Each triplicate value should stand in the +/- 2 x %CV of
            % the mean values of the triplicate with a QC flag=1
            % %CV is considered as the mean %CV obtained on a given type of
            % filter (GFF/10um) at T0 and at TF.
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %T0 and >0&<200 filters (GFF)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Identify all values obtained with >0&<200 filters at T0
            c1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'T0');

            %Identify as group the dil and wsw values
            [d1,D1]=findgroups(T1.dilution(c1));
                        
            %Create a new vector to store the mean values
            chl_avg=nan(max(d1),1);
            %Create a new vector to store the stdev values
            chl_std=nan(max(d1),1);
            %Create a new vector to store the %CV values
            chl_cv=nan(max(d1),1);

            for m1=1:max(d1)
                %Identify the values of the given group with a
                %iode_quality_flag of 1
                c2=c1 & strcmp(T1.dilution,D1(m1)) & T1.iode_quality_flag==1;
                chl_avg(m1)=mean(T1.chl(c2));
                chl_std(m1)=std(T1.chl(c2));
                chl_cv(m1)=chl_std(m1)/chl_avg(m1);
            end

            %Define the lower limit values
            chl_llim=chl_avg-2*mean(chl_cv)*chl_avg;
            %Define the upper limit values
            chl_ulim=chl_avg+2*mean(chl_cv)*chl_avg;

            %Check if the values for each type of treatment are in the
            %confidence interval llim<x<ulim. If they are not, values are
            %flag as questionable/suspect
            for m1=1:max(d1)
                c2=c1 & strcmp(T1.dilution,D1(m1)) & T1.iode_quality_flag==1;
                for M1=1:length(c2)
                    if c2(M1)==1
                        if (T1.chl(M1)<chl_llim(m1)) || (T1.chl(M1)>chl_ulim(m1))
                            T1.iode_quality_flag(M1)=3;
                        end
                    end
                end
            end

            clear chl_avg chl_std chl_llim chl_ulim chl_cv

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %T0 and >10&<200 filters (10um)
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            %Identify all values obtained with >10&<200 filters at T0
            c1=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'T0');

            %Identify as group the dil and wsw values
            [d1,D1]=findgroups(T1.dilution(c1));
                        
            %Create a new vector to store the mean values
            chl_avg=nan(max(d1),1);
            %Create a new vector to store the stdev values
            chl_std=nan(max(d1),1);
            %Create a new vector to store the %CV values
            chl_cv=nan(max(d1),1);

            for m1=1:max(d1)
                %Identify the values of the given group with a
                %iode_quality_flag of 1
                c2=c1 & strcmp(T1.dilution,D1(m1)) & T1.iode_quality_flag==1;
                chl_avg(m1)=mean(T1.chl(c2));
                chl_std(m1)=std(T1.chl(c2));
                chl_cv(m1)=chl_std(m1)/chl_avg(m1);
            end

            %Define the lower limit values
            chl_llim=chl_avg-2*mean(chl_cv)*chl_avg;
            %Define the upper limit values
            chl_ulim=chl_avg+2*mean(chl_cv)*chl_avg;

            %Check if the values for each type of treatment are in the
            %confidence interval llim<x<ulim. If they are not, values are
            %flag as questionable/suspect
            for m1=1:max(d1)
                c2=c1 & strcmp(T1.dilution,D1(m1)) & T1.iode_quality_flag==1;
                for M1=1:length(c2)
                    if c2(M1)==1
                        if (T1.chl(M1)<chl_llim(m1)) || (T1.chl(M1)>chl_ulim(m1))
                            T1.iode_quality_flag(M1)=3;
                        end
                    end
                end
            end
            clear chl_avg chl_std chl_llim chl_ulim chl_cv

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %TF and >0&<200 filters (GFF)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Identify all values obtained with >0&<200 filters at TF
            c1=b2 & strcmp(T1.filter_size,'>0&<200') & strcmp(T1.T0_TF,'TF');

            %Identify as group based on dilution, nutrient_treatment, light
            %level and replicate_bottle
            [d1,D1,D2,D3,D4]=findgroups(T1.dilution(c1),T1.nutrient_treatment(c1),T1.light_level(c1),T1.replicate_bottle(c1));
                        
            %Create a new vector to store the mean values
            chl_avg=nan(max(d1),1);
            %Create a new vector to store the stdev values
            chl_std=nan(max(d1),1);
            %Create a new vector to store the %CV values
            chl_cv=nan(max(d1),1);

            for m1=1:max(d1)
                %Identify the values of the given group with a
                %iode_quality_flag of 1
                c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
                    & strcmp(T1.light_level,D3(m1)) & strcmp(T1.replicate_bottle,D4(m1))...
                    & T1.iode_quality_flag==1;
                chl_avg(m1)=mean(T1.chl(c2));
                chl_std(m1)=std(T1.chl(c2));
                chl_cv(m1)=chl_std(m1)/chl_avg(m1);
            end

            %Define the lower limit values
            chl_llim=chl_avg-2*mean(chl_cv)*chl_avg;
            %Define the upper limit values
            chl_ulim=chl_avg+2*mean(chl_cv)*chl_avg;

            %Check if the values for each type of treatment are in the
            %confidence interval llim<x<ulim. If they are not, values are
            %flag as questionable/suspect
            for m1=1:max(d1)
                c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
                    & strcmp(T1.light_level,D3(m1)) & strcmp(T1.replicate_bottle,D4(m1))...
                    & T1.iode_quality_flag==1;
                for M1=1:length(c2)
                    if c2(M1)==1
                        if (T1.chl(M1)<chl_llim(m1)) || (T1.chl(M1)>chl_ulim(m1))
                            T1.iode_quality_flag(M1)=3;
                        end
                    end
                end
            end

            clear chl_avg chl_std chl_llim chl_ulim chl_cv

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %TF and >10&<200 filters (10um)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Identify all values obtained with >10&<200 filters at TF
            c1=b2 & strcmp(T1.filter_size,'>10&<200') & strcmp(T1.T0_TF,'TF');

            %Identify as group based on dilution, nutrient_treatment, light
            %level and replicate_bottle
            [d1,D1,D2,D3,D4]=findgroups(T1.dilution(c1),T1.nutrient_treatment(c1),T1.light_level(c1),T1.replicate_bottle(c1));
                        
            %Create a new vector to store the mean values
            chl_avg=nan(max(d1),1);
            %Create a new vector to store the stdev values
            chl_std=nan(max(d1),1);
            %Create a new vector to store the %CV values
            chl_cv=nan(max(d1),1);

            for m1=1:max(d1)
                %Identify the values of the given group with a
                %iode_quality_flag of 1
                c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
                    & strcmp(T1.light_level,D3(m1)) & strcmp(T1.replicate_bottle,D4(m1))...
                    & T1.iode_quality_flag==1;
                chl_avg(m1)=mean(T1.chl(c2));
                chl_std(m1)=std(T1.chl(c2));
                chl_cv(m1)=chl_std(m1)/chl_avg(m1);
            end

            %Define the lower limit values
            chl_llim=chl_avg-2*mean(chl_cv)*chl_avg;
            %Define the upper limit values
            chl_ulim=chl_avg+2*mean(chl_cv)*chl_avg;

            %Check if the values for each type of treatment are in the
            %confidence interval llim<x<ulim. If they are not, values are
            %flag as questionable/suspect
            for m1=1:max(d1)
                c2=c1 & strcmp(T1.dilution,D1(m1)) & strcmp(T1.nutrient_treatment,D2(m1)) ...
                    & strcmp(T1.light_level,D3(m1)) & strcmp(T1.replicate_bottle,D4(m1))...
                    & T1.iode_quality_flag==1;
                for M1=1:length(c2)
                    if c2(M1)==1
                        if (T1.chl(M1)<chl_llim(m1)) || (T1.chl(M1)>chl_ulim(m1))
                            T1.iode_quality_flag(M1)=3;
                        end
                    end
                end
            end

            clear chl_avg chl_std chl_llim chl_ulim chl_cv

            
        end


    end

    %Make sure cast and niskin are in text format in the table
    T1.cast=strcat(T1.cast,"'");
    T1.niskin=strcat(T1.niskin,"'");

    %Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    %cruise
    newname=strrep(list(n1).name,'fofa','chl-conc');%Replace raw by clean
    newtablename=strcat(rep2,newname);%New tablename and path
    writetable(T1,newtablename)

end


