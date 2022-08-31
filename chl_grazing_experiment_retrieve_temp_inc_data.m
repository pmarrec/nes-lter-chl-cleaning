%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matlab script for retrieving the mean/std temperature in the tanks during
%dilution experiments.
%
%For each cast/depth, get the mean/std temperature recorded by the HOBO
%data loggers in the different tanks.
%Because of the lack of records in terms of which bottles in which tanks,
%there are up to 3 CRUISE_temp_inc_X.csv files for each cruise.
% The mean/std in each tank are retrieved and the values with the lowest
% temperature difference will be kept and saved manually.
%From EN687, we are going to keep track of the tanks where the incubations
%were made.
%
%
%Input: CRUSIE-chl-grazing-experiments-ctd-raw.csv files and
%CRUISE_temp_inc_X.csv files
%
%Outputs: CRUISE-chla-grazing-experiments-temp-inc-raw.csv files.
%
%Written by Pierre Marrec
%
%pmarrec@uri.edu
%
%8/30/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\My Drive\NES-LTER_Chla_Cleaning_Rates_Computation\';
%Set the directory where the input raw data are
rep11 = strcat(rep,'chl-grazing-experiment-ctd-raw\');
rep12 = strcat(rep,'chl-grazing-experiment-temp-inc-data\');
%Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-temp-inc-raw\');
%URL of the REST-API
RESTAPI='https://nes-lter-data.whoi.edu/api/ctd/';
%Set the weboptions for downloading files as table
options = weboptions('ContentType', 'table');
%DateTime format
iso8601format = 'yyyy-mm-dd hh:MM:ss';

% Find all the *.csv files
ext='*.csv';
chemin = fullfile(rep11,ext);
list = dir(chemin);%List all files of interest in the directory

for n1=1:numel(list)%skip ar66 because the data are not yet in REST-API
    %load the .csv file of the corresponding cruise
    tablename=strcat(rep11,list(n1).name);
    T1=readtable(tablename);
    cruise=split(list(n1).name,'-');
    CRUISE=cruise{1,1};


    if n1==1%For AR66, don't do anything for the moment, CTD data not accessible in the REST-API yet
        %Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
        %cruise
        newname=strrep(list(n1).name,'ctd-raw','temp-inc-raw');%Replace raw by clean
        newtablename=strcat(rep2,newname);%New tablename and path
        writetable(T1,newtablename)
    else

        %load the CRUISE_temp_inc_X.csv files for the corresponding cruise
        % Find all the *.csv files
        ext='*.csv';
        rep121=strcat(rep12,CRUISE);
        chemin1 = fullfile(rep121,ext);
        list1 = dir(chemin1);%List all files of interest in the directory
        tablename1=strcat(rep121,'\',list1(1).name);
        T21=readtable(tablename1);
        if numel(list1)==2%if 2 tanks
            tablename2=strcat(rep121,'\',list1(2).name);
            T22=readtable(tablename2);
        elseif numel(list1)==3%if 3 tanks
            tablename2=strcat(rep121,'\',list1(2).name);
            T22=readtable(tablename2);
            tablename3=strcat(rep121,'\',list1(3).name);
            T23=readtable(tablename3);
        end

        % Erase the extra " ' " in T1.cast and T1.niskin
        T1.cast=erase(T1.cast,"'");
        T1.niskin=erase(T1.niskin,"'");

        % Identify each unique cast
        a1=unique(T1.cast);

        for n2=1:length(a1)

            b1=strcmp(T1.cast,a1(n2));
            % for each cast, identify the unique sampling depth
            a2=unique(T1.niskin(b1));

            for n3=1:length(a2)
                %Get the index of all the values from a given depth
                b2=b1 & strcmp(T1.niskin,a2(n3));

                %Get the Start/End DateTime of incubations
                B2 = find(b2, 1, 'first');%find the first occurence of b2=1
                START=datenum(T1.date_time_utc_start(B2),'yyyy-mm-dd hh:MM:ss');
                END=datenum(T1.date_time_utc_end(B2),'yyyy-mm-dd hh:MM:ss');

                %Create variables to store the temp mean/std of each tank.
                TEMP_INm=nan(3,1);
                TEMP_INstd=nan(3,1);

                %Find the matching Start/End DateTime in the HOBO files, and
                %get the mean/std Temp in the tanks at this time
                % Up to 3 mean/std Temp if 3 tanks

                %Tank 1
                m11=abs(START-datenum(T21.DateTimeHOBO_UTC));
                m21=abs(END-datenum(T21.DateTimeHOBO_UTC));
                t11=find(m11==min(m11));
                t21=find(m21==min(m21));
                TEMP_INm(1,1)=mean(T21.Temp_Inc(t11:t21),'omitnan');
                TEMP_INstd(1,1)=std(T21.Temp_Inc(t11:t21),'omitnan');

                %if only 2 tanks
                if numel(list1)==2
                    m12=abs(START-datenum(T22.DateTimeHOBO_UTC));
                    m22=abs(END-datenum(T22.DateTimeHOBO_UTC));
                    t12=find(m12==min(m12));
                    t22=find(m22==min(m22));
                    TEMP_INm(2,1)=mean(T22.Temp_Inc(t12:t22),'omitnan');
                    TEMP_INstd(2,1)=std(T22.Temp_Inc(t12:t22),'omitnan');
                elseif numel(list1)==3%if 3 tanks
                    %Tank 2
                    m12=abs(START-datenum(T22.DateTimeHOBO_UTC));
                    m22=abs(END-datenum(T22.DateTimeHOBO_UTC));
                    t12=find(m12==min(m12));
                    t22=find(m22==min(m22));
                    TEMP_INm(2,1)=mean(T22.Temp_Inc(t12:t22),'omitnan');
                    TEMP_INstd(2,1)=std(T22.Temp_Inc(t12:t22),'omitnan');
                    %Tank 3
                    m13=abs(START-datenum(T23.DateTimeHOBO_UTC));
                    m23=abs(END-datenum(T23.DateTimeHOBO_UTC));
                    t13=find(m13==min(m13));
                    t23=find(m23==min(m23));
                    TEMP_INm(3,1)=mean(T23.Temp_Inc(t13:t23),'omitnan');
                    TEMP_INstd(3,1)=std(T23.Temp_Inc(t13:t23),'omitnan');

                end

                %Identify the tank with the lowest temperature difference
                %comapred to the in-situ temperature. This tank will be
                %considered as the one where the incubation was performed and
                %the corresponding mean/std temp will be assigned to this
                %experiment
                Tdiff=abs(TEMP_INm-T1.temperature_sampling(B2));
                [MIN,Idx]=min(Tdiff);
                T1.incubation_tank(b2)=Idx;
                T1.temperature_incubation_avg(b2)=TEMP_INm(Idx);
                T1.temperature_incubation_std(b2)=TEMP_INstd(Idx);


            end


        end

        %Make sure cast and niskin are in text format in the table
        T1.cast=strcat(T1.cast,"'");
        T1.niskin=strcat(T1.niskin,"'");

        %Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
        %cruise
        newname=strrep(list(n1).name,'ctd-raw','temp-inc-raw');%Replace raw by clean
        newtablename=strcat(rep2,newname);%New tablename and path
        writetable(T1,newtablename)

    end

end