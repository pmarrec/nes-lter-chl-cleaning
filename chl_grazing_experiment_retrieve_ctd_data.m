%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for retrieving the DateTime, Temperature, Latitude and Longitude of
% sampling from the CTD bottles.csv files.
%
% Input: CRUSIE-chl-grazing-experiments-raw.csv files
%
% Outputs: CRUISE-chla-grazing-experiments-ctd-raw.csv files.
%
% Written by Pierre Marrec
%
% pmarrec@uri.edu
%
% 3/7/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars, clc, close all

%Set the directory where we work
rep = 'C:\Users\pierr\Desktop\PostDoc_URI_Desktop\NES-LTER\EDI_Growth_Grazing\';
%Set the directory where the input raw data are
rep1 = strcat(rep,'chl-grazing-experiment-raw\');
%Set the directory where the output clean data are
rep2 = strcat(rep,'chl-grazing-experiment-ctd-raw\');
%URL of the REST-API
RESTAPI='https://nes-lter-data.whoi.edu/api/ctd/';
%Set the weboptions for downloading files as table
options = weboptions('ContentType', 'table', 'Timeout', 30);
%DateTime format
iso8601format = 'yyyy-mm-dd hh:MM:ss';

%Find all the *.csv files
ext='*.csv';
chemin = fullfile(rep1,ext);
list = dir(chemin);%List all files of interest in the directory

for n1=1:numel(list)

    %load the .csv file of the corresponding cruise
    tablename=strcat(rep1,list(n1).name);
    T1=readtable(tablename);
    T1.date_time_utc_sampling=string(T1.date_time_utc_sampling);%Convert into string to store the datetime value from the CTD
    cruise=split(list(n1).name,'-');
    CRUISE=cruise{1,1};


    %Download the bottles.csv file for the corresponding cruise
    tablename1 = strcat(RESTAPI,CRUISE,'/bottles.csv');
    T2 = webread(tablename1, options);

    % Erase the extra " ' " in T1.cast and T1.niskin
    T1.cast=erase(T1.cast,"'");
    T1.niskin=erase(T1.niskin,"'");
    T1.niskin_second_cast=erase(T1.niskin_second_cast,"'");

    % Identify each unique cast
    a1=unique(T1.cast);

    for n2=1:length(a1)

        b1=strcmp(T1.cast,a1(n2));
        %Get the cast number from the string
        B1=regexp(a1{n2},'\d*','Match')';

        % for each cast, identify the unique sampling depth
        a2=unique(T1.niskin(b1));
        a3=unique(T1.niskin_second_cast(b1));

        for n3=1:length(a2)
            %Get the index of all the values from a given depth
            b2=b1 & strcmp(T1.niskin,a2(n3));

            %Get the niskin bottle number from the string
            B2 = regexp(a2{n3},'\d*','Match')';

            %Create variable to store the data for each bottle.
            %Max number of bottle = 7
            DateTime=nan(7,1);
            Latitude=nan(7,1);
            Longitude=nan(7,1);
            Depth=nan(7,1);
            Temperature=nan(7,1);

            %Get the 1st (or unique) cast number
            B11=str2double(B1{1,1});

            %1st Bottle
            Bottle11=str2double(B2{1,1});
            A11=find(T2.cast==B11 & T2.niskin==Bottle11);
            DateTime(1,1)=datenum(T2.date(A11),iso8601format);
            Latitude(1,1)=T2.latitude(A11);
            Longitude(1,1)=T2.longitude(A11);
            Depth(1,1)=T2.depsm(A11);
            Temperature(1,1)=T2.t090c(A11);


            %2nd Bottle
            if length(B2)>1%Check if at least 2 bottles
                Bottle12=str2double(B2{2,1});
                A12=find(T2.cast==B11 & T2.niskin==Bottle12);
                DateTime(2,1)=datenum(T2.date(A12),iso8601format);
                Latitude(2,1)=T2.latitude(A12);
                Longitude(2,1)=T2.longitude(A12);
                Depth(2,1)=T2.depsm(A12);
                Temperature(2,1)=T2.t090c(A12);
            else
            end

            %3rd Bottle
            if length(B2)>2%Check if at least 3 bottles
                Bottle13=str2double(B2{3,1});
                A13=find(T2.cast==B11 & T2.niskin==Bottle13);
                DateTime(3,1)=datenum(T2.date(A13),iso8601format);
                Latitude(3,1)=T2.latitude(A13);
                Longitude(3,1)=T2.longitude(A13);
                Depth(3,1)=T2.depsm(A13);
                Temperature(3,1)=T2.t090c(A13);
            else
            end

            %4th Bottle
            if length(B2)>3%Check if at least 4 bottles
                Bottle14=str2double(B2{4,1});
                A14=find(T2.cast==B11 & T2.niskin==Bottle14);
                DateTime(4,1)=datenum(T2.date(A14),iso8601format);
                Latitude(4,1)=T2.latitude(A14);
                Longitude(4,1)=T2.longitude(A14);
                Depth(4,1)=T2.depsm(A14);
                Temperature(4,1)=T2.t090c(A14);
            else
            end

            %5th Bottle
            if length(B2)>4%Check if at least 5 bottles
                Bottle15=str2double(B2{5,1});
                A15=find(T2.cast==B11 & T2.niskin==Bottle15);
                DateTime(5,1)=datenum(T2.date(A15),iso8601format);
                Latitude(5,1)=T2.latitude(A15);
                Longitude(5,1)=T2.longitude(A15);
                Depth(5,1)=T2.depsm(A15);
                Temperature(5,1)=T2.t090c(A15);
            else
            end


            %Check if 2 cast were made for a given station
            if length(B1)>1
                %Get the 2nd cast number
                B12=str2double(B1{2,1});
                %Get the niskin_other bottle number from the string
                B3 = regexp(a3{n3},'\d*','Match')';

                %1st Bottle during the 2nd cast
                Bottle21=str2double(B3{1,1});
                A21=find(T2.cast==B12 & T2.niskin==Bottle21);
                DateTime(6,1)=datenum(T2.date(A21),iso8601format);
                Latitude(6,1)=T2.latitude(A21);
                Longitude(6,1)=T2.longitude(A21);
                Depth(6,1)=T2.depsm(A21);
                Temperature(6,1)=T2.t090c(A21);

                %2nd Bottle other
                if length(B3)>1%Check if at least 2 bottles during the 2nd cast
                    Bottle22=str2double(B3{2,1});
                    A22=find(T2.cast==B12 & T2.niskin==Bottle22);
                    DateTime(7,1)=datenum(T2.date(A22),iso8601format);
                    Latitude(7,1)=T2.latitude(A22);
                    Longitude(7,1)=T2.longitude(A22);
                    Depth(7,1)=T2.depsm(A22);
                    Temperature(7,1)=T2.t090c(A22);
                else
                end

            else
            end

            T1.date_time_utc_sampling(b2)=datestr(mean(DateTime,'omitnan'),'yyyy-mm-dd HH:MM:ss+00:00');
            T1.latitude(b2)=mean(Latitude,'omitnan');
            T1.longitude(b2)=mean(Longitude,'omitnan');
            T1.depth(b2)=mean(Depth,'omitnan');
            T1.temperature_sampling(b2)=mean(Temperature,'omitnan');

            clear Temperature Depth Latitude Longitude DateTime

        end


    end

    %Make sure cast and niskin are in text format in the table
    T1.cast=strcat(T1.cast,"'");
    T1.niskin=strcat(T1.niskin,"'");
    T1.niskin_second_cast=strcat(T1.niskin_second_cast,"'");

    %Save the new CRUISE-chla-grazing-experiments-clean.csv files for each
    %cruise
    newname=strrep(list(n1).name,'raw','ctd-raw');%Replace raw by clean
    newtablename=strcat(rep2,newname);%New tablename and path
    writetable(T1,newtablename)



end

