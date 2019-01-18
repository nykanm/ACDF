function [] = csv_concat(loadfolder,DataDir)
% This function concatenates all your participants' toollogs into one
% toollog for each participant saved in CompactTools folder. To use, type
% csv_concat(loadfolder) where loadfolder is '1M' or '2J' or '3R' or '4C'

% Use this file to concatenate all your csv files from the OSSim into one
% big csv file for all tools used.

% Then type "csv2struct" in your command line and load the file "tools.csv"
% Change this to the directory where your ToolData files are:
loaddir = [DataDir loadfolder];
savedir = ['/Users/nykan/Documents/McGill Grad/Matlab/ACDF/Sample Data/'];


%% csv_concat
subdir = dir(loaddir);

for i=3:length(subdir)
    disp(i-2);
    foldername = [loaddir '/' subdir(i).name];
    subdir2 = dir(foldername);
    if length(subdir2)==3
        foldername2 = [foldername '/' subdir2(3).name '/ToolLog'];
    elseif length(subdir2)==4
        foldername2 = [foldername '/' subdir2(4).name '/ToolLog'];
    end
        subdir3 = dir(foldername2);
    for ii=3:length(subdir3)
        disp(ii-2);
        filename = [foldername2 '/' subdir3(ii).name];
        csv{ii-2} = readtable(filename);
    end
        allcsv = vertcat(csv{1},csv{2});
        for o=3:length(csv)
            allcsv = vertcat(allcsv,csv{o});
        end

        allcsv = sortrows(allcsv,'TimeSinceStart');
        
        writetable(allcsv,[savedir loadfolder '_Compact/' subdir(i).name '.csv']);
        clear csv;
        clear allcsv;
    
end





