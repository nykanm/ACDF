% Use this file to concatenate all your csv files from the OSSim into one
% big csv file for all tools used.

% Then type "csv2struct" in your command line and load the file "tools.csv"

% Change this to the directory where your ToolData files are:
loaddir = '/Users/nykan/Documents/McGill Grad/Matlab/ACDF/Data/ACDF85/ACDF_20181212105246/ToolLog/';
% Change this to the directory where you want to save tools.csv
savedir = '/Users/nykan/Documents/McGill Grad/Matlab/ACDF/Data/ACDF85/ACDF_20181212105246/ToolLog/';

subdir = dir(loaddir);

for i=3:length(subdir)
    filename = [loaddir '/' subdir(i).name];
    csv{i-2} = readtable(filename);
end

allcsv = vertcat(csv{1},csv{2});
for o=3:length(csv)
    allcsv = vertcat(allcsv,csv{o});
end

allcsv = sortrows(allcsv,'TimeSinceStart');

writetable(allcsv,[savedir 'tools.csv']);

