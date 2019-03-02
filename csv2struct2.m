function [RawData] = csv2struct2(loadfolder,DataDir)
% Function converts csv files into a structure. 
% csv2struct2(loadfolder) where loadfolder is '1M_Compact', '2J_Compact',
% '3S_Compact', or '4C_Compact'
loaddir = [DataDir loadfolder];

subdir = dir(loaddir);

if strcmp(subdir(3).name,'.DS_Store')==1
    n=4;
else
    n=3;
end

for i=n:length(subdir) 
    disp([num2str(i),' of ',num2str(length(subdir))]);
    filepath = [subdir(i).folder '/' subdir(i).name];

    ins = readtable(filepath);

    Instr = table2struct(ins,'ToScalar',true);

    fields = fieldnames(Instr);

    for idx=1:length(fieldnames(Instr))
        if n==3
            RawData(i-2).(fields{idx}) = vertcat(Instr.(fields{idx}));
        elseif n==4
            RawData(i-3).(fields{idx}) = vertcat(Instr.(fields{idx}));
        end
    end
end
old = cd;
cd(DataDir);

if strcmp(loadfolder,'1M_Compact')==1
    RawData_1M = RawData;
    disp('Saving 1M....');
    %save(['Data_' loadfolder '.mat'],'RawData_1M')
elseif strcmp(loadfolder,'2J_Compact')==1
    RawData_2J = RawData;
    disp('Saving 2J....');
   % save(['Data_' loadfolder '.mat'],'RawData_2J')
elseif strcmp(loadfolder,'3S_Compact')==1
    RawData_3S = RawData;
    disp('Saving 3S....');
   % save(['Data_' loadfolder '.mat'],'RawData_3S')
elseif strcmp(loadfolder,'4F_Compact')==1
    RawData_4F = RawData;
    disp('Saving 4F....');
   % save(['Data_' loadfolder '.mat'],'RawData_4F')
elseif strcmp(loadfolder,'5C_Compact')==1
    RawData_5C = RawData;
    disp('Saving 5C....');
   % save(['Data_' loadfolder '.mat'],'RawData_5C')
end

cd(old);