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
    filepath = [subdir(i).folder '\' subdir(i).name];

    ins = readtable(filepath);

    Instr = table2struct(ins,'ToScalar',true);

    fields = fieldnames(Instr);

    for idx=1:length(fieldnames(Instr))
        disp(idx);
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
    save(['Data_' loadfolder '.mat'],'RawData_1M')
elseif strcmp(loadfolder,'2J_Compact')==1
    RawData_2J = RawData;
    save(['Data_' loadfolder '.mat'],'RawData_2J')
elseif strcmp(loadfolder,'3S_Compact')==1
    RawData_3S = RawData;
    save(['Data_' loadfolder '.mat'],'RawData_3S')
elseif strcmp(loadfolder,'4C_Compact')==1
    RawData_4C = RawData;
    save(['Data_' loadfolder '.mat'],'RawData_4C')
end

cd(old);