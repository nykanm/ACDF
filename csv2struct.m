function [Instrument] = csv2struct()
% Function converts csv files into a structure. 
% Use this following format: [name_of_structure] = csv2struct().
[filename,pathname,~] = uigetfile('*.csv');
filepath = [pathname filename];

ins = readtable(filepath);

Instr = table2struct(ins,'ToScalar',true);

fields = fieldnames(Instr);
    
for idx=1:length(fieldnames(Instr))
    Instrument.(fields{idx}) = vertcat(Instr.(fields{idx}));
end
