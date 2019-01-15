function csv2mat()

[filename,pathname,~] = uigetfile('*.csv');
[~,name,ext] = fileparts(filename);
filepath = [pathname filename];

[num,txt,~] = xlsread(filepath);

startrow = 3;
rowname = matlab.lang.makeValidName(txt(startrow-1,:), 'replacementstyle', 'delete');


for i = 1:length(rowname)-1
    if contains(rowname{i},'ToolUsed')
        eval([rowname{i} '=txt(startrow:end,i);']);
    else
        eval([rowname{i} '=num(:,i);']);
    end
end

mkdir('mat_data');
savename = ['mat_data/' name '.mat'];
save(savename, rowname{1:end-1});

end