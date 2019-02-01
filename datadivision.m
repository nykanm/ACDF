%% function [sub,sig_table,sig_vars,nov,exp,novmat,expmat] = metrics_ACDF(data_nov,data_exp,sig_level)
% 1. Run csv_concat() to concatenate all tool files into one for
% '1M','2J', '3S', '4C'
% 2. Run [Instrument] = csv2struct; to convert into structure.

%% Data Conversions (only do this once)
% Take all the data from Ossim and put them into a folder called Sample
% Data. Create a folder for each group such as 1M, 2J, 3S, 4C
% Change this to the directory where your raw data is stored:
DataDir = 'F:\OSSIM STUDY DATA\';
 csv_concat('1M',DataDir);
 csv_concat('2J',DataDir);
 csv_concat('3S',DataDir);
 csv_concat('4C',DataDir);
% This creates the .mat files for all groups in the same folder as the raw
% data.
csv2struct2('1M_Compact',DataDir);
csv2struct2('2J_Compact',DataDir);
csv2struct2('3S_Compact',DataDir);
csv2struct2('4C_Compact',DataDir);

%% Load your data
DataDir = 'F:\OSSIM STUDY DATA\';

load([DataDir 'Data_1M_Compact.mat']);
load([DataDir 'Data_2J_Compact.mat']);
load([DataDir 'Data_3S_Compact.mat']);
load([DataDir 'Data_4C_Compact.mat']);

%% Time Blocks
% set timB equal to the time block in seconds
timB = 10;

for sj=1:length(data)
    disp(['Subject: ',num2str(sj),' of ',num2str(length(data))]);
    y=data(sj);

%% Find time duplicates
timeDupInd = find([0.1;diff(y.TimeSinceStart)]); 
fieldNam = fieldnames(y);
for i = 1:numel(fieldNam)
    x.(fieldNam{i}) = y.(fieldNam{i})(timeDupInd);
end


ToolUsed.BoneCurette = find(contains(x.ToolUsed,'Bone Curette'));

% Check for randomly changing tools
Chang = find(diff(ToolUsed.BoneCurette)>1);
ChangIdx = ToolUsed.BoneCurette(Chang);
Chang = flip(Chang);
for i = 1:length(ChangIdx)
    if contains(x.ToolUsed(ChangIdx(1)-10),'Bone Curette') || contains(x.ToolUsed(ChangIdx(1)+10),'Bone Curette')
    else
        ToolUsed.BoneCurette(Chang(i))=[];
    end
end
clear Chang

%% Find the start and end of each 10sec block

if any(ToolUsed.BoneCurette)>0
% For 10sec blocks
% TimeStart = min(x.TimeSinceStart(ToolUsed.BoneCurette));
% TimeEnd = max(x.TimeSinceStart(ToolUsed.BoneCurette));
% Diff = TimeEnd-TimeStart;
% groupStartIdx(1) = find(x.TimeSinceStart==TimeStart);
% correction=0;
%     for i=0:timB:Diff
%         if i==0
%             TimeLoop = TimeStart+timB;
%             groupEnd = find(x.TimeSinceStart>TimeLoop-1 & x.TimeSinceStart<TimeLoop+1);
%             groupEndIdx(1) = groupEnd(1);
%         else
%             groupStartIdx(i/timB+1-correction) = groupEndIdx(i/timB-correction);
%             TimeLoop = x.TimeSinceStart(groupEndIdx(i/timB-correction)) + timB;
%             groupEnd = find(x.TimeSinceStart>TimeLoop-1 & x.TimeSinceStart<TimeLoop+1);
%             if any(groupEnd)>0
%                 groupEndIdx(i/timB+1-correction) = groupEnd(1);
%             else
%                 groupStartIdx(i/timB+1-correction) = [];
%                 correction = correction + timB;
%             end
%         end
%     end
%     if length(groupStartIdx)>length(groupEndIdx)
%     groupStartIdx(end) = [];
%     end
% % -----------------------    
% For blocks that increase by 10sec every itteration
TimeStart = min(x.TimeSinceStart(ToolUsed.BoneCurette));
TimeEnd = max(x.TimeSinceStart(ToolUsed.BoneCurette));
Diff = TimeEnd-TimeStart;
groupStartIdx(1) = find(x.TimeSinceStart==TimeStart);
correction=0;
rep=Diff/timB;rep=round(rep);
    for i=1:rep
        if i==1
            TimeLoop = TimeStart+timB;
            groupEnd = find(x.TimeSinceStart>TimeLoop-1 & x.TimeSinceStart<TimeLoop+1);
            groupEndIdx(1) = groupEnd(1);
        else
            TimeLoop = TimeStart+(timB*i);
            groupEnd = find(x.TimeSinceStart>TimeLoop-1 & x.TimeSinceStart<TimeLoop+1);
            if any(groupEnd)>0
                groupEndIdx(i) = groupEnd(1);
            else
            end
        end
    end
for l = 1:length(groupEndIdx)
    jumps=linspace(groupStartIdx(1),groupEndIdx(l),(groupEndIdx(l)-groupStartIdx(1))+1);    
    for i = 1:numel(fieldNam)
        M1(sj).Block_Curette(l).(fieldNam{i}) = x.(fieldNam{i})(jumps');
        datTab{sj} = struct2table(M1(sj).Block_Curette);
    end
end
%% Extract Data in the blocks
% for l = 1:length(groupStartIdx)
%     jumps=linspace(groupStartIdx(l),groupEndIdx(l),(groupEndIdx(l)-groupStartIdx(l))+1);    
%     for i = 1:numel(fieldNam)
%         M1(sj).Block_Curette(l).(fieldNam{i}) = x.(fieldNam{i})(jumps');
%         datTab{sj} = struct2table(M1(sj).Block_Curette);
%     end
% end
% -----------
if sj==length(data)
    AllData = [datTab{1}];
    if sj>1
    for ss=2:sj
        AllData=[AllData;datTab{ss}];
    end
    end
end
clear groupEndIdx;
end
end




%% Metrics
for d=2
    
    if d==1
        %data=RawData_1M;
        data=table2struct(AllData_1M);
        disp(['Med Students Total: ',num2str(length(data))]);
    elseif d==2
        %data=RawData_2J;
        data=table2struct(AllData_2J);
        disp(['Junior Residents Total: ',num2str(length(data))]);
    elseif d==3
        %data=RawData_3S;
        data=table2struct(AllData_3S);
        disp(['Senior Residents Total: ',num2str(length(data))]);
    elseif d==4
        %data=RawData_4C;
        data=table2struct(AllData_4C);
        disp(['Consultant Total: ',num2str(length(data))]);
    end
for sj=1:length(data)
    tic
    disp(['Subject: ',num2str(sj),' of ',num2str(length(data))]);
    y=data(sj);
if length(y.TimeSinceStart)>0

%% Find time duplicates
timeDupInd = find([0.1;diff(y.TimeSinceStart)]); 
fieldNam = fieldnames(y);
for i = 1:numel(fieldNam)
    x.(fieldNam{i}) = y.(fieldNam{i})(timeDupInd);
end

%% Define some FINDS
% when in contact with structure
Contact.C4 = find(x.ContactVoxelsC4Vertebra);
Contact.C5 = find(x.ContactVoxelsC5Vertebra);
Contact.C4C5DiscAnnulus = find(x.ContactVoxelsC4C5DiscAnnulus);
Contact.C4C5DiscNucleus = find(x.ContactVoxelsC4C5DiscNucleus);
Contact.PllLeft = find(x.ContactVoxelsPllLeftBeam);
Contact.PllRight = find(x.ContactVoxelsPllRightBeam);
Contact.SC = find(x.ContactVoxelsSpinalCordNerves);
Contact.LVA = find(x.ContactVoxelsLeftVertebralArtery);
Contact.RVA = find(x.ContactVoxelsRightVertebralArtery);
% when cutting a structure
Cutting.C4 = find(x.CutVoxelsC4Vertebra);
Cutting.C5 = find(x.CutVoxelsC5Vertebra);
Cutting.C4C5DiscAnnulus = find(x.CutVoxelsC4C5DiscAnnulus);
Cutting.C4C5DiscNucleus = find(x.CutVoxelsC4C5DiscNucleus);
Cutting.PllLeft = find(x.CutVoxelsPllLeftBeam);
Cutting.PllRight = find(x.CutVoxelsPllRightBeam);
Cutting.SC = find(x.CutVoxelsSpinalCordNerves);
Cutting.LVA = find(x.CutVoxelsLeftVertebralArtery);
Cutting.RVA = find(x.CutVoxelsRightVertebralArtery);
% when each tool is used
ToolUsed.Scalpel = find(contains(x.ToolUsed,'Scalpel'));
ToolUsed.BoneCurette = find(contains(x.ToolUsed,'Bone Curette'));
ToolUsed.Rongeur2mm = find(contains(x.ToolUsed,'Pituitary Rongeur 2mm'));
ToolUsed.DiscRongeur = find(contains(x.ToolUsed,'Disc Rongeur'));
ToolUsed.Burr = find(contains(x.ToolUsed,'Burr'));
ToolUsed.NerveHook = find(contains(x.ToolUsed,'Nerve Hook'));
ToolUsed.Kerrison1mm = find(contains(x.ToolUsed,'Kerrison 1mm'));

%% Conversions
% Convert rows when in contact to 1 instead of number of voxels
Contact_con1.C4 = x.ContactVoxelsC4Vertebra; Contact_con1.C4(Contact.C4) = 1;
Contact_con1.C5 = x.ContactVoxelsC5Vertebra; Contact_con1.C5(Contact.C5) = 1;
Contact_con1.C4C5DiscAnnulus = x.ContactVoxelsC4C5DiscAnnulus; Contact_con1.C4C5DiscAnnulus(Contact.C4C5DiscAnnulus) = 1;
Contact_con1.C4C5DiscNucleus = x.ContactVoxelsC4C5DiscNucleus; Contact_con1.C4C5DiscNucleus(Contact.C4C5DiscNucleus) = 1;
Contact_con1.PllLeft = x.ContactVoxelsPllLeftBeam; Contact_con1.PllLeft(Contact.PllLeft) = 1;
Contact_con1.PllRight = x.ContactVoxelsPllRightBeam; Contact_con1.PllRight(Contact.PllRight) = 1;
Contact_con1.SC = x.ContactVoxelsSpinalCordNerves; Contact_con1.SC(Contact.SC) = 1;
Contact_con1.LVA = x.ContactVoxelsLeftVertebralArtery; Contact_con1.LVA(Contact.LVA) = 1;
Contact_con1.RVA = x.ContactVoxelsRightVertebralArtery; Contact_con1.RVA(Contact.RVA) = 1;
% Convert rows when cutting to 1 instead of number of voxels
anat_cutvox = {'CutVoxelsC4Vertebra','CutVoxelsC5Vertebra','CutVoxelsC4C5DiscAnnulus',...
    'CutVoxelsC4C5DiscNucleus','CutVoxelsPllLeftBeam','CutVoxelsPllRightBeam',...
    'CutVoxelsSpinalCordNerves','CutVoxelsLeftVertebralArtery','CutVoxelsRightVertebralArtery'};
Cutting_con1.C4 = x.CutVoxelsC4Vertebra; Cutting_con1.C4(Cutting.C4) = 1;
Cutting_con1.C5 = x.CutVoxelsC5Vertebra; Cutting_con1.C5(Cutting.C5) = 1;
Cutting_con1.C4C5DiscAnnulus = x.CutVoxelsC4C5DiscAnnulus; Cutting_con1.C4C5DiscAnnulus(Cutting.C4C5DiscAnnulus) = 1;
Cutting_con1.C4C5DiscNucleus = x.CutVoxelsC4C5DiscNucleus; Cutting_con1.C4C5DiscNucleus(Cutting.C4C5DiscNucleus) = 1;
Cutting_con1.PllLeft = x.CutVoxelsPllLeftBeam; Cutting_con1.PllLeft(Cutting.PllLeft) = 1;
Cutting_con1.PllRight = x.CutVoxelsPllRightBeam; Cutting_con1.PllRight(Cutting.PllRight) = 1;
Cutting_con1.SC = x.CutVoxelsSpinalCordNerves; Cutting_con1.SC(Cutting.SC) = 1;
Cutting_con1.LVA = x.CutVoxelsLeftVertebralArtery; Cutting_con1.LVA(Cutting.LVA) = 1;
Cutting_con1.RVA = x.CutVoxelsRightVertebralArtery; Cutting_con1.RVA(Cutting.RVA) = 1;

%% Define anat and tools
%tools = {'Scalpel','BoneCurette','Rongeur2mm','DiscRongeur','Burr','NerveHook','Kerrison1mm'};
tools = {'BoneCurette'};
anat = {'C4','C5','C4C5DiscAnnulus','C4C5DiscNucleus','PllLeft','PllRight','SC','LVA','RVA'};

%% Number of touches (instenses)

% Overall
for a = 1:length(anat)
    Metrics(sj).(strcat((anat{a}),'Overall','Contact_num')) = numel(find(diff(Contact_con1.(anat{a}))>0));
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Contact_1.(anat{a}).(tools{t}) = Contact_con1.(anat{a})(ToolUsed.(tools{t}));  
        Metrics(sj).(strcat((anat{a}),(tools{t}),'Contact_num')) = numel(find(diff(Contact_1.(anat{a}).(tools{t}))));
    end
end

%% Number of cuts (instenses)

% Overall
for a = 1:length(anat)
    Metrics(sj).(strcat((anat{a}),'Overall','Cutting_num')) = numel(find(diff(Cutting_con1.(anat{a}))>0));
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Cutting_1.(anat{a}).(tools{t}) = Cutting_con1.(anat{a})(ToolUsed.(tools{t}));  
        Metrics(sj).(strcat((anat{a}),(tools{t}),'Cutting_num')) = numel(find(diff(Cutting_1.(anat{a}).(tools{t}))));
    end
end

%% Number of voxels cut

% Overall
for a = 1:length(anat_cutvox)
    Metrics(sj).(strcat((anat{a}),'Overall','Cut_voxels')) = sum(x.(anat_cutvox{a}));
end

% With each tool
for a = 1:length(anat_cutvox)
    for t = 1:length(tools)
        Metrics(sj).(strcat((anat{a}),(tools{t}),'Cut_voxels')) = sum(x.(anat_cutvox{a})(ToolUsed.(tools{t})));  
    end
end

%% Average force
anat_forc_x = {'AverageForceC4Vertebra_X','AverageForceC5Vertebra_X','AverageForceC4C5DiscAnnulus_X',...
    'AverageForceC4C5DiscNucleus_X','AverageForcePllLeftBeam_X','AverageForcePllRightBeam_X',...
    'AverageForceSpinalCordNerves_X','AverageForceLeftVertebralArtery_X','AverageForceRightVertebralArtery_X'};
anat_forc_y = {'AverageForceC4Vertebra_Y','AverageForceC5Vertebra_Y','AverageForceC4C5DiscAnnulus_Y',...
    'AverageForceC4C5DiscNucleus_Y','AverageForcePllLeftBeam_Y','AverageForcePllRightBeam_Y',...
    'AverageForceSpinalCordNerves_Y','AverageForceLeftVertebralArtery_Y','AverageForceRightVertebralArtery_Y'};
anat_forc_z = {'AverageForceC4Vertebra_Z','AverageForceC5Vertebra_Z','AverageForceC4C5DiscAnnulus_Z',...
    'AverageForceC4C5DiscNucleus_Z','AverageForcePllLeftBeam_Z','AverageForcePllRightBeam_Z',...
    'AverageForceSpinalCordNerves_Z','AverageForceLeftVertebralArtery_Z','AverageForceRightVertebralArtery_Z'};

% Root sum square of all forces
for a = 1:length(anat)
    anat_forc = [x.(anat_forc_x{a}) x.(anat_forc_y{a}) x.(anat_forc_z{a})];
    for i = 1:length(anat_forc)
        Rssq.(anat{a})(i,:) = rssq(anat_forc(i,:));
    end
end

% Overall
for a = 1:length(anat)
    Rssq_mean.(anat{a}) = Rssq.(anat{a});
    Rssq_mean.(anat{a})(find(Rssq_mean.(anat{a})==0)) = []; %remove when forces = 0 as this affects the average
    if any(Rssq_mean.(anat{a}))>0
        Metrics(sj).(strcat((anat{a}),'Overall','Force_mean')) = mean(Rssq.(anat{a}));
    else
        Metrics(sj).(strcat((anat{a}),'Overall','Force_mean')) = 0;
    end
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Forc_tool = Rssq.(anat{a})(ToolUsed.(tools{t}));
        Forc_tool(find(Forc_tool==0)) = [];
        if any(Forc_tool)>0
            Metrics(sj).(strcat((anat{a}),(tools{t}),'Force_mean')) = mean(Forc_tool);
        else
            Metrics(sj).(strcat((anat{a}),(tools{t}),'Force_mean')) = 0;
        end
        clear Forc_tool
    end
end

%% Amount of time spent in contact with structure with each tool
% Overall
for a = 1:length(anat)
    if any(Contact.(anat{a}))>0
   %     Metrics(sj).(strcat((anat{a}),'Overall','Contact_time')) = ireg_time(Contact.(anat{a}),x.TimeSinceStart);
    else
   %     Metrics(sj).(strcat((anat{a}),'Overall','Contact_time')) = 0;
    end
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
  %          Metrics(sj).(strcat((anat{a}),(tools{t}),'Contact_time')) = ireg_time(inter_tool_anat,x.TimeSinceStart);
        else 
   %         Metrics(sj).(strcat((anat{a}),(tools{t}),'Contact_time')) = 0;
        end
    end
end

%% Velocity while in contact with structures for each tool
% check this with new data

% Velocity while contact with each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            [~,...
                ~,...
                ~,...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanOverall'))]...
                = ireg_vel_mean(inter_tool_anat,x.VirtualToolTipPosition_X,...
                x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z,...
                x.TimeSinceStart);
            %Metrics(sj).VelocityWhileContact_mean.(tools{t}).(anat{a}) = Rssq_vel(Contact.(anat{a}));
        else
            %Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanX')) = 0;
            %Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanY')) = 0;
            %Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanZ')) = 0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanOverall')) = 0;
        end
    end
end


%% Acceleration while in contact with structure for each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if length(inter_tool_anat)>2
            [Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanX')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanY')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanZ')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxX')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxY')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxZ')),...
                ~,...
                ~,...
                ~,...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numX')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numY')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numZ')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanX')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanY')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanZ')),...
                ~,...
                ~,...
                ~]...
                = ireg_acc_mean(inter_tool_anat,x.VirtualToolTipPosition_X,...
                x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z,...
                x.TimeSinceStart);
        else
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanX')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanY')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanZ')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxX')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxY')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxZ')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minX')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minY')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minZ')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numX')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numY')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numZ')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanX')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanY')) = 0;
                Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanZ')) = 0;
%                 Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanX')) = 0;
%                 Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanY')) = 0;
%                 Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanZ')) = 0;
        end
    end
end


%% Total Path Length for Tool while in contact with structure
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            [Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_X')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Y')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Z')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Overall')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_X')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Y')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Z')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Overall')),...
                Metrics(sj).(strcat((anat{a}),(tools{t}),'NumberStrokesWhileContact'))]=ireg_TTPL(inter_tool_anat,x.VirtualToolTipPosition_X,x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z);
        else
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_X'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Y'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Z'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Overall'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_X'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Y'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Z'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Overall'))=0;
            Metrics(sj).(strcat((anat{a}),(tools{t}),'NumberStrokesWhileContact'))=0;
        end
    end
end
if d==1
    tab1 = struct2table(Metrics(sj));
    Metrics_1M(sj,:) = tab1;
elseif d==2
    tab2 = struct2table(Metrics(sj));
    Metrics_2J(sj,:) = tab2;
elseif d==3
    tab3 = struct2table(Metrics(sj));
    Metrics_3S(sj,:) = tab3;
elseif d==4
    tab4 = struct2table(Metrics(sj));
    Metrics_4C(sj,:) = tab4;
end
clear Metrics
end
toc
end
end

%% Variance Thresholding
Metrics_matss = table2array(Metrics_4C);
for i=1:size(Metrics_matss,2)
    varb(i) = var(Metrics_matss(:,i));
end
avg_varb5 = mean(varb);


%% Preparation for Neural Nets
% Convert all tables to a 
Metrics_mat_1M = table2array(Metrics_1M);
Metrics_mat_2J = table2array(Metrics_2J);
Metrics_mat_3S = table2array(Metrics_3S);
Metrics_mat_4C = table2array(Metrics_4C);
All = [Metrics_mat_1M;Metrics_mat_2J;Metrics_mat_3S;Metrics_mat_4C];
All_normR=normalize(All);
All_norm=All_normR(:,all(~isnan(All_normR)));
% Create labels
num_1M=size(Metrics_1M,1);
num_2J=size(Metrics_2J,1);
num_3S=size(Metrics_3S,1);
num_4C=size(Metrics_4C,1);
total=num_1M+num_2J+num_3S+num_4C;

ANN_labs = repelem(1,num_1M);
ANN_labs = [ANN_labs repelem(0,total-num_1M)];
ANN_labs = [ANN_labs;repelem([0 1 0],[num_1M num_2J total-num_1M-num_2J])];
ANN_labs = [ANN_labs;repelem([0 1 0],[num_1M+num_2J num_3S num_4C])];
ANN_labs = [ANN_labs;repelem([0 1],[total-num_4C num_4C])];
ANN_labels=ANN_labs';


%% Feature Selection
labels = repelem([1 2 3 4],[num_1M num_2J num_3S num_4C])';
[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(All_norm,labels);
sig_mets = find(inmodel);
Best = All_norm(:,sig_mets);


%% Neural Net
inputs = Best';
targets = ANN_labels';itt=1;warning off
for layers=2:10
    for pp=1:50
hiddenLayerSize = layers;
net = patternnet(hiddenLayerSize);

net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 30/100;
% trnn=[1:24 41:93 109:166];
% ttt=[25:40 94:108 167:172];

    % Parameter Adjustments
    net.trainFcn='trainbr';
    net.trainParam.epochs=1000;
    net.trainParam.showWindow=false;
    
%     net.divideFcn = 'divideind';
% net.divideParam.trainInd = trnn;
% net.divideParam.testInd=ttt;


    % Train the Network
    [net,tr] = train(net,inputs,targets);

    % Test the Network
    outputs = net(inputs);
    errors = gsubtract(targets,outputs);
    performance = perform(net,targets,outputs);

    % View the Network
    %view(net)

    % Confusion Matrix
    %figure, plotconfusion(targets,outputs)

    % Confusion Matrix for test group
    tInd = tr.testInd;
    tstOutputs = net(inputs(:,tInd));
    tstTargets = targets(:,tInd);
    %figure, plotconfusion(tstTargets,tstOutputs);

    % Extract accuracy
    [c,cm]=confusion(tstTargets,tstOutputs);
    tst_acc=1-c;
    [ct,cm]=confusion(targets,outputs);
    tr_acc=1-ct;
        fprintf('Itt. %i, Layers %i, Epochs %i, Training Acc %2.2f, Testing Acc %2.2f \n',...
            pp,layers,...
            net.trainParam.epochs,...
            round(tr_acc,3,'significant'),round(tst_acc,3,'significant'));
         if tst_acc>0.8
             %figure, plotconfusion(tstTargets,tstOutputs);
             %break
         end
        itt=itt+1;
    end
end
