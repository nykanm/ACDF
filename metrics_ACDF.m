%% function [sub,sig_table,sig_vars,nov,exp,novmat,expmat] = metrics_ACDF(data_nov,data_exp,sig_level)
% 1. Run csv_concat() to concatenate all tool files into one for
% '1M','2J', '3S', '4C'
% 2. Run [Instrument] = csv2struct; to convert into structure.

%% Data Conversions (only do this once)
% Take all the data from Ossim and put them into a folder called Sample
% Data. Create a folder for each group such as 1M, 2J, 3S, 4C
% Change this to the directory where your raw data is stored:
DataDir = '/Volumes/Seagate/OSSIM STUDY DATA/Corrected_2/';
csv_concat('1M',DataDir);
csv_concat('2J',DataDir);
csv_concat('3S',DataDir);
csv_concat('4F',DataDir);
csv_concat('5C',DataDir);
% This creates the .mat files for all groups in the same folder as the raw
% data.
csv2struct2('1M_Compact',DataDir);
csv2struct2('2J_Compact',DataDir);
csv2struct2('3S_Compact',DataDir);
csv2struct2('4F_Compact',DataDir);
csv2struct2('5C_Compact',DataDir);

%% Load your data
DataDir = '/Users/nykan/Documents/McGill Grad/Matlab/ACDF/Sample Data/';
DataDir = '/Volumes/Seagate/OSSIM STUDY DATA/';

load([DataDir 'Data_1M_Compact.mat']);
load([DataDir 'Data_2J_Compact.mat']);
load([DataDir 'Data_3S_Compact.mat']);
load([DataDir 'Data_4C_Compact.mat']);


%%
% Cutting Disc: step=1;
% Removing Disc: step=2;
% Removing Osteophytes: step=3;
step=2;
for d=2:5
    tic
    if d==1
        data=Data_1M_Compact;
        %data=table2struct(AllData_1M);
        disp(['Med Students Total: ',num2str(length(data))]);
    elseif d==2
        data=Data_2J_Compact;
        %data=table2struct(AllData_2J);
        disp(['Junior Residents Total: ',num2str(length(data))]);
    elseif d==3
        data=Data_3S_Compact;
        %data=table2struct(AllData_3S);
        disp(['Senior Residents Total: ',num2str(length(data))]);
    elseif d==4
        data=Data_4F_Compact;
        %data=table2struct(AllData_4F);
        disp(['Fellows Total: ',num2str(length(data))]);
     elseif d==5
        data=Data_5C_Compact;
        %data=table2struct(AllData_5C);
        disp(['Consultants Total: ',num2str(length(data))]);
    end
for sj=1:length(data)
    disp(['Subject: ',num2str(sj),' of ',num2str(length(data))]);
    y=data(sj);
if length(y.TimeSinceStart)>0
%% Find time duplicates
timeDupInd = find([0.1;diff(y.TimeSinceStart)]); 
fieldNam = fieldnames(y);
for i = 1:numel(fieldNam)
    z.(fieldNam{i}) = y.(fieldNam{i})(timeDupInd);
end

%% Restrict the Step/tool used
% when each tool is used
ToolUsed.Scalpel = find(contains(z.ToolUsed,'Scalpel'));
ToolUsed.BoneCurette = find(contains(z.ToolUsed,'Bone Curette'));
ToolUsed.Rongeur2mm = find(contains(z.ToolUsed,'Pituitary Rongeur 2mm'));
ToolUsed.DiscRongeur = find(contains(z.ToolUsed,'Disc Rongeur'));
ToolUsed.Burr = find(contains(z.ToolUsed,'Burr'));
%ToolUsed.NerveHook = find(contains(x.ToolUsed,'Nerve Hook'));
%ToolUsed.Kerrison1mm = find(contains(x.ToolUsed,'Kerrison 1mm'));
SurgStep.Cutting = ToolUsed.Scalpel;
SurgStep.RemoveDisc = cat(1,ToolUsed.BoneCurette,ToolUsed.DiscRongeur,ToolUsed.Rongeur2mm);
SurgStep.Osteo = ToolUsed.Burr;
toolfield = fieldnames(SurgStep);
for i = 1:numel(fieldNam)
    x.(fieldNam{i}) = z.(fieldNam{i})(SurgStep.(toolfield{step}));
end
clear ToolUsed
ToolUsed.Scalpel = find(contains(x.ToolUsed,'Scalpel'));
ToolUsed.BoneCurette = find(contains(x.ToolUsed,'Bone Curette'));
ToolUsed.Rongeur2mm = find(contains(x.ToolUsed,'Pituitary Rongeur 2mm'));
ToolUsed.DiscRongeur = find(contains(x.ToolUsed,'Disc Rongeur'));
ToolUsed.Burr = find(contains(x.ToolUsed,'Burr'));
%% Define anat and tools
if step==1
    tools = {'Scalpel'};
    %anat = {'C4C5DiscAnnulus','C4C5DiscNucleus','PllLeft','PllRight','SC','LVA','RVA'};
elseif step==2
    tools = {'BoneCurette','Rongeur2mm','DiscRongeur'};
   % anat = {'C4C5DiscAnnulus','C4C5DiscNucleus','PllLeft','PllRight','SC','LVA','RVA'};
elseif step==3
    tools = {'Burr'};
   % anat = {'C4','C5','PllLeft','PllRight','SC','LVA','RVA'};
end
%tools = {'Scalpel','BoneCurette','Rongeur2mm','DiscRongeur','Burr'};
%tools = {'BoneCurette'};
anat = {'C4','C5','C4C5DiscAnnulus','C4C5DiscNucleus','PllLeft','PllRight','SC','LVA','RVA'};

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

%% Number of touches (instenses) (EFFICIENCY) (CHECKED)

% Overall
for a = 1:length(anat)
    Metrics(sj).(strcat((anat{a}),'Overall','Contact_num')) = numel(find(diff(Contact_con1.(anat{a}))>0));
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Contact_1.(anat{a}).(tools{t}) = Contact_con1.(anat{a})(ToolUsed.(tools{t}));  
        Metrics(sj).(strcat((anat{a}),(tools{t}),'Contact_num')) = numel(find(diff(Contact_1.(anat{a}).(tools{t}))>0));
    end
end

%% Number of cuts (instenses) (EFFICIENCY) (CHECKED)

% Overall
for a = 1:length(anat)
    Metrics(sj).(strcat((anat{a}),'Overall','Cutting_num')) = numel(find(diff(Cutting_con1.(anat{a}))>0));
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Cutting_1.(anat{a}).(tools{t}) = Cutting_con1.(anat{a})(ToolUsed.(tools{t}));  
        Metrics(sj).(strcat((anat{a}),(tools{t}),'Cutting_num')) = numel(find(diff(Cutting_1.(anat{a}).(tools{t}))>0));
    end
end

%% Number of voxels cut (QUALITY OR SAFETY) (CHECKED)

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

%% Average force (SAFETY) (CHECKED)
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
        Metrics(sj).(strcat((anat{a}),'Overall','Force_mean')) = mean(Rssq_mean.(anat{a}));
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

%% Max force (SAFETY) (CHECKED)

% Overall
for a = 1:length(anat)
    Rssq_mean.(anat{a}) = Rssq.(anat{a});
    if any(Rssq_mean.(anat{a}))>0
        Metrics(sj).(strcat((anat{a}),'Overall','Force_max')) = max(Rssq.(anat{a}));
    else
        Metrics(sj).(strcat((anat{a}),'Overall','Force_max')) = 0;
    end
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Forc_tool = Rssq.(anat{a})(ToolUsed.(tools{t}));
        Forc_tool(find(Forc_tool==0)) = [];
        if any(Forc_tool)>0
            Metrics(sj).(strcat((anat{a}),(tools{t}),'Force_max')) = max(Forc_tool);
        else
            Metrics(sj).(strcat((anat{a}),(tools{t}),'Force_max')) = 0;
        end
        clear Forc_tool
    end
end


%% Amount of time spent in contact with structure with each tool (EFFICIENCY) (CHECKED)
% Overall
for a = 1:length(anat)
    if any(Contact.(anat{a}))>0
        Metrics(sj).(strcat((anat{a}),'Overall','Contact_time')) = ireg_time(Contact.(anat{a}),x.TimeSinceStart);
    else
        Metrics(sj).(strcat((anat{a}),'Overall','Contact_time')) = 0;
    end
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            Metrics(sj).(strcat((anat{a}),(tools{t}),'Contact_time')) = ireg_time(inter_tool_anat,x.TimeSinceStart);
        else 
            Metrics(sj).(strcat((anat{a}),(tools{t}),'Contact_time')) = 0;
        end
    end
end


%% Velocity while in contact with structures for each tool (MOTIONS)
% check this with new data
if step==2
    anat_vel = {'C4C5DiscAnnulus','C4C5DiscNucleus'};
end
% Velocity while contact with each tool
for a = 1:length(anat_vel)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat_vel{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            [~,...
                ~,...
                ~,...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'VelocityWhileContact_mean'))]...
                = ireg_vel_mean(inter_tool_anat,x.VirtualToolTipPosition_X,...
                x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z,...
                x.TimeSinceStart);
            %Metrics(sj).VelocityWhileContact_mean.(tools{t}).(anat{a}) = Rssq_vel(Contact.(anat{a}));
        else
            %Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanX')) = 0;
            %Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanY')) = 0;
            %Metrics(sj).(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanZ')) = 0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'VelocityWhileContact_mean')) = 0;
        end
    end
end


%% Acceleration while in contact with structure for each tool (MOTIONS)
for a = 1:length(anat_vel)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat_vel{a}),ToolUsed.(tools{t}));
        if length(inter_tool_anat)>2
            [Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_meanX')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_meanY')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_meanZ')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_maxX')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_maxY')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_maxZ')),...
                ~,...
                ~,...
                ~,...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_numX')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_numY')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_numZ')),...
                ~,...
                ~,...
                ~,...
                ~,...
                ~,...
                ~]...
                = ireg_acc_mean(inter_tool_anat,x.VirtualToolTipPosition_X,...
                x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z,...
                x.TimeSinceStart);
        else
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_meanX')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_meanY')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_meanZ')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_maxX')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_maxY')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_maxZ')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minX')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minY')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minZ')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_numX')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_numY')) = 0;
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'AccelerationWhileContact_numZ')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanX')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanY')) = 0;
                %Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanZ')) = 0;
%                 Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanX')) = 0;
%                 Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanY')) = 0;
%                 Metrics(sj).(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanZ')) = 0;
        end
    end
end


%% Total Path Length for Tool while in contact with structure (EFFICIENCY)
for a = 1:length(anat_vel)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat_vel{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            [Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_X')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_Y')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_Z')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_XYZ')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_X')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_Y')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_Z')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_XYZ')),...
                Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'NumberStrokesWhileContact'))]=ireg_TTPL(inter_tool_anat,x.VirtualToolTipPosition_X,x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z);
        else
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_X'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_Y'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_Z'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLWhileContact_XYZ'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_X'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_Y'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_Z'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'TTPLPerStrokeWhileContact_XYZ'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'NumberStrokesWhileContact'))=0;
        end
    end
end

%% Choice of instrument (COGNITIVE)
for t=1:length(tools)
    if length(ToolUsed.(tools{t}))>2
        Metrics(sj).(strcat('ToolChoice_',(tools{t})))=1;
    else 
        Metrics(sj).(strcat('ToolChoice_',(tools{t})))=0;
    end
end
    
%% Angles (MOTION)
quat = [x.ToolQuaternion_QX x.ToolQuaternion_QY x.ToolQuaternion_QZ x.ToolQuaternion_QW];
[roll, pitch, yaw]=quat2angle(quat,'XYZ');
for t=1:length(tools)
    for a=1:length(anat_vel)
%     if t==5 
%         an=1:2;
%     else
%         an=3:4;
%     end
%    for a=an(1):an(end)
        inter_tool_anat = intersect(Contact.(anat_vel{a}),ToolUsed.(tools{t}));
        if length(ToolUsed.(tools{t}))>2
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'RollWhileContact_mean')) = mean(roll(inter_tool_anat));
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'PitchWhileContact_mean')) = mean(pitch(inter_tool_anat));
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'YawWhileContact_mean')) = mean(yaw(inter_tool_anat));
        else
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'RollWhileContact_mean')) = 0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'PitchWhileContact_mean')) = 0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'YawWhileContact_mean')) = 0;
        end
%    end
    end
end

%% Angular Velocity (MOTION)

for t=1:length(tools)
    for a=1:length(anat_vel)
%     if t==5 
%         an=1:2;
%     else
%         an=3:4;
%     end
%     for a=an(1):an(end)
        inter_tool_anat = intersect(Contact.(anat_vel{a}),ToolUsed.(tools{t}));
        if length(ToolUsed.(tools{t}))>2
        [Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'RollAngVelWhileContact_mean')),...
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'PitchAngVelWhileContact_mean')),...
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'YawAngVelWhileContact_mean')),~] = ireg_vel_mean(inter_tool_anat,roll,pitch,yaw,x.TimeSinceStart);
        else
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'RollAngVelWhileContact_mean'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'PitchAngVelWhileContact_mean'))=0;
            Metrics(sj).(strcat((anat_vel{a}),(tools{t}),'YawAngVelWhileContact_mean'))=0;
        end
%     end
    end
end

%%
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
    Metrics_4F(sj,:) = tab4;
elseif d==5
    tab5 = struct2table(Metrics(sj));
    Metrics_5C(sj,:) = tab5;
end
clear Metrics
end
end

%% Stats and Data Organization
% if d==1
%     nov = struct2table(Metrics(sj));
%     novmat = table2array(nov);
%     % repalce nan with the mean of the column
%     n = nanmean(novmat);
%     nn = isnan(novmat);
%     ii = sum(nn) < 3;
%     z = novmat(:,ii);
%     z(nn(:,ii)) = nonzeros(bsxfun(@times,nn(:,ii),n(ii)));
%     novmat(:,ii) = z;
%     novtbl = array2table(novmat);
%     nov(:,:) = novtbl;
% elseif d==2
%     exp = struct2table(Metrics(sj));
%     expmat = table2array(exp);
%     % repalce nan with the mean of the column
%     e = nanmean(expmat);
%     ee = isnan(expmat);
%     iie = sum(ee) < 3;
%     ze = expmat(:,iie);
%     ze(ee(:,iie)) = nonzeros(bsxfun(@times,ee(:,iie),e(iie)));
%     expmat(:,iie) = ze;
%     size(expmat)
%     exptbl = array2table(expmat);
%     exp(:,:) = exptbl;
% 
% [sub,sig_table,sig_vars] = metrics_lamin_stats(novmat,nov,expmat,exp,0.05);
% 
% toc
% end
% clear Metrics;
end

%% Remove irrelevant metrics
if step==1
    % remove metrics that say OVERALL, not necessary since only one
    % instrument
    toRem11 = find(contains(Metrics_1M.Properties.VariableNames,'Overall'));
    Metrics_1M(:,toRem11)=[];
    Metrics_2J(:,toRem11)=[];
    Metrics_3S(:,toRem11)=[];
    Metrics_4F(:,toRem11)=[];
    Metrics_5C(:,toRem11)=[];
    % remove metrics with scalpel on C4 or C5, since not relevant
    toRem12 = find(contains(Metrics_1M.Properties.VariableNames,'C4Scalpel'));
    Metrics_1M(:,toRem12)=[];
    Metrics_2J(:,toRem12)=[];
    Metrics_3S(:,toRem12)=[];
    Metrics_4F(:,toRem12)=[];
    Metrics_5C(:,toRem12)=[];    
    toRem13 = find(contains(Metrics_1M.Properties.VariableNames,'C5Scalpel'));
    Metrics_1M(:,toRem13)=[];
    Metrics_2J(:,toRem13)=[];
    Metrics_3S(:,toRem13)=[];
    Metrics_4F(:,toRem13)=[];
    Metrics_5C(:,toRem13)=[];     
end

%% Preparation for Neural Nets
% Convert all tables to a 
%Metrics_mat_1M = table2array(Metrics_1M);
Metrics_mat_2J = table2array(Metrics_2J);
Metrics_mat_3S = table2array(Metrics_3S);
Metrics_mat_4F = table2array(Metrics_4F);
Metrics_mat_5C = table2array(Metrics_5C);

All = [Metrics_mat_2J;Metrics_mat_3S;Metrics_mat_4F;Metrics_mat_5C];
% remove all 0 columns
emptyIdx=find(all(All==0));
All(:,emptyIdx)=[];
names=Metrics_2J;
names(:,emptyIdx)=[];
All_normR=normalize(All);
All_norm=All_normR(:,all(~isnan(All_normR)));
names = names(:,all(~isnan(All_normR)));
choices=find(contains(names.Properties.VariableNames,'ToolChoice'));
% Create labels
num_1M=0;
num_2J=size(Metrics_2J,1);
num_3S=size(Metrics_3S,1);
num_4F=size(Metrics_4F,1);
num_5C=size(Metrics_5C,1);

total=num_1M+num_2J+num_3S+num_4F+num_5C;

ANN_labs = repelem(1,num_2J);
ANN_labs = [ANN_labs repelem(0,total-num_2J)];
ANN_labs = [ANN_labs;repelem([0 1 0],[num_2J num_3S num_4F+num_5C])];
ANN_labs = [ANN_labs;repelem([0 0 1],[num_2J num_3S num_4F+num_5C])];
%ANN_labs = [ANN_labs;repelem([0 1 0],[total-num_4F-num_5C num_4F num_5C])];
%ANN_labs = [ANN_labs;repelem([0 1],[total-num_5C-num_4F num_5C+num_4F])];
ANN_labels=ANN_labs';


targets = ANN_labels';itt=1;warning off;

%% Feature selection
% Define training and testing set 
%trnn=[1 3:5 7 9:11 13 14 16:19 21 23:25 27];
% ttt=[2 6 8 12 15 20 22 26];
% no med students; step 1
trnn=[2 3:5 7 8 10:13 15 17:19 21];
ttt=[1 6 9 14 16 20]; 

% no med students; step 2
trnn=[2 3:5 7 8 10:13 15 17:19 21];
ttt=[1 6 9 14 16 20]; 

%num2_1 = 4;
%num2_1 = 0;
num2_2 = 5;
num2_3 = 4;
num2_4 = 6;

All_norm_train = All_norm(trnn,:);
labels2 = repelem([2 3 4],[num2_2 num2_3 num2_4])';
[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(All_norm_train,labels2,'penter',0.05);
sig_mets = find(inmodel);
sig_mets2=[sig_mets choices];
Best = All_norm(:,sig_mets2);
inputs = Best';

%name of metrics
Sig_names = names(1,sig_mets2);

%% Neural Net
clear train_accuracy test_accuracy outt
for layers=16
    for pp=1:10
hiddenLayerSize = layers;
net = patternnet(hiddenLayerSize);

% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 30/100;



%%
    % Parameter Adjustments
    net.trainFcn='trainbr';
    net.trainParam.epochs=1000;
    net.trainParam.showWindow=false;
    net.layers{1}.transferFcn='tansig';
    net.layers{2}.transferFcn='tansig';
    net.divideFcn = 'divideind';
 net.divideParam.trainInd = trnn;
 net.divideParam.testInd=ttt;
 %net.trainParam.lr=0.01;
 %net.trainParam.mc=0.75;
 %mu = pp*0.01;
 net.trainParam.Mu=0.01;
 net.trainParam.mu_dec=0.95;


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
         if tst_acc>0.60 & tr_acc>0.80
           % figure, plotconfusion(targets,outputs);
           % figure, plotconfusion(tstTargets,tstOutputs);
             %fprintf('WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOW \n');
            % break;break;break;
         end
        itt=itt+1;
    end
    train_accuracy(layers)=round(tr_acc,3,'significant');
    test_accuracy(layers)=round(tst_acc,3,'significant');
end
%outt=[[1:15]',train_accuracy',test_accuracy'];

%% Calculate metric importance

Weights1 = abs(net.IW{1});Weights2 = abs(net.LW{2});

for num_ins=1:16
    for k=1:16
        impJ(k)=Weights1(k,num_ins)*Weights2(1,k);
        impS(k)=Weights1(k,num_ins)*Weights2(2,k);
        impP(k)=Weights1(k,num_ins)*Weights2(3,k);
    end
    ImpJ(num_ins,1)=sum(impJ);
    ImpS(num_ins,1)=sum(impS);
    ImpP(num_ins,1)=sum(impP);
end

tot_J = sum(ImpJ);
tot_S = sum(ImpS);
tot_P = sum(ImpP);


for num_ins=1:16
    rel_impJ(num_ins,1) = ImpJ(num_ins,1)/tot_J * 100;
    rel_impS(num_ins,1) = ImpS(num_ins,1)/tot_S * 100;
    rel_impP(num_ins,1) = ImpP(num_ins,1)/tot_P * 100;

end

Feat_importance_J = table(ImpJ,rel_impJ,Sig_names.Properties.VariableNames');
Feat_importance_J = sortrows(Feat_importance_J,2,'descend');
Feat_importance_S = table(ImpS,rel_impS,Sig_names.Properties.VariableNames');
Feat_importance_S = sortrows(Feat_importance_S,2,'descend');
Feat_importance_P = table(ImpP,rel_impP,Sig_names.Properties.VariableNames');
Feat_importance_P = sortrows(Feat_importance_P,2,'descend');

%% Bar Graph
Metrics_mat_2J = table2array(Metrics_2J);
Metrics_mat_3S = table2array(Metrics_3S);
Metrics_mat_4F = table2array(Metrics_4F);
Metrics_mat_5C = table2array(Metrics_5C);

Metrics_mat_2J(:,emptyIdx)=[];
Metrics_mat_2J = Metrics_mat_2J(:,all(~isnan(All_normR)));
Metrics_mat_3S(:,emptyIdx)=[];
Metrics_mat_3S = Metrics_mat_3S(:,all(~isnan(All_normR)));
Metrics_mat_4F(:,emptyIdx)=[];
Metrics_mat_4F = Metrics_mat_4F(:,all(~isnan(All_normR)));
Metrics_mat_5C(:,emptyIdx)=[];
Metrics_mat_5C = Metrics_mat_5C(:,all(~isnan(All_normR)));


plot_2J = Metrics_mat_2J(:,sig_mets2);
plot_3S = Metrics_mat_3S(:,sig_mets2);
merge = [Metrics_mat_4F;Metrics_mat_5C];
plot_4C = merge(:,sig_mets2);
mean_plot_2J=mean(plot_2J);
mean_plot_3S=mean(plot_3S);
mean_plot_4C=mean(plot_4C);

bar_mean_plot = [mean_plot_2J(1) mean_plot_3S(1) mean_plot_4C(1)];
for i=2:16
    bar_mean_plot = [bar_mean_plot;mean_plot_2J(i) mean_plot_3S(i) mean_plot_4C(i)];
end
for i=1:16
    figure(1);
    subplot(2,8,i);
    bar(bar_mean_plot(i,[1:3]));
end

