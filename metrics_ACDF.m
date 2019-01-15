% 1. Run csv_concat to concatenate all tool files into one.
% 2. Run [Instrument] = csv2struct; to convert into structure.

for d=1:1
    tic
%     if d==1
%         data=data_nov;
%         disp(['Novices Total: ',num2str(length(data))]);
%     elseif d==2
%         data=data_exp;
%         disp(['Experts Total: ',num2str(length(data))]);
%     end
for sj=1:length(data)
    disp(['Subject: ',num2str(sj),' of ',num2str(length(data))]);
    y=data(sj);

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
tools = {'Scalpel','BoneCurette','Rongeur2mm','DiscRongeur','Burr','NerveHook','Kerrison1mm'};
anat = {'C4','C5','C4C5DiscAnnulus','C4C5DiscNucleus','PllLeft','PllRight','SC','LVA','RVA'};

%% Number of touches (instenses)

% Overall
for a = 1:length(anat)
    Metrics.(strcat((anat{a}),'Overall','Contact_num')) = numel(find(diff(Contact_con1.(anat{a}))>0));
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Contact_1.(anat{a}).(tools{t}) = Contact_con1.(anat{a})(ToolUsed.(tools{t}));  
        Metrics.(strcat((anat{a}),(tools{t}),'Contact_num')) = numel(find(diff(Contact_1.(anat{a}).(tools{t}))));
    end
end

%% Number of cuts (instenses)

% Overall
for a = 1:length(anat)
    Metrics.(strcat((anat{a}),'Overall','Cutting_num')) = numel(find(diff(Cutting_con1.(anat{a}))>0));
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Cutting_1.(anat{a}).(tools{t}) = Cutting_con1.(anat{a})(ToolUsed.(tools{t}));  
        Metrics.(strcat((anat{a}),(tools{t}),'Cutting_num')) = numel(find(diff(Cutting_1.(anat{a}).(tools{t}))));
    end
end

%% Number of voxels cut

% Overall
for a = 1:length(anat_cutvox)
    Metrics.(strcat((anat{a}),'Overall','Cut_voxels')) = sum(x.(anat_cutvox{a}));
end

% With each tool
for a = 1:length(anat_cutvox)
    for t = 1:length(tools)
        Metrics.(strcat((anat{a}),(tools{t}),'Cut_voxels')) = sum(x.(anat_cutvox{a})(ToolUsed.(tools{t})));  
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
        Metrics.(strcat((anat{a}),'Overall','Force_mean')) = mean(Rssq.(anat{a}));
    else
        Metrics.(strcat((anat{a}),'Overall','Force_mean')) = 0;
    end
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        Forc_tool = Rssq.(anat{a})(ToolUsed.(tools{t}));
        Forc_tool(find(Forc_tool==0)) = [];
        if any(Forc_tool)>0
            Metrics.(strcat((anat{a}),(tools{t}),'Force_mean')) = mean(Forc_tool);
        else
            Metrics.(strcat((anat{a}),(tools{t}),'Force_mean')) = 0;
        end
        clear Forc_tool
    end
end

%% Amount of time spent in contact with structure with each tool
% Overall
for a = 1:length(anat)
    if any(Contact.(anat{a}))>0
        Metrics.(strcat((anat{a}),'Overall','Contact_time')) = ireg_time(Contact.(anat{a}),x.TimeSinceStart);
    else
        Metrics.(strcat((anat{a}),'Overall','Contact_time')) = 0;
    end
end

% With each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            Metrics.(strcat((anat{a}),(tools{t}),'Contact_time')) = ireg_time(inter_tool_anat,x.TimeSinceStart);
        else 
            Metrics.(strcat((anat{a}),(tools{t}),'Contact_time')) = 0;
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
            [Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanX')),...
                Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanY')),...
                Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanZ')),...
                Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanOverall'))]...
                = ireg_vel_mean(inter_tool_anat,x.VirtualToolTipPosition_X,...
                x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z,...
                x.TimeSinceStart);
            %Metrics.VelocityWhileContact_mean.(tools{t}).(anat{a}) = Rssq_vel(Contact.(anat{a}));
        else
            Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanX')) = 0;
            Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanY')) = 0;
            Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanZ')) = 0;
            Metrics.(strcat((anat{a}),(tools{t}),'VelocityWhileContact_meanOverall')) = 0;
        end
    end
end


%% Acceleration while in contact with structure for each tool
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            [Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanX')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanY')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanZ')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxX')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxY')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxZ')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minX')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minY')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minZ')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numX')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numY')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numZ')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanX')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanY')),...
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanZ')),...
                ~,...
                ~,...
                ~]...
                = ireg_acc_mean(inter_tool_anat,x.VirtualToolTipPosition_X,...
                x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z,...
                x.TimeSinceStart);
        else
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanX')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanY')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_meanZ')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxX')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxY')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_maxZ')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minX')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minY')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_minZ')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numX')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numY')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_numZ')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanX')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanY')) = 0;
                Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_p2pmeanZ')) = 0;
%                 Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanX')) = 0;
%                 Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanY')) = 0;
%                 Metrics.(strcat((anat{a}),(tools{t}),'AccelerationWhileContact_regulmeanZ')) = 0;
        end
    end
end


%% Total Path Length for Tool while in contact with structure
for a = 1:length(anat)
    for t = 1:length(tools)
        inter_tool_anat = intersect(Contact.(anat{a}),ToolUsed.(tools{t}));
        if any(inter_tool_anat)>0
            [Metrics.(strcat((anat{a}),(tools{t}),'TTPLWhileContact_X')),...
                Metrics.(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Y')),...
                Metrics.(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Z')),...
                Metrics.(strcat((anat{a}),(tools{t}),'TTPLWhileContact_Overall')),...
                Metrics.(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_X')),...
                Metrics.(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Y')),...
                Metrics.(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Z')),...
                Metrics.(strcat((anat{a}),(tools{t}),'TTPLPerStrokeWhileContact_Overall')),...
                Metrics.(strcat((anat{a}),(tools{t}),'NumberStrokesWhileContact'))]=ireg_TTPL(inter_tool_anat,x.VirtualToolTipPosition_X,x.VirtualToolTipPosition_Y,x.VirtualToolTipPosition_Z);
        end
    end
end

end
end