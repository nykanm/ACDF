function [ttplX,ttplY,ttplZ,ttpl_overall,...
    ttplX_stroke,ttplY_stroke,ttplZ_stroke,ttpl_overall_stroke,...
    numRm] = ireg_TTPL(condition,posX,posY,posZ)
%% [ttplX, ttplY, ttplZ,numRm] = ireg_vel(condition,posX,posY,posZ,time)
% condition must be a vector from a find function (e.g. whether you are touching the
 %ttpl along each axis, numRm is the number of blocks that you have
  one = zeros(max(condition)+1,1);
  one(condition)=1; 
  % -1 are starts, 1 are ends
  bin = diff(one); 
    if any(condition)==1 % checks that condition is not empty
        % this adds a start if the action begins at time 0
        startpos=find(bin==1);endpos=find(bin==-1);
        if any(startpos)==0
            if condition(1)==1
                bin(1)=1;
                %bin(1)=-1; % for datadivision.m
            end
        else
            if endpos(1)==1
                bin(1)=0;
            elseif startpos(1) > endpos(1)
                if condition(1)==1
                    bin(1)=1;
                    %bin(1)=-1; % for datadivision.m
                end
            end
        end

        starts = find(bin==1);
        ends = find(bin==-1);
        
        numRm = length(starts);

        for idx=1:length(starts)
            pos_x = posX(starts(idx):ends(idx));
            pos_y = posY(starts(idx):ends(idx));
            pos_z = posZ(starts(idx):ends(idx));
            
            pos_xyz = [pos_x pos_y pos_z];
            
            % ttpl along each axis   
            ttpl_x{idx} = abs(diff(pos_x));
            ttpl_y{idx} = abs(diff(pos_y));
            ttpl_z{idx} = abs(diff(pos_z));
            
            % ttpl overall
            for i=1:size(pos_xyz,1)-1
                ttpl_xyz_i(i) = norm(pos_xyz(i,:) - pos_xyz(i+1,:));
            end
            ttpl_xyz{idx} = ttpl_xyz_i';
            clear ttpl_xyz_i;
        end

        % concatonate all the vectors
        ttpl_x_fin = cat(1,ttpl_x{:});
        ttpl_y_fin = cat(1,ttpl_y{:});
        ttpl_z_fin = cat(1,ttpl_z{:});
        ttpl_xyz_fin = cat(1,ttpl_xyz{:});

        ttplX = sum(ttpl_x_fin);
        ttplY = sum(ttpl_y_fin);
        ttplZ = sum(ttpl_z_fin);
        ttpl_overall = sum(ttpl_xyz_fin);
        
        ttplX_stroke = ttplX/numRm;
        ttplY_stroke = ttplY/numRm;
        ttplZ_stroke = ttplZ/numRm;
        ttpl_overall_stroke = ttpl_overall/numRm;
        
        
    else
        ttplX = NaN;
        ttplY = NaN;
        ttplZ = NaN;
        ttplX_stroke = NaN;
        ttplY_stroke = NaN;
        ttplZ_stroke = NaN;
        ttpl_overall = NaN;
        ttpl_overall_stroke = NaN;
        
        numRm = NaN;
        
    end

    
    
    
    
    