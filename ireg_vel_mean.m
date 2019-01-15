function [velX_mean,velY_mean,velZ_mean,velOverall_mean] = ireg_vel_mean(condition,posX,posY,posZ,time)
%% [velX_mean,velY_mean,velZ_mean,velOverall_mean] = ireg_vel(condition,posX,posY,posZ,time)
% condition must be a vector from a find function (e.g. whether you are touching the
 % structure or not, or whether you are using a tool or not
  one = zeros(max(condition)+1,1);
  one(condition)=1; 
  % -1 are starts, 1 are ends
  bin = diff(one); 
    if any(condition)==1 % checks that condition is not empty
        % this adds a start if the action begins at time 0
        if condition(1)==1
            bin(1)=1;
        end

        starts = find(bin==1);
        ends = find(bin==-1);

        for idx=1:length(starts)
            pos_x = posX(starts(idx):ends(idx));
            pos_y = posY(starts(idx):ends(idx));
            pos_z = posZ(starts(idx):ends(idx));
            times = time(starts(idx):ends(idx));

            vel_x{idx}=gradient(pos_x,times);
            vel_y{idx}=gradient(pos_y,times);
            vel_z{idx}=gradient(pos_z,times);
        end

        % concatonate all the vectors
        vel_x_fin = cat(1,vel_x{:});
        vel_y_fin = cat(1,vel_y{:});
        vel_z_fin = cat(1,vel_z{:});
        vel_xyz = [vel_x_fin vel_y_fin vel_z_fin];

        if any(vel_xyz)==1
            for i=1:size(vel_xyz,1)
                vel_total(i,:) = rssq(vel_xyz(i,:));
            end

            velX_mean = mean(abs(vel_x_fin));
            velY_mean = mean(abs(vel_y_fin));
            velZ_mean = mean(abs(vel_z_fin));

            velOverall_mean = mean(vel_total);
        else
            velX_mean = NaN;
            velY_mean = NaN;
            velZ_mean = NaN;

            velOverall_mean = NaN;
        end
    else
        velX_mean = NaN;
        velY_mean = NaN;
        velZ_mean = NaN;

        velOverall_mean = NaN;
    end

    
    
    
    
    