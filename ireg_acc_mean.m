function [accX_mean,accY_mean,accZ_mean,...
    accX_max,accY_max,accZ_max,...
    accX_min,accY_min,accZ_min,...
    acc_num_x,acc_num_y,acc_num_z,...
    acc_x_p2p_mean,acc_y_p2p_mean,acc_z_p2p_mean,...
    acc_regul_x_mean,acc_regul_y_mean,acc_regul_z_mean] = ireg_acc_mean(condition,posX,posY,posZ,time)
%% [accX_mean,accY_mean,accZ_mean,accX_max,accY_max,accZ_max,accX_min,accY_min,accZ_min,acc_num_x,acc_num_y,acc_num_z,acc_x_p2p_mean,acc_y_p2p_mean,acc_z_p2p_mean,acc_regul_x_mean,acc_regul_y_mean,acc_regul_z_mean] = ireg_vel(condition,posX,posY,posZ,time)
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
            
            speed_x{idx}=abs(vel_x{idx});
            speed_y{idx}=abs(vel_y{idx});
            speed_z{idx}=abs(vel_z{idx});
            
            acc_x_abs{idx}=gradient(speed_x{idx},times);
            acc_y_abs{idx}=gradient(speed_y{idx},times);
            acc_z_abs{idx}=gradient(speed_z{idx},times);
            
            if length(acc_x_abs{idx}) > 2 
                [~,loc_x]=findpeaks(acc_x_abs{idx});
                [~,loc_y]=findpeaks(acc_y_abs{idx});
                [~,loc_z]=findpeaks(acc_z_abs{idx});

                acc_regul_x_abs{idx}=diff(times(loc_x));
                acc_regul_y_abs{idx}=diff(times(loc_y));
                acc_regul_z_abs{idx}=diff(times(loc_z));
            else
                acc_regul_x_abs{idx}=NaN;
                acc_regul_y_abs{idx}=NaN;
                acc_regul_z_abs{idx}=NaN;
            end
                

            acc_x{idx}=gradient(vel_x{idx},times);
            acc_y{idx}=gradient(vel_y{idx},times);
            acc_z{idx}=gradient(vel_z{idx},times);

        end

        % concatonate all the vectors
        acc_x_abs_fin = cat(1,acc_x_abs{:});
        acc_y_abs_fin = cat(1,acc_y_abs{:});
        acc_z_abs_fin = cat(1,acc_z_abs{:});
        acc_xyz = [acc_x_abs_fin acc_y_abs_fin acc_z_abs_fin];

        acc_x_fin = cat(1,acc_x{:});
        acc_y_fin = cat(1,acc_y{:});
        acc_z_fin = cat(1,acc_z{:});
        
        acc_regul_x_fin = cat(1,acc_regul_x_abs{:});
        acc_regul_y_fin = cat(1,acc_regul_y_abs{:});
        acc_regul_z_fin = cat(1,acc_regul_z_abs{:});
 
        
        if any(acc_xyz)==1
            % average acceleration
            accX_mean = mean(acc_x_abs_fin);
            accY_mean = mean(acc_y_abs_fin);
            accZ_mean = mean(acc_z_abs_fin);
            % max acceleration
            accX_max = max(acc_x_abs_fin);
            accY_max = max(acc_y_abs_fin);
            accZ_max = max(acc_z_abs_fin);
            % max deceleration
            accX_min = min(acc_x_abs_fin);
            accY_min = min(acc_y_abs_fin);
            accZ_min = min(acc_z_abs_fin);
            
            acc_x_positives = acc_x_fin > 0;
            acc_y_positives = acc_y_fin > 0;
            acc_z_positives = acc_z_fin > 0;
            
            changes_x = xor(acc_x_positives(1:end-1),acc_x_positives(2:end));
            changes_y = xor(acc_y_positives(1:end-1),acc_y_positives(2:end));
            changes_z = xor(acc_z_positives(1:end-1),acc_z_positives(2:end));
            
            % sum of positive and negative peaks
            acc_num_x = sum(changes_x);
            acc_num_y = sum(changes_y);
            acc_num_z = sum(changes_z);

            row_x = find(changes_x == 1);
            row_y = find(changes_y == 1);
            row_z = find(changes_z == 1);

            acc_Mov_x=diff(acc_x_fin(row_x));
            acc_Mov_y=diff(acc_x_fin(row_y));
            acc_Mov_z=diff(acc_x_fin(row_z));
            
            % mean peak 2 peak difference
            acc_x_p2p_mean=mean(abs(acc_Mov_x));
            acc_y_p2p_mean=mean(abs(acc_Mov_y));
            acc_z_p2p_mean=mean(abs(acc_Mov_z));
            
            % mean regularity
            acc_regul_x_mean=mean(acc_regul_x_fin);
            acc_regul_y_mean=mean(acc_regul_y_fin);
            acc_regul_z_mean=mean(acc_regul_z_fin);

        else
            accX_mean = NaN;
            accY_mean = NaN;
            accZ_mean = NaN;
            
            accX_max = NaN;
            accY_max = NaN;
            accZ_max = NaN;

            accX_min = NaN;
            accY_min = NaN;
            accZ_min = NaN;

            acc_num_x = NaN;
            acc_num_y = NaN;
            acc_num_z = NaN;
            
            acc_x_p2p_mean = NaN;
            acc_y_p2p_mean = NaN;
            acc_z_p2p_mean = NaN;
            
            acc_regul_x_mean=NaN;
            acc_regul_y_mean=NaN;
            acc_regul_z_mean=NaN;
            
        end
    else
         accX_mean = NaN;
         accY_mean = NaN;
         accZ_mean = NaN;
            
         accX_max = NaN;
         accY_max = NaN;
         accZ_max = NaN;

         accX_min = NaN;
         accY_min = NaN;
         accZ_min = NaN;

         acc_num_x = NaN;
         acc_num_y = NaN;
         acc_num_z = NaN;
            
         acc_x_p2p_mean = NaN;
         acc_y_p2p_mean = NaN;
         acc_z_p2p_mean = NaN;
         
         acc_regul_x_mean=NaN;
            acc_regul_y_mean=NaN;
            acc_regul_z_mean=NaN;
    end

    