function [out] = ireg_time(pos, time)
    % Use this function to sum up the time that the user spends doing a certain action. 
    % Useful if time intervals are irregular.
    % pos: positions of where the action is happening (e.g.
    % find(x.LeftSwitch==1)
    % time: time vector (e.g. x.Time_s_)
    one = zeros(max(pos)+1,1);
    one(pos)=1; 
    bin = diff(one); % -1 are starts, 1 are ends
    % this adds a start if the action begins at time 0
    if any(pos)==1 % checks that condition is not empty
        % this adds a start if the action begins at time 0
        startpos=find(bin==1);endpos=find(bin==-1);
        if any(startpos)==0
            if pos(1)==1
                bin(1)=1;
                %bin(1)=-1; % for datadivision.m
            end
        else
            if endpos(1)==1
                bin(1)=0;
            elseif startpos(1) > endpos(1)
                if pos(1)==1
                    bin(1)=1;
                    %bin(1)=-1; % for datadivision.m
                end
            end
        end

    time_start = time(find(bin==1));
    time_end = time(find(bin==-1));
    time_diff = time_end-time_start;
    out = sum(time_diff);
    end