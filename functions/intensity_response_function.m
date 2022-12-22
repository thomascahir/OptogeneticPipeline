function [stim_struct] = intensity_response_function( ...
    temp_name, raw_data_spikes, data_stim, channels, irf_current_exp)

%%%%%%%%%%%% EXTRACT TIMINGS FROM FILE NAMES
timdif = 0;

if startsWith(temp_name, "SQ1S") % set stim time (yes i know its inelegant)
    stimtime = 0.006;
    stimgap = 2.1;
    if contains(irf_current_exp, "WT") == true
        timdif = 0;
    end
elseif startsWith(temp_name, "SQ3S")
    stimtime = 0.010;
    stimgap = 2.1;
    if contains(irf_current_exp, "WT") == true
        timdif = 0;
    end
elseif startsWith(temp_name, "SQ5S")
    stimtime = 0.005;
    stimgap = 2.1;
    if contains(irf_current_exp, "WT") == true
        timdif = 95;
    end
elseif startsWith(temp_name, "SQ10S")
    stimtime = 0.01;
    stimgap = 2.1;
    if contains(irf_current_exp, "WT") == true
        timdif = 137;
    end
elseif startsWith(temp_name, "SQ30S")
    stimtime = 0.03;
    stimgap = 2.1;
    if contains(irf_current_exp, "WT") == true
        timdif = 180;
    end
elseif startsWith(temp_name, "SQ50S")
    stimtime = 0.05;
    stimgap = 2.1;
    if contains(irf_current_exp, "WT") == true
        timdif = 222;
    end
elseif startsWith(temp_name, "SQ100S")
    stimtime = 0.1;
    stimgap = 2.1;
    if contains(irf_current_exp, "WT") == true
        timdif = 265;
    end
end


%%%%%%%%%%%% GET INTENSITIES INTO ARRAY / TABLE
stim_struct = struct(); %Init stim data struct
unit_struct = struct(); %Init unit data struct
%raw_data_spikes = raw_data_spikes()
%%%%%%%%%% Setup arrays
Amplitude = [];
OnsetTime = [];
SpikeCount = [];

for n = 1 : length(channels)
    Amplitude = [];
    amp = 12;
    OnsetTime = [];
    ons = 1;
    SpikeCount = [];

    data_spikes = [raw_data_spikes(raw_data_spikes.Channel == channels(n),:)];
    data_spikes{:,3} = (data_spikes{:,3}-timdif);

    units = max(data_spikes.Unit);
    name = string(channels(n));
    name = string("Channel_"+name);
    idx = any(data_stim.Stimulus(:,:)>0,2); %index and extra only times where stim>0
    data_stim = data_stim(idx,:); %create temp data with only stimulus active]
    
    %%%%%%%%%% Setup variables for extraction


    %%%%%%%%%% Setup variables for extraction
    for n2 = 1 : units
    stimstart = data_stim.Time(1,:); % gets first stim trial
    %if contains(experiment, "WT") 
    %    stimstart = stimstart + addtime   
    %end
    stim_n = stimstart;
    stim_n2 = stim_n + stimtime ;%set high range for trial
    rowmin = 0;
    rowmax = 0;
    Stimulus = [0];
    
    %%%%%%%%%%% EXTRACT STIMULUS AMPLITUDES AND TIMES
    
        data_spikes = data_spikes(data_spikes.Unit == n2,:); %%% Borked?
        unitname = string(n2);
        unitname = string("Unit_"+unitname);
        while stim_n < max(data_stim.Time(:,1))
            try
            rowmin = rowmax;
            rowmax = rowmin + sum((data_stim.Time)>=stim_n & ...
                data_stim.Time<=stim_n2);
            rowmin = rowmin + 1;
        
            stim = max(data_stim.Stimulus(rowmin:rowmax)); 
            % gets stim max (assumed amp for stim)
            
            spikes = sum(data_spikes.Timestamp > stim_n & ...
                data_spikes.Timestamp < stim_n2);
        
            Amplitude = [Amplitude;stim];
            OnsetTime = [OnsetTime;stim_n];
            SpikeCount = [SpikeCount;spikes];
        
            if stim_n < max(data_stim.Time(:,1))
                stim_n = data_stim.Time(rowmax+1); %new n
        
            end
            stim_n2 = stim_n + stimtime; %set high range for trial
            catch
                stim_n = stim_n2*10000; %just to stop errors, spaghetti code i know
            end
        end
        stimtable = table(Amplitude,OnsetTime,SpikeCount);
        %stimtable = renamevars(stimtable, ...
            %["Var1", "Var2", "Var3"], ...
            %["Amplitude", "OnsetTime", "SpikeCount"]);
        if temp_name == "SQ30MS_50REP"
        stimtable = stimtable(1:50,:);
        else
        stimtable = stimtable(1:12,:);
        end
        unit_struct.(unitname) = stimtable;
        ons = amp*n2+1; % to iterate through start
    end
    %%%%%%%%%% Build intensity response data table
    
    stim_struct.(name) = unit_struct;
end


