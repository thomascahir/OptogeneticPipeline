function [stim_struct] = intensity_response_function_combined( ...
    temp_name, raw_data_spikes, data_stim, channels)

% SINO times. ~3second stim, ~2.5 to 2.6s between stims
% SQ Times. Variable stim time, ~2.0 second time between stims
% plan: could either 1) extra what the timings should be based on file
% name e.g if contains SQ10 = 10ms, 1HZ=3 etc, or maybe have it entered as
% a variable option prompt

%%%%%%%%%%%% EXTRACT TIMINGS FROM FILE NAMES

if startsWith(temp_name, "SQ1S") % set stim time (yes i know its inelegant)
    stimtime = 0.006;
    stimgap = 2.1;
elseif startsWith(temp_name, "SQ3S")
    stimtime = 0.010;
    stimgap = 2.1;
elseif startsWith(temp_name, "SQ5S")
    stimtime = 0.005;
    stimgap = 2.1;
elseif startsWith(temp_name, "SQ10S")
    stimtime = 0.01;
    stimgap = 2.1;
elseif startsWith(temp_name, "SQ30S")
    stimtime = 0.03;
    stimgap = 2.1;
elseif startsWith(temp_name, "SQ50S")
    stimtime = 0.05;
    stimgap = 2.1;
elseif startsWith(temp_name, "SQ100S")
    stimtime = 0.1;
    stimgap = 2.1;
elseif startsWith(temp_name, "x1HZ")
    stimtime = 3.0;
    stimgap = 2.1; %might be incorrect for sino recordings
else
    stimtime = 3;
    stimgap = 2.1;
end

%%%%%%%%%%%% GET INTENSITIES INTO ARRAY / TABLE
stim_struct = struct();
%data_stim = DataStimulus.SQ3S_stimulus;
for n = 1 : length(channels)
    data_spikes = [raw_data_spikes(raw_data_spikes.Channel == channels(n),:)];
    name = string(channels(n));
    name = string("Channel_"+name);
    idx = any(data_stim.Stimulus(:,:)>0,2); %index and extra only times where stim>0
    data_stim = data_stim(idx,:); %create temp data with only stimulus active]
    
    %%%%%%%%%% Setup variables for extraction
    stimstart = data_stim.Time(1,:); % gets first stim trial
    stim_n = stimstart;
    stim_n2 = stim_n + stimtime; %set high range for trial
    rowmin = 0;
    rowmax = 0;
    Stimulus = [0];
    
    %%%%%%%%%% Setup arrays
    Amplitude = [];
    OnsetTime = [];
    SpikeCount = [];
    
    %%%%%%%%%%% EXTRACT STIMULUS AMPLITUDES AND TIMES
    while stim_n < max(data_stim.Time(:,1))
        try
        rowmin = rowmax;
        rowmax = rowmin + sum(data_stim.Time>=stim_n & data_stim.Time<=stim_n2);
        rowmin = rowmin + 1;
    
        stim = max(data_stim.Stimulus(rowmin:rowmax)); % gets stim max (assumed amp for stim)
        
        spikes = sum(data_spikes.Timestamp > stim_n & data_spikes.Timestamp < stim_n2);
    
        Amplitude = [Amplitude;stim];
        OnsetTime = [OnsetTime;stim_n];
        SpikeCount = [SpikeCount;spikes];
    
        if stim_n < max(data_stim.Time(:,1))
            stim_n = data_stim.Time(rowmax+1);
    
        end
        stim_n2 = stim_n + stimtime; %set high range for trial
        catch
            stim_n = stim_n2*100; %just to stop errors, spaghetti code i know
        end
    end
    
    %%%%%%%%%% Build intensity response data table
    
    stimtable = table(Amplitude, OnsetTime, SpikeCount);
    stim_struct.(name) = stimtable;
end


