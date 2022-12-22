function [stimulus, trace] = datastimimport(file, path)
    %%%% Run Intan Read File function. Imports stim and trace data
    [stimulus, trace] = read_Intan_RHD2000_stimulus(file, path);
    import_activeonly = 1;
    import_downsampletrace = 1;
    resample = 20000; %how much to downsample by.

    %%%% Process stimulus data
    stimulus = stimulus-0.68; % removes sqaure pulse floor 
    stimulus(stimulus<0.01) = 0; % removes noise values
    stimulus(:,2) = 1:size(stimulus); %adds time value to each row
    stimulus(:,2) = stimulus(:,2)./(resample); % convert into seconds (downsamp 20kHZ)
    stimulus = fliplr(stimulus); % flip columns so 1=time, 2=stimulus times
    stimulus = array2table(stimulus); %convert to table
    stimulus = renamevars(stimulus, ...
        ["stimulus1","stimulus2"], ...
        ["Time","Stimulus"]); 
    if import_activeonly >=1 
        idx = any(stimulus.Stimulus(:,:)>0,2); %index and extra only times where stim>0
        stimulus = stimulus(idx,:); %create temp data with only stimulus active]
    end

    %%%% Process partial/mixed data (SQ1-100 files)
    %% Warning: This is a bespoke workaround to get these specific files to work.
    if contains(path, "WT") == true
        if contains(file, "SQ100S")
            a = "Found the WT SQ100S"
            idx = any(stimulus.Time(:,:)<26,2); %index and extra only times where stim>0
            stimulus = stimulus(idx,:); %create temp data with only stimulus active]
        end
        if contains(file, "SQ50S")
            a = "Found the WT SQ50S"
            idx = any(stimulus.Time(:,:)<26,2); %index and extra only times where stim>0
            stimulus = stimulus(idx,:); %create temp data with only stimulus active]
        end
    end

    %%%% Process trace data
    trace(:,1) = 1:size(trace); %adds time value, but does delete chanel 1
    if import_downsampletrace >= 1
        trace(:,1) = trace(:,1)./(resample);
        trace=downsample(trace, 200);
    end
    trace = array2table(trace); %convert to table
    %trace = renamevars(trace, trace1 ,"Time"); 
    allVars = 2:width(trace);
    newNames = append("Channel ",string(allVars));
    trace = renamevars(trace,allVars,newNames);
end


%Z1=transpose(board_adc_stimulus)
%Z3=Z1(Z1<=0.05)=0
%Z3(Z3>=0.05)=Z3
%Z3(:,2) = 1:size(Z3);
%Z4=downsample(Z3, 1000)