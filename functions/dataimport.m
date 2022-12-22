function [data, waveform_avg, waveforms] = dataimport(input, importwaveform)
    if importwaveform > 0
        data = readtable(input);
        data = rmmissing(data); % CULLS EMPTY SPACES
        data = renamevars(data, ...
        ["Var1","Var2","Var3"], ...
        ["Channel","Unit","Timestamp"]);
        data = sortrows(data,["Channel" "Unit" "Timestamp"]);
    else
        data = readtable(input,'Sheet',1,'Range','A1:C65536');
        data = rmmissing(data); % CULLS EMPTY SPACES
        data = sortrows(data,["Channel" "Unit" "Timestamp"]);
    end
    rows = width(data(1,:));
    rows = rows-6; % to ignore non-wavelength columns (Above and PC0,1,2)
    waveforms = data(:,7:rows); %extract waveforms
    waveforms = table2array(waveforms); %convert to array
    waveform_avg = mean(waveforms); %convert to mean waveform

    clearvars -except data waveform_avg waveforms
end