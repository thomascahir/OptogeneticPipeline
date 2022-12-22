function [out] = rasterdata(Data, DataIRF, DataStimulus, rast_current_exp, rast_current_stim, rast_file, FileIDX, rast1)

if startsWith(rast_file, "SQ1S") % set stim time (yes i know its inelegant)
    stimtime = 0.006;
    stimgap = 2.1;
elseif startsWith(rast_file, "SQ3S")
    stimtime = 0.010;
    stimgap = 2.1;
elseif startsWith(rast_file, "SQ5S")
    stimtime = 0.005;
    stimgap = 2.1;
elseif startsWith(rast_file, "SQ10S")
    stimtime = 0.01;
    stimgap = 2.1;
elseif startsWith(rast_file, "SQ30S")
    stimtime = 0.03;
    stimgap = 2.1;
elseif startsWith(rast_file, "SQ50S")
    stimtime = 0.05;
    stimgap = 2.1;
elseif startsWith(rast_file, "SQ100S")
    stimtime = 0.1;
    stimgap = 2.1;
elseif startsWith(rast_file, "SQ30MS")
    stimtime = 3.0;
    stimgap = 2.1; %might be incorrect for sino recordings
else
    stimtime = 3;
    stimgap = 2.1;
end
Rasttable = table();
%%%%%%%%%%%% GET INTENSITIES INTO ARRAY / TABLE
rast_file = string(FileIDX(rast1));   
disp('RAST - ' + rast_file + ' processing...')
rast_spikes = Data.(rast_current_exp).(rast_current_stim).(rast_file)(:,1:3);
%rast_stimfile = DataStimulus.(rast_current_exp).(rast_current_stim).(rast_stimfile);

%%%%% Gets onset time (in a very roundabout, inellgent way)
rast_trial = string(fieldnames(DataIRF.(rast_current_exp).(rast_current_stim)));   
rast_trial = rast_trial(rast1);

rast_channels = string(fieldnames(DataIRF.(rast_current_exp).(rast_current_stim).(rast_trial)));
n = length(rast_channels)-1;

rast_unit = DataIRF.(rast_current_exp).(rast_current_stim).(rast_trial).(rast_channels(n));
rast_unit = struct2table(rast_unit);
rast_onset = rast_unit.Unit_1.OnsetTime;


%rast_onset = rast_unit.Unit_0.OnsetTime;


for n2 = 1 : length(rast_onset)
idx = (rast_spikes{:,3} >= rast_onset(n2)) & (rast_spikes{:,3} <= rast_onset(n2)+(stimtime*3));
rasttable = rast_spikes(idx,3);
rasttable = rasttable.Timestamp - (rast_onset(n2));
rasttable = table(rasttable);
Rasttable = [Rasttable; rasttable];
end
pre_out = Rasttable;
pre_out = sortrows(pre_out);
pre_out = renamevars(pre_out,"rastyable","Spikes");
%pre_out(:,2) = 1:size(pre_out); %adds time value to each row
out = pre_out;
end % Function End