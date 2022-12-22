%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   Master file for Optogenetics project                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The optogenetics project master is designed to import and organise spike
% sorted data generated from Plexon Offline Sorter 3.3.5. The program will
% also import and sort stimulus timing data from Intan .rhd files. 
% With the selection menu you can then run several functions to generate 
% PSTHs, Raw Traces, Response Time Course Functions, Raster Plots (All TBD) 
% and generate associated plots. 
% The "Data" struct stores all relavent .csvs and experimental data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Initilisation Settings                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; %clear command window
%clear; %% clear workspace

disp('-- Optogenetics Master 2022 --');
warning('off', 'MATLAB:Axes:NegativeLimitsInLogAxis');
warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
set(0, 'DefaultTextInterpreter', 'None')
set(0, 'DefaultLegendInterpreter', 'None')        % Don't interpret legend
set(0, 'DefaultFigureWindowStyle', 'normal');     % Don't plot in the dock.
set(0, 'DefaultAxesXMinorTick', 'on');     % Show minor tick marks
set(0, 'DefaultAxesYMinorTick', 'on');     

workspace; % ensure workshop panel is shown
mainpath = pwd; % ensure mainpath is working directory
addpath '.\functions'

%%%%%%%%%%%%%%%%%%%%%% Sets up basic import variables 
irf = 0; %placehold init variabl
DataIDXinit = 0;

%%%%%%%%%%%%%%%%%%%%%% Sets up containers 
FileIDX = {}; %setup FileIDX
DataIRF = struct();

%%%%%%%%%%%%%%%%%%%%%% PROGRAM SETTINGS

experiment = "Amogus"; % default name of experiment

extlen = 4; %length of file extension, just a hack until it can be autodet
export = 1; %export data
importwaveform = 1; %choose whether to import waveform or not.
import_rhd = 1; % Importing RHDS can make file very large, but are needed.
resample = 100; % amount to downsample trace by (file too large else)

%%%%%%%%%%%%%%%%%%%%%% Intan Channels. 
% Add 1 to name. E.g channel A0 is 1, A26 = 26. B0 = 33, B2=35 etc up to 64
channels = [17,18,19, 20, 23, 29, 30, 31, 32, 33, 34, 35, 36];
channels_wt = [19, 20, 23, 29, 30, 31, 32, 33, 34, 35, 36];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                         CHOICE 1 - What to Load                         %
%                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
txt = 9;
repeat = 0;
while txt > 4
    disp('- Data Import Menu -');
    prompt =...
        "[1] - Generate new database \n" + ...
        "[2] - Load database \n" + ...
        "[3] - Use Loaded Data \n" + ...
        "[4] - Close program \n" + ...
        ">> ";
    txt = input(prompt);
    if isempty(txt) %%% sends back to menu 
        txt = 5;
    end    
    if txt == 4 %%% End without saving, messy by easy
        error('Program End')
    end  
    clc;
end
choice = txt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Choose what to import                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1
    disp('How many experiments are you importing / processing?')
    prompt =...
    "WARNING: You must place experiment(s) within a subfolder \n" + ...
    "- Folder Structure: TOPFOLDER->EXPERIMENT(S)->STIMTYPE->Data/Subfolders"...
    +" \n"+ "- If you have no/one stim type, just make one stimtype folder"... 
    +"(e.g Light/Yellow Light/Base or other name.) \n" + ...
    "[1] - One (select its top level folder) \n" + ...
    "[2] - Multiple (select top level folder contaning experiments) " + ...
    "\n" + ">> ";
    exp_num = input(prompt);
    if isempty(exp_num)
        exp_num = '1';
    end
    if exp_num > 2
        exp_num = '2';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Select & Index Experimental Folder              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1
disp('Indexing Files & Folders. Select top level experiments folder');
% The main dataindex function
[FolderIndex, FolderCount, experiment_list, StimList] = dataindex(); 
% Set top level folder
TopLevelFolder = char(FolderIndex(1,1));
% Index every file
FullFileList = dir(fullfile(TopLevelFolder, '**\*.xls'));  
disp('Indexing Done');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Set Experiment Name                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1 && exp_num == 1 %% If you chose only 1
    %experiment_list = experiment_list(1);
    %disp('Building Database: '+string(experiment_list(1)))
end

if choice == 1 && exp_num >= 2 %% If you chose only 1 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Import & Organise Experimental Data             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 1
disp('Starting Data Import...');
Data = {}; %%%% Setup Data struct
DataWaveform = {}; %%%  Setup Waveform struct
DataStimulus = {}; %%%% Setup Data
DataIDX = struct();

%%%%%%%%%% FOR # EXPERIMENTS
for n_exp = 1 : length(experiment_list) % Sets import 
temp_experiment = convertCharsToStrings(genvarname(experiment_list(n_exp)));
ExperimentFolder = char(FolderIndex(1))+"\"+char(experiment_list(n_exp));
%FileList = dir(fullfile(ExperimentFolder, '**\*.xls'));  %FileList
disp('Building Database: '+string(experiment_list(n_exp)))
%%%%%%%%%% FOR # STIMS 
for n_stimtypes = 1 : length(StimList) %assumes equal accross all exps
    temp_stim = string(StimList(n_stimtypes));
    StimFolder = ExperimentFolder+"\"+string(StimList(n_stimtypes));
    FileList = dir(fullfile(StimFolder, '**\*.xls'));  %FileList
    %create a combined struct for later use
    Data.Combined.(temp_stim) = struct();
    DataStimulus.Combined.(temp_stim) = struct();
    disp('-> Building: '+temp_stim)
    disp("Starting XLS Import");
%%%%%%%%%% FOR # FILES
for n_files = 1 : length(FileList) % Iterate FileList from indexing
    temp_field = FileList(n_files,1); %create temp var to store targ file
    temp_name = convertCharsToStrings(temp_field.name(1:end-extlen)); %-.xls
    disp("Importing... "+temp_name);
    temp_field = convertCharsToStrings(temp_field.folder)+'\'...
        +convertCharsToStrings(temp_field.name); %format field name
    %data import func
    [temp_data, waveform_avg, waveforms] = dataimport(temp_field, importwaveform); 
    disp("Imported "+temp_name);
    %place data in struct
    Data.(temp_experiment).(temp_stim).(genvarname(temp_name)) = temp_data;
    %place waveformdata in struct (first one seems uneeded)
    %DataWaveform.(temp_experiment).(temp_stim).( ...
    %    genvarname(temp_name)) = waveforms; 
    DataWaveform.(temp_experiment).(temp_stim).( ...
        genvarname(temp_name+'_avg')) = waveform_avg; 
end %repeat until all files are placed in struct

if import_rhd >= 1
disp('.RHD Stimulus Import in progress...');
for n2 = 1 : length(FileList) % Iterate FileList from indexing
    try
    temp_field = FileList(n2,1); %create temp var to store targ file
    temp_file = convertCharsToStrings(temp_field.name(1:end-extlen)); %-.xls
    temp_trace_nm = temp_file;
    temp_name = temp_file+"_stimulus";
    temp_file = temp_file+'.rhd';
    temp_path = temp_field.folder+"\";
    [temp_data, trace] = datastimimport(temp_file, temp_path); %stim func
    disp("Imported "+temp_name);
    DataStimulus.(temp_experiment).(temp_stim).(genvarname(temp_name)) = temp_data; %place data in struct
    DataStimulus.(temp_experiment).(temp_stim).(genvarname(temp_trace_nm+'_trace')) = trace; %place data 
    catch
        disp("Warning: RHD file missing, too large to import or other error. Skipping...")         
    end
end %repeat until all files are placed in struct
end
%%%%%%%%%%%%% DataIDX files
    if DataIDXinit == 1
        DataIDX = [DataIDX; FileList];
    elseif DataIDXinit == 0
        DataIDX = FileList;
        DataIDXinit = 1;
    end
end %%% End of stim type import

end %%%% End of multi-experiment import

disp('Importing Done');

end %%% End of import
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Load Data                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 2 %Import DB and reconstruct workspace as if newely generated
    [file,path] = uigetfile({'*.mat'},'Select Data', ...
        '..\OptoPipeline\data\databases\');
   disp('Loading '+string(file)+'...')
    filepath = append(path,file);
    Data=load(filepath);
    DataStimulus=Data.db.DataStimulus;
    DataIRF=Data.db.DataIRF;
    DataWaveform=Data.db.DataWaveform;
    experiment_list = Data.db.experiment_list;
    FolderIndex=Data.db.FolderIndex;
    FileList=Data.db.FileList;
    TopLevelFolder=Data.db.TopLevelFolder;
    StimList=Data.db.StimList;
    Data=Data.db.Data;
    %Data=Data.data;
    disp('Data Import Done');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Use loaded data                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice == 3 %Import DB and reconstruct workspace as if newely generated
    disp("Using loaded data...")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                         CHOICE 2 - What to generate             
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice2 = 0;
while choice2 < 7
    if import_rhd <1
        disp("Warning! RHD data was not imported. Trace[2],"+...
        "Stacked Trace[3] and split PSTH[4] cannot be generated " +...
        "without stimulus data.")
    end
    
disp('- Data Analysis Menu -');
prompt =...
    "[1] - Stimulus Response \n" + ...
    "[2] - Single Channel Trace w/ stimulus \n" + ...
    "[3] - Stacked Multi-Channel Trace \n" + ...
    "[4] - Split Unit PSTH \n" + ...
    "[5] - Raster Plot - TBD \n" + ...
    "[6] - Response Intensity Boxplot \n" + ...
    "[7] - Save and Exit \n" + ...
    ">> ";
choice2 = input(prompt);
if isempty(choice2)
    choice2 = '7'; %automically close program if nothing selected
end
if choice2 >= 7
    outputProcessedData(Data, DataStimulus, FolderIndex, FileList, ...
        TopLevelFolder, experiment, DataIRF, DataWaveform, experiment_list, ...
        StimList) 
end
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Generate Intensity Reponse Function             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% REQUIRES SPIKE DATA
if choice2 == 1 %
    irf = 1;
    disp('Intense Response Function processing done.')
    analysis_type = 1;
    DataIRF = struct();

    %%%%%%%%%% Select specific protocol
    if analysis_type == 3
        FileIDX = fieldnames(Data);
        disp("Select Experiment (1 = First in list, 2 = second etc...")
        disp(FileIDX)
        prompt = ">>";
        irf_exp = string(FileIDX(input(prompt))); 
        
        FileIDX = fieldnames(Data.(irf_exp));
        disp("Select Stimtype (1 = First in list, 2 = second etc...")
        disp(FileIDX)
        prompt = ">>";
        irf_stim = string(FileIDX(input(prompt))); 
        
        FileIDX = fieldnames(Data.(irf_exp).(irf_stim));
        disp("Select Protocol (1 = First in list, 2 = second etc...")
        disp(FileIDX)
        prompt = ">>";
        irf_protocol = string(FileIDX(input(prompt))); 
    end
         
    %if protocol > length(FileIDX)
    %protocol = length(FileIDX);
    %end
    %protocol_nm = FileIDX(protocol);
    
    %%%%%%%%%% Select plot type
    if analysis_type == 3
    disp("Choose intensity response function type...")
    prompt =...
        "[1] - Plot Saturated Only \n" + ...
        "[2] - Plot Unaturated Only \n" + ...
        "[3] - Plot Both! \n" + ...
        ">> ";
    cutsat = input(prompt);
    if cutsat == 1
        output_file = (irf_protocol+" Intensity Response Function - Saturated")+'.jpg';
    end
    if cutsat == 2
        output_file = (irf_protocol+" Intensity Response Function - Unaturated")+'.jpg';
    end
    if isempty(cutsat) || cutsat == 0 || cutsat == 3
        cutsat = '3';
        output_file = (irf_protocol+" Intensity Response Function - Combined")+'.jpg';
    end
    end
    %%%%%%%%%%%%%%%%%%%% ITERATE EXPERIMENTS
    for n_exp = 1 : length(experiment_list) 
        irf_current_exp = string(experiment_list(n_exp));
        disp('Processing '+irf_current_exp)
    %%%%%%%%%%%%%%%%%%%% ITERATE STIMTYPES
    for n_stimtypes = 1 : length(StimList)
        irf_current_stim = string(StimList(n_stimtypes));
        FileIDX = fieldnames(Data.(irf_current_exp).(irf_current_stim));
        disp('Processing '+irf_current_stim)
    %%%%%%%%%%%%%%%%%%%% ITERATE FILES   
        disp("Processing per unit...")
    for irf1 = 1 : length(FileIDX) % Iterate FileList from stim timings
        try
        irf_file = string(FileIDX(irf1));   
        disp('IRF - ' + irf_file + ' processing...')
        irf_spikes = Data.(irf_current_exp).(irf_current_stim).(irf_file)(:,1:3);
        irf_stimfile = irf_file + "_stimulus";
        irf_stimfile = DataStimulus.(irf_current_exp).(irf_current_stim).(irf_stimfile);
        [stimtable] = intensity_response_function( ... % IRF Function
            irf_file, irf_spikes, irf_stimfile, channels, irf_current_exp);
        DataIRF.(irf_current_exp).(irf_current_stim).(genvarname(irf_file)) = stimtable; % Add data to struct
           disp('IRF - ' + irf_file + ' processed.')
        catch
           disp('IRF - ' + irf_file + ' failed to process. Error in file or stimfile missing.')
        end
    end %repeat for all files

    disp("Processing combined...")
    for irf2 = 1 : length(FileIDX) % Iterate FileList from stim timings
        try
        irf_file = string(FileIDX(irf2)); 
        disp('IRF - ' + irf_file + ' processing...')
        irf_spikes = Data.(irf_current_exp).(irf_current_stim).(irf_file)(:,1:3);
        irf_stimfile = irf_file + "_stimulus";
        irf_stimfile = DataStimulus.(irf_current_exp).(irf_current_stim).(irf_stimfile);
        [stimtable] = intensity_response_function_combined( ... % IRF 2
            irf_file, irf_spikes, irf_stimfile, channels);
        DataIRF.(irf_current_exp).(irf_current_stim).(irf_file).Combined = stimtable; % Add to struct
           disp('IRF - ' + irf_file + ' processed.')
        catch
           disp('IRF - ' + irf_file + ' failed to process. Error in file or stimfile missing.')
        end 
    end %repeat for all files
    end
    end
    disp('Processing done.')

    %%%%%%%%%%%%%%%% Plot averages protocol                                                                           
    if analysis_type >= 2                                                                                             
    y2 = [];                                                                                                          
    labels = [];                                                                                                      
    for protocol = 1 : length(FileIDX)                                                                                
        [x, y] = irf_plot(DataIRF, FileIDX, channels, ...                                                                 
        cutsat, protocol);                                                                                                
        avg_file = string(FileIDX(protocol));                                                                             
        y = mean(y, 2);                                                                                                   
        y2 = [y2, y];                                                                                                     
        labels = [labels, avg_file];                                                                                      
    end                                                                                                               
    ymean = mean(y2, 2);                                                                                              
    title('Summary of Intensity Response Function')                                                                   
    xlabel('Stimulus Intensity'); % adds x and y labels                                                               
    ylabel('Spike Count (mean)');                                                                                       

    %%%%%%%%%%%%%%%% Main Plot                                                                                        
    hold on                                                                                                           
                                                                                                                      
    if analysis_type == 1                                                                                             
    yyaxis left                                                                                                       
    plot(x,y2, 'LineWidth', 1.5)                                                                                      
    yyaxis right                                                                                                      
    plot(x, ymean, 'LineWidth', 3.0)                                                                                  
    ylim([0, 25])                                                                                                     
    %errorbar(ymean)                                                                                                  
    end                                                                                                               
    if analysis_type >= 2                                                                                             
    plot(x,y2, 'LineWidth', 2.0)                                                                                      
    legend(labels)                                                                                                    
    end                                                                                                               
    xlim([0.2, 1.6]);                                                                                                 
    grid on                                                                                                           
    hold off                                                                                                          
                                                                                                                      
    disp('Exporting plot to folder...');                                                                              
                                                                                                                      
    export = fullfile('figures', output_file');                                                                       
    exportgraphics(gcf, export, 'ContentType', 'vector')                                                              
                                                                                                                                                                
    if choice2 == 6                                                                                                   
        outputProcessedData(Data, DataStimulus, DataIndex, FileList, ...                                                  
        TopLevelFolder, experiment, DataIRF)                                                                              
    end                                                                                                                                                                                                                             
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Generate Raw Trace                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% REQUIRES RHD DATA
if choice2 == 2 %
    targetchannel = "Channel 35";
    %title('Channel Trace')
    xlabel('Time (seconds)'); % adds x and y labels
    %xlim([0.2, 1.6]);
    x=DataStimulus.DTA5.NOFILTER.SQ100S_stimulus.Time;
    x2=DataStimulus.DTA5.NOFILTER.SQ100S_trace.trace1;
    y1=DataStimulus.DTA5.NOFILTER.SQ100S_trace.(targetchannel);
    y2=DataStimulus.DTA5.NOFILTER.SQ100S_stimulus.Stimulus;

    yyaxis left
    plot(x2,y1)
    xlim([21 35])
    ylabel('Microvolts (µV)')

    yyaxis right;
    bar(x, y2);
    xlim([21 35])
    ylabel('Stimulus Intensity');
    ylim([0, 12]);
    %set(gca, 'YTick', []);
   export = fullfile('figures', 'Stimulus Trace.jpg');
   exportgraphics(gcf, export, 'ContentType', 'vector')
    disp('Generating Raw Trace')
   
   if choice2 == 6
   outputProcessedData(Data, DataStimulus, FolderIndex, FileList, ...
    TopLevelFolder, experiment, DataIRF) 
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Stacked Multi Channel Trace                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% REQUIRES RHD DATA
if choice2 == 3 %
    stackedplottable = table();
    for stacks = 1 : length(channels)
        targetchannel = channels(stacks);
        targetchannel = "Channel "+targetchannel;
        temp_stack_tbl = table(DataStimulus.DTA5.NOFILTER.SQ100S_trace.(targetchannel));
        temp_stack_tbl = renamevars(temp_stack_tbl,['Var1'],[targetchannel]);
        stackedplottable = [stackedplottable, temp_stack_tbl];
    end
    temp_stack_tbl = table(DataStimulus.DTA5.NOFILTER.SQ100S_stimulus.Stimulus);
    temp_stack_tbl = renamevars(temp_stack_tbl,'Var1','Stimulus');
    stackedplottable = [stackedplottable, temp_stack_tbl];
    disp('Generating Stacked Channel Plot')
    stackedplot(stackedplottable,"Title","Multichannel Trace")

   export = fullfile('figures', 'Stacked Multichannel Trace.jpg');
   exportgraphics(gcf, export, 'ContentType', 'vector')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Split Unit PSTH Graph                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% REQUIRES SPIKE AND RHD DATA
if choice2 == 4 %
   targetchannel = "Channel 36";
   disp('Generating Unit-PSTH')
   %target = char(FileIDX(2));
   histdata = [Data.DTA5.NOFILTER.SQ100S(Data.DTA5.NOFILTER.SQ100S.Channel == 36,:)];

    x=DataStimulus.DTA5.NOFILTER.SQ50S_stimulus.Time;
    x2=DataStimulus.DTA5.NOFILTER.SQ100S_stimulus.Time;
    y1=DataStimulus.DTA5.NOFILTER.SQ50S_trace.(targetchannel);
    y3=DataStimulus.DTA5.NOFILTER.SQ100S_trace.(targetchannel);
    y2=DataStimulus.DTA5.NOFILTER.SQ50S_stimulus.Stimulus;

   histdata_u1 = histdata(histdata.Unit == 1,:);
   histdata_u2 = histdata(histdata.Unit == 2,:);
   
   t = tiledlayout(2,2);
   t.Padding = 'compact';
   t.TileSpacing = 'compact';
   title(t,'Channel Unit PSTH')

   nexttile
  % ax1 = nexttile;
   %histogram(histdata_u1.Timestamp, 210);
   plot(x, y1)

   nexttile
   plot(DataWaveform.(target).SQ50S_avg)

   nexttile
   plot(x2, y3)
   nexttile
   plot(DataWaveform.(target).SQ100S_avg)

   export = fullfile('figures', 'Split Unit PSTH.jpg');
   exportgraphics(gcf, export, 'ContentType', 'vector')
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Rasterplot                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% REQUIRES SPIKE DATA & IRF TIMING (OR GET CODE FROM)
if choice2 == 5 %

    spike_struct = struct();
    %%%%%%%% stimtimes
    raw_spikes = Data.DTA5.NOFILTER.SQ30MS_50REP;
    spikes = [raw_spikes(raw_spikes.Channel == 46,:)];
    spikes = spikes(spikes.Unit == 1,:);
    %spikes = spikes(:,"Timestamp");
    trial_list = [];

    t1 = [1.9101;1.9399]; %0.0298 = 0.03 (599 rows, 1.9785 gap)
    t2 = [3.9185;3.9494]; % 0.03 (601 rows) 
    t3 = [5.9285;5.9584];
    t4 = [7.9629;7.9919];
    t5 = [9.9683;9.9983]; % 2400
    t6 = [11.9885;12.0184];
    t7 = [14.02195;14.0519];
    t8 = [16.0284;16.0583];%x,4800
    t9 = [18.0720;18.1019]; %4801, 5400
    t0 = [20.0786;20.1085];  %5401, 6001

    t11 = [22.0984;22.1284]; %0.0298 = 0.03 (599 rows, 1.9785 gap)
    t12 = [24.122;24.1519]; % 0.03 (601 rows) 
    t13 = [26.1386;26.1685];
    t14 = [28.182;28.2119];
    t15 = [30.1885;30.2185]; % 2400
    t16 = [32.2085;32.2384];
    t17 = [34.2319;34.2618];
    t18 = [36.2486;36.2786];%x,4800
    t19 = [38.2925;38.3224]; %4801, 5400
    t20 = [40.2989;40.3288];  %5401, 6001


    t = table(t1, t2, t3, t4, t5, t6, t7, t8, t9, t0, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20);

    for n = 1 : width(t)
    name = string(n);
    min = table2array(t(1,n));
    max = table2array(t(2,n));
    idx = spikes.Timestamp > (min-0.1) & spikes.Timestamp < (max+0.07);
    spiketable = spikes(idx,:);
    spiketable = spiketable.Timestamp-min;
    spiketable(:,2)=1;
    spike_struct.(genvarname('Trial '+name)) = spiketable;
    trial_name = ('Trial'+name);
    trial_list = [trial_list; trial_name];
    end

   DataRAST = spike_struct; % Add data to struct
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Response Intensity Boxplot                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% REQUIRES SPIKE DATA
if choice2 == 6 %
   
   disp('Not yet implemented')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Clean up                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                        Endings                                          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                        Local Functions (ye shall not pass)              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_out, y_out, labels] = irf_plot(DataIRF, FileIDX, channels, ...
    cutsat, protocol)

y_out = [];
prelabel = [];
postlabel = [];

for chan_count = 1 : length(channels)
    targetchannel = channels(chan_count);
    targetchannel = "Channel_"+targetchannel;
    channelunit = fieldnames(DataIRF.(string(FileIDX(protocol,1))).( ...
        targetchannel));

    %%% Set X range based on stimulus
    x_out = DataIRF.(string(FileIDX(protocol,1))).(targetchannel).( ...
        string(channelunit(1,1))).("Amplitude")(7:12);
    
    %%% Set Y plots based on units in channel
    for unitcount = 1 : length(channelunit)
        targetunit = channelunit(unitcount,:);
        y = DataIRF.(string(FileIDX(protocol,1))).(targetchannel).( ...
            string(targetunit)).SpikeCount(7:12);
        if max(y)*0.8 > y(length(y)) && cutsat == 1
            y_out = [y, y_out];
            prelabel = [prelabel; (targetchannel+" "+channelunit)];
        elseif max(y)*0.8 < y(length(y)) && cutsat == 2
            y_out = [y, y_out];
            prelabel = [prelabel; (targetchannel+" "+channelunit)]; 
        elseif cutsat >= 3
            y_out = [y, y_out];
            prelabel = [prelabel; (targetchannel+" "+channelunit)];  
        end
    end
    labels = [postlabel, prelabel] ;
    % = [postlabel, channelunit];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% TBD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputProcessedData(Data, DataStimulus, FolderIndex, FileList,...              
TopLevelFolder, experiment, DataIRF, DataWaveform, experiment_list, StimList)    
disp('Do you wish to save a database and code backup?');  
prompt =...
    "[1] - Save Database & Sourcecode \n" + ...
    "[2] - Save Sourcecode only \n" + ...
    "[3] - Exit without saving \n" + ...
    ">>";
saveexit = input(prompt);
if saveexit == 1 || saveexit == 2
disp('Saving '+string(experiment))

% Output Data and backup of script                                                   
estr = sprintf('%s.m', mfilename);                                          
FID = fopen(estr);                                                          
sourceCode = textscan(FID, '%s', 'delimiter', '\n');                        
fclose(FID);                                                                
                                                                            
timeStamp=datestr(now);                                                     
db.name = experiment;                                                       
db.Data = Data;
db.DataStimulus = DataStimulus;
db.DataWaveform = DataWaveform;
db.DataIRF = DataIRF;
db.experiment_list = experiment_list;
db.FolderIndex = FolderIndex;
db.FileList = FileList;
db.StimList = StimList;
db.TopLevelFolder = TopLevelFolder;
                                                                                                                                         
data.timeStamp=timeStamp;                                                   
data.sourceCode=char(sourceCode{1}); % Source code as plain text                                                  
data.mfilename= mfilename('fullpath');                                     
                                                                           
fprintf('\n%s\n',repmat(5786,88,1)); % This symbol                         
dataDir = fullfile('data', 'backups');                                     
dataDir_db = fullfile('data', 'databases'); 

%%%%% Create Backup & Data Directory (if none exists)                      
if ~exist(dataDir, 'dir')                                                  
fprintf('\n %s Creating directory: %s\n', 10171, dataDir)                 
mkdir(dataDir)                                                             
end                                                                         
if ~exist(dataDir_db, 'dir')                                                
fprintf('\n %s Creating directory: %s\n', 10171, dataDir_db)                
mkdir(dataDir_db)                                                           
end                   

%%%% Create Backup File                                                     
outName = ['backup_' (datestr(now, 'yyyy-mmmmdd-HHAM'))];                   
fprintf('\n %s Saving: %s\n\n', 10171, outName)  

if saveexit == 1
outNameDB = strcat(experiment+"_database");                                     
fprintf('\n %s Saving: %s\n\n', 10171, outNameDB)                           
save(fullfile(dataDir_db, outNameDB), 'db', '-v7.3') 
end
save(fullfile(dataDir, outName), 'data', '-v7.3')                                    


disp('Output Finished.'); 

else
disp('Cleaning up variables...');  
end
clearvars -except FileIDX Data DataStimulus DataIRF FileIndex FileList...
    TopLevelFolder channels FolderIndex x y2 DataWaveform choice2 experiment...
    DataIDX StimList experiment_list

disp('Program End.');                                      
end                                                                          
       