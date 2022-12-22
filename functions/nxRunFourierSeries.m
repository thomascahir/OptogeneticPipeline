
% nx_plot.nxRunFourierSeries
% 
% Calculate standard f0/f1/f2 for drifting grating data
% 
% Version History
%   V2.0   Calvin Eiber   06-Jul-17   Adapted from old nx_plot function

function [passes, blocks, options] = nxRunFourierSeries(passes, blocks, options)

% Set default options for calculation
options = setDef(options,'fold_to_cycle',1);
options = setDef(options,'bin_per_cycle',32);
options = setDef(options,'return_complex',false); % export complex-valued f0,f1,f2
options = setDef(options,'return_folded',false);  % export the folded PSTH? 

if strncmpi(options.time_unit,'s',1), tUnit = 1;    % seconds
else                                  tUnit = 1e-4; % EXPO sample rate
end


% nPasses = size(passes,1);
% nBlocks = size(blocks,1);

% Compute f0,f1,f2 for each pass
% we need to start by getting the folded PSTH 

TF_field = {'tf','Surface_driftrate','TempFreq'}; % possible names for temporal frequency
for ii = 1:length(TF_field) % Find the right one for this dataset
    if any(strcmp(blocks.Properties.VariableNames,TF_field{ii}))
        TF_field = TF_field{ii}; break
    end
end

if iscell(TF_field), error('Cannot do Fourier analysis without TemporalFrequency data'), end
block_TF = blocks.(TF_field);
pass_TF = block_TF(passes.Cond);

% This is done without loops through liberal usage of "cellfun" and friends

% filter times for spiketimes in stimulus block only
duration = num2cell(passes.EndTime - passes.StartTime); % Pass duration
times = cellget(@(t,d) t(t > 0 & t <= d), passes.SpikeTimes, duration);
tFold = num2cell(options.fold_to_cycle ./ pass_TF./tUnit); % duration to which to fold
times = cellget(@(t,w) mod(t,w)./w, times, tFold); % Folded spike phases (0-1)

nBins  = options.fold_to_cycle * options.bin_per_cycle;
psth_P = linspace(0,1,nBins+1); % PSTH phase bins (0-1)
psth_Y = cellget(@(t) hist(t,psth_P), times);
psth_Y = cellget(@(p) [p(1)+p(end) p(2:end-1)], psth_Y); % correct first and last bin

% Need to convert from counts to spike-rate in imp/s
nCycles = cellget(@(d,w) ceil(d/w), duration, tFold); % How many cycles or part-cycles? 
nCycles = cellget(@(n,w) bsxfun(@plus,psth_P(1:end-1)'*ones(1,n),0:(n-1))*w, nCycles, tFold); % expand out
nCycles = cellget(@(n,d) sum(n <= d,2)', nCycles, duration); % How many bins total in pass duration?
psth_Y = cellget(@(p,n,w) p./n./w./tUnit*nBins, psth_Y, nCycles, tFold);

% Now that we have an appropriately folded and scaled PSTH, compute FFT...
psth_F = cat(1,psth_Y{:}); % Flatten to matrix from cell array
psth_F = fft(psth_F,[],2)./nBins; % Compute FFT on matrix

pass.complex_f0 = psth_F(:,1); 
pass.complex_f1 = psth_F(:,1 + options.fold_to_cycle ); 
pass.complex_f2 = psth_F(:,1 + options.fold_to_cycle*2 ); 

pass.f0 = abs( pass.complex_f0 ); 
pass.f1 = abs( pass.complex_f1 ); 
pass.f2 = abs( pass.complex_f2 ); 

pass.f1phase = angle( pass.complex_f1 ); 
pass.f2phase = angle( pass.complex_f2 ); 

% Convert phase to degrees and add 90 to make it
% sine-phase normalised. [expo does this a funny way]
pass.f1phase =  pass.f1phase.*180/pi + 90;
pass.f2phase =  pass.f2phase.*180/pi + 90;

% Pass data is used for fitting, so always return it. 
% If the user didn't ask for it, it will get stripped off later.
if options.return_folded
    passes.folded_PSTH = psth_Y;
end
if options.return_complex
    passes.f0 = pass.complex_f0;
    passes.f1 = pass.complex_f1;
    passes.f2 = pass.complex_f2;
else
    passes.f0 = pass.f0;
    passes.f1 = pass.f1;
    passes.f2 = pass.f2;    
end
passes.phf1 = pass.f1phase;
passes.phf2 = pass.f2phase;


if options.return_folded % export average and sem of PSTH
    getPSTH = @(c) mean(cat(1,psth_Y{passes.Cond == c}),1);
    blocks.folded_PSTH = arrayfun(getPSTH,blocks.Cond,'UniformOutput',false);
    
    getPSTH = @(c) std(cat(1,psth_Y{passes.Cond == c}),[],1) ./ ...
                           sqrt(sum(passes.Cond == c));
    blocks.folded_PSTH_sem = arrayfun(getPSTH,blocks.Cond,'UniformOutput',false);
end

% The approach for averaging is:
% block_f1 = mean(abs(pass_f1))
% block_f1sem = std(abs(pass_f1)) / sqrt(len(...));
% block_f1phase = angle(mean(pass_f1))

the = @(f,c) pass.(f)(passes.Cond == c); % lambda to fetch field (e.g. f0) for a given cond
getBlockMeans = @(f) arrayfun(@(b) nanmean(the(f,b)), blocks.Cond);
getBlockSEM   = @(f) arrayfun(@(b) nanstd(the(f,b)) ./ sqrt(length(the(f,b))), blocks.Cond);
getBlockPhase = @(f) arrayfun(@(b) angle(nanmean(the(f,b))), blocks.Cond);

% Again, avoid loops by using arrayfun for f0/f1/f2. 
if options.return_complex
    blocks.f0 = getBlockMeans('complex_f0');
    blocks.f1 = getBlockMeans('complex_f1');
    blocks.f2 = getBlockMeans('complex_f2');
else
    blocks.f0 = getBlockMeans('f0');
    blocks.f1 = getBlockMeans('f1');
    blocks.f2 = getBlockMeans('f2');
end

blocks.f0sem = getBlockSEM('f0');
blocks.f1sem = getBlockSEM('f1');
blocks.f2sem = getBlockSEM('f2');
blocks.phf1 = getBlockPhase('complex_f1');
blocks.phf2 = getBlockPhase('complex_f2');

% Convert phase to degrees and add 90 to make it
% sine-phase normalised. [expo does this a funny way]
blocks.phf1 =  blocks.phf1.*180/pi +90;
blocks.phf2 =  blocks.phf2.*180/pi +90;

return




function foo = cellget(varargin)
% cellget(varargin) is an alias for cellfun(varargin,'UniformOutput',false)
% CDE 5-Apr-17 because I'm tired of writing "UniformOutput,false" heaps.
foo = cellfun(varargin{:},'UniformOutput',false);

function O = setDef(O,field,val)
if ~isfield(O,field), O.(field) = val; end
