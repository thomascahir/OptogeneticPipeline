%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fftdat,nh,harmdist] = vgbfft_2022(epoch, bins_per_cy, rate_per_bin, ~, nhist)

if ~exist('epoch', 'var') % Set up example data if called without arguments
    epoch = load_sample_epoch;
    bins_per_cy = 64;
    rate_per_bin = 0.0859;
    nhist = 1;
end

showFlag = 1;
h = epoch;
n = bins_per_cy;
r = rate_per_bin;
h = h./(n/2);					% change spike rates to account 
%                               % for number of bins in fft calculation
nh = reshape(h,n,length(h)/n);	% This gives Peak-to-trough amplitude of sine
Y=fft(nh,n); 	% FFT
% First three harmonics
Y0=Y(1,:);
Y1=Y(2,:);
Y2=Y(3,:);
Y3=Y(4,:);
Y4=Y(5,:);
Y5=Y(6,:);

% Amp, phase
amp0=Y0/2;
ph1 = angle(Y1);
amp1 = abs(Y1);
ph2 = angle(Y2);
amp2 = abs(Y2);

if nargout == 3
    % For harmonic distortion
    amp3 = abs(Y3);
    amp4 = abs(Y4);
    amp5 = abs(Y5);
    harmdist = sqrt(sum([(amp2/amp1)^2 (amp3/amp1)^2 (amp4/amp1)^2 (amp5/amp1)^2]));
end


hnum=(1:nhist)/1000;
fftdat = [hnum',amp0',amp1',ph1',amp2',ph2'];

end

function epoch = load_sample_epoch

  epoch = [11.6364
   11.6364
   58.1818
   34.9091
   23.2727
   23.2727
   11.6364
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
         0
   23.2727
         0
   11.6364
   11.6364
   11.6364
   23.2727
   23.2727
   69.8182
   23.2727
   23.2727
   69.8182
   93.0909
   69.8182
   58.1818
   34.9091
   23.2727
   23.2727
   58.1818
   11.6364
   46.5455
   46.5455
   34.9091
   23.2727
         0
   34.9091
   46.5455]'; % Note transpose to create column vector

end