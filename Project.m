clc
clear 
%% Reading data

data = importdata('BreathingPeriod_6s.txt');
%data = importdata('TestAbbraGiulia.txt');
sampling_frequency = 405;
%remove time intervals from the matrix (we use sample frequency) and first
%and last 5 seconds
data(:,4) = [];
data = data(sampling_frequency*5 : size(data,1) - sampling_frequency*5 ,:);
x = data(:,1);
y = data(:,2);
z = data(:,3);

length = size(data,1);
n = (0:1:length-1);

plot(n, data);

%% BREATHING

% FROM x AXIS
% PreProcessing:
% 4th order butterworth lp filter cutoff freq 0.5 Hz

cutOffFrequency = 0.5;
[b,a] = butter(4, cutOffFrequency/ (sampling_frequency/2));
%freqz(b,a)

xFiltered = filter(b,a,x);

%plot(n,xFiltered);

%% MINIMUN DETECTION

%Separate al least 3 seconds

Mins = islocalmin(xFiltered,'MinSeparation',3*sampling_frequency); 

plot(n,xFiltered,n(Mins),xFiltered(Mins),'r*')

CurrentBreathT = 0;
TotalBreathT = 0;
tot = 0;
start = find(Mins);
for i = start(1) + 1:start(end)
    CurrentBreathT = CurrentBreathT + 1/sampling_frequency;
    if Mins(i) ~= 0
        tot = tot + 1;
        TotalBreathT = TotalBreathT + CurrentBreathT;
        CurrentBreathT = 0;
        continue;
    end
end

BreathingPeriod = TotalBreathT / tot;

fprintf("Breathing rate estimation with Min Detection = %d \n", 1/BreathingPeriod);
%% Breathing rate estimation with PSD

t = 0:1/sampling_frequency:max(n)/sampling_frequency;

axis([0 1 -100 2]);
plomb(x,t,'psd');

[pxx, f] = plomb(x,t,'psd');
pxx = 10*log10(pxx);
[pk, f] = findpeaks(pxx,f);
[M,I] = max(pk);
fprintf("Breathing rate estimation with PSD = %d \n", f(I));


%% Heart Rate Threshold
[zeros,poles,gain] = butter(4,[5 25]/sampling_frequency);
[mg,ph] = zp2sos(zeros,poles,gain);
zFiltered = filtfilt(mg,ph,z);

%rectification?

zFiltered = abs(zFiltered);

%?

findpeaks(zFiltered);

[pk,~] = findpeaks(zFiltered, 'MinPeakDistance', 0.1*sampling_frequency);
dim = size(pk);

threshold = sum(pk) / (dim(1));   %avg peak

[peaks,f0] = findpeaks(zFiltered, 'MinPeakDistance', 0.1*sampling_frequency, ... 
                    'MinPeakHeight', threshold);

                
                
windowSize = 100; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% for each peack
%   compute the average in a window of 100 samples centered in the peak
%   compute cross correlation with original signal
%   if the max corr is more than 100ms lagged, discard else continue
% end


r = filter( b, a, zFiltered );

[c,lags] = xcorr(zFiltered,r);

% subplot(2,1,1);
% plot(r);
% subplot(2,1,2);
% plot(zFiltered);

% plot(zFiltered)
% hold on
% plot(f0, peaks, 'o');
% hold off

% Find avg beat T

e0 = 0;
totT = 0;
beatT = 0;
numB = size(f0); 

for i = 1:numB
    el = f0(i);
    beatT = ((el - e0) * 1/sampling_frequency);
    totT = totT + beatT;
    e0 = el;
end

bpm = 60/(totT/numB(1));

fprintf("Hart rate estimation = %d bpm\n", bpm);

% 1 battito : avgBeatT = x battiti : 60 sec

%% HeartRate Hilbertian transform
[zeros,poles,gain] = butter(4,[1 45]/sampling_frequency);
[mg,ph] = zp2sos(zeros,poles,gain);

xFiltered = filtfilt(mg,ph,x);

yFiltered = filtfilt(mg,ph,y);

zFiltered = filtfilt(mg,ph,z);

d = size(xFiltered);

% subplot(3,1,1);
% plot(xFiltered, 'g-');
% 
% subplot(3,1,2);
% plot(yFiltered, 'r-');
% 
% subplot(3,1,3);
% plot(yFiltered);

sn = (xFiltered.^2 + yFiltered.^2 + zFiltered.^2).^0.5;

shatn = hilbert(zFiltered);

At = (real(shatn).^2 + imag(shatn).^2).^0.5; % magnitude of its envelope curve

% plot(At)

ftrans = fft(At);

bpftrans = bandpass(ftrans, [0.5 3], sampling_frequency);

res = ifft(bpftrans);

plot(res);

% 
% subplot(2,1,1);
% plot(res)
% 
% subplot(2,1,2);
% plot(real(res))
% hold on
% plot(imag(res))
% hold off


%[zeros,poles,gain] = butter(10,[1 25]/sampling_frequency);
%[mg,ph] = zp2sos(zeros,poles,gain);
%yFiltered = filtfilt(mg,ph,y); 
%%

cose = bandpass(fft(zFiltered), [0.5 3], sampling_frequency);
coseres = ifft(cose);

plot(real(coseres))
hold on
plot(imag(coseres))
hold off

%% THreshold pt2

[zeros,poles,gain] = butter(5,[5 35]/sampling_frequency);
[mg,ph] = zp2sos(zeros,poles,gain);
zFiltered = filtfilt(mg,ph,z);

%1

HBI = sampling_frequency;
Seg_len = HBI/4;

peaks = nan(ceil(max(n)/Seg_len), 1);
pos = nan(ceil(max(n)/Seg_len), 1);

%2

prevI = 1;
nextI = Seg_len;
plot(zFiltered);
hold on
for i = 1 : ceil(max(n)/Seg_len) - 1
    [pk,f] = max(zFiltered(prevI : nextI ,:));
    f = f + (i-1)*Seg_len;
    peaks(i) = pk;
    pos(i) = f;
    prevI = prevI + Seg_len;
    nextI = nextI + Seg_len;
    plot(f,pk,'o');
end
hold off



%%
[pk,f] = max(zFiltered);

% plot(n,zFiltered,f,pk,'o')











