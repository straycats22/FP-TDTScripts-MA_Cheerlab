%% Housekeeping
% Clear workspace and close existing figures.
% Add data
clear all; clc;
data = TDTbin2mat('WIN78-FR1-6_WIN79-FR1-6', 'TYPE', {'epocs', 'scalars', 'streams'});

%% Remove Last reward Value  
data.epocs.bCu_.onset(end)=[];
data.epocs.bCu_.offset(end)=[];
data.epocs.bCu_.data(end)=[];
data.epocs.bCu_.onset(1)=[];
data.epocs.bCu_.offset(1)=[];
data.epocs.bCu_.data(1)=[];

% Extract rewarded nosepokes
%data.epocs.bRP_.onset = setdiff(data.epocs.bRP_.onset , data.epocs.bCu_.onset);

%% Setup the variables for the data you want to extract
% We will extract two different stream stores surrounding the 'PtAB' epoch 
% event. We are interested in a specific event code for the shock onset.

REF_EPOC = 'bCu/'; % Stimulation event to center on
% STIM_EPOC = aCu/'; %Stimulation pulses
STREAM_STORE1 = 'V2_B' % name of the 405 store
STREAM_STORE2 = 'B2_B' % name of the 470 store
TRANGE = [-5 35]; %window size [start time relative to epoc onset, entire duration]
% ARANGE = [1 1];

BASELINE_PER = [-5 -1]; % baseline period before stim
ARTIFACT405 = Inf;% variable created for artifact removal for 405 store
ARTIFACT470 = Inf;% variable created for artifact removal for 470 store
% Use TDTfilter to extract data around our epoc event
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event. Use the 'VALUES' parameter to specify allowed values of
% the REF_EPOC to extract.  For stream events, the chunks of data are 
% stored in cell arrays structured as data.streams.(STREAM_STORE1).filtered
data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);  
%data = TDTfilter(data, 'Stim', 'VALUES', ARANGE);

%%
% Optionally remove artifacts. If any waveform is above ARTIFACT level, or
% below -ARTIFACT level, remove it from the data set.
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT405), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT405), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);

art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT470), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT470), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
good2 = ~art1 & ~art2;
data.streams.(STREAM_STORE2).filtered = data.streams.(STREAM_STORE2).filtered(good2);

numArtifacts = sum(~good) + sum(~good2);

%%
% Applying a time filter to a uniformly sampled signal means that the
% length of each segment could vary by one sample.  Let's find the minimum
% length so we can trim the excess off before calculating the mean.
minLength1 = min(cellfun('prodofsize', data.streams.(STREAM_STORE1).filtered));
minLength2 = min(cellfun('prodofsize', data.streams.(STREAM_STORE2).filtered));
data.streams.(STREAM_STORE1).filtered = cellfun(@(x) x(1:minLength1), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false);
data.streams.(STREAM_STORE2).filtered = cellfun(@(x) x(1:minLength2), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false);

allSignals = cell2mat(data.streams.(STREAM_STORE1).filtered');

% downsample 10x and average 405 signal
N = 10;
F405 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength1 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
meanSignal1 = mean(F405);
stdSignal1 = std(double(F405))/sqrt(size(F405,1));
dcSignal1 = mean(meanSignal1);

% downsample 10x and average 465 signal
allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
F465 = zeros(size(allSignals(:,1:N:end-N+1)));
for ii = 1:size(allSignals,1)
    F465(ii,:) = arrayfun(@(i) mean(allSignals(ii,i:i+N-1)),1:N:length(allSignals)-N+1);
end
minLength2 = size(F465,2);

% Create mean signal, standard error of signal, and DC offset of 465 signal
meanSignal2 = mean(F465);
stdSignal2 = std(double(F465))/sqrt(size(F465,1));
dcSignal2 = mean(meanSignal2);

%% Plot Epoch Averaged Response
% Create the time vector for each stream store
ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;

% Subtract DC offset to get signals on top of one another
meanSignal1 = meanSignal1 - dcSignal1;
meanSignal2 = meanSignal2 - dcSignal2;

% Plot the 405 and 465 average signals
figure(1);
subplot(3,1,1);
plot(ts1, meanSignal1, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 3); hold on;
plot(ts2, meanSignal2, 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); 

% Plot vertical line at epoch onset, time = 0
line([0 0], [(min(F465(:) - dcSignal2)), ((max(F465(:)) - dcSignal2))], 'Color', [.7 .7 .7], 'LineStyle','-', 'LineWidth', 3)

% Create the standard error bands for the 405 signal
XX = [ts1, fliplr(ts1)];
YY = [meanSignal1 + stdSignal1, fliplr(meanSignal1 - stdSignal1)];

% Plot filled standard error bands.
h = fill(XX, YY, 'g');
set(h, 'facealpha',.25,'edgecolor','none')

% Repeat for 465
XX = [ts2, fliplr(ts2)];
YY = [meanSignal2 + stdSignal2, fliplr(meanSignal2 - stdSignal2)];
h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('mV', 'FontSize', 12)
title(sprintf('Nosepoke',numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
set(gcf, 'Position',[100, 100, 800, 500])

% Heat Map based on z score of 405 fit subtracted 465
% Fitting 405 channel onto 465 channel to detrend signal bleaching
% Scale and fit data
% Algorithm sourced from Tom Davidson's Github:
% https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m

bls = polyfit(F465(1:end), F405(1:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_dF_all = F465 - Y_fit_all;

zall = zeros(size(Y_dF_all));
for i = 1:size(Y_dF_all,1)
    ind = ts2(1,:) < BASELINE_PER(2) & ts2(1,:) > BASELINE_PER(1);
    zb = mean(Y_dF_all(i,ind)); % baseline period mean (-10sec to -6sec)
    zsd = std(Y_dF_all(i,ind)); % baseline period stdev
    zall(i,:)=(Y_dF_all(i,:) - zb)/zsd; % Z score per bin
end

% Standard error of the z-score
zerror = std(zall)/sqrt(size(zall,1));

% Plot heat map
subplot(3,1,2);
imagesc(ts2, 1, zall);
colormap('jet'); % c1 = colorbar; 
title(sprintf('Z-Score Heat Map', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14);
ylabel('Trials', 'FontSize', 12);

% Fill band values for second subplot. Doing here to scale onset bar
% correctly
XX = [ts2, fliplr(ts2)];
YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];

subplot(3,1,3)
plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)

h = fill(XX, YY, 'r');
set(h, 'facealpha',.25,'edgecolor','none')

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('Z-score', 'FontSize', 12)
title(sprintf('470 nm Z-Score: Nosepoke', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
%c2 = colorbar;
%%
figure(2)
plot(ts2, zall)
line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('Z-score', 'FontSize', 12)
title(sprintf('470 nm Z-Score: Nosepoke', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 14)
% findpeaks(zall, ts2, 'MinPeakProminence',20)
% [events, eventstime] = findpeaks(zall, ts2, 'MinPeakProminence',20);

%%
figure(3)
subplot(2,3,[1,2,4,5])
plot(ts2, mean(zall), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth', 3); hold on;
line([0 0], [min(YY*1.5), max(YY*1.5)], 'Color', [.7 .7 .7], 'LineWidth', 2)

h = fill(XX, YY, 'b');
set(h, 'facealpha',.25,'edgecolor','none')  

% Finish up the plot
axis tight
xlabel('Time, s','FontSize',18)
ylabel('Z-score +/- SEM', 'FontSize', 18)
title(sprintf('GrabDA NAc DA4 Ext', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 18)
box off

subplot(2,3,6);
imagesc(ts2, 1, zall);
colormap('jet'); colorbar; 
title(sprintf('Z-Score/Trial', numel(data.streams.(STREAM_STORE1).filtered), numArtifacts),'FontSize', 16);
xlabel('Time, s', 'FontSize', 12);
ylabel('Trial', 'FontSize', 12);

% Fill band values for second subplot. Doing here to scale onset bar
% correctly
XX = [ts2, fliplr(ts2)];
YY = [mean(zall)-zerror, fliplr(mean(zall)+zerror)];
Zmean = mean(zall);


%% Find max values, AUC and N� of spikes
 
% Find max value of Z-scores
Mepoch = max(zall,[],2);
M = mean(Mepoch)
medianM = median(Mepoch)

% Find AUC in TRANGE of Z-scores
AUC(1,1)=trapz(mean(zall(:,ts2(1,:) < -1 & ts2(1,:) > -5)));
AUC(1,2)=trapz(mean(zall(:,ts2(1,:) > 0 & ts2(1,:) < 2)));
AUC(1,3)=trapz(mean(zall(:,ts2(1,:) > 2 & ts2(1,:) < 30)))
% AUC(1,4)=trapz(mean(zall(:,ts2(1,:) > 15 & ts2(1,:) < 20)))
% AUC(1,5)=trapz(mean(zall(:,ts2(1,:) > 20 & ts2(1,:) < 25)))
% AUC(1,6)=trapz(mean(zall(:,ts2(1,:) > 25 & ts2(1,:) < 30)))

% Find number of spikes arround each epoch
for M1 = 1:size(zall,1)
    [pks,loc] = findpeaks(zall(M1,:),'MinPeakProminence',5);
    P{M1} = [pks; loc];
end
for M1 = 1:size(P,2)
    NrPks(M1) = size(P{M1},2);
end
AvgPeaks = mean(NrPks)

%% Create a fitted curve of the averaged TTL-limited signal 

[xData, yData] = prepareCurveData( ts2, Zmean );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
opts.SmoothingParam = 0.997;
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'Zmean', 'Fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
xlabel( 'Time s', 'Interpreter', 'none' );
ylabel( 'F/F0', 'Interpreter', 'none' );
grid on

% Prepare arrays to export raw data
coef=coeffvalues(fitresult);
fittedZ = coef.coefs(:,4);

