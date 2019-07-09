%% Housekeeping
% Clear workspace and close existing figures. Add data directories to Matlab
% path.
close all; clear all; clc;
data = TDTbin2mat('C:\Users\Tibo\Documents\MATLAB\Photometry\Data\GrabDA Sucrose Pellet\DA3_FR5-1', 'TYPE', {'epocs', 'scalars', 'streams'});
%[MAINPATH,name,ext] = fileparts(cd); % C:\TDT\Synapse\Tanks\
%DATAPATH = fullfile(MAINPATH, 'DanData'); % C:\TDT\Synapse\Tanks\DanData
%[DanPATH,name,ext] = fileparts(MAINPATH); % \TDTMatlabSDK
%addpath(genpath(DanPATH));


%%
% Set up the varibles for the data you want to extract. We will extract
% channel 3 from the eNe1 snippet data store, created by the PCA Sorting
% gizmo, and use our PulseGen epoc event ('PC0/') as our stimulus onset.
REF_EPOC = 'Rew/';
STREAM_STORE = 'x470A';
CHANNEL = 3;
ARTIFACT = inf; % optionally set an artifact rejection level
TRANGE = [-5, 10]; % window size [start time relative to epoc onset, window duration]

%% Use TDTfilter to extract data around our epoc event.
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event.  For stream events, the chunks of data are stored in cell
% arrays.
data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);

%%
% Optionally remove artifacts.
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE).filtered = data.streams.(STREAM_STORE).filtered(good);
numArtifacts = sum(~good);

%%
% Applying a time filter to a uniformly sampled signal means that the
% length of each segment could vary by one sample.  Let's find the minimum
% length so we can trim the excess off before calculating the mean.
minLength = min(cellfun('prodofsize', data.streams.(STREAM_STORE).filtered));
data.streams.(STREAM_STORE).filtered = cellfun(@(x) x(1:minLength), data.streams.(STREAM_STORE).filtered, 'UniformOutput',false);

%%
% Find the average signal.
allSignals = cell2mat(data.streams.(STREAM_STORE).filtered');
meanSignal = mean(allSignals);
stdSignal = std(allSignals);

%% Ready to plot
% Create the time vector.
ts = TRANGE(1) + (1:minLength) / data.streams.(STREAM_STORE).fs;

%%
% Plot all the signals as gray.
plot(ts, allSignals','Color', [.85 .85 .85]); hold on;

%%
% Plot vertical line at time=0.
line([0 0], [min(allSignals(:)), max(allSignals(:))], 'Color', 'r', 'LineStyle','-', 'LineWidth', 3)

%%
% Plot the average signal.
plot(ts, meanSignal, 'b', 'LineWidth', 3)

%%
% Plot the standard error of the mean.
plot(ts, (meanSignal)+(stdSignal/sqrt(size(allSignals,1))), 'b--', ts, (meanSignal)-(stdSignal/sqrt(size(allSignals,1))), 'b--');

%%
% Finish up the plot
axis tight
xlabel('Time, s','FontSize',12)
ylabel('mV', 'FontSize', 12)
title(sprintf('%s %d Trials (%d Artifacts Removed)', STREAM_STORE, numel(data.streams.(STREAM_STORE).filtered), numArtifacts))
set(gcf, 'Position',[100, 100, 800, 500])