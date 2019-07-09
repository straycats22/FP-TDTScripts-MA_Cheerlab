%% Fiber Photometry Epoch Averaging Example
%
% <html>
% This example goes through fiber photometry analysis using techniques <br>
% such as data smoothing, bleach detrending, and z-score analysis. <br>
% The epoch averaging was done using TDTfilter. <br><br>
% Author Contributions: <br>
% TDT, David Root, and the Morales Lab contributed to the writing and/or conceptualization of the code. <br>
% The signal processing pipeline was inspired by the workflow developed by <a href="https://doi.org/10.1016/j.celrep.2017.10.066">David Barker et al. (2017)</a> for the Morales Lab. <br>
% The data used in the example were provided by David Root. <br><br>
% Author Information: <br>
% David H. Root <br>
% Assistant Professor <br>
% Department of Psychology & Neuroscience <br>
% University of Colorado, Boulder <br>
% Lab Website: <a href="https://www.root-lab.org">https://www.root-lab.org</a> <br>
% david.root@colorado.edu <br><br>
% About the authors: <br>
% The Root lab and Morales lab investigate the neurobiology of reward, aversion, addiction, and depression. <br>
% <br> TDT edits all user submissions in coordination with the contributing
% author(s) prior to publishing.
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add data directories to Matlab
% path.
close all; clear all; clc;
%[MAINPATH,name,ext] = fileparts(cd); % C:\TDT\Synapse\Tanks\
%DATAPATH = fullfile(MAINPATH, 'DanData'); % C:\TDT\Synapse\Tanks\DanData
%[DanPATH,name,ext] = fileparts(MAINPATH); % \TDTMatlabSDK
%addpath(genpath(DanPATH));

%% Importing the Data
% This example assumes you downloaded our example data sets
% (<http://www.tdt.com/support/examples/TDTExampleData.zip link>) and extracted
% it into the \TDTMatlabSDK\Examples\ directory. To import your own data, replace
% |BLOCKPATH| with the path to your own data block.
%
% In Synapse, you can find the block path in the database. Go to Menu > History. 
% Find your block, then Right-Click > Copy path to clipboard.

addpath('C:\Users\Tibo\Documents\MATLAB\Photometry\Data\GrabDA Sucrose Pellet\DA3_FR5-1');
addpath('C:\TDT\TDTbin2mat');

%% Setup the variables for the data you want to extract
% We will extract two different stream stores surrounding the 'PtAB' epoch 
% event. We are interested in a specific event code for the shock onset.

REF_EPOC = 'Gate'; % Stimulation event to center on
STIM_EPOC = 'Stim'; %Stimulation pulses
STREAM_STORE1 = 'x470A' % name of the 470 store
STREAM_STORE2 = 'x405A' % name of the 405 store
TRANGE = [-5 25]; %window size [start time relative to epoc onset, entire duration]
ARANGE = [1 1];
BASELINE_PER = [-5 -1]; % baseline period before stim
ARTIFACT = Inf; % variable created for artifact removal

% Now read the specified data from our block into a Matlab structure.
data = TDTbin2mat('C:\Users\AF239\OneDrive\Desktop\Ach Photo Data\Raw Data\DanTest_Experiment-190108\AFTest8', 'TYPE', {'epocs', 'scalars', 'streams'});
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
art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE1).filtered, 'UniformOutput',false));
good = ~art1 & ~art2;
data.streams.(STREAM_STORE1).filtered = data.streams.(STREAM_STORE1).filtered(good);

art1 = ~cellfun('isempty', cellfun(@(x) x(x>ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
art2 = ~cellfun('isempty', cellfun(@(x) x(x<-ARTIFACT), data.streams.(STREAM_STORE2).filtered, 'UniformOutput',false));
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

allSignals1 = cell2mat(data.streams.(STREAM_STORE1).filtered');
allSignals2 = cell2mat(data.streams.(STREAM_STORE2).filtered');

%NEED TO CORRECT BELOW
% Delete row if 1st value is less than 300 or 100
allSignals1(allSignals1(:, 1)<= 300, :)= [];
allSignals2(allSignals2(:, 1)<= 100, :)= [];

N = length(art1);
t = minLength1;

% Create the time vector for each stream store
ts1 = TRANGE(1) + (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
ts2 = TRANGE(1) + (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;

%average 470 signal
F470 = zeros(size(allSignals1(:,1:N:end-N+1)));
for ii = 1:size(allSignals1,1)
    F470(ii,:) = arrayfun(@(i) mean(allSignals1(ii,i:i+N-1)),1:N:length(allSignals1)-N+1);
end
minLength1 = size(F470,2);

%average 405 signal
F405 = zeros(size(allSignals2(:,1:N:end-N+1)));
for ii = 1:size(allSignals2,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals2(ii,i:i+N-1)),1:N:length(allSignals2)-N+1);
end
minLength1 = size(F405,2);

% Create mean signal and DC offset of 470 signal
meanSignal1 = mean(F470);
dcSignal1 = mean(meanSignal1);

% Create mean signal and DC offset of 405 signal
meanSignal2 = mean(F405);
dcSignal2 = mean(meanSignal2);

% Subtract DC offset to get signals on top of one another
meanSignal1 = meanSignal1 - dcSignal1;
meanSignal2 = meanSignal2 - dcSignal2;

% Plot the 470 and 405 signals


a = 1:minLength1;
b = 1:minLength2;

%find value at stim onset
c1 = allSignals1(1:ii,5085);
c2 = allSignals2(1:ii,5085);

%fine mean of each row
c3 = allSignals1;
meanrow_c3 = mean(c3')';
c4 = allSignals2;
meanrow_c4 = mean(c4')';

%create matrices to subtract value at stim onset for 470 signal
m1 = allSignals1(1:ii,1:t);
p1 = repelem (c1,1,[t]);

%create matrices to subtract value at stim onset for 405 signal
m2 = allSignals2(1:ii,1:t);
p2 = repelem (c2,1,[t]);

norm1 = m1 - p1;
norm2 = m2 - p2;

%create matrices to subtract value at stim onset for 470 signal
p3 =  repelem (meanrow_c3,1,[t]);
p4 =  repelem (meanrow_c4,1,[t]);

dff1 = 100*((m1 - p3)./p3);
dff2 = 100*((m2 - p4)./p4);


figure;
subplot(2,1,1);
plot(ts1,norm1); 
%axis tight
xlabel('Time, s','FontSize',12)
ylabel('df', 'FontSize', 12)
title('470 nm')
legend('60', '30', '10', '50', '50', '40', '60')
%set(gcf, 'Position',[100, 100, 800, 500]); hold on;

subplot(2,1,2);
plot(ts2,norm2);
%axis tight
xlabel('Time, s','FontSize',12)
ylabel('df/f', 'FontSize', 12)
title('405nm')
%set(gcf, 'Position',[100, 100, 800, 500])

figure;
subplot(2,1,1);
plot(ts1,dff1); 
%axis tight
xlabel('Time, s','FontSize',12)
ylabel('df/f', 'FontSize', 12)
title('470 nm')
legend('60', '30', '10', '50', '50', '40', '60')
%set(gcf, 'Position',[100, 100, 800, 500]); hold on;

subplot(2,1,2);
plot(ts2,dff2);
%axis tight
xlabel('Time, s','FontSize',12)
ylabel('df/f', 'FontSize', 12)
title('405nm')
%set(gcf, 'Position',[100, 100, 800, 500])
