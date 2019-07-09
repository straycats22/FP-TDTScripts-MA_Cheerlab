%% Stream Plot Example
%
% <html>
% Import continuous data into Matlab using TDTbin2mat <br>
% Plot a single channel of data with various filtering schemes <br>
% Good for first-pass visualization of streamed data <br>
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clear all; clc;
%[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
DATAPATH = ('D:\PhD\CheerLab\WIN78-FR1-7_WIN79-FR1-7'); % \TDTMatlabSDK\Examples\ExampleData
%[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
%addpath(genpath(SDKPATH));

%% Importing the Data
% This example assumes you downloaded our example data sets
% (<http://www.tdt.com/support/examples/TDTExampleData.zip link>) and extracted
% it into the \TDTMatlabSDK\Examples\ directory. To import your own data, replace
% |BLOCKPATH| with the path to your own data block.
%
% In Synapse, you can find the block path in the database. Go to Menu > History. 
% Find your block, then Right-Click > Copy path to clipboard.
BLOCKPATH = fullfile(DATAPATH);
addpath('C:\Users\migue\OneDrive - upf.edu\Documentos\MATLAB\TDTMatlabSDK\');
%%
% Now read channel 1 from all stream data into a Matlab structure called 'data'.
data = TDTbin2mat('D:\PhD\CheerLab\WIN78-FR1-7_WIN79-FR1-7', 'TYPE', {'epocs', 'scalars', 'streams'});

%% 
% And that's it! Your data is now in Matlab. The rest of the code is a
% simple plotting example.

%% Stream Store Plotting
% Let's create time vectors for each stream store for plotting in time.

STREAM_STORE2 = 'x470A'
STREAM_STORE1 = 'x405A'
ARTIFACT = Inf; % variable created for artifact removal
%time_LFP1 = (1:length(data.streams.LFP1.data))/data.streams.LFP1.fs;
%time_pNe1 = (1:length(data.streams.pNe1.data))/data.streams.pNe1.fs;


%allSignals1 = TDTfilter(data.streams.(STREAM_STORE1).filtered);
%allSignals2 = TDTfilter(data.streams.(STREAM_STORE2).filtered);

allSignals1 = data.streams.x470A.data;
allSignals2 = data.streams.x405A.data;

time_x470D = (1:length(data.streams.x470A.data))/data.streams.x470A.fs;


figure;
subplot(2,1,1)
plot(time_x470D, allSignals1,  'color',[0, 0.5, 0], 'LineWidth', 1.5) 
axis tight
title({'Unexp Rews','470'},'FontSize',14)
ylabel('mV')
% xlim([100 inf])

hold on;
subplot(2,1,2)
plot(time_x470D, allSignals2, 'color',[0.4940, 0.1840, 0.5560], 'LineWidth', 1.5) 
axis tight
title({'405'},'FontSize',14)
ylabel('mV')
xlabel('time (s)')
% xlim([100 inf]);

%%
% Optionally remove artifacts. If any waveform is above ARTIFACT level, or
% below -ARTIFACT level, remove it from the data set.

t = 40; % time threshold below which we will discard
ind = find(time_x470D>t,1); % find first index of when time crosses threshold
time_x470D = time_x470D(ind:end); % reformat vector to only include allowed time
data.streams.x470A.data = data.streams.x470A.data(ind:end);
data.streams.x405A.data = data.streams.x405A.data(ind:end);

allSignals1 = data.streams.x470A.data;
allSignals2 = data.streams.x405A.data;

% downsample 10x and average 405 signal
N = 10;
F405 = zeros(size(allSignals2(:,1:N:end-N+1)));
for ii = 1:size(allSignals2,1)
    F405(ii,:) = arrayfun(@(i) mean(allSignals2(ii,i:i+N-1)),1:N:length(allSignals2)-N+1);
end
minLength1 = size(F405,2);

% Create mean signal, standard error of signal, and DC offset of 405 signal
%meanSignal1 = mean(F405);
%stdSignal2 = std(double(F405))/sqrt(size(F405,1));
%dcSignal2 = mean(allSignals2);

% downsample 10x and average 465 signal
%allSignals = cell2mat(data.streams.(STREAM_STORE2).filtered');
F465 = zeros(size(allSignals1(:,1:N:end-N+1)));
for ii = 1:size(allSignals1,1)
    F465(ii,:) = arrayfun(@(i) mean(allSignals1(ii,i:i+N-1)),1:N:length(allSignals1)-N+1);
end
minLength2 = size(F465,2);

% Fitting 405 channel onto 465 channel to detrend signal bleaching
% Scale and fit data
% Algorithm sourced from Tom Davidson's Github:
% https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m

bls = polyfit(F465(2:end), F405(2:end), 1);
Y_fit_all = bls(1) .* F405 + bls(2);
Y_dF_all = F465 - Y_fit_all;

dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));

%zall = zeros(size(Y_dF_all));

detrend_F465 = detrend(dFF);
% F465detrend = F465 - detrend_F465;
% 
% B = F465(:,[500:end]);
% C = F405(:, [500:end]);
% D = Y_dF_all(:, [500:end]);
% E = detrend_F465(:, [500:end]);
% meanSignal2 = mean(F465);
% stdSignal1 = std(double(F465))/sqrt(size(F465,1));
% dcSignal1 = mean(allSignals1);
% dcSignal1 = mean(B);
% dcSignal2 = mean(C);
% dcSignal3 = mean(D);
% 
% dcSignal4 = mean(E);
% delta1 = F465 - dcSignal1;
% delta2 = (F405/dcSignal2)*100;
% delta3 = Y_dF_all - dcSignal3;
% delta4 = F465detrend - dcSignal4;
% delta5 = (detrend_F465/dcSignal4)*100;
% delta6 = delta5 - delta2;



% figure;
% ax1 = subplot(2,1,1)
% plot(time_x470D, data.streams.x470A.data(1,:)*1e6, 'color',[0, 0.5, 0], 'LineWidth', 1.5)
% axis tight
% title({'Unexp Rews','dLight'},'FontSize',14)
% xlim([100 inf]);
% 
% ymin = min(delta1(200:end));
% ymax = max(delta1(200:end)); hold on;


% Create the time vector for each stream store
 ts1 = (1:minLength1) / data.streams.(STREAM_STORE1).fs*N;
 ts2 = (1:minLength2) / data.streams.(STREAM_STORE2).fs*N;

% Create mean signal, standard error of signal, and DC offset of 465 signal


%ax2 = subplot(3,1,2);
%plot(time_LFP1, data.streams.LFP1.data(1,:)*1e6,'r'); 
%title('LFP Waveform','FontSize',14);

%ax3 = subplot(3,1,3);
%plot(time_pNe1, data.streams.pNe1.data(1,:),'k'); 
%title('Plot Decimated Spikes','FontSize',14);
y1 = min(detrend_F465);
y2 = max(detrend_F465);



%% 
% Plot the 405 and 465 average signals
figure;
% subplot(3,1,1)
% plot(ts1, delta6,  'Color',[0, 0.5, 0], 'LineWidth', 1.5) %[0.4660, 0.6740, 0.1880]
% axis tight
% title({'Unexp Rews','470'},'FontSize',14)
% ylabel('df/f')
% xlim([10 inf])
% ylim([y1 y2])
% hold on;
% subplot(2,1,1)
% plot(ts1, delta2, 'Color',[0.4940, 0.1840, 0.5560], 'LineWidth', 1.5) %[0.4660, 0.6740, 0.1880]
% axis tight
% title({'405'},'FontSize',14)
% ylabel('df/f')
% xlabel('time (s)')
% %ylim([ymin 15])
% line([242 242],[-9 6], 'Color', [.7 .7 .7], 'LineWidth', 2) %[min(y1), max(y2)],
% xlim([10 inf]);
% 
% subplot(2,1,2)
plot(ts1, detrend_F465, 'Color',[0, 0.5, 0], 'LineWidth', 1.5), hold on  %[0.4660, 0.6740, 0.1880]
axis tight
%findpeaks(detrend_F465, ts2, 'MinPeakProminence',20)
%[events, eventstime] = findpeaks(detrend_F465, ts2, 'MinPeakProminence',20);
title({'470 de-trend'},'FontSize',14)
ylabel('df/f')
xlabel('time (s)')
% line([242 242], [95 110],  'Color', [.7 .7 .7], 'LineWidth', 2) % [min(y1), max(y2)],
% xlim([10 inf]);



%figure;
%subplot(4,1,1)
%ax = [ax1];
%axis([ax1], 'tight')
%linkaxes(ax, 'x')
%xlim(ax(end), [10 inf]);
%xlabel(ax(end), 'Time (s)','FontSize',12)
%ylabel(ax(end), 'Amplitude (\muV)','FontSize',12);


% Enlarge figure.
set(gcf, 'Units', 'centimeters', 'OuterPosition', [10, 10, 20, 20]);

%% Epoc Events
% Generate continuous time series for epoc data using epoc timestamps

% StimSync epoc event
STIMSYNC = 'aCu_';
aCu_off = data.epocs.aCu_.offset;
aCu_on = data.epocs.aCu_.onset;
aCu_x = reshape(kron([aCu_on, aCu_off], [1, 1])', [], 1);

%%
% Make a time series waveform of epoc values and plot them.
sz = length(aCu_on);
d = data.epocs.aCu_.data';
aCu_y = reshape([zeros(1, sz); d; d; zeros(1, sz)], 1, []);
hold on; plot(aCu_x, 7*(aCu_y) - 3, 'g', 'LineWidth', 2);