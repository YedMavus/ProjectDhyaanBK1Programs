% GFP and Microstate Analysis for Subject 013AR (G2 Phase)
% Automatically detects electrodes using lfpInfo.mat
clc;
clear all;
close all;
%% 1. Path Configuration
subjectName = '013AR';
expDate = '280122';
protocolName = 'G2';
rootPath = '/Volumes/Suvam''s 1TB/Supratim_Ray_Files/Supratim Ray''s files - segmentedData';

% Path directly to the LFP folder where lfpInfo.mat and elec files reside
folderSegment = fullfile(rootPath, subjectName, 'EEG', expDate, protocolName, 'segmentedData', 'LFP');

%% 2. Load Metadata (lfpInfo.mat)
lfpInfoFile = fullfile(folderSegment, 'lfpInfo.mat');
if ~exist(lfpInfoFile, 'file')
    error('Could not find lfpInfo.mat at: %s', lfpInfoFile);
end

lfpInfo = load(lfpInfoFile);
% analogChannelsStored contains the list of all available electrodes
electrodes = lfpInfo.analogChannelsStored; 

numElectrodes = length(electrodes);
numTrials = 120;   % Total trials in G1/G2
numSamples = 2500; % 1.25s gray + 1.25s grating

% Initialize 3D Matrix: [Electrodes x Trials x Samples]
eegDataAll = zeros(numElectrodes, numTrials, numSamples);

%% 3. Load All Electrode Data
fprintf('Loading %d electrodes for %s...\n', numElectrodes, subjectName);
for i = 1:numElectrodes
    elecFileName = fullfile(folderSegment, ['elec' num2str(electrodes(i)) '.mat']);
    if exist(elecFileName, 'file')
        dataVars = load(elecFileName);
        % analogData is a 120x2500 matrix
        eegDataAll(i, :, :) = dataVars.analogData; 
    else
        warning('Missing file: %s', elecFileName);
    end
end

%% 4. Compute Global Field Potential (GFP) - Continuous
% Reshape from [70 x 120 x 2500] to [70 x 300000]
% This concatenates all 120 trials into one long signal
totalSamples = numTrials * numSamples;
eegContinuous = reshape(eegDataAll, numElectrodes, totalSamples);

% Compute GFP as the spatial standard deviation at every time point
% This preserves the sharp peaks of individual trials
gfp = std(eegContinuous, 1, 1); 

% Create a continuous time axis for the overall duration
% Total duration: 120 trials * 2.5s = 300 seconds
totalDuration = numTrials * (numSamples / (numSamples/2.5)); % 300s
timeAxis = linspace(0, totalDuration, totalSamples);
%% 5. Peak Identification (Handling Float Outputs)
% 1. Find all peaks
xmax_all = findpeaks(gfp'); 

if isstruct(xmax_all)
    allLocs = [xmax_all.loc]; 
    allPks = gfp(allLocs);
elseif ismatrix(xmax_all) && any(mod(xmax_all, 1) ~= 0)
    % If findpeaks returned amplitudes (floats), find their indices in gfp
    allPks = xmax_all;
    allLocs = zeros(size(allPks));
    for p = 1:length(allPks)
        [~, idx] = min(abs(gfp - allPks(p))); % Find closest matching index
        allLocs(p) = idx;
    end
else
    % If findpeaks returned integer indices
    allLocs = xmax_all;
    allPks = gfp(allLocs);
end

% 2. Calculate rejection threshold
amplitudeThreshold = prctile(allPks, 15);

% 3. Extract final filtered peaks
xmax_filtered = findpeaks(gfp', amplitudeThreshold);

if isstruct(xmax_filtered)
    finalLocs = [xmax_filtered.loc];
elseif ismatrix(xmax_filtered) && any(mod(xmax_filtered, 1) ~= 0)
    % Convert filtered float amplitudes to indices
    finalPks_tmp = xmax_filtered;
    finalLocs = zeros(size(finalPks_tmp));
    for p = 1:length(finalPks_tmp)
        [~, idx] = min(abs(gfp - finalPks_tmp(p)));
        finalLocs(p) = idx;
    end
else
    finalLocs = xmax_filtered;
end

finalPeaks = gfp(finalLocs);

fprintf('Detected %d peaks; kept %d after 15%% rejection.\n', length(allLocs), length(finalLocs));
%% 6. Microstate Clustering 
mapsAtPeaks = eegContinuous(:, finalLocs); 

% Perform K-means clustering on the extracted peaks
numClusters = 5;
fprintf('Clustering %d peaks into %d states...\n', length(finalLocs), numClusters);

[clusterIdx, microstateMaps] = kmeans(mapsAtPeaks', numClusters, ...
    'Distance', 'correlation', 'Replicates', 15, 'MaxIter',500);

%% 7. Visualization (Updated for 5 States)
figure('Color', 'w', 'Name', ['Microstate Analysis: ' subjectName]);

% Plot GFP Profile across the top row
subplot(2, numClusters, 1:numClusters);
plot(timeAxis, gfp, 'k', 'LineWidth', 1.2); hold on;
plot(timeAxis(finalLocs), finalPeaks, 'ro', 'MarkerSize', 2);
title(['GFP Peaks for Subject ' subjectName ' (Total Duration)']);
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
grid on;

% Plot Topographic Maps for the 5 Clusters in the bottom row
for c = 1:numClusters
    subplot(2, numClusters, numClusters + c);
    % Displays the electrode vector as a column
    imagesc(reshape(microstateMaps(c, :), [], 1)); 
    title(['State ' char(64+c)]);
    set(gca, 'XTick', [], 'YTick', []);
    colorbar;
end
%% 8. Save Microstate Results for Downstream PSD
% Define the save filename
saveFileName = fullfile(folderSegment, [subjectName '_' protocolName '_Microstates.mat']);

% Variables to save:
% microstateMaps: The 4 cluster centers (70 x 4)
% clusterIdx: The label (1-4) for every identified peak
% finalLocs: The time indices (samples) where those peaks occurred
% numElectrodes, numSamples, numTrials: Metadata for reconstruction

fprintf('Saving microstate results to: %s\n', saveFileName);
save(saveFileName, 'microstateMaps', 'clusterIdx', 'finalLocs', ...
    'subjectName', 'protocolName', 'numElectrodes', 'numTrials', 'numSamples');
