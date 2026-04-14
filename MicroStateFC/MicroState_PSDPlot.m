% PSD Analysis for Identified Microstates (Subject 013AR)
% This script reads the saved microstate file and computes PSD per state

clc; clear all; close all;

%% 1. Configuration and Paths
subjectName = '013AR';
expDate = '280122'; 
protocolName = 'G2';
rootPath = '/Volumes/Suvam''s 1TB/Supratim_Ray_Files/Supratim Ray''s files - segmentedData';

% Path to the saved microstate results
folderSegment = fullfile(rootPath, subjectName, 'EEG', expDate, protocolName, 'segmentedData', 'LFP');
microstateFile = fullfile(folderSegment, [subjectName '_' protocolName '_Microstates.mat']);

if ~exist(microstateFile, 'file')
    error('Microstate result file not found at: %s', microstateFile);
end

%% 2. Load Microstate Labels and Metadata
load(microstateFile); 
% Loaded: microstateMaps, clusterIdx, finalLocs, numElectrodes, numTrials, numSamples

totalSamples = numTrials * numSamples;
Fs = 1000; % Assuming 1kHz sampling (2500 samples / 2.5s)

%% 3. Reconstruct Continuous EEG Data
% Initialize 3D matrix for all 70 electrodes
eegDataAll = zeros(numElectrodes, numTrials, numSamples);

fprintf('Reloading EEG data for %d electrodes...\n', numElectrodes);
lfpInfo = load(fullfile(folderSegment, 'lfpInfo.mat'));
electrodes = lfpInfo.analogChannelsStored;

for i = 1:numElectrodes
    elecFileName = fullfile(folderSegment, ['elec' num2str(electrodes(i)) '.mat']);
    dataVars = load(elecFileName);
    % analogData is 120x2500
    eegDataAll(i, :, :) = dataVars.analogData; 
end

% Reshape to continuous stream (70 x 300,000) to match finalLocs
eegContinuous = reshape(eegDataAll, numElectrodes, totalSamples);

%% 4. Compute PSD per Microstate
numClusters = max(clusterIdx); 
figure('Color', 'w', 'Name', ['Microstate PSD: ' subjectName]);
colors = lines(numClusters);

% We analyze a window around each peak to get sufficient data for pwelch
windowSize = 200; % 200ms window around each peak

for c = 1:numClusters
    % Find indices where this specific microstate (1-5) was active
    stateIndices = finalLocs(clusterIdx == c);
    
    if isempty(stateIndices)
        continue;
    end
    
    % Collect EEG segments for this state across all 70 electrodes
    allSegments = [];
    for s = 1:length(stateIndices)
        idx = stateIndices(s);
        % Ensure window is within signal bounds
        if idx > windowSize/2 && idx < (totalSamples - windowSize/2)
            % Use the spatial average (Global signal) for the PSD
            winData = eegContinuous(:, idx - windowSize/2 : idx + windowSize/2);
            allSegments = [allSegments, mean(winData, 1)]; 
        end
    end
    
    % Compute PSD using Welch's method
    [pxx, f] = pwelch(allSegments, hanning(windowSize), [], windowSize, Fs);
    
    % Plot on a log scale (dB)
    plot(f, 10*log10(pxx), 'Color', colors(c,:), 'LineWidth', 1.5); hold on;
end

%% 5. Formatting the Plot
title(['Power Spectral Density per Microstate (' subjectName ')']);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legendStrings = arrayfun(@(x) sprintf('State %s', char(64+x)), 1:numClusters, 'UniformOutput', false);
legend(legendStrings);
xlim([1 100]); % Focus on Delta through Gamma ranges
grid on;

fprintf('PSD plotting complete for %d microstates.\n', numClusters);