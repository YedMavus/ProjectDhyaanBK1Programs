clc; clear all; close all;

%% 1. Path and Subject Configuration
subjectName = '013AR';
expDate = '280122'; 
protocolName = 'G2';
rootPath = '/Volumes/Suvam''s 1TB/Supratim_Ray_Files/Supratim Ray''s files - segmentedData';
folderSegment = fullfile(rootPath, subjectName, 'EEG', expDate, protocolName, 'segmentedData', 'LFP');

%% 2. Define Frequency Bands
% Define a struct array with the names and frequency ranges
bands(1).name = 'alpha'; bands(1).range = [8 13];
bands(2).name = 'gamma'; bands(2).range = [30 80];

%% 3. Load Metadata and Raw Data
lfpInfo = load(fullfile(folderSegment, 'lfpInfo.mat'));
electrodes = lfpInfo.analogChannelsStored; 
numElectrodes = length(electrodes);
numTrials = 120;   % 360 for M1, 120 for G2
numSamples = 2500; 
Fs = 1000; 

eegDataAll = zeros(numElectrodes, numTrials, numSamples);
fprintf('Loading data for %d electrodes...\n', numElectrodes);
for i = 1:numElectrodes
    elecFile = fullfile(folderSegment, ['elec' num2str(electrodes(i)) '.mat']);
    if exist(elecFile, 'file')
        dataVars = load(elecFile);
        eegDataAll(i, :, :) = dataVars.analogData; 
    end
end

% Flatten to continuous for filtering
totalSamples = numTrials * numSamples;
eegContinuousRaw = reshape(eegDataAll, numElectrodes, totalSamples);

%% 4. Iterate Through Bands
for b = 1:length(bands)
    currentBand = bands(b).name;
    lowFreq = bands(b).range(1);
    highFreq = bands(b).range(2);
    
    fprintf('\n--- Processing %s band (%d-%d Hz) ---\n', currentBand, lowFreq, highFreq);

    %% 4.1 Manual FFT Filtering
    eegContinuous = zeros(size(eegContinuousRaw));
    for i = 1:numElectrodes
        signal = eegContinuousRaw(i,:);
        L = length(signal);
        signal_fft = fft(signal);
        freqs = (0:L-1)*(Fs/L);
        filterMask = (freqs >= lowFreq & freqs <= highFreq) | ...
                     (freqs >= (Fs - highFreq) & freqs <= (Fs - lowFreq));
        eegContinuous(i,:) = real(ifft(signal_fft .* filterMask));
    end

    %% 4.2 Compute GFP
    gfp = std(eegContinuous, 1, 1); 
    totalDuration = numTrials * 2.5;
    timeAxis = linspace(0, totalDuration, totalSamples);

    %% 4.3 Peak Identification
    xmax_all = findpeaks(gfp'); 
    if isstruct(xmax_all)
        allLocs = [xmax_all.loc]; 
    elseif ismatrix(xmax_all) && any(mod(xmax_all, 1) ~= 0)
        allPks_tmp = xmax_all;
        allLocs = zeros(size(allPks_tmp));
        for p = 1:length(allPks_tmp)
            [~, idx] = min(abs(gfp - allPks_tmp(p)));
            allLocs(p) = idx;
        end
    else
        allLocs = xmax_all;
    end
    
    allPks = gfp(allLocs);
    amplitudeThreshold = prctile(allPks, 15);
    xmax_filtered = findpeaks(gfp', amplitudeThreshold);
    
    if isstruct(xmax_filtered)
        finalLocs = [xmax_filtered.loc];
    elseif ismatrix(xmax_filtered) && any(mod(xmax_filtered, 1) ~= 0)
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

    %% 4.4 Clustering
    mapsAtPeaks = eegContinuous(:, finalLocs); 
    numClusters = 5;
    [clusterIdx, microstateMaps] = kmeans(mapsAtPeaks', numClusters, ...
        'Distance', 'correlation', 'Replicates', 15, 'MaxIter', 500);

    %% 4.5 Save Results with Band Suffix
    saveFileName = fullfile(folderSegment, [subjectName '_' protocolName '_Microstates_' currentBand '.mat']);
    save(saveFileName, 'microstateMaps', 'clusterIdx', 'finalLocs', ...
        'subjectName', 'protocolName', 'numElectrodes', 'numTrials', 'numSamples', 'lowFreq', 'highFreq');
    fprintf('Saved: %s\n', saveFileName);
    
    %% 4.6 Optional Visualization per Band
    figure('Color', 'w', 'Name', [upper(currentBand) ' Analysis: ' subjectName]);
    subplot(2, numClusters, 1:numClusters);
    plot(timeAxis, gfp, 'k'); hold on;
    plot(timeAxis(finalLocs), finalPeaks, 'ro', 'MarkerSize', 2);
    title([upper(currentBand) ' Band GFP Peaks']);
    
    for c = 1:numClusters
        subplot(2, numClusters, numClusters + c);
        imagesc(reshape(microstateMaps(c, :), [], 1)); 
        title(['State ' char(64+c)]); set(gca, 'XTick', [], 'YTick', []);
    end
end