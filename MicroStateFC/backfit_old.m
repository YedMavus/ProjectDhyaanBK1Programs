clc; clear; close all;

%% ===== LOAD MICROSTATE MODEL =====
load('MicroStateFC/microstate_model.mat');
templates = microstate.templates;
K = microstate.K;

%% ===== CONFIG =====
rootPath = '/Users/suvamdey/Desktop/NSP/data/segmentedData';
subjects = {'003S','013AR'};
stages = {'M1','G1','M2','G2'};

Fs = 1000;
minSamples = round(0.03 * Fs); % 30 ms

%% ===== LOAD MONTAGE =====
% Ensure the Functions folder is in your path to access createMontageEEG
addpath(genpath('/Users/suvamdey/Desktop/NSP/Montages')); 

capType = 'actiCap64_UOL';

% 1. Generate the montage files if they don't exist
% This creates 'actiCap64_UOL.mat' in Montages/Layouts/actiCap64_UOL/
% createMontageEEG(capType); 

% 2. Load the generated chanlocs file
% The variable inside the file is named 'chanlocs'
montagePath = fullfile('/Users/suvamdey/Desktop/NSP/Montages/Layouts', capType, [capType '.mat']);
tmpMontage = load(montagePath);
chanlocs = tmpMontage.chanlocs;
%% ===== NORMALIZE + FIX POLARITY =====
for k = 1:K
    templates(:,k) = templates(:,k) / norm(templates(:,k));
    if mean(templates(:,k)) < 0
        templates(:,k) = -templates(:,k);
    end
end

%% ===== BACKFITTING =====
results = struct();
idx = 1;

for s = 1:length(subjects)

    subjectName = subjects{s};
    subjPath = fullfile(rootPath, subjectName, 'EEG');

    dates = dir(subjPath);
    dates = dates([dates.isdir] & ~startsWith({dates.name}, '.'));

    for d = 1:length(dates)

        expDate = dates(d).name;

        for st = 1:length(stages)

            stage = stages{st};
            folderSegment = fullfile(rootPath, subjectName, ...
                'EEG', expDate, stage, ...
                'segmentedData','LFP');

            if ~exist(folderSegment,'dir'), continue; end

            %% LOAD DATA
            lfpInfo = load(fullfile(folderSegment,'lfpInfo.mat'));
            electrodes = lfpInfo.analogChannelsStored;
            numElectrodes = length(electrodes);

            tmp = load(fullfile(folderSegment, ['elec' num2str(electrodes(1)) '.mat']));
            [numTrials, numSamples] = size(tmp.analogData);

            eegDataAll = zeros(numElectrodes, numTrials, numSamples);

            for i = 1:numElectrodes
                file = fullfile(folderSegment, ['elec' num2str(electrodes(i)) '.mat']);
                if exist(file,'file')
                    D = load(file);
                    eegDataAll(i,:,:) = D.analogData;
                end
            end

            %% RESHAPE
            data = reshape(eegDataAll, numElectrodes, numTrials*numSamples);

            %% FILTER + CAR
            data = bandpass(data', [0.5 45], Fs)';
            data = data - mean(data,1);

            %% NORMALIZE
            data = data ./ vecnorm(data);

            %% BACKFIT
            T = size(data,2);
            labels = zeros(1,T);

            for t = 1:T
                corrVals = abs(templates' * data(:,t));
                [~, labels(t)] = max(corrVals);
            end

            %% SMOOTH
            labels = smooth_microstates(labels, minSamples);

            %% STATS
            stats = compute_stats(labels, Fs, K);

            %% STORE
            results(idx).subject = subjectName;
            results(idx).stage = stage;
            results(idx).date = expDate;
            results(idx).labels = labels;
            results(idx).stats = stats;

            idx = idx + 1;
        end
    end
end

%% ===== SAVE =====
save('microstate_results.mat','results','templates','-v7.3');

%% ===== PLOT TEMPLATES =====
figure;
for k = 1:K
    subplot(2, ceil(K/2), k);
    topoplot(templates(:,k), chanlocs, 'electrodes','on');
    title(['Microstate ' num2str(k)]);
end

%% ===== PLOT DIFFERENCES =====
figure;
count = 1;
for i = 1:K
    for j = i+1:K
        subplot(K-1,K-1,count);
        topoplot(templates(:,i)-templates(:,j), chanlocs);
        title(['M' num2str(i) '-' num2str(j)]);
        count = count + 1;
    end
end

%% ===== PLOT EXAMPLE SEQUENCE =====
figure;
plot(results(1).labels);
ylim([1 K]);
xlabel('Time'); ylabel('State');
title('Microstate Sequence');

%% Define required function
function labels = smooth_microstates(labels, minSamples)
% smooth_microstates: Removes segments shorter than minSamples
% and replaces them with the surrounding state.

N = length(labels);
i = 1;
while i <= N
    % Find the length of the current segment
    start_idx = i;
    current_val = labels(i);
    while i <= N && labels(i) == current_val
        i = i + 1;
    end
    end_idx = i - 1;
    segment_len = end_idx - start_idx + 1;

    % If segment is too short, replace it
    if segment_len < minSamples
        if start_idx > 1
            % Replace with the state immediately preceding the segment
            labels(start_idx:end_idx) = labels(start_idx - 1);
        elseif end_idx < N
            % If it's the very first segment, replace with the following state
            labels(start_idx:end_idx) = labels(end_idx + 1);
        end
    end
end
end

function stats = compute_stats(labels, Fs, K)
T = length(labels);
stats = struct();
for k = 1:K
    idx = (labels == k);
    % Coverage: % of total time spent in this state
    stats(k).coverage = sum(idx) / T;

    % Occurrence: Number of times state appeared per second
    starts = diff([0, idx]) == 1;
    num_entries = sum(starts);
    stats(k).occurrence = num_entries / (T / Fs);

    % Duration: Average length of a segment in milliseconds
    if num_entries > 0
        stats(k).meanDuration = (sum(idx) / num_entries) * (1000 / Fs);
    else
        stats(k).meanDuration = 0;
    end
end
end