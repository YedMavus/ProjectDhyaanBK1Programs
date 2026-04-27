

clc; clear; close all;

%% ===== LOAD MICROSTATE MODEL =====
load('MicroStateFC/microstate_model_gamma.mat');
templates = microstate.templates;
K = microstate.K;

%% ===== CONFIG =====
rootPath = '/Users/suvamdey/Desktop/NSP/data/segmentedData';
subjects = {'003S','013AR'};
stages = {'EO1', 'EC1', 'G1', 'M1','G2','EO2', 'EC2', 'M2'};

Fs = 1000;
minSamples = round(0.03 * Fs); % 30 ms

%% ===== LOAD MONTAGE =====
addpath(genpath('/Users/suvamdey/Desktop/NSP/Montages')); 

capType = 'actiCap64_UOL';
montagePath = fullfile('/Users/suvamdey/Desktop/NSP/Montages/Layouts', capType, [capType '.mat']);
tmpMontage = load(montagePath);
chanlocs = tmpMontage.chanlocs;

%% ===== NORMALIZE + FIX POLARITY =====
templates = templates ./ sqrt(sum(templates.^2,1));

for k = 1:K
    if mean(templates(:,k)) < 0
        templates(:,k) = -templates(:,k);
    end
end

%% ===== PRECOMPUTE FILTER =====
[b,a] = butter(4, [0.5 45]/(Fs/2), 'bandpass');

%% ===== BACKFITTING =====
results = struct();
idx = 1;

for s = 1:length(subjects)
    
    subjectName = subjects{s};
    subjPath = fullfile(rootPath, subjectName, 'EEG');
    
    dates = dir(subjPath);
    dates = dates([dates.isdir] & ~startsWith({dates.name}, '.'));
    
    
        
    expDate = dates.name;
    
    for st = 1:length(stages)
        
        stage = stages{st};
        
        folderSegment = fullfile(rootPath, subjectName, ...
            'EEG', expDate, stage, ...
            'segmentedData','LFP');
        
        if ~exist(folderSegment,'dir'), continue; end
        
        %% ===== LOAD DATA =====
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
        
        %% ===== RESHAPE =====
        data = reshape(eegDataAll, numElectrodes, []);
        
        %% ===== FILTER + CAR =====
        data = filtfilt(b, a, data')';
        data = data - mean(data,1);
        
        %% ===== NORMALIZE =====
        data = data ./ sqrt(sum(data.^2,1));
        
        %% ===== BACKFIT (VECTORIZED) =====
        corrMat = abs(templates' * data);   % K x T
        [~, labels] = max(corrMat, [], 1);
        
        %% ===== SMOOTH =====
        labels = smooth_microstates(labels, minSamples);
        
        %% ===== STATS =====
        stats = compute_stats(labels, Fs, K);
        
        %% ===== STORE =====
        results(idx).subject = subjectName;
        results(idx).stage = stage;
        results(idx).date = expDate;
        results(idx).labels = labels;
        results(idx).stats = stats;
        
        idx = idx + 1;
    end

end

%% ===== SAVE =====
save('MicroStateFC/microstate_results.mat','results','templates','-v7.3');

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

%% ===== PLOT EXAMPLE =====
figure;
plot(results(1).labels);
ylim([1 K]);
xlabel('Time'); ylabel('State');
title('Microstate Sequence');

%% ===== FUNCTIONS =====

function labels = smooth_microstates(labels, minSamples)

N = length(labels);
i = 1;

while i <= N
    start_idx = i;
    current_val = labels(i);
    
    while i <= N && labels(i) == current_val
        i = i + 1;
    end
    
    end_idx = i - 1;
    segLen = end_idx - start_idx + 1;
    
    if segLen < minSamples
        if start_idx > 1 && end_idx < N
            % choose neighbor with longer segment
            prev_val = labels(start_idx-1);
            next_val = labels(end_idx+1);
            
            labels(start_idx:end_idx) = prev_val; % simple rule
        elseif start_idx > 1
            labels(start_idx:end_idx) = labels(start_idx-1);
        elseif end_idx < N
            labels(start_idx:end_idx) = labels(end_idx+1);
        end
    end
end
end


function stats = compute_stats(labels, Fs, K)

T = length(labels);

for k = 1:K
    
    idx = (labels == k);
    
    stats(k).coverage = sum(idx)/T;
    
    d = diff([0 idx 0]);
    startIdx = find(d==1);
    endIdx = find(d==-1)-1;
    
    durations = (endIdx - startIdx + 1)/Fs;
    
    if ~isempty(durations)
        stats(k).meanDuration = mean(durations)*1000;
        stats(k).occurrence = length(durations)/(T/Fs);
    else
        stats(k).meanDuration = 0;
        stats(k).occurrence = 0;
    end
end
end
