clc; clear; close all;

%% ===== CONFIG =====
rootPath = '/Users/suvamdey/Desktop/NSP/data/segmentedData';

subjects = {'003S','013AR'};
stages = {'M1','G1'};%{'M1','G1','M2','G2'};

Fs = 1000;
K = 5;

allMaps = [];

%% ===== LOOP OVER SUBJECTS + STAGES =====
for s = 1:length(subjects)
    
    subjectName = subjects{s};
    subjPath = fullfile(rootPath, subjectName, 'EEG');
    
    dates = dir(subjPath);
    dates = dates([dates.isdir] & ~startsWith({dates.name}, '.'));
    

        
    expDate = dates.name;
    
    for st = 1:length(stages)
        
        protocolName = stages{st};
        
        folderSegment = fullfile(rootPath, subjectName, ...
            'EEG', expDate, protocolName, ...
            'segmentedData', 'LFP');
        
        if ~exist(folderSegment,'dir'), continue; end
        
        %% ===== LOAD METADATA =====
        lfpInfo = load(fullfile(folderSegment,'lfpInfo.mat'));
        electrodes = lfpInfo.analogChannelsStored;
        
        numElectrodes = length(electrodes);
        
        % infer trials from one file
        testFile = fullfile(folderSegment, ['elec' num2str(electrodes(1)) '.mat']);
        tmp = load(testFile);
        
        [numTrials, numSamples] = size(tmp.analogData);
        
        %% ===== LOAD ALL ELECTRODES =====
        eegDataAll = zeros(numElectrodes, numTrials, numSamples);
        
        for i = 1:numElectrodes
            elecFile = fullfile(folderSegment, ...
                ['elec' num2str(electrodes(i)) '.mat']);
            
            if exist(elecFile,'file')
                dataVars = load(elecFile);
                eegDataAll(i,:,:) = dataVars.analogData;
            end
        end
        
        %% ===== RESHAPE TO CONTINUOUS =====
        totalSamples = numTrials * numSamples;
        data = reshape(eegDataAll, numElectrodes, totalSamples);
        
        %% ===== FILTER 0.5–45 Hz =====
        data = bandpass(data', [0.5 45], Fs)';
        
        %% ===== COMMON AVERAGE REFERENCE =====
        data = data - mean(data,1);
        
        %% ===== GFP =====
        GFP = std(data, 0, 1);
        
        %% ===== GFP PEAKS =====
        [~, locs] = findpeaks(GFP);
        
        if isempty(locs), continue; end
        
        % remove lowest 15%
        thr_low = prctile(GFP(locs),15);
        locs = locs(GFP(locs) > thr_low);
        
        % remove extreme peaks
        thr_high = mean(GFP) + 3*std(GFP);
        locs = locs(GFP(locs) < thr_high);
        
        if isempty(locs), continue; end
        
        maps = data(:, locs);
        
        %% ===== NORMALIZE MAPS =====
        maps = maps ./ vecnorm(maps);
        
        allMaps = [allMaps, maps];
        
    end
end

disp(['Total maps: ', num2str(size(allMaps,2))]);

[templates, labels, GEV] = microstate_kmeans(allMaps, K, 50);

disp(['Final GEV: ', num2str(GEV)]);

%% ===== SAVE MODEL =====
microstate = struct();

microstate.templates = templates;
microstate.labels_gfp = labels;   % labels only for GFP maps
microstate.GEV = GEV;
microstate.K = K;

microstate.nMaps = size(allMaps,2);
microstate.nChannels = size(allMaps,1);

microstate.description = 'Global microstate model from pooled GFP maps';

save('microstate_model.mat', 'microstate', '-v7.3');

%% Function Defns


function GEV = compute_GEV(maps, templates, labels)

GFP = std(maps,0,1);
GEV_num = 0;
GEV_den = sum(GFP.^2);

for t = 1:length(labels)
    k = labels(t);

    c = abs(corr(maps(:,t), templates(:,k)));
    GEV_num = GEV_num + (GFP(t)^2)*(c^2);
end

GEV = GEV_num / GEV_den;

end



function [templates, labels, GEV] = microstate_kmeans(maps, K, nrep)

if nargin < 3, nrep = 50; end

bestGEV = -inf;

for rep = 1:nrep

    % init
    idx = randperm(size(maps,2), K);
    templates = maps(:, idx);

    for iter = 1:100

        % assign (polarity invariant)
        corrMat = abs(templates' * maps);
        [~, labels] = max(corrMat, [], 1);

        % update
        for k = 1:K
            clusterMaps = maps(:, labels == k);

            if isempty(clusterMaps), continue; end

            [U,~,~] = svd(clusterMaps, 'econ');
            templates(:,k) = U(:,1);
        end
    end

    % compute GEV
    GEV = compute_GEV(maps, templates, labels);

    if GEV > bestGEV
        bestGEV = GEV;
        bestTemplates = templates;
        bestLabels = labels;
    end
end

templates = bestTemplates;
labels = bestLabels;
GEV = bestGEV;

end