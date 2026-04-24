clc; clear; close all;

%% ===== CONFIG =====
rootPath = '/Users/suvamdey/Desktop/NSP/data/segmentedData';

subjects = {'003S','013AR'};
stages = {'EO1', 'EC1', 'G1', 'M1','G1','EO2', 'EC2', 'M2'};

Fs = 1000;
kRange = 2:10;
nrep = 4;   % reduce for speed (paper uses 100)

% %% ===== FILTER DESIGN (FAST) =====
% d = designfilt('bandpassfir', ...
%     'FilterOrder', 500, ...
%     'CutoffFrequency1', 0.5, ...
%     'CutoffFrequency2', 45, ...
%     'SampleRate', Fs);
%% ===== FILTER DESIGN =====

% ---- Stage 1: Broadband cleaning (0.5–45 Hz) ----
[b1,a1] = butter(4, [0.5 45]/(Fs/2), 'bandpass');
band = 'gamma';
% ---- Stage 2: Band-specific ----
switch band

    case 'alpha'
        fRange = [8 13];

    case 'gamma'
        fRange = [30 45];

    case 'broadband'
        fRange = []; % skip second filter

    otherwise
        error('Invalid band');
end

if ~isempty(fRange)
    [b2,a2] = butter(4, fRange/(Fs/2), 'bandpass');
end

disp(['Band mode: ', band]);


%% ===== STORE SUBJECT TEMPLATES =====
subjectTemplates = {};

%% ===== LOOP =====
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
        % data = filtfilt(d, data')';
        %% ===== TWO-STAGE FILTERING =====

        % Stage 1
        data = filtfilt(b1, a1, data')';

        % CAR
        data = data - mean(data,1);

        % Stage 2 (optional)
        if ~isempty(fRange)
            data = filtfilt(b2, a2, data')';
        end
        % data = data - mean(data,1);
        
        %% ===== GFP =====
        GFP = std(data,0,1);
        [~, locs] = findpeaks(GFP);
        
        if isempty(locs), continue; end
        
        % thresholding
        locs = locs(GFP(locs) > prctile(GFP(locs),15));
        locs = locs(GFP(locs) < mean(GFP)+3*std(GFP));
        
        if isempty(locs), continue; end
        
        maps = data(:, locs);
        maps = maps ./ sqrt(sum(maps.^2,1));
        
        %% ===== STAGE 1: SUBJECT CLUSTERING =====
        bestGEV = -inf;
        bestTemplates = [];
        
        for k = kRange
            
            [T, labels, GEV] = microstate_kmeans(maps, k, nrep);
            
            if GEV > bestGEV
                bestGEV = GEV;
                bestTemplates = T;
                bestK = k;
            end
        end
        
        % store templates (equal contribution)
        subjectTemplates{end+1} = bestTemplates(:,1:bestK);
        
    end
end

%% ===== STAGE 2: GROUP CLUSTERING =====
allSubjMaps = cat(2, subjectTemplates{:});

Kfinal = 5; % or choose via KL
[groupTemplates, ~, ~] = microstate_kmeans(allSubjMaps, Kfinal, 50);

%% ===== SAVE =====
microstate.templates = groupTemplates;
microstate.K = Kfinal;

save('MicroStateFC/microstate_model_gamma.mat','microstate','-v7.3');

disp('DONE: 2-stage clustering complete');


%% Function Defns
function [templates, labels, GEV] = microstate_kmeans(maps, K, nrep)

bestGEV = -inf;

for rep = 1:nrep

    idx = randperm(size(maps,2), K);
    templates = maps(:, idx);

    for iter = 1:50

        % vectorized assignment
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

    % fast GEV
    GFP = std(maps,0,1);
    corrVals = abs(sum(templates(:,labels).*maps,1));
    GEV = sum((GFP.^2).*(corrVals.^2)) / sum(GFP.^2);

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