clc; clear; close all;

%% ============================================================
% FAST MICROSTATE TEMPLATE EXTRACTION
% Stage K=5 -> Subject final K=5
% Saves final subject templates only
%% ============================================================

rootPath = '/Users/anatapmitra/Desktop/NSP project/data/Segmented Data';
subjectListPath = '/Users/anatapmitra/Desktop/NSP project/BK1AllSubjectList.mat';
chanlocsPath = '/Users/anatapmitra/Desktop/actiCap64_UOL.mat';

outputRoot = '/Users/anatapmitra/Desktop/NSP project/MicroStateTemplates_fastK5_v2';

stages = {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2'};
bands = {'broadband', 'alpha', 'gamma'};

Fs = 1000;

Kstage = 5;
KsubjectFinal = 5;

nrepStage = 3;
nrepSubject = 20;

maxIterStage = 20;
maxIterSubject = 50;

maxGFPMapsPerStage = 8000;

rng(1234);

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

%% ===================== LOAD CHANLOCS =====================

tmpMontage = load(chanlocsPath);
chanlocs = tmpMontage.chanlocs;
numChannels = length(chanlocs);

fprintf('Loaded chanlocs: %d channels\n', numChannels);

%% ===================== SUBJECT SELECTION =====================

S = load(subjectListPath);

controlList = cellstr(string(S.controlList(:)));
meditatorList = cellstr(string(S.meditatorList(:)));

fixedControl = '003S';
fixedMeditator = '013AR';

controlPool = setdiff(controlList, {fixedControl}, 'stable');
meditatorPool = setdiff(meditatorList, {fixedMeditator}, 'stable');

selectedControls = [{fixedControl}; controlPool(randperm(length(controlPool), 4))];
selectedMeditators = [{fixedMeditator}; meditatorPool(randperm(length(meditatorPool), 4))];

subjects = [selectedControls; selectedMeditators];

groupLabels = [repmat({'control'}, length(selectedControls), 1); ...
               repmat({'meditator'}, length(selectedMeditators), 1)];

save(fullfile(outputRoot, 'selected_subjects.mat'), ...
    'selectedControls', 'selectedMeditators', 'subjects', 'groupLabels');

disp('Selected controls:');
disp(selectedControls);

disp('Selected meditators:');
disp(selectedMeditators);

%% ===================== MAIN LOOP =====================

summaryResults = struct();
summaryCounter = 1;

for b = 1:length(bands)

    bandName = bands{b};

    fprintf('\n====================================================\n');
    fprintf('BAND: %s\n', bandName);
    fprintf('====================================================\n');

    bandOutDir = fullfile(outputRoot, bandName);
    if ~exist(bandOutDir, 'dir')
        mkdir(bandOutDir);
    end

    for s = 1:length(subjects)

        subjectName = subjects{s};
        groupName = groupLabels{s};

        fprintf('\n--------------------------------------------\n');
        fprintf('Subject: %s | Group: %s | Band: %s\n', subjectName, groupName, bandName);
        fprintf('--------------------------------------------\n');

        subjectPath = fullfile(rootPath, subjectName, 'EEG');

        if ~exist(subjectPath, 'dir')
            warning('Missing subject folder: %s', subjectPath);
            continue;
        end

        dateDirs = dir(subjectPath);
        dateDirs = dateDirs([dateDirs.isdir] & ~startsWith({dateDirs.name}, '.'));

        if isempty(dateDirs)
            warning('No date folder for subject %s', subjectName);
            continue;
        end

        expDate = dateDirs(1).name;

        subjectOutDir = fullfile(bandOutDir, groupName, subjectName);
        figDir = fullfile(subjectOutDir, 'figures');

        if ~exist(subjectOutDir, 'dir')
            mkdir(subjectOutDir);
        end

        if ~exist(figDir, 'dir')
            mkdir(figDir);
        end

        subjectStageTemplates = [];
        stageInfo = struct();
        stageCounter = 1;

        %% ===================== STAGE LOOP =====================

        for st = 1:length(stages)

            stageName = stages{st};

            fprintf('\nStage: %s\n', stageName);

            lfpPath = fullfile(rootPath, subjectName, ...
                'EEG', expDate, stageName, ...
                'segmentedData', 'LFP');

            if ~exist(lfpPath, 'dir')
                warning('Missing LFP folder: %s', lfpPath);
                continue;
            end

            %% Load directly as channels x time
            data2D = load_stage_eeg_2d(lfpPath, numChannels);

            if isempty(data2D)
                warning('Could not load data for %s %s', subjectName, stageName);
                continue;
            end

            fprintf('Loaded data2D: [%d x %d]\n', size(data2D,1), size(data2D,2));

            %% Preprocess
            data2D = preprocess_eeg_for_microstates(data2D, Fs, bandName);

            %% GFP peak maps with cap
            [maps, GFP, locs] = extract_gfp_peak_maps_fast(data2D, maxGFPMapsPerStage);

            clear data2D;

            if isempty(maps)
                warning('No valid GFP maps for %s %s', subjectName, stageName);
                continue;
            end

            fprintf('GFP maps used: %d\n', size(maps,2));

            if size(maps,2) < Kstage
                warning('Too few maps for K=%d', Kstage);
                continue;
            end

            %% Fast modified k-means
            [stageTemplates, stageLabels, stageGEV] = microstate_kmeans_fast( ...
                maps, Kstage, nrepStage, maxIterStage);

            clear maps;

            fprintf('Stage K=%d | GEV=%.4f\n', Kstage, stageGEV);

            subjectStageTemplates = [subjectStageTemplates, stageTemplates];

            stageInfo(stageCounter).stage = stageName;
            stageInfo(stageCounter).K_stage = Kstage;
            stageInfo(stageCounter).GEV_stage = stageGEV;
            stageInfo(stageCounter).numGFPMapsUsed = length(locs);

            stageCounter = stageCounter + 1;

        end

        %% ===================== SUBJECT FINAL CLUSTERING =====================

        if isempty(subjectStageTemplates)
            warning('No valid stage templates for %s %s', subjectName, bandName);
            continue;
        end

        fprintf('\nSubject final clustering: input [%d x %d]\n', ...
            size(subjectStageTemplates,1), size(subjectStageTemplates,2));

        [subjectTemplates, subjectLabels, subjectGEV] = microstate_kmeans_fast( ...
            subjectStageTemplates, KsubjectFinal, nrepSubject, maxIterSubject);

        subjectTemplates = normalize_maps(subjectTemplates);
        templateCorr = abs(subjectTemplates' * subjectTemplates);

        %% ===================== SAVE =====================

        subjectMicrostate = struct();

        subjectMicrostate.subject = subjectName;
        subjectMicrostate.group = groupName;
        subjectMicrostate.band = bandName;

        subjectMicrostate.Kstage = Kstage;
        subjectMicrostate.K = KsubjectFinal;

        subjectMicrostate.templates = subjectTemplates;
        subjectMicrostate.templateCorr = templateCorr;
        subjectMicrostate.GEV_subject = subjectGEV;
        subjectMicrostate.subjectLabels = subjectLabels;
        subjectMicrostate.stageInfo = stageInfo;

        subjectMicrostate.Fs = Fs;
        subjectMicrostate.stages = stages;
        subjectMicrostate.expDate = expDate;
        subjectMicrostate.maxGFPMapsPerStage = maxGFPMapsPerStage;
        subjectMicrostate.nrepStage = nrepStage;
        subjectMicrostate.nrepSubject = nrepSubject;
        subjectMicrostate.maxIterStage = maxIterStage;
        subjectMicrostate.maxIterSubject = maxIterSubject;
        subjectMicrostate.preprocessing = get_preprocessing_description(bandName);

        savePath = fullfile(subjectOutDir, ...
            sprintf('%s_%s_finalSubjectTemplates_K5_FAST.mat', subjectName, bandName));

        save(savePath, 'subjectMicrostate', '-v7.3');

        fprintf('Saved: %s\n', savePath);

        %% ===================== PLOTS =====================

        plot_subject_microstate_templates(subjectTemplates, chanlocs, ...
            subjectName, groupName, bandName, figDir);

        plot_template_correlation_matrix(templateCorr, ...
            subjectName, groupName, bandName, figDir);

        summaryResults(summaryCounter).subject = subjectName;
        summaryResults(summaryCounter).group = groupName;
        summaryResults(summaryCounter).band = bandName;
        summaryResults(summaryCounter).GEV_subject = subjectGEV;
        summaryResults(summaryCounter).savePath = savePath;

        summaryCounter = summaryCounter + 1;

    end
end

summaryPath = fullfile(outputRoot, 'microstate_template_extraction_summary_FAST.mat');

save(summaryPath, 'summaryResults', 'subjects', 'groupLabels', ...
    'selectedControls', 'selectedMeditators', ...
    'bands', 'stages', 'Kstage', 'KsubjectFinal', ...
    'nrepStage', 'nrepSubject', 'maxGFPMapsPerStage', '-v7.3');

fprintf('\nDONE.\nSummary saved at:\n%s\n', summaryPath);


%% ============================================================
% FUNCTIONS
%% ============================================================

function data2D = load_stage_eeg_2d(lfpPath, numChannels)

    data2D = [];

    firstFile = fullfile(lfpPath, 'elec1.mat');

    if ~exist(firstFile, 'file')
        warning('elec1.mat missing: %s', firstFile);
        return;
    end

    tmp = load(firstFile);

    if ~isfield(tmp, 'analogData')
        warning('analogData missing in elec1.mat');
        return;
    end

    [numTrials, numSamples] = size(tmp.analogData);
    totalTime = numTrials * numSamples;

    data2D = nan(numChannels, totalTime);

    for ch = 1:numChannels

        elecFile = fullfile(lfpPath, sprintf('elec%d.mat', ch));

        if ~exist(elecFile, 'file')
            warning('Missing %s', elecFile);
            continue;
        end

        D = load(elecFile);

        if ~isfield(D, 'analogData')
            warning('analogData missing in %s', elecFile);
            continue;
        end

        X = D.analogData;

        if ~isequal(size(X), [numTrials, numSamples])
            warning('Size mismatch in %s', elecFile);
            continue;
        end

        data2D(ch,:) = reshape(X, 1, []);

    end

    validChannels = all(~isnan(data2D), 2);
    data2D = data2D(validChannels, :);

end


function data = preprocess_eeg_for_microstates(data, Fs, bandName)

    validTime = all(~isnan(data), 1);
    data = data(:, validTime);

    data = data - mean(data, 2);

    [b1, a1] = butter(4, [0.5 45] / (Fs/2), 'bandpass');
    data = filtfilt(b1, a1, data')';

    data = data - mean(data, 1);

    switch lower(bandName)

        case 'broadband'
            % nothing extra

        case 'alpha'
            [b2, a2] = butter(4, [8 13] / (Fs/2), 'bandpass');
            data = filtfilt(b2, a2, data')';

        case 'gamma'
            [b2, a2] = butter(4, [30 45] / (Fs/2), 'bandpass');
            data = filtfilt(b2, a2, data')';

        otherwise
            error('Unknown band: %s', bandName);
    end

    data = data - mean(data, 1);

end


function [maps, GFP, locs] = extract_gfp_peak_maps_fast(data, maxMaps)

    GFP = std(data, 0, 1);

    [~, locs] = findpeaks(GFP);

    if isempty(locs)
        maps = [];
        return;
    end

    locs = locs(GFP(locs) > prctile(GFP(locs), 15));
    locs = locs(GFP(locs) < mean(GFP) + 3 * std(GFP));

    if isempty(locs)
        maps = [];
        return;
    end

    % Keep strongest GFP peaks if too many.
    % Better than random because high-GFP maps are more reliable for microstates.
    if length(locs) > maxMaps
        [~, ord] = sort(GFP(locs), 'descend');
        keepIdx = ord(1:maxMaps);
        locs = locs(keepIdx);
    end

    maps = data(:, locs);
    maps = normalize_maps(maps);

end


function maps = normalize_maps(maps)

    denom = sqrt(sum(maps.^2, 1));
    denom(denom == 0) = eps;
    maps = maps ./ denom;

end


function [templates, labels, GEV] = microstate_kmeans_fast(maps, K, nrep, maxIter)

    maps = normalize_maps(maps);

    bestGEV = -inf;
    bestTemplates = [];
    bestLabels = [];

    N = size(maps, 2);

    for rep = 1:nrep

        idx = randperm(N, K);
        templates = maps(:, idx);
        templates = normalize_maps(templates);

        labels = zeros(1, N);

        for iter = 1:maxIter

            oldLabels = labels;

            % assignment: polarity-invariant spatial correlation
            corrMat = abs(templates' * maps);
            [~, labels] = max(corrMat, [], 1);

            % update: first principal component per cluster
            % faster than svd(clusterMaps), same dominant spatial direction
            for k = 1:K

                clusterMaps = maps(:, labels == k);

                if isempty(clusterMaps)
                    templates(:,k) = maps(:, randi(N));
                    continue;
                end

                C = clusterMaps * clusterMaps';
                [V, D] = eig(C, 'vector');
                [~, imax] = max(D);
                templates(:,k) = V(:, imax);

            end

            templates = normalize_maps(templates);

            if isequal(labels, oldLabels)
                break;
            end

        end

        corrVals = abs(sum(templates(:, labels) .* maps, 1));
        GFP_proxy = std(maps, 0, 1);

        GEV_rep = sum((GFP_proxy.^2) .* (corrVals.^2)) / sum(GFP_proxy.^2);

        fprintf('    rep %d/%d | iter %d | GEV %.4f\n', ...
            rep, nrep, iter, GEV_rep);

        if GEV_rep > bestGEV
            bestGEV = GEV_rep;
            bestTemplates = templates;
            bestLabels = labels;
        end

    end

    templates = bestTemplates;
    labels = bestLabels;
    GEV = bestGEV;

end


function desc = get_preprocessing_description(bandName)

    desc = struct();

    desc.step1 = 'Load elec1.mat to elec64.mat corresponding to chanlocs(1:64)';
    desc.step2 = 'Concatenate trials into channels x time';
    desc.step3 = 'Channel-wise demean';
    desc.step4 = 'Butterworth 4th-order 0.5-45 Hz bandpass';
    desc.step5 = 'Common average reference';
    desc.step6 = 'Optional second band-specific filter';
    desc.step7 = 'Common average reference again';
    desc.step8 = 'GFP peak extraction';
    desc.step9 = 'Keep strongest GFP peaks up to maxGFPMapsPerStage';
    desc.step10 = 'GFP map normalization';

    switch lower(bandName)
        case 'broadband'
            desc.bandFilter = 'Broadband only: 0.5-45 Hz';
        case 'alpha'
            desc.bandFilter = '0.5-45 Hz followed by 8-13 Hz';
        case 'gamma'
            desc.bandFilter = '0.5-45 Hz followed by 30-45 Hz';
    end

end


function plot_subject_microstate_templates(templates, chanlocs, subjectName, groupName, bandName, figDir)

    K = size(templates, 2);

    fig = figure('Color', 'w', 'Position', [100 100 1200 600]);

    clim = max(abs(templates(:)));

    for k = 1:K
        subplot(2, ceil(K/2), k);
        topoplot(templates(:,k), chanlocs, 'electrodes', 'on');
        title(sprintf('M%d', k));
        caxis([-clim clim]);
        colorbar;
    end

    sgtitle(sprintf('%s | %s | %s | Final 5 Microstates', ...
        subjectName, groupName, bandName), 'Interpreter', 'none');

    saveas(fig, fullfile(figDir, ...
        sprintf('%s_%s_final_microstate_topoplots_FAST.png', subjectName, bandName)));

    savefig(fig, fullfile(figDir, ...
        sprintf('%s_%s_final_microstate_topoplots_FAST.fig', subjectName, bandName)));

    close(fig);

end


function plot_template_correlation_matrix(corrMat, subjectName, groupName, bandName, figDir)

    fig = figure('Color', 'w', 'Position', [100 100 650 550]);

    imagesc(corrMat);
    axis square;
    colorbar;
    caxis([0 1]);

    xlabel('Microstate template');
    ylabel('Microstate template');

    title(sprintf('%s | %s | %s | abs(template correlation)', ...
        subjectName, groupName, bandName), 'Interpreter', 'none');

    xticks(1:size(corrMat,1));
    yticks(1:size(corrMat,1));

    for i = 1:size(corrMat,1)
        for j = 1:size(corrMat,2)
            text(j, i, sprintf('%.2f', corrMat(i,j)), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 10, ...
                'Color', 'k');
        end
    end

    saveas(fig, fullfile(figDir, ...
        sprintf('%s_%s_template_correlation_matrix_FAST.png', subjectName, bandName)));

    savefig(fig, fullfile(figDir, ...
        sprintf('%s_%s_template_correlation_matrix_FAST.fig', subjectName, bandName)));

    close(fig);

end