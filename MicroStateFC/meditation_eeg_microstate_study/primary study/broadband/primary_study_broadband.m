%% Code for study3
%% CODE FOR Broadband case: STUDY 2
% State-specific microstate signatures of meditators vs controls
% Full10 only: 5 controls vs 5 meditators

clc; clear; close all;

%% ===================== CONFIG =====================

runMode = 'full10';
bandName = 'broadband';

rootPath = '/Users/anatapmitra/Desktop/NSP project/data/Segmented Data';
templateRoot = '/Users/anatapmitra/Desktop/NSP project/MicroStateTemplates_fastK5_v2';
chanlocsPath = '/Users/anatapmitra/Desktop/actiCap64_UOL.mat';

outputRoot = fullfile('/Users/anatapmitra/Desktop/NSP project/Study2_StateSpecific_Broadband_Results', runMode);

Fs = 1000;
numChannels = 64;

allStages = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
K = 5;
microstateNames = {'A','B','C','D','E'};
features = {'Coverage','Duration','GlobalWPLI'};

windowLengthSec = 30;
windowLengthSamples = windowLengthSec * Fs;
minSamplesForWPLI = round(1.0 * Fs);

nPerm = 5000;
rng(1234);

figDir = fullfile(outputRoot, 'figures');
tableDir = fullfile(outputRoot, 'tables');
modelDir = fullfile(outputRoot, 'models');

if ~exist(outputRoot, 'dir'), mkdir(outputRoot); end
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(tableDir, 'dir'), mkdir(tableDir); end
if ~exist(modelDir, 'dir'), mkdir(modelDir); end

%% ===================== LOAD CHANLOCS =====================

tmpMontage = load(chanlocsPath);
chanlocs = tmpMontage.chanlocs;

%% ===================== LOAD FULL10 SUBJECTS =====================

S = load(fullfile(templateRoot, 'selected_subjects.mat'));

if isfield(S, 'subjects') && isfield(S, 'groupLabels')
    subjects = S.subjects(:);
    groupLabels = S.groupLabels(:);
else
    error('Could not find subjects/groupLabels in selected_subjects.mat');
end

fprintf('\nSubjects used for Study 2:\n');
disp(table(subjects, groupLabels));

%% ============================================================
% 1. LOAD SUBJECT TEMPLATES AND DERIVE COMMON BROADBAND TEMPLATES
%% ============================================================

fprintf('\nLoading saved broadband subject templates...\n');

allSubjectTemplates = [];

for s = 1:length(subjects)

    subjectName = subjects{s};
    groupName = groupLabels{s};

    templatePath = fullfile(templateRoot, bandName, groupName, subjectName, ...
        sprintf('%s_%s_finalSubjectTemplates_K5_FAST.mat', subjectName, bandName));

    if ~exist(templatePath, 'file')
        error('Missing template file:\n%s', templatePath);
    end

    T = load(templatePath);

    if ~isfield(T, 'subjectMicrostate')
        error('subjectMicrostate missing in %s', templatePath);
    end

    subjTemplates = T.subjectMicrostate.templates;

    if size(subjTemplates,1) ~= numChannels || size(subjTemplates,2) ~= K
        error('Unexpected template size in %s', templatePath);
    end

    allSubjectTemplates = [allSubjectTemplates, subjTemplates];

end

fprintf('All subject template matrix: [%d x %d]\n', ...
    size(allSubjectTemplates,1), size(allSubjectTemplates,2));

fprintf('Deriving common K=5 broadband templates...\n');

[commonTemplates, ~, commonGEV] = microstate_kmeans_fast(allSubjectTemplates, K, 100, 100);
commonTemplates = normalize_maps(commonTemplates);

commonModel = struct();
commonModel.band = bandName;
commonModel.runMode = runMode;
commonModel.study = 'Study2_StateSpecific_GroupSignature';
commonModel.subjects = subjects;
commonModel.groupLabels = groupLabels;
commonModel.templates = commonTemplates;
commonModel.K = K;
commonModel.GEV = commonGEV;

save(fullfile(modelDir, sprintf('commonTemplates_%s_study2.mat', bandName)), ...
    'commonModel', '-v7.3');

plot_common_templates(commonTemplates, chanlocs, microstateNames, ...
    sprintf('Study 2 | Common %s templates', bandName), ...
    fullfile(figDir, sprintf('P01_common_templates_%s_study2.png', bandName)));

%% ============================================================
% 2. BACKFIT + WINDOW-LEVEL FEATURE EXTRACTION
%% ============================================================

fprintf('\nBackfitting and extracting broadband features...\n');

featureRows = {};

for s = 1:length(subjects)

    subjectName = subjects{s};
    groupName = groupLabels{s};

    subjectPath = fullfile(rootPath, subjectName, 'EEG');
    dateDirs = dir(subjectPath);
    dateDirs = dateDirs([dateDirs.isdir] & ~startsWith({dateDirs.name}, '.'));

    if isempty(dateDirs)
        warning('No date folder for subject %s', subjectName);
        continue;
    end

    expDate = dateDirs(1).name;

    for st = 1:length(allStages)

        stageName = allStages{st};

        fprintf('\nSubject %s | %s | %s\n', subjectName, groupName, stageName);

        lfpPath = fullfile(rootPath, subjectName, ...
            'EEG', expDate, stageName, ...
            'segmentedData', 'LFP');

        if ~exist(lfpPath, 'dir')
            warning('Missing LFP folder: %s', lfpPath);
            continue;
        end

        data2D = load_stage_eeg_2d(lfpPath, numChannels);

        if isempty(data2D)
            warning('Could not load stage data.');
            continue;
        end

        data2D = preprocess_broadband(data2D, Fs);

        nTotal = size(data2D, 2);
        nWindows = floor(nTotal / windowLengthSamples);

        fprintf('Total samples: %d | Windows: %d\n', nTotal, nWindows);

        for w = 1:nWindows

            idx1 = (w-1)*windowLengthSamples + 1;
            idx2 = w*windowLengthSamples;

            winData = data2D(:, idx1:idx2);

            labels = backfit_microstates(winData, commonTemplates);

            [coverage, duration] = compute_temporal_features(labels, K, Fs);

            globalWPLI = compute_global_wpli_by_microstate(winData, labels, K, minSamplesForWPLI);

            for m = 1:K
                featureRows(end+1,:) = { ...
                    subjectName, groupName, bandName, stageName, w, ...
                    microstateNames{m}, m, ...
                    coverage(m), duration(m), globalWPLI(m) ...
                };
            end
        end

        clear data2D;

    end
end

featureTable = cell2table(featureRows, ...
    'VariableNames', {'Subject','Group','Band','Stage','Window', ...
                      'Microstate','MicrostateID', ...
                      'Coverage','Duration','GlobalWPLI'});

writetable(featureTable, fullfile(tableDir, sprintf('window_features_%s_study2.csv', bandName)));
save(fullfile(tableDir, sprintf('window_features_%s_study2.mat', bandName)), ...
    'featureTable', '-v7.3');

%% ============================================================
% 3. SUBJECT-LEVEL AGGREGATION
%% ============================================================

fprintf('\nAggregating windows to subject-level stage averages...\n');

subjectFeatureTable = aggregate_subject_level(featureTable, features);

writetable(subjectFeatureTable, fullfile(tableDir, sprintf('subject_level_features_%s_study2.csv', bandName)));
save(fullfile(tableDir, sprintf('subject_level_features_%s_study2.mat', bandName)), ...
    'subjectFeatureTable', '-v7.3');

%% ============================================================
% 4. STAGE-SPECIFIC GROUP COMPARISONS
%% ============================================================

fprintf('\nRunning stage-specific group permutation tests...\n');

groupRows = {};

for f = 1:length(features)

    featName = features{f};

    for st = 1:length(allStages)

        stageName = allStages{st};

        for m = 1:K

            msName = microstateNames{m};

            idxC = strcmp(subjectFeatureTable.Group, 'control') & ...
                   strcmp(subjectFeatureTable.Stage, stageName) & ...
                   subjectFeatureTable.MicrostateID == m;

            idxM = strcmp(subjectFeatureTable.Group, 'meditator') & ...
                   strcmp(subjectFeatureTable.Stage, stageName) & ...
                   subjectFeatureTable.MicrostateID == m;

            x = subjectFeatureTable.(featName)(idxC);
            y = subjectFeatureTable.(featName)(idxM);

            x = x(~isnan(x));
            y = y(~isnan(y));

            if length(x) < 2 || length(y) < 2
                meanC = mean(x, 'omitnan');
                meanM = mean(y, 'omitnan');
                delta = meanM - meanC;
                effectD = NaN;
                pval = NaN;
            else
                [pval, delta, effectD] = permutation_test_between_groups(x, y, nPerm);
                meanC = mean(x, 'omitnan');
                meanM = mean(y, 'omitnan');
            end

            groupRows(end+1,:) = {featName, stageName, msName, m, ...
                meanC, meanM, delta, effectD, pval, length(x), length(y)};
        end
    end
end

groupStatsTable = cell2table(groupRows, ...
    'VariableNames', {'Feature','Stage','Microstate','MicrostateID', ...
                      'ControlMean','MeditatorMean','Delta_MeditatorMinusControl', ...
                      'CohensD','PermP','N_Control','N_Meditator'});

writetable(groupStatsTable, fullfile(tableDir, sprintf('stagewise_group_tests_%s_study2.csv', bandName)));
save(fullfile(tableDir, sprintf('stagewise_group_tests_%s_study2.mat', bandName)), ...
    'groupStatsTable', '-v7.3');

%% ============================================================
% 5. LMM: STAGE * GROUP PER FEATURE × MICROSTATE
%% ============================================================

fprintf('\nRunning LMM: Response ~ Stage*Group + (1|Subject)...\n');

lmmRows = {};

for f = 1:length(features)

    featName = features{f};

    for m = 1:K

        msName = microstateNames{m};

        idx = subjectFeatureTable.MicrostateID == m;
        tmpT = subjectFeatureTable(idx, :);

        response = tmpT.(featName);
        valid = ~isnan(response);

        LMMTable = table();
        LMMTable.Response = response(valid);
        LMMTable.Subject = categorical(tmpT.Subject(valid));
        LMMTable.Group = categorical(tmpT.Group(valid));
        LMMTable.Stage = categorical(tmpT.Stage(valid));

        LMMTable.Stage = reordercats(LMMTable.Stage, allStages);
        LMMTable.Group = reordercats(LMMTable.Group, {'control','meditator'});

        try
            fprintf('Running LMM: %s | MS%s | N=%d\n', featName, msName, height(LMMTable));

            lme = fitlme(LMMTable, 'Response ~ Stage*Group + (1|Subject)');
            anovaT = anova(lme);

            if istable(anovaT)
                anovaOut = anovaT;
            else
                anovaOut = dataset2table(anovaT);
            end

            save(fullfile(modelDir, sprintf('LMM_%s_MS%s_%s_study2.mat', ...
                featName, msName, bandName)), ...
                'lme', 'anovaT', 'LMMTable');

            writetable(anovaOut, fullfile(tableDir, sprintf('LMM_ANOVA_%s_MS%s_%s_study2.csv', ...
                featName, msName, bandName)));

            termNames = string(anovaOut.Term);

            pStage = get_anova_p(anovaOut, 'Stage');
            pGroup = get_anova_p(anovaOut, 'Group');

            interactionRows = contains(termNames, 'Stage:Group') | contains(termNames, 'Group:Stage');
            if any(interactionRows)
                pInteraction = anovaOut.pValue(find(interactionRows, 1));
            else
                pInteraction = NaN;
            end

            lmmRows(end+1,:) = {featName, msName, m, height(LMMTable), ...
                pStage, pGroup, pInteraction};

        catch ME
            warning('LMM failed for %s MS%s: %s', featName, msName, ME.message);

            lmmRows(end+1,:) = {featName, msName, m, height(LMMTable), ...
                NaN, NaN, NaN};
        end
    end
end

lmmSummary = cell2table(lmmRows, ...
    'VariableNames', {'Feature','Microstate','MicrostateID','N', ...
                      'StageP','GroupP','StageGroupInteractionP'});

writetable(lmmSummary, fullfile(tableDir, sprintf('LMM_summary_%s_study2.csv', bandName)));
save(fullfile(tableDir, sprintf('LMM_summary_%s_study2.mat', bandName)), ...
    'lmmSummary', '-v7.3');

%% ============================================================
% 6. DOMINANCE PATTERN USING COVERAGE
%% ============================================================

fprintf('\nComputing microstate dominance patterns...\n');

dominanceTable = compute_dominance_table(subjectFeatureTable, allStages, microstateNames);

writetable(dominanceTable, fullfile(tableDir, sprintf('coverage_dominance_%s_study2.csv', bandName)));
save(fullfile(tableDir, sprintf('coverage_dominance_%s_study2.mat', bandName)), ...
    'dominanceTable', '-v7.3');

%% ============================================================
% 7. PLOTS: MINIMAL STUDY 2 FIGURE SET
%% ============================================================

fprintf('\nGenerating Study 2 plots...\n');

% P02-P04: heatmaps of group differences
for f = 1:length(features)
    featName = features{f};

    plot_group_delta_heatmap(groupStatsTable, featName, allStages, microstateNames, ...
        sprintf('Study 2 | %s | Meditator - Control', featName), ...
        fullfile(figDir, sprintf('P0%d_heatmap_delta_%s_%s.png', f+1, featName, bandName)));
end

% P05: dominance pattern table/heatmap
plot_dominance_pattern(dominanceTable, allStages, ...
    fullfile(figDir, sprintf('P05_dominance_pattern_%s.png', bandName)));

% P06: M1 coverage profile
plot_stage_profile(subjectFeatureTable, 'Coverage', 'M1', microstateNames, ...
    'Study 2 | M1 Coverage profile', ...
    fullfile(figDir, sprintf('P06_M1_coverage_profile_%s.png', bandName)));

% P07: EC average coverage profile
plot_composite_stage_profile(subjectFeatureTable, 'Coverage', {'EC1','EC2'}, microstateNames, ...
    'Study 2 | EC1/EC2 average Coverage profile', ...
    fullfile(figDir, sprintf('P07_EC_avg_coverage_profile_%s.png', bandName)));

% P08: gamma average GlobalWPLI profile
plot_composite_stage_profile(subjectFeatureTable, 'GlobalWPLI', {'G1','G2'}, microstateNames, ...
    'Study 2 | G1/G2 average GlobalWPLI profile', ...
    fullfile(figDir, sprintf('P08_G_avg_GlobalWPLI_profile_%s.png', bandName)));

fprintf('\nDONE Study 2 broadband analysis.\nResults saved at:\n%s\n', outputRoot);


%% ============================================================
% FUNCTIONS
%% ============================================================

function data2D = load_stage_eeg_2d(lfpPath, numChannels)

    data2D = [];

    firstFile = fullfile(lfpPath, 'elec1.mat');

    if ~exist(firstFile, 'file')
        return;
    end

    tmp = load(firstFile);

    if ~isfield(tmp, 'analogData')
        return;
    end

    [numTrials, numSamples] = size(tmp.analogData);
    totalTime = numTrials * numSamples;

    data2D = nan(numChannels, totalTime);

    for ch = 1:numChannels

        elecFile = fullfile(lfpPath, sprintf('elec%d.mat', ch));

        if ~exist(elecFile, 'file')
            continue;
        end

        D = load(elecFile);

        if ~isfield(D, 'analogData')
            continue;
        end

        X = D.analogData;

        if ~isequal(size(X), [numTrials, numSamples])
            continue;
        end

        data2D(ch,:) = reshape(X, 1, []);

    end

    validChannels = all(~isnan(data2D), 2);
    data2D = data2D(validChannels, :);

end


function data = preprocess_broadband(data, Fs)

    validTime = all(~isnan(data), 1);
    data = data(:, validTime);

    data = data - mean(data, 2);

    [b, a] = butter(4, [0.5 45] / (Fs/2), 'bandpass');
    data = filtfilt(b, a, data')';

    data = data - mean(data, 1);

end


function labels = backfit_microstates(data, templates)

    templates = normalize_maps(templates);

    data = data - mean(data, 1);

    denom = sqrt(sum(data.^2, 1));
    denom(denom == 0) = eps;
    dataNorm = data ./ denom;

    corrMat = abs(templates' * dataNorm);
    [~, labels] = max(corrMat, [], 1);

end


function [coverage, duration] = compute_temporal_features(labels, K, Fs)

    N = length(labels);

    coverage = nan(1, K);
    duration = nan(1, K);

    for k = 1:K

        coverage(k) = sum(labels == k) / N;

        runs = get_state_runs(labels, k);

        if isempty(runs)
            duration(k) = NaN;
        else
            runLengths = runs(:,2) - runs(:,1) + 1;
            duration(k) = mean(runLengths) / Fs;
        end

    end

end


function runs = get_state_runs(labels, stateID)

    mask = labels == stateID;
    d = diff([false, mask, false]);

    starts = find(d == 1);
    ends = find(d == -1) - 1;

    runs = [starts(:), ends(:)];

end


function globalWPLI = compute_global_wpli_by_microstate(data, labels, K, minSamples)

    nCh = size(data,1);
    globalWPLI = nan(1,K);

    analytic = hilbert(data')';
    phaseData = angle(analytic);

    upperMask = triu(true(nCh), 1);

    for k = 1:K

        idx = labels == k;

        if sum(idx) < minSamples
            globalWPLI(k) = NaN;
            continue;
        end

        phaseSeg = phaseData(:, idx);

        wpliMat = compute_wpli_matrix_from_phase(phaseSeg);

        vals = wpliMat(upperMask);
        globalWPLI(k) = mean(vals, 'omitnan');

    end

end


function wpliMat = compute_wpli_matrix_from_phase(phaseData)

    nCh = size(phaseData,1);
    wpliMat = nan(nCh, nCh);

    for i = 1:nCh

        phi_i = phaseData(i,:);

        for j = i+1:nCh

            phaseDiff = phi_i - phaseData(j,:);
            imags = sin(phaseDiff);

            num = abs(mean(imags));
            den = mean(abs(imags));

            if den == 0
                val = NaN;
            else
                val = num / den;
            end

            wpliMat(i,j) = val;
            wpliMat(j,i) = val;

        end
    end

    wpliMat(1:nCh+1:end) = NaN;

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

            corrMat = abs(templates' * maps);
            [~, labels] = max(corrMat, [], 1);

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


function maps = normalize_maps(maps)

    denom = sqrt(sum(maps.^2, 1));
    denom(denom == 0) = eps;
    maps = maps ./ denom;

end


function subjectFeatureTable = aggregate_subject_level(featureTable, features)

    keyVars = {'Subject','Group','Band','Stage','Microstate','MicrostateID'};

    G = groupsummary(featureTable, keyVars, 'mean', features);

    subjectFeatureTable = G(:, keyVars);

    for f = 1:length(features)
        oldName = ['mean_' features{f}];
        subjectFeatureTable.(features{f}) = G.(oldName);
    end

end


function [pval, delta, effectD] = permutation_test_between_groups(controlVals, meditatorVals, nPerm)

    x = controlVals(:);
    y = meditatorVals(:);

    x = x(~isnan(x));
    y = y(~isnan(y));

    delta = mean(y, 'omitnan') - mean(x, 'omitnan');

    pooledStd = sqrt(((length(x)-1)*var(x, 'omitnan') + ...
                      (length(y)-1)*var(y, 'omitnan')) / ...
                      (length(x) + length(y) - 2));

    if pooledStd == 0
        effectD = NaN;
    else
        effectD = delta / pooledStd;
    end

    combined = [x; y];
    nX = length(x);

    permDiffs = nan(nPerm, 1);

    for p = 1:nPerm
        idx = randperm(length(combined));

        xPerm = combined(idx(1:nX));
        yPerm = combined(idx(nX+1:end));

        permDiffs(p) = mean(yPerm, 'omitnan') - mean(xPerm, 'omitnan');
    end

    pval = mean(abs(permDiffs) >= abs(delta));

end


function p = get_anova_p(anovaOut, termName)

    p = NaN;
    termNames = string(anovaOut.Term);
    idx = strcmp(termNames, termName);

    if any(idx)
        p = anovaOut.pValue(find(idx, 1));
    end

end


function dominanceTable = compute_dominance_table(T, allStages, microstateNames)

    groups = {'control','meditator'};
    rows = {};

    for g = 1:length(groups)

        groupName = groups{g};

        for st = 1:length(allStages)

            stageName = allStages{st};

            covMeans = nan(1, length(microstateNames));

            for m = 1:length(microstateNames)

                idx = strcmp(T.Group, groupName) & ...
                      strcmp(T.Stage, stageName) & ...
                      T.MicrostateID == m;

                covMeans(m) = mean(T.Coverage(idx), 'omitnan');
            end

            [sortedVals, sortedIdx] = sort(covMeans, 'descend');

            domID = sortedIdx(1);
            secondID = sortedIdx(2);

            rows(end+1,:) = {groupName, stageName, ...
                microstateNames{domID}, domID, sortedVals(1), ...
                microstateNames{secondID}, secondID, sortedVals(2)};
        end
    end

    dominanceTable = cell2table(rows, ...
        'VariableNames', {'Group','Stage','DominantMS','DominantMSID','DominantCoverage', ...
                          'SecondMS','SecondMSID','SecondCoverage'});

end


function plot_common_templates(templates, chanlocs, microstateNames, figTitle, savePath)

    K = size(templates,2);

    fig = figure('Color','w','Position',[100 100 1200 600]);

    clim = max(abs(templates(:)));

    for k = 1:K
        subplot(2, ceil(K/2), k);
        topoplot(templates(:,k), chanlocs, 'electrodes','on');
        title(sprintf('MS-%s', microstateNames{k}));
        caxis([-clim clim]);
        colorbar;
    end

    sgtitle(figTitle, 'Interpreter','none');

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end


function plot_group_delta_heatmap(groupStatsTable, featName, stages, microstateNames, figTitle, savePath)

    Z = nan(length(stages), length(microstateNames));

    for st = 1:length(stages)
        for m = 1:length(microstateNames)

            idx = strcmp(groupStatsTable.Feature, featName) & ...
                  strcmp(groupStatsTable.Stage, stages{st}) & ...
                  groupStatsTable.MicrostateID == m;

            if any(idx)
                Z(st,m) = groupStatsTable.Delta_MeditatorMinusControl(idx);
            end
        end
    end

    fig = figure('Color','w','Position',[100 100 850 550]);

    imagesc(Z);
    colorbar;

    xticks(1:length(microstateNames));
    xticklabels(strcat("MS-", microstateNames));

    yticks(1:length(stages));
    yticklabels(stages);

    xlabel('Microstate');
    ylabel('Stage');
    title(figTitle, 'Interpreter','none');

    maxAbs = max(abs(Z(:)), [], 'omitnan');
    if ~isempty(maxAbs) && ~isnan(maxAbs) && maxAbs > 0
        caxis([-maxAbs maxAbs]);
    end

    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            if ~isnan(Z(i,j))
                text(j, i, sprintf('%.3f', Z(i,j)), ...
                    'HorizontalAlignment','center', ...
                    'FontSize',8, ...
                    'Color','k');
            end
        end
    end

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end


function plot_dominance_pattern(dominanceTable, stages, savePath)

    groups = {'control','meditator'};
    Z = nan(length(groups), length(stages));

    for g = 1:length(groups)
        for st = 1:length(stages)

            idx = strcmp(dominanceTable.Group, groups{g}) & ...
                  strcmp(dominanceTable.Stage, stages{st});

            if any(idx)
                Z(g,st) = dominanceTable.DominantMSID(idx);
            end
        end
    end

    fig = figure('Color','w','Position',[100 100 1000 300]);

    imagesc(Z);
    colormap(parula(5));
    colorbar;
    caxis([1 5]);

    xticks(1:length(stages));
    xticklabels(stages);

    yticks(1:length(groups));
    yticklabels(groups);

    title('Study 2 | Dominant microstate by group and stage');

    labels = {'A','B','C','D','E'};

    for i = 1:size(Z,1)
        for j = 1:size(Z,2)
            if ~isnan(Z(i,j))
                text(j, i, sprintf('MS-%s', labels{Z(i,j)}), ...
                    'HorizontalAlignment','center', ...
                    'FontWeight','bold', ...
                    'Color','k');
            end
        end
    end

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end


function plot_stage_profile(T, featName, stageName, microstateNames, figTitle, savePath)

    groups = {'control','meditator'};
    meanMat = nan(length(groups), length(microstateNames));
    semMat = nan(length(groups), length(microstateNames));

    for g = 1:length(groups)
        for m = 1:length(microstateNames)

            idx = strcmp(T.Group, groups{g}) & ...
                  strcmp(T.Stage, stageName) & ...
                  T.MicrostateID == m;

            vals = T.(featName)(idx);
            vals = vals(~isnan(vals));

            meanMat(g,m) = mean(vals, 'omitnan');
            semMat(g,m) = std(vals, 'omitnan') / sqrt(length(vals));
        end
    end

    fig = figure('Color','w','Position',[100 100 800 500]);
    hold on;

    x = 1:length(microstateNames);

    errorbar(x - 0.08, meanMat(1,:), semMat(1,:), '-o', ...
        'LineWidth',1.5, 'DisplayName','control');

    errorbar(x + 0.08, meanMat(2,:), semMat(2,:), '-o', ...
        'LineWidth',1.5, 'DisplayName','meditator');

    xticks(x);
    xticklabels(strcat("MS-", microstateNames));

    ylabel(featName);
    title(figTitle, 'Interpreter','none');
    legend('Location','best');
    grid on;

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end


function plot_composite_stage_profile(T, featName, stageList, microstateNames, figTitle, savePath)

    groups = {'control','meditator'};
    meanMat = nan(length(groups), length(microstateNames));
    semMat = nan(length(groups), length(microstateNames));

    subjectList = unique(T.Subject);

    compositeRows = [];

    for s = 1:length(subjectList)

        subjectName = subjectList{s};

        idxSub = strcmp(T.Subject, subjectName);
        groupName = unique(T.Group(idxSub));

        for m = 1:length(microstateNames)

            idx = strcmp(T.Subject, subjectName) & ...
                  ismember(T.Stage, stageList) & ...
                  T.MicrostateID == m;

            vals = T.(featName)(idx);
            compVal = mean(vals, 'omitnan');

            compositeRows = [compositeRows; ...
                {subjectName, groupName{1}, m, compVal}];
        end
    end

    C = cell2table(compositeRows, ...
        'VariableNames', {'Subject','Group','MicrostateID','Value'});

    for g = 1:length(groups)
        for m = 1:length(microstateNames)

            vals = C.Value(strcmp(C.Group, groups{g}) & C.MicrostateID == m);
            vals = vals(~isnan(vals));

            meanMat(g,m) = mean(vals, 'omitnan');
            semMat(g,m) = std(vals, 'omitnan') / sqrt(length(vals));
        end
    end

    fig = figure('Color','w','Position',[100 100 800 500]);
    hold on;

    x = 1:length(microstateNames);

    errorbar(x - 0.08, meanMat(1,:), semMat(1,:), '-o', ...
        'LineWidth',1.5, 'DisplayName','control');

    errorbar(x + 0.08, meanMat(2,:), semMat(2,:), '-o', ...
        'LineWidth',1.5, 'DisplayName','meditator');

    xticks(x);
    xticklabels(strcat("MS-", microstateNames));

    ylabel(featName);
    title(figTitle, 'Interpreter','none');
    legend('Location','best');
    grid on;

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end