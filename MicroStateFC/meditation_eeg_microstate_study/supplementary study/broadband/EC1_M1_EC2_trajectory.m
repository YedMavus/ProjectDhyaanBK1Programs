%% CODE FOR Broadband case: STUDY 1
clc; clear; close all;

%% ============================================================
% STUDY 1 — BROADBAND MICROSTATE + GLOBAL wPLI ANALYSIS
%
% Cases supported:
%   runMode = 'pilot2'  -> 003S control + 013AR meditator
%   runMode = 'full10'  -> selected 5 controls + 5 meditators
%
% Pipeline:
%   load saved subject templates
%   derive common 5 templates
%   backfit EEG windows using common templates
%   compute:
%       Coverage_A-E
%       Duration_A-E
%       GlobalWPLI_A-E
%   statistics:
%       pilot2: within-subject window permutation tests
%       full10: LMM Feature ~ Stage*Group + (1|Subject)
%   plots:
%       common maps
%       EC1-M1-EC2 trajectories
%       window boxplots
%       change-score plots
%% ============================================================

%% ===================== CONFIG =====================

% runMode = 'pilot2';   
runMode = 'full10';

bandName = 'broadband';

rootPath = '/Users/anatapmitra/Desktop/NSP project/data/Segmented Data';
templateRoot = '/Users/anatapmitra/Desktop/NSP project/MicroStateTemplates_fastK5_v2';
subjectListPath = '/Users/anatapmitra/Desktop/NSP project/BK1AllSubjectList.mat';
chanlocsPath = '/Users/anatapmitra/Desktop/actiCap64_UOL.mat';

outputRoot = fullfile('/Users/anatapmitra/Desktop/NSP project/Study1_Broadband_Results', runMode);

Fs = 1000;
numChannels = 64;

allStages = {'EO1','EC1','G1','M1','G2','EO2','EC2','M2'};
mainStages = {'EC1','M1','EC2','M2'};
trajectoryStages = {'EC1','M1','EC2'};

K = 5;
microstateNames = {'A','B','C','D','E'};

windowLengthSec = 30;
windowLengthSamples = windowLengthSec * Fs;

minSamplesForWPLI = round(1.0 * Fs);  % at least 1 sec samples/state/window

nPerm = 5000;

rng(1234);

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

figDir = fullfile(outputRoot, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

tableDir = fullfile(outputRoot, 'tables');
if ~exist(tableDir, 'dir')
    mkdir(tableDir);
end

modelDir = fullfile(outputRoot, 'models');
if ~exist(modelDir, 'dir')
    mkdir(modelDir);
end

%% ===================== LOAD CHANLOCS =====================

tmpMontage = load(chanlocsPath);
chanlocs = tmpMontage.chanlocs;

%% ===================== SELECT SUBJECTS =====================

switch runMode

    case 'pilot2'
        subjects = {'003S'; '013AR'};
        groupLabels = {'control'; 'meditator'};

    case 'full10'
        S = load(fullfile(templateRoot, 'selected_subjects.mat'));

        if isfield(S, 'subjects') && isfield(S, 'groupLabels')
            subjects = S.subjects(:);
            groupLabels = S.groupLabels(:);
        else
            error('Could not find subjects/groupLabels in selected_subjects.mat');
        end

    otherwise
        error('Invalid runMode');
end

fprintf('\nSubjects used:\n');
disp(table(subjects, groupLabels));

%% ============================================================
% 1. LOAD SUBJECT TEMPLATES AND DERIVE COMMON TEMPLATES
%% ============================================================

fprintf('\nLoading saved subject templates...\n');

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
commonModel.subjects = subjects;
commonModel.groupLabels = groupLabels;
commonModel.templates = commonTemplates;
commonModel.K = K;
commonModel.GEV = commonGEV;

save(fullfile(modelDir, sprintf('commonTemplates_%s_%s.mat', bandName, runMode)), ...
    'commonModel', '-v7.3');

plot_common_templates(commonTemplates, chanlocs, microstateNames, ...
    sprintf('Common %s templates | %s', bandName, runMode), ...
    fullfile(figDir, sprintf('common_templates_%s_%s.png', bandName, runMode)));

%% ============================================================
% 2. BACKFIT + FEATURE EXTRACTION
%% ============================================================

fprintf('\nBackfitting and extracting features...\n');

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

featureCsvPath = fullfile(tableDir, sprintf('features_%s_%s.csv', bandName, runMode));
writetable(featureTable, featureCsvPath);

save(fullfile(tableDir, sprintf('features_%s_%s.mat', bandName, runMode)), ...
    'featureTable', '-v7.3');

fprintf('\nSaved feature table:\n%s\n', featureCsvPath);

%% ============================================================
% 3. STATISTICS
%% ============================================================

fprintf('\nRunning statistics...\n');

features = {'Coverage','Duration','GlobalWPLI'};
comparisons = {
    'EC1','M1';
    'EC1','EC2';
    'M1','M2'
};

if strcmp(runMode, 'pilot2')

    permRows = {};

    for s = 1:length(subjects)

        subjectName = subjects{s};
        groupName = groupLabels{s};

        for f = 1:length(features)

            featName = features{f};

            for m = 1:K

                msName = microstateNames{m};

                for c = 1:size(comparisons,1)

                    st1 = comparisons{c,1};
                    st2 = comparisons{c,2};

                    x = featureTable.(featName)( ...
                        strcmp(featureTable.Subject, subjectName) & ...
                        strcmp(featureTable.Stage, st1) & ...
                        featureTable.MicrostateID == m);

                    y = featureTable.(featName)( ...
                        strcmp(featureTable.Subject, subjectName) & ...
                        strcmp(featureTable.Stage, st2) & ...
                        featureTable.MicrostateID == m);

                    x = x(~isnan(x));
                    y = y(~isnan(y));

                    if length(x) < 2 || length(y) < 2
                        pval = NaN;
                        obsDiff = NaN;
                        effectD = NaN;
                    else
                        [pval, obsDiff, effectD] = permutation_test_equalized(x, y, nPerm);
                    end

                    permRows(end+1,:) = {subjectName, groupName, featName, msName, m, ...
                        st1, st2, obsDiff, effectD, pval, length(x), length(y)};

                end
            end
        end
    end

    permTable = cell2table(permRows, ...
        'VariableNames', {'Subject','Group','Feature','Microstate','MicrostateID', ...
                          'Stage1','Stage2','MeanDiff_Stage2MinusStage1', ...
                          'CohensD','PermP','N1','N2'});

    writetable(permTable, fullfile(tableDir, sprintf('pilot2_permutation_tests_%s.csv', bandName)));

    save(fullfile(tableDir, sprintf('pilot2_permutation_tests_%s.mat', bandName)), ...
        'permTable', '-v7.3');

%% change here 
elseif strcmp(runMode, 'full10')

    lmmRows = {};

    for f = 1:length(features)

        featName = features{f};

        for m = 1:K

            msName = microstateNames{m};

            idx = ismember(featureTable.Stage, mainStages) & ...
                  featureTable.MicrostateID == m;

            tmpT = featureTable(idx, :);

            response = tmpT.(featName);

            valid = ~isnan(response);

            LMMTable = table();

            LMMTable.Response = response(valid);
            LMMTable.Subject  = categorical(tmpT.Subject(valid));
            LMMTable.Group    = categorical(tmpT.Group(valid));
            LMMTable.Stage    = categorical(tmpT.Stage(valid));

            % Optional: force EC1 as baseline stage and control as baseline group
            LMMTable.Stage = reordercats(LMMTable.Stage, {'EC1','M1','EC2','M2'});
            LMMTable.Group = reordercats(LMMTable.Group, {'control','meditator'});

            try

                fprintf('Running LMM: %s | MS%s | N=%d\n', ...
                    featName, msName, height(LMMTable));

                lme = fitlme(LMMTable, ...
                    'Response ~ Stage*Group + (1|Subject)');

                anovaT = anova(lme);

                save(fullfile(modelDir, sprintf('LMM_%s_MS%s_%s.mat', ...
                    featName, msName, bandName)), ...
                    'lme', 'anovaT', 'LMMTable');

                % anova(lme) may return dataset/table depending on MATLAB version
                if istable(anovaT)
                    anovaOut = anovaT;
                else
                    anovaOut = dataset2table(anovaT);
                end

                writetable(anovaOut, fullfile(tableDir, sprintf('LMM_ANOVA_%s_MS%s_%s.csv', ...
                    featName, msName, bandName)));

                termNames = string(anovaOut.Term);

                interactionRows = contains(termNames, 'Stage:Group') | ...
                                  contains(termNames, 'Group:Stage');

                if any(interactionRows)
                    pInteraction = anovaOut.pValue(find(interactionRows, 1));
                else
                    pInteraction = NaN;
                end

                lmmRows(end+1,:) = {featName, msName, m, ...
                    height(LMMTable), pInteraction};

            catch ME

                warning('LMM failed for %s MS%s: %s', ...
                    featName, msName, ME.message);

                lmmRows(end+1,:) = {featName, msName, m, ...
                    height(LMMTable), NaN};

            end
        end
    end

    lmmSummary = cell2table(lmmRows, ...
        'VariableNames', {'Feature','Microstate','MicrostateID','N','StageGroupInteractionP'});

    writetable(lmmSummary, fullfile(tableDir, sprintf('LMM_summary_%s.csv', bandName)));

    save(fullfile(tableDir, sprintf('LMM_summary_%s.mat', bandName)), ...
        'lmmSummary', '-v7.3');

end

%% ============================================================
% 4. CHANGE SCORES
%% ============================================================

fprintf('\nComputing change scores...\n');

changeRows = {};

changeDefs = {
    'M1_minus_EC1', 'M1', 'EC1';
    'EC2_minus_EC1', 'EC2', 'EC1'
};

for s = 1:length(subjects)

    subjectName = subjects{s};
    groupName = groupLabels{s};

    for f = 1:length(features)

        featName = features{f};

        for m = 1:K

            msName = microstateNames{m};

            for cd = 1:size(changeDefs,1)

                changeName = changeDefs{cd,1};
                stageA = changeDefs{cd,2};
                stageB = changeDefs{cd,3};

                a = featureTable.(featName)( ...
                    strcmp(featureTable.Subject, subjectName) & ...
                    strcmp(featureTable.Stage, stageA) & ...
                    featureTable.MicrostateID == m);

                b = featureTable.(featName)( ...
                    strcmp(featureTable.Subject, subjectName) & ...
                    strcmp(featureTable.Stage, stageB) & ...
                    featureTable.MicrostateID == m);

                delta = mean(a, 'omitnan') - mean(b, 'omitnan');

                changeRows(end+1,:) = {subjectName, groupName, featName, msName, m, ...
                    changeName, stageA, stageB, delta};

            end
        end
    end
end

changeTable = cell2table(changeRows, ...
    'VariableNames', {'Subject','Group','Feature','Microstate','MicrostateID', ...
                      'ChangeName','StageA','StageB','Delta'});

writetable(changeTable, fullfile(tableDir, sprintf('change_scores_%s_%s.csv', bandName, runMode)));

save(fullfile(tableDir, sprintf('change_scores_%s_%s.mat', bandName, runMode)), ...
    'changeTable', '-v7.3');

%% ============================================================
% 5. PLOTS
%% ============================================================

fprintf('\nGenerating plots...\n');

for f = 1:length(features)

    featName = features{f};

    for m = 1:K

        msName = microstateNames{m};

        plot_trajectory(featureTable, featName, m, msName, trajectoryStages, ...
            sprintf('%s | MS-%s | EC1-M1-EC2 trajectory', featName, msName), ...
            fullfile(figDir, sprintf('trajectory_%s_MS%s_%s.png', featName, msName, bandName)));

        plot_window_boxplot(featureTable, featName, m, msName, trajectoryStages, ...
            sprintf('%s | MS-%s | window boxplot', featName, msName), ...
            fullfile(figDir, sprintf('boxplot_%s_MS%s_%s.png', featName, msName, bandName)));

        plot_change_scores(changeTable, featName, m, msName, ...
            sprintf('%s | MS-%s | change scores', featName, msName), ...
            fullfile(figDir, sprintf('changescore_%s_MS%s_%s.png', featName, msName, bandName)));

    end
end

fprintf('\nDONE Study 1 broadband analysis.\nResults saved at:\n%s\n', outputRoot);


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


function [pval, obsDiff, effectD] = permutation_test_equalized(x, y, nPerm)

    x = x(:);
    y = y(:);

    n = min(length(x), length(y));

    if length(x) > n
        x = x(randperm(length(x), n));
    end

    if length(y) > n
        y = y(randperm(length(y), n));
    end

    obsDiff = mean(y, 'omitnan') - mean(x, 'omitnan');

    pooledStd = std([x; y], 'omitnan');

    if pooledStd == 0
        effectD = NaN;
    else
        effectD = obsDiff / pooledStd;
    end

    combined = [x; y];
    n1 = length(x);

    permDiffs = nan(nPerm,1);

    for p = 1:nPerm

        idx = randperm(length(combined));

        xPerm = combined(idx(1:n1));
        yPerm = combined(idx(n1+1:end));

        permDiffs(p) = mean(yPerm, 'omitnan') - mean(xPerm, 'omitnan');

    end

    pval = mean(abs(permDiffs) >= abs(obsDiff));

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


function plot_trajectory(T, featName, msID, msName, stages, figTitle, savePath)

    subT = T(T.MicrostateID == msID & ismember(T.Stage, stages), :);

    G = groupsummary(subT, {'Group','Subject','Stage'}, 'mean', featName);

    fig = figure('Color','w','Position',[100 100 800 500]);
    hold on;

    groups = unique(G.Group);

    for g = 1:length(groups)

        groupName = groups{g};
        Tg = G(strcmp(G.Group, groupName), :);
        subjects = unique(Tg.Subject);

        for s = 1:length(subjects)

            Ts = Tg(strcmp(Tg.Subject, subjects{s}), :);

            y = nan(1,length(stages));

            for st = 1:length(stages)
                idx = strcmp(Ts.Stage, stages{st});
                if any(idx)
                    y(st) = Ts.(['mean_' featName])(idx);
                end
            end

            plot(1:length(stages), y, '-o', 'DisplayName', sprintf('%s-%s', groupName, subjects{s}));
        end
    end

    xticks(1:length(stages));
    xticklabels(stages);

    ylabel(featName);
    title(sprintf('%s | MS-%s', figTitle, msName), 'Interpreter','none');
    legend('Location','bestoutside');
    grid on;

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end


function plot_window_boxplot(T, featName, msID, msName, stages, figTitle, savePath)

    subT = T(T.MicrostateID == msID & ismember(T.Stage, stages), :);

    fig = figure('Color','w','Position',[100 100 900 500]);

    boxplot(subT.(featName), strcat(subT.Group, "_", subT.Stage));

    ylabel(featName);
    title(sprintf('%s | MS-%s', figTitle, msName), 'Interpreter','none');
    grid on;

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end


function plot_change_scores(changeTable, featName, msID, msName, figTitle, savePath)

    subT = changeTable(strcmp(changeTable.Feature, featName) & ...
                       changeTable.MicrostateID == msID, :);

    fig = figure('Color','w','Position',[100 100 900 500]);
    hold on;

    cats = unique(subT.ChangeName, 'stable');

    xBase = 1:length(cats);

    for c = 1:length(cats)

        Tc = subT(strcmp(subT.ChangeName, cats{c}), :);

        for r = 1:height(Tc)

            if strcmp(Tc.Group{r}, 'control')
                x = xBase(c) - 0.12;
            else
                x = xBase(c) + 0.12;
            end

            scatter(x, Tc.Delta(r), 80, 'filled');

            text(x + 0.03, Tc.Delta(r), Tc.Subject{r}, ...
                'FontSize', 8, 'Interpreter','none');

        end
    end

    yline(0, '--');

    xticks(xBase);
    xticklabels(cats);
    ylabel(sprintf('Delta %s', featName));

    title(sprintf('%s | MS-%s', figTitle, msName), 'Interpreter','none');
    grid on;

    saveas(fig, savePath);
    savefig(fig, replace(savePath, '.png', '.fig'));
    close(fig);

end