clc; clear; close all;

%% =========================================================
% Full microstate + FC feature extraction pipeline
%
% FINAL REDUCED VERSION:
%   1) Load raw EEG/LFP elec*.mat files
%   2) Compute full-signal spectral features using mtmfft + fourier:
%        PSD, COH, PPC, PLI, WPLI
%   3) For each band (alpha/gamma):
%        - load saved microstate templates
%        - recreate same band-limited signal used for extraction
%        - backfit templates to every sample
%        - smooth microstate labels
%        - compute temporal microstate features
%        - compute Hilbert-based PLI/WPLI:
%             full-signal and microstate-specific
%   4) Save outputs
%
% NO mtmconvol / NO time-frequency FC / NO microstate-specific COH/PPC
%
% Requires FieldTrip on MATLAB path.
%% =========================================================

%% -----------------------------
% USER SETTINGS
%% -----------------------------
cfgMain = struct();

% Input data root
cfgMain.rootPath = '/Users/anatapmitra/Desktop/NSP project/data/Segmented Data';

% Output root
cfgMain.outputBase = '/Users/anatapmitra/Desktop/NSP project/extracted_features';

% Subject / date / stages
cfgMain.subjectName = '013AR';
cfgMain.expDate     = '280122';
cfgMain.protocols   = {'M1', 'G2'};

% Bands for microstate backfitting / Hilbert FC
cfgMain.bands(1).name  = 'alpha';
cfgMain.bands(1).range = [8 13];
cfgMain.bands(2).name  = 'gamma';
cfgMain.bands(2).range = [30 80];

% Microstate smoothing
cfgMain.minSegmentMs = 30;

% Hilbert computation
cfgMain.hilbertDownsampleFactor = 10;   % increase to 1 for denser estimation

% Full-signal mtmfft settings
cfgMain.mtmfft.foilim    = [0 100];
cfgMain.mtmfft.taper     = 'dpss';
cfgMain.mtmfft.tapsmofrq = 4;          % +-4 Hz smoothing
cfgMain.mtmfft.pad       = 'nextpow2';

% Save raw backfitted labels
cfgMain.saveRawBackfit = true;

%% -----------------------------
% MAIN LOOP OVER STAGES
%% -----------------------------
for p = 1:numel(cfgMain.protocols)
    protocolName = cfgMain.protocols{p};

    fprintf('\n========================================================\n');
    fprintf('Processing subject=%s | stage=%s\n', cfgMain.subjectName, protocolName);
    fprintf('========================================================\n');

    folderSegment = fullfile(cfgMain.rootPath, cfgMain.subjectName, ...
        'EEG', cfgMain.expDate, protocolName, 'segmentedData', 'LFP');

    if ~exist(folderSegment, 'dir')
        warning('Folder not found: %s. Skipping.', folderSegment);
        continue;
    end

    %% -----------------------------
    % Load raw stage data
    %% -----------------------------
    [eegDataAll, electrodes, channelLabels, Fs, numTrials, numSamples, lfpInfo] = ...
        load_stage_lfp(folderSegment);

    numElectrodes = numel(electrodes);
    totalSamples  = numTrials * numSamples;
    totalDuration = totalSamples / Fs;

    fprintf('Loaded stage data: C=%d, trials=%d, samples/trial=%d, Fs=%g Hz\n', ...
        numElectrodes, numTrials, numSamples, Fs);

    eegContinuousRaw = reshape(eegDataAll, numElectrodes, totalSamples);

    % Build FieldTrip raw structure from raw stage data
    data_ft = build_ft_raw(eegDataAll, channelLabels, Fs);

    % Channel pair definitions
    [pairI, pairJ, labelcmb] = make_channel_pairs(channelLabels);

    %% -----------------------------
    % Create stage output folder
    %% -----------------------------
    outStageFolder = fullfile(cfgMain.outputBase, cfgMain.subjectName, protocolName);
    if ~exist(outStageFolder, 'dir')
        mkdir(outStageFolder);
    end

    %% =========================================================
    % A) FULL-SIGNAL SPECTRAL FEATURES (mtmfft + fourier)
    %% =========================================================
    fprintf('Computing full-signal mtmfft spectral features...\n');

    spectral_full = compute_full_signal_spectral_features(data_ft, pairI, pairJ, cfgMain.mtmfft);

    % stage-level metadata
    metadata_stage = struct();
    metadata_stage.subjectName   = cfgMain.subjectName;
    metadata_stage.expDate       = cfgMain.expDate;
    metadata_stage.protocolName  = protocolName;
    metadata_stage.numElectrodes = numElectrodes;
    metadata_stage.numTrials     = numTrials;
    metadata_stage.numSamples    = numSamples;
    metadata_stage.Fs            = Fs;
    metadata_stage.totalDuration = totalDuration;
    metadata_stage.electrodes    = electrodes(:);
    metadata_stage.channelLabels = channelLabels(:);
    metadata_stage.labelcmb      = labelcmb;
    metadata_stage.pairI         = pairI;
    metadata_stage.pairJ         = pairJ;

    save(fullfile(outStageFolder, 'full_signal_mtmfft_spectral_fc.mat'), ...
        'metadata_stage', 'spectral_full', '-v7.3');

    %% =========================================================
    % B) LOOP OVER BANDS (alpha / gamma)
    %% =========================================================
    for b = 1:numel(cfgMain.bands)
        bandName  = cfgMain.bands(b).name;
        bandRange = cfgMain.bands(b).range;

        fprintf('\n---- Band: %s [%d %d] Hz ----\n', bandName, bandRange(1), bandRange(2));

        % Load saved microstate file
        microstateFile = fullfile(folderSegment, ...
            sprintf('%s_%s_Microstates_%s.mat', cfgMain.subjectName, protocolName, bandName));

        if ~exist(microstateFile, 'file')
            warning('Microstate file not found: %s. Skipping band.', microstateFile);
            continue;
        end

        Sms = load(microstateFile);
        microstateMaps = Sms.microstateMaps;   % M x C
        M = size(microstateMaps, 1);

        % Create output folder for this band
        outBandFolder = fullfile(outStageFolder, bandName);
        if ~exist(outBandFolder, 'dir')
            mkdir(outBandFolder);
        end

        %% =========================================================
        % 1) RECREATE BAND-LIMITED SIGNAL
        %% =========================================================
        eegContinuousBand = fft_bandpass_continuous(eegContinuousRaw, Fs, bandRange(1), bandRange(2));
        eegBandTrials     = reshape(eegContinuousBand, numElectrodes, numTrials, numSamples);

        %% =========================================================
        % 2) BACKFIT MICROSTATE LABELS
        %% =========================================================
        fprintf('Backfitting microstates for %s...\n', bandName);

        z_cont_raw  = backfit_microstates_continuous(eegContinuousBand, microstateMaps);
        z_trial_raw = reshape(z_cont_raw, numTrials, numSamples);

        % Smooth if needed
        minSegSamples  = max(1, round(cfgMain.minSegmentMs * Fs / 1000));
        z_trial_smooth = smooth_microstate_labels(z_trial_raw, minSegSamples);
        z_cont_smooth  = reshape(z_trial_smooth, 1, []);

        %% =========================================================
        % 3) MICROSTATE FEATURES
        %% =========================================================
        fprintf('Computing microstate features...\n');
        ms_features = compute_microstate_features(z_trial_smooth, eegBandTrials, microstateMaps, Fs);

        %% =========================================================
        % 4) HILBERT-BASED PLI/WPLI
        %% =========================================================
        fprintf('Computing Hilbert-based PLI/WPLI...\n');
        hilbert_fc = compute_hilbert_fc( ...
            eegBandTrials, z_trial_smooth, pairI, pairJ, M, cfgMain.hilbertDownsampleFactor);

        %% =========================================================
        % 5) SAVE BAND OUTPUTS
        %% =========================================================
        metadata_band = metadata_stage;
        metadata_band.bandName  = bandName;
        metadata_band.bandRange = bandRange;

        % ---- backfit + microstate features ----
        backfitFile = fullfile(outBandFolder, 'microstate_backfit_and_features.mat');

        if cfgMain.saveRawBackfit
            z_trial_backfit_raw = z_trial_raw;
        end
        z_trial_backfit_smooth = z_trial_smooth;
        z_cont_backfit_smooth  = z_cont_smooth;

        save(backfitFile, ...
            'metadata_band', ...
            'microstateMaps', ...
            'ms_features', ...
            'z_trial_backfit_smooth', ...
            'z_cont_backfit_smooth', ...
            '-v7.3');

        if cfgMain.saveRawBackfit
            save(backfitFile, 'z_trial_backfit_raw', '-append');
        end

        % ---- Hilbert FC ----
        hilbertFile = fullfile(outBandFolder, 'fc_hilbert_bandlimited.mat');
        save(hilbertFile, ...
            'metadata_band', ...
            'microstateMaps', ...
            'hilbert_fc', ...
            '-v7.3');

        fprintf('Saved band outputs to: %s\n', outBandFolder);
    end
end

fprintf('\nAll done.\n');

%% =========================================================
% LOCAL FUNCTIONS
%% =========================================================

function [eegDataAll, electrodes, channelLabels, Fs, numTrials, numSamples, lfpInfo] = load_stage_lfp(folderSegment)

    lfpInfo = load(fullfile(folderSegment, 'lfpInfo.mat'));

    electrodes = lfpInfo.analogChannelsStored(:)';
    numElectrodes = numel(electrodes);

    % infer Fs if possible
    if isfield(lfpInfo, 'Fs')
        Fs = lfpInfo.Fs;
    elseif isfield(lfpInfo, 'fs')
        Fs = lfpInfo.fs;
    elseif isfield(lfpInfo, 'samplingRate')
        Fs = lfpInfo.samplingRate;
    else
        Fs = 1000;
    end

    % infer labels if possible
    if isfield(lfpInfo, 'analogChannelLabels') && numel(lfpInfo.analogChannelLabels) == numElectrodes
        channelLabels = lfpInfo.analogChannelLabels(:);
    elseif isfield(lfpInfo, 'channelLabels') && numel(lfpInfo.channelLabels) == numElectrodes
        channelLabels = lfpInfo.channelLabels(:);
    else
        channelLabels = arrayfun(@(x) sprintf('elec%d', x), electrodes, 'UniformOutput', false)';
    end

    % infer trial/sample counts from first file
    firstFound = false;
    for i = 1:numElectrodes
        elecFile = fullfile(folderSegment, sprintf('elec%d.mat', electrodes(i)));
        if exist(elecFile, 'file')
            tmp = load(elecFile, 'analogData');
            [numTrials, numSamples] = size(tmp.analogData);
            firstFound = true;
            break;
        end
    end
    assert(firstFound, 'No elec*.mat files found in %s', folderSegment);

    eegDataAll = zeros(numElectrodes, numTrials, numSamples, 'double');

    fprintf('Loading %d electrode files...\n', numElectrodes);
    for i = 1:numElectrodes
        elecFile = fullfile(folderSegment, sprintf('elec%d.mat', electrodes(i)));
        if ~exist(elecFile, 'file')
            error('Missing file: %s', elecFile);
        end
        tmp = load(elecFile);
        assert(isfield(tmp, 'analogData'), 'analogData missing in %s', elecFile);
        assert(isequal(size(tmp.analogData), [numTrials, numSamples]), ...
            'Unexpected analogData size in %s', elecFile);
        eegDataAll(i,:,:) = tmp.analogData;
    end
end

function data_ft = build_ft_raw(eegDataAll, channelLabels, Fs)
    % eegDataAll: C x T x S
    [~, T, S] = size(eegDataAll);

    data_ft = [];
    data_ft.label   = channelLabels(:);
    data_ft.fsample = Fs;
    data_ft.trial   = cell(1, T);
    data_ft.time    = cell(1, T);
    data_ft.sampleinfo = zeros(T, 2);

    tvec = (0:S-1) / Fs;
    for tr = 1:T
        data_ft.trial{tr} = squeeze(eegDataAll(:, tr, :)); % C x S
        data_ft.time{tr}  = tvec;
        begsamp = (tr-1)*S + 1;
        endsamp = tr*S;
        data_ft.sampleinfo(tr,:) = [begsamp endsamp];
    end
end

function [pairI, pairJ, labelcmb] = make_channel_pairs(channelLabels)
    C = numel(channelLabels);
    pairs = nchoosek(1:C, 2);
    pairI = pairs(:,1);
    pairJ = pairs(:,2);

    labelcmb = cell(size(pairs,1), 2);
    for p = 1:size(pairs,1)
        labelcmb{p,1} = channelLabels{pairI(p)};
        labelcmb{p,2} = channelLabels{pairJ(p)};
    end
end

function eegBandCont = fft_bandpass_continuous(eegContinuousRaw, Fs, lowFreq, highFreq)
    % Same filtering style as original microstate extraction
    [C, L] = size(eegContinuousRaw);
    eegBandCont = zeros(C, L);

    freqs = (0:L-1) * (Fs/L);
    filterMask = (freqs >= lowFreq & freqs <= highFreq) | ...
                 (freqs >= (Fs - highFreq) & freqs <= (Fs - lowFreq));

    for ch = 1:C
        signal = eegContinuousRaw(ch,:);
        signal_fft = fft(signal);
        eegBandCont(ch,:) = real(ifft(signal_fft .* filterMask));
    end
end

function z_cont = backfit_microstates_continuous(eegContinuousBand, microstateMaps)
    % eegContinuousBand: C x N
    % microstateMaps   : M x C
    [~, N] = size(eegContinuousBand);
    M = size(microstateMaps, 1);

    z_cont = zeros(1, N);

    % normalize templates
    A = microstateMaps';
    A = A ./ max(vecnorm(A, 2, 1), eps);

    for n = 1:N
        x = eegContinuousBand(:,n);
        nx = norm(x);
        if nx < eps
            z_cont(n) = 1;
            continue;
        end
        x = x / nx;
        sims = abs(A' * x);   % polarity-invariant
        [~, z_cont(n)] = max(sims);
    end
end

function z_trial_smooth = smooth_microstate_labels(z_trial_raw, minSegSamples)
    [T, ~] = size(z_trial_raw);
    z_trial_smooth = z_trial_raw;

    for tr = 1:T
        z = z_trial_raw(tr,:);
        changed = true;
        iter = 0;

        while changed && iter < 20
            iter = iter + 1;
            changed = false;

            [runStart, runEnd, runLabel] = run_length_encoding(z);
            runLen = runEnd - runStart + 1;

            for r = 1:numel(runLabel)
                if runLen(r) < minSegSamples
                    changed = true;

                    leftLabel  = [];
                    rightLabel = [];

                    if r > 1
                        leftLabel = runLabel(r-1);
                    end
                    if r < numel(runLabel)
                        rightLabel = runLabel(r+1);
                    end

                    if isempty(leftLabel) && isempty(rightLabel)
                        % do nothing
                    elseif isempty(leftLabel)
                        z(runStart(r):runEnd(r)) = rightLabel;
                    elseif isempty(rightLabel)
                        z(runStart(r):runEnd(r)) = leftLabel;
                    else
                        leftLen  = runLen(r-1);
                        rightLen = runLen(r+1);
                        if leftLen >= rightLen
                            z(runStart(r):runEnd(r)) = leftLabel;
                        else
                            z(runStart(r):runEnd(r)) = rightLabel;
                        end
                    end
                end
            end
        end

        z_trial_smooth(tr,:) = z;
    end
end

function [runStart, runEnd, runLabel] = run_length_encoding(z)
    d = [true, diff(z) ~= 0];
    runStart = find(d);
    runEnd   = [runStart(2:end)-1, numel(z)];
    runLabel = z(runStart);
end

function ms = compute_microstate_features(z_trial, eegBandTrials, microstateMaps, Fs)
    % z_trial       : T x S
    % eegBandTrials : C x T x S
    % microstateMaps: M x C

    [T, S] = size(z_trial);
    M = size(microstateMaps, 1);

    totalSamples = T * S;
    totalTimeSec = totalSamples / Fs;

    Coverage_m = zeros(M,1);
    Duration_m = nan(M,1);
    Occurrence_m = zeros(M,1);
    DwellStd_m = nan(M,1);
    numSegments_m = zeros(M,1);
    dwellSecondsCell = cell(M,1);

    TransitionCount = zeros(M, M);

    for tr = 1:T
        z = z_trial(tr,:);

        for m = 1:M
            Coverage_m(m) = Coverage_m(m) + sum(z == m);
        end

        [runStart, runEnd, runLabel] = run_length_encoding(z);
        runLenSec = (runEnd - runStart + 1) / Fs;

        for r = 1:numel(runLabel)
            m = runLabel(r);
            dwellSecondsCell{m}(end+1,1) = runLenSec(r); %#ok<AGROW>
            numSegments_m(m) = numSegments_m(m) + 1;
        end

        for n = 1:(S-1)
            i = z(n);
            j = z(n+1);
            TransitionCount(i,j) = TransitionCount(i,j) + 1;
        end
    end

    Coverage_m = Coverage_m / totalSamples;
    Occurrence_m = numSegments_m / totalTimeSec;

    for m = 1:M
        if ~isempty(dwellSecondsCell{m})
            Duration_m(m) = mean(dwellSecondsCell{m});
            DwellStd_m(m) = std(dwellSecondsCell{m});
        end
    end

    TransitionProb = zeros(M, M);
    for i = 1:M
        rowSum = sum(TransitionCount(i,:));
        if rowSum > 0
            TransitionProb(i,:) = TransitionCount(i,:) / rowSum;
        end
    end

    p_state = Coverage_m(:)';
    p_state = p_state(p_state > 0);
    LabelEntropy = -sum(p_state .* log2(p_state));

    rowMass = sum(TransitionCount, 2);
    rowMass = rowMass / max(sum(rowMass), eps);
    TransitionEntropy = 0;
    for i = 1:M
        pi_row = TransitionProb(i,:);
        pi_row = pi_row(pi_row > 0);
        if ~isempty(pi_row)
            TransitionEntropy = TransitionEntropy + rowMass(i) * (-sum(pi_row .* log2(pi_row)));
        end
    end

    % GEV
    Xcont = reshape(eegBandTrials, size(eegBandTrials,1), []);
    zcont = reshape(z_trial, 1, []);
    gfp = std(Xcont, 1, 1);

    A = microstateMaps';
    A = A ./ max(vecnorm(A, 2, 1), eps);

    gev_num = zeros(M,1);
    gev_den = sum(gfp.^2);

    for n = 1:size(Xcont,2)
        m = zcont(n);
        x = Xcont(:,n);
        nx = norm(x);
        if nx < eps
            continue;
        end
        x = x / nx;
        corrVal = A(:,m)' * x;
        gev_num(m) = gev_num(m) + (gfp(n)^2) * (corrVal^2);
    end

    GEV_m = gev_num / max(gev_den, eps);
    TotalGEV = sum(GEV_m);

    ms = struct();
    ms.Duration_m        = Duration_m;
    ms.Coverage_m        = Coverage_m;
    ms.Occurrence_m      = Occurrence_m;
    ms.GEV_m             = GEV_m;
    ms.TotalGEV          = TotalGEV;
    ms.DwellStd_m        = DwellStd_m;
    ms.NumSegments_m     = numSegments_m;
    ms.TransitionCount   = TransitionCount;
    ms.TransitionProb    = TransitionProb;
    ms.LabelEntropy      = LabelEntropy;
    ms.TransitionEntropy = TransitionEntropy;
end

function hilbert_fc = compute_hilbert_fc(eegBandTrials, z_trial, pairI, pairJ, M, dsFactor)
    % eegBandTrials: C x T x S
    % z_trial      : T x S
    [C, T, S] = size(eegBandTrials);

    obsAll = [];
    zAll   = [];

    for tr = 1:T
        X = squeeze(eegBandTrials(:, tr, :));   % C x S
        U = hilbert(X.').';                     % C x S analytic

        idx = 1:dsFactor:S;
        obsAll = [obsAll; U(:, idx).']; %#ok<AGROW>  % Nobs x C
        zAll   = [zAll; z_trial(tr, idx).']; %#ok<AGROW>
    end

    [PLI_full, WPLI_full] = compute_hilbert_pair_metrics(obsAll, pairI, pairJ, C);

    PLI_ms  = nan(C, C, M);
    WPLI_ms = nan(C, C, M);

    for m = 1:M
        mask = (zAll == m);
        [PLI_ms(:,:,m), WPLI_ms(:,:,m)] = compute_hilbert_pair_metrics(obsAll(mask,:), pairI, pairJ, C);
    end

    hilbert_fc = struct();
    hilbert_fc.PLI_full  = PLI_full;
    hilbert_fc.WPLI_full = WPLI_full;
    hilbert_fc.PLI_ms    = PLI_ms;
    hilbert_fc.WPLI_ms   = WPLI_ms;

    hilbert_fc.nodeStrength_PLI_full  = compute_node_strength(PLI_full);
    hilbert_fc.nodeStrength_WPLI_full = compute_node_strength(WPLI_full);

    hilbert_fc.nodeStrength_PLI_ms  = zeros(C, M);
    hilbert_fc.nodeStrength_WPLI_ms = zeros(C, M);
    for m = 1:M
        hilbert_fc.nodeStrength_PLI_ms(:,m)  = compute_node_strength(PLI_ms(:,:,m));
        hilbert_fc.nodeStrength_WPLI_ms(:,m) = compute_node_strength(WPLI_ms(:,:,m));
    end
end

function [PLI_mat, WPLI_mat] = compute_hilbert_pair_metrics(obsAll, pairI, pairJ, C)
    P = numel(pairI);

    PLI_vec  = nan(P,1);
    WPLI_vec = nan(P,1);

    if isempty(obsAll)
        PLI_mat  = nan(C, C);
        WPLI_mat = nan(C, C);
        return;
    end

    chunkSize = 64;
    for st = 1:chunkSize:P
        en = min(st + chunkSize - 1, P);
        idx = st:en;

        crossObs = obsAll(:, pairI(idx)) .* conj(obsAll(:, pairJ(idx)));  % Nobs x chunk
        imObs    = imag(crossObs);

        PLI_vec(idx) = abs(mean(sign(imObs), 1)).';

        den = sum(abs(imObs), 1);
        num = abs(sum(imObs, 1));
        tmp = nan(1, numel(idx));
        valid = den > 0;
        tmp(valid) = num(valid) ./ den(valid);
        WPLI_vec(idx) = tmp(:);
    end

    PLI_mat  = pairvec_to_symmetric(PLI_vec,  pairI, pairJ, C);
    WPLI_mat = pairvec_to_symmetric(WPLI_vec, pairI, pairJ, C);
end

function spectral_full = compute_full_signal_spectral_features(data_ft, pairI, pairJ, mtmCfg)
    % Full-signal spectral features using mtmfft + fourier
    % Returns:
    %   powspctrm_full : C x F
    %   COH_full       : C x C x F
    %   PPC_full       : C x C x F
    %   PLI_full       : C x C x F
    %   WPLI_full      : C x C x F
    %   node strengths for all

    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.output      = 'fourier';
    cfg.taper       = mtmCfg.taper;
    cfg.tapsmofrq   = mtmCfg.tapsmofrq;
    cfg.foilim      = mtmCfg.foilim;
    cfg.keeptrials  = 'yes';
    cfg.keeptapers  = 'yes';
    cfg.pad         = mtmCfg.pad;
    cfg.channel     = 'all';

    freq = ft_freqanalysis(cfg, data_ft);

    % expected dimord: rpttap_chan_freq
    Fcoef = freq.fourierspctrm;
    freqVec = freq.freq(:)';

    [Nobs, C, Fnum] = size(Fcoef); %#ok<ASGLU>
    P = numel(pairI);

    % Full PSD
    powspctrm_full = squeeze(mean(abs(Fcoef).^2, 1));  % C x F

    % Pairwise metrics
    COH_full  = nan(C, C, Fnum);
    PPC_full  = nan(C, C, Fnum);
    PLI_full  = nan(C, C, Fnum);
    WPLI_full = nan(C, C, Fnum);

    chunkSize = 64;

    for fi = 1:Fnum
        Xf = squeeze(Fcoef(:,:,fi));   % Nobs x C
        pwr = squeeze(powspctrm_full(:,fi)); %#ok<NASGU>

        coh_vec  = nan(P,1);
        ppc_vec  = nan(P,1);
        pli_vec  = nan(P,1);
        wpli_vec = nan(P,1);

        for st = 1:chunkSize:P
            en = min(st + chunkSize - 1, P);
            idx = st:en;

            Xi = Xf(:, pairI(idx));    % Nobs x chunk
            Xj = Xf(:, pairJ(idx));    % Nobs x chunk

            Sxy_obs = Xi .* conj(Xj);  % Nobs x chunk
            Sxx_obs = abs(Xi).^2;
            Syy_obs = abs(Xj).^2;

            % Coherence
            Sxx_mean = mean(Sxx_obs, 1);
            Syy_mean = mean(Syy_obs, 1);
            Sxy_mean = mean(Sxy_obs, 1);
            coh_vec(idx) = (abs(Sxy_mean) ./ sqrt(max(Sxx_mean .* Syy_mean, eps))).';

            % PLI
            imSxy = imag(Sxy_obs);
            pli_vec(idx) = abs(mean(sign(imSxy), 1)).';

            % WPLI
            num = abs(sum(imSxy, 1));
            den = sum(abs(imSxy), 1);
            tmp = nan(1, numel(idx));
            valid = den > 0;
            tmp(valid) = num(valid) ./ den(valid);
            wpli_vec(idx) = tmp(:);

            % PPC
            magSxy = abs(Sxy_obs);
            unit_phase = zeros(size(Sxy_obs));
            valid_phase = magSxy > 0;
            unit_phase(valid_phase) = Sxy_obs(valid_phase) ./ magSxy(valid_phase);

            Nvalid = sum(valid_phase, 1);
            ptmp = nan(1, numel(idx));
            ok = Nvalid >= 2;
            s = sum(unit_phase(:,ok), 1);
            ptmp(ok) = (abs(s).^2 - Nvalid(ok)) ./ (Nvalid(ok) .* (Nvalid(ok)-1));
            ppc_vec(idx) = ptmp(:);
        end

        COH_full(:,:,fi)  = pairvec_to_symmetric(coh_vec,  pairI, pairJ, C);
        PPC_full(:,:,fi)  = pairvec_to_symmetric(ppc_vec,  pairI, pairJ, C);
        PLI_full(:,:,fi)  = pairvec_to_symmetric(pli_vec,  pairI, pairJ, C);
        WPLI_full(:,:,fi) = pairvec_to_symmetric(wpli_vec, pairI, pairJ, C);
    end

    spectral_full = struct();
    spectral_full.freqVec        = freqVec;
    spectral_full.powspctrm_full = powspctrm_full;   % C x F
    spectral_full.COH_full       = COH_full;         % C x C x F
    spectral_full.PPC_full       = PPC_full;         % C x C x F
    spectral_full.PLI_full       = PLI_full;         % C x C x F
    spectral_full.WPLI_full      = WPLI_full;        % C x C x F

    spectral_full.nodeStrength_COH_full  = compute_node_strength_3d(COH_full);   % C x F
    spectral_full.nodeStrength_PPC_full  = compute_node_strength_3d(PPC_full);
    spectral_full.nodeStrength_PLI_full  = compute_node_strength_3d(PLI_full);
    spectral_full.nodeStrength_WPLI_full = compute_node_strength_3d(WPLI_full);
end

function mat = pairvec_to_symmetric(vec, pairI, pairJ, C)
    mat = nan(C, C);
    mat(1:C+1:end) = 0;
    for p = 1:numel(vec)
        i = pairI(p);
        j = pairJ(p);
        mat(i,j) = vec(p);
        mat(j,i) = vec(p);
    end
end

function ns = compute_node_strength(mat)
    mat = mat - diag(diag(mat));
    ns = nansum(mat, 2);
end

function ns3 = compute_node_strength_3d(mat3)
    [C, ~, Fnum] = size(mat3);
    ns3 = nan(C, Fnum);
    for fi = 1:Fnum
        A = mat3(:,:,fi);
        A = A - diag(diag(A));
        ns3(:,fi) = nansum(A, 2);
    end
end
