

basePath   = "/Volumes/Suvam's 1TB/Supratim_Ray_Files/Supratim Ray's files - segmentedData/";
outputBase = "/Users/iacv/Desktop/NSP_2026/";

experiments = {'G1','G2','M1','M2'};

subjects = dir(basePath);
subjects = subjects([subjects.isdir] & ~startsWith({subjects.name}, '.'));

for s = 1:length(subjects)
    
    subjName = subjects(s).name;
    fprintf('Processing Subject: %s\n', subjName);
    
    for e = 1:length(experiments)
        
        expName = experiments{e};
        
        lfpPath = fullfile(basePath, subjName, ...
            'EEG','221021',expName,'segmentedData','LFP');
        
        if ~exist(lfpPath,'dir')
            continue
        end
        
        % Create output folder
        outPath = fullfile(outputBase, subjName, expName);
        if ~exist(outPath,'dir')
            mkdir(outPath);
        end
        
        for elec = 1:70
            
            fileName = fullfile(lfpPath, sprintf('elec%d.mat', elec));
            
            if ~exist(fileName,'file')
                continue
            end
            
            S = load(fileName,'analogData');

            fprintf('Loaded')
            % 120x2500 -> 1x300000 (preserve time order)
            singleSignal = reshape(S.analogData.',1,[]);
            
            save(fullfile(outPath, sprintf('elec%d_1D.mat', elec)), ...
                 'singleSignal','-v7.3');
        end
        
        fprintf('   Done Experiment: %s\n', expName);
        
    end
end

fprintf('All done.\n');