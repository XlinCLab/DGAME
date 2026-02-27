function F_preproc_EEG(subject_ids, subject_dirs, experiment_root, ...
    experiment_outdir, ica_outdir, matlab_root, ...
    dgame_version, removed_electrodes, channels_to_remove_per_subj)

blocks = {'11','12','21','22'};

%% ------------------------------------------------------------------------
% EEGLAB
% -------------------------------------------------------------------------
cd(matlab_root);
addpath('./eeglab2025.1.0');
eeglab;

chanlocs = fullfile(matlab_root, 'eeglab2025.1.0', 'plugins', 'dipfit', 'standard_BESA', 'standard-10-5-cap385.elp');

%% ========================================================================
for s = 1:length(subject_ids)
    clear EEG EEG_raw TMP
    subj = subject_ids{s};
    subject_xdf_dir = string(subject_dirs{s});
    raw_set_before_filtering = [subj,'_raw_before_filtering.set'];
    pre_ica_file = [subj,'_director_preICA.set'];
    post_ica_file = [subj,'_director_postIC.set'];
    cleaned_set = [subj,'_director_cleaned.set'];
    all_ics_set = [subj,'_director_allICs.set'];

    outpath = fullfile(experiment_outdir,'eeg',subj);
    subj_ica_outdir = fullfile(ica_outdir,subj);
    if ~isfolder(outpath); mkdir(outpath); end
    if ~isfolder(subj_ica_outdir); mkdir(subj_ica_outdir); end

    %% =====================================================================
    % ====================== LOAD & MERGE BLOCKS ===========================
    %% =====================================================================

    for b = 1:length(blocks)
        block = string(blocks{b});

        % ----------- LOAD XDF ------------------------------
        xdfFile = fullfile(subject_xdf_dir, 'Director', "dgame" + string(dgame_version) + "_" + subj + "_Director_" + block + ".xdf");
        xdfFile = char(xdfFile);
        streams = load_xdf(xdfFile);

        % --- Find EEG stream ---
        eeg_stream = [];
        for i = 1:length(streams)
            if isfield(streams{i}.info,'type') && strcmpi(streams{i}.info.type,'EEG')
                eeg_stream = streams{i};
                break;
            end
        end

        % --- Extract data ---
        data = double(eeg_stream.time_series);
        srate = str2double(eeg_stream.info.nominal_srate);

        % --- Remove non-EEG channels ---
        labels = cell(1, length(eeg_stream.info.desc.channels.channel));
        for ch = 1:length(eeg_stream.info.desc.channels.channel)
            labels{ch} = eeg_stream.info.desc.channels.channel{ch}.label;
        end
        remove_labels = {'ACC128','ACC129','ACC130','Packet Counter','TRIGGER'};
        keep_idx = ~ismember(labels, remove_labels);
        data = data(keep_idx, :);
        labels = labels(keep_idx);

        % --- Apply MoBILAB scaling ---
        SCALE_FACTOR = 104.1178; % empirically determined 1.041177792474590e+02
        data = data / SCALE_FACTOR;

        % --- Create EEGLAB struct ---
        tmp_EEG = pop_importdata( ...
            'data', data, ...
            'srate', srate, ...
            'nbchan', size(data,1), ...
            'pnts', size(data,2), ...
            'xmin', 0);

        % --- Assign channel labels ---
        for ch = 1:tmp_EEG.nbchan
            tmp_EEG.chanlocs(ch).labels = labels{ch};
        end

        tmp_EEG.ref = 'common';
        tmp_EEG = eeg_checkset(tmp_EEG);
        
        %% -------------------- LOAD WORD EVENTS ---------------------------
        trialtime_filename = sprintf('%s_words2erp_%s_trialtime.csv',subj,block);
        event_file = fullfile(experiment_outdir,'audio',subj,trialtime_filename);
        event_file = char(event_file);

        words = table2struct(readtable(event_file));
        [words.saccAmpl] = deal([]);
        [words.fix_at] = deal([]);
        [words.latency] = deal([]);
        [words.trial_time_char] = words.trial_time;
        [words.trial_time] = deal([]);

        for ev = 1:length(words)
            words(ev).latency = words(ev).time * tmp_EEG.srate;
            words(ev).duration = words(ev).tmax - words(ev).time;
            words(ev).type = words(ev).pos;
            if strcmp(words(ev).trial_time_char, 'NA') == 0
                words(ev).trial_time = str2num(words(ev).trial_time_char);
            end
        end

        %% -------------------- LOAD FIXATIONS ------------------------------
        fix_filename = sprintf('fixations_times_%s_trials.csv',block);
        fix_file = fullfile(experiment_root,'preproc','eyetracking','fixations',subj,fix_filename);
        fix_file = char(fix_file);

        fixations = table2struct(readtable(fix_file));

        fixation_events = [];
        [fixation_events.id] = deal([]);
        [fixation_events.time] = deal([]);
        [fixation_events.text] = deal([]);
        [fixation_events.tmax] = deal([]);
        [fixation_events.freqClass] = deal([]);
        [fixation_events.pattern] = deal([]);
        [fixation_events.set] = deal([]);
        [fixation_events.position] = deal([]);
        [fixation_events.condition] = deal([]);
        [fixation_events.condition_code] = deal([]);
        [fixation_events.pos] = deal([]);
        [fixation_events.surface] = deal([]);
        [fixation_events.surface_end] = deal([]);
        [fixation_events.surface_competitor] = deal([]);
        [fixation_events.targetA_surface] = deal([]);
        [fixation_events.targetB_surface] = deal([]);
        [fixation_events.fillerA_surface] = deal([]);
        [fixation_events.fillerB_surface] = deal([]);
        [fixation_events.compA_surface] = deal([]);
        [fixation_events.compB_surface] = deal([]);
        [fixation_events.target_location] = deal([]);
        [fixation_events.p_target] = deal([]);
        [fixation_events.p_target_pre] = deal([]);
        [fixation_events.trial] = deal([]);
        [fixation_events.trial_time] = deal([]);
        [fixation_events.saccAmpl] = deal([]);
        [fixation_events.fix_at] = deal([]);
        [fixation_events.latency] = deal([]);
        [fixation_events.duration] = deal([]);
        [fixation_events.type] = deal([]);
        [fixation_events.trial_time] = deal([]);
        [fixation_events.trial_time_char] = deal([]);
       
        for ev = 1:length(fixations)
            fixation_events(ev).type = 'fixation';
            fixation_events(ev).latency = fixations(ev).time * 500;
            fixation_events(ev).duration = fixations(ev).duration;
            fixation_events(ev).fix_at = fixations(ev).fix_at;
            fixation_events(ev).trial_time = fixations(ev).trial_time;
            if fixations(ev).saccAmpl > 0
                fixation_events(ev).saccAmpl = fixations(ev).saccAmpl;
            end
        end

        % Merge word and fixation events
        tmp_EEG.event = [words; fixation_events'];
        tmp_EEG = eeg_checkset(tmp_EEG,'eventconsistency');

        %% -------------------- RESAMPLE BEFORE MERGE ----------------------
        tmp_EEG = pop_resample(tmp_EEG,250);
        tmp_EEG = pop_chanedit(tmp_EEG,'lookup',chanlocs);

        if exist('EEG','var')
            EEG = pop_mergeset(EEG,tmp_EEG);
            EEG = eeg_checkset(EEG);
        else
            EEG = tmp_EEG;
        end
    end
    % Interim save
    pop_saveset(EEG,'filename',raw_set_before_filtering,'filepath',outpath);

    EEG_raw = EEG;

    %% =====================================================================
    % ============================ CLEANING ================================
    %% =====================================================================
    % Pre-clean data
    EEG = pop_eegfiltnew(EEG,'locutoff',2);
    % Exclude channels that were removed to fit ET glasses
    EEG = pop_select(EEG,'nochannel',removed_electrodes);
    % Exclude any other specified channels, e.g. due to broken electrodes or excessive noise
    channels_to_remove = channels_to_remove_per_subj{s};
    if ~isempty(channels_to_remove)
        fprintf('Subject %s: Removing channels: %s\n', subj, ...
            strjoin(cellfun(@num2str, channels_to_remove, 'UniformOutput', false), ', '));
        EEG = pop_select(EEG,'nochannel',channels_to_remove);
    end

    EEG = pop_clean_rawdata(EEG,...
        'FlatlineCriterion',5,...
        'ChannelCriterion',0.8,...
        'LineNoiseCriterion',4,...
        'Highpass','off',...
        'BurstCriterion','off',...
        'WindowCriterion','off',...
        'BurstRejection','off',...
        'Distance','Euclidian');
    EEG = pop_rejchan(EEG,...
        'elec',[1:EEG.nbchan],...
        'threshold',2,...
        'norm','on',...
        'measure','kurt');

    pools = gcp('nocreate');
    cpus = feature('numCores');
    if size(pools) == 0
        pool = parpool(cpus);
    end
        
    % Low-pass filter at 100 Hz
    EEG = pop_eegfiltnew(EEG,...
        'hicutoff',100,...
        'plotfreqz',0);

    % Run cleanline to remove line noise
    EEG = pop_cleanline(EEG,...
        'bandwidth',2,...
        'linefreqs',50,...
        'chanlist',[1:EEG.nbchan],...
        'computepower',1,...
        'linefreqs',[50],...
        'newversion',1,...
        'normSpectrum',0,...
        'p',0.01,...
        'pad',2,...
        'plotfigures',0,...
        'scanforlines',0,...
        'sigtype',...
        'Channels',...
        'taperbandwidth',2,...
        'tau',100,...
        'verb',1,...
        'winsize',4,...
        'winstep',1);
    EEG = pop_chanedit(EEG, 'lookup',chanlocs);
    EEG = pop_clean_rawdata(EEG,...
        'FlatlineCriterion','off',...
        'ChannelCriterion','off',...
        'LineNoiseCriterion','off',...
        'Highpass','off',...
        'BurstCriterion',10,...
        'WindowCriterion','off',...
        'BurstRejection','off',...
        'Distance','Euclidian');

    %% ==================== INTERPOLATE MISSING CHANNELS ===========================
    all_chans = {EEG_raw.chanlocs.labels}';
    kept_chans = {EEG.chanlocs.labels}';
    bad_chans = setdiff(all_chans,kept_chans);
    EEG.bad_chans = bad_chans;
    EEG = pop_interp(EEG,EEG_raw.chanlocs,'spherical');

    %% ===================== AVERAGE REFERENCE =============================
    % Re-reference to average, adding and again removing the initial reference
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select(EEG,'nochannel',{'initialReference'});
    EEG = pop_select(EEG,'nochannel',bad_chans);
    EEG = eeg_checkset(EEG);

    % Save interim pre-ICA file
    EEG = pop_saveset(EEG,'filename',pre_ica_file,'filepath',outpath);

    %% =====================================================================
    % =============================== ICA ==================================
    %% =====================================================================
    % Define parameters for amica
    dataRank = sum(eig(cov(double(EEG.data'))) > 1E-6);
    numprocs= 1;        % number of nodes
    max_threads= 10;    % number of threads
    num_models= 1;      % number of models of mixture ICA
    max_iter= 20;       % max number of learning steps % run amica #20 for testing, set to 2000 at least
    num_rej = 10;       % # of rejections of unlikely data
    EEG = pop_resample(EEG,100);
    % Run the actual decomposition
    runamica15(EEG.data,...
        'num_chans',EEG.nbchan,...
        'outdir',subj_ica_outdir,...
        'pcakeep',dataRank,...
        'num_models', num_models,...
        'do_reject',1,...
        'numrej', num_rej,...
        'rejsig', 2.5,...
        'rejint', 1,...
        'max_iter',max_iter,...
        'numprocs',numprocs,...
        'max_threads',max_threads);

    mod = loadmodout15(subj_ica_outdir);
    EEG.icaweights = mod.W;
    EEG.icasphere  = mod.S;
    EEG = eeg_checkset(EEG, 'ica');
    EEG = pop_saveset(EEG, 'filename',post_ica_file,'filepath',outpath);

    % Clean data using ICs
    % Write ICA info to temporary structure to transfer to uncleaned data
    TMP.icasphere = EEG.icasphere;
    TMP.icaweights = EEG.icaweights;
    EEG_raw.bad_chans = bad_chans;

    % Use "raw" data again for different filtering etc.
    EEG = EEG_raw;
    EEG = pop_select(EEG,'nochannel',EEG.bad_chans);
    EEG = pop_eegfiltnew(EEG, 0.3,[]);% 'ftype', 'highpass', 'wtype', 'hamming', 'forder', 826, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
    EEG = pop_eegfiltnew(EEG,[],20);
    EEG = pop_chanedit(EEG, 'lookup',chanlocs);

    % ASR again
    EEG = pop_clean_rawdata(EEG,...
        'FlatlineCriterion','off',...
        'ChannelCriterion','off',...
        'LineNoiseCriterion','off',...
        'Highpass','off',...
        'BurstCriterion',10,...
        'WindowCriterion','off',...
        'BurstRejection','off',...
        'Distance','Euclidian');
    EEG = pop_chanedit(EEG, 'lookup',chanlocs);

    % Interpolate channels
    EEG = pop_interp(EEG,EEG_raw.chanlocs);
    EEG = eeg_checkset(EEG);
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(EEG.nbchan).labels = 'initialReference';

    % Average reference
    EEG = pop_reref(EEG, []);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
    EEG = pop_select(EEG,'nochannel',bad_chans);
    EEG.badchans = bad_chans;

    % Transfer ICA weights and sphere from TMP structure    
    EEG.icaweights = TMP.icaweights;
    EEG.icasphere = TMP.icasphere;
    EEG = eeg_checkset(EEG,'ica');
    
    %% =====================================================================
    % =========================== IC LABEL =================================
    %% =====================================================================
    % Classify components using IC label
    EEG = pop_iclabel(EEG, 'default');

    % Mark brain and non-brain components in a vector of booleans
    classifications = EEG.etc.ic_classification.ICLabel.classifications;
    [classes ~] = find(transpose(EEG.etc.ic_classification.ICLabel.classifications == max(EEG.etc.ic_classification.ICLabel.classifications,[],2)));
    ICs_keep = zeros(length(classes),1);
    summed_scores_to_keep = zeros(size(EEG.icawinv,2),1);
    classes_to_keep = 1;
    for class_to_keep = classes_to_keep
        ICs_keep(classes == class_to_keep) = 1;
        summed_scores_to_keep = summed_scores_to_keep + EEG.etc.ic_classification.ICLabel.classifications(:,1);
    end
    EEG.ICs_throw = find(~ICs_keep);
    EEG.ICs_keep = find(ICs_keep);

    % Interim save
    EEG = pop_saveset(EEG,'filename',all_ics_set,'filepath',outpath);

    % Remove non-brain components
    EEG = pop_subcomp(EEG,[EEG.ICs_throw]);

    %% ================= FINAL INTERPOLATION ===============================
    % Interpolate missing channels again
    EEG = pop_interp(EEG,EEG_raw.chanlocs);

    % Save set
    pop_saveset(EEG,'filename',cleaned_set,'filepath',outpath);

end
end