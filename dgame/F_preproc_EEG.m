function F_preproc_EEG(subject_dirs, experiment_root, matlab_root)

% Get script directory
this_file = mfilename('fullpath');
script_dir = fileparts(this_file);

% Setup Python environment
python_path = fullfile(script_dir, '..', 'venv', 'bin', 'python');
pyenv('Version', python_path);

% Add script directory to Python sys.path
if count(py.sys.path, script_dir) == 0
    insert(py.sys.path, int32(0), script_dir);
end

% Load constants from dgame.constants.py
% (convert Python set -> Python list -> MATLAB cell)
blocks = cell(py.list(py.constants.BLOCK_IDS));


% Mount dependencies / toolboxes
cd(matlab_root);
addpath('./eeglab2021.1');
addpath(genpath(['./eeglab2021.1/plugins/cleanline2.00/external/bcilab_partial/']));
addpath('./mobilab');
eeglab;

ica_outdir = fullfile(experiment_root, 'preproc/icatmp/');
ica_outdir = char(ica_outdir);
if ~exist(ica_outdir, 'dir')
    mkdir(ica_outdir);
end
chanlocs = './eeglab2021.1/plugins/dipfit4.3/standard_BESA/standard-10-5-cap385.elp';


subject_ids = subject_dirs.keys;
for s = 1:length(subject_ids)
    clear EEG
    T = 0;
    subj = subject_ids{s};
    subject_xdf_dir = subject_dirs(subj);
    subject_xdf_dir = string(subject_xdf_dir{1});
    raw_set_before_filtering = [subj,'_raw_before_filtering.set'];
    pre_ica_file = [subj,'_director_preICA.set'];
    post_ica_file = [subj,'_director_postIC.set'];
    cleaned_set = [subj,'_director_cleaned.set'];
    all_ics_set = [subj,'_director_allICs.set'];

    

    outfile = [subj,'_director_allBlocks.set'];
    
    bad_chans_out = [subject_xdf_dir, '_bad_chans.csv'];

    %% load data
    for b = 1:length(blocks)
        % Convert Python integer to MATLAB double and then to MATLAB string
        block = string(double(blocks{b}));
        event_file = [experiment_root,'preproc/audio/',subj,'/',subj,'_words2erp_',block,'_trialtime.csv'];
        fixations_file = [experiment_root,'fixations/',subj,'_fixations_times_',block,'_trials.csv'];
        xdfFile = fullfile(subject_xdf_dir, 'Director', "dgame2_" + subj + "_Director_" + block + ".xdf");
        xdfFile = char(xdfFile);
        mobipath = fullfile(subject_xdf_dir, 'Director', "dgame2_" + subj + "_Director_" + block + "_MoBI/");
        mobipath = char(mobipath);
        infile = [subj,'_director2.set'];
     
            %load the data
            if ~isfolder(mobipath)
                mobilab.allStreams = dataSourceXDF(xdfFile,mobipath);
            else
                mobilab.allStreams = dataSourceMoBI(mobipath);
            end
            exportIndex = mobilab.allStreams.getItemIndexFromItemClass('eeg');
            tmp = mobilab.allStreams.export2eeglab([exportIndex]);

        %% read word data

        tmp.event = table2struct(readtable(event_file));

        %add fields inside fixations but not words to have matching field names
        [tmp.event.saccAmpl] = deal([]);
        [tmp.event.fix_at] = deal([]);
        [tmp.event.latency] = deal([]);
        [tmp.event.trial_time_char] = tmp.event.trial_time;
        [tmp.event.trial_time] = deal([]);

        for ev=1:length(tmp.event)
            tmp.event(ev).latency = tmp.event(ev).time*tmp.srate;
            tmp.event(ev).duration = tmp.event(ev).tmax-tmp.event(ev).time;
            tmp.event(ev).type = tmp.event(ev).pos;
            if strcmp(tmp.event(ev).trial_time_char, 'NA') == 0
                tmp.event(ev).trial_time = str2num(tmp.event(ev).trial_time_char);
            end
        end
        
        %% read fixation data
        fixations = table2struct(readtable(fixations_file));
        for fix = 1:length(fixations)
            fixations(fix).type = 'fixation';
            fixations(fix).latency = fixations(fix).time * 500;
        end
            
        %create fields so word and fixation events have the same fields
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
            fixation_events(ev).latency = fixations(ev).latency;
            fixation_events(ev).duration = fixations(ev).duration;
            fixation_events(ev).type = fixations(ev).type;
            fixation_events(ev).fix_at = fixations(ev).fix_at;
            fixation_events(ev).trial_time = fixations(ev).trial_time;
            if fixations(ev).saccAmpl > 0
                fixation_events(ev).saccAmpl = fixations(ev).saccAmpl;
            end
        end
        
        tmp.event = [tmp.event;fixation_events'];
        tmp = eeg_checkset(tmp,'eventconsistency');
            
%% resample and merge
        tmp = pop_resample(tmp, 250);
% rename and exclude some channels for some subjects, because we had broken electrodes
        if strcmp('20',subj) == 1 | strcmp('21',subj) | strcmp('22',subj) == 1 | strcmp('23',subj) == 1|strcmp('24',subj) == 1|strcmp('25',subj) == 1 
           if strcmp('22',subj) == 1
              tmp=pop_chanedit(tmp, 'changefield',{118,'labels','FC6'},'changefield',{120,'labels','FC4'});
              elseif strcmp('23',subj) == 1|strcmp('24',subj) == 1 | strcmp('25',subj) == 1 
              tmp=pop_chanedit(tmp, 'changefield',{118,'labels','FC4'},'changefield',{120,'labels','FC6'});
           end
           tmp = pop_select(tmp,'nochannel',{'FC6','FC4'});
        end
        tmp=pop_chanedit(tmp, 'lookup',chanlocs);
        if exist('EEG','var')
            EEG = pop_mergeset(EEG,tmp);
            EEG = eeg_checkset(EEG);
        else
            EEG = tmp;
        end
    end
%interims save
    pop_saveset(EEG,'filename',raw_set_before_filtering,'filepath',outpath); 
end

    EEG_raw = EEG;
%% pre clean data
    EEG = pop_eegfiltnew(EEG, 'locutoff',2);
%exlude the channels that were removed to fit the ET glasses
    EEG = pop_select(EEG,'nochannel',{'FTT10h','FTT9h','FFT10h','FFT9h'});

%remove some more manually selected noisy channels for specific subjects.
    if strcmp('10',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'P10','TPP10h'});
    elseif strcmp('02',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'F10','F9','AF7'});
    elseif strcmp('01',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'TP7','TPP9h','TTP7h'});
    elseif strcmp('17',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'F10','T7','T8'});
    elseif strcmp('19',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'T7','FTT7h'});
    elseif strcmp('24',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'Iz','F10','T7'});
    elseif strcmp('07',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'F9'});
    elseif strcmp('27',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'Iz','PO8'});
    elseif strcmp('07',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'T7','T8','FTT8h','TTP7h'});
    elseif strcmp('08',subj) == 1
        EEG = pop_select(EEG,'nochannel',{'FC6'});
    end
%use clean_raw_data
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
    EEG = pop_rejchan(EEG, 'elec',[1:EEG.nbchan] ,'threshold',2,'norm','on','measure','kurt');
    
    pools = gcp('nocreate');
    cpus = feature('numCores');
    if size(pools) == 0
        pool = parpool(cpus);
    end
%filter at 100 Hz low-pass
    EEG = pop_eegfiltnew(EEG, 'hicutoff',100,'plotfreqz',0);
%run cleanline to remove line noise
    EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan] ,'computepower',1,'linefreqs',[50] ,'newversion',1,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);
    EEG=pop_chanedit(EEG, 'lookup',chanlocs);
    burstCriterion = 10;
    burstRejection = 'off';
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',burstCriterion,'WindowCriterion','off','BurstRejection',burstRejection,'Distance','Euclidian');

%%find and interpolate missing channels using the EEG_raw copy from above
    all_chans = {EEG_raw.chanlocs.labels}';
    same_chans = intersect({EEG.chanlocs.labels}',{EEG_raw.chanlocs.labels}');
    bad_chans=setdiff(all_chans,same_chans);
    EEG.bad_chans = bad_chans;
    %interpolate   
    EEG = pop_interp(EEG,EEG_raw.chanlocs,'spherical');
%% re-reference to average, adding and again removing the initial Reference    
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
    EEG = pop_select(EEG,'nochannel',bad_chans);
    EEG = eeg_checkset( EEG );

%save interims file
    EEG = pop_saveset(EEG,'filename',pre_ica_file,'filepath',outpath);

%% run amica
    % Run ICA
    % define parameters
    dataRank = sum(eig(cov(double(EEG.data'))) > 1E-6); 
    numprocs= 1;       % # of nodes
    max_threads= 10;   % # of threads
    num_models= 1;     % # of models of mixture ICA
    max_iter= 20;    % max number of learning steps % run amica #20 for testing, set to 2000 at least
    num_rej = 10;      % # of rejections of unlikely data
    ica_outname = [subj,'_ica/'];

    EEG = pop_resample(EEG,100);
% run the actual decomposition
    runamica15(EEG.data,'num_chans', EEG.nbchan,'outdir',ica_outdir,'pcakeep',dataRank,'num_models', num_models,'do_reject', 1,'numrej', num_rej,'rejsig', 2.5,'rejint', 1,'max_iter',max_iter,'numprocs',numprocs,'max_threads',max_threads);
    EEG.etc.amica  = loadmodout15(ica_outdir);
    EEG.icaweights = EEG.etc.amica.W;
    EEG.icasphere  = EEG.etc.amica.S;
    EEG = eeg_checkset(EEG, 'ica');
    EEG = pop_saveset(EEG, 'filename',post_ica_file,'filepath',outpath);
%% clean data using ICs
    all_chans = {EEG_raw.chanlocs.labels}';
    same_chans = intersect({EEG.chanlocs.labels}',{EEG_raw.chanlocs.labels}');
    bad_chans=setdiff(all_chans,same_chans);
    %write ica info to temporary structure to transfer to uncleaned data
    TMP.icasphere = EEG.icasphere;
    TMP.icaweights = EEG.icaweights;
    EEG_raw.bad_chans = bad_chans;   
%use "raw" data again for different filtering etc.
    EEG = EEG_raw;
    EEG = pop_select(EEG,'nochannel',EEG.bad_chans);
    EEG = pop_eegfiltnew(EEG, 0.3,[]);% 'ftype', 'highpass', 'wtype', 'hamming', 'forder', 826, 'minphase', 0, 'usefftfilt', 0, 'plotfresp', 0, 'causal', 0);
    EEG = pop_eegfiltnew(EEG,[],20);
    EEG=pop_chanedit(EEG, 'lookup',chanlocs);
% asr again
    burstCriterion = 10;
    burstRejection = 'off';
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off','BurstCriterion',burstCriterion,'WindowCriterion','off','BurstRejection',burstRejection,'Distance','Euclidian');
    EEG=pop_chanedit(EEG, 'lookup',chanlocs);
%interpolate channels
    EEG = pop_interp(EEG,EEG_raw.chanlocs);
    EEG = eeg_checkset(EEG);
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(EEG.nbchan).labels = 'initialReference';
%average reference
    EEG = pop_reref(EEG, []);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
    EEG = pop_select(EEG,'nochannel',bad_chans);
    EEG.badchans = bad_chans;
%transfer ICA weights and sphere from TMP        
    EEG.icaweights = TMP.icaweights;
    EEG.icasphere = TMP.icasphere;
    EEG = eeg_checkset(EEG,'ica');
% classify components using iclabel
    EEG = pop_iclabel(EEG, 'default');
%mark brain and non brain components in a vector of booleans
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
%interims save
    EEG = pop_saveset(EEG,'filename',all_ics_set,'filepath',outpath);
%remove non-brain components
    EEG = pop_subcomp(EEG,[EEG.ICs_throw]);
%interpolate missing channels again
    EEG = pop_interp(EEG,EEG_raw.chanlocs);
%save set
    EEG = pop_saveset(EEG,'filename',cleaned_set,'filepath',outpath);
end  
