%% optimized_glm_erp_reconstruction.m
% Script to reconstruct the ERP as predicted by the GLM
% Input: output of deconvolution_ERP.m (ufresult struct)
% Output: per-channel CSVs of predicted ERPs for "noun" and "fixation" events
%takes 10-15 min per single subject, approx. its parallelized, so depending
%on your N of workers, it can basically finish in that time.

clear all; close all;
%--- Parallel setup ---
pools = gcp('nocreate');
cpus = feature('numCores');
if isempty(pools)
    parpool(cpus);
end

eeglab nogui;

%--- Parameters ---
% Use current folder as study root
studyroot = pwd;
subjects  = {'02'};        % cell array of subject IDs

parfor s = 1:numel(subjects)
    subj     = subjects{s};
    eegpath  = fullfile(studyroot,'eeg');
    inpath   = fullfile(eegpath,'unfold_out');
    outpathS = fullfile(inpath);
    if ~exist(outpathS,'dir'), mkdir(outpathS); end

    %--- Load data ---
    unfold_file = fullfile(inpath, sprintf('%s_ufresult.mat', subj));
    D = load(unfold_file);
    U = D.ufresult;
    EEG = pop_loadset('filename',sprintf('%s_director_cleaned.set',subj), 'filepath',eegpath);

    %--- Precompute stats ---
    mean_tg      = mean([EEG.event(strcmp({EEG.event.type},'N')).trial_time], 'omitnan'); % mean target fixation time of nouns
    mean_saccAmp = mean([EEG.event(strcmp({EEG.event.type},'fixation')).saccAmpl], 'omitnan'); % mean saccadic amplitude of fixations

    %--- Common vectors ---
    ntime = numel(U.times);
    time  = U.times' * 1000;

    % Noun condition times & means
    con_idx     = strcmp({EEG.event.condition},'conflict') & strcmp({EEG.event.type},'N');
    ncon_idx    = strcmp({EEG.event.condition},'no_conflict') & strcmp({EEG.event.type},'N');
    %trial target gaze per condition
    tg_con      = unique([EEG.event(con_idx).trial_time]);  tg_con(isnan(tg_con)) = [];
    tg_nocon    = unique([EEG.event(ncon_idx).trial_time]); tg_nocon(isnan(tg_nocon)) = [];
    %mean trial number per condition
    trial_con   = mean([EEG.event(con_idx).trial], 'omitnan');
    trial_nocon = mean([EEG.event(ncon_idx).trial], 'omitnan');

    %--- Noun ERPs ---
    for ch = 1:numel(U.chanlocs)
        chanlabel = U.chanlocs(ch).labels;
        betaMat   = squeeze(U.beta(ch,:,:));   % [ntime x nParams]
        nParams   = size(betaMat,2);

        % Combine noun trials
        allTimesN = [tg_con, tg_nocon];
        isConN    = [true(1,numel(tg_con)), false(1,numel(tg_nocon))];
        nTrialsN  = numel(allTimesN);

        % Build contrast matrix
        Cn = zeros(nParams, nTrialsN);
        for t = 1:nTrialsN
            v = zeros(nParams,1);
            if isConN(t)
                v(end-4) = 1;   % conflict regressor
                v(end-2) = trial_con;
                v(end-1) = allTimesN(t);
                v(end)   = allTimesN(t);
            else
                v(end-4) = 1;   % always on
                v(end-3) = 0;   % no_conflict off
                v(end-2) = trial_nocon;
                v(end-1) = allTimesN(t);
                v(end)   = 0;
            end
            Cn(:,t) = v;
        end
        ERPn = betaMat * Cn;

        % Preallocate outputs
        nRows = ntime * nTrialsN;
        T_time = zeros(nRows,1);
        T_event= strings(nRows,1);
        T_cond = strings(nRows,1);
        T_mtf  = zeros(nRows,1);
        T_chan = strings(nRows,1);
        T_data = zeros(nRows,1);
        T_subj = strings(nRows,1);

        % Fill
        idx0 = 1;
        for t = 1:nTrialsN
            rows = idx0:(idx0+ntime-1);
            T_time(rows)  = time;
            T_event(rows) = "noun";
            if isConN(t)
                T_cond(rows) = "conflict";
            else
                T_cond(rows) = "no_conflict";
            end
            T_mtf(rows)   = allTimesN(t);
            T_chan(rows)  = chanlabel;
            T_data(rows)  = ERPn(:,t);
            T_subj(rows)  = subj;
            idx0 = idx0 + ntime;
        end

        nounTable = table(T_time, T_event, T_cond, T_mtf, T_chan, T_data, T_subj, ...
            'VariableNames',{'time','event','condition','mean_target_fixation','channel','data','subject'});
        writetable(nounTable, fullfile(outpathS,sprintf('%s_%s_unfold_N.csv',subj,chanlabel)), 'QuoteStrings', true);
    end

    %--- Fixation ERPs ---
    fixTimes  = unique(round([EEG.event(strcmp({EEG.event.type},'fixation')).trial_time],1));
    fixLabels = ["elsewhere","other","target"];
    conVals   = [true,false];
    nF        = numel(fixTimes)*numel(conVals)*numel(fixLabels);

    for ch = 1:numel(U.chanlocs)
        chanlabel = U.chanlocs(ch).labels;
        betaMat   = squeeze(U.beta(ch,:,:));
        nParams   = size(betaMat,2);

        % Build fixation contrast matrix & metadata
        Cfix = zeros(nParams, nF);
        condM = strings(nF,1);
        fixM  = strings(nF,1);
        timeM = zeros(nF,1);
        idxC = 1;
        for ft = fixTimes
            for cv = conVals
                for fl = fixLabels
                    v = zeros(nParams,1);
                    v(3) = 1;         % 'n' term
                    v(4) = cv;       % conflict
                    if fl=="other"
                        v(5) = 1;
                    elseif fl=="target"
                        v(6) = 1;
                    end
                    v(7) = ft;
                    v(8) = cv*(fl=="other");
                    v(9) = cv*(fl=="target");
                    v(10)= cv*ft;
                    v(11)= (fl=="other")*ft;
                    v(12)= (fl=="target")*ft;
                    v(13)= cv*((fl=="other")*ft);
                    v(14)= cv*((fl=="target")*ft);
                    v(15:18)=mean_saccAmp;
                    Cfix(:,idxC) = v;

                    % metadata
                    if cv
                        condM(idxC) = "conflict";
                    else
                        condM(idxC) = "no_conflict";
                    end
                    fixM(idxC)  = fl;
                    timeM(idxC)= ft;
                    idxC = idxC + 1;
                end
            end
        end
        ERPf = betaMat * Cfix;

        % Preallocate
        nRows = ntime * nF;
        T_time   = zeros(nRows,1);
        T_event  = strings(nRows,1);
        T_cond   = strings(nRows,1);
        T_fixat  = strings(nRows,1);
        T_trtime = zeros(nRows,1);
        T_chan   = strings(nRows,1);
        T_data   = zeros(nRows,1);
        T_subj   = strings(nRows,1);

        % Fill
        idx0 = 1;
        for t = 1:nF
            rows = idx0:(idx0+ntime-1);
            T_time(rows)   = time;
            T_event(rows)  = "fixation";
            T_cond(rows)   = condM(t);
            T_fixat(rows)  = fixM(t);
            T_trtime(rows) = timeM(t);
            T_chan(rows)   = chanlabel;
            T_data(rows)   = ERPf(:,t);
            T_subj(rows)   = subj;
            idx0 = idx0 + ntime;
        end

        fixTable = table(T_time, T_event, T_cond, T_fixat, T_trtime, T_chan, T_data, T_subj, ...
            'VariableNames',{'time','event','condition','fix_at','trial_time','channel','data','subject'});
        writetable(fixTable, fullfile(outpathS,sprintf('%s_%s_unfold_FIX.csv',subj,chanlabel)), 'QuoteStrings', true);
    end
end

