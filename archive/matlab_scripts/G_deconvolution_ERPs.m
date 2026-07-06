function G_deconvolution_ERPs(subject_ids, subject_dirs, matlab_root)

% Mount dependencies / toolboxes
cd(matlab_root);
addpath('./eeglab2025.1.0');
eeglab;

%% main loop over subjects
for s = 1:length(subject_ids)
    subj = subject_ids{s};
    subject_dir = subject_dirs{s};
    cleaned_set = [subj,'_director_cleaned.set'];
    cleaned_set_out = [subj,'_director_cleaned_w_unfold.set'];
    outpath = fullfile(subject_dir, 'unfold_out');
    unfold_matfile = fullfile(outpath, sprintf('%s_ufresult.mat', subj));
    unfold_outfile_beta = [outpath,subj,'_beta_dc.csv'];
    EEG = pop_loadset('filename',cleaned_set,'filepath',subject_dir);
%% setup the fixations for analysis
%backup event
[EEG.event_old] = EEG.event;


%loop through events and find fixations to add condition information retrieved from the entry of the noun of the corresponding trial
    for ev = 1:length(EEG.event)
        if strcmp(EEG.event(ev).type,'fixation')==1 % if its a fixation
        %summarize non-target fixations
            if strcmp(EEG.event(ev).type, 'fixation') == 1 & strcmp(EEG.event(ev).fix_at,'target') == 0 & strcmp(EEG.event(ev).fix_at,'other') == 0
               EEG.event(ev).fix_at = 'elsewhere';
            end
            if EEG.event(ev).trial_time < 0 % if it has trial time value smaller 0, i.e. occurs before noun onset
                if EEG.event(ev).trial_time >= -3.5 % if the trial time is greater or equal -3.5 relative to noun onset
                    for find_idx = 1:200 %go forward a few events and see if you find a noun, stop at the first one (200 is arbitrary so we are sure to get all nouns in all trials (e.g. with many many many fixations in a trial)
                         if find_idx + ev < length(EEG.event) % for events close to beginning/end of recording
                            if strcmp(EEG.event(ev+find_idx).type,'N') == 1 
                                EEG.event(ev).condition = EEG.event(ev+find_idx).condition; 
                                EEG.event(ev).set = EEG.event(ev+find_idx).set;
                                EEG.event(ev).pattern = EEG.event(ev+find_idx).pattern;
                                break;
                            end
                        end
                    end
                end
            %do the same for fixations after noun onset
            elseif EEG.event(ev).trial_time > 0
                if EEG.event(ev).trial_time <=3.5
                    for find_idx = 1:200 %go back a few events and see if you find a noun, stop at the first one (200 is arbitrary so we are sure to get all nouns in all trials (e.g. with many many many fixations in a trial)
                        if ev-find_idx > 0
                            if strcmp(EEG.event(ev-find_idx).type,'N') == 1 
                                EEG.event(ev).condition = EEG.event(ev-find_idx).condition;
                                EEG.event(ev).set = EEG.event(ev-find_idx).set;
                                EEG.event(ev).pattern = EEG.event(ev-find_idx).pattern;
                                break;
                            end
                        end
                    end
                end
            end
            %set all fixation without a condition to a different type        
            if strcmp(EEG.event(ev).type,'fixation')==1 & strcmp(EEG.event(ev).condition, 'conflict') == 0 & strcmp(EEG.event(ev).condition, 'no_conflict') == 0
                EEG.event(ev).type = ['other_fixation'];
                EEG.event(ev).condition_old = EEG.event(ev).condition;
                EEG.event(ev).condition = ['NA'];
            end 
        end
    end
    
    for ev = 1:length(EEG.event)
        if strcmp(EEG.event(ev).type, 'fixation') == 1 & strcmp(EEG.event(ev).fix_at,'target') == 0 & strcmp(EEG.event(ev).fix_at,'other') == 0
            EEG.event(ev).fix_at = 'elsewhere';
       
        end

    end
    
%do deconvolution
    init_unfold
    cfgDesign = [];
    cfgDesign.eventtypes = {{'prev'},{'next'},{'fixation'},{'D'},{'N'}};
    cfgDesign.codingschema = 'reference';
    % We use intercept-only formula for "other" words (=all words but the critical ones) because we are only interested in the overlap for now
    cfgDesign.formula = {'y ~ 1','y ~ 1','y ~ 1+cat(condition)*cat(fix_at)*trial_time+spl(saccAmpl,5)','y ~ 1+cat(condition)+trial','y ~ 1+cat(condition)*trial_time+trial'};

    cfgDesign.categorical = {'condition',{'no_conflict','conflict'}}; % level no-conflict is reference
    EEG = uf_designmat(EEG,cfgDesign);
    EEG = uf_imputeMissing(EEG,'method','drop');
    cfgTimeexpand = [];
    cfgTimeexpand.timelimits = [-0.5,1.5];
    EEG = uf_timeexpandDesignmat(EEG,cfgTimeexpand);
%detect and exclude artifacts (150 uV step size)
    winrej = uf_continuousArtifactDetect(EEG,'amplitudeThreshold',150);
    EEG = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));
    EEG= uf_glmfit(EEG,'method','lsmr','channel',[1:EEG.nbchan]);%,'ica',1,'channel',1:length(EEG.icaact(:,2)));
    
    ufresult= uf_condense(EEG);
    %optional: plot paramters
    %uf_plotParam(ufresult,'channel',18);
    save(unfold_matfile, 'ufresult');
    uftable = uf_unfold2csv(ufresult,'filename',unfold_outfile_beta,'deconv',1);
   pop_saveset(EEG,'filename',cleaned_set_out,'filepath',outpath);

end
