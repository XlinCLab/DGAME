function A_export_audio_and_et_times(subject_ids, subject_dirs, experiment_root, times_outdir, matlab_root, dgame_version)

blocks = {'11','12','21','22'};

% Mount dependencies / toolboxes
cd(matlab_root);
addpath('./eeglab2021.1');
eeglab;

for s = 1:length(subject_ids)
    subject = subject_ids{s};
    subject_xdf_dir = subject_dirs{s};
    subject_xdf_dir = string(subject_xdf_dir);
    
    for b = 1:length(blocks)
        block = string(blocks{b});
        inpath = fullfile(subject_xdf_dir, 'Director');
        xdfFile = fullfile(inpath,  "dgame" + string(dgame_version) + "_" + subject + "_Director_" + block + ".xdf");
        xdfFile = char(xdfFile);
        mobipath = fullfile(inpath, "dgame" + string(dgame_version) + "_" + subject + "_Director_" + block + "_MoBI");
        mobipath = char(mobipath);
        outpath_audio = fullfile(experiment_root, 'recordings/audio/', subject);  % TODO consider writing to experiment outdir instead
        outpath_times = fullfile(times_outdir, subject);
        director_outfile = fullfile(outpath_audio, subject + "_director_" + block + ".wav");
        decke_outfile = fullfile(outpath_audio, subject + "_decke_" + block + ".wav");
        if ~exist(outpath_times, 'dir')
            mkdir(outpath_times);
        end
        if ~exist(outpath_audio, 'dir')
            mkdir outpath;
        end

        % Read and save the first and last timestamp of ET recording
        if ~exist(xdfFile, 'file')
            error('xdfFile does not exist: %s', xdfFile);
        end
        tmpXDF = load_xdf(xdfFile,'HandleClockSynchronization',false);
        for items = 1:length(tmpXDF)
            if strcmp(tmpXDF{items}.info.name, 'pupil_capture') == 1
                times = [str2num(tmpXDF{items}.info.first_timestamp);str2num(tmpXDF{items}.info.last_timestamp)];
            end
        end
        dlmwrite(char(fullfile(outpath_times, subject + "_timestamps_max-min_" + block + ".csv")), times, 'precision', '%.7f');
        times = [];
        tmpXDF = [];

        % Load data with mobilab
        % Export ET times
        [ET] = pop_loadxdf(xdfFile, 'streamname', 'pupil_capture');
        
        times = ET.times/1000;
        dlmwrite(char(fullfile(outpath_times, subject + "_times_" + block + ".csv")), times, 'precision', '%.6f');
        times = [];

        % Extract audio, normalize and export to wav (normalization is needed or the signal will get clipped = garbage)
        [audio] = pop_loadxdf(xdfFile, 'streamname', 'audio');

        % Normalize audio
        audio_decke = audio.data(1,:)/max(abs(audio.data(1,:)));
        audio_director = audio.data(2,:)/max(abs(audio.data(2,:)));


        % Export to WAV
        audiowrite(director_outfile,audio_director,22050);
        audiowrite(decke_outfile,audio_decke,22050);

    end
end
