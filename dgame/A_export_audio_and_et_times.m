function A_export_audio_and_et_times(subject_dirs, experiment_root, matlab_root)

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

% Load constants from dgame.constants.py\
% (convert Python set -> Python list -> MATLAB cell)
blocks = cell(py.list(py.constants.BLOCK_IDS));

% Mount dependencies / toolboxes
cd(matlab_root);
addpath('./eeglab2021.1');
addpath('./mobilab');
eeglab;
chanlocs = './eeglab2021.1/plugins/dipfit4.3/standard_BESA/standard-10-5-cap385.elp';

subject_ids = subject_dirs.keys;
for s = 1:length(subject_ids)
    subject = subject_ids{s};
    subject_xdf_dir = subject_dirs(subject);
    subject_xdf_dir = string(subject_xdf_dir{1});
    
    for b = 1:length(blocks)
        % Convert Python integer to MATLAB double and then to MATLAB string
        block = string(double(blocks{b}));
        inpath = fullfile(subject_xdf_dir, 'Director');
        xdfFile = fullfile(inpath,  "dgame2_" + subject + "_Director_" + block + ".xdf");
        xdfFile = char(xdfFile);
        mobipath = fullfile(inpath, "dgame2_" + subject + "_Director_" + block + "_MoBI");
        mobipath = char(mobipath);
        outpath = fullfile(experiment_root, 'recordings/audios/', subject);
        outpath_times = fullfile(experiment_root, 'preproc/helper_files/');
        director_outfile = fullfile(outpath, subject + "_director_" + block + ".wav");
        decke_outfile = fullfile(outpath, subject + "_decke_" + block + ".wav");
        if ~isfolder(outpath)
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
        if ~isfolder(mobipath)
            mobilab.allStreams = dataSourceXDF(xdfFile,mobipath);
        else
            mobilab.allStreams = dataSourceMoBI(mobipath);
        end
        exportIndex = mobilab.allStreams.getItemIndexFromItemClass('eeg');
        indexET = mobilab.allStreams.getItemIndexFromItemName('pupil_capture_ieeg-rec-lap-2');
        audioIndex = mobilab.allStreams.getItemIndexFromItemName('audio_xlinc-recording');

        % Export ET times
        ET = mobilab.allStreams.export2eeglab([indexET]);
        
        times = ET.times/1000;
        dlmwrite(char(fullfile(outpath_times, subject + "_times_" + block + ".csv")), times, 'precision', '%.6f');
        times = [];

        % Extract audio, normalize and export to wav (normalization is needed or the signal will get clipped = garbage)
        audio = mobilab.allStreams.export2eeglab([audioIndex]);

        % Normalize audio
        audio_decke = audio.data(1,:)/max(abs(audio.data(1,:)));
        audio_director = audio.data(2,:)/max(abs(audio.data(2,:)));


        % Export to WAV
        audiowrite(director_outfile,audio_director,22050);
        audiowrite(decke_outfile,audio_decke,22050);

    end
end
