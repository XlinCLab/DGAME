clear all;
eeglab;
studyroot = './';
xdfpath = [studyroot,'recordings/xdf/'];

chanlocs = './eeglab2021.1/plugins/dipfit4.3/standard_BESA/standard-10-5-cap385.elp';

subjects = {'02'};
blocks = {'11','12','21','22'};

for s = 1:length(subjects)
    subject = subjects{s};
    blocks = {'11','12','21','22'};
    
    for b = 1:length(blocks)
        block = blocks{b};
        inpath = [xdfpath,subject,'/Director/'];
        xdfFile = [inpath,'dgame2_',subject,'_Director_',block,'.xdf'];
        mobipath = [inpath,'dgame2_',subject,'_Director_',block,'_MoBI/'];
        outpath = [studyroot,'recordings/audios/',subject,'/'];
        outpath_times = [studyroot,'preproc/helper_files/'];
        director_outfile = [outpath,subject,'_director_',block,'.wav'];
        decke_outfile = [outpath,subject,'_decke_',block,'.wav'];
        if ~isfolder(outpath)
            mkdir outpath;
        end
        

%% read and save the first and last timestamp of ET recording
        tmpXDF = load_xdf(xdfFile,'HandleClockSynchronization',false);
        for items = 1:length(tmpXDF)
            if strcmp(tmpXDF{items}.info.name, 'pupil_capture') == 1
                times = [str2num(tmpXDF{items}.info.first_timestamp);str2num(tmpXDF{items}.info.last_timestamp)];
            end
        end
        dlmwrite([outpath_times,subject,'_timestamps_max-min_',block,'.csv'],times,'precision','%.7f');
        times = [];
        tmpXDF = [];

%% load data with mobilab

        if ~isfolder(mobipath)
            mobilab.allStreams = dataSourceXDF(xdfFile,mobipath);
        else
            mobilab.allStreams = dataSourceMoBI(mobipath);
        end
        exportIndex = mobilab.allStreams.getItemIndexFromItemClass('eeg');
        indexET = mobilab.allStreams.getItemIndexFromItemName('pupil_capture_ieeg-rec-lap-2');
        audioIndex = mobilab.allStreams.getItemIndexFromItemName('audio_xlinc-recording');

%% export ET times for R

        ET = mobilab.allStreams.export2eeglab([indexET]);
        
        times = ET.times/1000;
        dlmwrite([outpath_times,subject,'_times_',block,'.csv'],times,'precision','%.6f');
        times = [];

%% extract audio, normalize and export to wav (normalization is needed or the signal will get clipped = garbage)

        audio = mobilab.allStreams.export2eeglab([audioIndex]);

        %normalize audio
        audio_decke = audio.data(1,:)/max(abs(audio.data(1,:)));
        audio_director = audio.data(2,:)/max(abs(audio.data(2,:)));


        %export to WAV
        audiowrite(director_outfile,audio_director,22050);
        audiowrite(decke_outfile,audio_decke,22050);

    end
end
