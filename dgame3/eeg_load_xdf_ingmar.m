function raw = eeg_load_xdf(filename, varargin)
% Import an XDF file from disk
%
% EEG = eeg_load_xdf(Filename, Options...)
%
% Options:
%   'streamname'  : import only the first stream with the given name
%   'streamtype'  : import only the first stream with the given content type (default: 'EEG')
%   'hostname'    : import only the first stream that also matches this hostname (optional)
%   'effective_rate' : if true, use the effective sampling rate instead of the nominal one
%   'exclude_markerstreams' : cell array of marker stream names to exclude

% -------------------------------------------------------------------------
% parse arguments (added 'hostname')
% -------------------------------------------------------------------------
args = hlp_varargin2struct(varargin, ...
    'streamname','', ...
    'streamtype','EEG', ...
    'hostname','', ...                    % <-- added
    'effective_rate',false, ...
    'exclude_markerstreams',{});

% ensure xdf folder is on path
addpath(fullfile(fileparts(mfilename('fullpath')), 'xdf'));

% -------------------------------------------------------------------------
% load file
% -------------------------------------------------------------------------
streams = load_xdf(filename);

% -------------------------------------------------------------------------
% select data stream
% -------------------------------------------------------------------------
found = false;
for s = 1:length(streams)
    info = streams{s}.info;

    % check optional fields safely
    hasName     = isfield(info,'name');
    hasType     = isfield(info,'type');
    hasHostname = isfield(info,'hostname');

    % check matching criteria
    nameOK     = isempty(args.streamname) || (hasName     && strcmp(info.name, args.streamname));
    typeOK     = isempty(args.streamtype) || (hasType     && strcmp(info.type, args.streamtype));
    hostOK     = isempty(args.hostname)   || (hasHostname && strcmp(info.hostname, args.hostname));

    if nameOK && typeOK && hostOK
        stream = streams{s};
        found = true;
        break;
    end
end

if ~found
    if ~isempty(args.hostname)
        error('No stream found matching name="%s", type="%s", hostname="%s".', ...
            args.streamname, args.streamtype, args.hostname);
    elseif ~isempty(args.streamname)
        error('No stream found with the name "%s".', args.streamname);
    else
        error('No stream found with the type "%s".', args.streamtype);
    end
end

% -------------------------------------------------------------------------
% construct EEGLAB structure
% -------------------------------------------------------------------------
raw = eeg_emptyset;

raw.data = stream.time_series;
[raw.nbchan, raw.pnts, raw.trials] = size(raw.data);
[raw.filepath, fname, fext] = fileparts(filename);
raw.filename = [fname fext];

if args.effective_rate && isfinite(stream.info.effective_srate) && stream.info.effective_srate>0
    raw.srate = stream.info.effective_srate;
else
    raw.srate = str2double(stream.info.nominal_srate);
end
raw.xmin = 0;
raw.xmax = (raw.pnts-1)/raw.srate;

% -------------------------------------------------------------------------
% channel locations and events (unchanged)
% -------------------------------------------------------------------------
% [keep your existing chanlocs/event code exactly as before]
% -------------------------------------------------------------------------
% (everything below this point stays identical)
