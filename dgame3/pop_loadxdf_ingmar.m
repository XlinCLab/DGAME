function [EEG, command] = pop_loadxdf(filename, varargin)

command = '';
filepath = '';
EEG = [];

if nargin < 1

    % ask user
    [filename, filepath] = uigetfile('*.xdf;*.xdfz', 'Choose an XDF file -- pop_loadxdf()');
    drawnow;
    if filename == 0, return; end;

    % popup window parameters
    % -----------------------
    uigeom = { [1 0.5] [1 0.5] [1 0.5] [1 0.5] 0.13};
    uilist = { ...
        { 'style' 'text' 'string' 'Stream name to import:' } ...
        { 'style' 'edit' 'string' '' } ...
        { 'style' 'text' 'string' 'Stream type to import:' } ...
        { 'style' 'edit' 'string' 'EEG' } ...
        { 'style' 'text' 'string' 'Hostname to import:' } ...        % <-- added line
        { 'style' 'edit' 'string' '' } ...                           % <-- added line
        { 'style' 'text' 'string' 'Exclude marker stream(s):' } ...
        { 'style' 'edit' 'string' '{}' } ...
        {} };

    result = inputgui(uigeom, uilist, 'pophelp(''pop_loadxdf'')', 'Load an XDF file');
    if isempty(result), return; end;

    % decode parameters
    % -----------------
    options = [];
    if ~isempty(result{1})
        options = [options, ', ''streamname'', ''' result{1} ''''];
    end
    if ~isempty(result{2})
        options = [options, ', ''streamtype'', ''' result{2} ''''];
    end
    if ~isempty(result{3})
        options = [options, ', ''hostname'', ''' result{3} ''''];    % <-- added line
    end
    if ~isempty(result{4})
        options = [options, ', ''exclude_markerstreams'', ' result{4} ''];
    end

else
    options = vararg2str(varargin);
end

% load data
% ----------
if exist('filepath','var')
    fullFileName = sprintf('%s%s', filepath, filename);
else
    fullFileName = filename;
end

fprintf('Now importing...');
if nargin > 0    
    EEG = eeg_load_xdf(fullFileName, varargin{:});
else
    eval(['EEG = eeg_load_xdf(fullFileName' options ');']);
end
fprintf('done.\n');

EEG = eeg_checkset(EEG);

if length(options) > 2
    command = sprintf('EEG = pop_loadxdf(''%s''%s);', fullFileName, options);
else
    command = sprintf('EEG = pop_loadxdf(''%s'');', fullFileName);
end