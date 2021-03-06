function varargout = onelayergui6(varargin)
%This QCM Analysis Program was created by Lauren Sturdy based on an 
%original version by Kenneth Shull from the Shull Research Group at Northwestern University.
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @onelayergui6_OpeningFcn, ...
    'gui_OutputFcn',  @onelayergui6_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before onelayergui6 is made visible.
function onelayergui6_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
warning('off','MATLAB:Axes:NegativeDataInLogAxis') %Supress warning about negative data
warning('off','optimlib:lsqlin:AlgAndLargeScaleConflict')
set(hObject,'toolbar','auto'); %Remove default toolbar

%Sets constants for use in the program.
handles.constants.f1 = 5e6;
handles.constants.zq = 8.84e6;
handles.constants.error.f = [44 0 120 0 273 0 500]; %Check 7harm
handles.constants.error.g = [11 0 22 0 14 0 10]; %Check 7harm
handles.contourtable.Data = [0 0.2; 0 90; 1 1.2;-0.2 0];

handles.constants.nhvals = {[1,3,1] [1,3,3] [1,5,1] [1,5,5] [3,5,3] [3,5,5]...
    [1,7,1] [1,7,7] [3,7,3] [3,7,7] [5,7,5] [5,7,7] [1,3,5] [1,5,3] [3,5,1]};
handles.constants.label = {'1:3,1' '1:3,3' '1:5,1' '1:5,5' '3:5,3' '3:5,5',...
    '1:7,1', '1:7,7', '3:7,3', '3:7,7', '5:7,5', '5:7,7', '1:3,5', '1:5,3', '3:5,1'};
handles.saveplot.plot = 0; %default to not
resettablesize(hObject,handles)
writetableheaders(hObject, handles, 'h3'); %Sets default size on opening to html h3
set(gcf,'Pointer','arrow'); %Sometimes the cursor gets stuck on as crosshairs

% --- Outputs from this function are returned to the command line.
function varargout = onelayergui6_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function changeqcmfile_Callback(hObject, eventdata, handles, qcmfile, qcmpath)
% This function imports data from a file, stores it (including raw f and g
% data, spectra, and error ranges if available) in handles, and displays plots.

if get(handles.issim, 'value') == 1 %Turn off simulation mode if on
    set(handles.issim, 'value', 0) 
    issim_Callback(hObject, eventdata, handles)
end

set(handles.opencond, 'visible', 'on') %sets the "show conductance" button to visible.
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default') %Clears the status bar

if ~exist('qcmfile','var')
    %Opens dialog to select file to load. It will go to the last used folder by
    %default.
    if isfield(handles, 'din') && isfield(handles.din, 'qcmpath') %Tries to use the folder for the current data
        try
            [qcmfile, qcmpath] = uigetfile({'*.xls;*.xlsx;*.fre;*.mat'}, 'Select the data file', handles.din.qcmpath);
        catch Err
            if strcmp(Err.identifier, 'MATLAB:load:couldNotReadFile')
                [qcmfile, qcmpath]= uigetfile({'*.xls;*.xlsx;*.fre;*.mat'}, 'Select the data file');
            else
                rethrow(Err)
            end
        end
    else
        [qcmfile, qcmpath]=uigetfile({'*.xls;*.xlsx;*.fre;*.mat'}, 'Select the data file'); %Starts in default folder
    end
end

% If someone cancels out of the uigetfile dialog, qcmfile and qcmpath will
% both be 0. This checks if that is the case.
if ~qcmfile
    set(handles.statusbar, 'string', 'You did not select a new file. No data has been changed', 'BackgroundColor', 'yellow')    
    return
end

set(handles.qcmfile,'string', ['Loading ' qcmfile]);
drawnow; %means that the above change will be visible.

%clears the data currently stored in handles.din, as that is all dependent
%on the data file.
if isfield(handles, 'din')
    rmfield(handles, 'din');
end

if isfield(handles, 'saveplot') %Checks if there is previously saved data stored
    rmfield(handles, 'saveplot'); %Clears data stored in saveplot
    handles.saveplot.plot = 0; %Indicates there's no data in saveplot
end

% Determines the file type and the base file name.
[~, filebase, filetypetxt] = fileparts(qcmfile);
name = {'' 'f1' 'g1' 'f3' 'g3' 'f5' 'g5' 'f7' 'g7'}; %Sets names for later.

if strcmp(qcmfile(end-7:end-4),'data') %_data.mat
    % Applied if the data has already been processed and just wants to be
    % viewed. In this mode, no new data file is saved and conductance data
    % is not available. But the data can be viewed.
    handles.filetype = 4;
    load([qcmpath, qcmfile]); %Loads the data in the file
    
    try %Rename the variables to the common names
        qcmt = time;
        cleandelf = delf;
        cleandelg = delg;
    catch
        warning('The data file is of an old type and does not have the proper fields')
        return
    end
    
    try %Loads reference and error data
        for i = 1:nhtoi(length(error.f))
            fref=['f',num2str(itonh(i))];
            gref=['g',num2str(itonh(i))];
            handles.offsettable.Data(i,1)={[commanumber(offset.(fref))...
                setstr(177) num2str(error.f(itonh(i)),'%8.0f')]};
            handles.offsettable.Data(i,2)={[commanumber(offset.(gref))...
                setstr(177) num2str(error.g(itonh(i)),'%8.0f')]};
        end
        if length(error.f) < 7
            handles.offsettable.Data(4,1:2) = {'-'};
            if length(error.f) < 5
                handles.offsettable.Data(3,1:2) = {'-'};
            end
        end
        handles.din.bare.error.f = error.f;
        handles.din.bare.error.g = error.g;
    catch Err
        set(handles.statusbar, 'string', 'Doesn''t contain reference information. Default error values will be used', 'BackgroundColor', 'yellow')
        for n = 2:9
            offset.(name{n}) = 0;
        end
        handles.offsettable.Data(1:4,1:2) = {'-'};
        handles.din.bare.offset = offset;
    end
    
    datasize = size(delf);
    if datasize(2) < 7 %checks number of columns of data (harmonics)
        set(handles.n7, 'value', 0)
        if datasize(2) < 5
            set(handles.n5, 'value', 0)
        end
        handles.activenh = getactivenh(handles);
    end
    
    %This code reads in which harmonic combinations were saved with the
    %data previouly and resets the harmonics to reflect that.
    if exist('nhvals', 'var')
        for i = 1:length(nhvals)
            set(handles.(['nhset' num2str(i) 'select']), 'value', 1); %Checks the checkbox
            %The file saves the actual harmonic combination, not the mnumber.
            %So this identifies which mnumber is associated with that
            %combination. This also works if the mnumbers change in the future
            %to provide a guaranteed reference, which is why I don't feel like
            %doing something simpler.
            p = cellfun(@(x) pdist2(x, [nhvals{i}])==0, handles.constants.nhvals);
            mnumber = find(p);
            set(handles.(['nhset' num2str(i) 'type']), 'value', mnumber);
        end
        %Unchecks any values not calculated for (and doesn't change the
        %harmonic combination)
        if length(nhvals) < 6
            for i = length(nhvals)+1:6
                set(handles.(['nhset' num2str(i) 'select']), 'value', 0);
            end
        end
    end
    
    handles.overwrite = 0; %Disables saving of a new _data file later on.
    set(handles.loadsavedplot, 'visible', 'on') %changes button options
    set(handles.savecleaned, 'visible', 'off')
    set(handles.loadcleaned, 'visible', 'off')
else %not _data.mat files
    % The program can read .xls or .xlsx (older data--Garret's primarily), .fre
    % (QTZ raw data), and .mat (new QCM program). This program assums standard
    % output data for each of those types of files.
    handles.overwrite = 1; %Update _data file each time calculations are done
    set(handles.loadsavedplot, 'visible', 'off') %Changes button options
    set(handles.savecleaned, 'visible', 'on')
    set(handles.loadcleaned, 'visible', 'on')
    if strcmp(filetypetxt,'.fre') %.fre files. Not adapted to 7th harmonic.
        handles.filetype = 1;
        %.fre files have the reference frequencies at the top. These are read in
        %first.
        fid=fopen([qcmpath, qcmfile]); %Opens .fre file for reading
        qcmheader = textscan(fid,'%s %s %s %s %s %s %s',1); %gets top line
        fclose(fid); %closes file to save memory (and because it isn't necessary)
        
        %This provides the reference data. First it checks if there is a _bare
        %file.
        filepathbase = [qcmpath filebase];
        [reference error.f error.g] = getreferencedata(hObject, handles, filepathbase);%, 1); %reference = 0 if unsuccessful
        if length(reference) > 1 %Uses _bare file if possible
            for n=2:7 %Assumes only 3 harmonics were measured. A fair assumption.
                offset.(name{n}) = reference(n);
            end
        else %Uses header values
            for n=2:7
                off = char(qcmheader{n});
                %This removes the Hz from the end of the cell and makes it a number
                offset.(name{n}) = str2num(regexprep(off,'Hz',''));
            end
        end
        
        %Reads in the actual QCM data (below the first line). Whether it is
        % shifted or absolute will be dealt with later.
        qcmdata=dlmread([qcmpath, qcmfile],'\t',1,0);
        
        %The rest of this tries to "clean up" the data and weed out erroneous
        %values before displaying them.
        qcmsubdataf = qcmdata(:,2:2:7); %breaks the data into chunks for easier processing
        qcmsubdatag = qcmdata(:,3:2:7);
        
        % Now the program has a quirk that sometimes it puts the
        % reference data instead of the point if no data is collected. To
        % determine if this is the case, use mode to check for repeated values
        [num, occur] = mode(qcmsubdatag);
        %if the number of occurences is greater than 1, it is probably the
        %reference value. This loop sets such values to NaN.
        for i = 1:3
            if occur(i) > 2
                qcmsubdataf(qcmsubdatag == num(i)) = NaN;
                qcmsubdatag(qcmsubdatag == num(i)) = NaN;
            end
        end
        
        % Another source of erroneous values is when the 15Mhz and 25MHz peaks
        % are actually the large peak at 21 (for Hera QTZ data). This removes
        % those as well, and assums that no real peaks will be between 16MHz and 24MHz.
        qcmsubdatag(qcmsubdataf < 24000000 & qcmsubdataf > 16000000) = NaN;
        qcmsubdataf(qcmsubdataf < 24000000 & qcmsubdataf > 16000000) = NaN;
        
        % Now that some extraneous values have been removed, the two arrays are
        % put back together.
        qcmcleandata(:,[2 4 6]) = qcmsubdataf;
        qcmcleandata(:,[3 5 7]) = qcmsubdatag;
        
    elseif strcmp(filetypetxt, '.mat') %.mat files (Josh's QCM program)
        handles.filetype = 2;
        load([qcmpath, qcmfile]) %Load the mat file into memory.
        
        % Since the file is saved with a million rows and they shouldn't all be
        % populated, this determines the extent of the data by checking for the
        % index of the largest time and assuming that is the highest index.
        [~, idx] = max(abs_freq(:,1));
        qcmdata = abs_freq(1:idx,1:9); %Extracts data for times and harmonics of interest
        %There is at least one file which starts off with some time NaNs's,
        %which need to be removed before they cause trouble later. I don't care
        %about any rows that begin with NaN
        nantime = isnan(abs_freq(1:idx, 1));
        qcmdata = qcmdata(~nantime, 1:9);
        clear abs_freq %It's a large array and may slow things down
        
        qcmsubdataf = qcmdata(:,2:2:9); %Splits for simplicity
        qcmsubdatag = qcmdata(:,3:2:9);
        
        qcmfits=chisq_values(~nantime, 1:9); %removes NaNs from this too
        qcmsubfitsf = qcmfits(:,2:2:9);
        
        % Sets to NaN anything with a bad fit
%         qcmsubdataf(qcmsubfitsf>1e-2) = NaN;
%         qcmsubdatag(qcmsubfitsf>1e-2) = NaN;
        
        %Assigns "clean" data to its own variables.
        qcmcleandata(:,[2 4 6 8]) = qcmsubdataf;
        qcmcleandata(:,[3 5 7 9]) = qcmsubdatag;
        
        %Checks for a _bare file to get reference data from, if not, uses
        %freq_shift_ref.
        [reference error.f error.g] = getreferencedata(hObject, handles, [qcmpath filebase]);
        if length(reference) > 1 %True if there wasn't an error
            for n=2:9
                if ~isnan(reference(n))
                    offset.(name{n}) = reference(n);
                else
                    offset.(name{n}) = NaN;
                    set(handles.(['n' num2str(floor(n/2)*2-1)]), 'value', 0)
                end
            end
        else %If there was an error in getting the reference, for isntance, there is no _bare file.
            reference = freq_shift_ref; %Use reference data from the mat file
            for n = 1:4
                offset.(name{n*2}) = reference(1,n);
                offset.(name{n*2+1}) = reference(2,n);
            end
            set(handles.statusbar, 'string', 'A reference data file was not found. Using data from .mat file', 'BackgroundColor', 'red')
        end
        
    elseif strcmp(filetypetxt,'.xls') || strcmp(filetype,'.xlsx') %Not corrected for 7th harmonic
        %Since I don't use excel files, this part is not updated, but will work
        % with old excel files, for instance, Garrett's data.
        handles.filetype = 3;
        tempdata=xlsread([qcmpath, qcmfile]);
        % set offsets to zero for  now - can change this if necessary
        % This assumes that all of the data is in delta form already
        clear qcmdata
        for n = 2:15
            offset.(name{n}) = 0;
        end
        [nrows,~]=size(tempdata);
        tempdata(isnan(tempdata))=0;
        % if all values are zero, then there's not really any data in this row
        harmlength =  min(size(tempdata));
        index=1;
        for i=1:nrows
            rowcheck=sum(tempdata(i,2:harmlength));
            if rowcheck~=0
                qcmdata(index,1:harmlength)=tempdata(i,1:harmlength);
                index=index+1;
            end
        end
        qcmcleandata = qcmdata;
        error.f = handles.constants.error.f; %Default error
        error.g = handles.constants.error.g; %Default error        
    end
      
    %Checks if there is 5th and 7th harmonic data 
    [nrows ncols]=size(qcmdata);
    if isempty(find(qcmdata(:,9) ~= 0)) %No data for 7th harmonic
        set(handles.n7, 'value', 0)
        if isempty(find(qcmdata(:,7) ~= 0)) %No data for 5th harmonic
            set(handles.n5, 'value', 0)
        end
    end
    handles.activenh = getactivenh(handles); %Updates active harmonics
    nhvals = handles.activenh.on;
    
    %This part writes the reference values for the nhvals in use to the gui
    for nh = 1:2:7
        if any(nhvals==nh)
            fref=['f',num2str(nh)];
            gref=['g',num2str(nh)];
            handles.offsettable.Data(nhtoi(nh),1)={[commanumber(offset.(fref)) setstr(177) num2str(error.f(nh),'%8.0f')]};
            handles.offsettable.Data(nhtoi(nh),2)={[commanumber(offset.(gref)) setstr(177) num2str(error.g(nh),'%8.0f')]};
        else
            handles.offsettable.Data(nhtoi(nh),1)={['0' setstr(177) '0']};
            handles.offsettable.Data(nhtoi(nh),2)={['0' setstr(177) '0']};
        end
    end
    
    %Extracts out the time data
    qcmt(1:nrows)=qcmdata(1:nrows,1);
    [qcmt,indexm,~]=unique(qcmt); %Not sure why this is here. Could possibly be removed?
    
    % Checks if the values are shifted or absolute. This check assumes the
    % frequency shift at 5MHz is less than 1MHz. If the first frequency value
    % is less than 4MHz, it assumes the data is already given in shifts and
    % removes the offset.
    if qcmdata(2,2)<4000000
        for n = 2:max(nhvals)*2+1
            offset.(name{n}) = 0;
        end
    end
    
    %Shifts the data if necessary and saves it.
    for nh=nhvals
        fref=['f',num2str(nh)];
        gref=['g',num2str(nh)];
        delf(:,nh)=(qcmdata(indexm,nh+1)-offset.(fref));
        delg(:,nh)=(qcmdata(indexm,nh+2)-offset.(gref));
        absf(:,nh)=(qcmcleandata(indexm,nh+1));
        absg(:,nh)=(qcmcleandata(indexm,nh+2));
        cleandelf(:,nh)=(qcmcleandata(indexm,nh+1)-offset.(fref));
        cleandelg(:,nh)=(qcmcleandata(indexm,nh+2)-offset.(gref));
    end
    
    %Only applies threshhold filters if the sample isn't a bare crystal.
    %This removes data points for which df/n > 100,000. Shifts shouldn't be
    %that big.
    if get(handles.isbare, 'value') == 0
        for nh = nhvals
            fieldname = ['nh' num2str(nh)];
            filter.(fieldname) = abs(cleandelf(:,nh)./nh) > 1e5;
            cleandelf(filter.(fieldname), nh) = NaN;
            cleandelg(filter.(fieldname), nh) = NaN;
        end
    end
end
%Saves all of the processed data into the handles structure.
handles.din.qcmt = qcmt;
handles.din.delf = delf;
handles.din.delg = delg;
handles.din.absf = absf;
handles.din.absg = absg;
handles.din.cleandelf = cleandelf;
handles.din.cleandelg = cleandelg;
handles.din.qcmpath = qcmpath;
handles.din.qcmfile = qcmfile;
handles.din.filebase = filebase;

try
    handles.din.bare.error.f = error.f;
    handles.din.bare.error.g = error.g;
    handles.din.bare.offset = offset;
catch
    set(handles.statusbar, 'string', 'Using default error', 'BackgroundColor', 'red')
    handles.din.bare.error = handles.constants.error;
end

%Plots the imported data to the lower plots
plotraw(hObject, eventdata, handles);

%Imports the conductance data, if any
conductance = getconductance(hObject, handles);
handles.cond = conductance;

%Changes the data file to not be "loading" anymore.
set(handles.qcmfile,'string',qcmfile);
set(handles.statusbar, 'string', 'File successfully loaded', 'BackgroundColor', 'green')
if ~strcmp(qcmfile(end-7:end-4),'data')
    handles = loadcleaned_Callback(hObject, eventdata, handles);
end
guidata(hObject, handles); %Saves new handles.

function [reference errorf errorg] = getreferencedata(hObject, handles, filepathbase);
% geterrorrange tries to open a '_bare' file of the same type as the main
% file and calculates the range of values for each quantity in it and saves
% these as the error range values.
defaulterrorf = handles.constants.error.f;
defaulterrorg = handles.constants.error.g;
% Base filename of potential bare file. 
if get(handles.isbare, 'value')
    filename = filepathbase;
    errorf = [0 0 0 0 0 0 0];
    errorg = [0 0 0 0 0 0 0];
    reference = [0 0 0 0 0 0 0 0 0];
else
    filename = [filepathbase '_bare'];
    % Method for calculating for .fre files
    if handles.filetype == 1 && exist([filename '.fre'], 'file') == 2
        barefile = [filename '.fre'];
        rawdata = importdata(barefile);
        reference = mean(rawdata.data(:,1:9));
        stdev = std(rawdata.data(:,1:9), 'omitnan');
        errorf = [stdev(2) 0 stdev(4) 0 stdev(6) 0 stdev(8)];
        errorg = [stdev(3) 0 stdev(5) 0 stdev(7) 0 stdev(9)];
       
    elseif handles.filetype == 2  % Method for .mat files
        if exist([filename '.mat'], 'file') == 2
            if exist([filename '_data.mat'], 'file') == 2
                load([filename '_data.mat'], 'absf', 'absg'); %loads absf and absg (among others)
                maxlength = min(size(absf));
                freqmean = nanmean(absf(:,1:maxlength));
                dissmean = nanmean(absg(:,1:maxlength));
                if maxlength == 5
                    reference = [0 freqmean(1) dissmean(1) freqmean(3) dissmean(3) freqmean(5) dissmean(5) NaN NaN];
                elseif maxlength == 7
                    reference = [0 freqmean(1) dissmean(1) freqmean(3) dissmean(3) freqmean(5) dissmean(5) freqmean(7) dissmean(7)];
                else
                    reference = 0;
                end
                errorf = std(absf, 'omitnan');
                errorg = std(absg, 'omitnan');
            else %use unedited bare file
                barefile = [filename '.mat'];
                load(barefile, 'abs_freq')
                [~, idx] = max(abs_freq(:,1));
                reference = nanmean(abs_freq(1:idx,1:9));
                stdev = std(abs_freq(1:idx,1:9), 'omitnan');
                errorf = [stdev(2) 0 stdev(4) 0 stdev(6) 0 stdev(8)];
                errorg = [stdev(3) 0 stdev(5) 0 stdev(7) 0 stdev(9)];
                set(handles.statusbar, 'string', 'Bare file may not have been looked at', 'BackgroundColor', 'yellow')
                drawnow
            end
        elseif exist([filename '.fre'], 'file') == 2 %If fre file exists
            set(handles.statusbar, 'string', 'Using .fre bare file', 'BackgroundColor', 'yellow')
            drawnow
            try
                barefile = [filename '.fre'];
                rawdata = importdata(barefile);
                stdev = std(rawdata.data(:,1:9), 'omitnan');
                errorf = [stdev(2) 0 stdev(4) 0 stdev(6) 0 stdev(8)];
                errorg = [stdev(3) 0 stdev(5) 0 stdev(7) 0 stdev(9)];
                reference = mean(rawdata.data(:,1:9));
            catch
                set(handles.statusbar, 'string', 'Default error values will be used since no bare data was found.', 'BackgroundColor', 'yellow')
                errorf = defaulterrorf;
                errorg = defaulterrorg;
                reference = 0;
            end
        else
            [barefile, barepath] = uigetfile({'*.fre;*.mat'}, 'Select the bare data file', filepathbase);
            try    
            load([barepath barefile(1:end-4) '_data.mat'], 'absf', 'absg'); %loads delf and delg (among others)
                maxlength = min(size(absf));
                freqmean = nanmean(absf(:,1:maxlength));
                dissmean = nanmean(absg(:,1:maxlength));
                if maxlength == 5
                    reference = [0 freqmean(1) dissmean(1) freqmean(3) dissmean(3) freqmean(5) dissmean(5) NaN NaN];
                elseif maxlength == 7
                    reference = [0 freqmean(1) dissmean(1) freqmean(3) dissmean(3) freqmean(5) dissmean(5) freqmean(7) dissmean(7)];
                else
                    reference = 0;
                end
                errorf = std(absf, 'omitnan');
                errorg = std(absg, 'omitnan');
            catch %use unedited bare file
                load([barepath barefile], 'abs_freq')
                [~, idx] = max(abs_freq(:,1)); %Finds the last timepoint with data
                reference = nanmean(abs_freq(1:idx,1:9));
                stdev = std(abs_freq(1:idx,1:9), 'omitnan');
                errorf = [stdev(2) 0 stdev(4) 0 stdev(6) 0 stdev(8)];
                errorg = [stdev(3) 0 stdev(5) 0 stdev(7) 0 stdev(9)];
                set(handles.statusbar, 'string', 'Bare file may not have been looked at', 'BackgroundColor', 'yellow')
                drawnow
            end

        end
    elseif handles.filetype == 3 % Method for excel files
        set(handles.statusbar, 'string', 'Default error values will be used since there isn''t a procedure for Excel files.', 'BackgroundColor', 'yellow')
        errorf = defaulterrorf;
        errorg = defaulterrorg;
        reference = 0;
    else % If file not found, use defaults
        set(handles.statusbar, 'string', 'Default error values will be used since no bare data was found.', 'BackgroundColor', 'yellow')
        errorf = defaulterrorf;
        errorg = defaulterrorg;
        reference = 0;
    end
    minerror = [1 0 1 0 1 0 1];
    errorg(errorg<1) = minerror(errorg<1);
end

function Solve_Callback(hObject,eventdata,handles)
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')
% Uses the top selected harmonic combination to solve the currently loaded
% point.
mnumber = getnhvals(handles);
[~, handles] = findsolution(hObject,eventdata,handles,1,mnumber(1));
guidata(hObject, handles); %Saves new handles.

function [outputs, handles] = findsolution(hObject, eventdata, handles, time, mnumber)
nhset = handles.constants.nhvals{mnumber};
handles.activenh = getactivenh(handles);
nhon = handles.activenh.on;
%Reads the values of df and dg out of the saved data shown on the screen
df(nhon) = handles.deldatadata(nhtoi(nhon),1);
dg(nhon) = handles.deldatadata(nhtoi(nhon),3);

refG = itonh(get(handles.modulustoplot, 'value'));
[outputs, handles] = solvedfstar(df, dg, mnumber, hObject, handles, refG);

if outputs.error ~= 1 %Possible reasons for solving incorrectly
    handles.propertytable.Data(1,5)={'not solved'};
    if outputs.error == 2.5 % Some data NaN
        set(handles.statusbar, 'string', ['Did not solve at time ' ...
            num2str(time) ' using ratio ' num2str(nhset(1)) ':' ...
            num2str(nhset(2)) ',' num2str(nhset(3)) ...
            ' because some data was NaN'], 'BackgroundColor', 'yellow')
    elseif outputs.error == 3.5 % Phi > 90
        set(handles.statusbar, 'string', ['Did not solve at time ' ...
            num2str(time) ' using ratio ' num2str(nhset(1)) ':' ...
            num2str(nhset(2)) ',' num2str(nhset(3)) ...
            ' because phi out of range'], 'BackgroundColor', 'yellow')
    elseif outputs.error == 1.5 % Misc error related to fsolve
        set(handles.statusbar, 'string', ['Did not solve at time ' ...
            num2str(time) ' using ratio ' num2str(nhset(1)) ':' ...
            num2str(nhset(2)) ',' num2str(nhset(3)) ...
            ' because unable to solve due to error in fsolve'],...
            'BackgroundColor', 'yellow')
    else % Other errors from fsolve--0, 2, etc.
        set(handles.statusbar, 'string', ['Did not solve at time ' ...
            num2str(time) ' using ratio ' num2str(nhset(1)) ':' ...
            num2str(nhset(2)) ',' num2str(nhset(3)) ' with error ' ...
            num2str(outputs.error)], 'BackgroundColor', 'yellow')
    end
else
    % update the data table on the gui with the updated values (but only if
    % it solved correctly)
    handles.propertytable.Data(1,1)={num2str(outputs.drho,'%8.4f')};
    handles.propertytable.Data(1,2)={num2str(outputs.grho(refG),'%8.3e')};
    handles.propertytable.Data(1,3)={num2str(outputs.phi,'%8.2f')};
    handles.propertytable.Data(1,5)={[num2str(nhset(1)) ':' num2str(nhset(2)) ',' num2str(nhset(3))]};
    
    [contourplots, handles] = redrawtable(hObject, handles, mnumber, outputs);
    assignin('caller','contourplots',contourplots)
    if get(handles.calcerror, 'value') == 1 && length(unique(nhset(1:3))) < 3 %No error calc yet for 1:3,5 etc.
        [drhoe grhoe phie] = finderror(hObject, handles, mnumber);
        outputs.grhoe = grhoe;
        outputs.drhoe = drhoe;
        outputs.phie = phie;
    else
        outputs.grhoe = NaN;
        outputs.drhoe = NaN;
        outputs.phie = NaN;
        handles.propertytable.Data(2,1:3)={'-','-','-'};
    end
end

outputs.time = time;
guidata(hObject, handles); %Saves new handles.



function Sim_Callback(hObject, eventdata, handles)
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')
% Does a theoretical calculation for the values currently entered for mass,
% modulus, density, and phi.
f1=handles.constants.f1;

% if get(handles.etabutton,'value')
%     etarho(1)=1000*str2double(handles.propertytable.Data(1,4));
%     grho(1)=etarho*2*pi*f1;
%     handles.propertytable.Data(1,2)={num2str(0.001*grho(1),'%6.1e')};
% else
    grho(1) = 1000*str2double(handles.propertytable.Data(1,2));
    etarho(1)=grho(1)/(2*pi*f1);
    handles.propertytable.Data(1,4)={num2str(0.001*etarho(1),'%5.2e')};
% end

mnumber = getnhvals(handles);
[contourplots, handles] = redrawtable(hObject, handles, mnumber(1));
handles.propertytable.Data(1,5)={['sim']};
handles.contourplots = contourplots;
if get(handles.calcerror, 'value')
    finderror(hObject, handles, mnumber(1))
else
    handles.propertytable.Data(2,1:3) = {'-'};
end
guidata(hObject, handles);

function F = fstar(x,nh,harmratio,dissratio, refn)
% Compares ideal and measured values for fstar for use in the solve
% function. The ideal output of the function is [0 0]. x(1) is the phase
% angle and x(2) is d/lambda.
phi = x(1);  % phase angle in degrees
drefn = x(2);
d(1)=drefn*(nh(1)/refn)^(1-phi/180); % d/lambda at nh(1)
d(2)=drefn*(nh(2)/refn)^(1-phi/180); % d/lambda at nh(2)
d(3)=drefn*(nh(3)/refn)^(1-phi/180); % d/lambda at nh(3)

for n=1:3 
    delf(n)=delfstar(d(n),phi);
end
harmratiocalc=real(delf(2))/real(delf(1));
dissratiocalc=imag(delf(3))/real(delf(3));
F=[dissratio-dissratiocalc;
    harmratio-harmratiocalc]; %difference in ratios.

function F=delfstar(d,phi)  % input phase angle is in degrees
% Calculates delfstar (rhs of delf/delfsn equation) with input of d/lambda
% and phi
F=-(1/((2*pi*d)*(1-1i*tand(phi/2))))* ...
    tan(2*pi*d*(1-1i*tand(phi/2)));

function F=sauerbrey(n,drho)
% Calculates the sauerbry shift based on the harmonic and areal density.
F=2*n*5e6^2*drho/8.84e6;

function conductance = getconductance(hObject, handles)
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')
%The rest of the file deals with things like plotting and finding the
%conductance curves.
%This finds all of the files in the same folder as the .fre file that have
%the extension .spc, which is associated with the conductance curves.

if handles.filetype == 2
    %if mat file, look first for a matching file, otherwise assume this is
    %a tempcorr file and look for the base name.
    try
        files = who('-file', [handles.din.qcmpath handles.din.qcmfile(1:end-4) '_raw_spectras.mat']);
        samefiles = 1:length(files)-1;
    catch Err
        try
            files = who('-file', [handles.din.qcmpath handles.din.qcmfile(1:end-12) '_raw_spectras.mat']);
        catch Err
            handles.cond = [];
            conductance = [];
            set(handles.opencond, 'visible', 'off')
            guidata(hObject,handles)
            return
        end
    end
    samefiles = [];
    for i = 1:length(files)
        if ~strcmp(files{i}, 'reference') && ~strcmp(files{i}, 'version')
            samefiles = [samefiles i];
        end
    end
elseif handles.filetype == 4 %Loading a _data.mat file (no spectra available)
    conductance = [];
    set(handles.opencond, 'visible', 'off')
    guidata(hObject,handles)
    return
else
    % if not a mat file, will be a bit more complicated. This opens the
    % folder and finds the files.
    files = dir(fullfile(handles.din.qcmpath,'*.spc')); %Find .spc files in folder
    numcond = max(size(files));
    %Troubleshoot not finding files.
    if sum(size(files)) == 1
        choice = questdlg(['The folder that you have selected does not contain any '...
            'conductance spectra. Either there are no conductance specra, or '...
            'they are in a different folder.'], 'Error finding spectra',...
            'Select Spectra from different folder', 'Ok', 'Ok');
        switch choice
            case 'Select Spectra from different folder'
                [~,PathName] = uigetfile('*.spc','Select the conductance files', handles.filename{1});
                handles.filename{1} = PathName;
                files = dir(fullfile(handles.filename{1},'*.spc'));
                if sum(size(files))
                    warndlg('Good try. But there still aren''t any condutance spectra')
                    conductance = [];
                    return
                end
            case 'Ok'
                handles.cond = [];
                conductance = [];
                set(handles.opencond, 'visible', 'off')
                guidata(hObject,handles)
                return
        end
    end
    %I want to know for each of the .spc files if it is associated with the
    %loaded file and not another one in the same folder, so this sorts through
    %finds the ones that have the same base filename.
    samefiles = [];
    
    for i = 1:numcond
        try
            %Looks for matching string. Assumes only appendage would be
            %'_bare'
            rightstring = strfind(files(i).name, handles.din.qcmfile(1:end-4));
            barestring = strfind(files(i).name, 'bare');
            if ~isempty(rightstring) && isempty(barestring)
                samefiles = [samefiles i];
            end
        catch Err
            %if ~strcmp(Err.identifier, 'MATLAB:nonLogicalConditional')
            rethrow(Err)
            %end
        end
    end
    
    if isempty(samefiles)
        choice = questdlg(['The folder that you have selected does not contain any '...
            'conductance spectra that match this file name.'], 'Error finding spectra',...
            'Select a different base name', 'Ok', 'Ok');
        switch choice
            case 'Select a different base name'
                [Filename,PathName] = uigetfile('*.fre','Select the base filename', handles.din.qcmpath);
                filebase = Filename(1:end-4);
                for i = 1:numcond
                    try
                        rightstring = strfind(files(i).name, filebase);
                        barestring = strfind(files(i).name, 'bare');
                        if rightstring && isempty(barestring)
                            samefiles = [samefiles i];
                        end
                    catch Err
                        if ~strcmp(Err.identifier, 'MATLAB:nonLogicalConditional')
                            rethrow Err
                        end
                    end
                end
            case 'Ok'
                handles.cond = [];
                handles.cond = [];
                conductance = [];
                set(handles.opencond, 'visible', 'off')
                guidata(hObject,handles)
                return
        end
    end
end

%Creates a matrix to hold the information about the conductance
%spectra--which time points and which harmonics the spectra are availabe
%for.
condfiles = zeros(length(samefiles),2);
%My filenames have a lot of underscores, but most probably don't. This
%hopefully determines how many underscores there are in the usual
%conductance filenames so the correct index can be called to extract the
%time information.

index = 1;
handles.cond.spectra = [];
for i = samefiles
    %Reads the filename and breaks it into chunks by the underscore
    %delimiter
    
    if handles.filetype == 2
        parsed = textscan(files{i},'%s','delimiter','_');
        condfiles(index,1) = str2double(regexprep([parsed{1}{length(parsed{1})-4}], 'dot', '.'));
    else
        parsed = textscan(files(i).name,'%s','delimiter','_');
        condfiles(index,1) = str2double([parsed{1}{length(parsed{1})-4}]);
    end
    %Saves the time index (which is 4th fron the end)
    %Gets and saves the number of the harmonic, 1-3
    harm = parsed{1}{end};
    condfiles(index,2) = str2double(harm(1));
    %Determines the associated index of the main data (freq, diss, etc.)
    [diff idx] = min(abs(handles.din.qcmt-condfiles(index,1)));
    if diff<0.05 %0.05 seems to be a good distinguisher between close and the same
        condfiles(index,3) = int32(idx);
    else
        condfiles(index,3) = 0;
    end
    
    if handles.filetype == 2
        handles.cond.filename{index} = files{i};
    else
        handles.cond.filename{index} = [handles.din.qcmpath files(i).name];
    end
    index = index + 1;
    
end
handles.cond.points = condfiles;
conductance = handles.cond;
guidata(hObject, handles);

function size = fontsize_Callback(hObject, eventdata, handles)  % change font size of gui
% Changes size of headers and fonts throughout
switch get(handles.fontsize,'Value')
    case 1
        fontsize=10;
        size = 'h4';
    case 2
        fontsize=11;
        size = 'h4';
    case 3
        fontsize=12;
        size = 'h3';
    case 4
        fontsize=13;
        size = 'h3';
    case 5
        fontsize=14;
        size = 'h2';
    case 6
        fontsize=15;
        size = 'h2';
    case 7
        fontsize=16;
        size = 'h2';
    case 8
        fontsize=17;
        size = 'h2';
    case 9
        fontsize=18;
        size = 'h1';
    case 10
        fontsize=19;
        size = 'h1';
    case 11
        fontsize=20;
        size = 'h1';
    otherwise
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize) %Changes fontsize
set(findall(gcf,'-property','FontSize'),'FontUnits','points')
handles = writetableheaders(hObject, handles, size); %Chanes headers (which are html)
resettablesize(hObject,handles) %Changes size of table to account for new fontsize
guidata(hObject, handles)

function savebutton_Callback(hObject, eventdata, handles)
%saves the fig file
guidata(hObject,handles)
hgsave('onelayergui6.fig')

function plotraw(hObject, eventdata, handles)
%Brings data into the function
try
    qcmt = handles.din.qcmt;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        set(handles.statusbar, 'string', 'Unable to replot data. Please reload file', 'BackgroundColor', 'red')
        return
    end
end
% Plots filtered or unfiltered data depending on radio button
if get(handles.cleaned, 'value')
    delf = handles.din.cleandelf;
    delg = handles.din.cleandelg;
else
    delf = handles.din.delf;
    delg = handles.din.delg;
end

% now we plot the data
% Set defaults, symbols, colors, etc.
cla(handles.axes2)
cla(handles.axes1)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinemarkersize',12)
set(0,'defaultlinelinewidth',1.5)
symbols{1}='r+-';
symbols{3}='go-';
symbols{5}='bs-';
symbols{7}='kd-';
colors{1}=[1,0,0];
colors{3}=[0,0.5,0];
colors{5}=[0,0,1];
colors{7}=[0,0,0];
legends{1}='n=1';
legends{3}='n=3';
legends{5}='n=5';
legends{7}='n=7';
legendtext={};

% Sometimes handles.activenh.on doesn't exist, so this tries to fix that.
% I'm pretty sure that I have tried an "if field" check to solve the
% problem, but it didn't work for some reason.
try
    handles.activenh.on;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        handles.activenh = getactivenh(handles);
    else
        rethrow(Err)
    end
end

try %Sometimes these values aren't saved--the try fixes it if not
    xlabeltext = handles.plotting.xlabeltext;
    timefactor = handles.plotting.timefactor;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        handles = xlabeltext_Callback(hObject, eventdata, handles);
        xlabeltext = handles.plotting.xlabeltext;
        timefactor = handles.plotting.timefactor;
    else
        rethrow(Err)
    end
end

for nh=handles.activenh.on
    plot(handles.axes1,qcmt./timefactor,delf(:,nh)/nh,symbols{nh},'color',colors{nh})
    plot(handles.axes2,qcmt./timefactor,delg(:,nh),symbols{nh}, 'color',colors{nh})
    legendtext=[legendtext,legends{nh}];
    hold(handles.axes1,'on')
    hold(handles.axes2,'on')
end

xlabel(handles.axes1, xlabeltext)
ylabel(handles.axes1,'{\Delta}f/n (Hz)')
legend(handles.axes1, legendtext,'location','best')

xlabel(handles.axes2, xlabeltext)
ylabel(handles.axes2,'\Delta\Gamma (Hz)')
legend(handles.axes2, legendtext,'location','best')

set(handles.axes1,'xlimmode','auto')
set(handles.axes1,'ylimmode','auto')
set(handles.axes2,'xlimmode','auto')
set(handles.axes2,'ylimmode','auto')

if get(handles.log,'value') %Sets time axis scale
    set([handles.axes1 handles.axes2],'xscale','log')
else
    set([handles.axes1 handles.axes2],'xscale','linear')
end

if get(handles.loggamma,'value') %Sets gamma axis scale
    set(handles.axes2,'yscale','log')
else
    set(handles.axes2,'yscale','linear')
end

linkaxes([handles.axes1, handles.axes2], 'x') %Links x-axes
guidata(hObject,handles);

function [str]=commanumber(num)
num=round(num);
str = num2str(num);
fin = length(str);
for i = fin-2:-3:2
    str(i+1:end+1) = str(i:end);
    str(i) = ',';
end
newfin=length(str);
% avoid -, output at beginning of string
if strcmp(str(1),'-') && strcmp(str(2),',')
    str=[str(1) str(3:newfin)];
end

function cursor_Callback(hObject, eventdata, handles)
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')
handles.activenh = getactivenh(handles);
checksolveharms(hObject, handles);
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
k=0; % keeps track of the number of points
try
    struct2var(handles.din)
catch Err
    if strcmp('MATLAB:nonExistentField', Err.identifier)
        warning = warndlg('The data was not saved correctly. Please reload the file')
        uiwait(warning)
        return
    end
end

nhvals = handles.constants.nhvals;
mnumber = getnhvals(handles);

% This doesn't work if there are lots of NaNs, so the following fills in
% the holes. Which actually doesn't seem a good idea in some cases, but
% here goes.
cleandelf(isnan(cleandelf)) = interp1(find(~isnan(cleandelf)), cleandelf(~isnan(cleandelf)), find(isnan(cleandelf)), 'PCHIP');
cleandelg(isnan(cleandelg)) = interp1(find(~isnan(cleandelg)), cleandelg(~isnan(cleandelg)), find(isnan(cleandelg)), 'PCHIP');

outputs={};
 
while but == 1
    k = k+1;
    handles.datapoints = k;
   [time,~,but] = ginput(1);
    try
        timefactor = handles.plotting.timefactor;
    catch
        handles = xlabeltext_Callback(hObject, eventdata, handles);
        timefactor = handles.plotting.timefactor;
    end
    for nh=handles.activenh.on
        df = interp1(qcmt,cleandelf(:,nh),time*timefactor);
        dg = interp1(qcmt,cleandelg(:,nh),time*timefactor);
        handles.deldatadata(nhtoi(nh),1) = df;
        handles.deldatadata(nhtoi(nh),3) = dg;
        handles.deldatatable.Data(nhtoi(nh),1) = {commanumber(handles.deldatadata(nhtoi(nh),1))};
        handles.deldatatable.Data(nhtoi(nh),3) = {commanumber(handles.deldatadata(nhtoi(nh),3))};
    end
    for nh=handles.activenh.off
        handles.deldatadata(nhtoi(nh),1) = NaN;
        handles.deldatadata(nhtoi(nh),3) = NaN;
        handles.deldatatable.Data(nhtoi(nh),1) = {' -'};
        handles.deldatatable.Data(nhtoi(nh),3) = {' -'};
    end
    
    if get(handles.autosolve,'value')
        for m = mnumber
            [outputs{k,m}, handles] = findsolution(hObject, eventdata, handles, time, m);
        end
    end
    drawnow;
end

handles.uitoggletool1.Enable = 'on'; %Turns zoom in back on
handles.uitoggletool2.Enable = 'on'; %Turns zoom out back on
handles.uitoggletool3.Enable = 'on'; %Turns pan back on
handles.uitoggletool4.Enable = 'on'; %Turns datatip back on

if k > 5 %Sets a threshold for asking about saving the data to the _data file
    ask = 1;
else
    ask = 0;
end
outputstructure = reformatdata(handles, outputs, ask);
handles.dout = outputstructure;
if exist('contourplots', 'var') %Should be created during findsolution if desired
    handles.contourplots = contourplots;
end
guidata(hObject,handles);
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')


function contourplots = plotcontours(handles, outputs, mnumber)
maingui = gcf;
try
    if get(handles.issim, 'value')
        comb = handles.constants.nhvals;
        n1 = comb{mnumber}(1);
        n2 = comb{mnumber}(2);
        n3 = comb{mnumber}(3);
        delf(1:2:5) = handles.deldatadata(1:3,2);
        delg(1:2:5) = handles.deldatadata(1:3,4);
        df1 = delf(n1);
        df2 = delf(n2);
        dg3 = delg(n3);
        df3 = delf(n3);
        drho = str2num(handles.propertytable.Data{1,1});
        grho = str2num(handles.propertytable.Data{1,2});
        phiout = str2num(handles.propertytable.Data{1,3});
        d1out = handles.deldatadata(1,8);
    else
        comb = handles.constants.nhvals;
        n1 = comb{mnumber}(1);
        n2 = comb{mnumber}(2);
        n3 = comb{mnumber}(3);
        delf = outputs.df;
        delg = outputs.dg;
        df1 = delf(n1);
        df2 = delf(n2);
        dg3 = delg(n3);
        df3 = delf(n3);
        drho = outputs.drho;
        grho = outputs.grho;
        phiout = outputs.phi;
        d1out = outputs.d;
    end
catch
    comb = handles.constants.nhvals;
    n1 = comb{mnumber}(1);
    n2 = comb{mnumber}(2);
    n3 = comb{mnumber}(3);
    delf(1:2:5) = handles.deldatadata(1:3,1);
    delg(1:2:5) = handles.deldatadata(1:3,3);
    df1 = delf(n1);
    df2 = delf(n2);
    dg3 = delg(n3);
    df3 = delf(n3);
end

dlofn=@(n,d1,phi) d1*n^(1-phi/180); %d/lambda
Dn=@(n,d1,phi) 2*pi*dlofn(n,d1,phi)*(1-1i*tand(phi/2));
dfstardn = @(Dn) -tan(Dn)/Dn;
dfstar = @(n, d1, phi) dfstardn(Dn(n, d1, phi));
dlof1 = @(n, dn, phi) dn*n^(phi/180-1); %Takes higher order dl and turns it into d1.
frh = @(dn, n, phi) real(dfstar(n2, dlof1(n, dn, phi), phi))/real(dfstar(n1, dlof1(n, dn, phi), phi));
frd = @(dn, n, phi) imag(dfstar(n3, dlof1(n, dn, phi), phi))/real(dfstar(n3, dlof1(n, dn, phi), phi));

%%  Make the contour plots
freqerror = handles.din.bare.error.f;
disserror = handles.din.bare.error.g;

% The following assumes that n3 = n1 or n2, so this exits if that isn't the
% case.
if n3 ~= n1 && n3 ~= n2
    set(handles.statusbar, 'string', 'Unable to make contour plot', 'BackgroundColor', 'yellow')  
    warning = warndlg('There isn''t a function to calculate the error in the contour plot for this harmonic combination. Please select one where n3 = n1 or n2.', 'Contour Error', 'modal');
    uiwait(warning);
    contourplots = [];
    return
end

shift = 2; %Sets the range over which the derivative is calculated
shifts = -shift:shift:shift;

%Sets up a matrix of all of the combinations that need to be looked at.
indecies = [unique(perms([1 2 2]), 'rows'); unique(perms([2 2 3]), 'rows')];
for i = 1:6
    k = indecies(i,1);
    l = indecies(i,2);
    m = indecies(i,3);
    df(n1) = df1 + shifts(k);
    df(n2) = df2 + shifts(l);
    dg(n3) = dg3 + shifts(m);
    
    dissratiosol(k, l, m) = dg(n3)./df(n3); %df(n3) will be equal to either df(n1) or df(n2)
    harmratiosol(k, l, m) = (n1/n2)*df(n2)./df(n1);
end

[yrd, xrd, zrd] = gradient(dissratiosol, shift); %Calculates the gradient
[yrh, xrh, zrh] = gradient(harmratiosol, shift);

errrh = [xrh(2,2,2) yrh(2,2,2) zrh(2,2,2)];
errrd = [xrd(2,2,2) yrd(2,2,2) zrd(2,2,2)];

shifts = [freqerror(n1) freqerror(n2) disserror(n3)];

rhe = sqrt(sum((errrh.*shifts).^2)); %Turns the gradient into an error based on the derfault error values
rde = sqrt(sum((errrd.*shifts).^2));

edissratio = delg(n3)./delf(n3);
eharmratio = (n1/n2)*delf(n2)./delf(n1);

dissrrange=[(1-rde)*edissratio, (1+rde)*edissratio]; % this is the range of dissipation ratios to consider (n=3)
harmrrange=[(1-rhe)*eharmratio, (1+rhe)*eharmratio]; % range of harmonic ratios to consider (delf(3)/delf(1))

dmin=handles.contourtable.Data(1,1);
dmax=handles.contourtable.Data(1,2);
phimin=handles.contourtable.Data(2,1);
phimax=handles.contourtable.Data(2,2);

n=75; % resolution of map
dplot = linspace(dmin,dmax,n);
phiplot = linspace(phimin,phimax,n);

for i = 1:n %Calculates the values of rh and rd over the contour plot
    for j = 1:n
        harmratio(j, i) = frh(dplot(i), n3, phiplot(j));
        dissratio(j, i) = frd(dplot(i), n3, phiplot(j));
    end
end

scrsize = get(0,'ScreenSize');
left=scrsize(1);
bot=scrsize(2);
width=scrsize(3);
height=scrsize(4);
contourplots.fig = figure('outerPosition',[width/2,height/2,width/2,height/2.3]);
contourplots.fig.Name = 'contourfig';
if get(handles.issim, 'value')
    comb = handles.constants.nhvals;
    n1 = comb{mnumber}(1);
    n2 = comb{mnumber}(2);
    n3 = comb{mnumber}(3);
    delf(1:2:5) = handles.deldatadata(1:3,2);
    delg(1:2:5) = handles.deldatadata(1:3,4);
    df1 = delf(n1);
    df2 = delf(n2);
    dg3 = delg(n3);
    df3 = delf(n3);
    drho = str2num(handles.propertytable.Data{1,1});
    grho = str2num(handles.propertytable.Data{1,2});
    phiout = str2num(handles.propertytable.Data{1,3});
    d1out = handles.deldatadata(1,8);
else
    comb = handles.constants.nhvals;
    n1 = comb{mnumber}(1);
    n2 = comb{mnumber}(2);
    n3 = comb{mnumber}(3);
    delf = outputs.df;
    delg = outputs.dg;
    df1 = delf(n1);
    df2 = delf(n2);
    dg3 = delg(n3);
    df3 = delf(n3);
    drho = outputs.drho;
    grho = outputs.grho;
    phiout = outputs.phi;
    d1out = outputs.d;
end
harmplotlim = [handles.contourtable.Data(3,1), handles.contourtable.Data(3,2)];
contourplots.harmplot = subplot(1,2,1);
contourf(dplot,phiplot,harmratio,linspace(harmplotlim(1), harmplotlim(2),256),'edgecolor','none');
hold on
contour(dplot,phiplot,dissratio,dissrrange,'edgecolor','black','linewidth',3,'linestyle','--')
contour(dplot,phiplot,harmratio,harmrrange,'edgecolor','black','linewidth',3)
colormap(jet(256*2));
%caxis([-.5,1.5]);
caxis(harmplotlim)
colorbar
xlabel(['d/\lambda_' num2str(n3)])
ylabel('\phi')
titlestring=strcat('r_h=', num2str(eharmratio,'%6.4f'), ': (', handles.constants.label(mnumber),')');
title(titlestring)

dissplotlim = [handles.contourtable.Data(4,1), handles.contourtable.Data(4,2)];
contourplots.dissplot = subplot(1,2,2);
contourf(dplot,phiplot,dissratio,linspace(dissplotlim(1), dissplotlim(2),256),'edgecolor','none');
hold on
contour(dplot,phiplot,dissratio,dissrrange,'edgecolor','black','linewidth',3,'linestyle','--')
contour(dplot,phiplot,harmratio,harmrrange,'edgecolor','black','linewidth',3)
colormap(jet(256));
caxis(dissplotlim);
colorbar
xlabel(['d/\lambda_' num2str(n3)])
ylabel('\phi')
titlestring=strcat('r_d=', num2str(edissratio,'%6.4f'), ': (', handles.constants.label(mnumber),')');
title(titlestring)

plot(contourplots.dissplot,d1out(n3), phiout, 'k+', 'markersize', 12, 'linewidth', 2)
plot(contourplots.harmplot,d1out(n3), phiout, 'k+', 'markersize', 12, 'linewidth', 2)

linkaxes([contourplots.dissplot, contourplots.harmplot],'xy')
set(0,'CurrentFigure', maingui); %Return focus (and crosshairs) to main gui window.

function handles = writecontourmap
dl = linspace(0,1,250);
phi = linspace(0,90,250);

for i = 1:length(dl)
    for j = 1:length(phi)
        dfstar(i,j) = delfstar(dl(i), phi(j));
    end
end

min_c = -3;
max_c = 3;
min_b = 0;
max_b = 1;

scrsize = get(0,'ScreenSize');
left=scrsize(1);
bot=scrsize(2);
width=scrsize(3);
height=scrsize(4);
handles.responsemap = figure('outerPosition',[width/2,height/2,width/2,height/2.3]);
set(gcf,'defaultaxesbox','on')

handles.plot1 = subplot(1,2,1);
hold on
[C h1] = contourf(dl, phi, real(dfstar)',[logspace(1,1.24,250)-15]);
xlabel('d/\lambda_n')
ylabel('\phi_n')
colorbar
title('{\Delta}f_n/{\Delta}f_{sn}')
caxis([min_c, max_c])
set(h1,'LineStyle','none');
%set(findall(gcf,'-property','FontSize'),'fontsize',20)

handles.plot2 = subplot(1,2,2);
hold on
[C h1] = contourf(dl, phi, imag(dfstar)',[0 logspace(-3,.7,250)]);
%[C h1] = contourf(dl, phi, imag(napprox)',[linspace(0,4,10)]);
xlabel('d/\lambda_n')
ylabel('\phi_n')
colorbar
title('\Delta\Gamma_n/{\Delta}f_{sn}')
caxis([min_b, max_b])
set(h1,'LineStyle','none');
%set(findall(gcf,'-property','FontSize'),'fontsize',20)
try
    load('onelayerguicolormaps.mat', 'mycmap')
    colormap(mycmap)
catch Err
    if strcmp(Err.identifier,'MATLAB:load:couldNotReadFile')
        warning('The custom colormap file ''custcolormaps.mat'' was not found. For optimal viewing, please load the file.')
    end
end
linkaxes([handles.plot1 handles.plot2])
savefig('contourmap.fig')

function solveall_Callback(hObject, eventdata, handles)
checksolveharms(hObject, handles);
handles.activenh=getactivenh(handles);
mnumber = getnhvals(handles); %Determine harmonics to solve for
nhvals = handles.constants.nhvals; %Make the variable shorter to type
set(handles.plotsolutioncontours, 'value', 0); %turn off plot contours for solve all

struct2var(handles.din);
% Solves only for the points in the window.
limits = handles.axes1.XLim;
[~, withinidx] = find(qcmt > limits(1) & qcmt < limits(2));
handles.datapoints = length(withinidx);

for k=1:length(withinidx); %Goes through the indexes of the relevant points
    idx = withinidx(k);
    for nh=handles.activenh.on
        handles.deldatadata(nhtoi(nh),1) = cleandelf(idx,nh);
        handles.deldatadata(nhtoi(nh),3) = cleandelg(idx,nh);
    end
    mnumber = getnhvals(handles);
    [contourplots, handles] = redrawtable(hObject, handles, mnumber);
    for m = mnumber
        [outputs{k,m}, handles] = findsolution(hObject, eventdata, handles, qcmt(idx), m);
        if outputs{k,m}.d(1) < 0 || outputs{k,m}.phi>90 || outputs{k,m}.error == 0 || outputs{k,m}.error == -2
            %             answer = questdlg({['Failed to find solution at time ' num2str(qcmt(k)) ...
            %                 ' and ratio ' num2str(nhvals{m}(1)) ':' num2str(nhvals{m}(2))...
            %                 ',' num2str(nhvals{m}(3))],[],['Do you want to change the input values?']},...
            %                 'Yes', 'No')
            %             switch answer
            %                 case 'Yes'
            %                     keystroke = 0;
            %
            %                     while keystroke ~= 'w';
            %                        waitforbuttonpress;
            %                        keystroke = get(gcf, 'CurrentCharacter')
            %                     end
            %                     outputs{k,m}=findsolution(hObject, eventdata, handles,qcmt(k),nhvals{m});
            %                 case 'No'
            %                 case 'Cancel'
            %                     return
            %             end
        end
    end
    if mod(k,100)==0
        outputstructure = reformatdata(handles, outputs, 0);
        handles.dout = outputstructure;
        handles = plotvalues(hObject,handles);
        set(handles.statusbar, 'string', [num2str(k/handles.datapoints*100, 2) '% done'], 'BackgroundColor', 'default')
        drawnow %Update the visual on the gui
    end
end
outputstructure = reformatdata(handles, outputs, 1);
handles.dout = outputstructure;
handles = plotvalues(hObject,handles);
handles.contourplos = contourplots;
guidata(hObject,handles);

function n1_Callback(hObject,eventdata,handles)
handles.activenh = getactivenh(handles);
guidata(hObject,handles)

function n3_Callback(hObject,eventdata,handles)
handles.activenh = getactivenh(handles);
guidata(hObject,handles)

function n5_Callback(hObject,eventdata,handles)
handles.activenh = getactivenh(handles);
guidata(hObject,handles)

function n7_Callback(hObject, eventdata, handles)
handles.activenh = getactivenh(handles);
guidata(hObject,handles)

function activenh = getactivenh(handles)
% Determines the harmonic buttons that are selected
activenh.on=[];
activenh.off=[];
for nh=[1,3,5,7]
    if get(handles.(['n',num2str(nh)]),'value')
        activenh.on=[activenh.on,nh];
    else
        activenh.off=[activenh.off,nh];
    end
end

function mnumber = getnhvals(handles)
% Determines which harmonics should be calculated
mnumber = [];
for i = 1:6
    % Reads in the checkmarks and records the numbers of the checked ones
    if get(handles.(['nhset',num2str(i),'select']),'value') == 1 %If box is checked
        m = get(handles.(['nhset',num2str(i),'type']), 'value'); %Get harmonic requested
        mnumber = [mnumber m];
    end
end
mnumber = sort(unique(mnumber)); % Remove duplicates and put in sensible order
% Uses the first box if nothing is selected.
if isempty(mnumber)
    warning('Using default harmonic since none was selected')
    set(handles.nhset1select, 'value', 1)
    mnumber = get(handles.nhset1type, 'value');
end

function struct2var(s)
if nargin < 1
    error('struct2var:invalidInput','No input structure found')
elseif nargin > 1
    error('struct2var:invalidInput','Too many inputs')
elseif ~isstruct(s)
    error('struct2var:invalidInput','Input must be a structure data type')
end

%[~,c] = size(s);
names = fieldnames(s);

for i=1:length(names)
    assignin('caller',names{i},s.(names{i}))
end

function log_Callback(hObject, eventdata, handles)
if get(handles.log,'value')
    set([handles.axes1 handles.axes2],'xscale','log')
    axistype = 'log';
else
    set([handles.axes1 handles.axes2],'xscale','linear')
    axistype = 'linear';
end
set(findall([getappdata(0, 'outputplot') getappdata(0, 'inputplot') handles.axes1 handles.axes2],...
    'type','axes'), 'XScale', axistype, 'Xlim', [-Inf Inf])

function CloseRequestFcn(hObject, eventdata, handles)
% %hObject is the handle of the object generating the callback (the source of the event)
% %evnt is the The event data structure (can be empty for some callbacks)
delete(hObject)
try
detete(figure(995))
end
% selection = questdlg('Do you want to save the file?',...
%     'Close Request Function',...
%     'Yes','No','Cancel','Yes');
% switch selection,
%     case 'Yes',
%         handles=guidata(hObject);
%         if isfield(handles,'din')
%             guidata(hObject,handles)
%             hgsave('onelayergui6.fig')
%             delete(hObject)
%         else
%             fprintf(1,'handles.din does not exist')
%             return
%         end
%     case 'No'
%         delete(hObject)
%     case 'Cancel'
%         return
% end

function plot_Callback(hObject, eventdata, handles)
set(0,'defaultaxesfontsize',10)
set(0,'defaultlinemarkersize',8)
try
    set(0, 'defaulterrorbarmarkersize', 8)
    set(0,'defaulterrorbarlinewidth',1)
catch
end
set(0,'defaultlinelinewidth',1)
handles = plotvalues(hObject,handles);
guidata(hObject,handles)

function Clear_Callback(hObject, eventdata, handles)
handles.n=0;
handles.dout.time={};
handles.dout.df={};
handles.dout.dg={};
handles.dout.drho={};
handles.dout.grho={};
handles.dout.phi={};
handles.dout.dfcalc={};
handles.dout.dgcalc={};
guidata(hObject,handles)

function buildcondfile_Callback(hObject, eventdata, handles)
% Pressing this button causes the gui to go through all of the spectra
% associated with the file and save them in the _cond file for later
% reading by condfig.

points = handles.cond.points;
filename = handles.cond.filename;
qcmt = handles.din.qcmt;

% Load raw spectra files if it is a mat file.
if handles.filetype == 2
    rawspectras = load([handles.din.qcmpath handles.din.filebase '_raw_spectras.mat']);
end

%Check for an existing _cond file, and determine if it need to be updated
if exist([handles.din.qcmpath handles.din.filebase '_cond.mat'], 'file') == 2
    cond = load([handles.din.qcmpath handles.din.filebase '_cond.mat'])
    if max(cond.time) == max(qcmt)
        set(handles.statusbar, 'string', 'The conductance file appears to be up to date', 'BackgroundColor', 'green')
        return
    else
        maximptime = max(cond.time); %So only need to add new times
    end
else
    maximptime = 0; %Set to 0 if none found (so it will add all)
    set(handles.statusbar, 'string', 'The cond file will have to be built from scratch. This could take a while.', 'BackgroundColor', 'default')
end

pointstoadd = find(points(:,1)>maximptime); %Finds spectra from later times
for i = pointstoadd'
    if points(i,3) ~= 0
        % Import spectra and save to array
        harm = itonh(points(i,2)); %harmonic in 1,3,5
        idx = points(i,3); %index in main data
        time = points(i,1);
        cond.time(idx) = time;
        if handles.filetype == 1
            condstruct = importdata(char(filename(i)));
            cond.spectra{idx,harm} = condstruct.data;
        elseif handles.filetype == 2
            cond.spectra{idx,harm} = rawspectras.(char(filename(i)));
        end
    end
end

save([handles.din.qcmpath handles.din.filebase '_cond.mat'], '-struct', 'cond');
set(handles.statusbar, 'string', 'New spectra added to file!', 'BackgroundColor', 'green')

function outstruct = reformatdata(handles, outputs, ask)
doerror = get(handles.calcerror, 'value');
% each value of k corresponds to one time point selected from the plot
% variables with 'p' in the name are the ones we use to create the plots
try
    dfp=handles.din.cleandelf;
    dgp=handles.din.cleandelg;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        set(handles.statusbar, 'string', 'Unable to replot data. Please reload file', 'BackgroundColor', 'red')
        changeqcmfile_Callback(hObject, 1, handles, handles.din.qcmfile, handles.din.qcmpath);
    end
end

mnumber = getnhvals(handles);

timeinp=handles.din.qcmt;
try %Try statement in case no data is loaded
    todo = size(outputs);
catch Err
    if strcmp(Err.identifier,'MATLAB:nonExistentField')
        if isfield(handles.saveplot, 'grhop')
            loadsavedplot_Callback(hObject, 1, handles)
            return
        else
            set(handles.statusbar, 'string', 'No data to plot. Please solve for values before plotting.', 'BackgroundColor', 'yellow')
            return
        end
    else
        rethrow(Err)
    end
end
refG = itonh(get(handles.modulustoplot, 'value'));
for k=1:todo(1)
    for m = mnumber
        struct2var(outputs{k,m})
        outstruct.grhop(k,m)=grho(refG);
        outstruct.drhop(k,m)=drho(1);
        outstruct.phip(k,m)=phi(1);
        outstruct.timep(k,m)=time;
        if doerror
            try
                outstruct.grhoep(k,m) = grhoe;
                outstruct.drhoep(k,m) = drhoe;
                outstruct.phiep(k,m) = phie;
            catch Err
                if strcmp(Err.identifier, 'MATLAB:UndefinedFunction')
                    outstruct.grhoep(k,m) = NaN;
                    outstruct.drhoep(k,m) = NaN;
                    outstruct.phiep(k,m) = NaN;
                else
                    rethrow(Err)
                end
            end
        else % If no error data was collected
            outstruct.grhoep(k,m) = NaN;
            outstruct.drhoep(k,m) = NaN;
            outstruct.phiep(k,m) = NaN;
        end
        
        for nh=handles.activenh.on
            outstruct.dfcalcp(k,m,nh)=dfcalc(nh);
            outstruct.dgcalcp(k,m,nh)=dgcalc(nh);
        end
    end
end
outstruct.nhvals = {handles.constants.nhvals{mnumber}};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves data to _data file, including raw data.
if get(handles.autosolve,'value') && ask == 1 %If something was actually solved for
    choice = questdlg(['Do you want to save these results?'], 'Save plotted results',...
        'Yes', 'No', 'Yes');
    switch choice
        case 'Yes'
            append = get(handles.qcmfile,'string');
            append = [handles.din.qcmpath append(1:end-4)];
            
            data = struct('dgcalcp', outstruct.dgcalcp, 'dfcalcp', outstruct.dfcalcp, 'dgp', dgp, 'dfp', dfp,...
                'timep', outstruct.timep, 'grhop', outstruct.grhop, 'drhop', ...
                outstruct.drhop, 'phip', outstruct.phip, 'grhoep', outstruct.grhoep,...
                'phiep', outstruct.phiep, 'drhoep', outstruct.drhoep,...
                'delf', handles.din.cleandelf, 'delg', handles.din.cleandelg,...
                'absf', handles.din.absf, 'absg', handles.din.absg,...
                'time', handles.din.qcmt, 'error', handles.din.bare.error,...
                'offset', handles.din.bare.offset, 'nhvals', {outstruct.nhvals}, 'refG', refG);
            save([append '_data.mat'], '-struct', 'data');
        case 'No'
            return
    end
end

function handles = plotvalues(hObject, handles, inputfigure, outputfigure, inputaxes, outputaxes)
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')
mnumber = getnhvals(handles);
handles.activenh = getactivenh(handles);

try %Sometimes these values aren't saved--the try fixes it if not
    xlabeltext = handles.plotting.xlabeltext;
    timefactor = handles.plotting.timefactor;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        handles = xlabeltext_Callback(hObject, mnumber, handles); %So I've included mnumber here instead of eventdata as a dummy variable because I don't use eventdata. 
        xlabeltext = handles.plotting.xlabeltext;
        timefactor = handles.plotting.timefactor;
    else
        rethrow(Err)
    end
end

doerror = get(handles.calcerror, 'value');
if handles.saveplot.plot == 1
    grhop = handles.saveplot.grhop;
    drhop = handles.saveplot.drhop;
    phip = handles.saveplot.phip;
    timep = handles.saveplot.timep;
    try
        grhoep = handles.saveplot.grhoep;
        drhoep = handles.saveplot.drhoep;
        phiep = handles.saveplot.phiep;
    catch
        doerror = 0;
    end
    dfcalcp = handles.saveplot.dfcalcp;
    dgcalcp = handles.saveplot.dgcalcp;
    dfp = handles.saveplot.delf;
    dgp = handles.saveplot.delg;
    timeinp = handles.saveplot.time;
else
    grhop = handles.dout.grhop;
    drhop = handles.dout.drhop;
    phip = handles.dout.phip;
    timep = handles.dout.timep;
    try
        grhoep = handles.dout.grhoep;
        drhoep = handles.dout.drhoep;
        phiep = handles.dout.phiep;
    catch
        doerror = 0;
    end
    dfcalcp = handles.dout.dfcalcp;
    dgcalcp = handles.dout.dgcalcp;
    dfp = handles.din.cleandelf;
    dgp = handles.din.cleandelg;
    timeinp = handles.din.qcmt;
end
% now we compare simulated and actual frequency and dissipation shifts
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinemarkersize',12)
try
    set(0, 'defaulterrorbarmarkersize', 12)
    set(0,'defaulterrorbarlinewidth',1.5)
catch
end
set(0,'defaultlinelinewidth',1.5)
% delete all previous plots (not necessary, but keeps things a bit cleaner--also doesn't work in 2015a)
% allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% delete(allPlots);
try
    if isvalid(getappdata(0, 'outputplot')) %Close previous plots--works in 2015a
        close(getappdata(0, 'outputplot'));
        close(getappdata(0, 'inputplot'));
    end
catch
    % If it doesn't work, that's ok. They just won't get closed. But I
    % don't want an error.
end

% determine the screensize so that the two figures are split vertically on
% the screen
scrsize = get(0,'ScreenSize');
left=scrsize(1);
bot=scrsize(2);
width=scrsize(3);
height=scrsize(4);

inputplot=figure('outerPosition',[width/2,height/2,width/2,height/2.3]);
set(gcf,'defaultaxesbox','on')
setappdata(0, 'inputplot', inputplot);
inputaxes(1)=subplot(1,2,1);
inputaxes(2)=subplot(1,2,2);

outputplot=figure('outerPosition',[width/2,bot,width/2,height/2.3]);
set(gcf,'defaultaxesbox','on')
setappdata(0, 'outputplot', outputplot)
outputaxes(1)=subplot(1,3,1);
outputaxes(2)=subplot(1,3,2);
outputaxes(3)=subplot(1,3,3);

symbols{1}='+';
symbols{2}='o';
symbols{3}='s';
symbols{4}='x';
symbols{5}='^';
symbols{6}='d';
colors{1}=[0, 0, 0];
colors{2}=[0, 0, 0];
colors{3}=[0, 0, 0];
colors{4}=[1, 0, 0];
colors{5}=[1, 0, 0];
colors{6}=[1, 0, 0];
linestyles = {'-' '' '-.' '' '--' '' ':'};
datalinecolor{1}=[1,0,0];
datalinecolor{3}=[0,0.5,0];
datalinecolor{5}=[0,0,1];
datalinecolor{7}=[0,0,0];

figure(inputplot)
set(gcf,'defaultaxesbox','on')
axes(inputaxes(1))
hold on
% first plot the measured frequency shifts for the three harmonics
index=0;
legendentries=[];

for nh=handles.activenh.on
    index=index+1;
    toplot = ~isnan(dfp(:,nh));
    try
        plots(index)=plot(timeinp(toplot)./timefactor,dfp(toplot,nh)/nh,linestyles{nh},'color',datalinecolor{nh});
        legendtext{index}=['n=',num2str(nh)];
    catch Err
        if strcmp(Err.identifier, 'MATLAB:badRectangle')
            set(handles.statusbar, 'string', 'No data to plot for the 5th harmoic', 'BackgroundColor', 'yellow')
            handles.activenh.on = handles.activenh.on(handles.activenh.on~=5);
            index = index-1;
        end
    end
end

harmcount = 1;
% now we plot the calculated values
for m = mnumber
    legendentries=[legendentries,index+1];
    for nh=handles.activenh.on
        index=index+1;
        plots(index)=plot(timep(:,m)./timefactor,dfcalcp(:,m,nh)/nh,symbols{harmcount},'color',datalinecolor{nh});
        legendtext{index}=handles.constants.label{m};
    end
    harmcount = harmcount + 1;
end
xlabel(xlabeltext)
ylabel('\Deltaf_{n}/n (Hz)')
legendin(1)=legend(plots(legendentries),legendtext(legendentries),'location','best');

axes(inputaxes(2))
% repeat the plot for the dissipation
% first plot the measured data for the three harmonics
index=0;
for nh=handles.activenh.on
    index=index+1;
    toplot = ~isnan(dfp(:,nh));
    plots(index)=semilogy(timeinp(toplot)./timefactor,dgp(toplot,nh),linestyles{nh},'color',datalinecolor{nh});
    legendtext{index}=['n=',num2str(nh)];
    hold on
end

harmnumb = 1;
% now we plot the calculated values
for m = mnumber
    for nh=handles.activenh.on
        index=index+1;
        plots(index)=semilogy(timep(:,m)./timefactor,dgcalcp(:,m,nh),symbols{harmnumb},'color',datalinecolor{nh});
        legendtext{index}=handles.constants.label{m};
    end
    harmnumb = harmnumb + 1;
end
if get(handles.loggamma, 'value')
    set(inputaxes(2), 'yscale', 'log');
else
    set(inputaxes(2), 'yscale', 'linear');
end

xlabel(xlabeltext)
ylabel('\Delta\Gamma_{n} (Hz)')
legendin(2)=legend(plots(legendentries),legendtext(legendentries),'location','best');
% set(legendin(2),'edgecolor',[1 1 1])

%  Now plot the property figure
symbols{1}='+';
symbols{2}='o';
symbols{3}='s';
symbols{4}='x';
symbols{5}='^';
symbols{6}='d';
colors{1}=[1, 0, 0];
colors{2}=[0, 0.5, 0];
colors{3}=[0, 0, 1];
colors{4}=[0, 0, 0];
colors{5}=[1, 0, 0];
colors{6}=[0, 0.5, 0];

figure(outputplot)
set(gcf,'defaultaxesbox','on')
axes(outputaxes(1))
hold on
harmnumb = 1;
if doerror
    for m = mnumber
        errorbar(timep(:,m)./timefactor,drhop(:,m),drhoep(:,m),symbols{harmnumb},'color',colors{harmnumb})
        harmnumb = harmnumb + 1;
    end
else
    for m = mnumber
        plot(timep(:,m)./timefactor,drhop(:,m),symbols{harmnumb},'color',colors{harmnumb})
        harmnumb = harmnumb + 1;
    end
end
legendout(1)=legend(handles.constants.label{mnumber},'location','best');
xlabel(xlabeltext)
ylabel('d\rho (g/m^{2})')

axes(outputaxes(2))
hold on
harmnumb = 1;
if doerror
    for m = mnumber
        errorbar(timep(:,m)./timefactor,grhop(:,m),grhoep(:,m),symbols{harmnumb},'color',colors{harmnumb})
        harmnumb = harmnumb + 1;
    end
else
    for m = mnumber
        plot(timep(:,m)./timefactor,grhop(:,m),symbols{harmnumb},'color',colors{harmnumb})
        harmnumb = harmnumb + 1;
    end
end
legendout(2)=legend(handles.constants.label{mnumber}, 'location','best');
xlabel(xlabeltext)
ylabel(['|G_' num2str(itonh(get(handles.modulustoplot, 'value'))) '^{*}|\rho (Pa-g/cm^{3})'])
if get(handles.logG, 'value')
    set(outputaxes(2), 'yscale', 'log');
else
    set(outputaxes(2), 'yscale', 'linear');
end

axes(outputaxes(3))
hold on
harmnumb = 1;
if doerror
    for m = mnumber
        errorbar(timep(:,m)./timefactor,phip(:,m),phiep(:,m),symbols{harmnumb},'color',colors{harmnumb})
        harmnumb = harmnumb + 1;
    end
else
    for m = mnumber
        plot(timep(:,m)./timefactor,phip(:,m),symbols{harmnumb},'color',colors{harmnumb})
        harmnumb = harmnumb + 1;
    end
end
legendout(3)=legend(handles.constants.label{mnumber},'location','best');
xlabel(xlabeltext)
ylabel('\phi (deg.)')

if get(handles.log, 'value') == 1
    set([inputaxes, outputaxes], 'xscale', 'log')
end

% link all the x axes together, so we don't have to change all the plots
% individually.
% linkaxes uses the limits of the first plot given, in this case,
% outputaxes.
linkaxes([outputaxes(1:3) inputaxes(1:2)],'x')
set(outputaxes(1), 'xlim', [-Inf Inf])

function write_Callback(hObject,eventdata,handles)
writeplots(hObject,handles)

function writeplots(hObject,handles)
% print the output files
set(handles.inputplot,'paperposition',[0 0 8 3.2])
set(handles.inputplot,'papersize',[8 3.2])
set(handles.outputplot,'paperposition',[0 0 12 3.2])
set(handles.outputplot,'papersize',[12 3.2])
set(handles.inputplot,'Resizefcn','')
set(handles.outputplot,'Resizefcn','')

print(handles.outputplot,[handles.din.qcmpath, handles.din.filebase,'_out.eps'],'-depsc2')
saveas(handles.outputplot,[handles.din.qcmpath,handles.din.filebase,'_out.fig'])

print(handles.inputplot, [handles.din.qcmpath,handles.din.filebase,'_in.eps'],'-depsc2')
saveas(handles.inputplot,[handles.din.qcmpath,handles.din.filebase,'_in.fig'])

function cleaned_Callback(hObject, eventdata, handles)
% The value of this function is called to determine things, but when it is
% changed the only thing necessary is to replot the data. Switching it to
% the cleaned data is taken care of by using "get" in the plotraw function.
plotraw(hObject, eventdata, handles)

function savecleaned_Callback(hObject, eventdata, handles)
% Changes variables in the _data file, which contains the "clean" values.
m = matfile([handles.din.qcmpath handles.din.filebase '_data.mat'],'Writable',true);
m.delf = handles.din.cleandelf;
m.delg = handles.din.cleandelg;
m.absf = handles.din.absf;
m.absg = handles.din.absg;
m.time = handles.din.qcmt;
m.offset = handles.din.bare.offset;
m.error = handles.din.bare.error;

function handles = loadcleaned_Callback(hObject, eventdata, handles)
% Loads the "clean" variables from the _data file (delf and delg).
try
    load([handles.din.qcmpath handles.din.filebase '_data.mat'], 'absf', 'absg')
catch Err
    if strcmp(Err.identifier, 'MATLAB:load:couldNotReadFile')
        set(handles.statusbar, 'string', 'There was no previously saved file to load', 'BackgroundColor', 'yellow')
        return
    elseif strcmp(Err.identifier, 'MATLAB:nonExistentField')
        set(handles.statusbar, 'string', 'The previous file seems to have been lost. Please load a new one.', 'BackgroundColor', 'red')
        return
    else
        rethrow(Err)
    end
end
%If there is new data since the last save, it should be included. So this
%starts with all of the data as it is currently

% Gets the length of the saved data
[datarows datacolumns] = size(absf); %Gets the size of the saved data
%Overwrites the times for which there is previously saved data with that
%data.
if length(handles.din.qcmt) < datarows %If the previously saved data is longer than the just loaded data
    datarows = length(handles.qcmt);
    set(handles.statusbar, 'string', 'Warning: Previously saved data is longer than current data', 'BackgroundColor', 'yellow')
end
handles.din.absf(1:datarows,1:datacolumns) = absf; %Overwrites with saved data
handles.din.absg(1:datarows,1:datacolumns) = absg; %Overwrites with saved data

for nh=1:2:datacolumns
    fref=['f',num2str(nh)];
    gref=['g',num2str(nh)];
    handles.din.cleandelf(:,nh)=(handles.din.absf(:,nh)-handles.din.bare.offset.(fref));
    handles.din.cleandelg(:,nh)=(handles.din.absg(:,nh)-handles.din.bare.offset.(gref));
end

guidata(hObject,handles); %Save changes to handles structure
plotraw(hObject, eventdata, handles) %Plot updated data

function opencond_Callback(hObject, eventdata, handles)
% This function opens the second program, condfig, which displays
% conductance data.
try
    if isempty(handles.cond)
        set(handles.statusbar, 'string', 'You cannot view the spectra since there aren''t any', 'BackgroundColor', 'red')
        return
    end
catch
    set(handles.statusbar, 'string', 'cond variable doesn''t exist', 'BackgroundColor', 'red')
end
set(handles.statusbar, 'string', 'Loading conductance data...', 'BackgroundColor', 'default')
drawnow
% Updates (or creates) _cond file for use by condfig.
buildcondfile_Callback(hObject, eventdata, handles);

% To facilitate getting conductance data for a specific point,
% datacursormode is turned on, with a special function that
% sends the time to condfig.
datacursormode on
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn, hObject, handles})

% This part actually runs condfig. When the program is closed,
% any changes to the cleaned figure file is saved. Which means
% that this function "stays open" until the window is closed.
handles.din = condfig6('onelayerguiwithcond', handles.figure);
guidata(hObject,handles);
datacursormode off
%The potentially changed data is plotted.
plotraw(hObject, eventdata, handles)
set(handles.statusbar, 'string', 'The data has been transferred', 'BackgroundColor', 'green')

function output_txt = myupdatefcn(~,event_obj,~, handles)
% This is the function that runs when datacursormode is employed. The
% output output-txt is what appears in the box.

%Determines output box--this is the usual function of datacursor, modified
%to know what the x axis actually is.
pos = get(event_obj,'Position');
if handles.plotting.timefactor == 1 %Show time and f/g values of the selected point
    output_txt = {['Time: ', num2str(pos(1),'%6.2f'), ' min'],...
    ['y: ',num2str(pos(2),'%6.2f')]};
else %If the time is not in minutes, show the time in minutes and whatever the units are
    output_txt = {['Time: ', num2str(pos(1),'%6.2f'), ' ', handles.plotting.unit],...
    ['Time: ', num2str(pos(1)*handles.plotting.timefactor,'%6.2f'), ' min'],...
    ['y: ',num2str(pos(2),'%6.2f')]};
end

% The following determines the time requested and sends it to the
% condfig window if that is the only other window open. Otherwise I'm not
% sure how to find the handle to the right window.

% Finds handles to all the windows
getting = get(0, 'children');


%If there are only two figure windows open, the second one is the condfig
%window. This then writes the time from above to the time box in condfig
% and runs the "view_Callback" function.

for i = 1:length(getting) %Look for the figure window which is the condfig window
    if strcmp(getting(i).Name, 'condfig')
        condfighandle = getting(i);
        break
    end
end

condhandles = guidata(getting(2));
if isfield(condhandles, 'currentpoint')
    set(condhandles.currentpoint, 'string', num2str(pos(1)*handles.plotting.timefactor, '%0.3f'))
    %This "imports" the function into this gui
    showcurves = get(condhandles.view, 'Callback');
    %Now I can call the function using the arguments in showcurves,
    %which can be seen if you view showcurves--note that it treats
    %handles as a function of hObject, so there are actually only two
    %inputs needed, not three.
    showcurves(condfighandle, event_obj)
end

datacursormode off

function checksolveharms(hObject, handles)
% Makes sure that the harmonics that are selected to be tried are ones for
% which the data was imported. Prohibits a [155] calculation if n=5 isn't
% checked.
try
    nhsel = handles.activenh.on;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        handles.activenh = getactivenh(handles);
        nhsel = handles.activenh.on;
    end
end
%This determines the required harmonics for each nhvals combination. It's
%written as a cellfun so that it doesn't need to be changed if
%handles.constants.nhvals is changed.
combinations = cellfun(@(x) sort(unique(x)), handles.constants.nhvals, 'UniformOutput', false);
if length(nhsel) < 4
    for i = 1:6
        if get(handles.(['nhset' num2str(i) 'select']), 'value') == 1
            mnumber = get(handles.(['nhset' num2str(i) 'type']), 'value');
            overlap = sum(ismember(nhsel, combinations{mnumber}));
            if overlap ~= length(combinations{mnumber}) %If one of the values is missing
                set(handles.(['nhset' num2str(i) 'select']), 'value', 0)
            end
        end
    end
end

function [drhoe grhoe phie] = finderror(hObject, handles, mnumber)
nh = handles.constants.nhvals{mnumber};
% Calculates the error in a measurement.
try
    nhs=handles.activenh.on;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        activenh = getactivenh(handles);
        nhs=activenh.on;
    else
        rethrow(Err)
    end
end

refG = itonh(get(handles.modulustoplot, 'value'));

if get(handles.issim, 'value') %Use the editable "calc" values
    dfi(nhs) = handles.deldatadata(ceil(nhs/2),2);
    dgi(nhs) = handles.deldatadata(ceil(nhs/2),4);
else %Use the values stored in the function
    dfi(nhs) = handles.deldatadata(ceil(nhs/2),1);
    dgi(nhs) = handles.deldatadata(ceil(nhs/2),3);
end

errorsf = handles.din.bare.error.f;
errorsg = handles.din.bare.error.g;

shift = 2;
shifts = -shift:shift:shift;

%This creates a matrix of all of the shifts necessary
indecies = [unique(perms([1 2 2]), 'rows'); unique(perms([2 2 3]), 'rows')];
for i = 1:6
    k = indecies(i,1);
    l = indecies(i,2);
    m = indecies(i,3);
    df(nh(1)) = dfi(nh(1)) + shifts(k);
    df(nh(2)) = dfi(nh(2)) + shifts(l);
    dg(nh(3)) = dgi(nh(3)) + shifts(m);
    
    outputs = solvedfstar(df, dg, mnumber, hObject, handles, refG);
    phis = outputs.phi;
    drhos = outputs.drho;
    grhos = outputs.grho(refG); %Unit conversions have already been done
    if outputs.error == 1 && grhos < 10^10 && phis > 0
        phi(k, l, m) = phis;
        drho(k, l, m) = drhos;
        grho(k, l, m) = grhos;
    else
        phi(k, l, m) = NaN;
        drho(k, l, m) = NaN;
        grho(k, l, m) = NaN;
    end
end

[yphi, xphi, zphi] = gradient(phi, shift);
[ydrho, xdrho, zdrho] = gradient(drho, shift);
[ygrho, xgrho, zgrho] = gradient(grho, shift);

errphi = [xphi(2,2,2) yphi(2,2,2) zphi(2,2,2)];
errdrho = [xdrho(2,2,2) ydrho(2,2,2) zdrho(2,2,2)];
errgrho = [xgrho(2,2,2) ygrho(2,2,2) zgrho(2,2,2)];
error = [errphi; errdrho; errgrho];

shifts = [errorsf(nh(1)) errorsf(nh(2)) errorsg(nh(3))];

phie = sqrt(sum((error(1,:).*shifts).^2));
drhoe = sqrt(sum((error(2,:).*shifts).^2));
grhoe = sqrt(sum((error(3,:).*shifts).^2));

pphi = phie/str2num(handles.propertytable.Data{1,3})*100;
pdrho = drhoe/str2num(handles.propertytable.Data{1,1})*100;
pgrho = grhoe/str2num(handles.propertytable.Data{1,2})*100;
handles.propertytable.Data{2,1} = num2str(pdrho, '%8.2f');
handles.propertytable.Data{2,2} = num2str(pgrho, '%8.2f');
handles.propertytable.Data{2,3} = num2str(pphi, '%8.2f');

function [outputs, handles] = solvedfstar(dfi, dgi, mnumber, hObject, handles, refn)
%SolveQCM takes the frequency and dissipation inputs, along with the
%reference angles and guesses about the properties of the film, to
%calculate the actual properties of the film.
nhi = handles.constants.nhvals{mnumber};

%This if statement determines if one of the points used in the evaluation
%was not collected (turned into NaN). If there is a NaN, this turned that
%data point into a whole set of NaN's so it won't plot.
if sum(isnan(dfi(nhi(1:3)))) > 0 %That is, if one of the values is NaN
    outputs = blankoutput();
    outputs.error = 2.5;
    handles.propertytable.Data(1,5)={'not solved'};
    
else
    f1=handles.constants.f1;
    zq=handles.constants.zq;
    dissratio=dgi(nhi(3))/dfi(nhi(3));
    harmratio=nhi(1)*dfi(nhi(2))/(nhi(2)*dfi(nhi(1)));
    % take initial guesses from simulation parameters
    phi=str2double(handles.propertytable.Data(1,3));
    grho=1000*str2double(handles.propertytable.Data(1,2)); %Pa kg/m^3
    
    drho = 0.001*str2double(handles.propertytable.Data(1,1)); %kg/m^2
    lrho = lambdarhof(refn, refn, grho, phi); %for n=refn, kg/m^2
    x0(1) = phi;
    x0(2) = drho/lrho; %d/lambda
    
    ftosolve = @(x) fstar(x, nhi, harmratio, dissratio, refn)*10;
    try
        options = optimset('Display','off');
        [x,fval, exitflag] = fsolve(ftosolve,x0,options); % Call solver
    catch MExc
        if isnan(harmratio) || isnan(dissratio) % Happens if df or dg is 0
            exitflag = 2.5;
        else
            % Sometimes the fields aren't properly loaded. This tries again.
            if strcmp(MExc.identifier,  'optim:trustnleqn:UsrObjUndefAtX0')
                redrawtable(hObject, handles, mnumber)
                x0(1) = str2double(handles.propertytable.Data(1,3));
                x0(2) = str2double(handles.deldatatable.Data{1,8});
                options = optimset('Display','off');
                try
                    [x,fval, exitflag] = fsolve(ftosolve,x0,options); %Try to call solver again.
                catch Err
                    %This catch attempts to fix d/lambda values that are out of
                    %whack because the mass has increased too far.
                    if strcmp(MExc.identifier,  'optim:trustnleqn:UsrObjUndefAtX0')
                        if x0(2) == 0 || isnan(x0(2))
                            x0(2) = .15;
                        end
                        try
                            [x,fval, exitflag] = fsolve(ftosolve,x0,options);
                        catch Err3
                            exitflag = 1.5;
                        end
                    else
                        rethrow(Err)
                    end
                end
            else
                rethrow(MExc)
            end
        end
    end
        
    %If the function did not solve satisfactorily (exitflag ~= 1), set data
    %to NaNs also.
    if exitflag ~= 1
        outputs = blankoutput();
        outputs.error = exitflag;        
        % Also discard data if the phi value is impossible.
    elseif x(1) > 90
        outputs = blankoutput();    
        outputs.error = 3.5;
    else %Exit flag was 1 -- properly solved
        phi=x(1); %#ok<*SAGROW>
        dref=x(2); %d/lambda with respect to refG
        
        refnh = nhi(1); %This ref needs to be a harmonic for which values were used, 
        % but it may not be the same as dref
        drefnh = dlf(refn, refnh, dref, phi);
        drho = (dfi(refnh)/real(delfstar(drefnh,phi)))*zq/(2*refnh*f1^2); %calculate reference drho based on n2
        grho(refn) = ((drho/dref)*refn*f1*cosd(phi/2))^2;
              
        outputs = calculatevalues(handles, drho*1000, grho(refn)*0.001, phi, refn);
        
%         decaylength(1:2:7) = outputs.lrho(1:2:7)*(1/(2*pi))*cotd(phi/2); %Calculate the decay length
%         maxnh = max(nhi(1:3)); %Get maximum harmonic in calculation
%         outputs.drho(outputs.drho > decaylength(maxnh)*0.9) = NaN;
        
        % add values to the output structure
        outputs.df=dfi;
        outputs.dg=dgi;
        outputs.error = exitflag;
           
        assignin('base','solveouptus',outputs);
    end
end

function calcerror_Callback(hObject,eventdata,handles)
handles.propertytable.Data{2,1} = '-';
handles.propertytable.Data{2,2} = '-';
handles.propertytable.Data{2,3} = '-';

function maps_Callback(hObject, eventdata, handles)
try
    responsemap = openfig('contourmap.fig');
catch
    responsemap = writecontourmap(hObject, handles);
end

zoom(responsemap, 'reset')
dfp=handles.din.cleandelf;
dgp=handles.din.cleandelg;
timeinp=handles.din.qcmt;
todo = size(handles.dout);

% grhop = handles.saveplot.grhop;
% drhop = handles.saveplot.drhop;
% phi  = handles.saveplot.phip;
grhop = handles.dout.grhop;
drhop = handles.dout.drhop;
phi = handles.dout.phip;

for n = 1:2:5
    for k = 1:length(grhop)
        lrho = lambdarhof(3, n, grhop(k,1)*1000, phi(k,1));
        dp(k,n) = drhop(k,1)/lrho/1000;
        phip(k) = phi(k,1);
    end
end
% for k=1:todo(1)
%     %for m = mnumber
%     struct2var(handles.dout{k}) %data for specific harmonic combination
%     dp(k,:)=d(1:5);
%     phip(k)=phi(1);
% end

if isempty(find(dp>0))
    warndlg('All of the data is negative and cannot be plotted')
    return
end

datalinecolor{1}=[1,0,0];
datalinecolor{3}=[0,0.5,0];
datalinecolor{5}=[0,0,1];
plot1 = responsemap.Children(2);
plot2 = responsemap.Children(4);

for i = 1:2:5
    left(i) = plot(plot1, dp(:,i), phip, 'o', 'markersize', 4);
    right(i) = plot(plot2, dp(:,i), phip, 'o', 'markersize', 4);
    set([left(i), right(i)], 'color', datalinecolor{i});
end

% set(plot1, 'ylim', [max(floor(min(phip)/10)*10,0) min(ceil(max(phip)/10)*10,90)]);
% set(plot1, 'xlim', [max(floor(min(min(dp(:,1:2:5)))/.3)*.3,0) min(ceil(max(max(dp(:,1:2:5)))/.3)*.3,1)]);

function handles = xlabeltext_Callback(hObject, eventdata, handles)
%Updates the saved values for the axes and plotting.
switch get(handles.xlabeltext,'value')
    case 1
        handles.plotting.unit = 'min';
        handles.plotting.xlabeltext = 't (min.)';
        handles.plotting.timefactor = 1;
    case 2
        handles.plotting.unit = 'hr.';
        handles.plotting.xlabeltext = 't (hr.)';
        handles.plotting.timefactor = 60;
    case 3
        handles.plotting.unit = 'day';
        handles.plotting.xlabeltext = 't (day)';
        handles.plotting.timefactor = 1440;
    case 4
        handles.plotting.unit = 'month';
        handles.plotting.xlabeltext = 't (month)';
        handles.plotting.timefactor = 43830;
    case 5
        handles.plotting.unit = 'year';
        handles.plotting.xlabeltext = 't (year)';
        handles.plotting.timefactor = 525969; % Sidereal year. I had to pick one.
end
plotraw(hObject, eventdata, handles); %Replots frequency and dissipation
% Replots calculated figure widows if present
try
    if isvalid(getappdata(0, 'outputplot'))
        close(getappdata(0, 'outputplot'));
        close(getappdata(0, 'inputplot'));
        plotvalues(hObject, handles);
    end
catch
end
guidata(hObject, handles);

function accesshandles_Callback(hObject, eventdata, handles)
keyboard;

function loadsavedplot_Callback(hObject, eventdata, handles)
handles.saveplot = load([handles.din.qcmpath handles.din.filebase '.mat']);
handles.saveplot.plot = 1;
plotvalues(hObject, handles);
handles.saveplot.plot = 0;
if isfield(handles.saveplot, 'refG') 
    set(handles.modulustoplot, 'value', ceil(handles.saveplot.refG/2))
else %Backwards compatibility from when refG was always 1
    set(handles.modulustoplot, 'value', 1)
end
guidata(hObject,handles);

function logG_Callback(hObject, eventdata, handles)
% Toggles the modulus axis type on the output plot if it is plotted.
if isvalid(getappdata(0, 'outputplot')) %If figure is currently plotted
    if get(handles.logG, 'value');
        axistype = 'log';
    else
        axistype = 'linear';
    end
    fig = getappdata(0, 'outputplot');
    fig.Children(4).YScale = axistype;
end

function loggamma_Callback(hObject, eventdata, handles)
    if get(handles.loggamma, 'value');
        axistype = 'log';
    else
        axistype = 'linear';
    end
handles.axes2.YScale = axistype;

function outputs = calculatevalues(handles, drho, grhoref, phi, refG)
drho = 0.001*drho;
grhoref = 1000*grhoref;

Zq = handles.constants.zq;
f1 = handles.constants.f1;

for n = [1,3,5,7]
    lrho(n)=lambdarhof(refG, n, grhoref, phi);
    deltarho(n)=lrho(n)*(1/(2*pi))*cotd(phi/2);
    grho(n)=grhof(refG, n, grhoref, phi);
    d(n)=drho/lrho(n);
    dfds(n)=real(delfstar(d(n),phi));
    f(n)=sauerbrey(n,drho)*dfds(n);
    g(n)=sauerbrey(n,drho)*imag(delfstar(d(n),phi));
    delfinf(n)=-(f1/(pi*Zq))*(grho(n))^0.5*sind(phi/2);
    delginf(n)=(f1/(pi*Zq))*(grho(n))^0.5*cosd(phi/2);
end

outputs.drho = drho*1000;
outputs.grho = grho*0.001;
outputs.phi = phi;
outputs.dfcalc = f;
outputs.dgcalc = g;
outputs.d = d;
outputs.lrho = lrho*1000;
outputs.deltarho = deltarho*1000;
outputs.dfds = dfds;
outputs.delfinf = delfinf;
outputs.delginf = delginf;

function [contourplots, handles] = redrawtable(hObject, handles, mnumber, outputs)
refG = itonh(get(handles.modulustoplot, 'value'));
phi = str2double(handles.propertytable.Data(1,3));
grhoref = str2double(handles.propertytable.Data(1,2));
drho = str2double(handles.propertytable.Data(1,1));

if ~exist('outputs', 'var')
    outputs = calculatevalues(handles, drho, grhoref, phi, refG);
end

handles.deldatadata(1:4,2) = outputs.dfcalc(1:2:7);
handles.deldatadata(1:4,4) = outputs.dgcalc(1:2:7);
handles.deldatadata(1:4,5) = outputs.dfds(1:2:7);
handles.deldatadata(1:4,6) = outputs.grho(1:2:7);
handles.deldatadata(1:4,7) = outputs.lrho(1:2:7);
handles.deldatadata(1:4,8) = outputs.d(1:2:7);
handles.deldatadata(1:4,9) = outputs.deltarho(1:2:7);
handles.deldatadata(1:4,10) = outputs.delfinf(1:2:7);
handles.deldatadata(1:4,11) = outputs.delginf(1:2:7);

formats={'%8.0f','%8.0f','%8.0f','%8.0f','%6.3f','%8.2e','%7.3f','%7.3f','%7.3f','%8.0f','%8.0f'};

for i = 1:4 % these are the rows
    for j=5:9; % these are the columns
        if isnan(handles.deldatadata(i,j))
            handles.deldatatable.Data(i,j) = {' -'};
        else
            handles.deldatatable.Data(i,j) = {num2str(handles.deldatadata(i,j),formats{j})};
        end
    end
    for j=[1,2,3,4,10,11]
        if isnan(handles.deldatadata(i,j)) || handles.deldatadata(i,j) == 0
            handles.deldatatable.Data(i,j) = {' -'};
        else
            handles.deldatatable.Data(i,j) = {commanumber(handles.deldatadata(i,j))};
        end
    end
end

contourplots = 0;
% Now plot contours if desired
if get(handles.plotsolutioncontours, 'value') == 1
    for q = mnumber
        contourplots = plotcontours(handles, outputs, q);
    end
end


function handles = writetableheaders(hObject,handles, size)
handles.propertytable.RowName={'value',...
    'error'};
handles.propertytable.ColumnName= {['<html><' size '>d&rho (g/m<sup>2</sup>)</' size '></html>'], ...
    ['<html><' size '>G<sub>' num2str(itonh((get(handles.modulustoplot, 'value'))))...
    '</sub>&rho<br/>(Pa-g/cm<sup>3</sup>)</' size '></html>'], ['<html><' size '>&phi (deg)</' size '></html>'],...
    ['<html><' size '>&eta&rho<br/>(Pa-s-g/cm<sup>3</sup>)</' size '></html>'],...
    ['<html><' size '>n<sub>1</sub>:n<sub>2</sub>,n<sub>3</sub></' size '?</html>']};
handles.deldatatable.RowName={'n=1','n=3','n=5','n=7'};
handles.deldatatable.ColumnName={['<html><' size '>&Delta f (Hz)</' size '></html> '],...
    ['<html><' size '>&Delta f<sub>calc</sub> (Hz)</' size '></html> '],...
    ['<html><' size '>&Delta &Gamma (Hz)</' size '></html> '],...
    ['<html><' size '>&Delta &Gamma<sub>calc</sub> (Hz)</' size '></' size '></html> '],...
    ['<html><' size '>&Delta f/&Delta f<sub>s</sub></' size '></html>'],...
    ['<html><' size '>G&rho (Pa-g/cm<sup>3</sup>)</' size '></html> )'],...
    ['<html><' size '>&lambda &rho (g/m<sup>2</sup>)</' size '></html> '],...
    ['<html><' size '>d/&lambda</' size '></html>'],...
    ['<html><' size '>&delta &rho (g/m<sup>2</sup>)</' size '></html>'],...
    ['<html><' size '>&Delta f<sub>&infin</sub> (Hz)</' size '></html>'],...
    ['<html><' size '>&Delta &Gamma<sub>&infin</sub> (Hz)</' size '></html>']};
if get(handles.issim, 'value')
    handles.offsettable.ColumnName={['<html><' size '>f<sub>bare</sub>(Hz) error</' size '></html>'],...
        ['<html><' size '>&Gamma<sub>bare</sub> (Hz) error</' size '></html>']};
else
    handles.offsettable.ColumnName={['<html><' size '>f<sub>bare</sub>(Hz)</' size '></html>'],...
        ['<html><' size '>&Gamma<sub>bare</sub> (Hz)</' size '></html>']};
end
handles.offsettable.RowName={'n=1','n=3','n=5','n=7'};
handles.contourtable.RowName={'d/lam', 'phi', 'rh', 'rd'};
handles.contourtable.ColumnName={['<html><' size '>min</' size '></html>'],...
    ['<html><' size '>max</' size '></html>']};

guidata(hObject, handles);

function resettablesize(hObject,handles)
%This function runs after resizing the figure to make sure the tables fit
%their contents
handles.propertytable.Position(3:4) = handles.propertytable.Extent(3:4);
handles.deldatatable.Position(3:4) = handles.deldatatable.Extent(3:4);
handles.offsettable.Position(3:4) = handles.offsettable.Extent(3:4);
guidata(hObject, handles);

% --- Executes when figure is resized.
function figure_SizeChangedFcn(hObject, eventdata, handles) %#ok<*INUSL>
resettablesize(hObject,handles)

function contourtable_CellEditCallback(hObject, eventdata, handles) %#ok<*DEFNU>
data = hObject.Data;
%Check that values are increasing--displays a warning if invalid
if data(1,1) >= data(1,2) || data(2,1) >= data(2,2) || data(3,1) >= data(3,2) || data(4,1) >= data(4,2)
    warndlg('Make sure that all pairs of values increase from left to right')
    return
end
%Update limits on all plots to reflect values in the table
currentfigures = get(0, 'children'); %Get all open figures
for i = 1:length(currentfigures)
    if strcmp(currentfigures(i).Name, 'contourfig') %4 is r_h, %2 is r_d
        currentfigures(i).Children(4).XLim = [data(1,1) data(1,2)];
        currentfigures(i).Children(4).YLim = [data(2,1) data(2,2)];
        currentfigures(i).Children(2).XLim = [data(1,1) data(1,2)];
        currentfigures(i).Children(2).YLim = [data(2,1) data(2,2)];
        currentfigures(i).Children(4).CLim = [data(3,1) data(3,2)];
        currentfigures(i).Children(2).CLim = [data(4,1) data(4,2)];
    end
end

function modulustoplot_Callback(hObject, eventdata, handles)
%Redraws table headers to reflect the new plotted modulus
writetableheaders(hObject,  handles, fontsize_Callback(hObject, eventdata, handles));
redrawtable(hObject, handles, 0)

function modulustoplot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function isbare_Callback(hObject, eventdata, handles)

function deldatatable_CellEditCallback(hObject, eventdata, handles)
%Updates the data name stored for deldatatable
idx = eventdata.Indices
handles.deldatadata(idx(1), idx(2)) = str2num(char(handles.deldatatable.Data(idx(1), idx(2))));
guidata(hObject, handles);

function issim_Callback(hObject, eventdata, handles)
size = fontsize_Callback(hObject, eventdata, handles);
if get(handles.issim, 'value')
    handles.offsettable.ColumnEditable = logical([1 1]); %Allow values to be input
    %Set offset table values to be the default error to give a reasonable
    %starting point (these can be changed manually)
    handles.din.bare.error.f = handles.constants.error.f;
    handles.din.bare.error.g = handles.constants.error.g;
    handles.offsettable.Data(1:4,1) = mat2cell(handles.constants.error.f(1:2:7), 1, [1 1 1 1]);
    handles.offsettable.Data(1:4,2) = mat2cell(handles.constants.error.g(1:2:7), 1, [1 1 1 1]);
    % Update header to reflect that it's just the error
    handles.offsettable.ColumnName={['<html><' size '>f<sub>bare</sub>(Hz) error</' size '></html>'],...
        ['<html><' size '>&Gamma<sub>bare</sub>(Hz)error</' size '></html>']};
    handles.deldatatable.Data(1:4,[1, 3]) = {'-'};
    cla(handles.axes2); cla(handles.axes1) %Clear plots
    set(handles.opencond, 'visible', 'off')
    set(handles.qcmfile, 'string', '')
else
    handles.offsettable.ColumnEditable = logical([0 0]); %Freeze table from input
    handles.offsettable.ColumnName={['<html><' size '>f<sub>bare</sub>(Hz)</' size '></html>'],...
    ['<html><' size '>&Gamma<sub>bare</sub>(Hz)</' size '></html>']};
    handles.deldatatable.Data(1:4,1:11) = {'-'}; %Clear all data from table
    handles.offsettable.Data(1:4,1:2) = {'-'}; %Clear all data from offset table
    set(handles.opencond, 'visible', 'on')
end
guidata(hObject, handles);

function offsettable_CellEditCallback(hObject, eventdata, handles)
%Updates stored error data to match values in the table
handles.din.bare.error.f(1:2:7) = [handles.offsettable.Data{1:4,1}];
handles.din.bare.error.g(1:2:7) = [handles.offsettable.Data{1:4,2}];
guidata(hObject, handles);

function i = nhtoi(nh)
%Converts [1 3 5] to [1 2 3]
i = ceil(nh/2);

function nh = itonh(i)
%Converts [1 2 3] to [1 3 5]
nh = i*2-1;

function outputs = blankoutput()
maxnh = 7;
outputs.df = NaN(1,maxnh);
outputs.dg = NaN(1,maxnh);
outputs.drho = NaN;
outputs.grho = NaN(1,maxnh);
outputs.phi = NaN;
outputs.grhoe = NaN;
outputs.drhoe = NaN;
outputs.phie = NaN;
outputs.dfcalc = NaN(1,maxnh);
outputs.dgcalc = NaN(1,maxnh);
outputs.d = NaN(1,maxnh);
outputs.lrho = NaN(1,maxnh);
outputs.deltarho = NaN(1,maxnh);
outputs.dfds = NaN(1,maxnh);
outputs.delfinf = NaN(1,maxnh);
outputs.delginf = NaN(1,maxnh);

function lrho = lambdarhof(refn, n, grho, phi)
if refn == n
    lrho = 1/(n*5e6)*(grho^0.5)/cosd(phi/2);
else
    lrho = 1/(n*5e6)*((grho*(n^(phi/90))/(refn^(phi/90)))^0.5)/cosd(phi/2);
end

function grho = grhof(refn, n, grhoref, phi)
grho = grhoref*(n^(phi/90))/(refn^(phi/90));

function dl = dlf(refn, refnh, dlin, phi);
dl = dlin*(refnh/refn)^(1-(phi/180));

function thickness_Callback(hObject, eventdata, handles)
%This function opens a pop up window which you can use to calculate the
%maximum (optimum) thickness of a film. The boxes auto-fill with the
%properties currently showing on the screen but can be changed.
close(figure(996)); figure(996); %close window if open
pos=get(figure(996),'position');
set(figure(996),'position',[pos(1)/1 pos(2)/1 pos(3)/1.7 pos(4)/1.5]); %Set window size
set(figure(996),'menubar','none','toolbar','none','numbertitle','off','name','Optimum thickness calculator');
laxis = axes('parent',figure(996),'units','normalized','position',[0 0 1 1],'visible','off');
%Create label and input box for the maximum dissipation
c1ledge = 0.15;
c2ledge = 0.6;
c1w = .4;
max_dg = uicontrol('style','edit','string','10','units','normalized',...
    'position',[c2ledge 0.85 .25 .1]);
% uicontrol('style','text','string','max DG (kHz)','units','normalized',...
%     'position',[c1ledge 0.82 c1w .1]);
text(c1ledge, 0.9, 'max \Delta\Gamma (kHz)');

%Create label and input box for the maximum thickness
max_drho = uicontrol('style','edit','string','10','units','normalized',...
    'position',[c2ledge 0.75 .25 .1]);
% uicontrol('style','text','string','max dp (g/m^2)','units','normalized',...
%     'position',[c1ledge 0.72 c1w .1]);
text(c1ledge, 0.8, 'max d\rho (g/m^2)');

%Create label and input box for the phase angle It defaults to the current
%phase angle in onelayergui
phi = uicontrol('style','edit','string',handles.propertytable.Data{1,3},'units','normalized',...
    'position',[c2ledge 0.65 .25 .1]);
% uicontrol('style','text','string','Phi','units','normalized',...
%     'position',[c1ledge 0.62 c1w .1]);
text(c1ledge, 0.7, '\phi');

%Create label and pull-down for the modulus reference harmonic It defaults
%to the value in onelayergui
ref_G = uicontrol('style','popup','string',{'<html>G<sub>1</sub></html>',...
    '<html>G<sub>3</sub></html>', '<html>G<sub>5</sub></html>',...
    '<html>G<sub>7</sub></html>'},'value', get(handles.modulustoplot, 'value'),...
    'units','normalized', 'position',[c2ledge 0.55 .25 .1]);
% uicontrol('style','text','string','Reference G','units','normalized',...
%     'position',[c1ledge 0.52 c1w .1]);
text(c1ledge, 0.6, 'reference G');

%Create label and input box for the modulus It defaults to the current
%modulus in onelayergui
grho = uicontrol('style','edit','string',handles.propertytable.Data{1,2},'units','normalized',...
    'position',[c2ledge 0.45 .25 .1]);
% uicontrol('style','text','string','|G*|p (Pa-g/cm^3)','units','normalized',...
%     'position',[c1ledge 0.42 c1w .1]);
text(c1ledge, 0.5, '|G*|p (Pa-g/cm^3)');

%Create check boxes for the different harmonics
nh1 = uicontrol('style','checkbox','string','5MHz','units','normalized',...
    'position',[.1 0.28 .25 .1]);
nh3 = uicontrol('style','checkbox','string','15MHz','units','normalized',...
    'value', 1, 'position',[.1+.19 0.28 .25 .1]);
nh5 = uicontrol('style','checkbox','string','25MHz','units','normalized',...
    'value', 1, 'position',[.1+.42 0.28 .25 .1]);
nh7 = uicontrol('style','checkbox','string','35MHz','units','normalized',...
    'position',[.1+.65 0.28 .25 .1]);

% calcmax = uicontrol('style','text','string','Max drho: <html>&Gamma;</html>','units','normalized',...
%     'position',[0.4 0.20 .3 .1]); %Displays the maximum thickness 
calcmax = text(0.50, 0.15, 'Max d\rho: ');

%Initializes button and callback for the Calculate button
uicontrol('style','push','string','Calculate','units','normalized',...
    'position',[c1ledge 0.10 .3 .1],'callback',{@maxthicknesscalc_Callback,...
    max_dg, max_drho, phi, ref_G, grho, nh1, nh3, nh5, nh7, calcmax});

set(findall(figure(996),'-property','FontSize'),'FontSize',12)

function maxthicknesscalc_Callback(hObject, ~, max_dg, max_drho, phi, ref_G, grho, nh1, nh3, nh5, nh7, calcmax)
%This function calculates the maximum thickness that a film can be before
%the dissipation of a particular harmonic reaches a critical
%threshhold--usually 10,000Hz, but this can be changed. This function does
%not take into account decay length limitations on the thickness for mass
%sensitivity and does not explore possibilities beyond the dissipation
%spike at d/lambda = 0.25.

% Get which harmonics the user is concerned about
nh = [1 3 5 7];
nh = nh(logical([get(nh1, 'value') get(nh3, 'value') get(nh5, 'value') get(nh7, 'value')]));
if isempty(nh)
    warndlg('You forgot to select harmonics to look at')
    return
end
grhoi = str2num(get(grho, 'string'))*1000; %Get modulus in the right units
phii = str2num(get(phi, 'string')); %Get phi
refn = itonh(get(ref_G, 'value')); %Get refn
maxdrho = str2num(get(max_drho, 'string'))*0.001; %Get maxdrho
maxdg = str2num(get(max_dg, 'string'))*1000; %Get max dissipation

% Returns if any values are either empty or otherwise not numbers
if isempty(grhoi) || isempty(phii)
    warndlg('Please give a valid value for the physical parameters')
    return
end

for n = nh
    lambdarhon = lambdarhof(refn, n, grhoi, phii); %Lambda rho for the harmonic of interest
    drho = 0.001*0.1; % Initialize for while loop starting at .1g/m^2
    delfstarc = 0;
    % Loops over increasing thicknesses to find when the dissipation
    % becomes to great
    while imag(delfstarc) < maxdg && drho < maxdrho
        drho = drho + 0.000025; % This sets the increment -- resolution is to the 0.025g/m^2
        dl = drho./lambdarhon; %Updates dl based on the new drho
        delfstarc = sauerbrey(n,drho)*delfstar(dl, phii); %Calculates the frequency and dissipation     
    end
    d(n) = (drho-0.000025)*1000;
end
set(calcmax, 'string', ['Max d\rho: ' num2str(min(d(d~=0))) 'g/m^2'])  %Displays the maximum thickness

function plotsolutioncontours_Callback(hObject, eventdata, handles)
%Runs when you toggle the plotcountours button. Opens window to edit the
%limits in the figures when it is turned on, and closes the window when
%turned off.
if get(handles.plotsolutioncontours, 'value') == 1
    close(figure(995)); figure(995); %close window if open
    pos=get(figure(995),'position');
    set(figure(995),'position',[pos(1)/2.2 pos(2)/2.2 pos(3)/2.7 pos(4)/2.4]); %Set window size
    set(figure(995),'menubar','none','toolbar','none','numbertitle','off','name','Contour limits');
    rname = {'d/l', 'phi', 'rh', 'rd'}; %Unable to get row html to work with row names
    % rname = {'<html><h3>d/&lambda;</h3></html>', '<html><h3>&phi;</h3></html>',...
    %     '<html><h3>r<sub>h</sub></h3></html>', '<html><h3>r<sub>d</sub></h3></html>'};
    cname = {['<html><h2>min</h2></html>'], ['<html><h2>max</h2></html>']};
    t = uitable(figure(995), 'Data', handles.contourtable.Data, ...
        'ColumnName', cname, 'RowName', rname, 'ColumnEditable', true, ...
        'CellEditCallback', {@contourtable_CellEditCallback, handles});
    %t.ColumnWidth{1} = 10;
    set(findall(figure(995),'-property','FontSize'),'FontSize',12)
    t.Position(3) = t.Extent(3);
    t.Position(4) = t.Extent(4);
elseif get(handles.plotsolutioncontours, 'value') == 0
    close(figure(995));
end
