function varargout = temperaturecorrect(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @temperaturecorrect_OpeningFcn, ...
    'gui_OutputFcn',  @temperaturecorrect_OutputFcn, ...
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

% --- Executes just before temperaturecorrect is made visible.
function temperaturecorrect_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.mdate.String='m: 11/3//15';
handles.figdate.String='fig: 9/1/15';
warning('off','MATLAB:Axes:NegativeDataInLogAxis') %Supress warning about negative data
set(hObject,'toolbar','auto'); %Remove default toolbar
% Update handles structure
linkaxes([handles.axes1, handles.axes2], 'off') % This seems a leftover from before
%Sets constants for use in the program.
handles.constants.f1 = 5e6;
handles.constants.zq = 8.84e6;
handles.constants.error.f = [44 0 120 0 273 0 500]; %Check 7harm
handles.constants.error.g = [11 0 22 0 14 0 10]; %Check 7harm

handles.constants.nhvals = {[1,3,1] [1,3,3] [1,5,1] [1,5,5] [3,5,3] [3,5,5]...
    [1,7,1] [1,7,7] [3,7,3] [3,7,7] [5,7,5] [5,7,7] [1,3,5] [1,5,3] [3,5,1]};
handles.constants.label = {'1:3,1' '1:3,3' '1:5,1' '1:5,5' '3:5,3' '3:5,5',...
    '1:7,1', '1:7,7', '3:7,3', '3:7,7', '5:7,5', '5:7,7', '1:3,5', '1:5,3', '3:5,1'};
handles.saveplot.plot = 0; %default to not
resettablesize(hObject,handles)
writetableheaders(hObject, handles, 'h3'); %Sets default size on opening to html h3
set(gcf,'Pointer','arrow'); %Sometimes the cursor gets stuck on (I think that was the reason for this)

% --- Outputs from this function are returned to the command line.
function varargout = temperaturecorrect_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function changeqcmfile_Callback(hObject, eventdata, handles, qcmfile, qcmpath)
% This function imports data from a file, stores it (including raw f and g
% data, spectra, and error ranges if available) in handles, and displays plots.

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

handles.restorefile = [qcmpath qcmfile];

% If someone cancels out of the uigetfile dialog, qcmfile and qcmpath will
% both be 0. This checks if that is the case.
if ~qcmfile
    set(handles.statusbar, 'string', 'You did not select a new file. No data has been changed', 'BackgroundColor', 'yellow')    
    return
end

set(handles.qcmfile,'string', ['Loading ' qcmfile])
drawnow; %means that the above change will be visible.

%clears the data currently stored in handles.din, as that is all dependent
%on the data file.
if isfield(handles, 'din')
    rmfield(handles, 'din');
end

if isfield(handles, 'saveplot') %Indicates that the program shouldn't plot saved data
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
        for nh = 1:2:max(length(error.f))
            fref=['f',num2str(nh)];
            gref=['g',num2str(nh)];
            handles.offsettable.Data(nhtoi(nh),1)={[commanumber(offset.(fref)) setstr(177) num2str(error.f(nh),'%8.0f')]};
            handles.offsettable.Data(nhtoi(nh),2)={[commanumber(offset.(gref)) setstr(177) num2str(error.g(nh),'%8.0f')]};
        end
        handles.din.bare.error.f = error.f;
        handles.din.bare.error.g = error.g;
    catch Err
        set(handles.statusbar, 'string', 'Doesn''t contain reference information. Default error values will be used', 'BackgroundColor', 'yellow')
        for n = 2:9
            offset.(name{n}) = 0;
        end
        handles.offsettable.Data(1:4,1:2) = {'-'};
        handles.din.bare.error.f = handles.constants.error.f; 
        handles.din.bare.error.g = handles.constants.error.g; 
    end
    
    handles.overwrite = 0; %Disables saving of a new _data file later on.
else 
    set(handles.statusbar, 'string', 'Please select a _data file', 'BackgroundColor', 'red')
    return
end

try
    m = matfile([handles.din.qcmpath handles.din.filebase(1:end-5) '_tempshift_data.mat'],'Writable',true);
    if isfield(m.temp, 'table') %Loads previously stored temperature data
        temp = m.temp;
        handles.datatemptable.Data = temp.table;
    end
catch Err
    Err
end
    
%Saves all of the processed data into the handles structure.
handles.din.qcmt = time;
handles.din.delf = delf;
handles.din.delg = delg;
handles.din.absf = absf;
handles.din.absg = absg;
handles.din.qcmpath = qcmpath;
handles.din.qcmfile = qcmfile;
handles.din.filebase = filebase;
handles.din.bare.offset = offset;
handles.din.bare.error.f = error.f; 
handles.din.bare.error.g = error.g;

clear qcmt delf delg absf absg temp
[qcmreffile, qcmrefpath] = uigetfile({'*.xls;*.xlsx;*.fre;*.mat'}, 'Select the temperature reference data file', handles.din.qcmpath);
load([qcmrefpath, qcmreffile]);
handles.dinref.qcmfile = [qcmrefpath, qcmreffile];
handles.dinref.qcmt = time;
handles.dinref.delf = delf;
handles.dinref.delg = delg;
handles.dinref.absf = absf;
handles.dinref.absg = absg;
if exist('corrections', 'var') && isfield(corrections, 'table') % Loads previously stored temp data
    handles.referencetemptable.Data = corrections.table;
    handles.reftemp = corrections.temp;
    handles.refshifts = corrections.refshifts;
end

%Plots the imported data to the lower plots
plotraw(hObject, eventdata, handles);

%Imports the conductance data, if any

%Changes the data file to not be "loading" anymore.
set(handles.qcmfile,'string',qcmfile);
set(handles.statusbar, 'string', 'File successfully loaded', 'BackgroundColor', 'green')
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
end

function Solve_Callback(hObject,eventdata,handles)
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')
% Uses the top selected harmonic combination to solve the currently loaded
% point.
mnumber = getnhvals(handles);
[~, handles] = findsolution(hObject,eventdata,handles,1,mnumber(1));
guidata(hObject, handles); %Saves new handles.

function size = fontsize_Callback(hObject, eventdata, handles)  % change font size of gui
% Changes size of headers and fonts throughout
switch get(handles.fontsize,'Value')
    case 1
        fontsize=10;
        size = 'h5';
    case 2
        fontsize=11;
        size = 'h5';
    case 3
        fontsize=12;
        size = 'h4';
    case 4
        fontsize=13;
        size = 'h4';
    case 5
        fontsize=14;
        size = 'h3';
    case 6
        fontsize=15;
        size = 'h3';
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
hgsave('onelayergui3.fig')

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

datat = handles.din.qcmt;
reft = handles.dinref.qcmt;
dataf = handles.din.delf;
refdelf = handles.dinref.delf;

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
[length harms] = size(refdelf);
for nh=1:2:harms
    plot(handles.axes1,reft./timefactor,refdelf(:,nh)/nh,symbols{nh},'color',colors{nh})
    plot(handles.axes2,datat./timefactor,dataf(:,nh)/nh,symbols{nh}, 'color',colors{nh})
    legendtext=[legendtext,legends{nh}];
    hold(handles.axes1,'on')
    hold(handles.axes2,'on')
end

xlabel(handles.axes1, xlabeltext)
ylabel(handles.axes1,'{\Delta}f/n reference (Hz)')
legend(handles.axes1, legendtext,'location','best')

xlabel(handles.axes2, xlabeltext)
ylabel(handles.axes2,'{\Delta}f data (Hz)')
legend(handles.axes2, legendtext,'location','best')

set(handles.axes1,'xlimmode','auto')
set(handles.axes1,'ylimmode','auto')
set(handles.axes2,'xlimmode','auto', 'yscale', 'lin')
set(handles.axes2,'ylimmode','auto')

if get(handles.log,'value') %Sets time axis scale
    set([handles.axes1 handles.axes2],'xscale','log')
else
    set([handles.axes1 handles.axes2],'xscale','linear')
end
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
outputstructure = reformatdata(handles, outputs);
handles.dout = outputstructure;
if exist('contourplots', 'var') %Should be created during findsolution if desired
    handles.contourplots = contourplots;
end
guidata(hObject,handles);
set(handles.statusbar, 'string', '', 'BackgroundColor', 'default')

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
selection = questdlg('Do you want to save the file?',...
    'Close Request Function',...
    'Yes','No','Cancel','Yes');
switch selection,
    case 'Yes',
        handles=guidata(hObject);
        if isfield(handles,'din')
            hgsave(hObject,'onelayergui3.fig')
            delete(hObject)
        else
            fprintf(1,'handles.din does not exist')
            return
        end
    case 'No'
        delete(hObject)
    case 'Cancel'
        return
end

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

function outstruct = reformatdata(handles, outputs)
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
if handles.overwrite == 1 %Disabled for _data files.
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
end

function cleaned_Callback(hObject, eventdata, handles)
% The value of this function is called to determine things, but when it is
% changed the only thing necessary is to replot the data. Switching it to
% the cleaned data is taken care of by using "get" in the plotraw function.
plotraw(hObject, eventdata, handles)

function savecleaned_Callback(hObject, eventdata, handles)
% Changes variables in the _data file, which contains the "clean" values.
m = matfile([handles.din.qcmpath handles.din.filebase(1:end-5) '_tempshift_data.mat'],'Writable',true);
m.delf = handles.din.delfshift;
m.delg = handles.din.delg;
m.absf = handles.din.absfshift;
m.absg = handles.din.absg;
m.time = handles.din.qcmt;
m.offset = handles.din.bare.offset;
m.error = handles.din.bare.error;
temp.table = handles.datatemptable.Data;
temp.reftemp = handles.reftemp;
temp.refshifts = handles.refshifts;
m.temp = temp;

mref = matfile(handles.dinref.qcmfile, 'Writable', true);
corrections.table = handles.referencetemptable.Data; %saves data from the table
corrections.temp = handles.reftemp; %Saves calculated temps
corrections.refshifts = handles.refshifts; %Saves calculated shift values
mref.corrections = corrections;
% mref.

function opencond_Callback(hObject, eventdata, handles)
% This function opens the second program, condfig, which displays
% conductance data.
if isempty(handles.cond)
    set(handles.statusbar, 'string', 'You cannot view the spectra since there aren''t any', 'BackgroundColor', 'red')
    return
end
set(handles.statusbar, 'string', 'Loading conductance data...', 'BackgroundColor', 'default')
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
handles.din = condfig5('onelayerguiwithcond', handles.figure);
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
if handles.plotting.timefactor == 1
    output_txt = {['Time: ', num2str(pos(1),'%f.2'), ' min'],...
    ['y: ',num2str(pos(2),'%f.1')]};
else
    output_txt = {['Time: ', num2str(pos(1),5), ' ', handles.plotting.unit],...
    ['Time: ', num2str(pos(1)*handles.plotting.timefactor,'%f.2'), ' min'],...
    ['y: ',num2str(pos(2),'%f.1')]};
end

% The following determines the time requested and sends it to the
% condfig window if that is the only other window open. Otherwise I'm not
% sure how to find the handle to the right window.

% Finds handles to all the windows
getting = get(0, 'children');

%If there are only two figure windows open, the second one is the condfig
%window. This then writes the time from above to the time box in condfig
% and runs the "view_Callback" function.
if max(size(getting)) == 2
    condfighandle = getting(2);
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
end
datacursormode off


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

function handles = writetableheaders(hObject,handles, size)
    handles.offsettable.ColumnName={['<html><' size '>f<sub>bare</sub>(Hz)</' size '></html>'],...
        ['<html><' size '>&Gamma<sub>bare</sub> (Hz)</' size '></html>']};
handles.offsettable.RowName={'n=1','n=3','n=5','n=7'};
guidata(hObject, handles);

function resettablesize(hObject,handles)
handles.offsettable.Position(3:4) = handles.offsettable.Extent(3:4);
guidata(hObject, handles);

% --- Executes when figure is resized.
function figure_SizeChangedFcn(hObject, eventdata, handles) %#ok<*INUSL>
resettablesize(hObject,handles)

function i = nhtoi(nh)
%Converts [1 3 5] to [1 2 3]
i = ceil(nh/2);

function nh = itonh(i)
%Converts [1 2 3] to [1 3 5]
nh = i*2-1;

function handles = ignoretimebox_Callback(hObject, eventdata, handles)
handles.ignoretime = str2num(get(handles.ignoretimebox, 'string'));
guidata(hObject, handles);

function handles = avgtimebox_Callback(hObject, eventdata, handles)
handles.avgtime = str2num(get(handles.avgtimebox, 'string'));
guidata(hObject, handles);

function calcref_Callback(hObject, eventdata, handles)
data = handles.referencetemptable.Data;
time = handles.dinref.qcmt;
delf = handles.dinref.delf;

%If end times are not given, take the start time of the next temperature as the end time
if isempty(data{1,3}) 
    for i = 1:length(data)-1
        data{i,3} = data{i+1,2};
    end
    data{length(data),3} = max(time); %final end time is the maximum time measured
end

for i = 1:length(data)
    tempi = data{i,1};
    starti = data{i,2};
    endi = data{i,3};
    if ~isa(tempi, 'double')
        tempi = str2num(tempi);
    end
    if ~isa(starti, 'double')
        tempi = str2num(starti);
    end
    if ~isa(endi, 'double')
        tempi = str2num(endi);
    end
    
    if isnan(tempi)
        data = data(1:i-1,:);
        if ~isa(data{i-1,3}, 'double')
            data{i-1,3} = str2num(data{i-1,3});
        end
        if isempty(data{i-1,3}) || isnan(data{i-1,3})
           data{i-1,3} = max(time); 
        end
        break
    end
    temp(i) = tempi;
    try
        start = starti+handles.ignoretime;
        begavg = endi-handles.avgtime;
    catch
        handles = ignoretimebox_Callback(hObject, eventdata, handles);
        handles = avgtimebox_Callback(hObject, eventdata, handles);
        start = starti+handles.ignoretime;
        begavg = endi-handles.avgtime;
    end
        
    last = endi;
    if start > begavg
        begavg = start;
    end
    higher = time > begavg;
    lower = time < last;
    refshifts{i} = nanmean(delf(logical(higher.*lower),:));    
end

[utemp, idx] = unique(temp); % Checks for repeated temperatures
% reassigns all calculated shifts to sorted temperatures. This step
% also removes the second appearance of a single temperature--this is
% accounted for later.
sortrefshifts(1:length(idx)) = refshifts(idx); %Sorts reference values by temperature from low to high

if length(utemp) ~= length(temp) %Checks if some temperatures appear twice   
    n = histc(temp,utemp); % counts how many times each temperature appears
    sametemps = utemp(n>1); % determines which temperatures appear more than once   
    % For each temperature which appears more than once, averages the two
    % appearances with each other.
    for i = 1:length(sametemps)
      idxsametemp{i} = find(temp==sametemps(i))
      k = 1; toavg = [];
      for j = idxsametemp{i}
          toavg(k,:) = refshifts{j};
          k = k+1;
      end
      sortrefshifts(find(utemp == sametemps(i))) = {mean(toavg)};
    end
end

% Plot the calculated values on the plot so the user can check that these
% values are correct.
plotraw(hObject, eventdata, handles) %Clear any previously plotted points
for j = 1:length(data)
    if ~isa(data{j,3}, 'double')
        data{j,3} = str2num(data{j,3});
    end
    if ~isa(data{j,1}, 'double')
        data{j,1} = str2num(data{j,1});
    end
    endtime(j) = max(find(time < data{j,3}));
    tempidx(j) = find(data{j,1} == utemp);
    for nh = 1:2:5
        plot(handles.axes1, time(endtime(j)), sortrefshifts{tempidx(j)}(nh)/nh, 'cx');
    end
end

handles.reftemp = utemp;
handles.refshifts = sortrefshifts;
set(handles.statusbar, 'string', 'Reference values calculated successfully', 'BackgroundColor', 'green')
guidata(hObject, handles);

function calcdata_Callback(hObject, eventdata, handles)
data = handles.datatemptable.Data;
time = handles.din.qcmt;
absfshift = handles.din.absf;


if isempty(data{1,3}) 
    for i = 1:length(data)-1
        data{i,3} = data{i+1,2};
    end
    data{length(data),3} = max(time); %final end time is the maximum time measured
end

for i = 1:length(data)
    if ~isa(data{i,1}, 'double')
        data{i,1} = str2num(data{i,1});
    end
    if ~isa(data{i,2}, 'double')
        data{i,2} = str2num(data{i,2});
    end
    if ~isa(data{i,3}, 'double')
        data{i,3} = str2num(data{i,3});
    end
    if isnan(data{i,1}) %Catches when the data in the table ends
        data = data(1:i-1,:);
        if isempty(data{i-1,3}) || isnan(data{i-1,3})
           data{i-1,3} = max(time); 
        end
        break
    end
    temp = data{i,1};
    start = data{i,2};
    try
        endignore = start+handles.ignoretime;
    catch
        handles = ignoretimebox_Callback(hObject, eventdata, handles);
        handles = avgtimebox_Callback(hObject, eventdata, handles);
        endignore = start+handles.ignoretime;
    end
    
    last = data{i,3};
    
    startignore = time > start;
    stopignore = time < endignore;
    higher = time > endignore;
    lower = time < last;
    
    ignoretimes = logical(startignore.*stopignore);
    absfshift(ignoretimes,1:2:5) = NaN; %remove all points within ignore time of a temperature change
    
    shifttimes = logical(higher.*lower);
    tempidx = find(temp == handles.reftemp); %Get the index of the temperature from the reference
    if isempty(tempidx)
        warndlg(['No reference data was found for ' num2str(temp) 'C. Data will be ignored.'])
        absfshift(shifttimes,1:2:5) = NaN;
    else
        absfshift(shifttimes,:) = handles.din.absf(shifttimes,:)-repmat(handles.refshifts{tempidx},[sum(shifttimes),1]);
    end
end
offsets = {'f1', '', 'f3', '', 'f5'};
for i = 1:2:5
    delfshift(:,i) = absfshift(:,i)-handles.din.bare.offset.(offsets{i});
end
handles.din.absfshift = absfshift;
handles.din.delfshift = delfshift;
set(handles.statusbar, 'string', 'Values shifted successfully', 'BackgroundColor', 'green')
handles.din.delf = delfshift;
plotraw(hObject, eventdata, handles)
guidata(hObject, handles);


function [attemp, temp] = getattemp(data, time, handles)
