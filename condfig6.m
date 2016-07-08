function varargout = condfig6(varargin)
%CONDFIG M-file for condfig6.fig
%      CONDFIG6, by itself, creates a new CONDFIG6 or raises the existing
%      singleton*.
%
%      H = CONDFIG6 returns the handle to a new CONDFIG6 or the handle to
%      the existing singleton*.
%
%      CONDFIG6('Property','Value',...) creates a new CONDFIG6 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CONDFIG_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CONDFIG6('CALLBACK') and CONDFIG6('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CONDFIG.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help condfig

% Last Modified by GUIDE v2.5 13-May-2016 12:48:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @condfig6_OpeningFcn, ...
                   'gui_OutputFcn',  @condfig6_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% --- Executes just before condfig is made visible.
function condfig6_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Is the changeme_main gui's handle is passed in varargin?
% if the name 'changeme_main' is found, and the next argument
% varargin{mainGuiInput+1} is a handle, assume we can open it.

dontOpen = false;
hObject.Name = 'condfig';
mainGuiInput = find(strcmp(varargin, 'onelayerguiwithcond'));
if (isempty(mainGuiInput)) ...
    || (length(varargin) <= mainGuiInput) ...
    || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.onelayerguiMain = varargin{mainGuiInput+1};
    
    % Obtain handles using GUIDATA with the caller's handle 
    mainHandles = guidata(handles.onelayerguiMain);
      
    % Position to be relative to parent:
    parentPosition = getpixelposition(handles.onelayerguiMain);
    currentPosition = get(hObject, 'Position');  
    % Set x to be directly in the middle, and y so that their tops align.
    newX = parentPosition(1) + (parentPosition(3)/2 - currentPosition(3)/2);
    %newY = parentPosition(2) + (parentPosition(4)/2 - currentPosition(4)/2);
    newY = parentPosition(2) + (parentPosition(4) - currentPosition(4));%*.95);
    newW = currentPosition(3);
    newH = currentPosition(4);
    
    set(hObject, 'Position', [newX, newY, newW, newH]);
end
handles.currenttimeidx = 1;
handles.main = mainHandles;
handles.direction = 1;
handles.refit = 0;
handles.refittypetoggle = 1;

% There should be a saved _cond file, so open it now.
handles.cond = load([handles.main.din.qcmpath handles.main.din.qcmfile(1:end-4) '_cond.mat']);

% Update handles structure
guidata(hObject, handles);

if dontOpen
    disp('-----------------------------------------------------');
    disp('Improper input arguments. Pass a property value pair')
    disp('whose name is "onelayerguiwithcond" and value is the handle')
    disp('to the onelayerguiwithcond figure, e.g:');
    disp('   x = onelayerguiwithcond()');
    disp('   condfig6(''onelayerguiwithcond'', x)');
    disp('-----------------------------------------------------');
else
    set(handles.currentpoint, 'string', '0')
       
    %plots the data in the main window
    plotraw(hObject, eventdata, handles)
    plotcond(hObject, handles)
    guidata(hObject,handles)
    %Not sure what is waiting, but this ensures that this window stays open
    uiwait(hObject);
end

% --- Outputs from this function are returned to the command line.
function varargout = condfig6_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This passes the changes to the cleaned up version of the data back to the
%main function.
varargout{1} = handles.main.din;
delete(hObject);

function condfig6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function condfig6_CloseRequestFcn(hObject, eventdata, handles)
if handles.refit
    time = handles.cond.time;
    spectra = handles.cond.spectra;
    fits = handles.cond.fits;
    save([handles.main.din.qcmpath handles.main.din.filebase '_cond.mat'], 'time', 'spectra', 'fits');
    disp('The spectra changes have been saved!')
end
uiresume(handles.figure);

function handles = removepoint(hObject, eventdata, handles, harmonic)
time = str2num(get(handles.currentpoint, 'string'));
if isfield(handles, 'currspectime')
    if time - handles.currspectime(1) <= 0.0051
        handles.main.din.cleandelf(handles.currenttimeidx, harmonic) = NaN;
        handles.main.din.cleandelg(handles.currenttimeidx, harmonic) = NaN;
        handles.main.din.absf(handles.currenttimeidx, harmonic) = NaN;
        handles.main.din.absg(handles.currenttimeidx, harmonic) = NaN;
    else
        warning('You haven''t updated the spectra')
        return 
    end
else
    warning('You haven''t updated the spectra')
    return    
end
set(handles.(['remove' num2str(harmonic)]), 'visible', 'off');
% guidata(hObject,handles)
plotraw(hObject, eventdata, handles)
plotcond(hObject, handles)

function remove1_Callback(hObject, eventdata, handles)
handles = removepoint(hObject, eventdata, handles, 1);
guidata(hObject,handles)

function remove3_Callback(hObject, eventdata, handles)
handles = removepoint(hObject, eventdata, handles, 3);
guidata(hObject,handles)

function remove5_Callback(hObject, eventdata, handles)
handles = removepoint(hObject, eventdata, handles, 5);
guidata(hObject,handles)

function remove7_Callback(hObject, eventdata, handles)
handles = removepoint(hObject, eventdata, handles, 7);
guidata(hObject,handles)

function plotcond(hObject, handles)
% Plots the conductance points onto the main gui.
cond = handles.main.cond.points;

%Queries current x and y limits so these can be left the same (if looking zoomed, for instance).
%ylimits1 = ylim(handles.main.axes1);
%ylimits2 = ylim(handles.main.axes2);
xlimits1 = xlim(handles.main.axes1);
xlimits2 = xlim(handles.main.axes2);

delf = handles.main.din.cleandelf;
delg = handles.main.din.cleandelg;
try
qcmt = handles.main.din.qcmt/handles.main.plotting.timefactor;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        switch get(handles.main.xlabeltext,'value')
            case 1
                handles.main.plotting.unit = 'min';
                handles.main.plotting.xlabeltext = 't (min.)';
                handles.main.plotting.timefactor = 1;
            case 2
                handles.main.plotting.unit = 'hr';
                handles.main.plotting.xlabeltext = 't (hr.)';
                handles.main.plotting.timefactor = 60;
            case 3
                handles.main.plotting.unit = 'day';
                handles.main.plotting.xlabeltext = 't (day)';
                handles.main.plotting.timefactor = 1440;
            case 4
                handles.main.plotting.unit = 'month';
                handles.main.plotting.xlabeltext = 't (month)';
                handles.main.plotting.timefactor = 43830; % month of 30.4 days.
            case 5
                handles.main.plotting.unit = 'year';
                handles.main.plotting.xlabeltext = 't (year)';
                handles.main.plotting.timefactor = 525969; % Sidereal year. I had to pick one.
        end
        xlabeltext = handles.main.plotting.xlabeltext;
        timefactor = handles.main.plotting.timefactor;
        qcmt = handles.main.din.qcmt/handles.main.plotting.timefactor;
    else
        rethrow(Err)
    end

end
%Define linestyle for the three different harmonics.
linestyle = [{'o-r'},{'o-b'},{'o-g'}];
%Goes through the three harmonics and plots them if desired. I couldn't
%figure out a way to do this such that the number was the index.

%Determines which spectra are present and plots them.
haveconductance = ~cellfun(@isempty, handles.cond.spectra); %logical 1 if spectra present
%Make this try/catch more robust--better check of how many harmonics there
%are
idx = handles.currenttimeidx;
if max(handles.main.activenh.on) == 3
    haveconductance = haveconductance(:,1:3); %eliminates higher harmonics even if present
    maxharm = 3;
    timearray = [qcmt(idx) qcmt(idx)];
elseif max(handles.main.activenh.on) == 5
    haveconductance = haveconductance(:,1:5); %eliminates higher harmonics even if present
    maxharm = 5;
    timearray = [qcmt(idx) qcmt(idx) qcmt(idx)];
else
    haveconductance = haveconductance(:,1:7); %eliminates higher harmonics even if present
    maxharm = 7;
    timearray = [qcmt(idx) qcmt(idx) qcmt(idx) qcmt(idx)];
end
qcmts = [qcmt' qcmt' qcmt' qcmt' qcmt' qcmt' qcmt'];
if maxharm == 7
    harms = repmat([1,0,3,0,5,0,7], length(qcmt), 1);
else
    harms = repmat([1,0,3,0,5], length(qcmt), 1);
end
plot(handles.main.axes1, qcmts(haveconductance(:,1:maxharm)), delf(haveconductance(:,1:maxharm))./harms(haveconductance(:,1:maxharm)), '.k');
plot(handles.main.axes2, qcmts(haveconductance(:,1:maxharm)), delg(haveconductance(:,1:maxharm)), '.k');

plot(handles.main.axes1, timearray, delf(idx, 1:2:maxharm)./[1:2:maxharm], '.c')
plot(handles.main.axes2, timearray, delg(idx, 1:2:maxharm), '.c')
%Sets ylim back to the original so very far off points aren't shown.
%ylim(handles.main.axes1,ylimits1);
%ylim(handles.main.axes2,ylimits2);
xlim(handles.main.axes1,xlimits1);
xlim(handles.main.axes2,xlimits2);
if get(handles.main.loggamma,'value')
    set(handles.main.axes2,'yscale','log')
else
    set(handles.main.axes2,'yscale','linear')
end

function plotraw(hObject, eventdata, handles)
% Plots the raw data. This function should be identical to the plotraw
% function in onelayergui3. Except that I want it to keep whatever the
% original zoom level was.

%Brings data into the function
try
    qcmt = handles.main.din.qcmt;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        disp('Unable to replot data. Please reload file')
        return
    end
end
% Plots filtered or unfiltered data depending on radio button
if get(handles.main.cleaned, 'value')
    delf = handles.main.din.cleandelf;
    delg = handles.main.din.cleandelg;
else
    delf = handles.main.din.delf;
    delg = handles.main.din.delg;
end

% now we plot the data
% Set defaults, symbols, colors, etc.
xlimits1 = xlim(handles.main.axes1);
xlimits2 = xlim(handles.main.axes2);
cla(handles.main.axes2)
cla(handles.main.axes1)
set(0,'defaultaxesfontsize',11)
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

try
    xlabeltext = handles.main.plotting.xlabeltext;
    timefactor = handles.main.plotting.timefactor;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        switch get(handles.main.xlabeltext,'value')
            case 1
                handles.main.plotting.unit = 'min';
                handles.main.plotting.xlabeltext = 't (min.)';
                handles.main.plotting.timefactor = 1;
            case 2
                handles.main.plotting.unit = 'hr';
                handles.main.plotting.xlabeltext = 't (hr.)';
                handles.main.plotting.timefactor = 60;
            case 3
                handles.main.plotting.unit = 'day';
                handles.main.plotting.xlabeltext = 't (day)';
                handles.main.plotting.timefactor = 1440;
            case 4
                handles.main.plotting.unit = 'month';
                handles.main.plotting.xlabeltext = 't (month)';
                handles.main.plotting.timefactor = 43830; % month of 30.4 days.
            case 5
                handles.main.plotting.unit = 'year';
                handles.main.plotting.xlabeltext = 't (year)';
                handles.main.plotting.timefactor = 525969; % Sidereal year. I had to pick one.
        end
        xlabeltext = handles.main.plotting.xlabeltext;
        timefactor = handles.main.plotting.timefactor;
    else
        rethrow(Err)
    end
end

for nh=handles.main.activenh.on
    plot(handles.main.axes1,qcmt./timefactor,delf(:,nh)/nh,symbols{nh},'color',colors{nh})
    plot(handles.main.axes2,qcmt./timefactor,delg(:,nh),symbols{nh}, 'color',colors{nh})
    legendtext=[legendtext,legends{nh}];
    hold(handles.main.axes1,'on')
    hold(handles.main.axes2,'on')
end

if get(handles.main.loggamma,'value')
    set(handles.main.axes2,'yscale','log')
else
    set(handles.main.axes2,'yscale','linear')
end
if get(handles.main.log,'value')
    set([handles.main.axes1 handles.main.axes2],'xscale','log')
else
    set([handles.main.axes1 handles.main.axes2],'xscale','linear')
end

xlabel(handles.main.axes1, xlabeltext)
ylabel(handles.main.axes1,'{\Delta}f/n (Hz)')
legend(handles.main.axes1, legendtext,'location','best')

xlabel(handles.main.axes2, xlabeltext)
ylabel(handles.main.axes2,'\Delta\Gamma (Hz)')
legend(handles.main.axes2, legendtext,'location','best')

xlim(handles.main.axes1,xlimits1);
xlim(handles.main.axes2,xlimits2);

% set(handles.main.axes1,'xlimmode','auto')
set(handles.main.axes1,'ylimmode','auto')
% set(handles.main.axes2,'xlimmode','auto')
set(handles.main.axes2,'ylimmode','auto')

guidata(hObject, handles);
linkaxes([handles.main.axes1, handles.main.axes2], 'x')


% --- Executes on button press in view.
function view_Callback(hObject, eventdata, handles)
time = str2num(get(handles.currentpoint, 'string'));

[~, righttimeidx] = min(abs(handles.cond.time-time));
righttime = handles.main.din.qcmt(righttimeidx);
set(handles.currentpoint, 'string', num2str(righttime, '%0.3f'));

%Goes through each harmonic and plots the spectra.
cla(handles.nh1axes)
cla(handles.nh3axes)
cla(handles.nh5axes)
cla(handles.nh7axes)
set(handles.nh1axes,'fontsize',10)
set(handles.nh3axes,'fontsize',10)
set(handles.nh5axes,'fontsize',10)
set(handles.nh7axes,'fontsize',10)
foundtypes = logical([1 1 1 1]);

duds = 0;
for i = 1:2:max(handles.main.activenh.on)
    data = handles.cond.spectra{righttimeidx, i};
    
    axes(handles.(['nh' num2str(i) 'axes']));
    cla reset
    hold on
    
    if ~isempty(data)
        if get(handles.susceptance, 'value')
            [dataplot, dataline1, dataline2]  = plotyy(handles.(['nh' num2str(i) 'axes']), data(:,1), data(:,2), data(:,1), data(:,3));
            [fitplot, fitline1, fitline2] = plotyy(handles.(['nh' num2str(i) 'axes']), data(:,1), data(:,4), data(:,1), data(:,5));
            set([dataline1 dataline2], {'color'}, {'r'; 'b'})
            set([fitline1 fitline2], 'color', 'k')
            set([dataplot fitplot], 'ycolor', 'k')
            set([fitplot dataplot], 'fontsize', 10)
            set(fitplot(2), 'ytick', [], 'xtick', [])
        else
            plot(handles.(['nh' num2str(i) 'axes']), data(:,1), data(:,2), 'r', 'LineWidth', 2);
            plot(handles.(['nh' num2str(i) 'axes']), data(:,1),data(:,4),  'k');
        end
        
        if isnan(handles.main.din.cleandelf(righttimeidx,i))
            set(handles.(['remove' num2str(i)]), 'Visible', 'off')
            if get(handles.goodfits, 'value') == 1
                cla
                duds = duds+1;
            end
        else
            set(handles.(['remove' num2str(i)]), 'Visible', 'on')
            set(handles.(['remove' num2str(i)]), 'ForegroundColor', [0 0 0])
        end      
    else
        set(handles.(['remove' num2str(i)]), 'Visible', 'off')
        duds = duds+1;
    end
end

if handles.main.din.qcmt(righttimeidx) < 100000
    handles.currspectime = round(handles.main.din.qcmt(righttimeidx)*1000)/1000;
else
    handles.currspectime = round(handles.main.din.qcmt(righttimeidx)*100)/100;
end
handles.currenttimeidx = righttimeidx;
plotcond(hObject, handles)
guidata(hObject,handles)

if get(handles.goodfits, 'value') == 1 && duds == 3
    if handles.direction == 1
        nextpoint_Callback(hObject, eventdata, handles)
    else
        prevpoint_Callback(hObject, eventdata, handles)
    end
end

% --- Executes on button press in nextpoint.
function nextpoint_Callback(hObject, eventdata, handles)
index = handles.currenttimeidx;

if length(handles.cond.time) > index
    set(handles.currentpoint, 'string', num2str(handles.cond.time(index + 1), '%0.3f'));
else
    warndlg('That''s the last point!')
    return
end

handles.direction = 1;
handles.currenttimeidx = index + 1;
view_Callback(hObject, eventdata, handles);
plotcond(hObject, handles);

% --- Executes on button press in prevpoint.
function prevpoint_Callback(hObject, eventdata, handles)
index = handles.currenttimeidx;

if index > 1
   set(handles.currentpoint, 'string', num2str(handles.cond.time(index - 1),'%0.3f'));
else
    warndlg('That''s the first point!')
    return
end
handles.direction = 0;
handles.currenttimeidx = index - 1;
view_Callback(hObject, eventdata, handles);
plotcond(hObject, handles);

% --- Executes on button press in goodfits.
function goodfits_Callback(hObject, eventdata, handles)
view_Callback(hObject, eventdata, handles);

function currentpoint_Callback(hObject, eventdata, handles)
view_Callback(hObject, eventdata, handles);

function figure_WindowKeyPressFcn(hObject, eventdata, handles)
if strcmp(eventdata.Key, 'leftarrow')
        prevpoint_Callback(hObject, eventdata, handles)
elseif strcmp(eventdata.Key, 'rightarrow')
        nextpoint_Callback(hObject, eventdata, handles)
end

function susceptance_Callback(hObject, eventdata, handles)
view_Callback(hObject, eventdata, handles)

function F_conductance = lfun4c(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
%p(5): Conductance offset value
%p(6): Susceptance offset value (unused in this function)
if length(p) == 6
    F_conductance = lorentzcond(p(1:4), x) + p(5);
elseif length(p) == 10
    F_conductance = lorentzcond(p(1:4),x) + lorentzcond(p(5:8), x) + p(9);
elseif length(p) == 14
    F_conductance = lorentzcond(p(1:4),x) + lorentzcond(p(5:8),x) + lorentzcond(p(9:12),x) + p(13);
else
    F_conductance = 0;
end

function F_both = lfun4cs(p,x)
F_both = [lfun4c(p,x) lfun4s(p,x)];

function [fitdata, spectra] = fitmultiplepeaks(spectra, condsus, fitrange)
if exist('fitrange', 'var')
    fullfreq = spectra(:,1);
    freq = spectra(fitrange,1);
    conductance = spectra(fitrange,2);
    susceptance = spectra(fitrange,3);
else
    fullfreq = spectra(:,1);
    freq = spectra(:,1);
    conductance = spectra(:,2);
    susceptance = spectra(:,3);
end
    
if mean(freq) < 24e6 && mean(freq) > 16e6
    fitdata = [0 0 0 0 0];
    disp('The peak is the one at 21')
    return
end

%Find the peaks of interest.
[peak_detect, index, numpeaks] = findrelavantpeaks(freq, conductance);
if numpeaks == 0
    fitdata = [0 0 0 0 0];
    disp('No peaks were found')
    return
end

%Takes the peaks of interest and fits them
[pmul, I, G_parameters] = fitindividualpeaks([conductance, susceptance], freq, peak_detect, index, condsus);
if isempty(G_parameters)
    disp('No solution was found')
    fitdata = [0 0 0 0 0];
    return
end

Imul = min(I{1}):max(I{end});
%To help the fitting, I applied a simple smoothing function over 5. I don't
%remember what happened if I didn't do this.
smoothcond = smooth(conductance,5);
smoothsus = smooth(susceptance,5);

[~, ~, G_parametersmul] = fit_spectra(pmul, freq, [smoothcond smoothsus], Imul, condsus);
for i = 1:numpeaks
    [~, ~, G_parametersmul] = fit_spectra(G_parametersmul, freq, [smoothcond smoothsus], Imul, condsus);
end
condspec = lfun4c(G_parametersmul, fullfreq);
susspec = lfun4s(G_parametersmul, fullfreq);

%This code can be used if you want to compare results from the left peak
%and the tallest peak. "imptpeak" used to be an input value.
% if strcmp(imptpeak, 'left')
%     fitdata = G_parametersmul(1:5);
% elseif strcmp(imptpeak, 'tallest')
%     [~, idx] = max(peak_detect);
%     fitdata = G_parametersmul(5*(idx)-4:5*idx);
% end

fitdata = [G_parametersmul(1:4) G_parametersmul(end-1) G_parametersmul(end)];
spectra(:,4) = condspec;
spectra(:,5) = susspec;


function [pmul, I, G_parameters] = fitindividualpeaks(data, freq, peak_detect, index, condsus)
numpeaks = length(index);
conductance = data(:,1);
susceptance = data(:,2);
phi=0;%Assume rotation angle is 0
factor_range_fit=3;

%Reduce the noise a bit in fitting.
smoothcond = smooth(conductance,5);

baseline = linspace(conductance(1), conductance(end), length(freq))';
coffset = mean(baseline); %conductance offset
soffset = mean([susceptance(1), susceptance(end)]); %susceptance offset
Gmax = peak_detect-baseline(index); %maximum in G relative to the baseline
f0 = freq(index); %Frequency at the identified peak
halfg = (Gmax)./2+baseline(index); %value of the half max

if numpeaks == 1
    %Assumes a symmetric peak and finds the frequency at the half max of
    %one side
    halfg_freq = freq(find(abs(halfg-conductance)==min(abs((halfg-conductance))),1));
    gamma0 = abs(halfg_freq-f0); %calculates the half width at the half max
    %uses the half width to calculate the range to fit over
    I{1}=find(freq>=(f0-gamma0*factor_range_fit)&freq<=f0 + gamma0*factor_range_fit);
    p{1} = [f0 gamma0 phi Gmax coffset soffset]; %starting values for the fit
elseif numpeaks == 2
    [~, trough] = min(smoothcond(index(1):index(2))); %finds the idx for the minimum between the two peaks
    trough = trough + index(1); %corrects for the fact that min only looked at part of the whole
    % This checks if the trough between the peaks is very shallow, and if it is, it
    % defines the index for the trough to be halfway between the two peaks.
    if abs(smoothcond(trough)-smoothcond(index(1)))<.1*(abs(smoothcond(index(1))-smoothcond(index(2))))
        trough = floor(mean(index));
    end
    halfg_freq(1) = freq(find(abs(halfg(1)-conductance(1:index(1)))==min(abs((halfg(1)-conductance(1:index(1))))),1));
    [~, idx] = min(abs(halfg(2)-conductance(index(2):end)));
    try
        halfg_freq(2) = freq(idx+index(2));
    catch Err
        if strcmp(Err.identifier, 'MATLAB:badsubscript')
            halfg_freq(2) = freq(end);
        else
            rethrow(err)
        end
    end
    for i = 1:numpeaks
        gamma0(i) = abs(halfg_freq(i)-f0(i));
        p{i} = [f0(i) gamma0(i) phi Gmax(i) coffset soffset];
    end
    I{1}=find(freq>=(f0(1)-gamma0(1)*factor_range_fit)&freq<=freq(trough));
    I{2}=find(freq>=freq(trough)&freq <= f0(2) + gamma0(2)*factor_range_fit);
elseif numpeaks == 3
    halfg_freq(1) = freq(find(abs(halfg(1)-conductance(1:index(1)))==min(abs((halfg(1)-conductance(1:index(1))))),1));
    for i = 1:numpeaks-1
        [~, trough(i)] = min(smoothcond(index(i):index(i+1)));
        trough(i) = trough(i)+ index(i);
    end
    for i = 2:numpeaks-1 %Ok, so when I wrote this I was clearly thinking
        %in terms that this could be expanded to more than three peaks, so
        %for now I'll leave it like this, even though the only value that
        %satisfies is 2.
        [val idx] = min(abs(halfg(i)-conductance(index(i):index(i+1))));
        halfg_freq(i) = freq(idx+index(i));
        I{i} = find(freq >=(freq(trough(i-1))) & freq <= freq(trough(i)));
    end
    [val idx] = min(abs(halfg(numpeaks)-conductance(index(numpeaks):end)));
    halfg_freq(numpeaks) = freq(idx+index(numpeaks)-1);
    for i = 1:numpeaks
        gamma0(i) = abs(halfg_freq(i)-f0(i));
        p{i} = [f0(i) gamma0(i) phi Gmax(i) coffset soffset];
    end
    %Define the ranges for the other two values, beginning and end.
    I{1}=find(freq>=(f0(1)-gamma0(1)*factor_range_fit)&freq<=freq(trough(1)));
    I{numpeaks} = find(freq >=(freq(trough(numpeaks-1))) & freq <= f0(numpeaks) + gamma0(numpeaks)*factor_range_fit);
end

assignin('base','p',p);
pmul = [];

for i = 1:numpeaks
    if p{i}(2) == 0
        p{i}(2) = 100;
    end
    [G_fit{i}, G_residual{i}, G_parameters{i}] = fit_spectra(p{i}, freq, [smoothcond susceptance], I{i}, condsus);
    [G_fit{i}, G_residual{i}, G_parameters{i}] = fit_spectra(G_parameters{i}, freq, [smoothcond susceptance], I{i}, condsus);
    [G_fit{i}, G_residual{i}, G_parameters{i}] = fit_spectra(G_parameters{i}, freq, [smoothcond susceptance], I{i}, condsus);
    pmul = [pmul G_parameters{i}(1:4)];
end
pmul = [pmul coffset soffset];

function [fitted_y,residual,parameters]=fit_spectra(x0,freq_data,y_data,I, condsus)%fit spectra to conductance curve
%This function takes the starting guess values ('guess_values'_) and fits a
%Lorentz curve to the the x_data and y_data. The variable 'guess_values'
%needs to be a 1x5 array. The designation foe each elements is as follows:
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
%Variables, 'lb' and 'ub', defines the lower and upper bound of each of the
%guess_paramters. Both 'lb' and 'ub' are 1x5 array.
lb(1:length(x0)) = -Inf; %Assigns the lower bound to the parameters to -Inf
ub(1:length(x0)) = Inf; %Assigns the upper bound to the parameters to Inf
ub(3:4:length(x0)) = 90; %Changes the phase angle upper bound to 90
options=optimset('display','off','tolfun',1e-10,'tolx',1e-10,'MaxIter',5000);

try
    if condsus == 1 %fit conductance and susceptance simultaneously
        [parameters resnorm residual exitflag]=lsqcurvefit(@lfun4cs,x0,freq_data(I),y_data(I,:),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
    elseif condsus == 2 %fit conductance and susceptance separately and average
        [parameterssus resnormsus residualsus exitflagsus]=lsqcurvefit(@lfun4s,x0,freq_data(I),y_data(I,2),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
        [parameterscond resnormcond residualcond exitflagcond]=lsqcurvefit(@lfun4c,x0,freq_data(I),y_data(I,1),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
        parameters = mean([parameterssus; parameterscond]);
        resnorm = mean([resnormsus; resnormcond]);
        residual = mean([residualsus; residualcond]);
        exitflag = [exitflagsus exitflagcond];
    else%fit conductance only
        [parameters resnorm residual exitflag]=lsqcurvefit(@lfun4c,x0,freq_data(I),y_data(I,1),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
    end
catch Err
    if strcmp(Err.identifier, 'optimlib:lsqncommon:ProblemNotHandled')
        parameters = [0 0 0 0 0; 0 0 0 0 0];
    else
        rethrow Err
    end
end
fitted_y=lfun4cs(parameters,freq_data);
residual=fitted_y-y_data;

function [peak_detect, index, numpeaks] = findrelavantpeaks(freq, conductance, varargin)
%So why I am splitting this out as its own function is that I want to be
%able to go back and try this step again if I don't like what it gave.
%There is a weird third peak that occasionally shows up when trying to fit
%a double peaked function, and it seems to have a negative phase angle
%(which I thought was illegal, but whatever)

%Since some of the plots have very sloped baselines that were interfering
%with my "qualifying" peaks determination, I decided to use a sloped
%baseline rather than a simple mean.
baseline = linspace(conductance(1), conductance(end), length(freq))';

%The smoothing function eliminates minor peaks due to noise.
smoothcond = smooth(smooth(conductance, 20),10);
[peak_detect,index]=findpeaks(smoothcond,'sortstr','descend');
%Eliminates where there is no peak, where the range is just noise, and
if isempty(peak_detect)
    disp('No peaks were found')
    peak_detect = [];
    index = [];
    numpeaks = 0;
    return
end
%The following determines the peaks that are reasonable to consider--those
%that are at least half the height of the main peak (that could be changed)

%The first part of this determines what the difference from the baseline
%is. The second part calculates the "height" of the tallest peak and sees
%how the others compare.
qualifying = abs((peak_detect-baseline(index))./baseline(index))>max(abs((peak_detect-baseline(index))./baseline(index))/3);
%This worked when I was using a simple mean as a baseline:
%qualifying = peak_detect>(max(peak_detect)-baseline)/3+baseline;

%Now include only those peak heights and indexes identifies as being tall
%enough.
peak_detect = peak_detect(qualifying);
index = index(qualifying);
if min(index) < 5
    disp('The peak is at the edge of the window. No tall peak found.')
    peak_detect = [];
    index = [];
    numpeaks = 0;
    return
end

%Right now they are in order of peak height, so if I use the varargin to
%limit the number, here is where to do it.
if ~isempty(varargin)
    peak_detect = peak_detect(1:varargin{1});
    index = index(1:varargin{1});
elseif length(index) > 3
    peak_detect = peak_detect(1:3);
    index = index(1:3);
end

%sorts the peaks in order from left to right--important for determining
%which are next to each other.
[index order] = sort(index);

%Sometimes the peaks identified are really close together--this happened at
%least once on a shoulder that was just becoming its own peak. So I don't
%want to consider as separate peaks peaks that are too close together.
peak_detect = peak_detect(order);
spacedpeaks = [logical(1) (diff(index)>10)'];
index = index(spacedpeaks);
peak_detect = peak_detect(spacedpeaks);
numpeaks = sum(spacedpeaks);

%Another problem is when one of the peaks found is actually in the trough
%between the two main peaks. Bizarre. But at the same time, I do not want
%to eliminate a third peak that may be between two taller ones, though I
%certainly hope that doesn't happen. On the other hand, why not make that
%assumption?
if numpeaks == 3
    absolutepeakheight = peak_detect - baseline(index);
    if absolutepeakheight(2)< absolutepeakheight(1) && absolutepeakheight(2) < absolutepeakheight(3)
        peak_detect = peak_detect([1,3]);
        index = index([1,3]);
        numpeaks = 2;
    end
end
    
function F_susceptance = lfun4s(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phase angle difference
%p(4): Gmax maximum conductance
%p(5): Conductance offset value %Unusd in this function
%p(6): Susceptance offset value
if length(p) == 6
    F_susceptance = lorentzsus(p(1:4), x) + p(6);
elseif length(p) == 10   
    F_susceptance = lorentzsus(p(1:4),x) + lorentzsus(p(5:8),x) + p(10);
elseif length(p) == 14
    F_susceptance = lorentzsus(p(1:4),x) + lorentzsus(p(5:8),x) + lorentzsus(p(9:12), x) + p(14);
end

function F_cond = lorentzcond(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
F_cond= p(4).*((((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*cosd(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*sind(p(3)));

function F_sus = lorentzsus(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
F_sus = -p(4).*(-(((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*sind(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*cosd(p(3)));

function refit(hObject, eventdata, handles, harmonic)
index = handles.currenttimeidx;
if isfield(handles, 'currspectime')
    plottime = handles.currspectime;
    if abs(plottime - str2num(get(handles.currentpoint, 'string'))) <= 0.0011
        spectra =  handles.cond.spectra{index, harmonic};
        %This part finds that section that is within the plot limits
        limits = handles.(['nh' num2str(harmonic) 'axes']).XLim;
        fitrange = find(spectra(:,1) > limits(1) & spectra(:,1) < limits(2));
       
        [fitdata newspectra] = fitmultiplepeaks(spectra, handles.refittypetoggle, fitrange);
        if fitdata(1)/harmonic < -3e6
            fitdata = [NaN NaN NaN NaN NaN];
        end
        handles.cond.spectra{index, harmonic} = newspectra;
        handles.main.din.cleandelf(index, harmonic) = fitdata(1) - handles.main.din.bare.offset.(['f' num2str(harmonic)]);
        handles.main.din.cleandelg(index, harmonic) = fitdata(2) - handles.main.din.bare.offset.(['g' num2str(harmonic)]);
        % Since this is data from a refitting, also change the complete data
        handles.main.din.delf(index, harmonic) = fitdata(1) - handles.main.din.bare.offset.(['f' num2str(harmonic)]);
        handles.main.din.delg(index, harmonic) = fitdata(2) - handles.main.din.bare.offset.(['g' num2str(harmonic)]);
        handles.main.din.absf(index, harmonic) = fitdata(1);
        handles.main.din.absg(index, harmonic) = fitdata(2);
        handles.cond.fits{index, harmonic} = fitdata;
    else
        disp('You haven''t updated the spectra')
        return 
    end
else
    disp('You haven''t updated the spectra')
    return    
end

set(handles.(['remove' num2str(harmonic)]), 'visible', 'on');
handles.refit = 1;
guidata(hObject,handles)
plotraw(hObject, eventdata, handles)
plotcond(hObject, handles)
view_Callback(hObject, eventdata, handles)

function refit1_Callback(hObject, eventdata, handles)
refit(hObject, eventdata, handles, 1)

function refit3_Callback(hObject, eventdata, handles)
refit(hObject, eventdata, handles, 3)

function refit5_Callback(hObject, eventdata, handles)
refit(hObject, eventdata, handles, 5)

function refit7_Callback(hObject, eventdata, handles)
refit(hObject, eventdata, handles, 7)

function currentpoint_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function close_Callback(hObject, eventdata, handles)
if handles.refit
    time = handles.cond.time;
    spectra = handles.cond.spectra;
    fits = handles.cond.fits;
    save([handles.main.din.qcmpath handles.main.din.filebase '_cond.mat'], 'time', 'spectra', 'fits');
    disp('The spectra changes have been saved!')
end
uiresume(handles.figure);

function accesshandles_Callback(hObject, eventdata, handles)
keyboard;

function starttime_Callback(hObject, eventdata, handles)

function starttime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function endtime_Callback(hObject, eventdata, handles)

function endtime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function executeaction_Callback(hObject, eventdata, handles)
% First, this needs to get the start and end times for the action
starttime = str2num(get(handles.starttime, 'string'));
endtime = str2num(get(handles.endtime, 'string'));
if isempty(starttime) || isempty(endtime)
    warning('Invalid times entered. Time must be a number')
end

% Retrieve the harmonics of interest
harmonics = [];
for i = 1:2:7
    if get(handles.(['harm' num2str(i*5) 'MHz']), 'value')
        harmonics = [harmonics i];
    end
end
% If there aren't any, say so.
if isempty(harmonics)
    warning('You do not have any harmonics selected')
    return
end

% Next, the times need to be correlated with the indexes.
time = handles.main.din.qcmt;
startidx = find(time>starttime, 1); %Finds first time in bounds
endidx = find(time > endtime, 1) - 1; %Finds last time in bounds
if isempty(endidx)
    endidx = length(time);
end

timeidx = startidx:endidx;
if get(handles.actionrefit, 'value')
    %Code for refitting
    for harmonic = harmonics
        calcpercent = 10;
        for index = timeidx
            spectra =  handles.cond.spectra{index, harmonic}; %retrieve spectra
            if ~isempty(spectra)
                [fitdata, newspectra] = fitmultiplepeaks(spectra, handles.refittypetoggle); %refit data
                handles.cond.spectra{index, harmonic} = newspectra;
                %Update both regular and clean data since this is a refitting
                handles.main.din.cleandelf(index, harmonic) = fitdata(1) - handles.main.din.bare.offset.(['f' num2str(harmonic)]);
                handles.main.din.cleandelg(index, harmonic) = fitdata(2) - handles.main.din.bare.offset.(['g' num2str(harmonic)]);
                % Since this is data from a refitting, also change the complete data
                handles.main.din.delf(index, harmonic) = fitdata(1) - handles.main.din.bare.offset.(['f' num2str(harmonic)]);
                handles.main.din.delg(index, harmonic) = fitdata(2) - handles.main.din.bare.offset.(['g' num2str(harmonic)]);
                handles.main.din.absf(index, harmonic) = fitdata(1);
                handles.main.din.absg(index, harmonic) = fitdata(2);
                handles.cond.fits{index, harmonic} = fitdata;
            end
            if index-startidx > calcpercent/100*length(timeidx)
                calcpercent = calcpercent + 10;
                set(handles.main.statusbar, 'string', [num2str(round((index-startidx)/(endidx - startidx) * 100)) ' percent done with ' num2str(harmonic*5) 'MHz'], 'BackgroundColor', 'green');
                drawnow;
            end
        end
        set(handles.main.statusbar, 'string',['100 percent done with ' num2str(harmonic*5) 'MHz'], 'BackgroundColor', 'green');
    end
    %
else % Since it is a button group, the only other option is to be removing.
    for harmonic = harmonics %for selected harmonics
        handles.main.din.cleandelf(startidx:endidx, harmonic) = NaN;
        handles.main.din.cleandelg(startidx:endidx, harmonic) = NaN;
        handles.main.din.absf(startidx:endidx, harmonic) = NaN;
        handles.main.din.absg(startidx:endidx, harmonic) = NaN;
    end
end

guidata(hObject,handles) % save data
plotraw(hObject, eventdata, handles) %replot
plotcond(hObject, handles)
view_Callback(hObject, eventdata, handles)

function removeall_Callback(hObject, eventdata, handles)
handles = removepoint(hObject, eventdata, handles, 1);
handles = removepoint(hObject, eventdata, handles, 3);
handles = removepoint(hObject, eventdata, handles, 5);
handles = removepoint(hObject, eventdata, handles, 7);
guidata(hObject,handles) % save data


function harm5MHz_Callback(hObject, eventdata, handles)
function harm15MHz_Callback(hObject, eventdata, handles)
function harm25MHz_Callback(hObject, eventdata, handles)
function harm35MHz_Callback(hObject, eventdata, handles)

function refitremove_SelectionChangeFcn(hObject, eventdata, handles)


% --- Executes when selected object is changed in refittype.
function refittype_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in refittype 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject, 'tag'), 'refitcond')
    handles.refittypetoggle = 0;
elseif strcmp(get(hObject, 'tag'), 'refitcondsus')
    handles.refittypetoggle = 1;
else %refitcondsusavg
    handles.refittypetoggle = 2;
end
guidata(hObject,handles) % save data


% --- Executes on button press in savefile.
function savefile_Callback(hObject, eventdata, handles)
% hObject    handle to savefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
time = str2num(get(handles.currentpoint, 'string'));

[~, righttimeidx] = min(abs(handles.cond.time-time));
righttime = handles.main.din.qcmt(righttimeidx);

foundtypes = logical([1 1 1 1]);
duds = 0;

h = figure;
titles = {'1st', '2nd', '3rd', '4th', '5th'};
choice = questdlg(['Do you want to include fits?'], 'Include fits?',...
    'Yes', 'No', 'Yes');
switch choice
    case 'Yes'
        fits = 1;
    case 'No'
        fits = 0;
end

for i = 1:2:max(handles.main.activenh.on)
    data = handles.cond.spectra{righttimeidx, i};
    subplot(1,3,ceil(i/2))
    hold on       
        if get(handles.susceptance, 'value')
            if fits == 1
                [dataplot, dataline1, dataline2]  = plotyy(data(:,1)/10^6, data(:,2), data(:,1)/10^6, data(:,3));
                [fitplot, fitline1, fitline2] = plotyy(data(:,1)/10^6, data(:,4), data(:,1)/10^6, data(:,5));
                set([dataline1 dataline2], {'color'}, {'r'; 'b'})
                set([fitline1 fitline2], 'color', 'k')
                set([dataplot fitplot], 'ycolor', 'k')
            else
                dataplot = plot(data(:,1)/10^6, data(:,2));
                fitplot = plot(data(:,1)/10^6, data(:,4));
            end
            set([fitplot dataplot], 'fontsize', 10)
            set(fitplot(2), 'ytick', [], 'xtick', [])
        else
            plot(data(:,1)/10^6, data(:,2), 'r', 'LineWidth', 2);
            if fits == 1
                plot(data(:,1)/10^6,data(:,4),  'k');
            end
        end
        xlabel('Frequency (MHz)')
        ylabel('Conductance (S)')
        if i > 2
            title([titles{i} ' harmonic'])
        else
            title('Fundamental (5 MHz)')
        end
        
end
figure(h);
set(findall(gcf,'-property','FontSize'),'fontsize',20)
set(h, 'Units','centimeters', 'Position',[0 0 40 10])
set(h, 'PaperUnits','centimeters', 'PaperPosition',[0 0 40 10])
%hgsave([file(5:end) '_' strtok(times{i}, '.')])
[~, savepath] = uigetfile({'*.eps;*.png'}, 'Select the folder to save data in', 'C:\Users\Lauren\Desktop\Dropbox\Research - Sturdy\Thesis');
if fits == 1
    print(h, '-depsc2', [savepath handles.main.din.filebase '_' get(handles.currentpoint, 'string') '_fits.eps'])
    print(h, '-dpng', [savepath handles.main.din.filebase '_' get(handles.currentpoint, 'string') '_fits.png'])
else
    print(h, '-depsc2', [savepath handles.main.din.filebase '_' get(handles.currentpoint, 'string') '.eps'])
    print(h, '-dpng', [savepath handles.main.din.filebase '_' get(handles.currentpoint, 'string') '.png'])
end
pause(1)
close(h)
