function varargout = condfig2(varargin)
%CONDFIG M-file for condfig2.fig
%      CONDFIG2, by itself, creates a new CONDFIG2 or raises the existing
%      singleton*.
%
%      H = CONDFIG2 returns the handle to a new CONDFIG2 or the handle to
%      the existing singleton*.
%
%      CONDFIG2('Property','Value',...) creates a new CONDFIG2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CONDFIG_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CONDFIG2('CALLBACK') and CONDFIG2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CONDFIG.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help condfig

% Last Modified by GUIDE v2.5 28-Oct-2014 10:59:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @condfig2_OpeningFcn, ...
                   'gui_OutputFcn',  @condfig2_OutputFcn, ...
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
function condfig2_OpeningFcn(hObject, eventdata, handles, varargin)
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
    newY = parentPosition(2) + (parentPosition(4) - currentPosition(4)*.95);
    newW = currentPosition(3);
    newH = currentPosition(4);
    
    set(hObject, 'Position', [newX, newY, newW, newH]);
end
handles.currenttimeidx = 1;
handles.main = mainHandles;
handles.direction = 1;
handles.refit = 0;

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
    disp('   condfig2(''onelayerguiwithcond'', x)');
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
function varargout = condfig2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This passes the changes to the cleaned up version of the data back to the
%main function.
varargout{1} = handles.main.din;
delete(hObject);

function condfig2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function condfig2_CloseRequestFcn(hObject, eventdata, handles)
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
    else
        warning('You haven''t updated the spectra')
        return 
    end
else
    warning('You haven''t updated the spectra')
    return    
end
set(handles.(['remove' num2str(harmonic)]), 'visible', 'off');
guidata(hObject,handles)
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
qcmt = handles.main.din.qcmt;
%Define linestyle for the three different harmonics.
linestyle = [{'o-r'},{'o-b'},{'o-g'}];
%Goes through the three harmonics and plots them if desired. I couldn't
%figure out a way to do this such that the number was the index.

%Determines which spectra are present and plots them.
haveconductance = ~cellfun(@isempty, handles.cond.spectra); %logical 1 if spectra present
qcmts = [qcmt'; qcmt'; qcmt'; qcmt'; qcmt'];
harms = repmat([1,0,3,0,5], length(qcmt), 1);
plot(handles.main.axes1, qcmts(haveconductance(:,1:2:end)), delf(haveconductance)./harms(haveconductance), '.k');
plot(handles.main.axes2, qcmts(haveconductance(:,1:2:end)), delg(haveconductance), '.k');
idx = handles.currenttimeidx;
plot(handles.main.axes1, [qcmt(idx) qcmt(idx) qcmt(idx)], delf(idx, 1:2:5)./[1 3 5], '.c')
plot(handles.main.axes2, [qcmt(idx) qcmt(idx) qcmt(idx)], delg(idx, 1:2:5), '.c')
%Sets ylim back to the original so very far off points aren't shown.
%ylim(handles.main.axes1,ylimits1);
%ylim(handles.main.axes2,ylimits2);
xlim(handles.main.axes1,xlimits1);
xlim(handles.main.axes2,xlimits2);

function plotraw(hObject, eventdata, handles)
% Plots the raw data. This function should be identical to the plotraw
% function in onelayergui2.
qcmt = handles.main.din.qcmt;
if get(handles.main.cleaned, 'value')
    delf = handles.main.din.cleandelf;
    delg = handles.main.din.cleandelg;
else
    delf = handles.main.din.delf;
    delg = handles.main.din.delg;
end

%ylimits1 = ylim(handles.main.axes1);
%ylimits2 = ylim(handles.main.axes2);
xlimits1 = xlim(handles.main.axes1);
xlimits2 = xlim(handles.main.axes2);
% now we plot the data
cla(handles.main.axes2)
cla(handles.main.axes1)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinemarkersize',12)
set(0,'defaultlinelinewidth',2)
colors{1}=[1,0,0];
colors{3}=[0,0.5,0];
colors{5}=[0,0,1];
symbols{1}='r+-';
symbols{3}='go-';
symbols{5}='bs-';
legends{1}='n=1';
legends{3}='n=3';
legends{5}='n=5';
legendtext={};

try
    handles.activenh.on;
catch ErMsg
    if strcmp(ErMsg.identifier, 'MATLAB:nonExistentField')
        activenh.on=[];
        activenh.off=[];
        for nh=[1,3,5]
            if get(handles.main.(['n',num2str(nh)]),'value')
                activenh.on=[activenh.on,nh];
            else
                activenh.off=[activenh.off,nh];
            end
        end
        handles.activenh = activenh;
    else
        rethrow(ErMsg)
    end
end

%This plots the values as in the main function, but applies hold on and removes the
%hittest property, which means that the curves cannot be selected by the
%datacursor. This comes in handy later so that only the conductance spectra
%points can be selected.
for nh=handles.activenh.on
    plot(handles.main.axes1,qcmt,delf(:,nh)/nh,symbols{nh}, 'color', colors{nh}, 'HitTest', 'off')
    plot(handles.main.axes2,qcmt,delg(:,nh),symbols{nh}, 'color', colors{nh}, 'HitTest', 'off')
    legendtext=[legendtext,legends{nh}];
    hold(handles.main.axes1,'on')
    hold(handles.main.axes2,'on')
end

set(handles.main.axes2,'yscale','log')
xlabel(handles.main.axes1,get(handles.main.xlabeltext,'string'))
ylabel(handles.main.axes1,'\Delta f/n (Hz)')
legend(handles.main.axes1, legendtext,'location','best')

xlabel(handles.main.axes2,get(handles.main.xlabeltext,'string'))
ylabel(handles.main.axes2,'\Delta\Gamma (Hz)')
legend(handles.main.axes2, legendtext,'location','best')

%ylim(handles.main.axes1,[-Inf Inf]);
%ylim(handles.main.axes2,[-Inf Inf]);
xlim(handles.main.axes1,xlimits1);
xlim(handles.main.axes2,xlimits2);
if get(handles.main.log,'value')
    set(handles.main.axes1,'xscale','log')
    set(handles.main.axes2,'xscale','log')
else
    set(handles.main.axes1,'xscale','linear')
    set(handles.main.axes2,'xscale','linear')
end

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
set(handles.nh1axes,'fontsize',12)
set(handles.nh3axes,'fontsize',12)
set(handles.nh5axes,'fontsize',12)
foundtypes = logical([1 1 1]);

duds = 0;
for i = 1:2:min(size(handles.cond.spectra))
    data = handles.cond.spectra{righttimeidx, i};
    
    axes(handles.(['nh' num2str(i) 'axes']));
    cla reset
    hold on
    
    if ~isempty(data)
        if get(handles.susceptance, 'value')
            [dataplot dataline1 dataline2]  = plotyy(handles.(['nh' num2str(i) 'axes']), data(:,1), data(:,2), data(:,1), data(:,3));
            [fitplot fitline1 fitline2] = plotyy(handles.(['nh' num2str(i) 'axes']), data(:,1), data(:,4), data(:,1), data(:,5));
            set([dataline1 dataline2], {'color'}, {'r'; 'b'})
            set([fitline1 fitline2], 'color', 'k')
            set([dataplot fitplot], 'ycolor', 'k')
            set([fitplot dataplot], 'fontsize', 12)
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
%p(5): Offset value
F_conductance= p(4).*((((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*cosd(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*sind(p(3)))+p(5);

function F_conductance_2 = lfun4c2(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
F_conductance_2 = p(4).*((((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*cosd(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*sind(p(3)))+p(5)+...
    p(4+5).*((((x.^2).*((2.*p(2+5)).^2))./(((((p(1+5)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2+5)).^2)))).*cosd(p(3+5))-((((p(1+5)).^2-x.^2)).*x.*(2.*p(2+5)))./...
    (((((p(1+5)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2+5)).^2))).*sind(p(3+5)))+p(5+5);

function F_conductance_2 = lfun4c3(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
F_conductance_2 = p(4).*((((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*cosd(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*sind(p(3)))+p(5)+...%end of 1st term
    p(4+5).*((((x.^2).*((2.*p(2+5)).^2))./(((((p(1+5)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2+5)).^2)))).*cosd(p(3+5))-((((p(1+5)).^2-x.^2)).*x.*(2.*p(2+5)))./...
    (((((p(1+5)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2+5)).^2))).*sind(p(3+5)))+p(5+5)+...%end of 2nd term
    p(4+10).*((((x.^2).*((2.*p(2+10)).^2))./(((((p(1+10)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2+10)).^2)))).*cosd(p(3+10))-((((p(1+10)).^2-x.^2)).*x.*(2.*p(2+10)))./...
    (((((p(1+10)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2+10)).^2))).*sind(p(3+10)))+p(5+10);%end of 3rd term

function [fitdata spectra] = fitmultiplepeaks(spectra)
freq = spectra(:,1);
conductance = spectra(:,2);
susceptance = spectra(:,3);
if mean(freq) < 24e6 && mean(freq) > 16e6
    fitdata = [NaN NaN NaN NaN NaN];
    spectra(:,4:5) = NaN;
    disp('The peak is the one at 21')
    return
end

%Find the peaks of interest.
[peak_detect, index, numpeaks] = findrelavantpeaks(freq, conductance);
if numpeaks == 0
    fitdata = [NaN NaN NaN NaN NaN];
    spectra(:,4:5) = NaN;
    disp('No peaks were found')
    return
end
[pmul, I, G_parameters] = fitindividualpeaks(conductance, freq, peak_detect, index);

if isempty(G_parameters)
    disp('No solution was found')
    fitdata = [NaN NaN NaN NaN NaN];
    spectra(:,4:5) = NaN;
    return
end

Imul = min(I{1}):max(I{end});
if length(Imul)<50
    disp('Insufficient points for good fit')
    fitdata = [NaN NaN NaN NaN NaN];
    spectra(:,4:5) = NaN;
    return
end

smoothcond = smooth(conductance,5);

if numpeaks == 1
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra(pmul, freq, smoothcond, Imul);
    plot(freq, lfun4c(G_parametersmul, freq), 'k')
    condspec = lfun4c(G_parametersmul, freq);
    susspec = lfun4s(G_parametersmul, freq);
elseif numpeaks == 2
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra2(pmul, freq, smoothcond, Imul);
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra2(G_parametersmul, freq, smoothcond, Imul);
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra2(G_parametersmul, freq, smoothcond, Imul);
    plot(freq, lfun4c2(G_parametersmul, freq), 'k')
    condspec = lfun4c2(G_parametersmul, freq);
    susspec = lfun4s2(G_parametersmul, freq);
elseif numpeaks == 3
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra3(pmul, freq, smoothcond, Imul);
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra3(G_parametersmul, freq, smoothcond, Imul);
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra3(G_parametersmul, freq, smoothcond, Imul);
    [G_fitmul, G_residualmul, G_parametersmul] = fit_spectra3(G_parametersmul, freq, smoothcond, Imul);
    plot(freq, lfun4c3(G_parametersmul, freq), 'k')
    condspec = lfun4c3(G_parametersmul, freq);
    susspec = lfun4s3(G_parametersmul, freq);
end

fitdata = G_parametersmul(1:5); %Data for leftmost peak

spectra(:,4) = condspec;
spectra(:,5) = susspec;


function [pmul, I, G_parameters] = fitindividualpeaks(conductance, freq, peak_detect, index) 
numpeaks = length(index);

phi=0;%Assume rotation angle is 0
offset=0;%Assume offset value is 0
factor_range_fit=3;

%Reduce the noise a bit in fitting.
smoothcond = smooth(conductance,5);

baseline = linspace(conductance(1), conductance(end), length(freq))';

Gmax = peak_detect;
f0 = freq(index);
halfg = (Gmax-baseline(index))./2+baseline(index);

if numpeaks == 1
    [~, idx] = min(abs(halfg-conductance));
    halfg_freq = freq(find(abs(halfg-conductance)==min(abs((halfg-conductance))),1));
    gamma0 = abs(halfg_freq-f0);
    I{1}=find(freq>=(f0-gamma0*factor_range_fit)&freq<=f0 + gamma0*factor_range_fit);
    p{1} = [f0 gamma0 phi Gmax offset];
    Imul = min(I{1}:max(I{end}));
elseif numpeaks == 2
    [~, trough] = min(smoothcond(index(1):index(2)));
    trough = trough + index(1);
    if abs(smoothcond(trough)-smoothcond(index(1)))<.1*(abs(smoothcond(index(1))-smoothcond(index(2))))
        trough = floor(mean(index));
    end
    [~, idx] = min(abs(halfg(1)-conductance(1:index(1))));
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
        p{i} = [f0(i) gamma0(i) phi Gmax(i) offset];
        pmul(5*i-4:5*i) = [f0(i) gamma0(i) phi Gmax(i) offset];
    end
    I{1}=find(freq>=(f0(1)-gamma0(1)*factor_range_fit)&freq<=freq(trough));
    I{2}=find(freq>=freq(trough)&freq <= f0(2) + gamma0(2)*factor_range_fit);
    Imul = min(I{1}:max(I{end}));
elseif numpeaks == 3
    [~, idx] = min(abs(halfg(1)-conductance(1:index(1))));
    halfg_freq(1) = freq(find(abs(halfg(1)-conductance(1:index(1)))==min(abs((halfg(1)-conductance(1:index(1))))),1));
    for i = 1:numpeaks-1
        [~, trough(i)] = min(smoothcond(index(i):index(i+1)));
        trough(i) = trough(i)+ index(i);
    end
   
    for i = 2:numpeaks-1
        [val idx] = min(abs(halfg(i)-conductance(index(i):index(i+1))));
        halfg_freq(i) = freq(idx+index(i));
        I{i} = find(freq >=(freq(trough(i-1))) & freq <= freq(trough(i)));
    end
    [val idx] = min(abs(halfg(numpeaks)-conductance(index(numpeaks):end)));
    halfg_freq(numpeaks) = freq(idx+index(numpeaks)-1);
    for i = 1:numpeaks
        gamma0(i) = abs(halfg_freq(i)-f0(i));
        p{i} = [f0(i) gamma0(i) phi Gmax(i) offset];
        pmul(5*i-4:5*i) = [f0(i) gamma0(i) phi Gmax(i) offset];
    end
    I{1}=find(freq>=(f0(1)-gamma0(1)*factor_range_fit)&freq<=freq(trough(1)));
    I{numpeaks} = find(freq >=(freq(trough(numpeaks-1))) & freq <= f0(numpeaks) + gamma0(numpeaks)*factor_range_fit);
% else
%     pmul = [0 0 0 0 0];
%     I = [];
%     G_parameters = [];
%     return
end

%assignin('base','p',p);
pmul = [];

for i = 1:numpeaks
    if p{i}(2) == 0
        p{i}(2) = 100;
    end
    [G_fit{i}, G_residual{i}, G_parameters{i}] = fit_spectra(p{i}, freq, smoothcond, I{i});
    [G_fit{i}, G_residual{i}, G_parameters{i}] = fit_spectra(G_parameters{i}, freq, smoothcond, I{i});
    [G_fit{i}, G_residual{i}, G_parameters{i}] = fit_spectra(G_parameters{i}, freq, smoothcond, I{i});
    pmul = [pmul G_parameters{i}];
end

function [fitted_y,residual,parameters]=fit_spectra(x0,freq_data,y_data,I,lb,ub)%fit spectra to conductance curve
%This function takes the starting guess values ('guess_values'_) and fits a
%Lorentz curve to the the x_data and y_data. The variable 'guess_values' 
%needs to be a 1x5 array. The designation for each elements is as follows:
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
%Variables, 'lb' and 'ub', defines the lower and upper bound of each of the
%guess_paramters. Both 'lb' and 'ub' are 1x5 array.
if nargin==4
    lb=[-inf -inf -inf -Inf -Inf];
    ub=[Inf Inf 90 Inf Inf];
end%if nargin==3
options=optimset('display','off','tolfun',1e-10,'tolx',1e-10,'MaxIter',5000);
try
    [parameters resnorm residual exitflag]=lsqcurvefit(@lfun4c,x0,freq_data(I),y_data(I),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
catch Err
    Err
    ub = [Inf Inf Inf Inf Inf];
    [parameters resnorm residual exitflag]=lsqcurvefit(@lfun4c,x0,freq_data(I),y_data(I),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
end
    fitted_y=lfun4c(parameters,freq_data);
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

smoothcond = smooth(smooth(conductance, 20),20);
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
%This worked when I was using a simple mean as a baseline.
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
        peak_detect = peak_detect([1,3])
        index = index([1,3])
        numpeaks = 2;
    end
end


plot(freq, conductance, 'b')
hold on
plot(freq(index), peak_detect, 'g*')
    
function [fitted_y,residual,parameters]=fit_spectra2(x0,freq_data,y_data,I,lb,ub)%fit spectra to conductance curve
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
if nargin==4
    lb=[-inf -inf -inf -Inf -Inf -inf -inf -inf -Inf -Inf];
    ub=[Inf Inf 90 Inf Inf Inf Inf 90 Inf Inf];
end%if nargin==3
options=optimset('display','off','tolfun',1e-10,'tolx',1e-10,'MaxIter',5000);
[parameters resnorm residual exitflag]=lsqcurvefit(@lfun4c2,x0,freq_data(I),y_data(I),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
fitted_y=lfun4c2(parameters,freq_data);
residual=fitted_y-y_data;

function [fitted_y,residual,parameters]=fit_spectra3(x0,freq_data,y_data,I,lb,ub)%fit spectra to conductance curve
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
if nargin==4
    lb=[-inf -inf -inf -Inf -Inf -inf -inf -inf -Inf -Inf -Inf -inf -inf -inf -Inf];
    ub=[Inf Inf 90 Inf Inf Inf Inf 90 Inf Inf Inf Inf 90 Inf Inf];
end%if nargin==3
options=optimset('display','off','tolfun',1e-10,'tolx',1e-10,'MaxIter',5000);
[parameters resnorm residual exitflag]=lsqcurvefit(@lfun4c3,x0,freq_data(I),y_data(I),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
fitted_y=lfun4c3(parameters,freq_data);
residual=fitted_y-y_data;
% disp('Conductance fitted parameters:');
% disp(parameters');
% exitflag

function [fitted_y,residual,parameters]=fit_spectra_sus(x0,freq_data,susceptance_data,I,lb,ub)%fit spectra to susceptance curve
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
if nargin==4
    lb=[0 0 -90 -inf -Inf];
    ub=[Inf Inf 90 Inf Inf];
end%if nargin==3
options=optimset('display','off','tolfun',1e-100,'tolx',1e-100,'MaxIter',1000);
[parameters resnorm residual exitflag]=lsqcurvefit(@lfun4s,x0,freq_data(I),susceptance_data(I),lb,ub,options);%use lsqcurvefit function to fit the spectra data to a Lorentz curve
fitted_y=lfun4s(parameters,freq_data);
residual=fitted_y-susceptance_data;
exitflag
% disp('Susceptance fitted parameters:');
% disp(parameters');

function F_susceptance_2 = lfun4s(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phase angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
F_susceptance_2 = -p(4).*(-(((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*sind(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*cosd(p(3)))+p(5);

function F_susceptance_2 = lfun4s2(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phase angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
F_susceptance_2 = -p(4).*(-(((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*sind(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*cosd(p(3)))+p(5)...
    -p(4+5).*(-(((x.^2).*((2.*p(2+5)).^2))./(((((p(1+5)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2+5)).^2)))).*sind(p(3+5))-((((p(1+5)).^2-x.^2)).*x.*(2.*p(2+5)))./...
    (((((p(1+5)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2+5)).^2))).*cosd(p(3+5)))+p(5+5);

function F_susceptance_2 = lfun4s3(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phase angle difference
%p(4): Gmax maximum conductance
%p(5): Offset value
F_susceptance_2 = -p(4).*(-(((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*sind(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*cosd(p(3)))+p(5)...
    -p(4+5).*(-(((x.^2).*((2.*p(2+5)).^2))./(((((p(1+5)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2+5)).^2)))).*sind(p(3+5))-((((p(1+5)).^2-x.^2)).*x.*(2.*p(2+5)))./...
    (((((p(1+5)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2+5)).^2))).*cosd(p(3+5)))+p(5+5)...
    -p(4+10).*(-(((x.^2).*((2.*p(2+10)).^2))./(((((p(1+10)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2+10)).^2)))).*sind(p(3+10))-((((p(1+10)).^2-x.^2)).*x.*(2.*p(2+10)))./...
    (((((p(1+10)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2+10)).^2))).*cosd(p(3+10)))+p(5+10);

function refit(hObject, eventdata, handles, harmonic)
index = handles.currenttimeidx;
plottime = handles.currspectime;
if isfield(handles, 'currspectime')
    if abs(plottime - str2num(get(handles.currentpoint, 'string'))) <= 0.0011
        spectra =  handles.cond.spectra{index, harmonic};
        [fitdata newspectra] = fitmultiplepeaks(spectra);
        if fitdata(1)/harmonic < -3e6
            fitdata = [NaN NaN NaN NaN NaN];
        end
        handles.cond.spectra{index, harmonic} = newspectra;
        handles.main.din.cleandelf(index, harmonic) = fitdata(1) - handles.main.offset.(['f' num2str(harmonic)]);
        handles.main.din.cleandelg(index, harmonic) = fitdata(2) - handles.main.offset.(['g' num2str(harmonic)]);
        % Since this is data from a refitting, also change the complete
        % data
        handles.main.din.delf(index, harmonic) = fitdata(1);
        handles.main.din.delg(index, harmonic) = fitdata(2);
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


function currentpoint_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
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
for i = 1:2:5
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
                [fitdata, newspectra] = fitmultiplepeaks(spectra); %refit data
                handles.cond.spectra{index, harmonic} = newspectra;
                %Update both regular and clean data since this is a refitting
                handles.main.din.cleandelf(index, harmonic) = fitdata(1) - handles.main.offset.(['f' num2str(harmonic)]);
                handles.main.din.cleandelg(index, harmonic) = fitdata(2) - handles.main.offset.(['g' num2str(harmonic)]);
                handles.main.din.delf(index, harmonic) = fitdata(1);
                handles.main.din.delg(index, harmonic) = fitdata(2);
                handles.cond.fits{index, harmonic} = fitdata;
            end
            if index-startidx > calcpercent/100*length(timeidx)
                calcpercent = calcpercent + 10;
                disp([num2str(round((index-startidx)/(endidx - startidx) * 100)) ' percent done with ' num2str(harmonic) 'MHz']);
            end
        end    
    end
    %
else % Since it is a button group, the only other option is to be removing.
    for harmonic = harmonics %for selected harmonics
        handles.main.din.cleandelf(startidx:endidx, harmonic) = NaN;
        handles.main.din.cleandelg(startidx:endidx, harmonic) = NaN;        
    end
end

guidata(hObject,handles) % save data
plotraw(hObject, eventdata, handles) %replot
plotcond(hObject, handles)
view_Callback(hObject, eventdata, handles)

function removeall_Callback(hObject, eventdata, handles)

function harm25MHz_Callback(hObject, eventdata, handles)

function harm15MHz_Callback(hObject, eventdata, handles)

function harm5MHz_Callback(hObject, eventdata, handles)

function refitremove_SelectionChangeFcn(hObject, eventdata, handles)
