function varargout = qcmplotter(varargin)
% QCMPLOTTER MATLAB code for qcmplotter.fig
%      QCMPLOTTER, by itself, creates a new QCMPLOTTER or raises the existing
%      singleton*.
%
%      H = QCMPLOTTER returns the handle to a new QCMPLOTTER or the handle to
%      the existing singleton*.
%
%      QCMPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QCMPLOTTER.M with the given input arguments.
%
%      QCMPLOTTER('Property','Value',...) creates a new QCMPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before qcmplotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to qcmplotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qcmplotter

% Last Modified by GUIDE v2.5 05-May-2014 13:34:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qcmplotter_OpeningFcn, ...
                   'gui_OutputFcn',  @qcmplotter_OutputFcn, ...
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

% --- Executes just before qcmplotter is made visible.
function qcmplotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to qcmplotter (see VARARGIN)

% Choose default command line output for qcmplotter
handles.output = hObject;
handles.folder = 'C:\Users\';
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = qcmplotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function importdata(hObject, eventdata, handles, num)
clear handles.data(num)
stringind = ['indicator' num2str(num)];
if ~isempty(get(handles.prefix, 'string'))
    prefix = get(handles.prefix, 'string');
name = get(handles.(['name' num2str(num)]), 'string');
filename = [prefix name];
else
    filename = get(handles.(['name' num2str(num)]), 'string');
end
    
try
    if isunix
        simplefile = [filename '.mat'];
        datafile = [filename '_data.mat'];
    else
        simplefile = [filename '.mat'];
        datafile = [filename '_data.mat'];
    end
newdata = load(datafile);
catch Err
    if strcmp(filename, prefix)
        set(handles.(stringind), 'backgroundcolor', 'default')
        return
    elseif strcmp(Err.identifier, 'MATLAB:load:couldNotReadFile')
        warndlg('The file could not be found. Please check that the filename is correct.')
        set(handles.(stringind), 'backgroundcolor', 'red')
        return
    end
end

if ~isfield(newdata, 'grhoep')
    warndlg('The saved files do not contain error data and cannot be processed.')
    set(handles.(stringind), 'backgroundcolor', 'red')
end
if ~isfield(newdata, 'delf')
    tempstruct = load(simplefile);
    newdata.delf = tempstruct.data.delf;
    newdata.delg = tempstruct.data.delg;
    newdata.time = tempstruct.data.time;
end
handles.data(num) = newdata;
stringind = ['indicator' num2str(num)];
stringcheck = ['checkbox' num2str(num)];
set(handles.(stringind), 'backgroundcolor', 'green')
set(handles.(stringcheck), 'value', 1)
guidata(hObject, handles);

function plotdata(hObject, eventdata, handles)
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinemarkersize',16)
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesbox','on')

if isfield(handles,'data')
    data = handles.data;
else
    warndlg('You haven''t read in any data yet!')
    return
end

plotting = whichtoplot(hObject, handles);
if isempty(plotting)
 warndlg('You forgot to select what to plot')
 return
end

mintime = str2num(get(handles.mintime, 'string'));
maxtime = str2num(get(handles.maxtime, 'string'));
numofpts = str2num(get(handles.numofpts, 'string'));

switch get(get(handles.timescale, 'SelectedObject'),'tag')
    case 'timemin'
        xlabeltext = 'Time (min.)';
        timecorr = 1;
    case 'timehr'
        xlabeltext = 'Time (hr.)';
        timecorr = 60;
    case 'timeday'
        xlabeltext = 'Time (days)';
        timecorr = 1440;
end
drhomin = str2num(get(handles.drholimmin, 'string'));
drhomax = str2num(get(handles.drholimmax, 'string'));
grhomin = str2num(get(handles.grholimmin, 'string'));
grhomax = str2num(get(handles.grholimmax, 'string'));
phimin = str2num(get(handles.philimmin, 'string'));
phimax = str2num(get(handles.philimmax, 'string'));
limits = {[drhomin drhomax] [grhomin grhomax] [phimin phimax]};

if get(handles.linlog, 'value')
    try
        times = linspace(mintime, maxtime, numofpts);
    catch Err
        if strcmp(Err.identifier, 'MATLAB:dimagree')
            warndlg('You must enter a start and end time')
        else
            return
        end
    end
else
    try
        times = logspace(log10(mintime+.01), log10(maxtime), numofpts);
    catch Err
        if strcmp(Err.identifier, 'MATLAB:dimagree')
            warndlg('You must enter a start and end time')
        else
            return
        end
    end
end


for j = plotting
    pointstoplot{j} = [];
    for i = times*timecorr
        [~, index] = min(abs(data(j).time - i));
        modulus = data(j).grhop;
        try
            while isnan(modulus(index,1)) && index<length(data(j).time)
                index = index + 1;
            end
            pointstoplot{j} = [pointstoplot{j} index];
        catch Err
        end
        
    end
end

symbols{1}='o';
symbols{2}='+';
symbols{4}='s';
symbols{3}='x';
symbols{5}='^';
symbols{6}='d';
symbols{7}='>';
symbols{8}='*';
symbols{9}='p';
symbols{10}='<';

outputaxes1=subplot(1,3,1);
hold on

for m=2
    for j = plotting
    errorbar(data(j).timep(pointstoplot{j},m)/timecorr,data(j).drhop(pointstoplot{j},m),data(j).drhoep(pointstoplot{j},m),symbols{j})
    end
end
xlabel(xlabeltext);
ylabel('d\rho (g/m^{2})');
ylim(limits{1});

outputaxes2=subplot(1,3,2);
hold on
for m=2
    for j = plotting
    errorbar(data(j).timep(pointstoplot{j},m)/timecorr,data(j).grhop(pointstoplot{j},m),data(j).grhoep(pointstoplot{j},m),symbols{j})
    end
end
xlabel(xlabeltext)
ylabel('|G^{*}|\rho (Pa-g/cm^{3})')
ylim(limits{2});

outputaxes3=subplot(1,3,3);
hold on
for m=2
    for j = plotting
    errorbar(data(j).timep(pointstoplot{j},m)/timecorr,data(j).phip(pointstoplot{j},m),data(j).phiep(pointstoplot{j},m),symbols{j})
    end
end

for i = plotting
    labels{i} = get(handles.(['name' num2str(i)]), 'string');
end

leg = legend(labels{plotting},'location','best');
set(leg, 'Interpreter', 'none');

xlabel(xlabeltext)
ylabel('\phi (deg.)')
ylim(limits{3});
p1 = get(outputaxes1,'Position');
p2 = get(outputaxes2,'Position');
p3 = get(outputaxes3,'Position');

p1(4) = p1(4)/1.5;
p2(4) = p2(4)/1.5;
p3(4) = p3(4)/1.5;
p1(2) = p1(2)*2;
p2(2) = p2(2)*2;
p3(2) = p3(2)*2;

set(outputaxes1,'position', p1/1.3)
set(outputaxes2, 'position', p2/1.3)
set(outputaxes3,'position', p3/1.3)

allaxes = [outputaxes1, outputaxes2, outputaxes3];
if get(handles.linlog, 'value') == 0
    set(allaxes, 'xscale', 'log')
end
if get(handles.logG, 'value') == 1
    set(outputaxes2, 'yscale', 'log')
end
set(leg, 'FontSize', 14)
%set(allaxes, 'Xtick', 0:6:36);
set(allaxes, 'Xlim', [mintime maxtime]);
hold off
handles.allaxes = allaxes;
if length(get(allaxes(1), 'xticklabel')) == 1
    set(allaxes, 'xtick', [10.^-2 10.^-1 10.^0 10.^1 10.^2 10.^3 10.^4 10.^5 10.^6 10.^7 10.^8 10.^9 10.^10]);
    if length(get(allaxes(1), 'xticklabel')) > 4
        set(allaxes, 'xtick', [10.^-1 10.^1 10.^3 10.^5 10.^7 10.^9]);
    end
end
guidata(hObject, handles)

function plotgphi(hObject, handles)
if isfield(handles,'data')
data = handles.data;
else
    warndlg('You haven''t read in any data yet!')
    return
end

h = figure;
hold on
plotting = whichtoplot(hObject, handles);

for j = plotting
    pointstoplot{j} = [];
    for i = 0:.25:20
        [~, index] = min(abs(data(j).phip - i));
        index = index(1);
        modulus = data(j).drhop;
        try
            while isnan(modulus(index,1)) && index<length(data(j).phip)
                index = index + 1;
            end
            pointstoplot{j} = [pointstoplot{j} index];
        catch Err
        end
        
    end
end

harmonics = 2;
symbols{1}='o';
symbols{2}='+';
symbols{4}='s';
symbols{3}='x';
symbols{5}='^';
symbols{6}='d';
symbols{7}='>';
symbols{8}='*';
symbols{9}='p';
symbols{10}='<';

for j = plotting
    lines(j) = plot(data(j).phip(pointstoplot{j},harmonics), data(j).grhop(pointstoplot{j},harmonics),symbols{j});
end
xlabel('\phi (deg.)')
ylabel('|G^{*}|\rho (Pa-g/cm^{3})')

for i = plotting
    labels{i} = get(handles.(['name' num2str(i)]), 'string');
end
set(gca, 'yscale', 'log')
leg = legend(labels{plotting},'location','northeast');
set(leg, 'Interpreter', 'none');
set(findall(gcf,'-property','FontSize'),'fontsize',20)
set(leg, 'FontSize', 14)

x = 1:.5:20;
for j = plotting
      [coeff stats] = polyfit(data(j).phip(pointstoplot{j},harmonics), log(data(j).grhop(pointstoplot{j}, harmonics)), 1) ;
   plot(x, exp((x.*coeff(1)+coeff(2)))); 
   slope(j) = coeff(1);
   intercept(j) = coeff(2);
end
slope
intercept

function saveimage_Callback(hObject, eventdata, handles)
newfig = figure('visible', 'off'); 
set(newfig, 'PaperUnits', 'inches', 'PaperSize', [15 5])
set(newfig, 'PaperPositionMode', 'manual')
set(newfig, 'PaperPosition', [0 0 15 5])
L1 = findobj(handles.figure1,'tag','legend');
allaxes = copyobj([handles.allaxes, L1], newfig);

set(allaxes(1), 'plotboxaspectratio', [1.1 1 1])
set(findall(gcf,'-property','FontSize'),'fontsize',16)
set(allaxes(2), 'fontsize', 10)
[FileName,PathName] = uiputfile('.eps');
print(newfig,'-depsc2',[PathName FileName])

function plotting = whichtoplot(hObject, handles)
%whichtoplot outputs an array with the length of the number of graphs
%desired, and values correxponding to the index of the desired data.
plotting = [];
%Loops through all of the checkbox callbacks to see which ones are checked.
for x = 1:6
    string = ['checkbox' num2str(x)];
    %If the checkbox is checked, the associated index is added to the
    %plotting arry.
    if get(handles.(string), 'value');
        plotting = [plotting x];
    end
end

function timescale_SelectionChangeFcn(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)

function axisscale_SelectionChangeFcn(hObject, eventdata, haneles)

function checkbox1_Callback(hObject, eventdata, handles)

function checkbox2_Callback(hObject, eventdata, handles)

function checkbox3_Callback(hObject, eventdata, handles)

function checkbox4_Callback(hObject, eventdata, handles)

function checkbox5_Callback(hObject, eventdata, handles)

function checkbox6_Callback(hObject, eventdata, handles)

function name1_Callback(hObject, eventdata, handles)
set(handles.indicator1, 'backgroundcolor', get(0,'defaultUicontrolBackgroundColor'))
importdata(hObject, eventdata, handles, 1)

function name1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function name2_Callback(hObject, eventdata, handles)
importdata(hObject, eventdata, handles, 2)


function name2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function name3_Callback(hObject, eventdata, handles)
importdata(hObject, eventdata, handles, 3)


function name3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function name4_Callback(hObject, eventdata, handles)
importdata(hObject, eventdata, handles, 4)


function name4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function name5_Callback(hObject, eventdata, handles)
importdata(hObject, eventdata, handles, 5)


function name5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function name6_Callback(hObject, eventdata, handles)
importdata(hObject, eventdata, handles, 6)


function name6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotbutton_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)

function numofpts_Callback(hObject, eventdata, handles)

function numofpts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function timehr_Callback(hObject, eventdata, handles)

function mintime_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)

function mintime_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxtime_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)
function maxtime_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function drholimmin_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)
function drholimmin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function drholimmax_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)
function drholimmax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function grholimmin_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)
function grholimmin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function grholimmax_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)
function grholimmax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function philimmax_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)
function philimmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function philimmin_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)
function philimmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function axesscale_SelectionChangeFcn(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)

function logG_Callback(hObject, eventdata, handles)
plotdata(hObject, eventdata, handles)

function plotgphi_Callback(hObject, eventdata, handles)
plotgphi(hObject, handles)

function prefix_Callback(hObject, eventdata, handles)

function prefix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
selection = questdlg('Do you want to save the file?',...
    'Close Request Function',...
    'Yes','No','Cancel','Yes');
switch selection,
    case 'Yes',
        handles=guidata(hObject);
        hgsave(hObject,'qcmplotter.fig')
            delete(hObject)
        
            return
      
    case 'No'
        delete(hObject)
    case 'Cancel'
        return
end
