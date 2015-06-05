function varargout = onelayergui3(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @onelayergui3_OpeningFcn, ...
    'gui_OutputFcn',  @onelayergui3_OutputFcn, ...
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

% --- Executes just before onelayergui3 is made visible.
function onelayergui3_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
date='6/25/14';
set(handles.mdate,'string',['m: ',date])
set(handles.harmonicadjustment,'string',num2str(1))
warning('off','MATLAB:Axes:NegativeDataInLogAxis')
set(hObject,'toolbar','figure');
% Update handles structure

%Sets constants for use in the program.
handles.constants.f1 = 5e6;
handles.constants.zq = 8.84e6;

handles.constants.nhvals = {[1,3,1,5] [1,3,3,5] [1,5,1,3] [1,5,5,3] [3,5,3,1] [3,5,5,1] [1,3,5,NaN] [1,5,3,NaN] [3,5,1,NaN]};
handles.constants.label = {'1:3,1' '1:3,3' '1:5,1' '1:5,5' '3:5,3' '3:5,5', '1:3,5', '1:5,3', '3:5,1'};
handles.saveplot.plot = 0; %default to not
guidata(hObject, handles);
set(gcf,'Pointer','arrow');

% --- Outputs from this function are returned to the command line.
function varargout = onelayergui3_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function changeqcmfile_Callback(hObject, eventdata, handles, qcmfile, qcmpath)
% This function imports data from a file, stores it (including raw f and g
% data, spectra, and error ranges if available) in handles, and displays plots. 

%sets the "show conductance" button to visible.
set(handles.opencond, 'visible', 'on')

if ~exist('qcmfile','var')
    %Opens dialog to select file to load. It will go to the last used folder by
    %default.
    if isfield(handles, 'din') && isfield(handles.din, 'qcmpath')
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
        [qcmfile, qcmpath]=uigetfile({'*.xls;*.xlsx;*.fre;*.mat'}, 'Select the data file');
    end
end

handles.restorefile = [qcmpath qcmfile];

% If someone cancels out of the uigetfile dialog, qcmfile and qcmpath will
% both be 0. This checks if that is the case.
if ~qcmfile
    warning('You did not select a new file. No data has been changed')
    return
end

set(handles.qcmfile,'string', ['Loading ' qcmfile])
drawnow; %means that the above change will be visible.

%clears the data currently stored in handles.din, as that is all dependent
%on the data file.
if isfield(handles, 'din')
    rmfield(handles, 'din');
end

if isfield(handles, 'saveplot')
    rmfield(handles, 'saveplot');
    handles.saveplot.plot = 0;
end

% Determines the file type and the base file name.
[~, filebase, filetypetxt] = fileparts(qcmfile);
name = {'' 'f1' 'g1' 'f3' 'g3' 'f5' 'g5'}; %Sets names for later.

if strcmp(qcmfile(end-7:end-4),'data')
    % Applied if the data has already been processed and just wants to be
    % viewed.
    handles.filetype = 4;
    load([qcmpath, qcmfile]);
    try
        qcmt = time;
        cleandelf = delf;
        cleandelg = delg;
    catch
        warning('The data file is of an old type and does not have the proper fields')
        return
    end
    try
        for nh = [1 3 5]
            fref=['f',num2str(nh)];
            gref=['g',num2str(nh)];
            set(handles.(['offset', fref]),'string',commanumber(offset.(fref)))
            set(handles.(['offset', gref]),'string',commanumber(offset.(gref)))
        end
        handles.din.bare.error.f = error.f;
        handles.din.bare.error.g = error.g;
    catch Err
        warning('Data file doesn''t contain offset or error information')
        disp('Default error values will be used')
        offset.f1=0;
        offset.f3=0;
        offset.f5=0;
        offset.g1=0;
        offset.g3=0;
        offset.g5=0;
        handles.din.bare.error.f = [44 0 120 0 273];
        handles.din.bare.error.g = [11 0 22 0 14];
    end
    datasize = size(dfcalcp);
    if datasize(3) < 5
        set(handles.n5, 'value', 0)
        handles.activenh = getactivenh(handles);
    end
    handles.overwrite = 0; %Disables saving of a new _data file later on.
    set(handles.loadsavedplot, 'visible', 'on')
    set(handles.savecleaned, 'visible', 'off')
else
    % The program can read .xls or .xlsx (older data--Garret's primarily), .fre
    % (QTZ raw data), and .mat (new QCM program). This program assums standard
    % output data for each of those types of files.
    handles.overwrite = 1;
    set(handles.loadsavedplot, 'visible', 'off')
    set(handles.savecleaned, 'visible', 'on')
    if strcmp(filetypetxt,'.fre')
        handles.filetype = 1;
        %.fre files have the reference frequencies at the top. These are read in
        %first.
        fid=fopen([qcmpath, qcmfile]); %Opens .fre file for reading
        qcmheader = textscan(fid,'%s %s %s %s %s %s %s',1); %gets top line
        fclose(fid); %closes file to save memory (and because it isn't necessary)
        
        %This provides the reference data. First it checks if there is a _bare
        %file.
        reference = getreference(handles, qcmpath, filebase, 1); %reference = 0 if unsuccessful
        if length(reference) > 1 %Uses _bare file if possible
            for n=2:7
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
        
    elseif strcmp(filetypetxt, '.mat')
        handles.filetype = 2;
        load([qcmpath, qcmfile]) %Load the mat file into memory.
        
        % Since the file is saved with a million rows and they shouldn't all be
        % populated, this determines the extent of the data by checking for the
        % index of the largest time and assuming that is the highest index.
        [~, idx] = max(abs_freq(:,1));
        qcmdata=abs_freq(1:idx,1:7); %Extracts data for times and harmonics of interest
        
        %There is at least one file which starts off with some time NaNs's,
        %which need to be removed before they cause trouble later. I don't care
        %about any rows that begin with NaN
        nantime = isnan(abs_freq(1:idx, 1));
        qcmdata = qcmdata(~nantime, 1:7);
        
        qcmsubdataf = qcmdata(:,2:2:7); %Splits for simplicity
        qcmsubdatag = qcmdata(:,3:2:7);
        
        qcmfits=chisq_values(nantime, 1:7); %removes NaNs from this too
        qcmsubfitsf = qcmfits(:,2:2:7);
        
        % Sets to NaN anything with a bad fit
        qcmsubdataf(qcmsubfitsf>1e-3) = NaN;
        qcmsubdatag(qcmsubfitsf>1e-3) = NaN;
        
        %Assigns "clean" data to its own file.
        qcmcleandata(:,[2 4 6]) = qcmsubdataf;
        qcmcleandata(:,[3 5 7]) = qcmsubdatag;
        
        %Checks for a _bare file to get reference data from, if not, uses
        %freq_shift_ref.
        reference = getreference(handles, qcmpath, filebase, 2);
        if length(reference) > 1
            for n=2:7
                offset.(name{n}) = reference(n);
            end
        else
            reference = freq_shift_ref;
            for n = 1:3
                offset.(name{n*2}) = reference(1,n);
                offset.(name{n*2+1}) = reference(2,n);
            end
        end
        
    elseif strcmp(filetypetxt,'.xls') || strcmp(filetype,'.xlsx')
        %Since I don't use excel files, this part is not updated, but will work
        % with old excel files, for instance, Garrett's data.
        handles.filetype = 3;
        tempdata=xlsread([qcmpath, qcmfile]);
        % set offsets to zero for  now - can change this if necessary
        clear qcmdata
        offset.f1=0;
        offset.f3=0;
        offset.f5=0;
        offset.g1=0;
        offset.g3=0;
        offset.g5=0;
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
    end
    
    %This function will calculate the error range from raw crystal data as long
    %as it is saved in a file of the same time with _bare appended to the base
    %filename. If there isn't one (or the file was an excel file), default
    %values will be used.
    filepathbase = [qcmpath filebase];
    [handles.din.bare.error.f handles.din.bare.error.g] = geterrorrange(hObject, handles, filepathbase);
    
    [nrows ncols]=size(qcmdata);
    if ncols<7 % account for possibility that there is not data from nh=5
        nhvals=[1,3];
        set(handles.n5, 'value', 0)
        handles.activenh = getactivenh(handles);
    else
        nhvals=[1,3,5];
    end
    
    %This part writes the reference values for the nhvals in use to the gui
    for nh=nhvals
        fref=['f',num2str(nh)];
        gref=['g',num2str(nh)];
        set(handles.(['offset', fref]),'string',commanumber(offset.(fref)))
        set(handles.(['offset', gref]),'string',commanumber(offset.(gref)))
    end
    
    %Extracts out the time data
    qcmt(1:nrows)=qcmdata(1:nrows,1);
    [qcmt,indexm,~]=unique(qcmt); %Not sure why this is here. Could possibly be removed?
    
    % Checks if the values are shifted or absolute. This check assumes the
    % frequency shift at 5MHz is less than 1MHz. If the first frequency value
    % is less than 4MHz, it assumes the data is already given in shifts and
    % removes the offset.
    if qcmdata(2,2)<4000000
        offset.f1=0;
        offset.f3=0;
        offset.f5=0;
        offset.g1=0;
        offset.g3=0;
        offset.g5=0;
    end
    
    %Shifts the data if necessary and saves it.
    for nh=nhvals
        fref=['f',num2str(nh)];
        gref=['g',num2str(nh)];
        delf(:,nh)=(qcmdata(indexm,nh+1)-offset.(fref));
        delg(:,nh)=(qcmdata(indexm,nh+2)-offset.(gref));
        cleandelf(:,nh)=(qcmcleandata(indexm,nh+1)-offset.(fref));
        cleandelg(:,nh)=(qcmcleandata(indexm,nh+2)-offset.(gref));
    end
    
    %Only applys threshhold filters if the sample isn't a bare crystal.
    if get(handles.isbare, 'value') == 0
        filter1 = abs(cleandelf(:,1)) > 1e5;
        filter3 = abs(cleandelf(:,3)./3) > 1e5;
        filter5 = abs(cleandelf(:,5)./5) > 1e5;
        % cleandelf([filter1 filter3 filter5], (1:2:5)) = NaN;
        % cleandelg([filter1 filter3 filter5], (1:2:5)) = NaN;
        cleandelf(filter1, 1) = NaN;
        cleandelf(filter3, 3) = NaN;
        cleandelf(filter5, 5) = NaN;
        cleandelg(filter1, 1) = NaN;
        cleandelg(filter3, 3) = NaN;
        cleandelg(filter5, 5) = NaN;
    end
end
%Saves all of the processed data into the handles structure.
handles.din.qcmt=qcmt;
handles.din.delf=delf;
handles.din.delg=delg;
handles.din.cleandelf = cleandelf;
handles.din.cleandelg = cleandelg;
handles.din.qcmpath=qcmpath;
handles.din.qcmfile=qcmfile;
handles.din.filebase=filebase;
handles.din.bare.offset = offset;

%Plots the imported data to the lower plots
plotraw(hObject, eventdata, handles);

%Imports the conductance data, if any
conductance = getconductance(hObject, handles);
handles.cond = conductance;

%Changes the data file to not be "loading" anymore.
set(handles.qcmfile,'string',qcmfile)

guidata(hObject, handles); %Saves new handles.

function [errorf errorg] = geterrorrange(hObject, handles, filepathbase);
% geterrorrange tries to open a '_bare' file of the same type as the main
% file and calculates the range of values for each quantity in it and saves
% these as the error range values.

% Base filename of potential bare file.
filename = [filepathbase '_bare'];

% Method for calculating for .fre files
if handles.filetype == 1 && exist([filename '.fre'], 'file') == 2
    barefile = [filename '.fre'];
    rawdata = importdata(barefile);
    ranges = range(rawdata.data(:,1:7));
    errorf = [ranges(2) 0 ranges(4) 0 ranges(6)];
    errorg = [ranges(3) 0 ranges(5) 0 ranges(7)];
% Method for .mat files
elseif handles.filetype == 2 && exist([filename '.mat'], 'file') == 2 || exist([filename '.fre'], 'file')
    try
        barefile = [filename '.mat'];
        load(barefile)
        [~, idx] = max(abs_freq(:,1));
        ranges = range(abs_freq(1:idx,1:7));
        errorf = [ranges(2) 0 ranges(4) 0 ranges(6)];
        errorg = [ranges(3) 0 ranges(5) 0 ranges(7)];
    catch Err
        barefile = [filename '.fre'];
        rawdata = importdata(barefile);
        ranges = range(rawdata.data(:,1:7));
        errorf = [ranges(2) 0 ranges(4) 0 ranges(6)];
        errorg = [ranges(3) 0 ranges(5) 0 ranges(7)];
    end
    % Method for excel files
elseif handles.filetype == 3
    disp('Default error values will be used since there isn''t a procedure for Excel files.')
    errorf = [44 0 120 0 273];
    errorg = [11 0 22 0 14];
% If file not found, use defaults
else
    disp('Default error values will be used since no bare data was found.')
    errorf = [44 0 120 0 273];
    errorg = [11 0 22 0 14];
end

function Solve_Callback(hObject,eventdata,handles)
% Uses the top selected harmonic combination to solve the currently loaded
% point.
try
    firstnhvals = handles.constants.nhvals(handles.mnumber);
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        handles.mnumber = getnhvals(handles);
        firstnhvals = handles.constants.nhvals(handles.mnumber);
    end
end
findsolution(hObject,eventdata,handles,1,firstnhvals{1});

function outputs = findsolution(hObject, eventdata, handles, time, nhset)
handles.activenh = getactivenh(handles);
%Reads the values of df and dg off the screen
for nh=handles.activenh.on
    fname=['f',num2str(nh),'exp'];
    gname=['g',num2str(nh),'exp'];
    df(nh)=str2num(get(handles.(fname),'string'));
    dg(nh)=str2num(get(handles.(gname),'string'));
end

%This if statement determines if one of the points used in the evaluation
%was not collected (turned into NaN). If there is a NaN, this turned that
%data point into a whole set of NaN's so it won't plot.
if sum(isnan(df(nhset(1:3)))) > 0 %That is, if one of the values is NaN
    outputs.df = nan(1,length(df));
    outputs.dg = nan(1,length(df));
    
    %NOTE: These may need to be set to the maximum nh value
    outputs.drho=nan(1,length(df));
    outputs.grho=nan(1,length(df));
    outputs.phi=NaN;
    outputs.dfcalc=nan(1,length(df));
    outputs.dgcalc=nan(1,length(df));
    outputs.d = nan(1,length(df));
    outputs.time=time;
    outputs.error = 1.5;
    set(handles.shownnhs, 'foregroundcolor', 'red')
    set(handles.shownnhs, 'string', 'not solved')
    disp(['Did not solve at time ' num2str(time) ' using ratio '...
            num2str(nhset(1)) ':' num2str(nhset(2)) ',' num2str(nhset(3))...
            ' because some data was NaN'])
else
    set(handles.shownnhs, 'foregroundcolor', 'black')
    f1=handles.constants.f1;
    zq=handles.constants.zq;
    % bulk phase properties hard coded for now
%     if get(handles.water,'value') && get(handles.drybase,'value')
%         rhobulk=1000;  % density of bulk phase
%         phibulk=90;  % bulk phase angle
%         phibulk=phibulk*pi/180; % convert to radians
%         etabulk=0.001; %  complex viscosity in radians
%         for n=handles.activenh.on
%             gstarbulk(n)=2*pi*f1*n*etabulk*exp(1i*phibulk);
%             zstarbulk(n)=(gstarbulk(n)*rhobulk)^0.5;
%             df(n)=df(n)-real(1i*f1*zstarbulk(n)/(pi*zq));
%             dg(n)=dg(n)-imag(1i*f1*zstarbulk(n)/(pi*zq));
%         end
%     end
    dissratio=dg(nhset(3))/df(nhset(3));
    harmratio=nhset(1)*df(nhset(2))/(nhset(2)*df(nhset(1)));
    % take initial guesses from simulation parameters
    x0(1) = str2num(get(handles.phi,'String'));
    dfield=['d',num2str(nhset(1)),'calc'];
    x0(2) = str2num(get(handles.(dfield),'String'));
    harmonicadjustment=str2num(get(handles.harmonicadjustment,'string'));
    % harmonic adjustment should always be 1, unless you are testing some
    % of the detailed behavior of the solutions
    
    ftosolve = @(x) fstar(x,nhset,harmratio,dissratio,harmonicadjustment);
    try
        options = optimset('Display','off');
        %[x,fval,exitflag] = fsolve(ftosolve,x0, options);
        [x,fval, exitflag,output, jacobian] = fsolve(ftosolve,x0,options); % Call solver
    catch MExc
        % Sometimes the fields aren't properly loaded. This tries again.
        if strcmp(MExc.identifier,  'optim:trustnleqn:UsrObjUndefAtX0')
            Sim_Callback(hObject, eventdata, handles)
            x0(1) = str2num(get(handles.phi,'String'));
            dfield=['d',num2str(nhset(1)),'calc'];
            x0(2) = str2num(get(handles.(dfield),'String'));
            options = optimset('Display','off');
            %[x,fval,exitflag] = fsolve(ftosolve,x0, options);
            try
            [x,fval, exitflag, output, jacobian] = fsolve(ftosolve,x0); %Try to call solver again.
            catch Err
                %This catch attempts to fix d/lambda values that are out of
                %whack because the mass has increased too far.
                if strcmp(MExc.identifier,  'optim:trustnleqn:UsrObjUndefAtX0')
                    if x0(2) == 0
                        x0(2) = .15;
                    end
                    try
                    [x,fval, exitflag, output, jacobian] = fsolve(ftosolve,x0);
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
    outputs.error = exitflag;
    
    %If the function did not solve satisfactorily (exitflag ~= 1), set data
    %to NaNs also.
    if exitflag ~= 1
        if length(df)==3
            outputs.df = [NaN NaN NaN];
            outputs.dg = [NaN NaN NaN];
        else
            outputs.df = [NaN NaN NaN NaN NaN];
            outputs.dg = [NaN NaN NaN NaN NaN];
        end
        outputs.drho=[NaN NaN NaN NaN NaN];
        outputs.grho=[NaN NaN NaN NaN NaN];
        outputs.phi=NaN;
        outputs.dfcalc=[NaN NaN NaN NaN NaN];
        outputs.dgcalc=[NaN NaN NaN NaN NaN];
        outputs.d = [NaN NaN NaN NaN NaN];
        outputs.time=time;
        
        set(handles.shownnhs, 'foregroundcolor', 'red')
        set(handles.shownnhs, 'string', 'not solved')
        disp(['Did not solve at time ' num2str(time) ' using ratio '...
            num2str(nhset(1)) ':' num2str(nhset(2)) ',' num2str(nhset(3))...
            ' with error ' num2str(exitflag)])
    % Also discard data if the phi value is impossible.
    elseif x(1) > 90
        if length(df)==3
            outputs.df = [NaN NaN NaN];
            outputs.dg = [NaN NaN NaN];
        else
            outputs.df = [NaN NaN NaN NaN NaN];
            outputs.dg = [NaN NaN NaN NaN NaN];
        end
        outputs.drho=[NaN NaN NaN NaN NaN];
        outputs.grho=[NaN NaN NaN NaN NaN];
        outputs.phi=NaN;
        outputs.dfcalc=[NaN NaN NaN NaN NaN];
        outputs.dgcalc=[NaN NaN NaN NaN NaN];
        outputs.d = [NaN NaN NaN NaN NaN];
        outputs.time=time;
        
        set(handles.shownnhs, 'foregroundcolor', 'red')
        set(handles.shownnhs, 'string', 'not solved')
        disp(['Did not solve at time ' num2str(time) ' using ratio '...
            num2str(nhset(1)) ':' num2str(nhset(2)) ',' num2str(nhset(3))...
            ' because phi out of range '])
        outputs.error = 2.5;
    else
        phi=x(1); %#ok<*SAGROW>
        dref=x(2);
        phir=phi*pi/180;  % phase angle in radians
        
        refnh = nhset(1);
        drho(1:2:5) = (df(refnh)/real(delfstar(dref,phi)))*zq/(2*refnh*f1^2); %calculate reference drho based on n2
        grho(refnh)=((drho(refnh)/dref)*refnh*f1*cos(phir/2))^2;
        for nh = 1:2:5
            d(nh)=dref*(nh/nhset(1))^(1-phi/180);
            grho(nh)=grho(nhset(1))*(nh/nhset(1))^(phi/90);
            dfcalc(nh)=sauerbrey(nh,drho(nh))*real(delfstar(d(nh),phi));
            dgcalc(nh)=sauerbrey(nh,drho(nh))*imag(delfstar(d(nh),phi));
            lambdarho(nh)=drho(nh)/(d(nh));
        end
%         for nh=[nhset(1),nhset(2)]
%             d(nh)=dref*(nh/nhset(1))^(1-phi/180);
%             drho(nh)=(df(nh)/real(delfstar(d(nh),phi)))*zq/(2*nh*f1^2);  %#ok<*AGROW>
%             grho(nh)=((drho(nh)/d(nh))*nh*f1*cos(phir/2))^2;
%             dfcalc(nh)=sauerbrey(nh,drho(nh))*real(delfstar(d(nh),phi));
%             dgcalc(nh)=sauerbrey(nh,drho(nh))*imag(delfstar(d(nh),phi));
%             lambdarho(nh)=drho(nh)/(d(nh));
%         end
%         % now we get the properties for the unused harmonic
%         for nh=nhset(4)
%             d(nh)=dref*(nh/nhset(1))^(1-phi/180);
%             drho(nh)=drho(nhset(1));  %#ok<*AGROW>
%             grho(nh)=grho(nhset(1))*(nh/nhset(1))^(phi/90);
%             dfcalc(nh)=sauerbrey(nh,drho(nh))*real(delfstar(d(nh),phi));
%             dgcalc(nh)=sauerbrey(nh,drho(nh))*imag(delfstar(d(nh),phi));
%             lambdarho(nh)=drho(nh)/(d(nh));
%         end        
        % adjust values to account for shift due to water immersion, if necessary
%         if get(handles.water,'value') && get(handles.drybase,'value')
%             for nh=handles.activenh.on
%                 dfcalc(nh)=dfcalc(nh)+real(1i*f1*zstarbulk(nh)/(pi*zq));
%                 dgcalc(nh)=dgcalc(nh)+imag(1i*f1*zstarbulk(nh)/(pi*zq));
%             end
%         end
        
        decaylength = sqrt(grho./(2.*pi.^2.*[5e6 NaN 15e6 NaN 25e6].^2)); %units of 
        maxnh = max(nhset(1:3)); %Get maximum harmonic in calculation
        drho(drho > decaylength(maxnh)*0.9) = NaN;
        
        if get(handles.calcerror, 'value') == 1 && length(unique(nhset(1:3))) < 3
            [drhoe grhoe phie] = finderror(hObject, handles, nhset);
            outputs.grhoe = grhoe;
            outputs.drhoe = drhoe;
            outputs.phie = phie;
        else
            set(handles.errordrho, 'string', '');
            set(handles.errorgrho, 'string', '');
            set(handles.errorphi, 'string', '');
        end
        
        % update the data table on the gui with the updated values
        drho=drho*1e3; % convert to g/m^2
        grho=grho*1e-3;  % convert to Pa-g/cm^3
        lambdarho=1e6*lambdarho;  % convert to microns
        set(handles.drho1,'string',num2str(drho(1),'%8.3f'))
        set(handles.grho1,'string',num2str(grho(1),'%8.2e'))
        set(handles.grho3,'string',num2str(grho(3),'%8.2e'))
        set(handles.grho5,'string',num2str(grho(5),'%8.2e'))
        set(handles.phicalc,'string',num2str(phi,'%8.2f'))
        set(handles.phi,'String',num2str(phi,'%8.2f'))
        set(handles.drho,'String',num2str(drho(1),'%8.2f'))
        set(handles.grho,'String',num2str(grho(1),'%6.2e'))
        set(handles.f1calc,'string',num2str(dfcalc(1),'%8.0f'))
        set(handles.f3calc,'string',num2str(dfcalc(3),'%8.0f'))
        set(handles.f5calc,'string',num2str(dfcalc(5),'%8.0f'))
        set(handles.f1fs,'string',num2str(real(delfstar(d(1),phi)),'%8.3f'))
        set(handles.f3fs,'string',num2str(real(delfstar(d(3),phi)),'%8.3f'))
        set(handles.f5fs,'string',num2str(real(delfstar(d(5),phi)),'%8.3f'))
        set(handles.g1calc,'string',num2str(dgcalc(1),'%8.0f'))
        set(handles.g3calc,'string',num2str(dgcalc(3),'%8.0f'))
        set(handles.g5calc,'string',num2str(dgcalc(5),'%8.0f'))
        set(handles.lambdarho1,'string',num2str(0.001*lambdarho(1),'%8.2f'))
        set(handles.lambdarho3,'string',num2str(0.001*lambdarho(3),'%8.2f'))
        set(handles.lambdarho5,'string',num2str(0.001*lambdarho(5),'%8.2f'))
        set(handles.d1calc,'string',num2str(d(1),'%8.4f'))
        set(handles.d3calc,'string',num2str(d(3),'%8.4f'))
        set(handles.d5calc,'string',num2str(d(5),'%8.4f'))
        
        set(handles.shownnhs, 'string', [num2str(nhset(1)) ':' num2str(nhset(2)) ',' num2str(nhset(3))])
        % add values to the output structure
        outputs.df=df;
        outputs.dg=dg;
        outputs.drho=drho;
        outputs.grho=grho;
        outputs.phi=phi;
        outputs.dfcalc=dfcalc;
        outputs.dgcalc=dgcalc;
        outputs.time=time;
        outputs.d = d;
        
        assignin('base','solveouptus',outputs)
    end
end

function Sim_Callback(hObject, eventdata, handles)
% Does a theoretical calculation for the values currently entered for mass,
% modulus, density, and phi.
f1=handles.constants.f1;
harmonicadjustment=str2num(get(handles.harmonicadjustment,'string'));
drho = 1e-3*str2num(get(handles.drho,'String'));
if get(handles.etabutton,'value')
    etarho(1)=1000*str2num(get(handles.etarho,'String'));
    grho(1)=eta1rho*2*pi*f1;
    set(handles.grho, 'string', num2str(grho(1),'%6.1e'))
    set(handles.grho, 'foregroundcolor', 'red')
    set(handles.etarho, 'foregroundcolor', 'blue')
else
    grho(1) = 1000*str2num(get(handles.grho,'String'));
    etarho(1)=grho(1)/(2*pi*f1);
    set(handles.etarho, 'string', num2str(etarho(1),'%8.3f'))
    set(handles.grho, 'foregroundcolor', 'blue')
    set(handles.etarho, 'foregroundcolor', 'red')
end
phi = str2num(get(handles.phi,'String'));
phir=phi*pi/180;
for n=[1,3,5]
    lambdarho(n)=n^(phi/180-1)*(1/5e6)*((grho(1))^0.5)/cos(phir/2);
    grho(n)=grho(1)*n^(phi/90);
    d(n)=drho/lambdarho(n);
    f(n)=sauerbrey(n,drho)*real(delfstar(d(n),phi));
    g(n)=sauerbrey(n,drho)*imag(delfstar(d(n),phi));
end
set(handles.f1calc,'string',num2str(f(1),'%8.0f'))
set(handles.f3calc,'string',num2str(f(3),'%8.0f'))
set(handles.f5calc,'string',num2str(f(5),'%8.0f'))
set(handles.f1fs,'string',num2str(real(delfstar(d(1),phi)),'%8.3f'))
set(handles.f3fs,'string',num2str(real(delfstar(d(3),phi)),'%8.3f'))
set(handles.f5fs,'string',num2str(real(delfstar(d(5),phi)),'%8.3f'))
set(handles.g1calc,'string',num2str(g(1),'%8.0f'))
set(handles.g3calc,'string',num2str(g(3),'%8.0f'))
set(handles.g5calc,'string',num2str(g(5),'%8.0f'))
set(handles.lambdarho1,'string',num2str(1e3*lambdarho(1),'%8.2f'))
set(handles.lambdarho3,'string',num2str(1e3*lambdarho(3),'%8.2f'))
set(handles.lambdarho5,'string',num2str(1e3*lambdarho(5),'%8.2f'))
set(handles.d1calc,'string',num2str(d(1),'%8.4f'))
set(handles.d3calc,'string',num2str(d(3),'%8.4f'))
set(handles.d5calc,'string',num2str(d(5),'%8.4f'))
set(handles.grho1,'string',num2str(0.001*grho(1),'%6.2e'))
set(handles.grho3,'string',num2str(0.001*grho(3),'%6.2e'))
set(handles.grho5,'string',num2str(0.001*grho(5),'%6.2e'))
set(handles.drho1,'string',num2str(drho*1000,'%6.3f'))

guidata(hObject, handles);

function F = fstar(x,nh,harmratio,dissratio,harmonicadjustment)
% Compares ideal and measured values for fstar for use in the solve
% function. The ideal output of the function is [0 0]. x(1) is the phase
% angle and x(2) is d/lambda.
phi=x(1);  % phase angle in degrees
d(1)=x(2);  % d/lambda at nh(1)
d(2)=d(1)*(nh(2)/nh(1))^(1-phi/180); % d/lambda at nh(2)
d(2)=d(2)*harmonicadjustment;
d(3)=d(1)*(nh(3)/nh(1))^(1-phi/180); % d/lambda at nh(3)
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
phir=phi*pi/180;
F=-(1/((2*pi*d)*(1-1i*tan(phir/2))))* ...
    tan(2*pi*d*(1-1i*tan(phir/2)));

function F=sauerbrey(n,drho)
% Calculates the sauerbry shift based on the harmonic and areal density.
F=2*n*5e6^2*drho/8.84e6;

function etabutton_Callback(hObject, eventdata, handles)
if get(handles.etabutton, 'value')
    set(handles.etarho, 'style', 'edit')
    set(handles.grho, 'style', 'text')
    set(handles.grho, 'foregroundcolor', 'red')
    set(handles.etarho, 'foregroundcolor', 'blue')
else
    set(handles.etarho, 'style', 'text')
    set(handles.grho, 'style', 'edit')
    set(handles.grho, 'foregroundcolor', 'blue')
    set(handles.etarho, 'foregroundcolor', 'red')
end

function conductance = getconductance(hObject, handles)
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
    samefiles = 1:length(files)-1;
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
            if rightstring && isempty(barestring)
                samefiles = [samefiles i];
            end
        catch Err
            %if ~strcmp(Err.identifier, 'MATLAB:nonLogicalConditional')
                rethrow Err
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
        if handles.filetype == 2 && ~strcmp(files{i}, 'reference') && ~strcmp(files{i}, 'version')
            parsed = textscan(files{i},'%s','delimiter','_');
            condfiles(index,1) = str2double(regexprep([parsed{1}{length(parsed{1})-4}], 'dot', '.'));
        elseif handles.filetype == 2 %if it is the "reference" or "version" file
            break
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

function fontsize_Callback(hObject, eventdata, handles)  % change font size of gui
switch get(handles.fontsize,'Value')
    case 1
        fontsize=10;
    case 2
        fontsize=11;
    case 3
        fontsize=12;
    case 4
        fontsize=13;
    case 5
        fontsize=14;
    case 6
        fontsize=15;
    case 7
        fontsize=16;
    case 8
        fontsize=17;
    case 9
        fontsize=18;
    case 10
        fontsize=19;
    case 11
        fontsize=20;
    otherwise
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
set(findall(gcf,'-property','FontSize'),'FontUnits','points')
%set([axes1, axes2], 'FontSize', fontsize)
guidata(hObject, handles)

function savebutton_Callback(hObject, eventdata, handles)
%saves the fig file
guidata(hObject,handles)
hgsave('onelayergui3.fig')

function reference = getreference(handles, qcmpath, filebase, filetype);
% Looks for reference data in bare file.
% If it is the bare file, set the reference values to 0. No shift.
if get(handles.isbare, 'value')
    reference = [0 0 0 0 0 0 0];
else
    try
        if filetype == 1 %fre file
            barefile = [qcmpath filebase '_bare.fre'];
            rawdata = importdata(barefile);
            reference = mean(rawdata.data(:,1:7));
        elseif filetype == 2 %mat file
            % There is the original bare file, but may also be a modified
            % one to correct for errors. Use the modified one if present.
            if exist([qcmpath filebase '_bare_data.mat']) == 2
                load([qcmpath filebase '_bare_data.mat']);
                freqmean = nanmean(delf(:,1:5));
                dissmean = nanmean(delg(:,1:5));
                reference = [ 0 freqmean(1) dissmean(1) freqmean(3) dissmean(3) freqmean(5) dissmean(5)];
            else
                barefile = [qcmpath filebase '_bare.mat'];
                load(barefile)
                [~, idx] = max(abs_freq(:,1));
                reference = nanmean(abs_freq(1:idx,1:7));
            end
        end
    catch Err
        %If it can't find a bare mat file, look for a bare fre file
        %instead.
        if filetype == 2 && strcmp(Err.identifier, 'MATLAB:load:couldNotReadFile')
            try
                barefile = [qcmpath filebase '_bare.fre'];
                rawdata = importdata(barefile);
                reference = mean(rawdata.data(:,1:7));
            catch
                reference = 0;
            end
            %If that doesn't work for whatever reason, set it to 0 to indicate failure
        else
            reference = 0;
        end
    end
end

function plotraw(hObject, eventdata, handles)
%Brings data into the function
try
    qcmt = handles.din.qcmt;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        disp('Unable to replot data. Please reload file')
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
colors{1}=[1,0,0];
colors{3}=[0,0.5,0];
colors{5}=[0,0,1];
legends{1}='n=1';
legends{3}='n=3';
legends{5}='n=5';
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

try
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

set(handles.axes2,'yscale','log')
%set(handles.axes2, 'yscale', 'lin')

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

if get(handles.log,'value')
    set([handles.axes1 handles.axes2],'xscale','log')
else
    set([handles.axes1 handles.axes2],'xscale','linear')
end

linkaxes([handles.axes1, handles.axes2], 'x')

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
handles.activenh = getactivenh(handles);
checksolveharms(hObject, handles);
disp('Left mouse button picks points.')
disp('Right mouse button picks last point.')
but = 1;
k=0; % keeps track of the number of points
try
struct2var(handles.din);
catch Err
    if strcmp('MATLAB:nonExistentField', Err.identifier)
       warning = warndlg('The data was not saved correctly. Please reload the file')
       uiwait(warning)
       changeqcmfile_Callback(hObject, eventdata, handles)     
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
    timefactor = handles.plotting.timefactor;
    for nh=handles.activenh.on
        nhstr=num2str(nh);
        df=interp1(qcmt,cleandelf(:,nh),time*timefactor);
        dg=interp1(qcmt,cleandelg(:,nh),time*timefactor);
        set(handles.(['f',nhstr,'exp']),'string',num2str(df,'%6.0f'))
        set(handles.(['g',nhstr,'exp']),'string',num2str(dg,'%6.0f'))
    end
    for nh=handles.activenh.off
        nhstr=num2str(nh);
        set(handles.(['f',nhstr,'exp']),'string','-')
        set(handles.(['g',nhstr,'exp']),'string','-')
    end
    if get(handles.autosolve,'value')
        for m = mnumber
            outputs{k,m} = findsolution(hObject, eventdata, handles, time, nhvals{m});
        end
        if get(handles.plotsolutioncontours, 'value') && outputs{k,m}.error == 1
            plotcontours(outputs{k,m}, m) %Make contour plot for last m value
            figure(onelayergui3) %Return figure focus to gui so cursor works properly
        end
    end
end
outputstructure = reformatdata(handles, outputs);
handles.dout = outputstructure;
guidata(hObject,handles);

function plotcontours(outputs, mnumber)
comb = {[1 3 1] [1 3 3] [1 5 1] [1 5 5] [3 5 3] [3 5 5] [1 3 5] [1 5 3] [3 5 1]};
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

dlofn=@(n,d1,phi) d1*n^(1-phi/180); %d/lambda
Dn=@(n,d1,phi)   2*pi*dlofn(n,d1,phi)*(1-1i*tand(phi/2));
dfstardn = @(Dn) -tan(Dn)/Dn; 
dfstar = @(n, d1, phi) dfstardn(Dn(n, d1, phi));
dlof1 = @(n, dn, phi) dn*n^(phi/180-1); %Takes higher order dl and turns it into d1.
frh = @(dn, n, phi) real(dfstar(n2, dlof1(n, dn, phi), phi))/real(dfstar(n1, dlof1(n, dn, phi), phi));
frd = @(dn, n, phi) imag(dfstar(n3, dlof1(n, dn, phi), phi))/real(dfstar(n3, dlof1(n, dn, phi), phi));

%%  Make the contour plots
    freqerror = [20 0 60 0 130];
    harmerror = [5 0 11 0 7];
    dissrtolerance = abs(((dg3-harmerror(n3))-(dg3+harmerror(n3)))/dg3);
    harmtolerance = abs(((df1-freqerror(n1))-(df1+freqerror(n1)))/df1);
   
    edissratio = delg(n3)./delf(n3);
    eharmratio = (n1/n2)*delf(n2)./delf(n1);
    
    dissrrange=[(1-dissrtolerance)*edissratio, (1+dissrtolerance)*edissratio]; % this is the range of dissipation ratios to consider (n=3)
    harmrrange=[(1-harmtolerance)*eharmratio, (1+harmtolerance)*eharmratio]; % range of harmonic ratios to consider (delf(3)/delf(1))

        n=75; % resolution of map
        dplot = linspace(0.0,0.23,n);  
        phiplot = linspace(0,90,n);
        for i = 1:n
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
        contourplotfig = figure('outerPosition',[width/2,height/2,width/2,height/2.3]);
         
        harmplotlim = [0, 3];
        harmplot = subplot(1,2,1);
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
        title('r_h')
        
        dissplotlim = [-2, 0];
        dissplot = subplot(1,2,2);
        contourf(dplot,phiplot,dissratio,linspace(dissplotlim(1), dissplotlim(2),256),'edgecolor','none');
        hold on
        contour(dplot,phiplot,dissratio,dissrrange,'edgecolor','black','linewidth',3,'linestyle','--')
        contour(dplot,phiplot,harmratio,harmrrange,'edgecolor','black','linewidth',3)
        colormap(jet(256));
        caxis(dissplotlim);
        colorbar
        xlabel(['d/\lambda_' num2str(n3)])
        ylabel('\phi')
        title('r_d')
        % plot(d3out(index),phiout(index),'k+','markersize',16,'color','red')
        % plot(harmplot,d3out(index),phiout(index),'k+','markersize',16,'color','red')
        
        plot(dissplot,d1out(n3), phiout, 'r+')
        plot(harmplot,d1out(n3), phiout, 'r+')
        
        linkaxes([dissplot, harmplot],'xy')



function responsemap = writecontourmap
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
responsemap = figure('outerPosition',[width/2,height/2,width/2,height/2.3]);

plot1 = subplot(1,2,1);
hold on
[C h1] = contourf(dl, phi, real(dfstar)',[logspace(1,1.24,250)-15]);
xlabel('d/\lambda_n')
ylabel('\phi_n')
colorbar
title('{\Delta}f_n/{\Delta}f_{sn}')
caxis([min_c, max_c])
set(h1,'LineStyle','none');
%set(findall(gcf,'-property','FontSize'),'fontsize',20)

plot2 = subplot(1,2,2);
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
linkaxes([plot1 plot2])
savefig('contourmap.fig')

function solveall_Callback(hObject, eventdata, handles)
checksolveharms(hObject, handles)
handles.activenh=getactivenh(handles);
mnumber = getnhvals(handles); %Determine harmonics to solve for
nhvals = handles.constants.nhvals; %Make the variable shorter to type

struct2var(handles.din);
[~,handles.datapoints]=size(qcmt);
for k=1:handles.datapoints;
    for nh=handles.activenh.on
        set(handles.(['f',num2str(nh),'exp']),'string',num2str(cleandelf(k,nh),'%6.0f'))
        set(handles.(['g',num2str(nh),'exp']),'string',num2str(cleandelg(k,nh),'%6.0f'))
    end
    for m = mnumber
        outputs{k,m}=findsolution(hObject, eventdata, handles,qcmt(k),nhvals{m});
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
    if mod(k,500)==0
        outputstructure = reformatdata(handles, outputs);
        handles.dout = outputstructure;
        handles = plotvalues(hObject,handles);
        disp([ num2str(k/handles.datapoints*100, 2) '% done'])
    end
end
outputstructure = reformatdata(handles, outputs)
handles.dout = outputstructure;
handles = plotvalues(hObject,handles);
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

function activenh = getactivenh(handles)
% Determines the harmonic buttons that are selected
activenh.on=[];
activenh.off=[];
for nh=[1,3,5]
    if get(handles.(['n',num2str(nh)]),'value')
        activenh.on=[activenh.on,nh];
    else
        activenh.off=[activenh.off,nh];
    end
end

function mnumber = getnhvals(handles)
% Determines which harmonics should be calculated

mnumber = [];
for m = 1:9
    % Reads in the checkmarks and records the numbers of the checked ones
    if get(handles.(['nhset',num2str(m),'select']),'value') == 1
        mnumber = [mnumber m];
    end
end

% Uses the first box if nothing is selected.
if isempty(mnumber)
    warning('Using default harmonic since none was selected')
    set(handles.nhset1select, 'value', 1)
    mnumber = 1;
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
if exist([handles.din.qcmpath handles.din.filebase '_cond.mat']) == 2
    cond = load([handles.din.qcmpath handles.din.filebase '_cond.mat'])
    if max(cond.time) == max(qcmt)
         disp('The conductance file appears to be up to date')
         return
    else
        maximptime = max(cond.time); %So only need to add new times
    end
else
    maximptime = 0; %Set to 0 if none found (so it will add all) 
    disp('The cond file will have to be built from scratch. This could take a while.')
end

pointstoadd = find(points(:,1)>maximptime); %Finds spectra from later times
for i = pointstoadd'
    if points(i,3) ~= 0
        % Import spectra and save to array
        harm = points(i,2)*2-1; %harmonic in 1,3,5
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
disp('New spectra added to file!')


function popoutplots_Callback(hObject, eventdata, handles)

handles.outputplot=outputplot;
handles.inputplot=inputplot;
guidata(hObject,handles);
plotvalues(hObject, handles);

function outstruct = reformatdata(handles, outputs)
doerror = get(handles.calcerror, 'value');
% each value of k corresponds to one time point selected from the plot
% variables with 'p' in the name are the ones we use to create the plots
try
    dfp=handles.din.cleandelf;
    dgp=handles.din.cleandelg;
catch Err
    if strcmp(Err.identifier, 'MATLAB:nonExistentField')
        disp('Unable to replot data. Please reload file')
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
            disp('No data to plot. Please solve for values before plotting.')
            return
        end
    else
        rethrow(Err)
    end
end

for k=1:todo(1)
    for m = mnumber
        struct2var(outputs{k,m})
        outstruct.grhop(k,m)=grho(1);
        outstruct.drhop(k,m)=drho(1);
        outstruct.phip(k,m)=phi(1);
        outstruct.timep(k,m)=time;
        if doerror
            try
                outstruct.grhoep(k,m) = grhoe(1);
                outstruct.drhoep(k,m) = drhoe(1);
                outstruct.phiep(k,m) = phie(1);
            catch Err
                if strcmp(Err.identifier, 'MATLAB:UndefinedFunction')
                    outstruct.grhoep(k,m) = NaN;
                    outstruct.drhoep(k,m) = NaN;
                    outstruct.phiep(k,m) = NaN;
                end
            end
        end
        for nh=handles.activenh.on
            outstruct.dfcalcp(k,m,nh)=dfcalc(nh);
            outstruct.dgcalcp(k,m,nh)=dgcalc(nh);
        end
    end
end

assignin('base','dfp',dfp); %What are these here for?
assignin('base','dgp',dgp);
assignin('base','grhop',outstruct.grhop);
assignin('base','drhop',outstruct.drhop);
assignin('base','phip',outstruct.phip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves data to _data file, including raw data.
if handles.overwrite == 1
    append = get(handles.qcmfile,'string');
    append = [handles.din.qcmpath append(1:end-4)];  
    if doerror
        data = struct('dgcalcp', outstruct.dgcalcp, 'dfcalcp', outstruct.dfcalcp, 'dgp', dgp, 'dfp', dfp,...
            'timep', outstruct.timep, 'grhop', outstruct.grhop, 'drhop', ...
            outstruct.drhop, 'phip', outstruct.phip, 'delg', handles.din.cleandelg,...
            'delf', handles.din.cleandelf, 'grhoep', outstruct.grhoep,...
            'phiep', outstruct.phiep, 'drhoep', outstruct.drhoep,...
            'time', handles.din.qcmt, 'error', handles.din.bare.error,...
            'offset', handles.din.bare.offset);
    else
        data = struct('dgcalcp', outstruct.dgcalcp, 'dfcalcp', outstruct.dfcalcp, 'dgp', dgp, 'dfp', dfp,...
            'timep', outstruct.timep, 'grhop', outstruct.grhop,...
            'drhop', outstruct.drhop, 'phip', outstruct.phip, 'delg',...
            handles.din.cleandelg, 'delf', handles.din.cleandelf,...
            'time', handles.din.qcmt, 'error', handles.din.bare.error,...
            'offset', handles.din.bare.offset);
    end
    save([append '_data.mat'], '-struct', 'data');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plotvalues(hObject, handles, inputfigure, outputfigure, inputaxes, outputaxes)
mnumber = getnhvals(handles);
handles.activenh = getactivenh(handles);

xlabeltext = handles.plotting.xlabeltext;
timefactor = handles.plotting.timefactor;

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
end

% determine the screensize so that the two figures are split vertically on
% the screen
scrsize = get(0,'ScreenSize');
left=scrsize(1);
bot=scrsize(2);
width=scrsize(3);
height=scrsize(4);

inputplot=figure('outerPosition',[width/2,height/2,width/2,height/2.3]);
setappdata(0, 'inputplot', inputplot);
inputaxes(1)=subplot(1,2,1);
inputaxes(2)=subplot(1,2,2);

outputplot=figure('outerPosition',[width/2,bot,width/2,height/2.3]);
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
symbols{7}='+';
symbols{8}='o';
symbols{9}='s';
colors{1}=[0, 0, 0];
colors{2}=[0, 0, 0];
colors{3}=[0, 0, 0];
colors{4}=[1, 0, 0];
colors{5}=[1, 0, 0];
colors{6}=[1, 0, 0];
colors{7}=[0, 0, 0];
colors{8}=[1, 0, 0];
colors{9}=[1, 0, 0];
linestyles = {'-' '' '-.' '' '--'};
datalinecolor{1}=[1,0,0];
datalinecolor{3}=[0,0.5,0];
datalinecolor{5}=[0,0,1];

figure(inputplot)
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
            disp('No data to plot for the 5th harmoic')
            handles.activenh.on = handles.activenh.on(handles.activenh.on~=5);
            index = index-1;
        end
    end
end

% now we plot the calculated values
for m = mnumber
    legendentries=[legendentries,index+1];
    for nh=handles.activenh.on
        index=index+1;
        plots(index)=plot(timep(:,m)./timefactor,dfcalcp(:,m,nh)/nh,symbols{m},'color','black');
        legendtext{index}=handles.constants.label{m};
    end
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

% now we plot the calculated values
for m = mnumber
    for nh=handles.activenh.on
        index=index+1;
        plots(index)=semilogy(timep(:,m)./timefactor,dgcalcp(:,m,nh),symbols{m},'color','black');
        legendtext{index}=handles.constants.label{m};
    end
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
axes(outputaxes(1))
hold on
if doerror
    for m = mnumber
        errorbar(timep(:,m)./timefactor,drhop(:,m),drhoep(:,m),symbols{m},'color',colors{m})    
    end
else
    for m = mnumber
        plot(timep(:,m)./timefactor,drhop(:,m),symbols{m},'color',colors{m})
    end
end
legendout(1)=legend(handles.constants.label{mnumber},'location','best');
xlabel(xlabeltext)
ylabel('d\rho (g/m^{2})')

axes(outputaxes(2))
hold on
if doerror
    for m = mnumber
        errorbar(timep(:,m)./timefactor,grhop(:,m),grhoep(:,m),symbols{m},'color',colors{m})
    end
else
    for m = mnumber
        plot(timep(:,m)./timefactor,grhop(:,m),symbols{m},'color',colors{m})
    end
end
legendout(2)=legend(handles.constants.label{mnumber}, 'location','best');
xlabel(xlabeltext)
ylabel('|G^{*}|\rho (Pa-g/cm^{3})')
if get(handles.logG, 'value')
    set(outputaxes(2), 'yscale', 'log');
else
    set(outputaxes(2), 'yscale', 'linear');
end

axes(outputaxes(3))
if doerror
    for m = mnumber
        errorbar(timep(:,m)./timefactor,phip(:,m),phiep(:,m),symbols{m},'color',colors{m})
        hold on
    end
else
    for m = mnumber
        plot(timep(:,m)./timefactor,phip(:,m),symbols{m},'color',colors{m})
        hold on
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

% we also need to expand all the legend boxes so that they look better in
% the printed vesion

% expandFact=1.5;  % this can be adjusted as appropriate
% for hL=[legendin(1:2) legendout(1:3)]
%     p = get(hL, 'position');
%     p(3) = p(3)*expandFact;
%     set(hL,'position', p)
%     % ht = findobj( get(hL,'children'), 'type', 'text');
%     set(gcf,'Resizefcn','')
% end

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
m.time = handles.din.qcmt;

function loadcleaned_Callback(hObject, eventdata, handles)
% Loads the "clean" variables from the _data file (delf and delg).
try
    load([handles.din.qcmpath handles.din.filebase '_data.mat'])
catch Err
    if strcmp(Err.identifier, 'MATLAB:load:couldNotReadFile')
        warndlg('There was no previously saved file to load');
        return
    elseif strcmp(Err.identifier, 'MATLAB:nonExistentField')
        disp('The previous file seems to have been lost. Please load a new one.')
        return
    else
        rethrow(Err)
    end
end
%If there is new data since the last save, it should be included. So this
%starts with all of the data as it is currently
handles.din.cleandelf = handles.din.delf;
handles.din.cleandelg = handles.din.delg;
%Gets the length of the saved data
curdata = length(time); %"time" is a variable from the saved data file
%Overwrites the times for which there is previously saved data with that
%data.
handles.din.cleandelf(1:curdata,:) = delf; 
handles.din.cleandelg(1:curdata,:) = delg;

guidata(hObject,handles); %Save changes to handles structure
plotraw(hObject, eventdata, handles) %Plot updated data

function opencond_Callback(hObject, eventdata, handles)
% This function opens the second program, condfig, which displays
% conductance data.
if isempty(handles.cond)
    warndlg('You cannot view the spectra since there aren''t any')
    return
end

% Updates (or creates) _cond file for use by condfig.
buildcondfile_Callback(hObject, eventdata, handles);

% To facilitate getting conductance data for a specific point,
% datacursormode is turned on, with a special function that
% sends the time to condfig.
datacursormode on
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,hObject, handles})

% This part actually runs condfig. When the program is closed, 
% any changes to the cleaned figure file is saved. Which means 
% that this function "stays open" until the window is closed.
handles.din = condfig2('onelayerguiwithcond', handles.figure);
guidata(hObject,handles);
datacursormode off
%The potentially changed data is plotted.
plotraw(hObject, eventdata, handles)
disp('The data has been transferred')

function output_txt = myupdatefcn(~,event_obj,~, handles)
% This is the function that runs when datacursormode is employed. The
% output output-txt is what appears in the box.

%Determines output box--this is the usual function of datacursor, modified
%to know what the x axis actually is.
pos = get(event_obj,'Position');
output_txt = {['Time: ',num2str(pos(1),5)],...
    ['y: ',num2str(pos(2),5)]};

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
        set(condhandles.currentpoint, 'string', num2str(pos(1), '%0.3f'))
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
combinations = {[1 3] [1 3] [1 5] [1 5] [3 5] [3 5]};
if length(nhsel) < 3
    for i = 1:6
        if get(handles.(['nhset' num2str(i) 'select']), 'value') == 1
            overlap = sum(ismember(nhsel, combinations{i}))
            if overlap ~= 2
                set(handles.(['nhset' num2str(i) 'select']), 'value', 0)
            end
        end
    end
end

function [drhoe grhoe phie] = finderror(hObject, handles, nh)
% Calculates the error in a measurement.
for nhs=handles.activenh.on
    fname=['f',num2str(nhs),'exp'];
    gname=['g',num2str(nhs),'exp'];
    dfi(nhs)=str2num(get(handles.(fname),'string'));
    dgi(nhs)=str2num(get(handles.(gname),'string'));
end

errorsf = handles.din.bare.error.f;
errorsg = handles.din.bare.error.g;

shift = 2;
shifts = -shift:shift:shift;

indecies = [unique(perms([1 2 2]), 'rows'); unique(perms([2 2 3]), 'rows')];
for i = 1:6
    k = indecies(i,1);
    l = indecies(i,2);
    m = indecies(i,3);
    df(nh(1)) = dfi(nh(1)) + shifts(k);
    df(nh(2)) = dfi(nh(2)) + shifts(l);
    dg(nh(3)) = dgi(nh(3)) + shifts(m);
    
    results = SolveQCM(df, dg, nh, hObject, handles);
    if results(4) == 1 && results(3) < 10^10 && results(1) > 0
        phi(k, l, m) = results(1);
        drho(k, l, m) = results(2);
        grho(k, l, m) = results(3);
    else
        phi(k, l, m) = NaN;
        drho(k, l, m) = NaN;
        grho(k, l, m) = NaN;
    end
end

[yphi xphi zphi] = gradient(phi, shift);
[ydrho xdrho zdrho] = gradient(drho, shift);
[ygrho xgrho zgrho] = gradient(grho, shift);

errphi = [xphi(2,2,2) yphi(2,2,2) zphi(2,2,2)];
errdrho = [xdrho(2,2,2) ydrho(2,2,2) zdrho(2,2,2)];
errgrho = [xgrho(2,2,2) ygrho(2,2,2) zgrho(2,2,2)];
error = [errphi; errdrho; errgrho];

shifts = [errorsf(nh(1)) errorsf(nh(2)) errorsg(nh(3))];

phie = sqrt(sum((error(1,:).*shifts).^2));
drhoe = sqrt(sum((error(2,:).*shifts).^2));
grhoe = sqrt(sum((error(3,:).*shifts).^2));

pphi = phie/str2num(get(handles.phi, 'string'))*100; 
pdrho = drhoe/str2num(get(handles.drho, 'string'))*100;
pgrho = grhoe/str2num(get(handles.grho, 'string'))*100;
set(handles.errordrho, 'string', num2str(pdrho));
set(handles.errorgrho, 'string', num2str(pgrho));
set(handles.errorphi, 'string', num2str(pphi));

function data = SolveQCM(dfi, dgi, nhi, hObject, handles)
%SolveQCM takes the frequency and dissipation inputs, along with the
%reference angles and guesses about the properties of the film, to
%calculate the actual properties of the film.

df = dfi;
dg = dgi;

f1=handles.constants.f1;
zq=handles.constants.zq;
dgn3dfn3=dg(nhi(3))/df(nhi(3));
dfn2dfn1=nhi(1)*df(nhi(2))/(nhi(2)*df(nhi(1)));
% take initial guesses from simulation parameters

%Tell the output not to display. I'll make separate errors based on the
%exitflags.
options = optimset('Display','off');

x0(1) = str2num(get(handles.phi,'String'));
dfield=['d',num2str(nhi(1)),'calc'];
x0(2) = str2num(get(handles.(dfield),'String'));
harmonicadjustment=str2num(get(handles.harmonicadjustment,'string'));
% if x0(1) > 90 || x0(2) > 1e10  || x0(2) == 0
%     answer = questdlg({['Failed to find solution at time ' num2str(qcmt(k)) ...
%         ' and ratio ' num2str(nhvals{m}(1)) ':' num2str(nhvals{m}(2))...
%         ',' num2str(nhvals{m}(3))],[],['Do you want to change the input values?']},...
%         'Yes', 'No')
%     switch answer
%         case 'Yes'
%             keystroke = 0;
%             while keystroke ~= 'w';
%                 waitforbuttonpress;
%                 keystroke = get(gcf, 'CurrentCharacter')
%             end
%             x0(1) = str2num(get(handles.phi,'String'));
%             dfield=['d',num2str(nhi(1)),'calc'];
%             x0(2) = str2num(get(handles.(dfield),'String'));
%         case 'No'
%         case 'Cancel'
%             return
%     end
% end

dissratio=dg(nhi(3))/df(nhi(3));
harmratio=nhi(1)*df(nhi(2))/(nhi(2)*df(nhi(1)));

if ~isnan(harmratio) && ~isnan(dissratio)
ftosolve = @(x) fstar(x,nhi,harmratio,dissratio,harmonicadjustment);
[x,fval,exitflag] = fsolve(ftosolve,x0,options); % Call solver
else
    exitflag = 0;
end
    
if exitflag ~= 1
    disp('Problem was not solved');
    x(1) = NaN;
    x(2) = NaN;
end
phi=x(1); 
dref=x(2);
phir=phi*pi/180;  % phase angle in radians

for nh=[nhi(1),nhi(2)]
    d(nh)=dref*(nh/nhi(1))^(1-phi/180);
    drho(nh)=(df(nh)/real(delfstar(d(nh),phi)))*zq/(2*nh*f1^2);  %#ok<*AGROW>
    grho(nh)=((drho(nh)/d(nh))*nh*f1*cos(phir/2))^2;
    dfcalc(nh)=sauerbrey(nh,drho(nh))*real(delfstar(d(nh),phi));
    dgcalc(nh)=sauerbrey(nh,drho(nh))*imag(delfstar(d(nh),phi));
    lambdarho(nh)=drho(nh)/(d(nh));
end

drho=drho*1e3; % convert to g/m^2
grho=grho*1e-3;  % convert to Pa-g/cm^3
data = [phi, drho(1), grho(1), exitflag];

function calcerror_Callback(hObject,eventdata,handles)
set(handles.errordrho, 'string', '');
set(handles.errorgrho, 'string', '');
set(handles.errorphi, 'string', '');

function maps_Callback(hObject, eventdata, handles)
try
    responsemap = openfig('contourmap.fig')
catch
    responsemap = writecontourmap
end
zoom(responsemap, 'reset')

dfp=handles.din.cleandelf;
dgp=handles.din.cleandelg;
timeinp=handles.din.qcmt;
todo = size(handles.dout);

for k=1:todo(1)
    %for m = mnumber
        struct2var(handles.dout{k}) %data for specific harmonic combination
        dp(k,:)=d(1:5);
        phip(k)=phi(1);
end

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

set(plot1, 'ylim', [max(floor(min(phip)/10)*10,0) min(ceil(max(phip)/10)*10,90)]);
set(plot1, 'xlim', [max(floor(min(min(dp(:,1:2:5)))/.3)*.3,0) min(ceil(max(max(dp(:,1:2:5)))/.3)*.3,1)]);

function autosolve_Callback(hObject, eventdata, handles)

function handles = xlabeltext_Callback(hObject, eventdata, handles)
%Updates the saved values for the axes and plotting.
switch get(handles.xlabeltext,'value')
    case 1
        handles.plotting.xlabeltext = 't (min.)';
        handles.plotting.timefactor = 1;
    case 2
        handles.plotting.xlabeltext = 't (hr.)';
        handles.plotting.timefactor = 60;
    case 3
        handles.plotting.xlabeltext = 't (day)';
        handles.plotting.timefactor = 1440;
    case 4
        handles.plotting.xlabeltext = 't (month)';
        handles.plotting.timefactor = 43830;
    case 5
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
    
function xlabeltext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nhset1_Callback(hObject, eventdata, handles)
function nhset1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nhset2_Callback(hObject, eventdata, handles)
function nhset2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nhset3_Callback(hObject, eventdata, handles)
function nhset3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nhset4_Callback(hObject, eventdata, handles)
function nhset4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nhset5_Callback(hObject, eventdata, handles)
function nhset5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nhset6_Callback(hObject, eventdata, handles)
function nhset6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nhset1select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset2select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset3select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset4select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset5select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset6select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset7select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset9select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function nhset8select_Callback(hObject, eventdata, handles)
handles.mnumber = getnhvals(handles);
guidata(hObject,handles);

function qcmfile_Callback(hObject, eventdata, handles)

function qcmfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function khz_Callback(hObject, eventdata, handles)

function fontsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function grho_Callback(hObject, eventdata, handles)

function grho_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function phi_Callback(hObject, eventdata, handles)

function phi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function drho_Callback(hObject, eventdata, handles)

function drho_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n1select_Callback(hObject, eventdata, handles)
function n3select_Callback(hObject, eventdata, handles)

function n5select_Callback(hObject, eventdata, handles)
function offsetf5_Callback(hObject, eventdata, handles)

function offsetf5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function offsetf3_Callback(hObject, eventdata, handles)

function offsetf3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function offsetf1_Callback(hObject, eventdata, handles)

function offsetf1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function offsetg1_Callback(hObject, eventdata, handles)

function offsetg1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function offsetg3_Callback(hObject, eventdata, handles)

function offsetg3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function offsetg5_Callback(hObject, eventdata, handles)

function offsetg5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function water_Callback(hObject, eventdata, handles)
function drybase_Callback(hObject, eventdata, handles)

function harmonicadjustment_Callback(hObject, eventdata, handles)
function harmonicadjustment_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function g5exp_Callback(hObject, eventdata, handles)
function g5exp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function g3exp_Callback(hObject, eventdata, handles)

function g3exp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function g1exp_Callback(hObject, eventdata, handles)

function g1exp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function f5exp_Callback(hObject, eventdata, handles)

function f5exp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function f3exp_Callback(hObject, eventdata, handles)

function f3exp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function f1exp_Callback(hObject, eventdata, handles)
function f1exp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function accesshandles_Callback(hObject, eventdata, handles)
keyboard;

function isbare_Callback(hObject, eventdata, handles)

function loadsavedplot_Callback(hObject, eventdata, handles)
handles.saveplot = load([handles.din.qcmpath handles.din.filebase '.mat']);
handles.saveplot.plot = 1;
plotvalues(hObject, handles);
handles.saveplot.plot = 0;
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

function plotsolutioncontours_Callback(hObject, eventdata, handles)

% --- Executes on button press in hideplots.
function hideplots_Callback(hObject, eventdata, handles)
% hObject    handle to hideplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.calcdataplotspanel, 'visible', 'off')