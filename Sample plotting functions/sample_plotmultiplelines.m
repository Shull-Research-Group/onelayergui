%% Inputs
%Ideally this should be all that you need to change to get the
%output that you want.

%The filenames should include the full file path except the extension.
%
n = 1;
fileroot{n} = 'LFS_LO_01';
labels{n} = '1';
n = n+1;
fileroot{n} = 'LFS_LO_02';
labels{n} = '2';
n = n+1;
fileroot{n} = 'LFS_LO_03';
labels{n} = '3';

savefilename = 'linseedoil'; %Filename (and location) to save the files

timescale = 'day';  %options are 'min' 'hr' or 'day'
timerange = [0 100];  %should be in units of timescale;

axisscale = 'log';  %x-axis scale
gaxisscale = 'log'; %y-axis of Grho plot
ploterror = 1;      %1 for yes, 0 for no
numofpts = 50;      %number of points to be plotted in the time range
harm = 1;           %numeric value from 1-6, in order 1:3,1, 1:3,3, 1:5,1, 1:5,5, 3:5,3, 3:5,5

drhorange = [1 6];
grhorange = [1e6 1e9];
phirange = [0 100];

legendlocation = 'best';
legendframe = 3;
colors = [];
symbols = [];
normalize = 'none';

%% Do not change below this line

ranges = {timerange, drhorange, grhorange, phirange};

[allaxes newfig] = plotmultiplelines(fileroot, labels, savefilename, timescale, ranges, axisscale, gaxisscale, ploterror, numofpts, harm, symbols, colors, legendlocation, legendframe, normalize);