%% Inputs
clear all
%Ideally this should be all that you need to change to get the
%output that you want.

%The filenames should include the full file path except the extension.
%For these files, assume they are in a subfolder of "Current_Projects", and
%put the folders in after this. 
file = 'LFS_LO_01';
savefilename = 'LFS_LO_01';

timescale = 'day';  %options are 'min' 'hr' or 'day'
timerange = [0 15];  %should be in units of timescale;

axisscale = 'log';  %x-axis scale
gaxisscale = 'log'; %y-axis of Grho plot
ploterror = 1;
numofpts = 50;      %points to be plotted in the time range
harm = [1 2]; %numeric value from 1-6, in order 1:3,1, 1:3,3, 1:5,1, 1:5,5, 3:5,3, 3:5,5
harmonics = [1 3 5]; %harmonics for which you want data plotted ([1 3 5])

drhorange = [1 5];
grhorange = [1e6 1e9];
phirange = [0 90];

delfrange = [-50 0];
delgrange = [-Inf Inf];

symbols = [];
colors = [];

legendlocation = {'best','southeast'}; %raw data, calculated data
legendframe = 2;


%% Do not change below this line
ranges = {timerange drhorange grhorange phirange delfrange delgrange};
allaxes = plotsinglesample(file, savefilename, timescale, ranges, axisscale, gaxisscale, ploterror, numofpts, harm, harmonics, symbols, colors, legendlocation, legendframe);
