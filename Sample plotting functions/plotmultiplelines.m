function [allaxes newfig] = plotmultiplelines(fileroot, labels, savefilename, timescale, ranges, axisscale, gaxisscale, ploterror, numofpts, harm, symbols, colors, legendlocation, legendframe, normalize)
%Variable arguments should be symbols, then colors, of an array the same
%size or larger than the number of files to be plotted.

timerange = ranges{1};
drhorange = ranges{2};
grhorange = ranges{3};
phirange = ranges{4};
%% Read data
numfiles = length(fileroot);

% This reads in each of the files, saving the relevant data to the "data"
% structure, where data(1) is the first file, data(2) the second, etc. 
for i = 1:numfiles
    %   simplefile{i} = [fileroot{i} '.mat'];
    
    datafile{i} = [fileroot{i} '_data.mat'];
    
    try
        data(i) = load(datafile{i});
    catch Err
        %This happens if the files were saved in different orders, and the
        %order of the structure names is different.
        if strcmp(Err.identifier, 'MATLAB:load:couldNotReadFile')
            datafile{i} = [fileroot{i} 'data.mat'];
            
            try
                data(i) = load(datafile{i});
            catch Err
                if strcmp(Err.identifier, 'MATLAB:heterogeneousStrucAssignment')
                    data2 = load(datafile{i});
                    fields = {'dgcalcp', 'dfcalcp', 'dgp', 'dfp', 'timep', 'grhop',...
                        'drhop', 'phip', 'delg', 'delf', 'grhoep', 'phiep', 'drhoep', 'time'};
                    for k = 1:length(fields)
                        data(i).(fields{k}) = data2.(fields{k});
                    end
                    
                else
                    rethrow(Err)
                end
            end
            
        elseif strcmp(Err.identifier, 'MATLAB:heterogeneousStrucAssignment')
            data2 = load(datafile{i});
            fields = {'dgcalcp', 'dfcalcp', 'dgp', 'dfp', 'timep', 'grhop',...
                'drhop', 'phip', 'delg', 'delf', 'grhoep', 'phiep', 'drhoep', 'time'};
            for k = 1:length(fields)
                data(i).(fields{k}) = data2.(fields{k});
            end
            
        else
            rethrow Err
        end
    end
    if ~isfield(data(i), 'grhoep') && ploterror == 1 
        warndlg(['The saved file ' fileroot{i} ' does not contain error data and cannot be processed.'])
    end
end

%% Set limits, find data points
set(0,'defaultaxesfontsize',16)
set(0,'defaultlinemarkersize',12)
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesbox','on')

switch timescale
    case 'min'
        xlabeltext = 'Time (min.)';
        timecorr = 1;
    case 'hr'
        xlabeltext = 'Time (hr.)';
        timecorr = 60;
    case 'day'
        xlabeltext = 'Time (days)';
        timecorr = 1440;
end

limits = {drhorange grhorange phirange};

switch axisscale
    case 'lin'
        times = linspace(timerange(1), timerange(2), numofpts);
    case 'log'
        times = logspace(log10(timerange(1)+.001), log10(timerange(2)), numofpts);       
end


for j = 1:numfiles
    pointstoplot{j} = [];
    for i = times*timecorr
        [~, index] = min(abs(data(j).timep(:,harm(1)) - i));
        modulus = data(j).grhop;
        try
            while isnan(modulus(index,harm(1))) && index<length(data(j).timep(:,harm(1)))
                index = index + 1;
            end
            pointstoplot{j} = [pointstoplot{j} index];
        catch Err
            keyboard
        end       
    end
end

%% Plot lines
figuredefaults

scrsize = get(0,'ScreenSize');
left=scrsize(1);
bot=scrsize(2);
width=scrsize(3);
height=scrsize(4);

%  first plot the property figure
newfig=figure('outerPosition',[width/2,bot,width/1.5,height/2.3]);

set(newfig, 'PaperUnits', 'inches', 'PaperSize', [14 4])
set(newfig, 'PaperPositionMode', 'manual')
set(newfig, 'PaperPosition', [0 0 14 4])

outputaxes1=subplot(1,3,1);
hold on

m = harm;

for j = 1:numfiles
    if strcmp(normalize, 'min')
        [minval minidx] = findpeaks(-data(j).drhop(pointstoplot{j},m));
        if isempty(minval)
            norm(j) = data(j).drhop(pointstoplot{j}(1),m); %Normalize to beginning
        else
            norm(j) = -minval(1);
        end
        disp(['Sample ' labels{j} ' normalized to ' num2str(norm)])
        drhoylabel = 'Normalized mass';
    elseif strcmp(normalize, 'max')
        [minval minidx] = findpeaks(data(j).drhop(pointstoplot{j},m));
        norm(j) = minval(1);
        drhoylabel = 'Normalized mass';
    elseif strcmp(normalize, 'beg')
        norm(j) = data(j).drhop(pointstoplot{j}(1),m);
        drhoylabel = 'Normalized mass';
    else
        drhoylabel = 'd\rho (g/m^2)';
        norm(j) = 1;
    end
    
    if ploterror == 1
        h = errorbar(data(j).timep(pointstoplot{j},m)/timecorr,data(j).drhop(pointstoplot{j},m)./norm(j),data(j).drhoep(pointstoplot{j},m)./norm(j),symbols{j},'color',colors{j});
        h.MarkerSize = markersize;
        h.LineWidth = linewidth;
    else
        plot(data(j).timep(pointstoplot{j},m)/timecorr,data(j).drhop(pointstoplot{j},m)./norm(j),symbols{j},'color',colors{j});
    end
end

if legendframe == 1
    leg = legend(labels,'location','best');
    set(leg, 'Interpreter', 'none');
end

xlabel(xlabeltext);
ylabel(drhoylabel);
ylim(limits{1});

outputaxes2=subplot(1,3,2);
hold on

for j = 1:numfiles
    if ploterror == 1
        h = errorbar(data(j).timep(pointstoplot{j},m)/timecorr,data(j).grhop(pointstoplot{j},m),data(j).grhoep(pointstoplot{j},m),symbols{j},'color',colors{j});
        h.MarkerSize = markersize;
        h.LineWidth = linewidth;
    else
        plot(data(j).timep(pointstoplot{j},m)/timecorr,data(j).grhop(pointstoplot{j},m),symbols{j},'color',colors{j});
    end
end

if legendframe == 2
    leg = legend(labels,'location','best');
    set(leg, 'Interpreter', 'none');
end

xlabel(xlabeltext)
try
    Gref = [];
    for i = 1:numfiles
        Gref = [data(i).refG Gref];
    end
    if length(unique(Gref)) > 1
        warndlg('The samples have different reference G values')
    else
        ylabel(['|G^{*}_' num2str(Gref(1)) '|\rho (Pa-g/cm^{3})'])
    end
catch
    ylabel(['|G^{*}|\rho (Pa-g/cm^{3})'])
end
ylim(limits{2});

outputaxes3=subplot(1,3,3);
hold on

for j = 1:numfiles
    if ploterror == 1
        h = errorbar(data(j).timep(pointstoplot{j},m)/timecorr,data(j).phip(pointstoplot{j},m),data(j).phiep(pointstoplot{j},m),symbols{j},'color',colors{j});
        h.MarkerSize = markersize;
        h.LineWidth = linewidth;
    else
        plot(data(j).timep(pointstoplot{j},m)/timecorr,data(j).phip(pointstoplot{j},m),symbols{j},'color',colors{j});
    end
end

if legendframe == 3
    leg = legend(labels,'location','best');
    set(leg, 'Interpreter', 'none');
end

xlabel(xlabeltext)
ylabel('\phi (deg.)')
ylim(limits{3});

% phi/grho plot
% h = figure;
% hold on
% 
% for j = 1:numfiles
%     plot(data(j).phip(pointstoplot{j},m), data(j).grhop(pointstoplot{j},m),symbols{j},'color',colors{j})
% end
% phigrho = gca;
% 
% xlabel('\phi (deg.)');
% set(phigrho, 'yscale', 'log')
% ylabel('|G^{*}|\rho (Pa-g/cm^{3})')
% 
% leg = legend(labels, 'location','southwest');
% set(leg, 'Interpreter', 'none');
% set(findall(gcf,'-property','FontSize'),'fontsize',20)
% set(leg, 'FontSize', 14, 'location', 'southeast')
% 
% print(h,'-depsc', [savefilename '_phigrho.eps'])
% print(h,'-dpng', [savefilename '_phigrho.png'])

%% Finalize axis and save
allaxes = [outputaxes1, outputaxes2, outputaxes3];

set(allaxes, 'xscale', axisscale)

set(outputaxes2, 'yscale', gaxisscale)

%set(allaxes, 'Xtick', 0:6:36);
set(allaxes, 'Xlim', timerange);

set(leg, 'FontSize', 14)
set(leg, 'Location', legendlocation);
hold off
if length(get(allaxes(1), 'xticklabel')) < 4
    if strcmp(axisscale, 'log')
        set(allaxes, 'xtick', [10.^-2 10.^-1 10.^0 10.^1 10.^2 10.^3 10.^4 10.^5 10.^6 10.^7 10.^8 10.^9 10.^10]);
    else
        %set(allaxes, 'xtick', [timerange(1):round((timerange(2)-timerange(1))/20)*5:timerange(2)]);
    end
elseif strcmp(axisscale, 'log')
    set(allaxes, 'xtick', [10.^-1 10.^1 10.^3 10.^5 10.^7 10.^9]);
end

print(newfig,'-depsc', [savefilename '.eps'])
%print(newfig,'-dpng', [savefilename '.png'])

