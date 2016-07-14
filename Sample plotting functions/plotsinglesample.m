function allaxes = plotsinglesample(file, savefilename, timescale, ranges, axisscale, gaxisscale, ploterror, numofpts, harm, harmonics, symbols, colors, legendlocation, legendframe)
timerange = ranges{1};
drhorange = ranges{2};
grhorange = ranges{3};
phirange = ranges{4};
delfrange = ranges{5};
delgrange = ranges{6};

try
    load([file 'data.mat']);
catch
    load([file '_data.mat']);
end

figuredefaults

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
inplimits = {delfrange delgrange};

% each value of k corresponds to one time point selected from the plot
% variables with 'p' in the name are the ones we use to create the plots

%times = linspace(1, 25671, 24);
if ~strcmp(numofpts, 'all')
    switch axisscale
        case 'lin'
            times = linspace(timerange(1), timerange(2), numofpts);
        case 'log'
            times = logspace(log10(timerange(1)+.01), log10(timerange(2)), numofpts);
    end
    
    pointstoplot = [];
    for i = times*timecorr
        [~, index] = min(abs(timep - i));
        while isnan(grhop(index(1),harm))
            index = index + 1;
            if index > length(drhop)
                break
            end
        end
        if index <= length(timep)
            pointstoplot = [pointstoplot index];
        end
    end
    
    %pointstoplot = pointstoplot <= length(timep);
else
    pointstoplot = 1:length(timep);
end

% pointstoplotraw = [];
% for i = times*timecorr
%     [~, index] = min(abs(time - i));
%     while isnan(delp(index,1)) 
%         index = index + 1;
%     end
%     pointstoplotraw = [pointstoplotraw index];
% end


labels = { '1:3,1'    '1:3,3'    '1:5,1'    '1:5,5'    '3:5,3'    '3:5,5'};
% determine the screensize so that the two figures are split vertically on
% the screen
scrsize = get(0,'ScreenSize');
left=scrsize(1);
bot=scrsize(2);
width=scrsize(3);
height=scrsize(4);

%  first plot the property figure

outputplot=figure('outerPosition',[width/2,bot,width/1.5,height/2.3]);
if isempty(symbols)
    symbols{1}='o';
    symbols{2}='+';
    symbols{4}='s';
    symbols{3}='x';
    symbols{5}='^';
    symbols{6}='d';
    colors{1}=[1, 0, 0];
    colors{2}=[0, 0.5, 0];
    colors{4}=[0, 0, 1];
    colors{3}=[0, 0, 0];
    colors{5}=[1, 0, 0];
    colors{6}=[0, 0.5, 0];
end

outputaxes1=subplot(1,3,1);
datasize = size(grhop);
%Just do the first one.
whichtoplot = harm;
for m=whichtoplot
    if ploterror
        errorbar(timep(pointstoplot,m)/timecorr,drhop(pointstoplot,m),drhoep(pointstoplot,m),symbols{m},'color',colors{m})
    else
        plot(timep(pointstoplot,m)/timecorr,drhop(pointstoplot,m),symbols{m},'color',colors{m})
    end
    hold on
end
if legendframe == 1
    leg = legend(labels{whichtoplot},'location',legendlocation{2});
end
xlabel(xlabeltext);
ylabel('d\rho (g/m^{2})');
ylim(limits{1});

outputaxes2=subplot(1,3,2);
for m=whichtoplot
    if ploterror
        errorbar(timep(pointstoplot,m)/timecorr,grhop(pointstoplot,m),grhoep(pointstoplot,m),symbols{m},'color',colors{m})
    else
        plot(timep(pointstoplot,m)/timecorr,grhop(pointstoplot,m),symbols{m},'color',colors{m})
    end
    hold on
end
if legendframe == 2
    leg = legend(labels{whichtoplot},'location',legendlocation{2});
end
xlabel(xlabeltext)
try
ylabel(['|G^{*}_' num2str(refG) '|\rho (Pa-g/cm^{3})'])
catch
    ylabel(['|G^{*}_1|\rho (Pa-g/cm^{3})'])
end
ylim(limits{2});


outputaxes3=subplot(1,3,3);
for m=whichtoplot
    if ploterror
    errorbar(timep(pointstoplot,m)/timecorr,phip(pointstoplot,m),phiep(pointstoplot,m),symbols{m},'color',colors{m})
    else
        plot(timep(pointstoplot,m)/timecorr,phip(pointstoplot,m),symbols{m},'color',colors{m})
    end
    hold on
end
if legendframe == 3
    leg = legend(labels{whichtoplot},'location',legendlocation{2});
end
%leg = legend(labels{whichtoplot},'location','best')
xlabel(xlabeltext)
ylabel('\phi (deg.)')
ylim(limits{3});

% now we compare simulated and actual frequency and dissipation shifts
set(0,'defaultlinemarkersize',4)
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
linestyles = {'o' '' 'x' '' '+'};
datalinecolor{1}=[1,0,0];
datalinecolor{3}=[0,0.5,0];
datalinecolor{5}=[0,0,1];
inputplot=figure('outerPosition',[width/2,height/2,width/2,height/2.3]);
inputaxes1=subplot(1,2,1);
% first plot the measured frequency shifts for the three harmonics
index=0;
legendentries=[];

try
for nh=harmonics
    index=index+1;
    legendentries=[legendentries,index];
    plots(index)=plot(time/timecorr,delf(:,nh)/nh/1000,linestyles{nh},'color',datalinecolor{nh});
    hold on
    legendtext{index}=['n=',num2str(nh)];
end
set(inputaxes1, 'ylim', inplimits{1})


% now we plot the calculated values
% for m=1:datasize(2)
%     legendentries=[legendentries,index+1];
%     for nh=[1 3 5]
%         index=index+1;
%         plots(index)=plot(timep(pointstoplot,m)/60,dfcalcp(pointstoplot,m,nh)/nh,symbols{m},'color',datalinecolor{nh});
%         legendtext{index}=labels{m};
%     end
% end
xlabel(xlabeltext);
ylabel('\Deltaf_{n}/n (kHz)');
leg = legend(plots(legendentries),legendtext(legendentries),'location',legendlocation{1});
%leg = legend(plots(legendentries),legendtext(legendentries),'location','best');
catch
end
inputaxes2=subplot(1,2,2);
% repeat the plot for the dissipation
% first plot the measured data for the three harmonics
index=0;
try
for nh=harmonics
    index=index+1;
    plots(index)=plot(time/timecorr,delg(:,nh),linestyles{nh},'color',datalinecolor{nh});
    legendtext{index}=['n=',num2str(nh)];
    
    hold on
end
catch
end

xlabel(xlabeltext);
ylabel('\Delta\Gamma_{n} (Hz)');

%leg = legend(plots(legendentries),legendtext(legendentries),'location',legendlocation{1});
%set(inputaxes2, 'ytick', logspace(0,6,7));
set(inputaxes2, 'ylim', inplimits{2})

symbols{1}='o';
symbols{2}='+';
symbols{4}='s';
symbols{3}='x';
symbols{5}='^';
symbols{6}='d';
colors{1}=[1, 0, 0];
colors{2}=[0, 0.5, 0];
colors{4}=[0, 0, 1];
colors{3}=[0, 0, 0];
colors{5}=[1, 0, 0];
colors{6}=[0, 0.5, 0];

set(outputplot, 'PaperUnits', 'inches', 'PaperSize', [14 4])
set(outputplot, 'PaperPositionMode', 'manual')
set(outputplot, 'PaperPosition', [0 0 14 4])
% set(outputplot, 'PaperUnits', 'inches');
% set(outputplot, 'PaperSize', [12 3.2]);
% set(outputplot, 'PaperPositionMode', 'manual');
% set(outputplot, 'PaperPosition', [0 0 12 3.2]);

set(inputplot, 'PaperUnits', 'inches');
set(inputplot, 'PaperSize', [8 3.2]);
set(inputplot, 'PaperPositionMode', 'manual');
set(inputplot, 'PaperPosition', [0 0 8 3.2]);

allaxes = [inputaxes1, inputaxes2, outputaxes1, outputaxes2, outputaxes3];
%set(allaxes, 'Xtick', 0:150:600);
%set(allaxes, 'Xlim', [0 480]);
%set([inputplot, outputplot],'PaperPositionMode','auto') 
set(findall(gcf,'-property','FontSize'),'fontsize',20)
set(leg,'FontSize',14);
set(allaxes, 'Xlim', timerange);
set(allaxes, 'xscale', axisscale);
set(outputaxes2, 'yscale', gaxisscale);

if strcmp(axisscale, 'log')
    set(allaxes, 'xscale', 'log')
    if length(get(allaxes(1), 'xticklabel')) < 4
        set(allaxes, 'xtick', [10.^-2 10.^-1 10.^0 10.^1 10.^2 10.^3 10.^4 10.^5 10.^6 10.^7 10.^8 10.^9 10.^10]);
    else
        set(allaxes, 'xtick', [10.^-1 10.^1 10.^3 10.^5 10.^7 10.^9]);
    end
end
    
set(inputaxes1,'xlim',get(outputaxes1,'xlim'));
set(inputaxes2,'xlim',get(outputaxes1,'xlim'));
set(leg, 'Location', legendlocation{1});
set(leg,'FontSize',14);

print(outputplot,'-depsc2',['../Figures/' savefilename '_calc.eps']);
print(outputplot,'-dpng',['../Figures/' savefilename '_calc.png']);
print(inputplot, '-depsc2', ['../Figures/' savefilename '_dfdg.eps']);
print(inputplot, '-dpng', ['../Figures/' savefilename '_dfdg.png']);

%%
% figure
% semilogy(phip(pointstoplot,2), grhop(pointstoplot,2));
% xlabel('\phi (deg.)')
% ylabel('|G^*|\rho')
% print(gcf, '-depsc2', [savefilename '_gvsphi.eps'])