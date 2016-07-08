dl = linspace(0,.25,250);
phi = linspace(0,90,250);

delfstar = @(d,phi) -(1/((2*pi*d)*(1-1i*tand(phi/2))))* ...
    tan(2*pi*d*(1-1i*tand(phi/2)));

for i = 1:length(dl)
    for j = 1:length(phi)
        dfstar(i,j) = delfstar(dl(i), phi(j));
    end
end

min_c = -3;
max_c = 0;
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

