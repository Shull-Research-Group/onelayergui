function optimumthickness(hObject, eventdata, handles)
%This function creates a figure which takes inputs to calculate the maximum
%thickness that a film can be before the dissipation becomes too high.
close(figure(996)); figure(996); %close window if open
pos=get(figure(996),'position');
set(figure(996),'position',[pos(1)/1 pos(2)/1 pos(3)/1.7 pos(4)/1.5]); %Set window size
set(figure(996),'menubar','none','toolbar','none','numbertitle','off','name','Optimum thickness calculator');
laxis = axes('parent',figure(996),'units','normalized','position',[0 0 1 1],'visible','off'); %Creates an axis for the text boxes
% The text boxes are necessary rather than uicontrols because they contain
% greek letters

c1ledge = 0.15; %location of label edges
c2ledge = 0.6; %location of box edges
c1w = .4; %Width of the left column

%Create label and input box for the maximum dissipation
max_dg = uicontrol('style','edit','string','10','units','normalized',...
    'position',[c2ledge 0.85 .25 .1]);
text(c1ledge, 0.9, 'max \Delta\Gamma (kHz)');

%Create label and input box for the maximum thickness
max_drho = uicontrol('style','edit','string','10','units','normalized',...
    'position',[c2ledge 0.75 .25 .1]);
text(c1ledge, 0.8, 'max d\rho (g/m^2)');

%Create label and input box for the phase angle It defaults to the current
%phase angle in onelayergui
phi = uicontrol('style','edit','string','45','units','normalized',...
    'position',[c2ledge 0.65 .25 .1]);
text(c1ledge, 0.7, '\phi');

%Create label and pull-down for the modulus reference harmonic It defaults
%to the value in onelayergui
ref_G = uicontrol('style','popup','string',{'<html>G<sub>1</sub></html>',...
    '<html>G<sub>3</sub></html>', '<html>G<sub>5</sub></html>',...
    '<html>G<sub>7</sub></html>'},'value', 2,...
    'units','normalized', 'position',[c2ledge 0.55 .25 .1]);
text(c1ledge, 0.6, 'reference G');

%Create label and input box for the modulus It defaults to the current
%modulus in onelayergui
grho = uicontrol('style','edit','string','1e8','units','normalized',...
    'position',[c2ledge 0.45 .25 .1]);
text(c1ledge, 0.5, '|G*|p (Pa-g/cm^3)');

%Create check boxes for the different harmonics
nh1 = uicontrol('style','checkbox','string','5MHz','units','normalized',...
    'position',[.1 0.28 .25 .1]);
nh3 = uicontrol('style','checkbox','string','15MHz','units','normalized',...
    'value', 1, 'position',[.1+.19 0.28 .25 .1]);
nh5 = uicontrol('style','checkbox','string','25MHz','units','normalized',...
    'value', 1, 'position',[.1+.42 0.28 .25 .1]);
nh7 = uicontrol('style','checkbox','string','35MHz','units','normalized',...
    'position',[.1+.65 0.28 .25 .1]);

calcmax = text(0.50, 0.15, 'Max d\rho: '); %Text box for the answer

%Initializes button and callback for the Calculate button
uicontrol('style','push','string','Calculate','units','normalized',...
    'position',[c1ledge 0.10 .3 .1],'callback',{@maxthicknesscalc_Callback,...
    max_dg, max_drho, phi, ref_G, grho, nh1, nh3, nh5, nh7, calcmax});

set(findall(figure(996),'-property','FontSize'),'FontSize',12)

function maxthicknesscalc_Callback(hObject, ~, max_dg, max_drho, phi, ref_G, grho, nh1, nh3, nh5, nh7, calcmax)
%This function calculates the maximum thickness that a film can be before
%the dissipation of a particular harmonic reaches a critical
%threshhold--usually 10,000Hz, but this can be changed. This function does
%not take into account decay length limitations on the thickness for mass
%sensitivity and does not explore possibilities beyond the dissipation
%spike at d/lambda = 0.25.

% Get which harmonics the user is concerned about
nh = [1 3 5 7];
nh = nh(logical([get(nh1, 'value') get(nh3, 'value') get(nh5, 'value') get(nh7, 'value')]));
if isempty(nh)
    warndlg('You forgot to select harmonics to look at')
    return
end
grhoi = str2num(get(grho, 'string'))*1000; %Get modulus in the right units
phii = str2num(get(phi, 'string')); %Get phi
refn = (get(ref_G, 'value'))*2-1; %Get refn
maxdrho = str2num(get(max_drho, 'string'))*0.001; %Get maxdrho
maxdg = str2num(get(max_dg, 'string'))*1000; %Get max dissipation

% Returns if any values are either empty or otherwise not numbers
if isempty(grhoi) || isempty(phii)
    warndlg('Please give a valid value for the physical parameters')
    return
end

for n = nh
    lambdarhon = lambdarhof(refn, n, grhoi, phii); %Lambda rho for the harmonic of interest
    drho = 0.001*0.1; % Initialize for while loop starting at .1g/m^2
    delfstarc = 0;
    % Loops over increasing thicknesses to find when the dissipation
    % becomes to great
    while imag(delfstarc) < maxdg && drho < maxdrho
        drho = drho + 0.000025; % This sets the increment -- resolution is to the 0.025g/m^2
        dl = drho./lambdarhon; %Updates dl based on the new drho
        delfstarc = sauerbrey(n,drho)*delfstar(dl, phii); %Calculates the frequency and dissipation     
    end
    d(n) = (drho-0.000025)*1000;
end
%Displays the maximum thickness
set(calcmax, 'string', ['Max d\rho: ' num2str(min(d(d~=0))) 'g/m^2']) 

function lrho = lambdarhof(refn, n, grho, phi)
if refn == n
    lrho = 1/(n*5e6)*(grho^0.5)/cosd(phi/2);
else
    lrho = 1/(n*5e6)*((grho*(n^(phi/90))/(refn^(phi/90)))^0.5)/cosd(phi/2);
end

function F=delfstar(d,phi)  % input phase angle is in degrees
% Calculates delfstar (rhs of delf/delfsn equation) with input of d/lambda
% and phi
F=-(1/((2*pi*d)*(1-1i*tand(phi/2))))* ...
    tan(2*pi*d*(1-1i*tand(phi/2)));

function F=sauerbrey(n,drho)
% Calculates the sauerbry shift based on the harmonic and areal density.
F=2*n*5e6^2*drho/8.84e6;