clear; clc; close all



%% FAMAPbI3 cell simulation

addpath("C:\Users\glozano\Nextcloud\G\teaching\máster UPO\24-25\energía\práctica sims\códigos\SCOptC-1.1\src");

%% Parameters
wl =350:2:900;  % Wavelength (nm)

theta.values = 0;  % Incident angles of the measurements (degrees)


%% Layer models of the structutre


% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index

% LiF
LiF_top.type = "file";  
LiF_top.filename = "LiF_nk.mat";   
LiF_top.D = 0;    % Thickness (nm)
LiF_top.active = "false";


% Al2O3
Al2O3.type = "file";  
Al2O3.filename = "Al2O3_ald_nk.mat";   
Al2O3.D = 0;    % Thickness (nm)
Al2O3.active = "false";

% FTO
FTO.type = "file";
FTO.filename ="FTO_nk.mat";
FTO.D = 0;
FTO.active = "false";

% TiO2_bulk
TiO2_bulk.type = "file";
TiO2_bulk.filename ="TiO2_bulk_nk.mat";
TiO2_bulk.D = 0;
TiO2_bulk.active = "false";

% TiO2_meso
TiO2_meso.type = "file";
TiO2_meso.filename ="TiO2_meso_nk.mat";
TiO2_meso.D = 0;
TiO2_meso.active = "active";

% ITO thin
ITO_top.type = "file";
ITO_top.filename ="ito_harsh_nk.mat";
ITO_top.D = 0;
ITO_top.active = "false";

% C60
C60.type = "file";  
C60.filename = "C60_nk.mat";   
C60.D = 0;    % Thickness (nm)
C60.active = "false";


% SnO2
SnO2.type = "file";  
SnO2.filename = "SnO2_nk.mat";   
SnO2.D = 0;    % Thickness (nm)
SnO2.active = "false";



% PVK
pvk.active = "true";

pvk.type = "file";
pvk.filename = "pvk_nk.mat"; 

pvk.D = 0;


% TaTM
TaTM.type = "file";  
TaTM.filename = "TaTm_nk.mat";   
TaTM.D = 0;    % Thickness (nm)
TaTM.active = "false";



% ITO thick
ITO_bot.type = "file";
ITO_bot.filename ="ito_harsh_nk.mat";
ITO_bot.D = 0;
ITO_bot.active = "false";

% Glass substrate
glass.type = "cnst";    % Known constant refractive index type
glass.n = 1.51;      % Refractive index
glass.D = 1e6;      % Layer thickness (nm)
glass.active = "false";


% LiF
LiF_bot.type = "file";  
LiF_bot.filename = "LiF_nk.mat";   
LiF_bot.D = 0;    % Thickness (nm)
LiF_bot.active = "false";

% spiro
spiro.type = "file";  
spiro.filename = "spiro_nk.mat";   
spiro.D = 0;    % Thickness (nm)
spiro.active = "false";

% Au
Au.type = "file";  
Au.filename = "Au_nk.mat";   
Au.D = 0;    % Thickness (nm)
Au.active = "false";


%% Layers

models = {air , glass , XXX ,air};




%% Fitting options (no fitting will take place as there isn't any Unknown type material)

foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)
foptions.backwards = "false";
foptions.zstep = 0.5;
foptions.plot = "true";
foptions.SQgap = 1.55;



%% Results

[models_out,N,D,results,foptions_out] = SCOptC(wl,theta,models,foptions);



%% Plot E2

h = cumsum(D)-glass.D;
E2 = squeeze(results.E2);
Pabs = squeeze(results.Pabs);
z =results.z(:)-glass.D;

figure;
pcolor(wl,z,E2(:,:)'); shading flat; colormap hot; colorbar
hold on

for k=1:length(h)-1
    line([wl(1),wl(end)],[h(k),h(k)],'LineWidth',2, 'LineStyle', '--', 'Color', 'white')
end

ylim([-200,1200])
xlabel("\lambda (nm)")
ylabel("z (nm)")
title("|E|^2/|E_0|^2")

figure;
pcolor(wl,z,Pabs(:,:)'); shading flat; colormap hot; colorbar
hold on

for k=1:length(h)-1
    line([wl(1),wl(end)],[h(k),h(k)],'LineWidth',2, 'LineStyle', '--', 'Color', 'white')
end

ylim([-200,1200])
%set(gca,"ColorScale","log")
xlabel("\lambda (nm)")
ylabel("z (nm)")
title("Absorption Density")

%



