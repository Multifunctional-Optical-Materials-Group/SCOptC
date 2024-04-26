clear; clc;



%% FAMAPbI3 cell simulation


%% Parameters
wl =400:5:800;  % Wavelength (nm)

theta.values = 0;  % Incident angles of the measurements (degrees)


%% Layer models of the structutre


% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index

% LiF
LiF_top.type = "file";  
LiF_top.filename = "LiF_nk.mat";   
LiF_top.D = 64;    % Thickness (nm)
LiF_top.active = "false";


% Al2O3
Al2O3.type = "file";  
Al2O3.filename = "Al2O3_ald_nk.mat";   
Al2O3.D = 24;    % Thickness (nm)
Al2O3.active = "false";

% ITO thin
ITO_top.type = "file";
ITO_top.filename ="ito_harsh_nk.mat";
ITO_top.D = 10.5;
ITO_top.active = "false";

% C60
C60.type = "file";  
C60.filename = "C60_nk.mat";   
C60.D = 12;    % Thickness (nm)
C60.active = "false";


% SnO2
SnO2.type = "file";  
SnO2.filename = "SnO2_nk.mat";   
SnO2.D = 35;    % Thickness (nm)
SnO2.active = "false";



% PVK
pvk.active = "true";

pvk.type = "Fh-N";
pvk.n0 = 2.2696;
pvk.Eg = 1.5205;
pvk.Ei = [1.5890 2.5300 3.0990];
pvk.Gi = [0.1012 0.6205 0.5663];
pvk.fi = [0.1347 0.0960 0.3251];

pvk.D = 500;


% TaTM
TaTM.type = "file";  
TaTM.filename = "TaTm_nk.mat";   
TaTM.D = 6;    % Thickness (nm)
TaTM.active = "false";



% ITO thick
ITO_bot.type = "file";
ITO_bot.filename ="ito_harsh_nk.mat";
ITO_bot.D = 40;
ITO_bot.active = "false";

% Glass substrate
glass.type = "cnst";    % Known constant refractive index type
glass.n = 1.51;      % Refractive index
glass.D = 1e6;      % Layer thickness (nm)
glass.active = "false";


% LiF
LiF_bot.type = "file";  
LiF_bot.filename = "LiF_nk.mat";   
LiF_bot.D = 80;    % Thickness (nm)
LiF_bot.active = "false";



models = {air , LiF_top, Al2O3, ITO_top, C60, SnO2, pvk, TaTM, ITO_bot, glass, LiF_bot ,air};




%% Fitting options (no fitting will take place as there isn't any Unknown type material)

foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)
foptions.backwards = "false";
foptions.zstep = 0.5;
foptions.plot = "true";
foptions.SQgap = 1.55;



%% Results

[models_out,N,D,results,foptions_out] = SCOptC(wl,theta,models,foptions);



%% Plot E2

h = cumsum(D);
E2 = squeeze(results.E2);
z =results.z(:);

figure;
pcolor(wl,z,E2(:,:)'); shading flat; colormap hot; colorbar
hold on

for k=1:length(h)-2
    line([wl(1),wl(end)],[h(k),h(k)],'LineWidth',2, 'LineStyle', '--', 'Color', 'white')
end

ylim([0,800])

%



