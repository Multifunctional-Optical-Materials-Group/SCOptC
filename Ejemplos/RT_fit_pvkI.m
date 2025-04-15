%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; 
% close all;

addpath("C:\Users\glozano\Nextcloud\G\teaching\máster UPO\24-25\energía\práctica sims\códigos\ReTraF-1.4\src");


%% ReTraF fitting test

plot_exp_data = 1; % Plot experimental data (1 if you want the plot)


%% Files

%  Size of delta & psi data must be (  length(wl_exp) ,  length(theta_exp) )
%
data_file = "pvk.mat";


%% Parameters
wl =350:10:900;  % Wavelength (nm)

theta.values = [6 30 50];  % Incident angles of the measurements (degrees)
theta.index  = [1,3];  % Index of the angles to use


%% Layer models of the structutre


%% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index

%% Al2O3
Al2O3.type = "U-cnst";  
Al2O3.l_n = 1.60;    
Al2O3.u_n = 1.85;
Al2O3.l_D = 25;    % Thickness (nm)
Al2O3.u_D = 65;

%% PVK
pvk.type = "U-Fh-N";
pvk.l_n0 = 0.5*1.99;
pvk.u_n0 = 1.5*1.99;

pvk.l_Eg = 1239.8/830;
pvk.u_Eg = 1239.8/830;

pvk.l_Ei = 0.5*[1.4939 2.8737 4.0264 ];
pvk.u_Ei = 1.5*[1.4939 2.8737 4.0264 ];

pvk.l_Gi = 0.5*[0.080 0.387 0.448  ];
pvk.u_Gi = 1.5*[0.080 0.387 0.448  ];
pvk.l_fi = 0.5*[0.149 0.078 0.056  ];
pvk.u_fi = 1.5*[0.149 0.078 0.056  ];

pvk.l_D = 40;
pvk.u_D = 100;
%pvk.D = 50; % 85nm sale bien


%% TaTM
TaTM.type = "file";  
TaTM.filename = "TaTm_nk.mat";   
TaTM.D = 5;    % Thickness (5 nm)



%% Glass substrate
glass.type = "cnst";    % Known constant refractive index type
glass.n = 1.51;      % Refractive index
glass.D = 1e6;      % Layer thickness (nm)

%%
models = {air, Al2O3, pvk, TaTM, glass ,air};

% models = fliplr(models);

%% Plot data
if plot_exp_data == 1
    data = load("pvk.mat");

    wl_exp = data.wl_exp;
    Rexp = 0.5*(data.RSample_S+data.RSample_P);
    Texp = 0.5*(data.TSample_S+data.TSample_P);

    figure;
    subplot(1,2,1)
    plot(wl_exp,Texp)
    xlabel("\lambda (nm)")
    ylabel("T")

    subplot(1,2,2)
    plot(wl_exp,Rexp)
    xlabel("\lambda (nm)")
    ylabel("R")
    return
end


%% Fitting options

foptions.method = "fmincon";  % "fmincon" or "genetic"
% foptions.method = "genetic";  % "fmincon" or "genetic"
foptions.scatt = false;
foptions.itermax = 50;      % Maximum number of iterations
foptions.parallel = true;   % Parallel evaluation of the objective function
foptions.popsize = 300;      % Population size (if "genetic" is used)
foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)



%% Fit

%data_file = [];

[layers_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions);

%%

figure;
subplot(1,2,1)
plot(wl,real(N(:,3)),'LineWidth',2)
xlabel("\lambda (nm)")
ylabel("n")
subplot(1,2,2)
plot(wl,imag(N(:,3)),'LineWidth',2)
xlabel("\lambda (nm)")
ylabel("k")
%