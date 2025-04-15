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
data_file = "transparent.mat";


%% Parameters
wl =400:10:800;  % Wavelength (nm)

theta.values = [6 30 50];  % Incident angles of the measurements (degrees)
theta.index  = [1,2,3];  % Index of the angles to use


%% Layer models of the structutre


%% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index



%% Material (nanophosphor)
Material.type = "U-cnst";  
Material.l_n = 1.2;    
Material.u_n = 1.6;
Material.l_D = 20;    % Thickness (nm)
Material.u_D = 100;



%% Glass substrate
glass.type = "cnst";    % Known constant refractive index type
glass.n = 1.46;      % Refractive index
glass.D = 1e6;      % Layer thickness (nm)

%%

models = {air, Material, glass ,air};



%% Fitting options

foptions.method = "fmincon";  % "fmincon" or "genetic"
%foptions.method = "genetic";  % "fmincon" or "genetic"
foptions.scatt = true;
foptions.itermax = 50;      % Maximum number of iterations
foptions.parallel = true;   % Parallel evaluation of the objective function
foptions.popsize = 300;      % Population size (if "genetic" is used)
foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)



%% Plot data
if plot_exp_data == 1
    data = load("transparent.mat");

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

%% Fit

%data_file = [];

[layers_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions);

%%

figure;
subplot(1,2,1)
plot(wl,real(N(:,2)),'LineWidth',2)
xlabel("\lambda (nm)")
ylabel("n")
subplot(1,2,2)
plot(wl,imag(N(:,2)),'LineWidth',2)
xlabel("\lambda (nm)")
ylabel("k")
%