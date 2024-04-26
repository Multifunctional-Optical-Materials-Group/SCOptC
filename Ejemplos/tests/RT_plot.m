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

%% ReTraF calculation test

%% Files
data_file   =[];

%% Parameters
wl =350:10:900;  % Wavelength (nm)

theta.values = [0];  % Incident angles of the measurements (degrees)
theta.index  = [1:length(theta.values)];  % Index of the angles to use


%% Layer models of the structutre

%% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index

%% Al2O3
Al2O3.type = "cnst";  
Al2O3.n = 1.64;
Al2O3.D = 30;

%% PVK
pvk.type = "Fh-N";
pvk.Eg   = 1.4937;
pvk.n0   = 2.0777;
pvk.fi   = [0.1337 0.0810 0.0551];
pvk.Ei   = [1.5041 2.4662 3.3123];
pvk.Gi   = [0.0905 0.4218 0.3844];
pvk.D    = 91.2763;

%% TaTM
TaTM.type = "file";  
TaTM.filename = "TaTm_nk.mat";   
TaTM.D = 5; 


%% Glass substrate
glass.type = "cnst";    % Known constant refractive index type
glass.n = 1.51;      % Refractive index
glass.D = 1e6;      % Layer thickness (nm)

%%
% layers = {air, Al2O3, pvk, TaTM, glass, air};
%layers = {air, Al2O3, Al2O3};
% layers = fliplr(layers);

%% Options
foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)
foptions.scatt = false;

%% Calculation

%data_file = [];

[layers_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,layers,data_file,foptions);

%%

% n1=air.n;
% n2=Al2O3.n;
% Rfres=((n1-n2)/(n1+n2))^2

figure;
subplot(1,2,1)
plot(wl,real(N(:,1)),'LineWidth',2)
xlabel("\lambda (nm)")
ylabel("n")
subplot(1,2,2)
plot(wl,imag(N(:,1)),'LineWidth',2)
xlabel("\lambda (nm)")
ylabel("k")
%