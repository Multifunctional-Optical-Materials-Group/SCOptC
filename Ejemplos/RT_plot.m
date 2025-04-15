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


%% ReTraF calculation test

%% Files
data_file   =[];

%% Parameters
wl =350:1:900;  % Wavelength (nm)

theta.values = [0];  % Incident angles of the measurements (degrees)
theta.index  = [1:length(theta.values)];  % Index of the angles to use


%% Layer models of the structutre

%% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index

%% Al2O3
Al2O3.type = "cnst";  
Al2O3.n = 1.64;
Al2O3.D = 300;

%% PVK
pvk.type = "Fh-N";
pvk.n0 = 2.2696;
pvk.Eg = 1.5205;
pvk.Ei = [1.5890 2.5300 3.0990];
pvk.Gi = [0.1012 0.6205 0.5663];
pvk.fi = [0.1347 0.0960 0.3251];
% pvk.type = "file";  
% pvk.filename = "pvk.mat";   
pvk.D    = 200;

%% TaTM
TaTM.type = "file";  
TaTM.filename = "TaTm_nk.mat";   
TaTM.D = 15; 

%% ITO
ITO.type = "file";
ITO.filename ="ito_harsh_nk.mat";
ITO.D = 50;

%% C60
C60.type = "file";  
C60.filename = "C60_nk.mat";   
C60.D = 5;    % Thickness (nm)

%% SnO2
SnO2.type = "file";  
SnO2.filename = "SnO2_nk.mat";   
SnO2.D = 35;    % Thickness (nm)

%% Au
Au.type = "file";  
Au.filename = "Au_nk.mat";   
Au.D = 100;    % Thickness (nm)


%% Glass substrate
glass.type = "cnst";    % Known constant refractive index type
glass.n = 1.51;      % Refractive index
glass.D = 1e6;      % Layer thickness (nm)

%%
layers = {air, glass, air}; % Prueba material transparente


%% Options
foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)
foptions.scatt = false;

%% Calculation

% data_file = [];

[layers_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,layers,data_file,foptions);

%%

n1=air.n;
n2=Al2O3.n;
Rfres=((n1-n2)/(n1+n2))^2

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
