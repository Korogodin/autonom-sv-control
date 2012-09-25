try 
    close(handle_fig_main); % Close old output form
end

close all
clear all
clc

globals;
addpath([pwd '/data/']);

mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
omega_e = 0.7292115e-4; % [rad/s] Earth's rotation rate
        
Tmod = 15*60*60;
dTmod = 300;
t = 0:dTmod:Tmod;
Nmod = length(t);

GPS_Const = CreateGPS();