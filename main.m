try 
    close(handle_fig_main); % Close old output form
end

close all
clear 
clc

globals;
addpath([pwd '/data/']);

mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
omega_e = 0.7292115e-4; % [rad/s] Earth's rotation rate
R_e = 6371e3;  % [m] Mean Earth radiius
        
Tmod = 15*60*60;
dTmod = 300;
t = CTime(0:dTmod:Tmod, 24, 09, 2012);
Nmod = length(t.t);

% global Omega1 Omega2 Omega3 cOmega1 cOmega2 cOmega3
% Omega1 = 0; Omega2 = 0; Omega3 = 0;
GLO_Const = LoadGLOConst([pwd '/data/Topcon_120925.agl'], t);
figure(1)
hold on
for j = 1:length(GLO_Const)
    plot3(GLO_Const(j).x, GLO_Const(j).y, GLO_Const(j).z, 'm')
    plot3(GLO_Const(j).x(1), GLO_Const(j).y(1), GLO_Const(j).z(1), '*g', 'MarkerSize', 10)
%     fprintf('%f\n', GLO_Const(j).Alm.Omega);
end
hold off


GPS_Const = LoadGPSConst([pwd '/data/Topcon_120925.agp'], t);
figure(1)
hold on
for j = 1:length(GPS_Const)
    plot3(GPS_Const(j).x, GPS_Const(j).y, GPS_Const(j).z, 'y')
    plot3(GPS_Const(j).x(1), GPS_Const(j).y(1), GPS_Const(j).z(1), '*b', 'MarkerSize', 10)
%     fprintf('%f\n', GPS_Const(j).Alm.Omega);
end
hold off

N_GPS = length(GPS_Const);
N_GLO = length(GLO_Const);

aa = 0;
h_min = 0e3;
h_max = 5.0e6;
dh = 100e3;
Na = length((R_e+h_min):dh:(R_e + h_max));

min_S = nan(1, Na);
max_S = nan(1, Na);
proc_S = nan(1, Na);

min_S_GLO = nan(1, Na);
max_S_GLO = nan(1, Na);
proc_S_GLO = nan(1, Na);

min_S_GPS = nan(1, Na);
max_S_GPS = nan(1, Na);
proc_S_GPS = nan(1, Na);

h_S = nan(1, Na);

for a = (R_e+h_min):dh:(R_e + h_max)
    aa = aa+1;

    fprintf('Calculation for LEO satt, h = %f\n', a-R_e);
    % Космос-2480 — российский разведывательный спутник типа Кобальт-М.
    % http://ru.wikipedia.org/wiki/Космос-2480
    % R_e + 259.5e3 = (1 + e)*a;
    % R_e + 196.4e3 = (1 - e)*a; => a = 6598950; e = 0.0048;
    e = 0.0048;
    % a = 6598950; 
    % a = a+750e3;
    Omega = deg2rad(10);
    omega = 0; M0 = 0; 
    i = deg2rad(81.4);
    LEO_SV =  SpaceVehicle_LEO('LEO', a, e, Omega, omega, i, M0, t);
    LEO_SV.calcKeplerOrbit(t);
    fprintf('OK\n');
    figure(1)
    hold on
    plot3(LEO_SV.x, LEO_SV.y, LEO_SV.z, 'r')
    plot3(LEO_SV.x(1), LEO_SV.y(1), LEO_SV.z(1), '*r', 'MarkerSize', 10)
    hold off

    Sucs_Counter_GLO = zeros(1, Nmod);
    Sucs_Counter_GPS = zeros(1, Nmod);
    for k = 1:Nmod

        LEO_SV.X = [LEO_SV.x(k); LEO_SV.y(k); LEO_SV.z(k)];
        LEO_SV.AntM = LEO_SV.X/LEO_SV.r(k);
        
        for j = 1:N_GPS
            GPS_Const(j).X = [GPS_Const(j).x(k); GPS_Const(j).y(k); GPS_Const(j).z(k)];
            GPS_Const(j).AntM = -GPS_Const(j).X /GPS_Const(j).r(k);
            SVSV = (GPS_Const(j).X - LEO_SV.X); SVSV = SVSV / norm(SVSV);
            AntAngleLEO = rad2deg(acos(LEO_SV.AntM' * SVSV));
            AntAngleGPS = rad2deg(acos(GPS_Const(j).AntM' * (-SVSV)));
            if (abs(AntAngleLEO) <= 85) && (abs(AntAngleGPS) <= 25)
                Sucs_Counter_GPS(k) = Sucs_Counter_GPS(k) + 1;
            end
        end
        
        for j = 1:N_GLO
            GLO_Const(j).X = [GLO_Const(j).x(k); GLO_Const(j).y(k); GLO_Const(j).z(k)];
            GLO_Const(j).AntM = -GLO_Const(j).X /GLO_Const(j).r(k);
            SVSV = (GLO_Const(j).X - LEO_SV.X); SVSV = SVSV / norm(SVSV);
            AntAngleLEO = rad2deg(acos(LEO_SV.AntM' * SVSV));
            AntAngleGLO= rad2deg(acos(GLO_Const(j).AntM' * (-SVSV)));
            if (abs(AntAngleLEO) <= 85) && (abs(AntAngleGLO) <= 25)
                Sucs_Counter_GLO(k) = Sucs_Counter_GLO(k) + 1;
            end        
        end
    end

    Sucs_Counter = Sucs_Counter_GLO + Sucs_Counter_GPS;

    figure(2)
    plot(t.t, Sucs_Counter, t.t, Sucs_Counter_GPS, t.t, Sucs_Counter_GLO)
    xlabel('t, sec');
    legend('Sum', 'GPS', 'GLO');
    
    min_S(aa) = min(Sucs_Counter);
    max_S(aa) = max(Sucs_Counter);
    proc_S(aa) = sum(Sucs_Counter>=4) / Nmod * 100;

    min_S_GLO(aa) = min(Sucs_Counter_GLO);
    max_S_GLO(aa) = max(Sucs_Counter_GLO);
    proc_S_GLO(aa) = sum(Sucs_Counter_GLO>=4) / Nmod * 100;

    min_S_GPS(aa) = min(Sucs_Counter_GPS);
    max_S_GPS(aa) = max(Sucs_Counter_GPS);
    proc_S_GPS(aa) = sum(Sucs_Counter_GPS>=4) / Nmod * 100;

    h_S(aa) = a - R_e;
    clear LEO_SV;
end

figure(3)
plot(h_S/1000, min_S, h_S/1000, min_S_GPS, h_S/1000, min_S_GLO)
xlabel('H, km')
ylabel('Min number of visible SVs');
legend('Sum', 'GPS', 'GLO');

figure(4)
plot(h_S/1000, max_S, h_S/1000, max_S_GPS, h_S/1000, max_S_GLO)
xlabel('H, km')
ylabel('Max number of visible SVs');
legend('Sum', 'GPS', 'GLO');

figure(5)
plot(h_S/1000, proc_S, h_S/1000, proc_S_GPS, h_S/1000, proc_S_GLO)
xlabel('H, km')
ylabel('The percentage of time during which it is possible positioning (N>=4)');
legend('Sum', 'GPS', 'GLO');