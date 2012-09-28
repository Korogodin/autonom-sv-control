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
t = CTime(0:dTmod:Tmod, 26, 09, 2012);
Nmod = length(t.t);

GPS_Const = LoadGPSConst([pwd '/data/Topcon_120925.agp'], t);
figure(1)
hold on
for j = 1:length(GPS_Const)
    plot3(GPS_Const(j).x, GPS_Const(j).y, GPS_Const(j).z, 'y')
    plot3(GPS_Const(j).x(1), GPS_Const(j).y(1), GPS_Const(j).z(1), '*b', 'MarkerSize', 10)
%     fprintf('%f\n', GPS_Const(j).Alm.Omega);
end
hold off

aa = 0;
for a = (R_e+0e3):1e3:(R_e + 3.5e6)
aa = aa+1;

fprintf('Calculation for LEO satt...');
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


Sucs_Counter = zeros(1, Nmod);
for k = 1:Nmod
    N_GPS = length(GPS_Const);
    N_GLO = 0;
    
    LEO_SV.X = [LEO_SV.x(k); LEO_SV.y(k); LEO_SV.z(k)];
    LEO_SV.AntM = LEO_SV.X/LEO_SV.r(k);
    for j = 1:N_GPS
        GPS_Const(j).X = [GPS_Const(j).x(k); GPS_Const(j).y(k); GPS_Const(j).z(k)];
        GPS_Const(j).AntM = -GPS_Const(j).X /GPS_Const(j).r(k);
        SVSV = (GPS_Const(j).X - LEO_SV.X); SVSV = SVSV / norm(SVSV);
        AntAngleLEO = rad2deg(acos(LEO_SV.AntM' * SVSV));
        AntAngleGPS = rad2deg(acos(GPS_Const(j).AntM' * (-SVSV)));
        if (AntAngleLEO < 80) && (AntAngleGPS < 20)
            Sucs_Counter(k) = Sucs_Counter(k) + 1;
        end
    end
    
end
figure(2)
plot(t.t, Sucs_Counter)

min_S(aa) = min(Sucs_Counter);
max_S(aa) = max(Sucs_Counter);
proc_S(aa) = sum(Sucs_Counter>=4) / Nmod * 100;
h_S(aa) = a - R_e;
clear LEO_SV;
end

figure(3)
plot(h_S/1000, min_S, h_S/1000, max_S)
xlabel('H, km')
ylabel('Min and Max number of visible GPS SV');

figure(4)
plot(h_S/1000, proc_S)
xlabel('H, km')
ylabel('The percentage of time during which it is possible positioning (N>=4)');