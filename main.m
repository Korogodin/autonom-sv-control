try 
    close(handle_fig_main); % Close old output form
end

close all
clear 
clc

globals;
figures_handles;
addpath([pwd '/data/']);

load RadPatt.mat % Radiation patterns for transmitter's and receiver's antennas

mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
omega_e = 0.7292115e-4; % [rad/s] Earth's rotation rate
R_e = 6371e3;  % [m] Mean Earth radiius
        
% Tmod = 48*60*60; % Model time
% dTmod = 600; % time step
Tmod = 30*60*60; % Model time
dTmod = 100; % time step
t = CTime(0:dTmod:Tmod, 24, 09, 2012);
Nmod = length(t.t);

% Sensivity for NV08C
SensLoopGPS = -185; % dBWt
SensSearchGPS = -167;
SensLoopGLO = -182;
SensSearchGLO = -168;

% Load GLONASS constellation from almanac
GLO_Const = LoadGLOConst([pwd '/data/Topcon_120925.agl'], t);
N_GLO = length(GLO_Const);
figure(hF_sat);
pos = get(hF_sat, 'Position'); 
pos(3:4) = pos(3:4)*2;
set(hF_sat, 'Position', pos);
hold on
lambda =  3e8 / 1.602e9;
TransPower = 18;
for j = 1:N_GLO
    GLO_Const(j).RL_L1 = CRadioLine(lambda, Gr_L1, Gt_GLOL1, TransPower, SensLoopGLO, SensSearchGLO);
    plot3(GLO_Const(j).x, GLO_Const(j).y, GLO_Const(j).z, 'b'); % m
    plot3(GLO_Const(j).x(end), GLO_Const(j).y(end), GLO_Const(j).z(end), '*b', 'MarkerSize', 12)
%     fprintf('%f\n', GLO_Const(j).Alm.Omega);
end
hold off


% Load GPS constellation from almanac
GPS_Const = LoadGPSConst([pwd '/data/Topcon_120925.agp'], t);
N_GPS =  length(GPS_Const);
figure(hF_sat)
hold on
lambda =  3e8 / 1.575e9;
TransPower = 17.5;
for j = 1:N_GPS
    GPS_Const(j).RL_L1 = CRadioLine(lambda, Gr_L1, Gt_GPSL1, TransPower, SensLoopGPS, SensSearchGPS);
    plot3(GPS_Const(j).x, GPS_Const(j).y, GPS_Const(j).z, 'g'); % y
    plot3(GPS_Const(j).x(end), GPS_Const(j).y(end), GPS_Const(j).z(end), '*g', 'MarkerSize', 12)
%     fprintf('%f\n', GPS_Const(j).Alm.Omega);
end
xlabel('X_0'); ylabel('Y_0'); zlabel('Z_0');
hold off


% Earth sphere
R_earth_pol = 6356777; % Polar radius of Earth
R_earth_equa = 6378160; % Equatorial radius of Earth
Earth_axe_angl = deg2rad(23); 
[x_e,y_e,z_e] = sphere(50);
x_e = R_earth_equa * x_e; y_e = R_earth_equa * y_e; z_e = R_earth_pol * z_e;
hF = 0;
hF = figure(hF + 1);
hold on; surface(x_e, y_e, z_e); hold off;

Type = 'LEO';
CNS = LoadCNS( Type, t );

Hvect = R_e:25e3:R_e+40000e3;

if length(Hvect) > 1
    a = R_e + Hvect;
    Na = length(a);
    Availability_GLO = nan(1, length(Hvect));
    Availability_GPS = nan(1, length(Hvect));
    Availability_Comm = nan(1, length(Hvect));
    
    GoodDOP_GLO = nan(1, length(Hvect));
    GoodDOP_GPS = nan(1, length(Hvect));
    GoodDOP_Comm = nan(1, length(Hvect));
    
    MaxPower_GLO = nan(1, length(Hvect));
    MaxPower_GPS = nan(1, length(Hvect));
else 
    Na = 1;
end

for ja = 1:Na
    
    
    if Na > 1
        fprintf('Calculate for h = %.0f km...\n', (Hvect(ja)-R_e)/1e3);
        CNS.a = Hvect(ja);
        CNS.calcKeplerOrbit(t);
    end
%         CNS.a =R_earth_equa + 4000e3;
%         CNS.calcKeplerOrbit(t);    
    
    figure(hF_sat)
    hold on
    plot3(CNS.x, CNS.y, CNS.z, 'r', 'LineWidth', 2)
    plot3(CNS.x(end), CNS.y(end), CNS.z(end), '*r', 'MarkerSize', 12)
    hold off
    
    % Radioline GLO-CNS
    fprintf('Calculate radioline for GLO...\n');
    for js = 1:N_GLO
        GLO_Const(js).RL_L1.CalcRL(t, GLO_Const(js), CNS);
    end
        figure(hF_P_GLO);
        cla
        spr = 'plot(t.t, GLO_Const(1).RL_L1.Power';
        spr_leg = ['legend(''GLONASS ' GLO_Const(1).Name ''''];
        for j = 2:N_GLO
            spr = [spr sprintf(', t.t, GLO_Const(%.0f).RL_L1.Power', j)];
            spr_leg = [spr_leg ', ''GLONASS ' GLO_Const(j).Name ''''];
        end
        spr = [spr ');'];
        spr_leg = [spr_leg ');'];
        eval(spr);
        eval(spr_leg);
        xlabel('t, s');
        ylabel('P_{ant}, dBWt')
        title('GLO')        

    
    
    % Radioline GPS-CNS
    fprintf('Calculate radioline for GPS...\n');
    for js = 1:N_GPS
        GPS_Const(js).RL_L1.CalcRL(t, GPS_Const(js), CNS);
    end
        figure(hF_P_GPS);
        cla
        spr = 'plot(t.t, GPS_Const(1).RL_L1.Power';    
        spr_leg = ['legend(''GPS ' GLO_Const(1).Name ''''];
        for j = 2:N_GPS
            spr = [spr sprintf(', t.t, GPS_Const(%.0f).RL_L1.Power', j)];
            spr_leg = [spr_leg ', ''GPS ' GPS_Const(j).Name ''''];
        end
        spr = [spr ');'];
        spr_leg = [spr_leg ');'];
        eval(spr);
        eval(spr_leg);
        xlabel('t, s');
        ylabel('P_{ant}, dBWt')
        title('GPS')    
    
    fprintf('Calculate DOP...\n');
    CNS.calcDOP(t, GLO_Const, GPS_Const);
        figure(hF_DOP)
        cla
        plot(t.t, CNS.DOP_GLO, t.t, CNS.DOP_GPS, t.t, CNS.DOP_Comm);
        xlabel('t, sec')
        ylabel('GDOP');
        legend('GLONASS', 'GPS', 'GLONASS+GPS');
    

    fprintf('Calculate number of satellites...\n');
    CNS.calcSatNums(t, GLO_Const, GPS_Const);
        figure(hF_Num)
        cla
        plot(t.t, CNS.Num_Sync_GLO, t.t, CNS.Num_Sync_GPS, t.t, CNS.Num_Sync_GLO+CNS.Num_Sync_GPS);
        xlabel('t, sec')
        ylabel('N');
        legend('GLONASS', 'GPS', 'GLONASS+GPS');
    
        figure(hF_MinIn4)
        cla
        plot(t.t, CNS.MinIn4_GLO, t.t, CNS.MinIn4_GPS);
        xlabel('t, sec')
        ylabel('Minimum P_{ant} in working four, dBWt');
        legend('GLONASS', 'GPS');
        
%         figure(hF_Angles)
%         cla
%         plot(t.t, GLO_Const(21).RL_L1.alpha, t.t, GLO_Const(21).RL_L1.beta, t.t, GPS_Const(17).RL_L1.alpha, t.t, GPS_Const(17).RL_L1.beta);
%         xlabel('t, sec')
%         ylabel('Angles to vision line');
%         legend('GLONASS GL21 \alpha', 'GLONASS GL21 \beta', 'GPS GP17 \alpha', 'GPS GP17 \beta');
%         
%         figure(hF_Gain)
%         cla
%         plot(t.t, GLO_Const(21).RL_L1.G_rec, t.t, GLO_Const(21).RL_L1.G_tr, t.t, GPS_Const(17).RL_L1.G_rec, t.t, GPS_Const(17).RL_L1.G_tr);
%         xlabel('t, sec')
%         ylabel('Antenna''s gain');
%         legend('GLONASS GL21 Gr', 'GLONASS GL1 Gt', 'GPS GP17 Gr', 'GPS GP17 Gt');        
%         
%         % Histogram of angles
%         x_hist = 0:180;
%         alphas = nan(1, N_GLO*Nmod);
%         betas = nan(1, N_GLO*Nmod);
%         for jg = 1:N_GLO
%             alphas((jg-1)*Nmod+1:jg*Nmod) = GLO_Const(jg).RL_L1.alpha;
%             betas((jg-1)*Nmod+1:jg*Nmod) = GLO_Const(jg).RL_L1.beta;
%         end        
%         figure(hF_Hist_alpha_GLO)
%         cla
%         hist(alphas, x_hist)
%         xlabel('\alpha_{GLO}, deg');
%         xlim([0 180]);
%         figure(hF_Hist_beta_GLO)
%         cla
%         hist(betas, x_hist)
%         xlabel('\beta_{GLO}, deg');
%         xlim([0 180]);
%         
%         alphas = nan(1, N_GPS*Nmod);
%         betas = nan(1, N_GPS*Nmod);
%         for jg = 1:N_GPS
%             alphas((jg-1)*Nmod+1:jg*Nmod) = GPS_Const(jg).RL_L1.alpha;
%             betas((jg-1)*Nmod+1:jg*Nmod) = GPS_Const(jg).RL_L1.beta;
%         end        
%         figure(hF_Hist_alpha_GPS)
%         cla
%         hist(alphas, x_hist)
%         xlabel('\alpha_{GPS}, deg');
%         xlim([0 180]);
%         figure(hF_Hist_beta_GPS)
%         cla
%         hist(betas, x_hist)
%         xlabel('\beta_{GPS}, deg');          
%         xlim([0 180]);
%         
%         beta_arr = 0:1:180;
%         EnergyFromBeta_GLO = zeros(1, 181);
%         EnergyFromBeta_GPS = zeros(1, 181);
%         for kb = 1:Nmod
%             for jg = 1:N_GLO
%                 if ~isnan(GLO_Const(jg).RL_L1.beta(kb))
%                     EnergyFromBeta_GLO(round(GLO_Const(jg).RL_L1.beta(kb))+1) = ...
%                        EnergyFromBeta_GLO(round(GLO_Const(jg).RL_L1.beta(kb))+1) + ...
%                        10^(GLO_Const(jg).RL_L1.Power_noNaN(kb)/10);
%                 end
%             end
%             for jg = 1:N_GPS
%                 if ~isnan(GPS_Const(jg).RL_L1.beta(kb))
%                     EnergyFromBeta_GPS(round(GPS_Const(jg).RL_L1.beta(kb))+1) = ...
%                        EnergyFromBeta_GPS(round(GPS_Const(jg).RL_L1.beta(kb))+1) + ...
%                        10^(GPS_Const(jg).RL_L1.Power_noNaN(kb)/10);
%                 end
%             end
%         end
%         EnergyFromBeta_GLO = 10*log10(EnergyFromBeta_GLO);
%         EnergyFromBeta_GPS = 10*log10(EnergyFromBeta_GPS);
%         figure(hF_EnergyFromBeta_GLO)
%         bar(beta_arr, EnergyFromBeta_GLO  + 220)
%         ylabel('Energy histogram for GLONASS, dB')
%         xlabel('beta, deg');
%         xlim([0 180]);
%         figure(hF_EnergyFromBeta_GPS)
%         bar(beta_arr, EnergyFromBeta_GPS + 220)
%         ylabel('Energy histogram for GPS, dB')
%         xlabel('beta, deg');
%         xlim([0 180]);
        
        Availability_GLO(ja) = sum(CNS.Num_Sync_GLO >= 4) / Nmod;
        Availability_GPS(ja) = sum(CNS.Num_Sync_GPS >= 4) / Nmod;
        Availability_Comm(ja) = sum( (CNS.Num_Sync_GPS >= 3)&(CNS.Num_Sync_GLO >= 2) |...
                                                   (CNS.Num_Sync_GPS >= 2)&(CNS.Num_Sync_GLO >= 3) | ...
                                                   (CNS.Num_Sync_GPS >= 4) | ...
                                                   (CNS.Num_Sync_GLO >= 4)) / Nmod;
       
        MaxPower_GLO(ja) = max(CNS.MinIn4_GLO);
        MaxPower_GPS(ja) = max(CNS.MinIn4_GPS);
    
        GoodDOP_GLO(ja) = sum(CNS.DOP_GLO <= 20) / Nmod;
        GoodDOP_GPS(ja) = sum(CNS.DOP_GPS <= 20) / Nmod;
        GoodDOP_Comm(ja) = sum(CNS.DOP_Comm <= 20) / Nmod;
end


if Na > 1
    Hvect = Hvect - R_e;
    hF = figure(hF + 1);
    plot(Hvect, Availability_GLO*100, Hvect, Availability_GPS*100, Hvect, Availability_Comm*100)
    xlabel('H, m');
    ylabel('Measurament availability percent');
    legend('GLONASS', 'GPS', 'GLONASS+GPS');
    
    hF = figure(hF + 1);
    plot(Hvect, MaxPower_GLO, Hvect, MaxPower_GPS)
    xlabel('H, m');
    ylabel('Max Signal Power for 4 satellites');
    legend('GLONASS', 'GPS');    
    
    hF = figure(hF + 1);
    plot(Hvect, GoodDOP_GLO*100, Hvect, GoodDOP_GPS*100, Hvect, GoodDOP_Comm*100)
    xlabel('H, m');
    ylabel('GDOP<=20 PerCent');
    legend('GLONASS', 'GPS', 'GLONASS+GPS');       
end