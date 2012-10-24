classdef SpaceVehicle_GLO < SpaceVehicle
    %SPACEVEHICLE_GLO Class of GLO sattellites
        
    properties
        Alm
        
        % Kepler
%         M
%         E
%         theta
%         Omega
        
    end
    
    methods
        function SV = SpaceVehicle_GLO(Name, Alm)
            SV = SV@SpaceVehicle(Name);
            SV.Alm = Alm;
        end
        
        function calcOrbit(SV, t)

            t_star = t.t - SV.Alm.t.t + (t - SV.Alm.t); 
            
            Nmod = length(t.t);
          
            % Mean anomaly
            global mu_earth
            if isempty(mu_earth)
                mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
            end
            
            i_mean = deg2rad(63); % Mean inclination for GLONASS
            T_mean = 43200; % [s] Mean draconian period for GLONASS 
            a_e = 6378136; % [m] Equatorial radius of Earth
            C20 = -1082625.75e-9; % second zone garmonic
            
            Tdr = T_mean + SV.Alm.dT; % draconian period
            i = i_mean + SV.Alm.di; % inclination
            e = SV.Alm.e; % Eccentricity
            omega = SV.Alm.omega; % Argument of perigee
            
            % Semi-major axis
            a_old = 0;
            a = ( (Tdr/2/pi)^2 * mu_earth )^(1/3); % semi-major axis
            while (a - a_old) > 1 % [m]
                a_old = a;
                p = a*(1-e^2); % Focal parameter
                Tosc = Tdr*1/(1 + 1.5*C20*(a_e/p)^2 * ...
                    ( (2 - 2.5*sin(i)^2)*(1-e^2)^(3/2)/(1+e*cos(omega))^2 + (1 + e*cos(-omega))^3 / (1 - e^2) ) );
                a = ( (Tosc/2/pi)^2 * mu_earth )^(1/3);
            end
            p = a*(1-e^2);
            
            
            % Omega and t_lambda for u=0
            W = fix(t_star(1) / Tdr); % Recalc to model time 
            tsr = SV.Alm.t.t + Tdr*W + SV.Alm.ddT*W.^2;
            
            global omega_e
            if isempty(omega_e)
                omega_e = 0.7292115e-4; % [rad/s] Earth's rotation rate
            end            
            
            Omega_h = 1.5 * C20 * 2*pi / Tdr * (a_e / a)^2 * cos(i) / (1 - e^2)^2;
            lak = SV.Alm.Omega + (Omega_h - omega_e)*(Tdr*W + SV.Alm.ddT*W.^2);
            
            S = SV.Alm.t.S0 + omega_e*(mod(tsr, 86400) - 10800); % -3 hour (Moscow)
            Omega0 = lak + S; 
            
            %Corrections
            
            %
            
            tk = t.t - tsr + (t - SV.Alm.t); % time from u=0 on t.t(1) turn
            
            M = 0 + sqrt(mu_earth ./ a.^3) .*tk; % if t == SV.Alm.t.tsr) then M=0; E=0; u=0.
%             M = 0 + sqrt(mu_earth ./ a.^3) .*(t.t - t.t(1)); % if t == SV.Alm.t.tsr) then M=0; E=0; u=0.
            
%             global Omega1 Omega2 Omega3 cOmega1 cOmega2 cOmega3
%             if Omega1 == 0
%                 Omega1 = Omega0;
%                 cOmega1 = 0;
%                 nO = cOmega1;
%             else
%                 if mod(abs(Omega0 - Omega1), 2*pi) < pi/6
%                     cOmega1 = cOmega1 + 1;
%                     nO = cOmega1;
%                 else
%                     if Omega2 == 0
%                         Omega2 = Omega0;
%                         cOmega2 = 0;
%                         nO = cOmega2;
%                     else
%                         if mod(abs(Omega0 - Omega2), 2*pi) < pi/6
%                             cOmega2 = cOmega2 + 1;
%                             nO = cOmega2;                        
%                         else
%                             if Omega3 == 0
%                                 Omega3 = Omega0;
%                                 cOmega3 = 0;
%                                 nO = cOmega3;
%                             else
%                                 if mod(abs(Omega0 - Omega3), 2*pi) < pi/6
%                                     cOmega3 = cOmega3 + 1;
%                                     nO = cOmega3;          
%                                 else
%                                     disp('Error');
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%             M = M + nO * 2*pi/8;
                        
            % Eccentricity anomaly
            E = nan(1, Nmod); 
            options_solve = optimset('Display','off');
            E(1) = fsolve(@(Ex)(Ex-e*sin(Ex) - M(1)), M(1), options_solve);
            for j  = 2:Nmod
                E(j) = fsolve(@(Ex)(Ex-e*sin(Ex)  - M(j)), E(j-1), options_solve); 
            end

            % True anomaly
            theta = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2)); 
            for j = 2:Nmod
                while abs(theta(j) - theta(j-1)) > pi
                    if (theta(j) - theta(j-1)) > pi
                        theta(j) = theta(j) - 2*pi;
                    elseif (theta(j) - theta(j-1)) < -pi
                        theta(j) = theta(j) + 2*pi;
                    end
                end
            end             
            
            Omega = Omega0 + Omega_h*tk;       
%             Omega = Omega0 + zeros(1, length(tk));
            
            munapi = sqrt(mu_earth ./ p);
            Vr = munapi.*e.*sin(theta); 
            Vu = munapi.*(1+e.*cos(theta));
            SV.r = p ./ (1+e.*cos(theta)); 
%             u = theta + omega; 
            u = theta;

            for j = 1:Nmod
                X = U3(-Omega(j))*U1(-i)*U3(-u(j))*[SV.r(j); 0; 0];
                SV.x(j) = X(1); SV.y(j) = X(2); SV.z(j) = X(3);
            end
            
            SV.Vx = Vr.*SV.x./SV.r - Vu.*(sin(u).*cos(Omega) + cos(u).*sin(Omega).*cos(i));
            SV.Vy = Vr.*SV.y./SV.r - Vu.*(sin(u).*sin(Omega) - cos(u).*cos(Omega).*cos(i));
            SV.Vz = Vr.*SV.z./SV.r + Vu.*cos(u).*sin(i);
            
        end
    end
    
    
    
end

