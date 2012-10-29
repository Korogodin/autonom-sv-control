classdef SpaceVehicle_GPS < SpaceVehicle
    %SPACEVEHICLE_GPS Class of GPS sattellites
        
    properties
        Alm
        RL_L1
        
        % Kepler
        M
        E
        theta
        Omega
    end
    
    methods
        function SV = SpaceVehicle_GPS(Name, Alm)
            SV = SV@SpaceVehicle(Name, 'GPS');
            SV.Alm = Alm;
        end
        
        function calcOrbit(SV, t)

            tk = t.t - SV.Alm.t.t + (t - SV.Alm.t); 
            
            Nmod = length(t.t);
          
            % Mean anomaly
            global mu_earth
            if isempty(mu_earth)
                mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
            end
            SV.M = SV.Alm.M0 + sqrt(mu_earth) ./ (SV.Alm.sqrta.^3) .*tk; 
            
            % Eccentricity anomaly
            SV.E = nan(1,Nmod); 
            options_solve = optimset('Display','off');
            SV.E(1) = fsolve(@(E)(E-SV.Alm.e*sin(E) - SV.M(1)), SV.Alm.M0, options_solve);
            for j  = 2:Nmod
                SV.E(j) = fsolve(@(E)(E-SV.Alm.e*sin(E)  - SV.M(j)), SV.E(j-1), options_solve); 
            end

            % True anomaly
            SV.theta = 2*atan(sqrt((1+SV.Alm.e)./(1-SV.Alm.e)).*tan(SV.E/2)); 
            for j = 2:Nmod
                while abs(SV.theta(j) - SV.theta(j-1)) > pi
                    if (SV.theta(j) - SV.theta(j-1)) > pi
                        SV.theta(j) = SV.theta(j) - 2*pi;
                    elseif (SV.theta(j) - SV.theta(j-1)) < -pi
                        SV.theta(j) = SV.theta(j) + 2*pi;
                    end
                end
            end             
            
            SV.Omega = SV.Alm.Omega + SV.Alm.dOmega*tk;             
            
            SV.Alm.a = SV.Alm.sqrta^2;
            SV.Alm.p = SV.Alm.a .* (1-SV.Alm.e.^2);
            munapi = sqrt(mu_earth ./ SV.Alm.p);
            Vr = munapi.*SV.Alm.e.*sin(SV.theta); 
            Vu = munapi.*(1+SV.Alm.e.*cos(SV.theta));
            SV.r = SV.Alm.p ./ (1+SV.Alm.e.*cos(SV.theta)); 
            u = SV.theta + SV.Alm.omega; 

            for j = 1:Nmod
                X = U3(-SV.Omega(j))*U1(-SV.Alm.i)*U3(-u(j))*[SV.r(j); 0; 0];
                SV.x(j) = X(1); SV.y(j) = X(2); SV.z(j) = X(3);
            end
            
            SV.Vx = Vr.*SV.x./SV.r - Vu.*(sin(u).*cos(SV.Omega) + cos(u).*sin(SV.Omega).*cos(SV.Alm.i));
            SV.Vy = Vr.*SV.y./SV.r - Vu.*(sin(u).*sin(SV.Omega) - cos(u).*cos(SV.Omega).*cos(SV.Alm.i));
            SV.Vz = Vr.*SV.z./SV.r + Vu.*cos(u).*sin(SV.Alm.i);
            
        end
    end
    
    
    
end

