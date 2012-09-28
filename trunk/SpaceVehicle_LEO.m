classdef SpaceVehicle_LEO < SpaceVehicle
    %SPACEVEHICLE_LEO Class of LEO sattellites
        
    properties
        p
        a
        e
        Omega
        omega
        i
        M0
        t0
        
        GPS_Power
    end
    
    methods
        function SV = SpaceVehicle_LEO(Name, a, e, Omega, omega, i, M0, t0)
            SV = SV@SpaceVehicle(Name);
            SV.a = a;
            SV.e = e;
            SV.Omega = Omega;
            SV.omega = omega;
            SV.i = i;
            SV.M0 = M0;
            SV.t0 = t0; SV.t0.t = SV.t0.t(1);
        end
        
        function calcKeplerOrbit(SV, t)

            Nmod = length(t.t);
            global mu_earth
            
            e = SV.e;
            a = SV.a;
            p = a .* (1-e.^2);
            i = SV.i;
            Omega = SV.Omega;
            omega = SV.omega;
            M0 = SV.M0;
            
            % Mean anomaly
            if isempty(mu_earth)
                mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
            end
            
            tk = t.t - SV.t0.t + (t - SV.t0);
            M = M0 + sqrt(mu_earth) ./ (sqrt(a).^3) .*tk;
            
            % Eccentricity anomaly
            E = nan(1,Nmod);
            options_solve = optimset('Display','off');
            E(1) = fsolve(@(E)(E-e*sin(E) - M(1)), M0, options_solve);
            for j  = 2:Nmod
                E(j) = fsolve(@(E)(E-e*sin(E)  - M(j)), E(j-1), options_solve);
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
            
            munapi = sqrt(mu_earth ./ p);
            Vr = munapi.*e.*sin(theta); 
            Vu = munapi.*(1+e.*cos(theta));
            SV.r = p ./ (1+e.*cos(theta)); 
            u = theta + omega; 

            for j = 1:Nmod
                X = U3(-Omega)*U1(-i)*U3(-u(j))*[SV.r(j); 0; 0];
                SV.x(j) = X(1); SV.y(j) = X(2); SV.z(j) = X(3);
            end
            
            SV.Vx = Vr.*SV.x./SV.r - Vu.*(sin(u).*cos(Omega) + cos(u).*sin(Omega).*cos(i));
            SV.Vy = Vr.*SV.y./SV.r - Vu.*(sin(u).*sin(Omega) - cos(u).*cos(Omega).*cos(i));
            SV.Vz = Vr.*SV.z./SV.r + Vu.*cos(u).*sin(i);

        end
        
    end
    
    
    
end

