classdef SpaceVehicle < handle
    %SpaceVehicle Class of space vehicles
   
    properties 
        
        % Cartesian
%         X = [0; 0; 0]; % Coordinates
%         V = [0; 0; 0]; % Velocities
%         vmDN = [0; 0; 0]; % Vector of antenna's axe (Antennas with axial symmetry only)
        x
        y
        z
        Vx
        Vy
        Vz
        
        % Kepler
        a
        p
        e
        Omega
        i
        omega       
        M
        E
        theta
        r
        
        Name
    end
    
    methods
        function SV = SpaceVehicle(a, e, Omega, dOmega, i, omega, M0, str_name)
            if nargin < 8
                error('SpaceVehicle:InvalidInitialization',...
                   'Input argumets must be: p, e, Omega, dOmega, i, omega, M0, name')
            end            
            
            global mu_earth Nmod t
            
            SV.Name = str_name;
            
            SV.a(1:Nmod) = a;
            SV.e(1:Nmod) = e;
            SV.p = SV.a .* (1-SV.e.^2);
            
            if (e == 0)
                SV.e(1:Nmod) = 0;
                SV.theta = sqrt(mu_earth ./ SV.a.^3).*t + M0;
            else % Calc true anomaly
                % Mean anomaly
                SV.M = M0 + sqrt(mu_earth ./ SV.a.^3).*t; 
                
                % Extr anomaly
                SV.E = nan(1,Nmod); 
                options_solve = optimset('Display','off');
                SV.E(1) = fsolve(@(E)(E-SV.e(1)*sin(E) - SV.M(1)), M0, options_solve);
                for j  = 2:length(t)
                    SV.E(j) = fsolve(@(E)(E-SV.e(j)*sin(E)  - SV.M(j)), SV.E(j-1), options_solve); 
                end

                % True anomaly
                SV.theta = 2*atan(sqrt((1+SV.e)./(1-SV.e)).*tan(SV.E/2)); 
                for j = 2:Nmod
                    while abs(SV.theta(j) - SV.theta(j-1)) > pi
                        if (SV.theta(j) - SV.theta(j-1)) > pi
                            SV.theta(j) = SV.theta(j) - 2*pi;
                        elseif (SV.theta(j) - SV.theta(j-1)) < -pi
                            SV.theta(j) = SV.theta(j) + 2*pi;
                        end
                    end
                end                
            end
            
            % Not pertubed motion yet
            
            SV.omega(1:Nmod) = omega;
            SV.i(1:Nmod) = i;
            SV.Omega(1:Nmod) = Omega + dOmega*t; 
                
            SV.x = nan(1,Nmod); SV.y = nan(1,Nmod); SV.z = nan(1,Nmod);
            SV.Vx = nan(1,Nmod); SV.Vy = nan(1,Nmod); SV.Vz = nan(1,Nmod);
            
            SV.calcCartesian;
        end
        
        function l = dist(obj1, obj2)
            % Calculate distance between two SV
            l = sqrt((obj1.X - obj2.X)'*(obj1.X - obj2.X));
            return;
        end
        
        function diff_r = dr(obj1, obj2)
            % Radius-vector between SVs (obj2 - this)
            diff_r = obj2.X - obj1.X;
            return;
        end

        function a = angl(obj1, obj2)
            % Calculate angle between antenna axe and vision line
            if obj1.dist(obj2) > 0
                a = real(acos( obj1.vmDN' * obj1.dr(obj2)/obj1.dist(obj2) ));
            else
                disp('SVs are in one point of space');
            end
            return;
        end        
        
        function calcCartesian(SV)
            
            global mu_earth
            p = SV.p;
            e = SV.e;
            theta = SV.theta;
            Omega = SV.Omega;
            omega = SV.omega;
            i = SV.i;
            
            munapi = sqrt(mu_earth ./ p);
            Vr = munapi.*e.*sin(theta); 
            Vu = munapi.*(1+e.*cos(theta));
            r = p ./ (1+e.*cos(theta)); SV.r = r;
            u = theta + omega; 

            for j = 1:length(theta)
                X = U3(-Omega(j))*U1(-i(j))*U3(-u(j))*[r(j); 0; 0];
                SV.x(j) = X(1); SV.y(j) = X(2); SV.z(j) = X(3);
            end
            
            SV.Vx = Vr.*SV.x./r - Vu.*(sin(u).*cos(Omega) + cos(u).*sin(Omega).*cos(i));
            SV.Vy = Vr.*SV.y./r - Vu.*(sin(u).*sin(Omega) - cos(u).*cos(Omega).*cos(i));
            SV.Vz = Vr.*SV.z./r + Vu.*cos(u).*sin(i);
        end
      
    end
    
end

