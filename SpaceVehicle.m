classdef SpaceVehicle
    %SpaceVehicle Class of space vehicles

    properties (Access=private)
        mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
        omega_e = 0.7292115e-4; % [rad/s] Earth's rotation rate
    end
    
    properties 
        
        % Cartesian
        X = [0; 0; 0]; % Coordinates
        V = [0; 0; 0]; % Velocities
        vmDN = [0; 0; 0]; % Vector of antenna's axe (Antennas with axial symmetry only)
        
        % Kepler
        p
        e
        Omega
        i
        omega        
        theta
        
        H_earth % dist to Earth's center
        
        Name
    end
    
    methods
        function SV = SpaceVehicle(p, e, Omega, i, omega, theta, str_name)
            if nargin < 7
                error('SpaceVehicle:InvalidInitialization',...
                   'Input argumets must be: p, e, Omega, i, omega, theta, name')
            end            
            
            globals;
            
            SV.Name = str_name;
            
            
            SV.X = [Omega;Omega;Omega];
            
            SV.H_earth = sqrt(SV.X'*SV.X);
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
                disp('SVs in one point of space');
            end
            return;
        end        
        
        function getCartesian(obj, j)
            
            mu_earth = obj.mu_earth;
            p = obj.p(j);
            e = obj.e(j);
            theta = obj.theta(j);
            Omega = obj.Omega(j);
            omega = obj.omega(j);
            i = obj.i(j);
            
            munapi = sqrt(mu_earth / p);
            Vr = munapi*e*sin(theta); 
            Vu = munapi*(1+e*cos(theta));
            r = p / (1+e*cos(theta));
            u = theta + omega; 

            obj.X = U3(-Omega)*U1(-i)*U3(-u)*[r; 0; 0];

            x = obj.X(1); y = obj.X(2); z = obj.X(3);
            
            obj.V(1) = Vr.*x./r - Vu.*(sin(u).*cos(Omega) + cos(u).*sin(Omega).*cos(i));
            obj.V(2) = Vr.*y./r - Vu.*(sin(u).*sin(Omega) - cos(u).*cos(Omega).*cos(i));
            obj.V(3) = Vr.*z./r + Vu.*cos(u).*sin(i);
        end
        
    end
    
end

