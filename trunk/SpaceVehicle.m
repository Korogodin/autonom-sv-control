classdef SpaceVehicle
    %SpaceVehicle Class of space vehicles
        
    properties 
        % Dekart
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
    end
    
    methods
        function SV = SpaceVehicle(t, p, e, Omega, i, omega)
            if nargin < 6
            error('SpaceVehicle:InvalidInitialization',...
               'Must provide an account number and initial balance')
            end            
            SV.X = [t;t;t];
            SV.H_earth = sqrt(SV.X'*SV.X);
            SV.vmDN = SV.X / SV.H_earth;
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
        
        
    end
    
end

