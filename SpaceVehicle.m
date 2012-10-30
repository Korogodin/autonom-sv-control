classdef SpaceVehicle < handle
    %SpaceVehicle Class of space vehicles
    
    properties
        
        % Cartesian, equatorial
        X % 3x1
        V % 3x1
        x % 1:K
        y % 1:K
        z % 1:K
        Vx % 1:K
        Vy % 1:K
        Vz % 1:K
      
        r % 1:K
        
        Name
        Type
        
        AntM % 3x1
    end
    
    methods
        function SV = SpaceVehicle(str_name, str_type)
            if nargin < 2
                error('SpaceVehicle:InvalidInitialization',...
                    'Input argumets must be: Name, Type')
            end
            
            SV.Name = str_name;
            SV.Type = str_type; % GLO, GPS, GEO, LEO, HElO and others
        end
        
        function touchX(SV1, k)
            SV1.X =  [SV1.x(k); SV1.y(k); SV1.z(k)];
        end
        
        function l = dist(SV1, SV2)
            % Calculate distance between two SV
            l = norm(SV2.X - SV1.X);
        end
        
        function diff_r = dr(obj1, obj2)
            % Radius-vector between SVs (obj2 - this)
            diff_r = obj2.X - obj1.X;
            return;
        end
        
        function CalcAntM(SV, k)
            % Calculate antenna direction for t.t(k) 
            % SV.X, SV.r(k) must be actual
            switch SV.Type
                case 'LEO' % Low Earth Orbit
                    SV.AntM = -SV.X/SV.r(k); % TEEEEEMPPP
                case {'GEO', 'MEO'} % Geostationary Earth Orbit or Middle Earth Orbit
                    SV.AntM = -SV.X/SV.r(k);
                case 'HElO' % High Elliptical Orbit
                    SV.AntM = -[SV.x(1); SV.y(1); SV.z(1)]/SV.r(1);
                case {'GPS', 'GLO'}
                    SV.AntM = -SV.X/SV.r(k);
                otherwise
                    disp('Unknown Type of satellite');
            end
        end
       
    end
    
end
