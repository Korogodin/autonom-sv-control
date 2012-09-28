classdef SpaceVehicle < handle
    %SpaceVehicle Class of space vehicles
    
    properties
        
        % Cartesian, equatorial
        X
        V
        x
        y
        z
        Vx
        Vy
        Vz
      
        r 
        
        Name
        
        AntM
    end
    
    methods
        function SV = SpaceVehicle(str_name)
            if nargin < 1
                error('SpaceVehicle:InvalidInitialization',...
                    'Input argumets must be: Name')
            end
            
            SV.Name = str_name;
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
        
    end
    
end

% function out = resize_to_t(arr, Nmod)
%     out = NaN;
%     if isempty(arr)
%         return
%     end
%     if isnan(arr) || isempty(arr)
%         return;
%     else
%         if length(arr) ~= Nmod
%             out = arr(1)*ones(1, Nmod);
%         else
%             out = arr;
%         end
%     end
% end
