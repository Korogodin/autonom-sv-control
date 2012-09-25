function GPS_Const = CreateGPS()
%CREATEGPS Create GPS constellation

incl = deg2rad(55); % [rad]
e = 0;
a = 26559.7e3; % [m]
p = a*(1-e^2); % [m]
omega = 0;

orb = importdata([pwd '/data/GPS_Orbits.csv'], ';');

fprintf('Begin GPS constellation creation...\n');

for j = 1:length(orb.data)
    GPS_Const(j) = SpaceVehicle(p, e, orb.data(j, 1), incl, omega, orb.data(j ,2), orb.textdata{j+1, 1});
    fprintf('\t%s\tOK\n', orb.textdata{j+1, 1});
end
fprintf('\tGPS constellation creation complete\n');
   

end

