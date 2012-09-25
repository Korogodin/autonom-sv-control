function GPS_Const = CreateGPS()
%CREATEGPS Create GPS constellation

% incl = deg2rad(55); % [rad]
% e = 0;
% a = 26559.7e3; % [m]
% p = a*(1-e^2); % [m]
% omega = 0;

% orb = importdata([pwd '/data/GPS_Orbits.csv'], ';');
% orb.data = deg2rad(orb.data);

AlmanacGPS;

fprintf('Begin GPS constellation creation...\n');

orb_old = -1; stri = 'plot3(NaN, NaN, NaN';
col_stri = 'ymcrgbbbbbbbbbbbbbbbbbbbbbbbbbbbb'; col_stri_i = 0;
col_stri = 'yyyyyyyyyyyyyyyyyyyyyyyyyyyyy'; col_stri_i = 0;

orb_old = deg2rad(AlmGPS(15, 7));
for j = 1:size(AlmGPS, 1)
    
    e = AlmGPS(j, 3);
    a = AlmGPS(j, 6);
    incl = deg2rad(AlmGPS(j, 4));
    omega =  deg2rad(AlmGPS(j, 8));
    Omega = deg2rad(AlmGPS(j, 7));
    dOmega = deg2rad(AlmGPS(j, 5));
    m0 = deg2rad(AlmGPS(j, 9));
    
    GPS_Const(j) = SpaceVehicle(a, e, Omega, dOmega, incl, omega, m0, num2str(AlmGPS(j, 1)));
    fprintf('\t%s\tOK\n', num2str(AlmGPS(j, 1)));
    if abs(Omega - orb_old) < deg2rad(10)
        col_stri_i = col_stri_i + 1;
        stri = sprintf('%s, GPS_Const(%.0f).x, GPS_Const(%.0f).y, GPS_Const(%.0f).z, ''%s''', stri, j, j, j, col_stri(col_stri_i));
        stri = sprintf('%s, GPS_Const(%.0f).x(1), GPS_Const(%.0f).y(1), GPS_Const(%.0f).z(1), ''b*''',stri, j, j, j);
%         orb_old = Omega;
    end

end
stri = [stri ')'];
eval(stri);
fprintf('\tGPS constellation creation complete\n');

end

