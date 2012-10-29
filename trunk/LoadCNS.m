function SV_CNS = LoadCNS( Type, t )
%LOADCNS Construct SV-consumer of Type

switch Type
    case 'LEO'
        % Космос-2480 — российский разведывательный спутник типа Кобальт-М.
        % http://ru.wikipedia.org/wiki/Космос-2480
        % R_e + 259.5e3 = (1 + e)*a;
        % R_e + 196.4e3 = (1 - e)*a; => a = 6598950; e = 0.0048;
        e = 0.0048;
        a = 6598950;
        Omega = deg2rad(10);
        omega = 0;
        i = deg2rad(81.4);
    case 'GEO'
        %Экспресс-МД1
        % R_e + 35804.5e3 = (1 + e)*a;
        % R_e + 35782.3e3 = (1 - e)*a;
        % a = 42164.4e3;
        % => e = 0.0048;
        a = 42164.4e3;
        e = 0.0048;
        Omega = deg2rad(80);
        omega = 0;
        i = 0;
    case 'HElO'
        % EO
        % Меридиан-4
        % R_e + 40000e3 = (1 + e)*a;
        % R_e + 500e3 = (1 - e)*a;
        % => a = 26621000
        % => e = 0.0048;
        e = 0.7419;
        a = 26621000;
        Omega = deg2rad(70);
        omega = deg2rad(280);
        i = deg2rad(62.8);
    otherwise
        disp('Unknown Type of CNS');
end
    M0 = 0;
    
    SV_CNS =  SpaceVehicle_CNS(Type, Type, a, e, Omega, omega, i, M0, t);
    SV_CNS.calcKeplerOrbit(t);

end

