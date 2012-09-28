function GPS_Const = LoadGPSConst( filename, t )
%LOADGPSCONST Load GPS constellation from almanac

fprintf('Load Almanac for GPS constellation...\n');

fid = fopen(filename);

N_GPS = 32;
AlmSVList = zeros(1, N_GPS);

j = 0;
while (~feof(fid))
    s1 = fgets(fid); 
    if (ischar(s1))
        s2 = fgets(fid); 
        if (ischar(s2))
            s3 = fgets(fid);
            if (ischar(s3))
                j = j + 1;
                [Alm(j) ok] = ParseAGP(s1, s2, s3);
                if (~isnan(Alm(j).PRN))
                    if (AlmSVList(Alm(j).PRN) == 0)
                       AlmSVList(Alm(j).PRN) = j;
                    else
                       if (Alm(j).t - Alm(AlmSVList(Alm(j).PRN)).t > 0) % Newer
                           AlmSVList(Alm(j).PRN) = j;
                       end
                    end
                end
                if ~ok
                    fclose(fid);
                    error('Error of reading of almanac AGP file');
                end
            end
        end
    end
end
fclose(fid);

k = 0;
for (j = 1:N_GPS)
    if AlmSVList(j) > 0
        k = k + 1;
        GPS_Const(k) = SpaceVehicle_GPS(['GP' num2str(Alm(AlmSVList(j)).PRN)], Alm(AlmSVList(j)));
        fprintf('\tSV%02.0f - almanac %s - orbit calculation for %s ...', ...
                Alm(AlmSVList(j)).PRN, Alm(AlmSVList(j)).t.Stri, t.Stri);
        GPS_Const(k).calcOrbit(t); 
        fprintf('OK\n');
    end
end

fprintf('Constellation of %.0f GPS SVs are loaded\n', k);

end

function [Alm ok] = ParseAGP(s1, s2, s3)
%PARSEAGP Parser for AGP-files
%See: ftp://ftp.glonass-iac.ru/MCC/FORMAT/Format.agp
%       http://www.navipedia.net/index.php/GPS_and_Galileo_Satellite_Coordinates_Computation

ok = 0;

% First string:
% 12 12 1234  12345  1234 123456  123456789101234567890
% 1 - число получения альманаха (Almanac reciept day)
% 2 - месяц получения альманаха (Almanac reciept month)
% 3 - год получения альманаха (Almanac reciept year)
% 4 - время получения альманаха от начала суток, с UTC (Time of almanac reciept from begining of day)
% 5 - неделя GPS (получения альманаха) (GPS week of --//--)
% 6 - время недели GPS, с (получения альманаха) (GPS sec of --//--)
% 7 - комментарий (приемник с которого получено, версия SW и т.д.) (Comment)
ns1 = str2num(s1);
sns1 = length(ns1);
if sns1 > 7
    Alm = 0;
    return;
end

for j = 7:-1:(sns1+1)
    ns1(j) = NaN;
end

Alm.RecieptDay = ns1(1);
Alm.RecieptMonth = ns1(2);
Alm.RecieptYear = ns1(3);
Alm.RecieptSecond = ns1(4);
Alm.RecieptGPSWeek = ns1(5);
Alm.RecieptGPSSecond = ns1(6);


% Second string
% 12 12345  1234  123456   12 12 1234  12345.123  1.12345678E-12  1.12345678E-12  1.12345678E-12
% 1 - номер PRN (PRN)
% 2 - обобщенный признак здоровья (Health)
% 3 - неделя GPS (альманаха) (GPS Week)
% 4 - время недели GPS, с (альманаха) (GPS sec of week, s)
% 5 - число (Day)
% 6 - месяц (Month)
% 7 - год (Year)
% 8 - время альманаха, с (Almanac time, s)
% 9 - поправка времени КА GPS относительно системного времени, с (Time shift SV<>System, s)
% 10- скорость поправки времени КА GPS относительно системного времени, с/с (Rate of Time shift SV<>System, s/s)
% 11- Om0 - скорость долготы узла, полуциклы/c (Rate of longtitude of ascending node, half-of-cicle/s)
ns2 = str2num(s2);
sns2 = length(ns2);
if sns2 > 11
    return;
end

for j = 11:-1:(sns2+1)
    ns2(j) = NaN;
end

Alm.PRN = ns2(1);
Alm.Health = ns2(2);
Alm.GPSWeek = ns2(3);
Alm.GPSSec = ns2(4);
% Alm.Day = ns2(5);
% Alm.Month = ns2(6);
% Alm.Year = ns2(7);
% Alm.AlmTime = ns2(8);
Alm.t = CTime(ns2(8), ns2(5), ns2(6), ns2(7));
Alm.TimeShift = ns2(9);
Alm.dTimeShift = ns2(10);
Alm.dOmega = ns2(11)*pi; % [rad/s]


%Third string
% 1.12345678E-12  1.12345678E-12  1.12345678E-12  1.12345678E-12  1.12345678E-12  1.12345678E-12
% 1 - Om0 - долгота узла, полуциклы (Omega: longtitude of ascending node, half-of-cicle) 
% 2 - I - наклонение, полуциклы (i: Inclination, half-of-cicle)
% 3 - w - аргумент перигея, полуциклы (omega: argument of perigee, half-of-cicle)
% 4 - E - эксцентриситет (Eccenticity)
% 5 - SQRT(A) - корень из большой полуоси, м^0.5 (root of semi-major axis)
% 6 - M0 - средняя аномалия, полуциклы (mean anomaly, half-of-cicle)
ns3 = str2num(s3);
sns3 = length(ns3);
if sns3 > 6
    return;
end

for j = 6:-1:(sns3+1)
    ns3(j) = NaN;
end

Alm.Omega = ns3(1)*pi; % [s]
Alm.i = ns3(2)*pi; % [s]
Alm.omega = ns3(3)*pi; % [s]
Alm.e = ns3(4);
Alm.sqrta = ns3(5);
Alm.M0 = ns3(6)*pi; % [s]

ok = 1;
end




