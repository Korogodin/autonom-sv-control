function GLO_Const = LoadGLOConst( filename, t )
%LOADGPSCONST Load GLONASS constellation from almanac

fprintf('Load Almanac for GLONASS constellation...\n');

fid = fopen(filename);

N_GLO = 24;
AlmSVList = zeros(1, N_GLO);

j = 0;
while (~feof(fid))
    s1 = fgets(fid); 
    if (ischar(s1))
        s2 = fgets(fid); 
        if (ischar(s2))
            s3 = fgets(fid);
            if (ischar(s3))
                j = j + 1;
                [Alm(j) ok] = ParseAGL(s1, s2, s3);
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
                    error('Error of reading of almanac AGL file');
                end
            end
        end
    end
end
fclose(fid);

k = 0;
for (j = 1:N_GLO)
    if AlmSVList(j) > 0
        k = k + 1;
        GLO_Const(k) = SpaceVehicle_GLO(['GL' num2str(Alm(AlmSVList(j)).PRN)], Alm(AlmSVList(j)));
        fprintf('\tSV%02.0f - almanac %s - orbit calculation for %s ...', ...
                Alm(AlmSVList(j)).PRN, Alm(AlmSVList(j)).t.Stri, t.Stri);
        GLO_Const(k).calcOrbit(t); 
        fprintf('OK\n');
    end
end

fprintf('Constellation of %.0f GLONASS SVs was loaded\n', k);

end

function [Alm ok] = ParseAGL(s1, s2, s3)
%PARSEAGP Parser for AGL-files
%See: ftp://ftp.glonass-iac.ru/MCC/FORMAT/Format.agl
%       ICD 5.1, pages 64+

ok = 0;

% First string:
% 12 12 1234  12345  123456789101234567890
% 1 - число получения альманаха (Almanac reciept day)
% 2 - месяц получения альманаха (Almanac reciept month)
% 3 - год получения альманаха (Almanac reciept year)
% 4 - время получения альманаха от начала суток, с UTC (Time of almanac reciept from begining of day, sec UTC)
% 5 - комментарий (приемник с которого получено, версия SW и т.д.) (Comment)

ns1 = str2num(s1);
sns1 = length(ns1);
if sns1 > 5
    Alm = 0;
    return;
end

for j = 5:-1:(sns1+1)
    ns1(j) = NaN;
end

Alm.RecieptDay = ns1(1);
Alm.RecieptMonth = ns1(2);
Alm.RecieptYear = ns1(3);
Alm.RecieptSecond = ns1(4);


% Second string
% 12 12  1  12 12 1234  1.123456789E-12  1.123456789E-12  1.123456789E-12  1.123456789E-12

% 1 - номер КА в группировке (Number of SV in constellation)
% 2 - номер частотного слота (-7 - 24) (Number of freq slot (-7 - 24))
% 3 - признак здоровья по альманаху (0 - 1) (Health)
% 4 - число (Day)
% 5 - месяц (Month)
% 6 - год (Year)
% 7 - время прохождения первого узла, на которое все дано, с (Time of transit of first node, s (base point of time - t0))
% 8 - поправка ГЛОНАСС-UTC, с (Shift GLO-UTC, sec)
% 9 - поправка GPS-ГЛОНАСС, с (Shift GPS-GLO, sec)
% 10 - поправка времени КА ГЛОНАСС относительно системного времени, с (Shift SV-SystemTime, sec)

ns2 = str2num(s2);
sns2 = length(ns2);
if sns2 > 10
    return;
end

for j = 10:-1:(sns2+1)
    ns2(j) = NaN;
end

Alm.PRN = ns2(1); % [num]
Alm.FreqSlot = ns2(2); % [num]
Alm.Health = ns2(3);  % [bool]
Alm.t = CTime(ns2(7), ns2(4), ns2(5), ns2(6)); % [s] [d] [m] [y]
Alm.TS_GLOUTC = ns2(8); % [s]
Alm.TS_GPSGLO = ns2(9); % [s]
Alm.TS_SVGLO = ns2(10); % [s]


%Third string
% 1.1234567E-12  1.1234567E-12  1.1234567E-12  1.1234567E-12  1.1234567E-12  1.1234567E-12
% 1 - Lam - долгота узла, полуциклы (Omega: longtitude of ascending node, half-of-cicle) 
% 2 - dI - коррекция наклонения, полуциклы (delta i: Correction for inclination, half-of-cicle)
% 3 - w - аргумент перигея, полуциклы (omega: argument of perigee, half-of-cicle)
% 4 - E - эксцентриситет  (Eccentricity)
% 5 - dT - поправка к драконическому периоду, с (correction for draconian period, s)
% 6 - dTT - поправка к драконическому периоду, с/виток (correction for correction for draconian period, s/cicle)

ns3 = str2num(s3);
sns3 = length(ns3);
if sns3 > 6
    return;
end

for j = 6:-1:(sns3+1)
    ns3(j) = NaN;
end

Alm.Omega = ns3(1)*pi; % [s]
Alm.di = ns3(2)*pi; % [s]
Alm.omega = ns3(3)*pi; % [s]
Alm.e = ns3(4);
Alm.dT = ns3(5); % [s]
Alm.ddT = ns3(6); % [s/cicle]

ok = 1;
end




