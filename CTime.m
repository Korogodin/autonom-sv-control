classdef CTime
    %CTIME Class of time scale
    
    properties
        Day
        Month
        Year
        DMY
        t
        Stri
        JDN % Julian Day Number
        MJDN
        ST0 % Sidereal time (for 00:00) (mean)
        S0
    end
    
    methods
        function t = CTime(tk, Day, Month, Year)
            t.t = tk;
            t.Day = Day;
            t.Month = Month;
            t.Year = Year;
            t.DMY = datevec([num2str(Day) '.' num2str(Month) '.' num2str(Year)], 'dd.mm.yyyy');
            t.Stri = sprintf('%02.0f.%02.0f.%.0f', Day, Month, Year);
            
            % Julian
            a = fix((14 - t.Month)/12);
            y = t.Year + 4800 - a;
            m = t.Month + 12*a - 3;
            t.JDN = t.Day + fix((153*m + 2)/5) + 365*y + fix(y/4) - fix(y/100) + fix(y/400) - 32045;
            t.MJDN = t.JDN - 2400001;
            
            %Sidereal time
            A1 = 24110.54841;
            A2 = 8640184.812;
            A3 = 0.093104;
            A4 = 0.0000062;
            T0 = (t.MJDN - 51544.5) / 36525;
            t.ST0 = A1 + A2 * T0 + A3 * T0 ^ 2 - A4 * T0 ^ 3; % sideral time, s
            t.S0 = deg2rad(mod(t.ST0 *15 / 3600, 360)); % -> rad
        end
        function res = minus(obj1, obj2)
            res = (datenum(obj1.DMY) - datenum(obj2.DMY))*60*60*24;
        end
    end
    
end

