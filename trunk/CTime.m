classdef CTime
    %CTIME Class of time scale
    
    properties
        Day
        Month
        Year
        DMY
        t
        Stri
    end
    
    methods
        function t = CTime(tk, Day, Month, Year)
            t.t = tk;
            t.Day = Day;
            t.Month = Month;
            t.Year = Year;
            t.DMY = datevec([num2str(Day) '.' num2str(Month) '.' num2str(Year)], 'dd.mm.yyyy');
            t.Stri = sprintf('%02.0f.%02.0f.%.0f', Day, Month, Year);
        end
        function res = minus(obj1, obj2)
            res = (datenum(obj1.DMY) - datenum(obj2.DMY))*60*60*24;
        end
    end
    
end

