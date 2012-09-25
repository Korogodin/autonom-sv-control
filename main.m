try 
    close(handle_fig_main); % Close old output form
end

close all
clear 
clc

globals;
t = 0:dTmod:Tmod;
Nmod = length(t);

GPS_Const = CreateGPS();