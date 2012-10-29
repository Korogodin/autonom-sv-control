classdef SpaceVehicle_CNS < SpaceVehicle
    %SPACEVEHICLE_CNS Class of sattellites-consumers of GNSS
        
    properties
        p
        a
        e
        Omega
        omega
        i
        M0
        t0
        
        % GDOP  (1:K)
        DOP_GPS
        DOP_GLO
        DOP_Comm
        DOP_GPS_noNaN
        DOP_GLO_noNaN
        DOP_Comm_noNaN

        % Min power in max 4 of satts, dBWt (1:K)
        MinIn4_GPS 
        MinIn4_GLO
        MinIn4_GPS_noNaN 
        MinIn4_GLO_noNaN 
        
        % Number of satellites with measurement
        Num_Sync_GLO % (1:K)
        Num_Sync_GPS

        
        
    end
    
    methods
        function SV = SpaceVehicle_CNS(Name, Type, a, e, Omega, omega, i, M0, t0)
            SV = SV@SpaceVehicle(Name, Type);
            SV.a = a;
            SV.e = e;
            SV.Omega = Omega;
            SV.omega = omega;
            SV.i = i;
            SV.M0 = M0;
            SV.t0 = t0; SV.t0.t = SV.t0.t(1);
        end
        
        function calcKeplerOrbit(SV, t)

            Nmod = length(t.t);
            global mu_earth
            
            e = SV.e;
            a = SV.a;
            p = a .* (1-e.^2);
            i = SV.i;
            Omega = SV.Omega;
            omega = SV.omega;
            M0 = SV.M0;
            
            % Mean anomaly
            if isempty(mu_earth)
                mu_earth = 3.9860044e14; % [m^3/s^2] Gravity constant
            end
            
            tk = t.t - SV.t0.t + (t - SV.t0);
            M = M0 + sqrt(mu_earth) ./ (sqrt(a).^3) .*tk;
            
            % Eccentricity anomaly
            E = nan(1,Nmod);
            options_solve = optimset('Display','off');
            E(1) = fsolve(@(E)(E-e*sin(E) - M(1)), M0, options_solve);
            for j  = 2:Nmod
                E(j) = fsolve(@(E)(E-e*sin(E)  - M(j)), E(j-1), options_solve);
            end
            
            % True anomaly
            theta = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
            for j = 2:Nmod
                while abs(theta(j) - theta(j-1)) > pi
                    if (theta(j) - theta(j-1)) > pi
                        theta(j) = theta(j) - 2*pi;
                    elseif (theta(j) - theta(j-1)) < -pi
                        theta(j) = theta(j) + 2*pi;
                    end
                end
            end
            
            munapi = sqrt(mu_earth ./ p);
            Vr = munapi.*e.*sin(theta); 
            Vu = munapi.*(1+e.*cos(theta));
            SV.r = p ./ (1+e.*cos(theta)); 
            u = theta + omega; 

            for j = 1:Nmod
                X = U3(-Omega)*U1(-i)*U3(-u(j))*[SV.r(j); 0; 0];
                SV.x(j) = X(1); SV.y(j) = X(2); SV.z(j) = X(3);
            end
            
            SV.Vx = Vr.*SV.x./SV.r - Vu.*(sin(u).*cos(Omega) + cos(u).*sin(Omega).*cos(i));
            SV.Vy = Vr.*SV.y./SV.r - Vu.*(sin(u).*sin(Omega) - cos(u).*cos(Omega).*cos(i));
            SV.Vz = Vr.*SV.z./SV.r + Vu.*cos(u).*sin(i);

        end
        
        function calcDOP(SV, t, GLO_Const, GPS_Const)
            % Execute after calculation of Redioline!
            
            Nmod = length(t.t);
            N_GLO = length(GLO_Const);
            N_GPS = length(GPS_Const);
            
            SV.DOP_GPS = nan(1, Nmod);
            SV.DOP_GLO= nan(1, Nmod);
            SV.DOP_Comm = nan(1, Nmod);

            SV.DOP_GPS_noNaN = nan(1, Nmod);
            SV.DOP_GLO_noNaN = nan(1, Nmod);
            SV.DOP_Comm_noNaN = nan(1, Nmod);
            

            globals;
            
            for k = 1:Nmod
                
                % Calc GDOP for GLONASS
                A = -ones(N_GLO, 4);
                jA = 0;
                for jg = 1:N_GLO
                    if (GLO_Const(jg).RL_L1.Power_noNaN(k) > SensLoopGLO)&&GLO_Const(jg).RL_L1.Sync(k) 
                        jA = jA + 1;
                        A(jA, 1:3) = GLO_Const(jg).RL_L1.dr' / GLO_Const(jg).RL_L1.R;
                    end
                end
                if jA >= 4
                    A = A(1:jA, 4);
                    SV.DOP_GLO(k) = sqrt(trace( (A'*A)^-1));
                    SV.DOP_GLO_noNaN(k) = SV.DOP_GLO(k);
                else
                    SV.DOP_GLO(k) = NaN;
                    SV.DOP_GLO_noNaN(k) = 80;
                end
                

                % Calc GDOP for GPS
                A = -ones(N_GPS, 4);
                jA = 0;
                for jg = 1:N_GPS
                    if (GPS_Const(jg).RL_L1.Power_noNaN(k) > SensLoopGPS)&&GPS_Const(jg).RL_L1.Sync(k)
                        jA = jA + 1;
                        A(jA, 1:3) = GPS_Const(jg).RL_L1.dr' / GPS_Const(jg).RL_L1.R;
                    end
                end
                if jA >= 4
                    A = A(1:jA, 4);
                    SV.DOP_GPS(k) = sqrt(trace( (A'*A)^-1));
                    SV.DOP_GPS_noNaN(k) = SV.DOP_GPS(k);
                else
                    SV.DOP_GPS(k) = NaN;
                    SV.DOP_GPS_noNaN(k) = 80;
                end          
                
                
                % Calc GDOP for GPS+GLO
                A = -ones(N_GPS+N_GLO, 4);
                A1 = -ones(N_GPS, 4);
                jA1 = 0;
                for jg = 1:N_GPS
                    if (GPS_Const(jg).RL_L1.Power_noNaN(k) > SensLoopGPS)&&GPS_Const(jg).RL_L1.Sync(k)
                        jA1 = jA1 + 1;
                        A1(jA, 1:3) = GPS_Const(jg).RL_L1.dr' / GPS_Const(jg).RL_L1.R;
                    end
                end
                A2 = -ones(N_GLO, 4);
                jA2 = 0;
                for jg = 1:N_GLO
                    if (GLO_Const(jg).RL_L1.Power_noNaN(k) > SensLoopGLO)&&GPS_Const(jg).RL_L1.Sync(k)
                        jA2 = jA2 + 1;
                        A(jA2, 1:3) = GLO_Const(jg).RL_L1.dr' / GLO_Const(jg).RL_L1.R;
                    end
                end         
                
                jA = 0;
                if jA1 >= 2
                    A(jA+1:jA+jA1) = A1(1:jA1, 4);
                    jA = jA + jA1;
                end
                if jA2 >= 2
                    A(jA+1:jA+jA2) = A2(1:jA2, 4);
                    jA = jA + jA2;
                end
                    
                 if (jA >= 5)||(jA1 >= 4)||(jA2 >= 4)
                    A = A(1:jA, 4);
                    SV.DOP_Comm(k) = sqrt(trace( (A'*A)^-1));
                    SV.DOP_Comm_noNaN(k) = SV.DOP_Comm(k);
                else
                    SV.DOP_Comm(k) = NaN;
                    SV.DOP_Comm_noNaN(k) = 80;
                end 
                
            end
                
        end
        
        function calcSatNums(SV, t, GLO_Const, GPS_Const)
            % Execute after calculation of Redioline!
            
            Nmod = length(t.t);
            N_GLO = length(GLO_Const);
            N_GPS = length(GPS_Const);
            
            SV.MinIn4_GPS = nan(1, Nmod);
            SV.MinIn4_GLO= nan(1, Nmod);            
            SV.MinIn4_GPS_noNaN = nan(1, Nmod);
            SV.MinIn4_GLO_noNaN = nan(1, Nmod);            
            
            % Number of satellites with measurements
            SV.Num_Sync_GLO = zeros(1, Nmod);
            SV.Num_Sync_GPS = zeros(1, Nmod);
            for k = 1:Nmod
                
                % For GLONASS
                PowMatr = nan(1, N_GLO);
                for jg = 1:N_GLO
                    if GLO_Const(jg).RL_L1.Sync(k) 
                        PowMatr(jg) = GLO_Const(jg).RL_L1.Power_noNaN(k);
                        SV.Num_Sync_GLO(k) = SV.Num_Sync_GLO(k) + 1;
                    else
                        PowMatr(jg) = -250;
                    end
                    
                end
                if N_GLO >= 4
                    for jj = 1:4
                        [val num] = max(PowMatr);
                        PowMatr(num) = -1000;
                    end
                    SV.MinIn4_GLO_noNaN(k) = val;
                    if val <= -250
                        SV.MinIn4_GLO(k) = NaN;
                    else
                        SV.MinIn4_GLO(k) = val;
                    end
                else
                    SV.MinIn4_GLO(k) = NaN;
                    SV.MinIn4_GLO_noNaN(k) = -250;
                end
                    
                % For GPS
                PowMatr = nan(1, N_GPS);
                for jg = 1:N_GPS
                    if GPS_Const(jg).RL_L1.Sync(k) 
                        PowMatr(jg) = GPS_Const(jg).RL_L1.Power_noNaN(k);
                        SV.Num_Sync_GPS(k) = SV.Num_Sync_GPS(k) + 1;
                    else
                        PowMatr(jg) = -250;
                    end                    
                end
                if N_GPS >= 4
                    for jj = 1:4
                        [val num] = max(PowMatr);
                        PowMatr(num) = -1000;
                    end
                    SV.MinIn4_GPS_noNaN(k) = val;
                    if val <= -250
                        SV.MinIn4_GPS(k) = NaN;
                    else
                        SV.MinIn4_GPS(k) = val;
                    end
                else
                    SV.MinIn4_GPS(k) = NaN;
                    SV.MinIn4_GPS_noNaN(k) = -250;
                end
                
            end
        end
        
    end
    
    
    
end

