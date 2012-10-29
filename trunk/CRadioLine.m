classdef CRadioLine < handle
    %CRADIOLINE Radioline between 2 satellites
    
    properties
        % SV1 - transmitter
        % SV2 - receiver
        
        R % length (k)
        dr % radius-vector SV1 -> SV2 = transmitter->receiver (k)
        ThroughEarth  % flag of Earth (1:K)
        lambda % length of wave
        TransPower % Power of transmitter, dBWt
        
        Gr  % Radiation pattern for receiver's antenna
        Gt % Radiation pattern of transmitter's antenna
        
        Power % power on receiver's antenna output (1:K)
        Power_noNaN % NaN -> -200 dBWt (1:K)
        
        SensSync % Sensivity, dBWt
        SensLoop
        
        Sync % Syncronization flag
        
        alpha % angle between SV1 antenna axe and vision line
        beta % angle between SV2 antenna axe and vision line
        
        G_rec % Gain in receiver antenna (1:K)
        G_tr % Gain in receiver antenna (1:K)
    end
    
    methods
        function RL = CRadioLine(lambda, Gr, Gt, TransPower, SensLoop, SensSync) % Constructor
            RL.lambda = lambda; % length of wave
            RL.Gr = Gr;  % Radiation pattern for receiver's antenna
            RL.Gt = Gt; % Radiation pattern of transmitter's antenna
            RL.TransPower = TransPower; % Power of transmitter, dBWt
            RL.SensLoop = SensLoop;
            RL.SensSync = SensSync;
        end
        
        function CalcRL(RL, t, SV1, SV2) % Calculate radioline
            % SV1 - must be transmitter
            % SV2 - must be receiver
            Nmod = length(t.t);

            % Memory allocation
            RL.ThroughEarth = nan(1, Nmod);
            RL.Power = nan(1, Nmod);
            RL.Power_noNaN = nan(1, Nmod);
            RL.alpha = nan(1, Nmod);
            RL.beta = nan(1, Nmod);
            RL.G_rec = nan(1, Nmod);
            RL.G_tr = nan(1, Nmod);     
            
            RL.Sync = zeros(1, Nmod);
            
            for k = 1:Nmod
                SV1.touchX(k); % set X for k
                SV2.touchX(k);
                RL.dr = SV1.dr(SV2);
                RL.R = SV1.dist(SV2);
                
                if RL.thrEarth(SV2, k) % if throuh planet
                    RL.Power(k) = NaN; % no power
                    RL.Power_noNaN(k) = -250;
                    RL.Sync(k) = 0;
                else
                    SV1.CalcAntM(k); % Calculate direction of antenna
                    SV2.CalcAntM(k);
                    
                    RL.Angles(SV1, SV2, k); % Calc angles between AntM and vision line
                    RL.ApproxGain(k); % Calculate gain for receiver's and transmitter's antennas
                    
                    RL.Power(k) = RL.TransPower + RL.G_tr(k) + 20*log10(RL.lambda / (4*pi*RL.R)) - 1 + RL.G_rec(k);
                    
                    if k > 1
                        RL.Sync(k) = RL.Sync(k - 1);
                    else
                        RL.Sync(1) = 0;
                    end
                    if isnan(RL.Power(k))
                        RL.Power_noNaN(k) = -250;
                        RL.Sync(k) = 0;
                    else
                        RL.Power_noNaN(k) = RL.Power(k);
                        if RL.Power(k) < RL.SensLoop
                            RL.Sync(k) = 0;
                        elseif RL.Power(k) >= RL.SensSync
                            RL.Sync(k) = 1;
                        end
                    end
                end
                
            end
        end
        
        function res = thrEarth(RL, SV2, k) % Pass a planet?
            % SV2.X, SV2.r(k), RL.dr, RL.R must be actual
            R_e = 6371e3;  

            gamma = acos(-RL.dr' * (-SV2.X) / RL.R / SV2.r(k));
            if (sin(gamma) * SV2.r(k) < R_e) && (gamma < pi/2)
                res = 1;
                RL.ThroughEarth(k) = 1;
                return;
            end
            res = 0;
            RL.ThroughEarth(k) = 0;
        end    
        
        function Angles(RL, SV1, SV2, k)
            % Calculate angle between antenna axe and vision line
            % SV1.AntM, SV2.AntM, RL.dr, RL.R must be actual
            
            if RL.R > 0
                dr1 = RL.dr / RL.R;
                RL.alpha(k) = rad2deg(acos( SV1.AntM' * dr1 ));
                RL.beta(k) = rad2deg(acos( SV2.AntM' * (-dr1) ));
            else
                disp('SVs are in one point of space');
            end
            return;
        end        
        
        function ApproxGain(RL, k)
            % Approx gain of antenna
            % RL.alpha(k), RL.beta(k) must be actual
            
            cA = ceil(RL.alpha(k));
            fA = fix(RL.alpha(k));
            dA = RL.alpha(k) - fA;
            RL.G_tr(k) = (RL.Gt(cA+1) - RL.Gt(fA+1))*dA + RL.Gt(fA+1);
          
            cB = ceil(RL.beta(k));
            fB = fix(RL.beta(k));
            dB = RL.beta(k) - fB;
            RL.G_rec(k) = (RL.Gr(cB+1) - RL.Gr(fB+1))*dB + RL.Gr(fB+1);

        end
    end
    
end

