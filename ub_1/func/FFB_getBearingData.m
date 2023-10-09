function bearingDataMatrix = FFB_getBearingData(bearingData)

% getBearingData function interpolates suitable table containing 
% dimensionless stiffness and damping coefficients with respect to 
% Sommerfeld number. Data tables are from DIN 31657-2 (Multilobe bearings)
% and DIN 31657-3 (Tilting pad bearings).
% 
% Syntax: bearingDataMatrix = getBearingData(bearingData)
% 
% where bearingData.Req.DIN             2 for multilobe bearings and 3 for
%                                       tilting pad bearings
%       bearingData.Req.Z               Number of lobes/pads
%       bearingData.Req.B_D             Geometric ratio of bearings leangth and
%                                       diameter
%       bearingData.Req.Omega           Lobe/pad angle
%       bearingData.Req.fiiP1           Force direction with respect to lobe/pad
%       bearingData.Req.h0max           Gap height ratio (only multilobe)
%       bearingData.Req.OmegaF_Omega    Pivot offset (only tilting pad)
%       bearingData.Req.dRB_CR          Profiling (only tilting pad)
% 
% bearingDataMatrix is n x 9 matrix containing dimensionless stiffness and 
% damping coefficients with respect to Sommerfeld number.
% 
% Here are basic structures as in parts 2 and 3 of DIN 31657
% 
% Principle here is to define can DIN 31657 table be directly used or is
% there need to interpolate with respect to desired characteristics from
% between different tables.

% Redefining geometric ratio to shorten some if estatements
B_D=bearingData.Req.B_D;

% Formatting variable for possible force direction (Omega) interpolation.
% NOTE: Only DIN3 (Tinting pad) tables provide possibility to interpolate
% with respect to force direction (Omega).
OmegaInterpolation=0;


%% DIN2 interpolation check (Multilobes)
if bearingData.Req.DIN == 2
    
    
    %% 2 lobes (Z) and lobe angle 150 deg (Omega) check
    if bearingData.Req.Z == 2 && bearingData.Req.Omega == 150
        
        
        %% 2 lobes, lobe angle 150 deg (Omega) and force direction 180 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 180
            
            % If requested h0max value is not supported in standard it is
            % being interpolated.
            % .Req.interp struct defines which characteristic value is
            % needed to be interpolated, in this case h0max if it is not
            % supported in standard.
            if sum(bearingData.Req.h0max == [3 5]) ~= 1
                bearingData.Req.interp.h0max=[3 5];
            end
            
            % From standard can be seen if Z=2 and Omega=180 that it
            % provides only one geometric ratio (B_D) so user defined .Req
            % struct is being manipulated to one provided in standard
            % regardless what h0max is.
            bearingData.Req.B_D=0.75;
            
        end % end of 2 lobes (Z), lobe angle 150 deg (Omega) and force 
        % direction 180 deg (fiiP1) check
        
        
        %% 2 lobes, lobe angle 150 deg (Omega) and force direction 240 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 240
            
            % If requested h0max value is not supported in standard it is
            % being interpolated.
            % .Req.interp struct defines which characteristic value is
            % needed to be interpolated, in this case h0max if it is not
            % supported in standard.
            if sum(bearingData.Req.h0max == [3 5]) ~= 1
                bearingData.Req.interp.h0max=[3 5];
            end
            
            % From standard can be seen if Z=2 and Omega=180 that it
            % provides only one geometric ratio (B_D) so user defined .Req
            % struct is being manipulated to one provided in standard
            % regardless what h0max is.
            bearingData.Req.B_D=0.75;
            
        end % end of 2 lobes (Z), lobe angle 150 deg (Omega) and force 
        % direction 240 deg (fiiP1) check
        
        
        %% 2 lobes (Z), lobe angle 150 deg (Omega) and force direction 270 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 270
            
            % User defined h0max value is being check if it is directly
            % usable from standard or is it necessary to interpolate.
            % Now standard provides three (3) different values and if
            % interpolation is needed here is defined closest values for
            % interpolation.
            if sum(bearingData.Req.h0max == [1 3 5]) ~= 1
                if bearingData.Req.h0max < 3
                    
                    % Generally speaking valiables of .Req.interp struct 
                    % over defines always user defined values from .Req.
                    % This is because .Req.interp expresses needed
                    % interpolation steps. This over defining is done in
                    % searchPerms function.
                    bearingData.Req.interp.h0max=[1 3];
                else
                    bearingData.Req.interp.h0max=[3 5];
                end
            end
            
            % In this case (2 lobes (Z), lobe angle 150 deg (Omega) and 
            % force direction 270 deg (fiiP1)) standard provides also
            % interpolation with respect to geometric ratio (B_D) and these
            % h0max and B_D interpolations can work together because
            % standard provides enough tables, nine (3*3=9) altogether.
            % Here is also defined closest tables for interpolation if
            % geometric ratio is not supported directly in standard.
            if sum(B_D == [0.5 0.75 1]) ~= 1
                if B_D < 0.75
                    bearingData.Req.interp.B_D=[0.5 0.75];
                else
                    bearingData.Req.interp.B_D=[0.75 1];
                end
            end
            
        end % end of 2 lobes (Z), lobe angle 150 deg (Omega) and force 
        % direction 270 deg (fiiP1) check
        
        
        %% 2 lobes (Z), lobe angle 150 deg (Omega) and force direction 300 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 300
            
            % If requested h0max value is not supported in standard it is
            % being interpolated.
            % .Req.interp struct defines which characteristic value is
            % needed to be interpolated, in this case h0max if it is not
            % supported in standard.
            if sum(bearingData.Req.h0max == [3 5]) ~= 1
                bearingData.Req.interp.h0max=[3 5];
            end
            
            % From standard can be seen if Z=2 and Omega=180 that it
            % provides only one geometric ratio (B_D) so user defined .Req
            % struct is being manipulated to one provided in standard
            % regardless what h0max is.
            bearingData.Req.B_D=0.75;
            
        end % end of 2 lobes (Z), lobe angle 150 deg (Omega) and force 
        % direction 300 deg (fiiP1) check
        
        
    %% 3 lobes (Z) and lobe angle 100 deg (Omega) check
    elseif bearingData.Req.Z == 3 && bearingData.Req.Omega == 100
        
        
        %% 3 lobes (Z), lobe angle 100 deg (Omega) and force direction 240 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 240
            
            % If requested h0max value is not supported in standard it is
            % being interpolated.
            % .Req.interp struct defines which characteristic value is
            % needed to be interpolated, in this case h0max if it is not
            % supported in standard.
            if sum(bearingData.Req.h0max == [3 5]) ~= 1
                bearingData.Req.interp.h0max=[3 5];
            end
            
            % From standard can be seen if Z=2 and Omega=180 that it
            % provides only one geometric ratio (B_D) so user defined .Req
            % struct is being manipulated to one provided in standard
            % regardless what h0max is.
            bearingData.Req.B_D=0.75;
            
        end % end of 3 lobes (Z), lobe angle 100 deg (Omega) and force 
        % direction 240 deg (fiiP1) check
        
        
        %% 3 lobes (Z), lobe angle 100 deg (Omega) and force direction 300 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 300
            
            % User defined h0max value is being check if it is dirextly
            % usable from standard or is it necessary to interpolate.
            % Now standard provides three (3) different values and if
            % interpolation is needed here is defined closest values for
            % interpolation.
            if sum(bearingData.Req.h0max == [1 3 5]) ~= 1
                if bearingData.Req.h0max < 3
                    
                    % Generally speaking valiables of .Req.interp struct 
                    % over defines always user defined values from .Req.
                    % This is because .Req.interp expresses needed
                    % interpolation steps. This over defining is done in
                    % searchPerms function.
                    bearingData.Req.interp.h0max=[1 3];
                else
                    bearingData.Req.interp.h0max=[3 5];
                end
            end
            
            % In this case (3 lobes (Z), lobe angle 100 deg (Omega) and 
            % force direction 300 deg (fiiP1)) standard provides also
            % interpolation with respect to geometric ratio (B_D) and these
            % h0max and B_D interpolations can work together because
            % standard provides enough tables, nine (3*3=9) altogether.
            % Here is also defined closest tables for interpolation if
            % geometric ratio is not supported directly in standard.
            if sum(B_D == [0.5 0.75 1]) ~= 1
                if B_D < 0.75
                    bearingData.Req.interp.B_D=[0.5 0.75];
                else
                    bearingData.Req.interp.B_D=[0.75 1];
                end
            end
            
        end % end of 3 lobes (Z), lobe angle 100 deg (Omega) and force 
        % direction 300 deg (fiiP1) check
        
        
    %% 4 lobes (Z) and lobe angle 70 deg (Omega) check
    elseif bearingData.Req.Z == 4 && bearingData.Req.Omega == 70
        
        
        %% 4 lobes (Z), lobe angle 70 deg (Omega) and force direction 270 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 270
            
            % If requested h0max value is not supported in standard it is
            % being interpolated.
            % .Req.interp struct defines which characteristic value is
            % needed to be interpolated, in this case h0max if it is not
            % supported in standard.
            if sum(bearingData.Req.h0max == [3 5]) ~= 1
                bearingData.Req.interp.h0max=[3 5];
            end
            
            % From standard can be seen if Z=2 and Omega=180 that it
            % provides only one geometric ratio (B_D) so user defined .Req
            % struct is being manipulated to one provided in standard
            % regardless what h0max is.
            bearingData.Req.B_D=0.75;
            
        end % end of 4 lobes (Z), lobe angle 70 deg (Omega) and force 
        % direction 270 deg (fiiP1) check
        
        
        %% 4 lobes (Z), lobe angle 70 deg (Omega) and force direction 270 deg (fiiP1) check
        if bearingData.Req.fiiP1 == 315
            
            % User defined h0max value is being check if it is dirextly
            % useable from standard or is it necessary to interpolate.
            % Now standard provides three (5) different values and if
            % interpolation is needed here is defined closest values for
            % interpolation.
            if sum(bearingData.Req.h0max == [1 2 3 4 5]) ~= 1
                if bearingData.Req.h0max < 2
                    
                    % Generally speaking valiables of .Req.interp struct 
                    % over defines always user defined values from .Req.
                    % This is because .Req.interp expresses needed
                    % interpolation steps. This over defining is done in
                    % searchPerms function.
                    bearingData.Req.interp.h0max=[1 2];
                elseif bearingData.Req.h0max < 3
                    bearingData.Req.interp.h0max=[2 3];
                elseif bearingData.Req.h0max < 4
                    bearingData.Req.interp.h0max=[3 4];
                else
                    bearingData.Req.interp.h0max=[4 5];
                end
            end
            
            % In this case (4 lobes (Z), lobe angle 70 deg (Omega) and 
            % force direction 315 deg (fiiP1)) standard provides also
            % interpolation with respect to geometric ratio (B_D) and these
            % h0max and B_D interpolations can work together because
            % standard provides enough tables, 15 (5*3=15) altogether.
            % Here is also defined closest tables for interpolation if
            % geometric ratio is not supported directly in standard.
            if sum(B_D == [0.5 0.75 1]) ~= 1
                if B_D < 0.75
                    bearingData.Req.interp.B_D=[0.5 0.75];
                else
                    bearingData.Req.interp.B_D=[0.75 1];
                end
            end
            
        end % end of 4 lobes (Z), lobe angle 70 deg (Omega) and force 
        % direction 270 deg (fiiP1) check
        
    end % end of 2, 3 and 4 lobes (Z) check
    
end % end of DIN2 interpolation check (Multilobes)


%% DIN3 interpolation check (Tilting pad)
% Here is also available pad angle (Omega) interpolation. Notice also every
% table in DIN 31657-3 is for 0.5 and 0.6 values of pivot offset 
% (OmegaF_Omega). Check more information from manual.
if bearingData.Req.DIN == 3
    
    
    %% Four pads (Z) check
    if bearingData.Req.Z == 4
        
        
        %% Four pads (Z) and pad angle (Omega) interpolation check
        if sum(bearingData.Req.Omega == [80 60]) ~= 1
            
            % If user defined pad angle (Omega) is not one provided
            % in stadard, the pad angle (Omega) interpolation is
            % activated.
            OmegaInterpolation=1;
            
            % Storing original user defined Omega value for further use
            bearingData.Req.Omega_orig=bearingData.Req.Omega;
            
            % Making .Req.interp struct for pad angle (Omega)
            % interpolation.
            bearingData.Req.interp.Omega=[80 60];
            
            
            %% Four pads (Z), pad angle (Omega) interpolation when Omega = 80 and force direction 45 deg (fiiP1) check
            if bearingData.Req.fiiP1 == 45
                
                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end
                
                % In this case (4 pads (Z), pad angle 80 deg (Omega) and 
                % force direction 45 deg (fiiP1)) standard provides also
                % interpolation with respect to geometric ratio (B_D) and these
                % dRB_CR and B_D interpolations can work together because
                % standard provides enough tables, nine (3*3=9) altogether.
                % Here is also defined closest tables for interpolation if
                % geometric ratio is not supported directly in standard.
                if sum(B_D == [0.5 0.75 1]) ~= 1
                    if B_D < 0.75
                        bearingData.Req.interp.B_D=[0.5 0.75];
                    else
                        bearingData.Req.interp.B_D=[0.75 1];
                    end
                end
                
                
            %% Four pads (Z), pad angle (Omega) interpolation when Omega = 80 and force direction 0 deg (fiiP1) check
            elseif bearingData.Req.fiiP1 == 0
                
                % If force direction (fiiP1) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.75;
                
            end % end of Four pads (Z), pad angle (Omega) interpolation 
            % when Omega = 80 and force direction 45 deg (fiiP1) check
            
            
            %% Four pads (Z), pad angle (Omega) interpolation first step when Omega = 80
            % Here is done first interpolation step. Pad angle (Omega)
            % interpolation needs two ordinary interpolation steps and
            % third interpolation step provides final interpolated table.
            bearingData.Req.Omega=80;
            
            % First interpolation step results
            bearingDataMatrix1 = FFB_OmegaInterp(bearingData);
            
            
            %% Four pads (Z), pad angle (Omega) interpolation when Omega = 60 and force direction 45 deg (fiiP1) check
            if bearingData.Req.fiiP1 == 45
                
                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end
                
                % In this case (4 pads (Z), pad angle 80 deg (Omega) and 
                % force direction 45 deg (fiiP1)) standard provides also
                % interpolation with respect to geometric ratio (B_D) and these
                % dRB_CR and B_D interpolations can work together because
                % standard provides enough tables, six (3*2=6) altogether.
                % Here is also defined closest tables for interpolation if
                % geometric ratio is not supported directly in standard.
                if sum(B_D == [0.5 0.75]) ~= 1
                    bearingData.Req.interp.B_D=[0.5 0.75];
                end
                
                
            %% Four pads (Z), pad angle (Omega) interpolation when Omega = 60 and force direction 0 deg (fiiP1) check
            elseif bearingData.Req.fiiP1 == 0
                
                % If force direction (fiiP1) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.5;
                
            end % end of Four pads (Z), pad angle (Omega) interpolation 
            % when Omega = 60 and force direction 45 deg (fiiP1) check
            
            
            %% Four pads (Z), pad angle (Omega) interpolation second step when Omega = 80
            % Here is done second interpolation step. Pad angle (Omega)
            % interpolation needs two ordinary interpolation steps and
            % third interpolation step provides final interpolated table.
            bearingData.Req.Omega=60;
            
            % Second interpolation step results
            bearingDataMatrix2 = FFB_OmegaInterp(bearingData);
            
        end % end of Four pads (Z) and pad angle (Omega) interpolation check
        
        
        %% Four pads (Z), pad angle (Omega) 80 check
        if bearingData.Req.Omega == 80
            
            
            %% Four pads (Z), pad angle (Omega) 80 and force direction 45 deg check
            if bearingData.Req.fiiP1 == 45
                
                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end
                
                % In this case (4 pads (Z), pad angle 80 deg (Omega) and 
                % force direction 45 deg (fiiP1)) standard provides also
                % interpolation with respect to geometric ratio (B_D) and these
                % dRB_CR and B_D interpolations can work together because
                % standard provides enough tables, nine (3*3=9) altogether.
                % Here is also defined closest tables for interpolation if
                % geometric ratio is not supported directly in standard.
                if sum(B_D == [0.5 0.75 1]) ~= 1
                    if B_D < 0.75
                        bearingData.Req.interp.B_D=[0.5 0.75];
                    else
                        bearingData.Req.interp.B_D=[0.75 1];
                    end
                end
            
            
            %% Four pads (Z), pad angle (Omega) 80 and force direction 0 deg check
            elseif bearingData.Req.fiiP1 == 0
                
                % If force direction (Omega) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.75;
                
            end % end of Four pads (Z), pad angle (Omega) 80 and force 
            % direction 0 and 45 deg check
            
        end % end of Four pads (Z), pad angle (Omega) 80 check
        
        
        %% Four pads (Z), pad angle (Omega) 60 check
        if bearingData.Req.Omega == 60
            
            %% Four pads (Z), pad angle (Omega) 60 and force direction 45 deg check
            if bearingData.Req.fiiP1 == 45
                
                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end
                
                % In this case (4 pads (Z), pad angle 60 deg (Omega) and 
                % force direction 45 deg (fiiP1)) standard provides also
                % interpolation with respect to geometric ratio (B_D) and these
                % dRB_CR and B_D interpolations can work together because
                % standard provides enough tables, six (3*2=6) altogether.
                % Here is also defined closest tables for interpolation if
                % geometric ratio is not supported directly in standard.
                if sum(B_D == [0.5 0.75]) ~= 1
                    bearingData.Req.interp.B_D=[0.5 0.75];
                end
            
            
            %% Four pads (Z), pad angle (Omega) 60 and force direction 0 deg check
            elseif bearingData.Req.fiiP1 == 0
                
                % If force direction (Omega) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.5;
                
            end % end of Four pads (Z), pad angle (Omega) 60 and force 
            % direction 0 and 45 deg check
            
        end % end of Four pads (Z), pad angle (Omega) 60 check
        
    end % end of Four pads (Z) check
    
    
    %% Five pads (Z) check
    if bearingData.Req.Z == 5
        
        
        %% Five pads (Z) and pad angle interpolation (Omega) check
        if sum(bearingData.Req.Omega == [60 45]) ~= 1
            
            % If user defined pad angle (Omega) is not one provided
            % in standard, the pad angle (Omega) interpolation is
            % activated.
            OmegaInterpolation=1;
            
            % Storing original user defined Omega value for further use
            bearingData.Req.Omega_orig=bearingData.Req.Omega;
            
            % Making .Req.interp struct for pad angle (Omega)
            % interpolation.
            bearingData.Req.interp.Omega=[60 45];
            
            
            %% Five pads (Z), pad angle (Omega) interpolation when Omega = 60 and force direction 36 deg (fiiP1) check
            if bearingData.Req.fiiP1 == 36
                
                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end
                
                % In this case (5 pads (Z), pad angle 60 deg (Omega) and 
                % force direction 36 deg (fiiP1)) standard provides also
                % interpolation with respect to geometric ratio (B_D) and these
                % dRB_CR and B_D interpolations can work together because
                % standard provides enough tables, six (3*2=6) altogether.
                % Here is also defined closest tables for interpolation if
                % geometric ratio is not supported directly in standard.
                if sum(B_D == [0.5 0.75]) ~= 1
                    bearingData.Req.interp.B_D=[0.5 0.75];
                end
            
                
            %% Five pads (Z), pad angle (Omega) interpolation when Omega = 60 and force direction 0 deg (fiiP1) check
            elseif bearingData.Req.fiiP1 == 0
                
                % If force direction (fiiP1) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.5;
                
            end % end of Five pads (Z), pad angle (Omega) interpolation 
            % when Omega = 60 and force direction 36 deg (fiiP1) check
           
            
            %% Five pads (Z), pad angle (Omega) interpolation first step when Omega = 60
            % Here is done first interpolation step. Pad angle (Omega)
            % interpolation needs two ordinary interpolation steps and
            % third interpolation step provides final interpolated table.
            bearingData.Req.Omega=60;
            
            % First interpolation step results
            bearingDataMatrix1 = FFB_OmegaInterp(bearingData);
            
            
            %% Five pads (Z), pad angle (Omega) interpolation when Omega = 45 and force direction 36 deg (fiiP1) check
            if bearingData.Req.fiiP1 == 36

                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end

                % From standard can be seen if Z=5 and Omega=45 that it
                % provides only one geometric ratio (B_D) so user defined .Req
                % struct is being manipulated to one provided in standard
                % regardless what dRB_CR is.
                bearingData.Req.B_D=0.5;
                
                
            %% Five pads (Z), pad angle (Omega) interpolation when Omega = 45 and force direction 0 deg (fiiP1) check
            elseif bearingData.Req.fiiP1 == 0

                % If force direction (fiiP1) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.5;

            end % end of Five pads (Z), pad angle (Omega) interpolation 
            % when Omega = 45 and force direction 0 and 36 deg (fiiP1) check
            
            
            %% Five pads (Z), pad angle (Omega) interpolation second step when Omega = 45
            % Here is done second interpolation step. Pad angle (Omega)
            % interpolation needs two ordinary interpolation steps and
            % third interpolation step provides final interpolated table.
            bearingData.Req.Omega=45;
            
            % Second interpolation step results
            bearingDataMatrix2 = FFB_OmegaInterp(bearingData);
            
            
        end % end of Five pads (Z) and pad angle interpolation (Omega) check
        
        
        %% Five pads (Z), pad angle (Omega) 60 check
        if bearingData.Req.Omega == 60
            
            
            %% Five pads (Z), pad angle (Omega) 60 and force direction 36 deg check
            if bearingData.Req.fiiP1 == 36
                
                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end
                
                % In this case (4 pads (Z), pad angle 60 deg (Omega) and 
                % force direction 36 deg (fiiP1)) standard provides also
                % interpolation with respect to geometric ratio (B_D) and these
                % dRB_CR and B_D interpolations can work together because
                % standard provides enough tables, six (3*2=6) altogether.
                % Here is also defined closest tables for interpolation if
                % geometric ratio is not supported directly in standard.
                if sum(B_D == [0.5 0.75]) ~= 1
                    bearingData.Req.interp.B_D=[0.5 0.75];
                end
                
                
            %% Five pads (Z), pad angle (Omega) 60 and force direction 0 deg check
            elseif bearingData.Req.fiiP1 == 0
                
                % If force direction (Omega) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.5;
                
            end % end of Five pads (Z), pad angle (Omega) 60 and force 
            % direction 0 and 36 deg check
            
        end % end of Five pads (Z), pad angle (Omega) 60 check
        
        
        %% Five pads (Z), pad angle (Omega) 45 check
        if bearingData.Req.Omega == 45
            
            
            %% Five pads (Z), pad angle (Omega) 45 and force direction 36 deg check
            if bearingData.Req.fiiP1 == 36
                
                % If requested dRB_CR value is not supported in standard it is
                % being interpolated.
                % .Req.interp struct defines which characteristic value is
                % needed to be interpolated, in this case dRB_CR if it is not
                % supported in standard.
                if sum(bearingData.Req.dRB_CR == [2 3 5]) ~= 1
                    if bearingData.Req.dRB_CR < 3
                        bearingData.Req.interp.dRB_CR=[2 3];
                    else
                        bearingData.Req.interp.dRB_CR=[3 5];
                    end
                end
                
                % From standard can be seen if Z=5 and Omega=45 that it
                % provides only one geometric ratio (B_D) so user defined .Req
                % struct is being manipulated to one provided in standard
                % regardless what dRB_CR is.
                bearingData.Req.B_D=0.5;
                
                
            %% Five pads (Z), pad angle (Omega) 45 and force direction 0 deg check
            elseif bearingData.Req.fiiP1 == 0
                
                % If force direction (Omega) 0 deg then standard provides only
                % one table and user defined values for dRB_CR and B_D are
                % being replaced with values provided in before mentioned
                % table.
                % 
                % NOTICE that this 0 deg force direction may cause some
                % error in results because of this automatic replacement of
                % user defined values.
                bearingData.Req.dRB_CR=3;
                bearingData.Req.B_D=0.5;
                
            end % end of Five pads (Z), pad angle (Omega) 45 and force 
            % direction 0 and 36 deg check
            
        end % end of Five pads (Z), pad angle (Omega) 45 check
        
    end % end of Five pads (Z) check
    
end % end of DIN3 interpolation check (Tilting pad)


%% Bearing data interpolation with pad angle (Omega) interpolation
if OmegaInterpolation == 1
    
    % NOTICE here is returned the original user defined pad angle (Omega)
    % value
    bearingData.Req.Omega=bearingData.Req.Omega_orig;
    
    % Here is done last pad angle (Omega) interpolation step meaning two
    % before provided tables are being here interpolated into one table.
    bearingDataMatrix = FFB_OmegaInterpCombine...
        (bearingData,bearingDataMatrix1,bearingDataMatrix2);

    
%% Bearing data interpolation without pad angle (Omega) interpolation
else
    
    % defining vectors for searching needed tables in interpolation
    combinations = FFB_searchPerms(bearingData);
    
    % here is done search for needed tables
    tables = FFB_searchTables(combinations);
    
    
    % defining interpolation vectors
    [tables,wc,wci] = FFB_createInterpVectors(combinations,tables,bearingData);
    
    
    % actual interpolation
    bearingDataMatrix = FFB_dataInterpolation(tables,wc,wci,bearingData.Req.DIN);
    
end

end



function bearingDataMatrix = FFB_OmegaInterp(bearingData)

% Making initial interpolation with respect to h0max/dRB_CR and B_D values
% for pad angle (Omega interpolation)

% defining vectors for searching needed tables in interpolation
combinations = FFB_searchPerms(bearingData);

% here is done search for needed tables
tables = FFB_searchTables(combinations);

% defining interpolation vectors
[tables,wc,wci] = FFB_createInterpVectors(combinations,tables,bearingData);

% actual interpolation
bearingDataMatrix = FFB_dataInterpolation(tables,wc,wci,bearingData.Req.DIN);

end



function bearingDataMatrix = FFB_OmegaInterpCombine...
    (bearingData,bearingDataMatrix1,bearingDataMatrix2)

% Combining two earlier interpolated tables by interpolating with respect to
% pad angle (Omega) hereby creating final table

% Interpolation values as bearingDataMatrix1 and bearingDataMatrix2
% matrices represents, also called as weight coefficients for interpolation
data.wc=bearingData.Req.interp.Omega;

% Original user defined pad angle (Omega) as interpolation coefficient,
% also desired weight coefficient
data.wci=bearingData.Req.Omega;

% Converting matrices into struct
data1.So=bearingDataMatrix1(:,1);
data1.c11=bearingDataMatrix1(:,2);
data1.c12=bearingDataMatrix1(:,3);
data1.c21=bearingDataMatrix1(:,4);
data1.c22=bearingDataMatrix1(:,5);
data1.d11=bearingDataMatrix1(:,6);
data1.d12=bearingDataMatrix1(:,7);
data1.d21=bearingDataMatrix1(:,8);
data1.d22=bearingDataMatrix1(:,9);

data2.So=bearingDataMatrix2(:,1);
data2.c11=bearingDataMatrix2(:,2);
data2.c12=bearingDataMatrix2(:,3);
data2.c21=bearingDataMatrix2(:,4);
data2.c22=bearingDataMatrix2(:,5);
data2.d11=bearingDataMatrix2(:,6);
data2.d12=bearingDataMatrix2(:,7);
data2.d21=bearingDataMatrix2(:,8);
data2.d22=bearingDataMatrix2(:,9);

Out = FFB_fittingToolCustomData(data1,data2,data);

% Converting struct to n x 9 matrix
bearingDataMatrix = [Out.So Out.c11 Out.c12 Out.c21 Out.c22...
    Out.d11 Out.d12 Out.d21 Out.d22];

end



function combinations = FFB_searchPerms(bearingData)

% Creating vectors for DIN table search tool. Search vector containing as
% follows [DIN,OmegaF_Omega,Z,Omega,fiiP1,h0max,dRB_CR,B_D]

%% Formatting variables
% If value = -1 means that variable in question is not included in table 
% search

DIN=-1;
OmegaF_Omega=-1;
Z=-1;
Omega=-1;
fiiP1=-1;
h0max=-1;
dRB_CR=-1;
B_D=-1;


%% Loading parameters
% Determining if needed to use user defined values or .Req.interp values 
% meaning variable in question needs to be interpolated. Values defined for
% use of table search are now converter to scalar variables.

% Defining standard type
DIN=bearingData.Req.DIN;

% Defining number of lobes/pads
Z=bearingData.Req.Z;

% Defining lobe/pad angle
Omega=bearingData.Req.Omega;

% Defining force direction
% try .. catch used for determining if variables exists, try catch because
% there were problems with exist() function possibly because of structural
% arrays

try
    
    % Trying to insert request value for interpolation, if it does not
    % exist it is causing error and going to catching point.
    fiiP1=bearingData.Req.interp.fiiP1;
    
catch
    
    % If request value for interpolation did not exist user defined value
    % is being used, meaning there is no need for interpolation.
    fiiP1=bearingData.Req.fiiP1;
    
end

% Defining geometric ratio
try
    B_D=bearingData.Req.interp.B_D;
catch
    B_D=bearingData.Req.B_D;
end

% Defining variables used only for multilobe bearings
if bearingData.Req.DIN == 2
    
    % Defining relative gap
    try
        h0max=bearingData.Req.interp.h0max;
    catch
        h0max=bearingData.Req.h0max;
    end
  
end

% Defining variables used only for tilting pad bearings
if bearingData.Req.DIN == 3
    
    % Defining pivot offset
    OmegaF_Omega=bearingData.Req.OmegaF_Omega;
    
    % Defining profiling
    try
        dRB_CR=bearingData.Req.interp.dRB_CR;
    catch
        dRB_CR=bearingData.Req.dRB_CR;
    end
    
end


%% Creating row vectors for bearing data table search
combinations=[];

for n1=1:length(DIN)
    
    for n2=1:length(OmegaF_Omega)
        
        for n3=1:length(Z)
            
            for n4=1:length(Omega)
                
                for n5=1:length(fiiP1)
                    
                    for n6=1:length(h0max)
                        
                        for n7=1:length(dRB_CR)
                            
                            for n8=1:length(B_D)
                                
                                combinations=[combinations
                                    [DIN(n1) OmegaF_Omega(n2) Z(n3) ...
                                    Omega(n4) fiiP1(n5) h0max(n6) ...
                                    dRB_CR(n7) B_D(n8)]];
                            end
                        end
                    end
                end
            end
        end
    end
end

end



function [tables,wc,wci] = FFB_createInterpVectors(combinations,tables,bearingData)

% Here is being created vectors of coefficients for interpolation. tables 
% vector includes all table numbers used in interpolation. wc vector 
% includes specified coefficients provided by DIN 31657 data tables. wci 
% vector includes user defined characteristic values.
% 
% Example: Tilting pad bearing (OmegaF_Omega=0.5, Z=4, Omega=60, B_D=0.625,
% fiiP1=45, dRB_CR=1.59)
% 
% tables=[21 27 23 29];       tables 21 and 23 represents geometric ratio (B_D) 0.5
%                             tables 27 and 29 represents geometric ratio (B_D) 0.75
%                             tables 21 and 27 represents profiling (dRB_CR) 2
%                             tables 23 and 29 represents profiling (dRB_CR) 3
% wc=[0.5 0.75 0.5 0.75 2 3]; coefficients wc(1:2) represents tables 21 and 27
%                             coefficients wc(3:4) represents tables 23 and 29
%                             coefficients wc(5:6) represents tables combinations of
%                             21 <-> 27 and 23 <-> 29
% wci=[0.625 0.625 1.59];     coefficients wci(1) and wci(2) are used for
% interpotion with respect to geometric ratio (B_D) and wci(3) for interpolation
% profiling (dRB_CD)
% 
%   [21 (wc=0.5)   27 (wc=0.75)]   [23 (wc=0.5)   29 (wc=0.75)]
%     \           /                  \            /
%   [wci=0.625 (wc=2)]             [wci=0.625 (wc=3)]
%                   \               /
%                 [wci=0.625, wci=1.59]
% 
% Length of vectors:
% If interpolation is needed to do with gap height ratio (h0max)/profiling (dRB_CR)
% and geometric ratio (B_D): tables 1 x 4, wc 1 x 6, wci 1 x 3
% If interpolation is needed to do with gap height ratio (h0max)/profiling (dRB_CR)
% or geometric ratio (B_D): tables 1 x 2, wc 1 x 2, wci 1 x 1
% If table search tool finds only one table using user defined values
% meaning table is directly usable: tables 1 x 1, wc empty, wci empty
% therefore no interpolation is needed

% Creating empty arrays for interpolation coefficients
wc=[];
wci=[];

% Interpolation using only gap height ratio (h0max)/profiling (dRB_CR) or 
% geometric ratio
if length(tables) == 2
    
    % Determines if using multilobe (2) or tilting pad (3) bearing tables
    if combinations(1,1) == 2
        
        % Interpolation with respect to geometric ratio (B_D)
        if combinations(1,6) == combinations(2,6)
            
            % Interpolation coefficients in order: geometric ratio (B_D) 
            % (2 coefficients)
            wc=[combinations(1,8) combinations(2,8)];
            
            % User defined interpolation coefficient (1 coefficient)
            wci=bearingData.Req.B_D;
        
        % Interpolation with respect to realtive gap (h0max)
        elseif combinations(1,8) == combinations(2,8)
            
            % Interpolation coefficients in order: gap height ratio (h0max) 
            % (2 coefficients)
            wc=[combinations(1,6) combinations(2,6)];
            
            % User defined interpolation coefficient (1 coefficient)
            wci=bearingData.Req.h0max;
            
        end
        
    % Determines if using multilobe (2) or tilting pad (3) bearing tables
    elseif combinations(1,1) == 3
        
        % Interpolation with respect to geometric ratio (B_D)
        if combinations(1,7) == combinations(2,7)
            
            % Interpolation coefficients in order: geometric ratio (B_D) 
            % (2 coefficients)
            wc=[combinations(1,8) combinations(2,8)];
            
            % User defined interpolation coefficient (1 coefficient)
            wci=bearingData.Req.B_D;
        
        % Interpolation with respect to profiling (dRB_CR)
        elseif combinations(1,8) == combinations(2,8)
            
            % Interpolation coefficients in order: profiling (dRB_CR) 
            % (2 coefficients)
            wc=[combinations(1,7) combinations(2,7)];
            
            % User defined interpolation coefficient (1 coefficient)
            wci=bearingData.Req.dRB_CR;
            
        end
    
    end % end of Determines if using multilobe (2) or tilting pad (3) 
    % bearing tables
    
end % end of Interpolation using only gap height ratio (h0max)/profiling (dRB_CR)
% or geometric ratio

% Interpolation using only gap height ratio (h0max)/profiling (dRB_CR) and 
% geometric ratio
if length(tables) == 4
    
    % Determines if using multilobe (2) or tilting pad (3) bearing tables
    if combinations(1,1) == 2
        
        % Interpolation coefficients in order: gap height ratio (h0max) 
        % (4 coefficients) and geometric ratio (B_D) (2 coefficients)
        wc=combinations(1:4,8)';
        wc=[wc combinations(1,6) combinations(3,6)];
        
        % User defined interpolation coefficients (3 coefficients)
        wci=[bearingData.Req.B_D...
            bearingData.Req.B_D bearingData.Req.h0max];
    
    % Determines if using multilobe (2) or tilting pad (3) bearing tables
    elseif combinations(1,1) == 3
        
        % Interpolation coefficients in order: profiling (dRB_CR) 
        % (4 coefficients) and geometric ratio (B_D) (2 coefficients)
        wc=combinations(1:4,8)';
        wc=[wc combinations(1,7) combinations(3,7)];
        
        % User defined interpolation coefficients (3 coefficients)
        wci=[bearingData.Req.B_D...
            bearingData.Req.B_D bearingData.Req.dRB_CR];
    
    end
    
end % end of Interpolation using only gap height ratio (h0max)/profiling (dRB_CR)
% and geometric ratio

end



function tables = FFB_searchTables(combinations)

% Searching suitable tables for interpolation using predefined row vectors

% Loading DIN 31657-2 and -3 tables
load FFB_BearingData.mat

tables=[];

% Going through all combined search rows
for j=1:size(combinations,1)
    
    % Loading one row at a time
    DIN=combinations(j,1);
    OmegaF_Omega=combinations(j,2);
    Z=combinations(j,3);
    Omega=combinations(j,4);
    fiiP1=combinations(j,5);
    h0max=combinations(j,6);
    dRB_CR=combinations(j,7);
    B_D=combinations(j,8);
    
    % Defining maximum number of tables
    if DIN==2, imax=43;
    elseif DIN==3, imax=56;
    end
    
    % Going through all tables of standard in question
    for i=1:imax
        
        % Loading every table section of standard for search algorithm using
        % eval() function
        if DIN==2, eval(['DINtable=bearings2.table' num2str(i) ';'])
        elseif DIN==3, eval(['DINtable=bearings3.table' num2str(i) ';'])
        end
        
        variables=0;
        hits=0;
        
        if Z > -1
            % If variable in question is not -1 meaning it is included in
            % search therefore number of search variables increases by one
            variables=variables+1;
            if Z==DINtable.Z
                % If table contains searched value hits variable increases
                % by one
                hits=hits+1;
            end,end
        
        if OmegaF_Omega > -1
            variables=variables+1;
            if OmegaF_Omega==DINtable.OmegaF_Omega
                hits=hits+1;
            end,end
        
        if Omega > -1
            variables=variables+1;
            if Omega==DINtable.Omega
                hits=hits+1;
            end,end
        
        if fiiP1 > -1
            variables=variables+1;
            if fiiP1==DINtable.fiiP1
                hits=hits+1;
            end,end
        
        if h0max > -1
            variables=variables+1;
            if h0max==DINtable.h0max
                hits=hits+1;
            end,end
        
        if B_D > -1
            variables=variables+1;
            if B_D==DINtable.B_D
                hits=hits+1;
            end,end
        
        if dRB_CR > -1
            variables=variables+1;
            if dRB_CR==DINtable.dRB_CR
                hits=hits+1;
            end,end
        
        % If number of all variables used in search is equal than founded
        % hits number of table is added into tables vector
        if hits==variables
            tables(end+1)=i;
        end
        
    end % end of Going through all tables of standard in question

end % end of Going through all combined search rows

end



function bearingDataMatrix = FFB_dataInterpolation(tables,wc,wci,DIN)

% This function controls two different type of interpolation functions for
% two (2) and four (4) table interpolation

% Loading DIN 31657 bearing data tables
load FFB_BearingData.mat
Input.bearings2=bearings2;
Input.bearings3=bearings3;

% Case where user defined values match directly in one table provided in
% standard, no interpolation is done
if length(tables) == 1
    
    eval(['Out.So=Input.bearings' num2str(DIN) '.table' num2str(tables) '.So;'])
    
    eval(['Out.c11=Input.bearings' num2str(DIN) '.table' num2str(tables) '.c11;'])
    eval(['Out.c12=Input.bearings' num2str(DIN) '.table' num2str(tables) '.c12;'])
    eval(['Out.c21=Input.bearings' num2str(DIN) '.table' num2str(tables) '.c21;'])
    eval(['Out.c22=Input.bearings' num2str(DIN) '.table' num2str(tables) '.c22;'])
    
    eval(['Out.d11=Input.bearings' num2str(DIN) '.table' num2str(tables) '.d11;'])
    eval(['Out.d12=Input.bearings' num2str(DIN) '.table' num2str(tables) '.d12;'])
    eval(['Out.d21=Input.bearings' num2str(DIN) '.table' num2str(tables) '.d21;'])
    eval(['Out.d22=Input.bearings' num2str(DIN) '.table' num2str(tables) '.d22;'])
    
% Case where interpolation is done using only geometric ratio or gap height
% ratio (h0max)/profiling (dRB_CR)
elseif length(tables) == 2
    
    % data struct for interpolation function
    data.DIN=DIN;
    data.tables=tables;
    data.wc=wc;
    data.wci=wci;
    Out=FFB_fittingTool(Input,data);
    
% Case where interpolation is done using gap height ratio (h0max)/profiling 
% (dRB_CR) and geometric ratio
elseif length(tables) == 4
    
    % data struct for interpolation function: first interpolation step is 
    % geometric ratio (B_D)
    data.DIN=DIN;
    data.tables=tables(1:2);
    data.wc=wc(1:2);
    data.wci=wci(1);
    Out1=FFB_fittingTool(Input,data);
    
    % data struct for interpolation function: second interpolation step is 
    % also geometric ratio (B_D)
    data.DIN=DIN;
    data.tables=tables(3:4);
    data.wc=wc(3:4);
    data.wci=wci(2);
    Out2=FFB_fittingTool(Input,data);
    
    % Third interpolation step: combining two previously created tables in
    % respect of geometric ratio (B_D). data struct for interpolation 
    % function: interpolation step is now gap height ratio (h0max) or profiling 
    % (dRB_CR) depending on bearing type
    data.wc=wc(5:6);
    data.wci=wci(3);
    Out=FFB_fittingToolCustomData(Out1,Out2,data);
    

end

% Converting Out struct into n x 9 matrix
bearingDataMatrix = [Out.So Out.c11 Out.c12 Out.c21 Out.c22...
    Out.d11 Out.d12 Out.d21 Out.d22];

end



function Output = FFB_fittingTool(Input,data)

% Interpolating two curves into one using specific coefficients and DIN
% 31657 bearing data tables


% Loading different values from specific bearing data table
eval(['x1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.So;'])
eval(['x2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.So;'])

% Output of Sommerfeld number of interpolated curve
Output.So=x1(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.c11;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.c11;'])

% Interpolating points of Somemrfeld number of second curve to respond 
% point of Sommerfeld number of first curve, without this step there can be
% differences because usually increment of Sommerfeld number between two
% curves is not constant.
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip'); % new interpolated c11 of second curve

% Empty array for results
Output.c11=[];

% Actual interpolation to desired coefficient is done by interpolating two
% curves with respect to their specific coefficients
for i=1:length(x1)
    
    % Interpolation using specific coefficients (data.wc) as x values and
    % stiffness/damping coefficient with respect to specific Sommerfeld
    % number in loop (x1) as y values
    Output.c11(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c11=Output.c11(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.c12;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.c12;'])

y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.c12=[];
for i=1:length(x1)
    Output.c12(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c12=Output.c12(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.c21;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.c21;'])

y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.c21=[];
for i=1:length(x1)
    Output.c21(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c21=Output.c21(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.c22;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.c22;'])

y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.c22=[];
for i=1:length(x1)
    Output.c22(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c22=Output.c22(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.d11;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.d11;'])

y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d11=[];
for i=1:length(x1)
    Output.d11(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d11=Output.d11(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.d12;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.d12;'])

y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d12=[];
for i=1:length(x1)
    Output.d12(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d12=Output.d12(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.d21;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.d21;'])

y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d21=[];
for i=1:length(x1)
    Output.d21(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d21=Output.d21(:);


eval(['y1=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(1)) '.d22;'])
eval(['y2=Input.bearings' num2str(data.DIN) '.table' num2str(data.tables(2)) '.d22;'])

y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d22=[];
for i=1:length(x1)
    Output.d22(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d22=Output.d22(:);

end



function Output = FFB_fittingToolCustomData(data1,data2,data)

% Interpolating two curves into one using specific coefficients and in this
% function also previously interpolated tables.

% Loading previously interpolated curves
x1=data1.So;
x2=data2.So;


% Output of Sommerfeld number of interpolated curve
Output.So=x1(:);


Output.c11=[];

y1=data1.c11;
y2=data2.c11;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip'); % new interpolated c11 of second curve

% Actual interpolation to desired coefficient is done by interpolating two
% curves with respect to their specific coefficients
for i=1:length(x1)
    
    % Interpolation using specific coefficients (data.wc) as x values and
    % stiffness/damping coefficient with respect to specific Sommerfeld
    % number in loop (x1) as y values
    Output.c11(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c11=Output.c11(:);


y1=data1.c12;
y2=data2.c12;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.c12=[];
for i=1:length(x1)
    Output.c12(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c12=Output.c12(:);


y1=data1.c21;
y2=data2.c21;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.c21=[];
for i=1:length(x1)
    Output.c21(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c21=Output.c21(:);


y1=data1.c22;
y2=data2.c22;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.c22=[];
for i=1:length(x1)
    Output.c22(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.c22=Output.c22(:);


y1=data1.d11;
y2=data2.d11;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d11=[];
for i=1:length(x1)
    Output.d11(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d11=Output.d11(:);


y1=data1.d12;
y2=data2.d12;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d12=[];
for i=1:length(x1)
    Output.d12(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d12=Output.d12(:);


y1=data1.d21;
y2=data2.d21;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d21=[];
for i=1:length(x1)
    Output.d21(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d21=Output.d21(:);


y1=data1.d22;
y2=data2.d22;
y2=interp1(x2,y2,x1(1:min([length(x1) length(x2)])),'pchip');

Output.d22=[];
for i=1:length(x1)
    Output.d22(end+1)=interp1([data.wc(1) data.wc(2)],[y1(i) y2(i)],data.wci,'pchip');
end
Output.d22=Output.d22(:);

end

