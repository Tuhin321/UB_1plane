% Simple Rotor with an overhung disk
% Rotor-Bearing Dynamic Analysis (RoBeDyn) Main Script
% April 2020
% ======================================================
% clearvars; close all;clc

function Features_matrix= robedynMain(Request,ub_magnphase)

disp('Start of RoBeDyn Analysis'); disp(' ');


% PREPROCESSING...........................................................
% Read Input Data
Inp=INDATA_fan(Request,ub_magnphase);

if strcmp(Request.show_plots, 'Yes')
    % Draw Model as a 3D WireFrame plot..........................................................
    f_RGeomPlotWireFr(Inp, [1 1 1 1 1], [0 1],{'-b',3});
    view(2); axis equal; axis vis3d; view3d rot %figure settings
    
    % % Free-free Modes..........................................................
    Request.FreeModes.NroFreqs=6;      % Number of frequencies to be solved and plotted
    Request.FreeModes.PlotModes='No'; % Plot mode shapes 'Yes' or 'No'
    % Mode plotting options
    %                           Mode#, ScaleFactor, linetype, linewidth
    Request.FreeModes.ModePlot={1, 1000, '-k',  2
                                2, 1000, '--k', 2
                                3, 12000, '-r',  2
                                4, 12000, '--r', 2  
                                5, 20000, '-g', 2
                                6, 20000, '--g', 2  
                                };
    Request.FreeModes.Options={'view','XY', 'originloc', [-0.3,0.0,0],'shading', 'on', 'transparency',0.4};   
    % Calculate Free-Free Modes
    % [FrqR,DampR,FiiR,Matrx]=f_FreeModes(Inp,Request)
%     [FrqR,DampR,FiiR,Matrx]=f_FreeModes(Inp,Request);
end

%% unbalance response
if strcmp(Request.UnbalanceResponse,'Yes')

    % Calculate and Plot Unbalance response of given nodes
    [Probe1_amp, Probe1_phase, probe2_amp, probe2_phase]=f_UBResponse(Inp,Request);   
end

Features_matrix=[Probe1_amp; Probe1_phase; probe2_amp; probe2_phase];

end


