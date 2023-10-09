function [bearingCoeff, AMBs]=bearingParameters(Inp, bearingConfig, Kx, Ki)
%% for bearing parameters ------------------------------------%
% nodes are taken directly, matrices are over ridden here
% Stiffness matrix (directions in global coordinate system, x=axial dir. yz=radial)
%     Kb=[kxx  kxy  kxz
%         kyx  kyy  kyz
%         kzx  kzy  kzz];

ki=2e6;
ci=7e4*[0   0   0
        0   1   0
        0   0   1];

switch bearingConfig
     case{'AB'}
        % bearing A (left side)
        Inp.Bearing(1).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(1).Cb= ci;

        % bearing B (right 1st)
        Inp.Bearing(2).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(2).Cb= ci;
    case{'ABD','ACD'}
        % bearing A (left side)
        Inp.Bearing(1).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(1).Cb= ci;

        % bearing B or C (right 1st)
        Inp.Bearing(2).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(2).Cb= ci;

        % bearing D (extended)
       Inp.Bearing(3).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(3).Cb= ci;
    case{'ABCD'}
        % bearing A (left side)
        Inp.Bearing(1).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(1).Cb= ci;

        % bearing B or C (right 1st)
        Inp.Bearing(2).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(2).Cb= ci;

        % bearing D (extended)
        Inp.Bearing(3).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(3).Cb= ci;
        
        % bearing 4 (right 3rd- EXTENDED)
        Inp.Bearing(4).Kb=[0   0     0
                           0   ki    0
                           0   0    ki];
        Inp.Bearing(4).Cb= ci;
    otherwise
        disp('invalid bearing configuration')
end
bearingCoeff=Inp.Bearing;

%% for AMB nodes and information

% Defining node locations for AMBs, sensors and impeller in the model

AMBs.Actuator(1).Inode = Inp.MassPoints(1,2); % impeller disc node location
AMBs.Actuator(1).Name = 'DE impeller';

AMBs.Actuator(2).Name = 'A End';
AMBs.Actuator(2).Inode = Inp.Bearing(1).Inode;    
% Inp.AMBs.Actuator(2).Kx = Kx;
% Inp.AMBs.Actuator(2).Ki = Ki;
AMBs.Actuator(2).Jnode=0; 

switch bearingConfig
    case{'AB'}
     AMBs.Actuator(3).Name = 'B End';
      AMBs.Actuator(3).Inode = Inp.Bearing(2).Inode; 
        % Inp.AMBs.Actuator(3).Kx = Kx;
        % Inp.AMBs.Actuator(3).Ki = Ki;
        AMBs.Actuator(3).Jnode=0; 
    case{'ABD','ACD'}
        switch bearingConfig
            case{'ABD'}
                AMBs.Actuator(3).Name = 'B End';
            case {'ACD'}
                AMBs.Actuator(3).Name = 'C End';
        end
        AMBs.Actuator(3).Inode = Inp.Bearing(2).Inode; 
        % Inp.AMBs.Actuator(3).Kx = Kx;
        % Inp.AMBs.Actuator(3).Ki = Ki;
        AMBs.Actuator(3).Jnode=0; 

        AMBs.Actuator(4).Name = 'D End';
        AMBs.Actuator(4).Inode = Inp.Bearing(3).Inode;   
        % Inp.AMBs.Actuator(4).Kx = Kx;
        % Inp.AMBs.Actuator(4).Ki = Ki;
        AMBs.Actuator(4).Jnode=0; 
        
    case{'ABCD'}
        
        AMBs.Actuator(3).Name = 'B End';
        AMBs.Actuator(3).Inode = Inp.Bearing(2).Inode; 
        % Inp.AMBs.Actuator(3).Kx = Kx;
        % Inp.AMBs.Actuator(3).Ki = Ki;
        AMBs.Actuator(3).Jnode=0; 

        
        AMBs.Actuator(4).Name = 'C End';
        AMBs.Actuator(4).Inode = Inp.Bearing(3).Inode; 
        % Inp.AMBs.Actuator(4).Kx = Kx;
        % Inp.AMBs.Actuator(4).Ki = Ki;
        AMBs.Actuator(4).Jnode=0; 
        
        AMBs.Actuator(5).Name = 'D End';
        AMBs.Actuator(5).Inode = Inp.Bearing(4).Inode;   
        % Inp.AMBs.Actuator(5).Kx = Kx;
        % Inp.AMBs.Actuator(5).Ki = Ki;
        AMBs.Actuator(5).Jnode=0; 
   otherwise
        disp('invalid bearing configuration')
end


% sensor nodes are common for all bearing configurations
AMBs.Sensors(1).Inode = Inp.sensorNodes(1);
AMBs.Sensors(1).name = 'A';

AMBs.Sensors(2).Inode = Inp.sensorNodes(2);
AMBs.Sensors(2).name = 'B';
% 
% AMBs.Sensors(3).Inode = Inp.sensorNodes(3);
% AMBs.Sensors(3).name = 'C';
% 
% AMBs.Sensors(4).Inode = Inp.sensorNodes(4); 
% AMBs.Sensors(4).name = 'D';