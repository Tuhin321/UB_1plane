function [xm, m,Jd,Jp,Mm,Matrx]=f_RotorMassProperty(Inp)
% This function calculates the rigid body mass properties of rotor in
% RoBeDyn 
% 
% Input: 
% Inp   structural array containing all RoBeDyn input data
%
% Output:
% xm = rotor's center of mass x-coordinate with respect to model origin
% m = Rotor's mass
% Jd = Diametral mass moment of inertia of the rotor with respect to center
%      of mass
% Jp = Polar mass moment of inertia of the rotor with respect to center
%      of mass
% Mm = Rigid body mass matrix (6x6)

%% DEVELOPMENT VERSION 11.11.2015
% NOT fully tested yet
% To be tested: 
% - Does the calculation fail if model contains support mass ??
% - Does it work with all models ?? 
%


%% Step 1
% Remove Constraints temporaly 
% Does not affect to Disp definition ???
Inp.Disp=[];

%% Step 2 
% Calculation of Rotor Matrices------------------------------------------
% These matrices are constant with respect to rotation speed
[Matrx]=f_Rotor_Matrix(Inp);
% Contents of structure Matrx:
% Matrx.Kc=Kc; % Stiffness Matrix
% Matrx.Cc=Cc; % Damping matrix
% Matrx.Gc=Gc; % Gyroscopic matrix (must be multiplied with Spin Speed)
% Matrx.Mc=Mc; % Mass matrix
% Matrx.FR=FR; % Vector of external applied forces
% Matrx.F2=F2; % Unbalance vector (see documentation)
% Matrx.F3=F3; % Unbalance vector (see documentation)
% Matrx.nnum2gdof=nnum2gdof; % Indexing matrix of unconstrained dofs
% Matrx.nnum2cgdof=nnum2cgdof; % Indexing matrix of constrained dofs

%% Development
% Here should be some routine which eliminates the support masses from
% the mass matrix Matrx.Mc
% Something that detects if a mass point in connected to beam element or
% not. In this I can see a problem if the masspoint is connected to shaft by
% springs. Should there be additional column in Node matrix, which tells if
% the node is rotating (1) or not-rotating (0)???
% OR just accept this flaw and be carefull when using it.

% THIS development created 30.6.2016 by JSo. Tested and verifed with one model
SupportMassGlobalDofs=[];
for j=1:size(Inp.Bearing,2) %all bearings in a model
    
    % It is assumed that IF there is bearing support defined, its node must
    % be defined as Jnode of bearing
    if Inp.Bearing(j).Jnode~=0
        
        %Find row of Bearing Jnode from node matrix
        NodeMatrRow=Inp.Node(:,1)==Inp.Bearing(j).Jnode;
        % remove node corresponding to NodeMatrxRow
        Inp.Node(NodeMatrRow,:)=[];
        
        % Add all global dofs of bearing Jnode into a list
        SupportMassGlobalDofs=[SupportMassGlobalDofs Matrx.nnum2gdof(Inp.Bearing(j).Jnode,:)];
        
    end
end
% Remove rows and colunms corresponding to support mass from the mass
% matrix
Matrx.Mc(:,SupportMassGlobalDofs)=[]; %Remove colums
Matrx.Mc(SupportMassGlobalDofs,:)=[]; %Remove rows

%% Creating rigid body modes manually (3 translations, 3 rotations)
%Create rigid body modes
for i=1:size(Inp.Node,1)
    if i==1; 
        v1x=[1 0 0 0 0 0]'; %translation x
        v1y=[0 1 0 0 0 0]'; %translation y
        v1z=[0 0 1 0 0 0]'; %translation z
        vrot_x=[0 0 0 1 0 0]'; %rotation x
        vrot_y=[0 0 -Inp.Node(i,2) 0 1 0]'; %rotation y %TARKISTA MERKIT!!!
        vrot_z=[0 Inp.Node(i,2) 0 0 0 1]'; %rotation z %TARKISTA MERKIT!!!
    else 
        v1x=[v1x     %translation x
            1; 0; 0; 0; 0; 0];
        v1y=[v1y      %translation y
            0; 1; 0; 0; 0; 0];
        v1z=[v1z      %translation z
            0; 0; 1; 0; 0; 0];
        vrot_x=[vrot_x     %rotation x
                0; 0; 0; 1; 0; 0];        
        vrot_y=[vrot_y     %rotation y
                0; 0; -Inp.Node(i,2); 0; 1; 0]; %TARKISTA MERKIT!!!
        vrot_z=[vrot_z     %rotation z
                0; Inp.Node(i,2); 0; 0; 0; 1]; %TARKISTA MERKIT!!!
    end;
end


% size(Matrx.Mc)
% size(vrot_z)
% size(v1y)

%% Evaluation of center of mass location 
% Center of mass (Lantto 1997)
% zm=abs(vrot_y'*Matrx.Mc*v1x/(v1x'*Matrx.Mc*v1x))
xm=abs(vrot_z'*Matrx.Mc*v1y/(v1y'*Matrx.Mc*v1y)); %Why abs ???
m=v1x'*Matrx.Mc*v1x; % Rotor's mass

% Create the rotational modes
v2x=vrot_x;
v2y=vrot_y+xm*v1z; %TARKISTA MERKIT!!!
v2z=vrot_z-xm*v1y; %TARKISTA MERKIT!!!


%Rigid body mode shape matrix
FiiR=[v1x v1y v1z v2x v2y v2z];

%% Calculation of rotor's mass properties
% Theoretical explanation lacks behind, but idea is to calculate the
% modal mass matrix using the rigid body modes. Rigid body modes are scaled
% so that each displacement has unity value (1). With these mode shapes
% the modal mass matrix should be diagonal and diagonal elements are as
% follows: Mm = diag(m,m,m,Jp,Jd,Jd). Mass moment's of inertia Jd and Jp 
% are with respect to center of mass.
% This seems to happen in case of simple Jeffcott rotor...

%Modal mass matrix
Mm=FiiR'*Matrx.Mc*FiiR;
m=Mm(1,1);
Jp=Mm(4,4);
Jd=Mm(5,5);

