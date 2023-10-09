function [K_bear,Disp,FB]=f_Bearing_Equilibrium(FR, Bearing, mm, Disp0, Vel0)
% Ball Bearing Equilibrium analysis
% 
% INPUT: 
% FR = external force FR=[FX FY FZ] in global coordinate system
% Bearing =  structural array containing the bearing input data and
% parameters
% mm = number of parameters to be iterated (same as DOFs)
%        (if mm=3 only displacemtns are solved, if mm=5 --> also rotations are iterated)
% Disp0 = initial values for dipslacements [ex0, ey0, ez0, gammax0, gammay0, thetaz0]
% Vel0 = initial values for velocities [vx, vy, vz, omegax, omegay, omegaz]
%
% OUTPUT
% Disp = bearing displacements in a global coordinate system
% K_bear = bearing stiffness matrrix in a global coordinate system
% FB = bearing forces in a global coordinate system

% Modifications:
% - 'HS Ball Bearing' type added by E. Kurvinen
% - ee=ee0-inv(Jac)*F0; replaced to ee=ee0-Jac\F0; by E. Kurvinen

%%
%BALL BEARING HIGH-SPEED FORCES INCLUDED
if strcmp(Bearing.type,'HS Ball Bearing')

% Add moments into force vector
FR=[FR'; [0 0 0]']; %Force vector Axial (X),Radial (Y),Radial (Z)

% External Force
% Transformation matrix
Ab=[Bearing.A   zeros(3,3)
    zeros(3,3)  Bearing.A];
% In bearing coordinate system
% XY = radial direction
% Z = axial direction
FR=Ab'*FR; %[FX FY FZ MX MY MZ] FX=Radial,FY=Radial,FZ=Axial

%Initialize forces (mm = 3 or 5) variables to be integrated
F0=FR(1:mm)';

% Newton-Raphson %iteration loop
    
    % Calculation of the bearing equilibrium------------------------------- 
    
    Displ=f_HSBB_Disp(F0,Bearing,Bearing.omega);
%     pause
%     disp('eka')
%     pause
    h = (1+abs(F0))*1e-5;
    Jac = zeros(length(F0)); %Creates the empty Jacobian matrix
%     return
    for jcol=1:mm
        Fh = F0; % Set all elements to F0 values
        % Perturb jcol element of Fh
        Fh(jcol) = Fh(jcol)+h(jcol); % one variable at a time
        % Calculate variated displacement, Disph
        Disph=f_HSBB_Disp(Fh,Bearing,Bearing.omega); 
%         disp('toka')
        % Form column of Jacobian by finite differences in displacement
        Jac(:,jcol) = (Disph-Displ)/h(jcol);
    end
    
    %Jac %INVERSE! flexibility matrix to stiffness matrix
    Jacinv=inv(Jac);

%     % Solve new forces
%     Fnew=F0'-Jacinv\Displ; % temporaly variable ee is used 
%     % Save new values of forces
%     if mm==3
%         F0=[Fnew' FR(4:6)']';
%     elseif mm==5
%         F0= [Fnew' FR(6)']';
%     end 
%     
% %     % Set convergence criteria
%     if ii==1; Re=norm(F0)*1e-8; end; %ANSYS uses factor 0.005
%     % Check for convergence and end if converged
%     if norm(F0)<Re
%         %disp(['Solution CONVERGED at iteration number ' num2str(ii)])
%         %disp(['Convergence Norm = ' num2str(norm(F0),'%0.5e') '    Criterion = ' num2str(Re,'%0.5e')] )
%         %disp(Disp') % display coordinates of the end tip
%         return
%     end
% clear Jac ii jcol 

% Calculation of the bearing stiffness matrix in the previously solved equilibrium point

% % Bearing forces in the global coordinate system
% FnewGlobal=Ab*Fnew;
% 
% % Save Bearing force vector
% FB_loc=[F0' 0]'; 
% FB=Ab*FB_loc; % Force vector in the global coordinate system

% Save bearing stiffness matrix
% minus sign gives sensible results (NOT IN TURBO)
K_bear_loc=1*[Jacinv     zeros(3,3)
              zeros(3,3)    zeros(3,3)     ];
          
% bearing stiffness matrrix in a global coordinate system
K_bear=Ab'*K_bear_loc*Ab;
% pause
end

%%
%BALL BEARING MODEL WITHOUT HIGH-SPEED FORCES
if strcmp(Bearing.type,'Ball Bearing')

% Add moments into force vector
FR=[FR'; [0 0 0]'];

% External Force
% Transformation matrix
Ab=[Bearing.A   zeros(3,3)
    zeros(3,3)  Bearing.A];
% In bearing coordinate system
% XY = radial direction
% Z = axial direction
FR=Ab'*FR;

% Initialize displacements
Disp= Disp0; %[ex0, ey0, ez0, gammax0, gammay0, thetaz0];

%Initialize displacements
ee0=Disp(1:mm)';

%Initialize velocities
Vel=Vel0; % [vx, vy, vz, omegax, omegay, omegaz];

% Newton-Raphson %iteration loop
for ii=1:100 
    
    % Calculation of the bearing equilibrium------------------------------- 
    
    % bearing force vector at given operation point 
    FB=f_DGBB_Forces(Disp, Vel, Bearing);
    F0=FB(1:mm)-FR(1:mm); %Residual force
    
    % Set size of perturbation for each variable
    h = (1+abs(Disp))*1e-5; 
    
    % loop over colunms of the Jacobian
    for jcol=1:mm
        
        Disph = Disp; % Set all elements to x0 values
        % Perturb jcol element of Disph
        Disph(jcol) = Disp(jcol)+h(jcol);
        
        % Calculate variated Force Fh
        FBh=f_DGBB_Forces(Disph, Vel, Bearing);
        Fh=FBh(1:mm)-FR(1:mm);
        
        % Form column of Jacobian by finite differences
        Jac(:,jcol) = (Fh-F0)/h(jcol);
    end

    % Solve new displacements
    ee=ee0-Jac\F0; % temporaly variable ee is used 
    ee0=ee;
    % Save new values
    if mm==3 
        Disp= [ee' Disp0(4:6)]';
    elseif mm==5
        Disp= [ee'  Disp0(6)]';
    end
    
    % Set convergence criteria
    if ii==1; Re=norm(F0)*1e-8; end; %ANSYS uses factor 0.005
    % Check for convergence and end if converged
    if norm(F0)<Re
        %disp(['Solution CONVERGED at iteration number ' num2str(ii)])
        %disp(['Convergence Norm = ' num2str(norm(F0),'%0.5e') '    Criterion = ' num2str(Re,'%0.5e')] )
        %disp(Disp') % display coordinates of the end tip
        break
    end
end
clear Jac ii jcol 

% Calculation of the bearing stiffness matrix in the previously solved equilibrium point

% bearing force vector at given operation point 
F0=f_DGBB_Forces(Disp, Vel, Bearing);

% Set size of perturbation for each variable
h = (1+abs(Disp))*1e-5; 

% loop over columns of the Jacobian
for jcol=1:5
    
    Disph = Disp; % Set all elements to x0 values
    % Perturb jcol element of Disph
    Disph(jcol) = Disp(jcol)+h(jcol);
 
    % Calculate variated Force Fh
    Fh=f_DGBB_Forces(Disph, Vel, Bearing);
    
    % Form column of Jacobian by finite differences
    Jac(:,jcol) = (Fh-F0)/h(jcol);
end

% Bearing displacements in the global coordinate system
Disp=Ab*Disp;

% Save Bearing force vector
FB_loc=[F0' 0]'; 
FB=Ab*FB_loc; % Force vector in the global coordinate system

% Save bearing stiffness matrix
% minus sign gives sensible results
K_bear_loc=-1*[Jac     zeros(5,1)
              zeros(1,5)    0     ];
          
% bearing stiffness matrix in a global coordinate system
K_bear=Ab'*K_bear_loc*Ab; 
end
