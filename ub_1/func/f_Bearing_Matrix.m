function [Kb,Cb]=f_Bearing_Matrix(Bearing,omega)
% Calculates bearing materix
% 
% - Choises: Ball Bearing, Journal Bearing, HS Ball Bearing and Spherical Roller Bearing
% - Calculated the matrix from only one bearing 

% Update 2.5.2013 by J. Sopanen
% 'User Defined Speed Dependent A' option --> rotation speed is defined as
% rpm in the INDATA file. See user manual for details.
%
% Modifications:
% - 'HS Ball Bearing' type added by E. Kurvinen
% - 17.9.2014 'Spherical Roller Bearing' type coded by B. Ghalamchi, added 
%   by E. Sikanen 

% Initializing bearing matrices
Kb=zeros(3,3);
Cb=zeros(3,3);

% Testing if constant transformation matrix works on all files
Bearing.A=[0  0 -1
    0  1  0
    1  0  0];


% First, checking the bearing type
if strcmp(Bearing.type,'Journal Bearing') % strcmp(Inp.Bearing(1).type,'Journal Bearing')
    
    c_s=zeros(3,3);
    d_s=zeros(3,3);
    
    % Sommerfeld-Number 
    %    So = (F*psi^2)/(B*D*eta*omega) 
    % psi - relative clearance,
    % B - bearing width, 
    % eta - effective viscosity,
    % omega - angular velocity of the shaft
    % 
    % c_ik* = c_ik (C_R So)/F_r 
    % d_ik* = d_ik (omega C_R So)/F_r
    % 
    % C_R - radial clearance
    % F_r - radial load
    
    % calculating bearing forces
    F_r=sqrt(Bearing.LoadVec(2)^2+Bearing.LoadVec(3)^2);
    
    % viscosity of the lubricant (DIN 31654 T2)
    eta_x=0.18e-3; % supplementary variable
    eta = eta_x*exp((159.56/(Bearing.lubr_temp+95)-0.181913)*log(Bearing.lubr_rho*Bearing.VG/(1e6*eta_x)));
    
    % Sommerfeld-Number 
    if omega==0
        omega=omega+1*2*pi/60; % If rotation speed is zero, it is changed to 1. Getting rid of possibility of dividing by zero
    end
    So = (F_r*Bearing.psi^2)/((Bearing.B*1e-3)*(Bearing.D*1e-3)*eta*omega); % transformation mm --> m 
     
    
    % interpolation of dimensioless coefficients-------
    intp_method='pchip';
    
    % dimensioless stiffness coefficients
    c_s(1,1)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+1),So,intp_method);
    c_s(1,2)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+2),So,intp_method);
    c_s(2,1)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+3),So,intp_method);
    c_s(2,2)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+4),So,intp_method);
    % dimensioless damping coefficients
    d_s(1,1)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+5),So,intp_method);
    d_s(1,2)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+6),So,intp_method);
    d_s(2,1)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+7),So,intp_method);
    d_s(2,2)=interp1(Bearing.coeff(:,1),Bearing.coeff(:,1+8),So,intp_method);
    
    
    % Radial Clearance
    C_r=0.5*Bearing.D*Bearing.psi*1e-3; %radial clearance [mm] --> [m]
    
    % stiffness and damping matrices
    Kb= c_s*F_r/(So*C_r );  % C_r is inserted mm, changing to m
    Cb= d_s*F_r/(omega*So*C_r );
    
    % Transformation 
    Kb=Bearing.A'*Kb*Bearing.A;
    Cb=Bearing.A'*Cb*Bearing.A;
end % if Journal Bearing

% If ball bearing
if strcmp(Bearing.type,'Ball Bearing') % 
    
    
    % Initial values of bearing for iteration
    ex0=2*Bearing.cd; ey0=2*Bearing.cd; ez0=Bearing.cd;
    gammax0=0.0; gammay0=0.00; thetaz0=0;
    Disp0= [ex0, ey0, ez0, gammax0, gammay0, thetaz0];
    vx=0; vy=0; vz=0; omegax=0.0; omegay=0; omegaz=0*omega;
    Vel0= [vx, vy, vz, omegax, omegay, omegaz];
    mm=3; % number of iteration variables
     
    % calculating bearing equilibrium and stiffness matrix
    % [K_bear1,U_bear1,FB1]=f_Bearing_Equilibrium(FR1, Bearing(1), mm, Disp0, Vel0);
    [Kb]=f_Bearing_Equilibrium(Bearing.LoadVec, Bearing, mm, Disp0, Vel0);
    
    % Scaling the stiffness matrix
    Kb(1:3,1:3)=Kb(1:3,1:3)*1000; % Transformation N/mm --> N/m
    % moment (is not used)
    Kb(4:6,4:6)=Kb(4:6,4:6)*0.001; % Transformation Nmm/rad --> Nm/rad
    
    % Grabbing only diagonal terms
    Kb=diag(diag(Kb));
    
    % Damping matrix based on Krämerin
    Cb=2.5e-5*Kb;

end % if Ball Bearing

if strcmp(Bearing.type,'HS Ball Bearing') % 
    
    
    % Initial values of bearing for iteration
    ex0=2*Bearing.cd; ey0=2*Bearing.cd; ez0=Bearing.cd;
    gammax0=0.0; gammay0=0.00; thetaz0=0;
    Disp0= [ex0, ey0, ez0, gammax0, gammay0, thetaz0];
    vx=0; vy=0; vz=0; omegax=0.0; omegay=0; omegaz=0*omega;
    Vel0= [vx, vy, vz, omegax, omegay, omegaz];
    mm=3; % number of iteration parameters
    % calculating bearing equilibrium and stiffness matrix
    Bearing.omega=omega;
    % [K_bear1,U_bear1,FB1]=f_Bearing_Equilibrium(FR1, Bearing(1), mm, Disp0, Vel0);
    [Kb]=f_Bearing_Equilibrium(Bearing.LoadVec, Bearing, mm, Disp0, Vel0);
    
%     % Scaling the stiffness matrix
%     Kb(1:3,1:3)=Kb(1:3,1:3)*1000; % Transformation N/mm --> N/m
%     % moment (is not used)
%     Kb(4:6,4:6)=Kb(4:6,4:6)*0.001; % Transformation Nmm/rad --> Nm/rad
    
    % Grabbing only diagonal terms
    Kb=diag(diag(Kb));
    
    % Damping matrix based on Krämerin
    Cb=2.5e-5*Kb;

end % if Ball Bearing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SRB START%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Bearing.type,'Spherical Roller Bearing') % 
    
    
    % Initial values of bearing for iteration
    ex0=2*Bearing.cd; ey0=2*Bearing.cd; ez0=Bearing.cd;
    gammax0=0.0; gammay0=0.00; thetaz0=0;
    Disp0= [ex0, ey0, ez0, gammax0, gammay0, thetaz0];
    vx=0; vy=0; vz=0; omegax=0.0; omegay=0; omegaz=0*omega;
    Vel0= [vx, vy, vz, omegax, omegay, omegaz];
    mm=3; % number of iteration parameters
     
    % calculating bearing equilibrium and stiffness
    % [K_bear1,U_bear1,FB1]=f_Bearing_Equilibrium(FR1, Bearing(1), mm, Disp0, Vel0);
    [Kb]=f_Bearing_Equilibrium(Bearing.LoadVec, Bearing, mm, Disp0, Vel0);
    
    % Scaling bearing stiffness matrix
    Kb(1:3,1:3)=Kb(1:3,1:3)*1000; % Transformation N/mm --> N/m
    % moment (is not used)
    Kb(4:6,4:6)=Kb(4:6,4:6)*0.001; % Transformation Nmm/rad --> Nm/rad
    
    % Grabbing only diagonal terms
    Kb=diag(diag(Kb));
    
    % Damping matrix based on Krämerin
    Cb=2.5e-5*Kb;

end % if SRB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SRB END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If bearing type: Bearing Matrix
if strcmp(Bearing.type,'Bearing Matrix') % 
    
    Kb=Bearing.Kb;
    Cb=Bearing.Cb; 
end % if Bearing Matrix

% Jos laakerin tyyppi 'User Defined Speed Dependent A'
if strcmp(Bearing.type,'User Defined Speed Dependent A') % 
    
    % Format:        kb_n=[Speed,   kxx,    kyy,    kzz];
    % Units:               rpm,   N/m,    N/m,    N/m    
    
    % interpolation of coefficients-------
    intp_method='pchip';
    
    % convert units
    rpm=omega*60/(2*pi);
    
    % Stiffness terms are interpolated according to angular velocity
    kxx=interp1(Bearing.kb_n(:,1),Bearing.kb_n(:,1+1),rpm,intp_method);
    kyy=interp1(Bearing.kb_n(:,1),Bearing.kb_n(:,1+2),rpm,intp_method);
    kzz=interp1(Bearing.kb_n(:,1),Bearing.kb_n(:,1+3),rpm,intp_method);
    
    % Damping terms are interpolated according to angular velocity
    cxx=interp1(Bearing.cb_n(:,1),Bearing.cb_n(:,1+1),rpm,intp_method);
    cyy=interp1(Bearing.cb_n(:,1),Bearing.cb_n(:,1+2),rpm,intp_method);
    czz=interp1(Bearing.cb_n(:,1),Bearing.cb_n(:,1+3),rpm,intp_method);
    
    % forming the matrix
    Kb=[kxx   0    0
        0    kyy   0
        0     0   kzz];
    Cb=[cxx   0    0
        0    cyy   0
        0     0   czz];
end % if 'User Defined Speed Dependent A'


% If bearing type 'User Defined Speed Dependent B'
if strcmp(Bearing.type,'User Defined Speed Dependent B') % 
    
    % Stiffness matrix
    % Format:        kb_n=[Speed,   kxx,      kxy,       kyx,   kyy];
    % Units:               rpm,    N/mm,     N/mm,      N/mm,   N/mm 
    % Damping matrix
    % Format:        cb_n=[Speed,   cxx,    cxy,      cyx,     cyy];
    % Units:               rpm,    Ns/mm,   Ns/mm,    Ns/mm,  Ns/mm   
    
    % interpolation of coefficients-------
    intp_method='pchip';
    
    % convert units
    rpm=omega*60/(2*pi);
    
    % Stiffness terms are interpolated according to angular velocity
    kxx=interp1(Bearing.kb_n(:,1),Bearing.kb_n(:,1+1),rpm,intp_method);
    kxy=interp1(Bearing.kb_n(:,1),Bearing.kb_n(:,1+2),rpm,intp_method);
    kyx=interp1(Bearing.kb_n(:,1),Bearing.kb_n(:,1+3),rpm,intp_method);
    kyy=interp1(Bearing.kb_n(:,1),Bearing.kb_n(:,1+4),rpm,intp_method);
    
    % Damping terms are interpolated according to angular velocity
    cxx=interp1(Bearing.cb_n(:,1),Bearing.cb_n(:,1+1),rpm,intp_method);
    cxy=interp1(Bearing.cb_n(:,1),Bearing.cb_n(:,1+2),rpm,intp_method);
    cyx=interp1(Bearing.cb_n(:,1),Bearing.cb_n(:,1+3),rpm,intp_method);
    cyy=interp1(Bearing.cb_n(:,1),Bearing.cb_n(:,1+4),rpm,intp_method);
    
    %  forming matrices and unit transformation N/mm --> N/m
    Kb=1e3*[kxx  kxy   0
        kyx  kyy   0
        0     0    0];
    Cb=1e3*[cxx  cxy    0
        cyx  cyy   0
        0     0   0];
    
    % Transformation
    Kb=Bearing.A'*Kb*Bearing.A;
    Cb=Bearing.A'*Cb*Bearing.A;  
    
end % if 'User Defined Speed Dependent B'
