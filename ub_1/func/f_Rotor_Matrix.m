function [Matrx]=f_Rotor_Matrix(Inp)
%%========================================================================
% [Matrx]=f_Rotor_Matrix(Inp)
% Function calculates rotor matrices. Rotor is modelled using 3D beam elements
% (Timoshenko) and using stiff disks. Masses of support and stiffness and damping
%  function are based on ANSYS 5.7 theory manual (BEAM4)
%
% OUTPUT: Matrx.
% -------
% Kg     = model global stiffness matrix
% Mg     = model global mass matrix
% Kc     = global stiffness matrix with constraints take into account
% Mc     = global stiffness matrix with constraints take into account
% FR     = vector for external forces
% nnum2gdof= node numbers and corresponding global DOFs
%            nnum2gdof(node,:)=[UX, UY, UZ, ROTX, ROTY, ROTZ]    
%                - global dof is e.g. 
%                  node 10, dof UX -->  gdof=nnum2gdof(10,1)            
% INPUT: Inp.
% ------
% Node    = node matrix         [ID, X, Y, Z]
% Elem    = element matrix      [ID, I, J, Mat, Real]
% MassPoints= point mass matrix [ID, node, mass, Jxx, Jyy, Jzz, L, dia]
% Real    = Real Constant       [ID, A, Ixx, Iyy, Izz, B, H, sheary, shearz]
% Mat     = material properties [ID, E, G, nuxy, rho]
% Disp    = constraints         [Node, Dof, Value] (Dof: 1=UX, 2=UY, 3=UZ, 4=ROTX, 5=ROTY, 6=ROTZ)
% Force   = node forces         [Node, Dof, Value]
% lumpm   = mass matrix type, 0=consistency, 1=lumped
% grav    = Gravity acceleration vector [gy, gz] eg. Inp.grav=[-9.81,0];
% 
% Written by Jussi Sopanen, LTKK/KoA, 2002
% ========================================================================

%-Modifications------------------------
% 24.2.2006 - Disk gyromatrix corrected, terms were in wrong order
%           - Renamed F2 --> Fs and F3 --> Fc
% 16.7.2013 - Added Jxx into beam element mass matrix calculation
% 29.8.2014 - Added gravity force vector calculation
%----------------------------------

%% Tallennetaan Inputit
Node=Inp.Node;
Elem=Inp.Elem;
MassPoints=Inp.MassPoints;
Real=Inp.Real;
Mat=Inp.Mat;
Force=Inp.Force;
lumpm=Inp.lumpm;
Disp=Inp.Disp;
UB=Inp.UB;
SpringDamper=Inp.SpringDamper;


% Calculating model DOFs-------------------------
doftmp=zeros( max(max(Elem(:,2:3))),1 ); % vector preallocation, length = number of DOFs

% Calculating number of nodes used in Elem matrix
% Going through all node numbers from min to max node
MinNode=max([1 min(min(Elem(:,2:3)))]);
MaxNode=max(max(Elem(:,2:3)));
for aa=MinNode:MaxNode
    for bb=2:3
        for cc=1:size(Elem,1)
            % If element node number is I or J, value will be 1 in the table
            if aa==Elem(cc,bb)
                doftmp(aa)=1;
            end
        end    
    end
end
% Number of active DOFs
DOFs=sum(doftmp)*6;
% str=['Number of DOFs=' int2str(DOFs)];
% disp(str)
clear doftmp aa bb cc MinNode MaxNode     
%----------------------------------------------------------


% Node matrix is sorted for indexing global DOFs
[Y,indx]=sort(Node(:,1)); % Y is node number, indx is corresponding index Node matrix

%% Going through elements and forming matrices
for ii=1:size(Elem,1)
    
    
    
    % Calculating global DOFs and connections of I and J nodes
    for jj=1:size(Node,1)
        % Getting row of I node from Node Matrix
        if Elem(ii,2)==Node(jj,1) %
            Iindx=jj; 
            
            % Start index if I node global DOFs
            % Y include node numbers sorted from highest t lowest
            % - if I node number is smallest, global DOF Igdof is 1
            % - if second smallest, global DOF Igdof is 7
            % - etc.
            for nn=1:size(Node,1)
                if Node(Iindx,1)==Y(nn)
                    Igdof=(nn-1)*6+1;
                    % node number and corresponding global dofs
                    %     nnum2gdof(node,:)=[UX, UY, UZ, ROTX, ROTY, ROTZ]    
                    %     -global dof is get by: 
                    %        node 10, dof UX -->  gdof=nnum2gdof(10,1)
                    nnum2gdof(Y(nn),:)=[(nn-1)*6+1 (nn-1)*6+2 (nn-1)*6+3 (nn-1)*6+4 (nn-1)*6+5 (nn-1)*6+6 ];
                end
            end
        end
        % Getting row of J node from Node matrix
        if Elem(ii,3)==Node(jj,1)
            Jindx=jj; 
            
            % start index of J node global dofs
            % Y include node numbers sorted from highest t lowest
            % - if I node number is smallest, global DOF Igdof is 1
            % - if second smallest, global DOF Igdof is 7
            % - etc.
            for nn=1:size(Node,1)
                if Node(Jindx,1)==Y(nn)
                    Jgdof=(nn-1)*6+1;
                    % node number and corresponding global dofs
                    %     nnum2gdof(node,:)=[UX, UY, UZ, ROTX, ROTY, ROTZ]         
                    nnum2gdof(Y(nn),:)=[(nn-1)*6+1 (nn-1)*6+2 (nn-1)*6+3 (nn-1)*6+4 (nn-1)*6+5 (nn-1)*6+6  ];
                end
            end    
        end
    end
    
    
    if ~Elem(ii,3) == 0  % If not mass element
        
        % Element length
        L=sqrt((Node(Jindx,2)-Node(Iindx,2) )^2+ (Node(Jindx,3)-Node(Iindx,3))^2+...
            (Node(Jindx,4)-Node(Iindx,4))^2);
        
        % Element node coordinates
        coord=[Node(Iindx,2) Node(Jindx,2) Node(Iindx,3) Node(Jindx,3) Node(Iindx,4) Node(Jindx,4)];
        
        % Element properties
        ind=find(Real(:,1)==Elem(ii,5)); % searching row corresponding to real constant
        A=Real(ind,2);    % cross-section area
        Izz=Real(ind,3);  % moment of inertia YY
        Iyy=Real(ind,4);  % moment of inertia ZZ
        theta=Real(ind,7);
        if Real(ind,8)==0; % If torsional stiffness is not definedy, Ixx is used to calculate ki
            Iv=Real(ind,9); 
        else 
            Iv=Real(ind,8); 
        end 
        Ixx=Real(ind,9);
        shearz=Real(ind,10); % shear
        sheary=Real(ind,11); % shear
        
        
        % Material properties
        ind=find(Mat(:,1)==Elem(ii,4)); % searching row corresponding to real constant
        E=Mat(ind,2);
        nuxy=Mat(ind,3);
        rho=Mat(ind,4);
        
        %3D beam element stiffness matrix
        %function [K3]=K3D(A,E,Ixx,Iyy,Izz,L,nuxy,sheary,shearz,Igdof,Jgdof,DOFs,theta,coord);
        [ki]=K3D(A,E,Iv,Iyy,Izz,L,nuxy,sheary,shearz,Igdof,Jgdof,DOFs,theta,coord);
        
        % Summing elements' stiffness matrices
        if ii==1
            Kg=ki;
        else
            Kg=Kg+ki;
        end
        
        %3D beam element mass matrix
        % function [M3]=Mass3D(A,E,rho,Ixx,Iyy,Izz,L,nuxy,sheary,shearz,lumpm,Igdof,Jgdof,DOFs,theta,coord);
        [mi]=Mass3D(A,E,rho,Ixx,Iyy,Izz,L,nuxy,sheary,shearz,lumpm,Igdof,Jgdof,DOFs,theta,coord);
        
        % Summing elements' mass matrices
        if ii==1
            Mg=mi;
        else
            Mg=Mg+mi;
        end
        
        %3D beam element gyroscopic matrix
        % function [G3]=Gyro3D(A,E,rho,I,L,nuxy,shear,Igdof,Jgdof,DOFs,coord)
        [gi]=Gyro3D(A,E,rho,Iyy,L,nuxy,sheary,Igdof,Jgdof,DOFs,coord);
        
        % Summing elements' gyro matrices
        if ii==1
            Gg=gi;
        else
            Gg=Gg+gi;
        end
        
        %% Gravity force vector calculation
        %Inp.grav=[gy,gz]
        if isfield(Inp,'grav') %check if gravity vector is defined in indata Inp.grav=[-9.81,0];
            %Global gravity force vector (complete model, all dofs)
            if ii==1 Fg=zeros(DOFs,1); end %DOFS
            % Calculate element gravity force vector
            fgi = rho*A*[0, 1/2*Inp.grav(1)*L, 1/2*Inp.grav(2)*L, 0, -1/12*Inp.grav(2)*L^2, 1/12*Inp.grav(1)*L^2, ...
                         0, 1/2*Inp.grav(1)*L, 1/2*Inp.grav(2)*L, 0,  1/12*Inp.grav(2)*L^2, -1/12*Inp.grav(1)*L^2]';
            % Add to global matrix
            Fg(Igdof:Igdof+5)=Fg(Igdof:Igdof+5)+fgi(1:6); 
            Fg(Jgdof:Jgdof+5)=Fg(Jgdof:Jgdof+5)+fgi(7:12); 
        end
        
    end % if
end

%% **************************************************************************
% Adding point mass values to mass mtrix and gyro matrix
%MassPoints=[ID, node, mass, Jxx,Jyy,Jzz,L,dia]
for ii=1:size(MassPoints,1)
    
    Node_m1=MassPoints(ii,2); % point mass node number
    h=MassPoints(ii,7); % offset
    % mass matrix terms of node
    mp=zeros(6,6);
    mp(1,1)=MassPoints(ii,3); 
    mp(2,2)=MassPoints(ii,3); 
    mp(3,3)=MassPoints(ii,3); 
    mp(4,4)=MassPoints(ii,4); % Jxx = Jp
    mp(5,5)=MassPoints(ii,5)+MassPoints(ii,3)*h^2; % Jyy=Jd
    mp(6,6)=MassPoints(ii,6)+MassPoints(ii,3)*h^2; % Jzz=Jd
    % offset terms
    mp(2,6)=MassPoints(ii,3)*h;
    mp(6,2)=MassPoints(ii,3)*h;
    mp(3,5)=-MassPoints(ii,3)*h;
    mp(5,3)=-MassPoints(ii,3)*h;
    
    % gyro matrix terms of node
    % NOTE! Rotor is assumed to be along x axis
    gp=zeros(6,6);
    gp(5,6)=+MassPoints(ii,4); % Omega*Jxx later on multiplied by Spin*
    gp(6,5)=-MassPoints(ii,4); % -Omega*Jxx later on multiplied by Spin*
    
    %% Add mass points to Gravity force vector
    %Inp.grav=[gy,gz]
    if isfield(Inp,'grav') % check if gravity vector is defined in indata Inp.grav=[-9.81,0];
        
        % Calculate element gravity force vector
        fgp = MassPoints(ii,3)*[0, Inp.grav(1), Inp.grav(2), 0, 0, 0];
        
        % Get global dof indexes
        gdof_y=nnum2gdof(Node_m1,2); %get index node translation y
        gdof_z=nnum2gdof(Node_m1,3); %get index node translation 3
        % Add to global gravity force vector
        Fg(gdof_y)=Fg(gdof_y)+fgp(2);
        Fg(gdof_z)=Fg(gdof_z)+fgp(3);
    end

    %% Adding terms to matrices
    for i=1:6 % rows 1-6
        for j=1:6 % columns 1-6
            
            % point mass global dofs
            gdof_r=nnum2gdof(Node_m1,i); % Getting global dof for row
            gdof_c=nnum2gdof(Node_m1,j); % Getting global dof for column
            % Mass matrix
            Mg(gdof_r,gdof_c)=Mg(gdof_r,gdof_c) + mp(i,j); % adding point mass values
            % Gyro matrix
            Gg(gdof_r,gdof_c)=Gg(gdof_r,gdof_c) + gp(i,j); % adding point mass values
            % Gravity force vector
            
        end
    end
end
clear i j gdof_r gdof_c mp gp Node_m1


%% Adding Spring-Damper parameters to stiffness and damping matrices
Cg=zeros(DOFs,DOFs); % Preallocating damping matrix
%   SpringDamper=[ID, Inode, Jnode, Type, Dir,  Value];
% Type = 1 --> spring
% Type = 2 --> damper
% Dir 1=X, 2=Y, 3=Z
% If J node is 0, pring is attached to ground
for ii=1:size(SpringDamper,1)
    
    NodeI=SpringDamper(ii,2); % spring I node
    NodeJ=SpringDamper(ii,3); % spring J node
    
    % spring global dofs
    Igdof=nnum2gdof(NodeI,SpringDamper(ii,5)); % Getting global dof 
    
    % Inserting stiffness and damping coefficient
    if SpringDamper(ii,4)==1 % if spring
        Kg(Igdof,Igdof)=Kg(Igdof,Igdof) + SpringDamper(ii,6);
    else
        Cg(Igdof,Igdof)=Cg(Igdof,Igdof) + SpringDamper(ii,6);
    end
    
    if ~NodeJ==0 % if other end is not attached to ground
        Jgdof=nnum2gdof(NodeJ,SpringDamper(ii,5)); % getting global dof 
        
        % Inserting stiffness or damping coefficient
        if SpringDamper(ii,4)==1 % if spring
            Kg(Jgdof,Jgdof)=Kg(Jgdof,Jgdof) + SpringDamper(ii,6);
            % cross-coupled terms
            Kg(Igdof,Jgdof)=Kg(Igdof,Jgdof) - SpringDamper(ii,6);
            Kg(Jgdof,Igdof)=Kg(Jgdof,Igdof) - SpringDamper(ii,6);
        else
            Cg(Jgdof,Jgdof)=Cg(Jgdof,Jgdof) + SpringDamper(ii,6);
            % cross-coupled terms
            Cg(Igdof,Jgdof)=Cg(Igdof,Jgdof) - SpringDamper(ii,6);
            Cg(Jgdof,Igdof)=Cg(Jgdof,Igdof) - SpringDamper(ii,6);
        end
        
    end
end
clear ii gdof NodeI NodeJ

%% Force vector-------------------------------------------------
% Forces
%-----Node, Dof, Value
% e.g. Force=[10,  1,   -1000];
FR=zeros(DOFs,1);
for ii=1:size(Force,1)
    gdof=nnum2gdof(Force(ii,1),Force(ii,2)); % getting global DOF
    FR(gdof)=Force(ii,3); % adding value to force vector
end
clear gdof

% Unbalance masses-------------------------------------------
% -- Node, value (kg*m), angle
% UB=[1, 0.001*0.5, 0];
F2=zeros(DOFs,1); F3=zeros(DOFs,1);
for ii=1:size(UB,1)
    gdof_Y=nnum2gdof(UB(ii,1),2); % getting global DOF Y
    gdof_Z=nnum2gdof(UB(ii,1),3); % getting global DOF Y
    % Unbalance forces
    % UB(ii,3) = phase angle alfa of unbalance, starting according to positive Y axis
    % and rotates clockwise (see theory manual)
    F2(gdof_Y)= -UB(ii,2)*sin(UB(ii,3)) ; % adding value to force vector
    F2(gdof_Z)= UB(ii,2)*cos(UB(ii,3)) ; % adding value to force vector
    F3(gdof_Y)= UB(ii,2)*cos(UB(ii,3))  ; % adding value to force vector
    F3(gdof_Z)= UB(ii,2)*sin(UB(ii,3)) ; % adding value to force vector
end




%% Taking into account the constraints -----------------------------------------------------

% Constraints
%----------Node,  Dof, Value
% e.g.  Disp=[5,    1,   0
%             5,    2,   0
%             5,    3,   0];

%--------------------------------------------------------------------------
% index matrix where constraints are taken into account
nnum2cgdof=nnum2gdof;
% Going through all costraint matrix loops
for ii=1:size(Disp,1)
    % node number and dof
    node=Disp(ii,1);
    dof=Disp(ii,2);
    % going through all index matrix cells
    for jj=1:size(nnum2gdof,1)
        for kk=1:6
            % if constrained dof is smaller, subtracting 1
            if nnum2gdof(node,dof) < nnum2gdof(jj,kk)
                nnum2cgdof(jj,kk) = max(0,nnum2cgdof(jj,kk)-1); % max, if not appear to be negative
            %  if constrained dof = gdof, adding 0 to matrix   
            elseif nnum2gdof(node,dof) == nnum2gdof(jj,kk)
                nnum2cgdof(jj,kk) = 0;
            end
        end
    end
end
clear ii jj node dof

% Taking constraints into account-------------------------------------------------
% saving constrained global dofs' numbers to vector
% Preallocating constrained matrices
Kc=Kg;  Mc=Mg;   Gc=Gg; Cc=Cg;

if ~isempty(Disp)
    for ii=1:size(Disp,1)
        gdof(ii)=nnum2gdof(Disp(ii,1),Disp(ii,2)); % Getting global dof
    end
    gdof=sort(gdof); % sorting from lowest to highest
    

    % Removing rows and columns corresponding to constrained dofs
    for ii=1:size(gdof,2)
        % Stiffness matrix
        Kc(:,gdof(ii))=[]; % removing column
        Kc(gdof(ii),:)=[]; % removing row
        % Mass matrix
        Mc(:,gdof(ii))=[]; % removing column
        Mc(gdof(ii),:)=[]; % removing row
        % Damping matrix
        Cc(:,gdof(ii))=[]; % removing column
        Cc(gdof(ii),:)=[]; % removing row
        % Gyroscopic matrix
        Gc(:,gdof(ii))=[]; % removing column
        Gc(gdof(ii),:)=[]; % removing row
        % Force vector
        FR(gdof(ii))=[]; % removing row        
        % Unbalance force vector
        F2(gdof(ii))=[]; % removing row 
        F3(gdof(ii))=[]; % removing row 
        % Gravity vector
        if isfield(Inp,'grav')
            Fg(gdof(ii))=[]; % delete row
        end
        % reducing one from global dofs
        gdof=gdof-1;
    end
end

%% UNDER CONSTRUCTION..................................................
% Apply Gear Constraints in Torsional Vibration analysis-------------------
% In torsional vibration modelling shafts are made equivalent with respect
% to 1st shaft by multiplying stiffnesses and mass moments of inertia with
% gear ratio n
% Example definition in Indata file
%-Gear constraints (Torsional analysis) --------------------%                
%   Gear=[ID, Inode, Jnode, GearRatio];
% Inp.Gear=[1,      3,     4,  2
%           2,      3,     7,  2
%           ];

%% Add ModalDamping if specified
if isfield(Inp,'ModalDamping')
    % Creating damping matrix using modal damping coefficients
    % Solution of eigenvalue problem
    [FiiN d]=eig(Kc,Mc);
    % Sorting matrices according to frequency
    [omega, idx]=sort(sqrt(diag(d)));
    FiiN=FiiN(:,idx); clear idx d;

    % Trying to remove rigid body modes
    % omega2=omega(find(omega>0.1));

    nrbm=min(find(omega>0.1))-1; % number of rigid body modes

    Ctmp=zeros(size(Cc));
    for i=1:size(Inp.ModalDamping,1)
        % setting damping ratios to diagonal
        mode=Inp.ModalDamping(i,1);
        Ctmp(mode+nrbm,mode+nrbm)=2*Inp.ModalDamping(i,2)*omega(mode+nrbm);
    end
    % Damping matrix using modal coefficients
    Cm=(FiiN')^(-1)*Ctmp*FiiN^(-1);
else
    Cm=zeros(size(Kc));
end




%% Creating structure for return parameters
% [Kg,Mg,Gg,Kc,Mc,Gc,FR,F2,F3,Fg,nnum2gdof,nnum2cgdof]

Matrx.Kc=Kc; % constrained stiffness matrix
Matrx.Cc=Cc+Cm; % constrained damping matrix
Matrx.Gc=Gc; % constrained gyro matrix (later on multiplied by rotational speed)
Matrx.Mc=Mc; % constrained mass matrix
Matrx.FR=FR; % External force vector
Matrx.F2=F2; % unbalance vector
Matrx.F3=F3; % unbalance vector
if isfield(Inp,'grav') Matrx.Fg=Fg; end % gravity force vector
Matrx.nnum2gdof=nnum2gdof; % index vector for unconstrained dofs
Matrx.nnum2cgdof=nnum2cgdof; % index vector for constrained dofs

%% LOCAL FUNCTIONS******************************************************************
%**************************************************************************

function [K3]=K3D(A,E,Ixx,Iyy,Izz,L,nuxy,sheary,shearz,Igdof,Jgdof,DOFs,theta,coord);
%==============================================================
% 3D beam element stiffness matrix:
%[K3]=K3D(A,E,Ixx,Iyy,Izz,L,nuxy,sheary,shearz,Igdof,Jgdof,DOFs,theta,coord);
% -----------------------------------------------------
% Function calculates global stiffness matrix for3D beam element (Timoshenko) 
%
%Parameters:
%=========== 
%
% - A     = Element cross-sectional area
% - E     = Material elastic modulus
% - Ixx   = Beam moment of inertia
% - Iyy   = Beam moment of inertia
% - Izz   = Beam moment of inertia
% - L     = Beam length
% - Lxy   = Beam projected length to xy plane
% - nuxy  = Poisson's vakio
% - sheary= shear strain coefficient
% - shearz= shear strain coefficient
% - Igdof = beam I node global dof number (=UX global)
% - Jgdof = beam J node global dof number (=UX global)
% - DOFs  = number of dofs
% - theta = beam orientation around x axis
% - coord = Element node coordinates
%
% Written by Jussi Sopanen, LTKK/KoA, 2002
% ===============================================================

G=E/(2*(1+nuxy)); % shear modulus
% If shearz =0 shear strain is not taken into account
if sheary==0
    fiiz=0;
else
    fiiz=12*E*Iyy/(G*(A*sheary)*L^2);
end

% If sheary =0 shear strain is not taken into account
if shearz==0
    fiiy=0;
else
    fiiy=12*E*Izz/(G*(A*shearz)*L^2);
end

% Stiffness matrix cells (ANSYS Theory)
J=Ixx;
AA=E*A/L;
BB=G*J/L;
Ay=12*E*Iyy/(L^3*(1+fiiz));
Az=12*E*Izz/(L^3*(1+fiiy));
By=-12*E*Iyy/(L^3*(1+fiiz));
Bz=-12*E*Izz/(L^3*(1+fiiy));
Cy=6*E*Iyy/(L^2*(1+fiiz));
Cz=6*E*Izz/(L^2*(1+fiiy));
Dy=-6*E*Iyy/(L^2*(1+fiiz));
Dz=-6*E*Izz/(L^2*(1+fiiy));
Ey=E*Iyy*(4+fiiz)/(L*(1+fiiz));
Ez=E*Izz*(4+fiiy)/(L*(1+fiiy));
Fy=E*Iyy*(2-fiiz)/(L*(1+fiiz));
Fz=E*Izz*(2-fiiy)/(L*(1+fiiy));


K=[  AA   0    0    0     0    0    -AA   0    0    0     0    0
     0    Az   0    0     0    Cz    0    Bz   0    0     0    Cz
     0    0    Ay   0     Dy   0     0    0    By   0     Dy   0
     0    0    0    BB    0    0     0    0    0   -BB    0    0
     0    0    Dy   0     Ey   0     0    0    Cy   0     Fy   0
     0    Cz   0    0     0    Ez    0    Dz   0    0     0    Fz
    -AA   0    0    0     0    0     AA   0    0    0     0    0
     0    Bz   0    0     0    Dz    0    Az   0    0     0    Dz
     0    0    By   0     Cy   0     0    0    Ay   0     Cy   0
     0    0    0   -BB    0    0     0    0    0    BB    0    0
     0    0    Dy   0     Fy   0     0    0    Cy   0     Ey   0
     0    Cz   0    0     0    Fz    0    Dz   0    0     0    Ez]; 

%Fixed 3.9.2014 the following terms by JSo
% K(9,11)=Cy; K(11,9)=Cy;
% K(8,12)=Dz; K(12,8)=Dz; 

% Transformation matrix
T=T3D(theta,L,coord);

% Global stiffness matrix (of element)
KG=T'*K*T;

% Global stiffness matrix (whole model)
K3=zeros(DOFs); %DOFS

% 6 = number of node dofs
for ii=1:6 % rows
    for jj=1:6 % columns
        
        % upper-left
        K3(Igdof+(ii-1),Igdof+(jj-1))=KG(ii,jj);
        %upper-right
        K3(Igdof+(ii-1),Jgdof+(jj-1))=KG(ii,jj+6);
        % lower-left
        K3(Jgdof+(ii-1),Igdof+(jj-1))=KG(ii+6,jj);
        % lower-right
        K3(Jgdof+(ii-1),Jgdof+(jj-1))=KG(ii+6,jj+6);
        
    end
end

function [M3]=Mass3D(A,E,rho,Ixx,Iyy,Izz,L,nuxy,sheary,shearz,lumpm,Igdof,Jgdof,DOFs,theta,coord);
%========================================================================
% 3D beam element mass matrix:
% [M3]=Mass3D(A,E,rho,Ixx,Iyy,Izz,L,nuxy,sheary,shearz,lumpm,Igdof,Jgdof,DOFs,theta,coord)
%------------------------------------------------------------------
% Function calculates global mass amtrix for 3D beam element (Timoshenko)
%
% Paramaters:
% -----------
%
% - A     = Element cross-sectional area
% - E     = Material elastic modulus
% - rho   = Material density
% - Ixx   = Polar area moment of inertia
% - Iyy   = area moment of inertia about y-axis
% - Izz   = area moment of inertia about z-axis
% - L     = Beam length
% - nuxy  = Poisson's ratio
% - sheary= shear strain coefficient
% - shearz= shear strain coefficient
% - lumpm = mass matrix type (1=lumped, 0=consistency)
% - Igdof = beam I node global dof number (=UX global)
% - Jgdof = beam J node global dof number (=UX global)
% - DOFs  = number of dofs
% - theta = beam orientation around x axis
% - coord = Element node coordinates
% 
% Written by Jussi Sopanen, LTKK/KoA, 2002
% ========================================================================


% cross-section radius of gyration
ry=sqrt(Iyy/A); % radius of gyration
rz=sqrt(Izz/A); % radius of gyration
% polar moment of inertia
% ----Update 16.7.2013
Jx=Ixx; 
if (Iyy==0 && Izz==0)
    disp('Cross-section properties Iyy and Izz are zero, only Torsion Analysis is possible')
end
if Ixx==0
    Jx=Iyy+Izz; % Use default value
    disp('Value not defined for Ixx, a default value Ixx=Iyy+Izz is used')
end
%------------------
% If shearz =0 shear strain is not taken into account
if sheary==0
    fiiz=0;
else
    G=E/(2*(1+nuxy)); % glide modulus% glide modulus
    fiiz=12*E*Iyy/(G*(A*sheary)*L^2);
end
% fiiz=0
% If sheary =0 shear strain is not taken into account
if shearz==0
    fiiy=0;
else
    G=E/(2*(1+nuxy)); % glide modulus
    fiiy=12*E*Izz/(G*(A*shearz)*L^2);
end
% fiiy=0
% forming mass matrix
if lumpm==0
    % Consistency mass mass matrix cells (ANSYS Theory)
    AAz=(13/35+7/10*fiiy+1/3*fiiy^2+6/5*(rz/L)^2)/(1+fiiy)^2;
    AAy=(13/35+7/10*fiiz+1/3*fiiz^2+6/5*(ry/L)^2)/(1+fiiz)^2;
    
    BBz=(9/70+3/10*fiiy+1/6*fiiy^2-6/5*(rz/L)^2)/(1+fiiy)^2;
    BBy=(9/70+3/10*fiiz+1/6*fiiz^2-6/5*(ry/L)^2)/(1+fiiz)^2;
    
    CCz=(11/210+11/120*fiiy+1/24*fiiy^2+(1/10-1/2*fiiy)*(rz/L)^2)*L/(1+fiiy)^2;
    CCy=(11/210+11/120*fiiz+1/24*fiiz^2+(1/10-1/2*fiiz)*(ry/L)^2)*L/(1+fiiz)^2;
    
    DDz=(13/420+3/40*fiiy+1/24*fiiy^2-(1/10-1/2*fiiy)*(rz/L)^2)*L/(1+fiiy)^2;
    DDy=(13/420+3/40*fiiz+1/24*fiiz^2-(1/10-1/2*fiiz)*(ry/L)^2)*L/(1+fiiz)^2;
    
    EEz=(1/105+1/60*fiiy+1/120*fiiy^2+(2/15+1/6*fiiy+1/3*fiiy^2)*(rz/L)^2)*L^2/(1+fiiy)^2;
    EEy=(1/105+1/60*fiiz+1/120*fiiz^2+(2/15+1/6*fiiz+1/3*fiiz^2)*(ry/L)^2)*L^2/(1+fiiz)^2;
    
    FFz=-(1/140+1/60*fiiy+1/120*fiiy^2+(1/30+1/6*fiiy-1/6*fiiy^2)*(rz/L)^2)*L^2/(1+fiiy)^2;
    FFy=-(1/140+1/60*fiiz+1/120*fiiz^2+(1/30+1/6*fiiz-1/6*fiiz^2)*(ry/L)^2)*L^2/(1+fiiz)^2;
    
    % consistency
    Mtmp=zeros(12,12);
    % Column 1 + symmetry terms
    Mtmp(1,1)=1/3;    
    Mtmp(7,1)=1/6;    Mtmp(1,7)=1/6;  
    % Column 2 + symmetry terms
    Mtmp(2,2)=AAz;    
    Mtmp(6,2)=CCz;    Mtmp(2,6)=CCz; 
    Mtmp(8,2)=BBz;    Mtmp(2,8)=BBz;
    Mtmp(12,2)=-DDz;  Mtmp(2,12)=-DDz;
    % Column 3 + symmetry terms
    Mtmp(3,3)=AAy;    
    Mtmp(5,3)=-CCy;   Mtmp(3,5)=-CCy; 
    Mtmp(9,3)=BBy;    Mtmp(3,9)=BBy;
    Mtmp(11,3)=DDy;   Mtmp(3,11)=DDy;
    % Column 4 + symmetry terms
    Mtmp(4,4)=Jx/(3*A);    
    Mtmp(10,4)=Jx/(6*A);   Mtmp(4,10)=Jx/(6*A);   
    % Column 5 + symmetry terms
    Mtmp(5,5)=EEy;    
    Mtmp(9,5)=-DDy;   Mtmp(5,9)=-DDy; 
    Mtmp(11,5)=FFy;   Mtmp(5,11)=FFy;
    % Column 6 + symmetry terms
    Mtmp(6,6)=EEz;    
    Mtmp(8,6)=DDz;    Mtmp(6,8)=DDz; 
    Mtmp(12,6)=FFz;   Mtmp(6,12)=FFz;
    % Column 7
    Mtmp(7,7)=1/3;  
    % Column 8 + symmetry terms
    Mtmp(8,8)=AAz;    
    Mtmp(12,8)=-CCz;   Mtmp(8,12)=-CCz;
    % Column 9 + symmetry terms
    Mtmp(9,9)=AAy;    
    Mtmp(11,9)=CCy;   Mtmp(9,11)=CCy;
    % Column 10
    Mtmp(10,10)=Jx/(3*A); 
    % Column 11
    Mtmp(11,11)=EEy; 
    % Column 12
    Mtmp(12,12)=EEz; 
    
    % Final matrix
    M=rho*A*L*Mtmp;
    
    
else
    
    % Lumped
    M=rho*A*L/2*[ 1   0   0   0   0   0   0   0   0   0   0   0 
                  0   1   0   0   0   0   0   0   0   0   0   0
                  0   0   1   0   0   0   0   0   0   0   0   0
                  0   0   0   0   0   0   0   0   0   0   0   0
                  0   0   0   0   0   0   0   0   0   0   0   0
                  0   0   0   0   0   0   0   0   0   0   0   0
                  0   0   0   0   0   0   1   0   0   0   0   0
                  0   0   0   0   0   0   0   1   0   0   0   0
                  0   0   0   0   0   0   0   0   1   0   0   0
                  0   0   0   0   0   0   0   0   0   0   0   0
                  0   0   0   0   0   0   0   0   0   0   0   0
                  0   0   0   0   0   0   0   0   0   0   0   0];
    
end

% Transformation matrix
T=T3D(theta,L,coord);

% Global stiffness matrix
MG=T'*M*T;

% Preallocating global mass matrix for whole model
M3=zeros(DOFs); %DOFS

% 6 = number of node dofs
for ii=1:6 %rows
    for jj=1:6 % columns
        
         % upper-left
         M3(Igdof+(ii-1),Igdof+(jj-1))=MG(ii,jj);
         % upper-right
         M3(Igdof+(ii-1),Jgdof+(jj-1))=MG(ii,jj+6);
         % lower-left
         M3(Jgdof+(ii-1),Igdof+(jj-1))=MG(ii+6,jj);
         % lower-right
         M3(Jgdof+(ii-1),Jgdof+(jj-1))=MG(ii+6,jj+6);

    end
end



function [G3]=Gyro3D(A,E,rho,I,L,nuxy,shear,Igdof,Jgdof,DOFs,coord);
%========================================================================
% 3D beam element gyroscopic unsymmetrical damping matrix:
% [M3]=Mass3D(A,E,rho,I,L,nuxy,shear,Igdof,Jgdof,DOFs,coord)
%------------------------------------------------------------------
%
% Parameters:
% -----------
%
% - Spin  = Angular velocity around local x axis
% - A     = Element cross-sectional area
% - E     = Material elastic modulus
% - rho   = Material density
% - I     = Beam moment of inertia
% - L     = Beam length
% - nuxy  = Poisson's ratio
% - shear = shear strain coefficient
% - Igdof = beam I node global dof number (=UX global)
% - Jgdof = beam J node global dof number (=UX global)
% - DOFs  = number of dofs
% - theta = beam orientation around x axis
% - coord = Element node coordinates
% 
% Written by Jussi Sopanen, LTY/IMVe, 2004
% ========================================================================


% cross-section radius of gyration
% Assuming that cross-section is circular i.e. Iyy=Izz and sheary=shearz
r=sqrt(I/A); % radius of gyration

% If shear =0 shear strain is not taken into account
if shear==0
    fii=0;
else
    G=E/(2*(1+nuxy)); % glide modulus
    fii=12*E*I/(G*(A*shear)*L^2);
end
% fii=0

% forming gyro matrix

% Consistency mass matrix cells  (ANSYS Theory PIPE16)
gg=(6/5*r^2)/(L^2*(1+fii)^2);
hh= -(1/10-1/2*fii)*r^2/(L*(1+fii)^2);
ii= (2/15+1/6*fii+1/3*fii^2)*r^2/(1+fii)^2;
jj=-(1/30+1/6*fii-1/6*fii^2)*r^2/(1+fii)^2;

% consistency
Gtmp=zeros(12,12);
% Column 2 + antisymmetry terms
Gtmp(3,2)=-gg;    Gtmp(2,3)=gg;
Gtmp(5,2)=-hh;    Gtmp(2,5)=hh; 
Gtmp(9,2)=gg;     Gtmp(2,9)=-gg; 
Gtmp(11,2)=-hh;   Gtmp(2,11)=hh;
% Column 3 + antisymmetry terms
Gtmp(6,3)=-hh;    Gtmp(3,6)=hh;
Gtmp(8,3)=-gg;    Gtmp(3,8)=gg; 
Gtmp(12,3)=-hh;   Gtmp(3,12)=hh; 
% Column 5 + antisymmetry terms
Gtmp(6,5)=-ii;    Gtmp(5,6)=ii;
Gtmp(8,5)=-hh;    Gtmp(5,8)=hh; 
Gtmp(12,5)=-jj;   Gtmp(5,12)=jj; 
% Column 6 + antisymmetry terms
Gtmp(9,6)=-hh;    Gtmp(6,9)=hh;
Gtmp(11,6)=jj;    Gtmp(6,11)=-jj; 
% Column 8 + antisymmetry terms
Gtmp(9,8)=-gg;    Gtmp(8,9)=gg;
Gtmp(11,8)=hh;    Gtmp(8,11)=-hh; 
% Column 9 + antisymmetry terms
Gtmp(12,9)=hh;    Gtmp(9,12)=-hh;
% Column 11 + antisymmetry terms
Gtmp(12,11)=-ii;    Gtmp(11,12)=ii;


% Final matrix
G1=2*rho*A*L*Gtmp; % at the end multiplied by angular velocity Spin = Omega ergo G1=2*Spin*rho*A*L*Gtmp;

theta=0;
% Transformation matrix
T=T3D(theta,L,coord);

% Global gyro matrix
GG=T'*G1*T;

% Preallocating global gyro matrix for whole model
G3=zeros(DOFs); %DOFS

% 6 = number of node dofs
for ii=1:6 % rows
    for jj=1:6 % columns
        
         % upper-left
         G3(Igdof+(ii-1),Igdof+(jj-1))=GG(ii,jj);
         % upper-right
         G3(Igdof+(ii-1),Jgdof+(jj-1))=GG(ii,jj+6);
         % lower-left
         G3(Jgdof+(ii-1),Igdof+(jj-1))=GG(ii+6,jj);
         % lower-right
         G3(Jgdof+(ii-1),Jgdof+(jj-1))=GG(ii+6,jj+6);

    end
end


function [T]=T3D(theta,L,coord)
% Function T3D calculates rational matrix for beam elements


T=zeros(12);
d=0.0001*L;

X1=coord(1);
X2=coord(2);
Y1=coord(3);
Y2=coord(4);
Z1=coord(5);
Z2=coord(6);

% Element projected length in xy plane
Lxy=sqrt(L^2-(Z2-Z1)^2);

if Lxy>d
    S1=(Y2-Y1)/Lxy;
    C1=(X2-X1)/Lxy;
else
    S1=0;
    C1=1;
end

S2=(Z2-Z1)/L;
S3=sin(theta);
C2=Lxy/L;
C3=cos(theta);

Tpart=[  C1*C2             S1*C2            S2
       (-C1*S2*S3-S1*C3) (-S1*S2*S3+C1*C3)  S3*C2
       (-C1*S2*C3+S1*S3) (-S1*S2*C3-C1*S3)  C3*C2];
   
T(1:3,1:3)=Tpart;
T(4:6,4:6)=Tpart;
T(7:9,7:9)=Tpart;
T(10:12,10:12)=Tpart;









