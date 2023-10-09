function [FB1,FB2]=f_Bearing_Radial_Loads(Inp,grav)
% Function calculates bearing forces for rotor supported by two bearings
% Rotor needs to be along X axis!!!!!!
% INPUT: 
% Node = nodematrix
% Node_B1 = node of bearing 1
% Node_B2 = node of bearing 2
% Mg = centralized mass matrix
% grav = gravitation vector
% nnum2gdof = index matrix

% Calculation of Rotor Matrices------------------------------------------
Inp.lumpm=1; % using lumped mass matrix
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


mg=diag(Matrx.Mc);
disp(' ')
disp(['Rotor mass : ' num2str(sum(mg)/3,4) ' kg'] )
% Force vector from gravitation-------------------------------------------------

Node=Inp.Node;
nnum2cgdof=Matrx.nnum2cgdof;

Fgrav=zeros(size(mg));
Mgrav=zeros(size(mg));
for ii=1:size(Node,1)
    
    % removing mass of bearing housings
    % if Node(ii,1) is equal to Jnode of bearing, it is neglected
    % from calculation
    if ~(Node(ii,1)==Inp.Bearing(1).Jnode | Node(ii,1)==Inp.Bearing(2).Jnode)
        
        % getting global DOFs
        gdofX=nnum2cgdof(Node(ii,1),1); 
        gdofY=nnum2cgdof(Node(ii,1),2); 
        gdofZ=nnum2cgdof(Node(ii,1),3); 
        
        % Mass matrix of the node
        mtmp=zeros(3);
        if ~gdofX==0; mtmp(1,1)=mg(gdofX); end;
        if ~gdofY==0; mtmp(2,2)=mg(gdofY); end;
        if ~gdofZ==0; mtmp(3,3)=mg(gdofZ); end;
        
        % Gravitation force of the node
        Ftmp=mtmp*grav;
        
        % Placing node forces into correct coordinates
        if ~gdofX==0; Fgrav(gdofX)=Ftmp(1);  end;
        if ~gdofY==0; Fgrav(gdofY)=Ftmp(2);  end;
        if ~gdofZ==0; Fgrav(gdofZ)=Ftmp(3);  end;
        
        % Moment caused by gravitation M=F*x
        if ~gdofX==0; Mgrav(gdofX)=Ftmp(1)*Node(ii,2); end;
        if ~gdofY==0; Mgrav(gdofY)=Ftmp(2)*Node(ii,2); end;
        if ~gdofX==0; Mgrav(gdofZ)=Ftmp(3)*Node(ii,2); end;
    end
end

Node_B2=Inp.Bearing(2).Inode;
Node_B1=Inp.Bearing(1).Inode;
% Bearing forces can be solved used following equations:
%    FB1+FB2+sum(Fgrav)=0                                    (1)
%    FB1*Node(Node_B1,2)+FB2*Node(Node_B2,2)+sum(Mgrav)=0;   (2)
% (1) ==> 
%    FB2=(-sum(Fgrav)-FB1)
% sij. (2) ==>
%    FB1*Node(Node_B1,2)+(-sum(Fgrav)-FB1)*Node(Node_B2,2)+sum(Mgrav)=0; ==>
%    FB1*(Node(Node_B1,2)-Node(Node_B2,2)) -sum(Fgrav)*Node(Node_B2,2) + sum(Mgrav)=0 ==>

FB1= (+sum(Fgrav)*Node(Node_B2,2) - sum(Mgrav))/(Node(Node_B1,2)-Node(Node_B2,2));
FB2=(-sum(Fgrav)-FB1);
