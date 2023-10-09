function [m]=f_ModalData(Inp)
% Function calculates modal data for the rotor model
[Matrx]=f_Rotor_Matrix(Inp);m=Matrx;% Calculate matrices
% Matrx structure include following fields:
% Matrx.K= Stiffness matrix Matrx.M= Mass matrix Matrx.G= Gyromatrix    
% Matrx.F2= unbalance       Matrx.F3= unbalance

%% Create rigid body modes
for i=1:Inp.RotorDOF/4                                      % size(Inp.Node,1)
    if i==1; 
        v1y=[1 0 0 0]';                                 % translation y
        v1z=[0 1 0 0]';                                 % translation z
        vrot_y=[0 -Inp.Node(i,2) 1 0]';                 % rotation y
        vrot_z=[Inp.Node(i,2) 0 0 1]';                  % rotation z
    else 
        v1y=[v1y;1; 0; 0; 0];                           % translation y
        v1z=[v1z; 0; 1; 0; 0];                          % translation z
        vrot_y=[vrot_y;0; -Inp.Node(i,2); 1; 0];        % rotation y
        vrot_z=[vrot_z;Inp.Node(i,2); 0; 0; 1];         % rotation z
    end;
end

zm=abs(vrot_z'*Matrx.Mc(1:Inp.RotorDOF,1:Inp.RotorDOF)*v1y/(v1y'*Matrx.Mc(1:Inp.RotorDOF,1:Inp.RotorDOF)*v1y));% Center of mass (Lantto 1997)
v2y=vrot_y+zm*v1z;v2z=vrot_z-zm*v1y;    % rotation modes
FiiR=[v1y v1z v2y v2z];                 % Rigid body modes
[FiiN,d]=eig(Matrx.Kc(1:Inp.RotorDOF,1:Inp.RotorDOF),Matrx.Mc(1:Inp.RotorDOF,1:Inp.RotorDOF));% Eigenvalues
[Freq,idx]=sort(sqrt(diag(d))/(2*pi));  % Sort according to frequency
for k=1:length(idx)                     % Sort modes according to frequency
    Fiiapu(:,k)=FiiN(:,idx(k));
end
FiiN=Fiiapu; clear Fiiapu idx k d
FiiNn=FiiN(:,5:4+Inp.Modes);% Select modes
%% Orthogonalization of modes and rotation according to main axes
% Solved eigenmodes in YX or ZX plane are not absolutely accurate
% Modes are mapped according to these YX or ZX planes. Quick and dirty method...
for k=1:Inp.Modes % going through the modes
    for i=1:4%size(Inp.Node,1) % going through the nodes
        % saving node displacements and rotations
        Ydisp=FiiNn(i*4-3,k);                   % y displacement
        Zdisp=FiiNn(i*4-2,k);                   % z displacement
        Yrot=FiiNn(i*4-1,k);                    % y rotation
        Zrot=FiiNn(i*4,k);                      % z rotation
        DispMag=sqrt( Zdisp^2 + Ydisp^2);       % magnitude of displacement
        RotMag=sqrt( Zrot^2 + Yrot^2);          % magnitude of displacement
        angle=atan2(Zdisp,Ydisp);               % phase angle of displacement -pi...pi
        angleR=atan2(Zrot,Yrot);                % phase angle of rotation -pi...pi
        % setting the displacements according to main axes
        if angle >= -pi/4 && angle <= pi/4      % closest displacement according to positive y axis
            FiiNn(i*4-3,k)=DispMag; 
            FiiNn(i*4-2,k)=0;
        elseif angle >= pi/4 && angle <= 3*pi/4 % closest displacement according to positive z axis
            FiiNn(i*4-3,k)=0;
            FiiNn(i*4-2,k)=DispMag;
        elseif angle >= -3*pi/4 && angle <= -pi/4 % closest displacement according to negative z axis
            FiiNn(i*4-3,k)=0;
            FiiNn(i*4-2,k)=-DispMag;
        elseif abs(angle) >= 3*pi/4             % closest displacement according to negative y axis
            FiiNn(i*4-3,k)=-DispMag;
            FiiNn(i*4-2,k)=0;
        end
        % asetetaan rotaatiot pääakseleiden suunt.
        if angleR >= -pi/4 && angleR <= pi/4    % closest rotation around positive y axis
            FiiNn(i*4-1,k)=RotMag; 
            FiiNn(i*4,k)=0;
        elseif angleR >= pi/4 && angleR <= 3*pi/4 % closest rotation around positive z axis
            FiiNn(i*4-1,k)=0;
            FiiNn(i*4,k)=RotMag;
        elseif angleR >= -3*pi/4 && angleR <= -pi/4 % closest rotation around negative z axis
            FiiNn(i*4-1,k)=0;
            FiiNn(i*4,k)=-RotMag;
        elseif abs(angleR) >= 3*pi/4            % closest rotation around negative y axis
            FiiNn(i*4-1,k)=-RotMag;
            FiiNn(i*4,k)=0;
        end           
    end
end
%% Modal matrices for output
m.ModeS=[FiiR FiiNn];           % Assemble mode shape matrix
m.Mm=m.ModeS'*Matrx.Mc(1:Inp.RotorDOF,1:Inp.RotorDOF)*m.ModeS;      m.fgm=m.ModeS'*m.Fg(1:Inp.RotorDOF);
m.Km=m.ModeS'*Matrx.Kc(1:Inp.RotorDOF,1:Inp.RotorDOF)*m.ModeS;      m.F2m=m.ModeS'*Matrx.F2(1:Inp.RotorDOF); 
m.Cm=m.ModeS'*(Matrx.Cc(1:Inp.RotorDOF,1:Inp.RotorDOF))*m.ModeS;    m.F3m=m.ModeS'*Matrx.F3(1:Inp.RotorDOF);
m.Gm=m.ModeS'*(Matrx.Gc(1:Inp.RotorDOF,1:Inp.RotorDOF))*m.ModeS; 
                 
