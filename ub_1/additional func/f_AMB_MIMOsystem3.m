function rotor = f_AMB_MIMOsystem3(Inp,varargin)
% Description: functions combines all the information avalilable and makes
% a reduced ss 'rotor' structure for simulations and controller synthesis
% 
%  PREFORMATTED
%  TEXT
% 
% Inputs:
%           Inp - structure with geometry of the rotor and information
%           abouts bearings and sensors
%           prm.Nmodes - how many modes to include in the final model (2 modes
%           = 1 flexible in x and y directions (by default = 4)
% Outputs:
%           rotor - structure with a information about reduced ss model of
%           the system
% if nargin>1
%     Nmodes = varargin{1};
% else
%     Nmodes = 4;
% end
% Coded by Alexadner Smirnov
%

% 2019.05.07: GYROSCOPIC EFFECT DOES NOT WORK IF NMODES OVER 2! CHECK GYRO MATRIX & "DIRTY FIX" SECTION BELOW FOR FAULTS!
pars = inputParser;
pars.KeepUnmatched = true;
pars.StructExpand = true;

addOptional(pars,'Nmodes',4, @(x) isnumeric(x) && isreal(x));
addOptional(pars,'ModeFrequencies',NaN, @(x) isnumeric(x) && isreal(x));
addOptional(pars,'ModeUncertaintyPercent',NaN, @(x) isnumeric(x) && isreal(x));
addOptional(pars,'GyroscopicScale',NaN, @(x) isnumeric(x) && isreal(x) && ismatrix(x));
addOptional(pars,'SensorModesScale',NaN, @(x) isnumeric(x) && isreal(x) && ismatrix(x));
addOptional(pars,'ActuatorModesScale',NaN, @(x) isnumeric(x) && isreal(x) && ismatrix(x));
parse(pars,varargin{:});
prm = pars.Results;

% Calculate matrices
[Matrx]=f_Rotor_Matrix(Inp); % Genetates full matrices (in global dofs)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrx structure include following fields:
% Matrx.K  = Stiffness matrix
% Matrx.M  = Mass matrix
% Matrx.G  = Gyromatrix
% Matrx.F2 = unbalance
% Matrx.F3 = unbalance

Kc = Matrx.Kc; 
Cc = Matrx.Cc;
%Kc2=Matrx.Kc; %Kc2 is the same as Kc when symmetric bearing

% Here are made manually rigid body modes and only for 4 dofs/node
% (Y,Z,theta_y,theta_z) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(Inp.Node,1)
    if i==1; 
        v1x     = [1 0 0 0]';                           % translation x
        v1y     = [0 1 0 0]';                           % translation y
        vrot_x  = [0 -Inp.Node(i,2) 1 0]';              % rotation x
        vrot_y  = [Inp.Node(i,2) 0 0 1]';               % rotation y
    else 
        v1x     = [v1x; 1; 0; 0; 0];                    % translation x
        v1y     = [v1y; 0; 1; 0; 0];                    % translation y
        vrot_x  = [vrot_x; 0; -Inp.Node(i,2); 1; 0];    % rotation x
        vrot_y  = [vrot_y; Inp.Node(i,2); 0; 0; 1];     % rotation y
    end;
end
    
% Center of mass (Lantto 1997)
% Mc = global mass matrix with constraints take into account
zm=abs(vrot_y'*Matrx.Mc*v1x/(v1x'*Matrx.Mc*v1x));

%rotation modes
v2x=vrot_x+zm*v1y;
v2y=vrot_y-zm*v1x;

%Rigid body modes
FiiR=[v1x v1y v2x v2y];

%% TEEMU: see also f_ModalData() & f_Campbell() ??
% Eigenvalues
[FiiN, d]=eig(Kc,Matrx.Mc);
% Matriisien jarjestely taajuuden mukaan.
[Freq, idx]=sort(sqrt(diag(d))/(2*pi));
FiiN=FiiN(:,idx); clear idx d;

% Muotojen maaran valinta (Nmuoto alinta ominaismuotoa)
% FiiNn=FiiN(:,1:prm.Nmodes);
FiiNn=FiiN(:,5:4+prm.Nmodes); % Select modes

FiiNn_max1=max(abs(FiiNn(1:4:end,:))); % Something wrong%%%%%%%%%%%%%%%%%%%%%%%%
FiiNn_max2=max(abs(FiiNn(2:4:end,:)));
FiiNn_max=max([FiiNn_max1; FiiNn_max2]);

for zz=1:length(FiiNn_max)     %Nb of flex. modes
    FiiNn(:,zz)=1/FiiNn_max(zz)*FiiNn(:,zz);
end

FiiNormN = [FiiR FiiNn];

% % Massanormeeraus This should be used, but this needs to be fixed%%%%%%%%%%%%%%%%%%%%%%%
% for k=1:prm.Nmodes
%     FiiNormN(:,k)=FiiNn(:,k)/sqrt(FiiNn(:,k)'*Matrx.Mc*FiiNn(:,k));
% end

% Moodaaliset matriisit

Mm=FiiNormN'*Matrx.Mc*FiiNormN;
%Km=FiiNormN'*Kc2*FiiNormN + FiiNormN'*(Kc-Kc2)*FiiNormN; % jälkimmäinen termi sisältää laakerin ristitermit
Km=FiiNormN'*Kc*FiiNormN;

% Gyromatriisi kerrotaan pyörimisnopeudella
Gm=FiiNormN'*Matrx.Gc*FiiNormN;
% vaimennusmatriisi + gyromatriisi
Cm=FiiNormN'*Cc*FiiNormN;

% Clean matricies from numerical errors
Mm = diag(diag(Mm));
Km = blkdiag(zeros(4),diag(diag(Km(5:end,5:end))));
Gm(abs(Gm)<1e-9) = 0; % This needs to checked, recommemded to not use before this is verified
% Dirty fix for the gyroscopic matrix (07.05.2019 OR DIRTY BUG GENERATOR?!?)
Gm = triu(abs(Gm),1)+tril(-abs(Gm),-1); % This should not be used, seems incorrect%%%%
if ~isnan(prm.GyroscopicScale) % This should be be needed %%%%%%%%%%%%%%%%%%%%%%
    Gm = prm.GyroscopicScale'*Gm*prm.GyroscopicScale;
end


% If 'Inp.ModalDamping' exists then Dm==Cm
% Introduce default damping factor
zeta_d2=0.002;
df_mat=(zeros(1,4));
if prm.Nmodes>0 % This should be checked, I dont know what this loop makes%%%%%%%%%%%%%%%
    for ii=1:prm.Nmodes/2;
        if all(~isnan(prm.ModeFrequencies)) && (numel(prm.ModeFrequencies)>=ii)
            w0=prm.ModeFrequencies(ii);
            w_old=2*pi*Freq(1+(2*(ii-1))+4);
            Km(1+(2*(ii-1))+4,1+(2*(ii-1))+4) = Km(1+(2*(ii-1))+4,1+(2*(ii-1))+4)*w0^2/w_old^2;
            Km(2+(2*(ii-1))+4,2+(2*(ii-1))+4) = Km(2+(2*(ii-1))+4,2+(2*(ii-1))+4)*w0^2/w_old^2;
        else
            w0=2*pi*Freq(1+(2*(ii-1))+4);
        end
        if isfield(Inp,'ModalDamping') % Acutally this is alrady included in the matrix calculations and Dm=Cm can be used
%             if 2*ii<numel(Inp.ModalDamping)/2
            if 2*ii<=numel(Inp.ModalDamping)/2
                df_mat=[df_mat 2*Inp.ModalDamping(2*ii,2)*[1/w0 1/w0]];
            else
                df_mat=[df_mat 2*zeta_d2*[1/w0 1/w0]]; %Use estimate
            end
        else
            df_mat=[df_mat 2*zeta_d2*[1/w0 1/w0]]; %Use estimate
        end

    end
end
Dm = diag(df_mat)*Km; % Dm = Cm can be used with 'Inp.ModalDamping' filled

if ~isnan(prm.ModeUncertaintyPercent)
    Kmu = umat(Km);
    for ii=1:numel(prm.ModeUncertaintyPercent)
        indx = 2*(ii-1)+4;
        if prm.ModeUncertaintyPercent(ii)>0
            Munc = ucomplex(['Flex' num2str(ii)],Km(indx+1,indx+1),'Percentage',...
                prm.ModeUncertaintyPercent(ii));
            Kmu(indx+1,indx+1) = Munc; % for x direction
            Kmu(indx+2,indx+2) = Munc; % for y direction
        end
    end
else
    Kmu = Km;
end

omega = 0; % NOTE this should be updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenvalue problem A*X=lambda*X
zeroX=zeros(size(Mm));  eyeX=eye(size(Mm));
Asy=[   zeroX     	eyeX
        -Mm\Km      -Mm\(Dm+omega*Gm)];
    
omega_u = ureal('speed',0.5*Inp.NominalSpeed,'Range',[0 Inp.NominalSpeed]);
% Asy_u=[ zeroX     	eyeX
%         -Mm\Km      -Mm\(Dm+omega_u*Gm)];
Asy_u=[ zeroX     	eyeX
        -Mm\Kmu     -Mm\(Dm+omega_u*Gm)];
% Asy_u=[ zeroX     	eyeX
%         -Mm\Km      -Mm\Dm-Mm\Gm*omega_u];
    
mat_dat.Ak = [	zeroX     	eyeX
                -Mm\Km      zeroX];
mat_dat.Ad = [	zeroX     	zeroX
                zeroX      -Mm\Dm];
mat_dat.Ag = [	zeroX     	zeroX
                zeroX      -Mm\Gm]; 
    
num_amb = 0; % Number of nodes related to actuators
num_act = 0; % Number of actuators 2=xy direction = 1 amb
num_sen = 0; % Number of sensors
inp_names = {}; % Input names for ss model
inp_grp = {};   % Input groups
out_grp = {};   % Output groups
out_names = {}; % Output names for ss model
nda = [];       % Actuator nodes
nds = [];       % Sensor nodes
Kx = [];        % Position stiffness
Ki = [];        % Current stiffness
ndi = [];       % Input nodes
ndo = [];       % Output nodes



for ii=1:length(Inp.AMBs.Actuator)
    if isfield(inp_grp,'amb')
        inp_grp.amb=[inp_grp.amb 2*ii-1 2*ii];
    else
        inp_grp.amb=[2*ii-1 2*ii];
    end
    num_amb = num_amb+1;
    num_act = num_act+2;
    nda = [nda Inp.AMBs.Actuator(ii).Inode]; 
    Kx = blkdiag(Kx,Inp.AMBs.Actuator(ii).Kx*eye(2));
    Ki = blkdiag(Ki,Inp.AMBs.Actuator(ii).Ki*eye(2));
    % If there are names copy them
    if isfield(Inp.AMBs.Actuator(ii),'name')
        inp_names{2*ii-1}=[Inp.AMBs.Actuator(ii).name ', x']; 
        inp_names{2*ii}=[Inp.AMBs.Actuator(ii).name ', y']; 
    else
        inp_names{2*ii-1}=[num2str(ii) 'x']; 
        inp_names{2*ii}=[num2str(ii) 'y']; 
    end
end

% Iterate througth the available sensors
if isfield(Inp.AMBs, 'Sensors')
    for ii=1:length(Inp.AMBs.Sensors)
        if isfield(out_grp,'amb')
            out_grp.amb=[out_grp.amb 2*ii-1 2*ii];
        else
            out_grp.amb=[2*ii-1 2*ii];
        end

        num_sen = num_sen+1;
        nds = [nds Inp.AMBs.Sensors(ii).Inode]; 
        if isfield(Inp.AMBs.Sensors(ii),'name')
            out_names{2*ii-1}=[Inp.AMBs.Sensors(ii).name ', x']; 
            out_names{2*ii}=[Inp.AMBs.Sensors(ii).name ', y']; 
        else
            out_names{2*ii-1}=[num2str(ii) 'x']; 
            out_names{2*ii}=[num2str(ii) 'y']; 
        end
    end
end

% Iterate througth the additonal inputs
if isfield(Inp.AMBs, 'Inputs')
    for ii=1:length(Inp.AMBs.Inputs)
        if isfield(inp_grp,'inp')
            inp_grp.inp=[inp_grp.inp 2*ii-1 2*ii];
        else
            inp_grp.inp=[2*ii-1 2*ii];
        end
        if isfield(inp_grp,'amb')
            inp_grp.amb = inp_grp.amb+2;
        end
        ndi = [ndi Inp.AMBs.Inputs(ii).Inode]; 
        if isfield(Inp.AMBs.Inputs(ii),'name')
            inp_names = [[Inp.AMBs.Inputs(ii).name ', x'] inp_names]; 
            inp_names = [[Inp.AMBs.Inputs(ii).name ', y'] inp_names]; 
        else
            inp_names = [['Out ' num2str(ii) 'x'] inp_names]; 
            inp_names = [['Out ' num2str(ii) 'y'] inp_names]; 
        end
    end
end

% Iterate througth the additonal inputs
if isfield(Inp.AMBs, 'Outputs')
    for ii=1:length(Inp.AMBs.Outputs)
        if isfield(out_grp,'out')
            out_grp.out=[out_grp.out 2*ii-1 2*ii];
        else
            out_grp.out=[2*ii-1 2*ii];
        end
        if isfield(out_grp,'amb')
            out_grp.amb = out_grp.amb+2;
        end
        ndo = [ndo Inp.AMBs.Outputs(ii).Inode]; 
        if isfield(Inp.AMBs.Outputs(ii),'name')
            out_names = [[Inp.AMBs.Outputs(ii).name ', x'] out_names]; 
            out_names = [[Inp.AMBs.Outputs(ii).name ', y'] out_names]; 
        else
            out_names = [['Out ' num2str(ii) 'x'] out_names];
            out_names = [['Out ' num2str(ii) 'y'] out_names];
        end
    end
end

Sai = sparse([(nda-1)*4+1 (nda-1)*4+2], ...
            [1:2:num_amb*2 2:2:num_amb*2],ones(1,num_amb*2),...
            size(FiiNormN,1),num_amb*2);
Keye = speye(num_act);
Keyem = FiiNormN'*Sai*Keye;
Keyem(abs(Keyem)<1e-9) = 0; % clean matrix

Ss = sparse([1:2:num_sen*2 2:2:num_sen*2],[(nds-1)*4+1 (nds-1)*4+2],...
            ones(1,num_sen*2), num_sen*2, size(FiiNormN,1));
Bamb = [zeros(4+prm.Nmodes,num_act); Mm\Keyem];
C0b = [(Sai')*FiiNormN zeros(2*num_amb,4+prm.Nmodes)];
C0b(abs(C0b)<1e-9) = 0; % clean matrix

if ~isnan(prm.ActuatorModesScale)
    indx_r = [1 1 2 2];
    indx_c = reshape([1:prm.Nmodes/2;1:prm.Nmodes/2],1,prm.Nmodes);
    ActuatorScale = prm.ActuatorModesScale(indx_r,indx_c);
    C0b(:,4+(1:prm.Nmodes)) = ActuatorScale(:,1:prm.Nmodes).*C0b(:,4+(1:prm.Nmodes));
    Bamb(end-prm.Nmodes+1:end,:) = ActuatorScale(:,1:prm.Nmodes)'.*Bamb(end-prm.Nmodes+1:end,:);
end


A0 = Asy + Bamb*Kx*C0b;
A0_u = Asy_u + Bamb*Kx*C0b;
% A0_u = mat_dat.Ak+mat_dat.Ad+omega_u*mat_dat.Ag+Bamb*Kx*C0b;

B0 = Bamb*Ki;

C0 = full([Ss*FiiNormN zeros(2*num_sen,4+prm.Nmodes)]);
C0(abs(C0)<1e-9) = 0; % clean matrix

if ~isnan(prm.SensorModesScale)
    indx_r = [1 1 2 2];
    indx_c = reshape([1:prm.Nmodes/2;1:prm.Nmodes/2],1,prm.Nmodes);
    SensorScale = prm.SensorModesScale(indx_r,indx_c);
    C0(:,4+(1:prm.Nmodes)) = SensorScale(:,1:prm.Nmodes).*C0(:,4+(1:prm.Nmodes));
end

D0 = zeros(2*num_sen,num_act);

mat_dat.Bb = Bamb;
mat_dat.Cb = C0b;
mat_dat.Cs = C0;

% Additonal inputs
if isfield(Inp.AMBs, 'Inputs')
    num_inp = length(Inp.AMBs.Inputs);
    Sinp = sparse([(ndi-1)*4+1 (ndi-1)*4+2], ...
                [1:2:num_inp*2 2:2:num_inp*2],ones(1,num_inp*2),...
                size(FiiNormN,1),num_inp*2);
    Binp = [zeros(4+prm.Nmodes,2*num_inp); Mm\(FiiNormN'*Sinp)];
end

% Additonal outputs
if isfield(Inp.AMBs, 'Outputs')
    num_out = length(Inp.AMBs.Outputs);
    Sout = sparse([(ndo-1)*4+1 (ndo-1)*4+2], ...
                [1:2:num_out*2 2:2:num_out*2],ones(1,num_out*2),...
                size(FiiNormN,1),num_out*2);
    Cout = [Sout'*FiiNormN zeros(2*num_out,4+prm.Nmodes)];
end

rotor.m.Mm = Mm;
rotor.m.Km = Km;
rotor.m.Dm = Dm;
rotor.m.Gm = Gm;
rotor.m.zm = zm;
rotor.m.ModeS = full(FiiNormN);
rotor.Sai = full(Sai);
rotor.infb.prm.Nmodes = prm.Nmodes;
rotor.infb.num_act = num_act;
rotor.infb.num_sen = num_sen;
    
if isempty(Inp.Disp)
    Matrx.Fg(1:3:end) = [];
end
rotor.fgv=FiiNormN'*Matrx.Fg;


if isfield(Inp.AMBs, 'Inputs')
    B0 = [full(Binp) B0];
    D0 = [D0 zeros(size(D0,1),size(Binp,2))];
end
if isfield(Inp.AMBs, 'Outputs')
    C0 = [full(Cout); C0];
    D0 = [D0; zeros(size(Cout,1),size(D0,2))];
end

rotor.C0b = full(C0b);
rotor.A0 = A0;
rotor.B0 = B0;
rotor.C0 = C0;
rotor.D0 = D0;
rotor.Bamb = Bamb;

% Make a ss representation
rotor.sys = ss(A0,B0,C0,D0);
rotor.sys.InputName = inp_names;
rotor.sys.OutputName = out_names;

rotor.sys.OutputGroup = out_grp;
rotor.sys.InputGroup = inp_grp;
% 
rotor.sysu = uss(A0_u,B0,C0,D0);
rotor.sysu.InputName = inp_names;
rotor.sysu.OutputName = out_names;

rotor.sysu.OutputGroup = out_grp;
rotor.sysu.InputGroup = inp_grp;

rotor.m.mat_dat = mat_dat;
% rotor = [];

end
