function [t,y,Matrx,Matrx_g]=f_TransientSRB(Inp,Request)

% coded by Behnam Ghalamchi
% f_Gravity
% F_system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Rotor Matrices
[Matrx]=f_Rotor_Matrix(Inp);
% Calculate rotor matrices using lumped mass matrices
Inp.lumpm=1;
[Matrx_g]=f_Rotor_Matrix(Inp);
%%%%% Calculation of Gravity Force
Matrx.Fgrav=Matrx.Fg;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Define Constrain number of degrees of freedom
% Matrx.Rdofs=length(Matrx.Kc(:,1));    %number of Rotor dofs
% Matrx.Sdofs=0;                        %number of Support dofs
Matrx.Tdofs=length(Matrx.Kc(:,1));     %number of total dofs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Total matrix include supports
% Matrx.M=[Matrx.Mc                zeros(Matrx.Rdofs,Matrx.Sdofs)
%        zeros(Matrx.Sdofs,Matrx.Rdofs)   eye(Matrx.Sdofs)*Inp.ms  ];
% 
% Matrx.C=[Inp.Omega*Matrx.Gc+Matrx.Cc   zeros(Matrx.Rdofs,Matrx.Sdofs)
%        zeros(Matrx.Sdofs,Matrx.Rdofs)    eye(Matrx.Sdofs)*Inp.Cs    ];
% 
% Matrx.K=[Matrx.Kc         zeros(Matrx.Rdofs,Matrx.Sdofs)
%        zeros(Matrx.Sdofs,Matrx.Rdofs)    eye(Matrx.Sdofs)*Inp.Ks    ];
   
%%%%%%%%%%%%%%% Eigenvalue solution
[FiiN d]=eig(Matrx.Kc,Matrx.Mc);
[Freq, idx]=sort(sqrt(diag(d))/(2*pi)); % Matrices arrangement of frequency.
FiiN=FiiN(:,idx); clear idx d; % Sort shapes according to frequency
FiiNn=FiiN(:,1:Request.Transient.NroModes); % Selection of the number of modes (N form of the lowest eigenmode)

%%%%%%%%%%% Normalaized mass%%%%%%%%%%%
for k=1:Request.Transient.NroModes
    FiiNormN(:,k)=FiiNn(:,k)/sqrt(FiiNn(:,k)'*Matrx.Mc*FiiNn(:,k));
end 
Matrx.FiiNormN=FiiNormN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transient Analysis
tstep=Inp.tstep;
tend=Inp.tend;
t=0:tstep:tend;
err=1e-4;
options = odeset('RelTol',err,'MaxStep',5e-4,'Stats','on');
len=2*Matrx.Tdofs;
y0(1,:)=zeros(2*Matrx.Tdofs,1); 

h = waitbar(0,'Computing Transient Response...');
Inp.h=h;
tic
% [t,y]=ode15s(@EOM_transient,t,y0,options,Matrx,tend,Inp,Request);
[t,y]=ode15s(@f_system,t,y0,options,Matrx,Inp,Request);
toc
delete(h)
full=len/2;
figure; plot(t,y(:,17)); grid on; title('transient-ode15s')
xlabel('time [sec]','Fontsize',12)
% hold on
% plot(t,y(:,45)); grid on; title('transient-ode15s')
ylabel('Y-Displacement [\mum]','Fontsize',12)
figure(4); plot(y(:,17),y(:,18)); grid on; title('orbit')
xlabel('z-Displacement [\mum] ','Fontsize',12)
ylabel('Y-Displacement [\mum]','Fontsize',12)
figure(5); plot(t,y(:,18)); grid on; title('transient-ode15s')
xlabel('time [sec]','Fontsize',12)
ylabel('Z-Displacement [\mum]','Fontsize',12)
end