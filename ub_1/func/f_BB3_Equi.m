function [Eq, beta_i,beta_o, Fi, Fo, Fc,L3,L4,deltai,deltao,Ki,Ko,EQtmp11,EQtmp22,lambdaa,Mg]=f_BB3_Equi(X,B,FR)
% Function calculates ball bearing equilibrium equations

% Eq = Vector of Equations [(3.100); (3.101); (3.104); (3.105);], Ref [1] 
%
% X= Vector of known parameters [L_3, L_4, delta_i, delta_o]
% delta = vector of unknown bearing displacements [delta_t; delta_r]
% B = Structure of Bearing parameters
% FR = Bearing load vector [F_t, F_r]

% Equilibrium equations for high speed ball bearing------------------

% initial values for parameters to be iterated
delta_t=X(4*B.n+1); % axial displacement
delta_r=X(4*B.n+2); % radial displacement

j=0;
beta_i=zeros(length(B.n));
beta_o=zeros(length(B.n));
for j=1:B.n
    
L_3=X(4*(j-1)+1,1);     % Fig. 3.11, p. 94
L_4=X(4*(j-1)+2,1);     % 
delta_i=X(4*(j-1)+3,1);
delta_o=X(4*(j-1)+4,1);

    % attitude angle
    psi=2*pi*(j-1)/B.n; 
    % Axial distance between inner and outer curvature centers
    L_1=B.D*sin(B.beta_f)+delta_t; % Eq. (3.93), p. 95
    % Radial distance between inner and outer curvature centers
    L_2=B.D*cos(B.beta_f)+delta_r*cos(psi); % Eq. (3.94), p. 95

    % Sines and cosines of the contact angles C=cos, S=sin
    C_beta_o=L_4/(B.d*(B.f_o-0.5)+delta_o);
    S_beta_o=L_3/(B.d*(B.f_o-0.5)+delta_o);
    C_beta_i=(L_2-L_4)/(B.d*(B.f_i-0.5)+delta_i);
    S_beta_i=(L_1-L_3)/(B.d*(B.f_i-0.5)+delta_i);
    
     beta_i(j)=acos(C_beta_i); % Inner ring contact angle
     beta_o(j)=acos(C_beta_o); % Outer ring contact angle
%     if C_beta_o > 1  %Removes negative angles
%     beta_o(j)=0;
%     C_beta_o=1;
%     end

    % Equations that will be solved iteratively
    % see fig 3.11
    Eq(4*(j-1)+1,1)=L_4^2+L_3^2-(B.d*(B.f_o-0.5)+delta_o)^2;  %Eq. 3.100
    Eq(4*(j-1)+2,1)=(L_2-L_4)^2+(L_1-L_3)^2-((B.f_i-0.5)*B.d+delta_i)^2; %Eq. 3.101    
    % contact stiffnesses
    [K_i,K_o]=f_Hertz_K(C_beta_i,C_beta_o,B);
    % Normal forces-----
    delta_i=max([0 delta_i]);
    delta_o=max([0 delta_o]);
    if delta_i < 0
        disp('delta_i lower than zero')
        pause
    end
    if delta_o < 0
        disp('delta_o lower than zero')
        pause
    end
    F_i=K_i*delta_i^(3/2); 
    F_o=K_o*delta_o^(3/2);
    %disp([delta_i, delta_o]);
    %sign(delta_i)*
    % Equilibrium of force in horizontal and vertical directions----
    
    % Race control check ++++
    % lambda=1 for outer race control
    % lambda=0 for inner race control
    [lambda]=f_race_control(F_i,F_o,C_beta_i,C_beta_o,B);
    
    % Centrifugal force ++++
    d_pe=B.d_e+2*L_4-2*B.d*(B.f_o-0.5)*cos(B.beta_f); % Eq. (3.107)
    
    % ball spin axis angle
    if lambda==0
        zeta=atan( (d_pe*S_beta_i)/(d_pe*C_beta_i-B.d) ) ; %Eq. (3.89), inner race control
    elseif lambda==1 
        zeta=atan( (d_pe*S_beta_o)/(d_pe*C_beta_o+B.d) ); %Eq. (3.90), outer race control
    lambda=0.5;
    end
    
    
    % cos and sin function cant be used in iteration, transformation
    % equation cos(x-y)=cos(x)*cos(y)+sin(x)*sin(y)
    %cos(beta_i-zeta)=
    C_beta_i_zeta=C_beta_i*cos(zeta)+S_beta_i*sin(zeta);
    % cos(beta_o-zeta)
    C_beta_o_zeta=C_beta_o*cos(zeta)+S_beta_o*sin(zeta);
    
    m=B.rho_b*4/3*pi*(B.d/2)^3; %mass of the ball
    omega_B=-B.Omega_i/(B.d*( C_beta_i_zeta/(d_pe-B.d*C_beta_i) + (C_beta_o_zeta/(d_pe+B.d*C_beta_o)) )); %spin speed of the ball around its own axis (3.81)
    omega_c=B.Omega_i/(1+ ( (d_pe+B.d*C_beta_o)/(d_pe-B.d*C_beta_i) *(C_beta_i_zeta/C_beta_o_zeta)  )); %cage speed (3.82)???
    
    F_c=1/2*m*d_pe*omega_c^2;
    
    % Gyroscopic moment ++++
    I_p=2/5*m*(B.d/2)^2; %Mass moment of the inertia of the ball
    M_g=I_p*omega_B*omega_c*sin(zeta);
    %M_g=0;
       
    % Equations that will be solved iteratively
    Eq(4*(j-1)+3,1)=F_o*S_beta_o-F_i*S_beta_i-(2*M_g/B.d)*(lambda*C_beta_o-(1-lambda)*C_beta_i);
    Eq(4*(j-1)+4,1)=F_o*C_beta_o-F_i*C_beta_i+(2*M_g/B.d)*(lambda*S_beta_o-(1-lambda)*S_beta_i)-F_c;
%
    %Alternative formulation
    if j==1; EQ_tmp1=0; EQ_tmp2=0; end;
   % EQ_tmp1=EQ_tmp1 + (F_i*S_beta_i+2*(1-lambda)*M_g*C_beta_i/B.d); %%%%%%%CHECK + sign in front of 2*M_g!!!!!!!!!!!!!!!!!!
   % EQ_tmp2=EQ_tmp2 + (F_i*C_beta_i-2*(1-lambda)*M_g*S_beta_i/B.d)*cos(psi); % Eq. (3.112)
    %(F_i*C_beta_i-2*M_g/B.d*(1-lambda)*S_beta_i)*cos(psi)
   
    
    %v1.0
%     EQ_tmp1=EQ_tmp1 + (F_i*S_beta_i-2*M_g/B.d*(1-lambda)*C_beta_i);  %%%%%%%CHECK + sign in front of 2*M_g!!!!!!!!!!!!!!!!!! IN THE BOOK is - sign!
%     EQ_tmp2=EQ_tmp2 + (F_i*C_beta_i-2*M_g/B.d*(1-lambda)*S_beta_i)*cos(psi); % Eq. (3.112)

    %v2.0
   EQ_tmp1=EQ_tmp1 + (F_i*S_beta_i-(M_g/B.d)*C_beta_i);  %%%%%%%CHECK + sign in front of 2*M_g!!!!!!!!!!!!!!!!!! IN THE BOOK is - sign!
   EQ_tmp2=EQ_tmp2 + (F_i*C_beta_i+(M_g/B.d)*S_beta_i)*cos(psi); % Eq. (3.112)

    
    Fi(j)=F_i;
    Fo(j)=F_o;
    Fc(j)=F_c;
    
    L3(j)=L_3;
    L4(j)=L_4;
    deltai(j)=delta_i;
    deltao(j)=delta_o;
    
    Ki(j)=K_i;
    Ko(j)=K_o;
    EQtmp11(j)=(F_i*S_beta_i-(M_g/B.d)*C_beta_i);
    EQtmp22(j)=(F_i*C_beta_i+(M_g/B.d)*S_beta_i)*cos(psi);
    lambdaa(j)=lambda;
    Mg(j)=M_g;
    


    
end


Eq(4*(j-1)+5,1)=FR(1)-EQ_tmp1;
Eq(4*(j-1)+6,1)=FR(2)-EQ_tmp2;
% tulos_axial(i)=Eq(57);
% tulos_radial(i)=Eq(58);

% % disp('eka')
% % FR(1)
% % EQ_tmp1
% % disp('ois kiva jos ois -300N')
% % axiaaliresiduaali=Eq(57)
% % 
% % disp('toka')
% % FR(2)
% % EQ_tmp2
% % disp('ois kiva jos ois 700N')
% % radiaaliresiduaali=Eq(58)
% % 
%     figure(4)
% plot(EQtmp11)
% hold on
% plot(EQtmp22,'r')
% pause
