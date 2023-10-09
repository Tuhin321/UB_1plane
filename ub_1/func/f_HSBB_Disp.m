function [Displ]=f_HSBB_Disp(FR,B,omega)

%WITH HIGH-SPEED EFFECTS
% ex=Disp(1); %radial
% ey=Disp(2); %radial
% ez=Disp(3); %axial
% Rotation speed in revolutions per minute
if B.showRPM == 1
omega_rpm=(omega*60)/(2*pi);
disp(['Currently running ', num2str(omega_rpm),' rpm'])
end
% Initial guess for the axial and radial displacement
% delta_t=B.delta_ti(B.tk); 
% delta_r=B.delta_r(B.rk);
% L_3=sin(B.beta_f)*B.d*(B.f_o-0.5);  % Fig. 3.11, p. 94
% L_4=cos(B.beta_f)*B.d*(B.f_o-0.5);  
% delta_i=7.0e-7;
% delta_o=48.0e-7;

% Initial guesses for the iteration
delta_t=B.delta_t; 
delta_r=B.delta_r;
L_3=B.L_3;  % Fig. 3.11, p. 94
L_4=B.L_4;  
delta_i=B.delta_i;
delta_o=B.delta_o;


% Vars=[L_3; L_4; delta_o; delta_i; delta_t; delta_r]
X0=[L_3; L_4; delta_i; delta_o];
% pause
X1=X0;
for ii=1:B.n-1
    X1=[X1;
        X0];
end
delta0 =[delta_t; delta_r];
X0=[X1
    delta0];

% Radial force X,Y,Z (x and y are radial forces, z is axial)
F_r=sqrt(FR(1)^2+FR(2)^2); % Radial force varied in x,y and z. => computed the axial and radial forces
F_a=FR(3);

theta2=atan2(FR(1),FR(2)); % 90 degrees 

FR_hs=[F_a; F_r]; % axial and radial

B.Omega_i=omega;

%%%%%%%X0(4*B.n+2) = radial
%%%%%%%X0(4*B.n+1) = axial
[X]=f_NR_Iter3v2('f_BB3_Equi',X0,B,FR_hs); % local

disp(['Axial Displacement  : ' num2str(X(4*B.n+1)*1e6,3) ' um'])
disp(['Radial Displacement (x) : ' num2str(X(4*B.n+2)*1e6*sin(theta2),3) ' um'])
disp(['Radial Displacement (y) : ' num2str(X(4*B.n+2)*1e6*cos(theta2),3) ' um'])

Displ=[X(4*B.n+2)*sin(theta2); X(4*B.n+2)*cos(theta2); X(4*B.n+1)];
% attitude angles of the balls stored in a vector 
% x=rad,y=rad,z=ax
% beta0=0:2*pi/Bearing.z:1.99*pi;
% beta0=beta0';

% calculating forces using final displacements
% Iterated displacements for every point in local coordinate system
[Eq, beta_i,beta_o,Finner,Fo,Fc,L3,L4,deltai,deltao,Ki,Ko,EQtmp1,EQtmp2,lambda,Mg]=f_BB3_Equi(X,B,FR_hs);



outercangle=beta_o*180/pi;
% innercangle=beta_i*180/pi;

c=1;
for kk=1:B.n
    c=c+1;
    if kk == B.n
        c=c-1;
    end
    deltacangle(kk)=outercangle(kk)+outercangle(c)/2;
end
% deltacangle
ratio=(sum(deltacangle)/length(deltacangle))/deltacangle(1);
% sum(deltacangle)
% length(deltacangle)
if  1.05 < ratio > 0.95 % if the ratio is higher 1.05 or 0.95
    disp('Contact angle is variating too much. Results are not correct! Make better initial guess or change the load step amount!')
    pause
end

if B.plotting == 1;
     for jj=1:B.n
    for ii=1:4  
    X_plot(ii,jj)=X(4*(jj-1)+ii);
    end
end

figure(3)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)

subplot(2,1,1)
plot(1:B.n,X_plot(1,:)*1e3,'b','linewidth',2)
hold on
plot(1:B.n,X_plot(2,:)*1e3,'r','linewidth',2)
legend('L_{3}','L_{4}',0)
ylabel('[mm]')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;

subplot(2,1,2)
plot(1:B.n,X_plot(3,:)*1e6,'b','linewidth',2)
hold on
plot(1:B.n,X_plot(4,:)*1e6,'r','linewidth',2)
legend('\delta_{i}','\delta_{o}',0 )
ylabel('Elastic Deformation [\mum]')
xlabel('Ball #')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)
% print -dtiff def_steel_1.tif
print -dtiff def_HC_1.tif
% pause

% pause
figure(4)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)
plot(1:B.n,beta_i*180/pi,'b','linewidth',2)
hold on
plot(1:B.n,beta_o*180/pi,'r','linewidth',2)
legend('\beta_{i}','\beta_{o}',0 )
ylabel('Contact Angle [\circ]')
xlabel('Ball #')
v=axis;
axis([1 B.n 0 30])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)
% print -dtiff angle_steel_1.tif
print -dtiff angle_HC_1.tif

figure(5)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)
plot(1:B.n,Finner,'b','linewidth',2)
hold on
plot(1:B.n,Fo,'r','linewidth',2)
plot(1:B.n,Fc,'g','linewidth',2)
legend('F_{i}','F_{o}','F_{c}',0 )
ylabel('Force [N]')
xlabel('Ball #')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)
% print -dtiff force_steel_1.tif
%print -dtiff force_HC_1.tif



figure(6)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)
plot(1:B.n,L4,'b','linewidth',2)
hold on
plot(1:B.n,L3,'r','linewidth',2)
plot(1:B.n,deltai,'g','linewidth',2)
plot(1:B.n,deltao,'g','linewidth',2)
legend('L4','L3','delta_i','delta_o',0 )
ylabel('Force [N]')
xlabel('Ball #')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)

figure(7)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)
plot(1:B.n,Ki,'b','linewidth',2)
hold on
plot(1:B.n,Ko,'r','linewidth',2)
legend('K_i','K_o',0 )
ylabel('Stiffness [N/m]')
xlabel('Ball #')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)

figure(8)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)
plot(1:B.n,EQtmp1,'b','linewidth',2)
hold on
plot(1:B.n,EQtmp2,'r','linewidth',2)
legend('EQtmp1','EQtmp2',0 )
ylabel('EQ_tmp [N]')
xlabel('Ball #')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)

figure(9)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)
plot(1:B.n,lambda,'b','linewidth',2)
hold on
%plot(1:B.n,EQtmp2,'r','linewidth',2)
legend('Lambda',0 )
ylabel('Value')
xlabel('Ball #')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)

figure(10)
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.2 0.2 0.7 0.7])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)
plot(1:B.n,Mg,'b','linewidth',2)
hold on
%plot(1:B.n,EQtmp2,'r','linewidth',2)
legend('Lambda',0 )
ylabel('Value')
xlabel('Ball #')
v=axis;
axis([1 B.n v(3) v(4)])
grid on;
str=[B.type ', F_{r}=' num2str(B.F_r) ' N, F_{a}='  num2str(B.F_a) ' N, ' num2str(B.Omega_i*60/(2*pi)) ' rpm'];
title(str, 'Fontsize', 13)



return
 else
 end

% %Z IS THE AXIAL DIRECTION, X,Y RADIAL
% FX=-sum(Finner.*cos(beta_i).*cos(beta')); %radial
% FY=-sum(Finner.*cos(beta_i).*sin(beta'));  %radial force
% FZ=-sum(Finner.*sin(beta_i));             %axial force
% % Momentit
% TX=-sum(Finner.*(R_inner+Bearing.d/2).*sin(beta_i).*sin(beta'));
% TY=-sum(Finner.*(R_inner+Bearing.d/2).*sin(beta_i).*(-cos(beta')));


