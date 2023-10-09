function [Forcevec]=f_DGBB_Forces(Disp, Vel, Bearing)
% Ball Bearing force calculation routine
% Written by Jussi Sopanen

% INPUT: 
% Disp= [ex, ey, ez, gammax, gammay, thetaz ]
% Vel= [vx, vy, vz, omegax, omegay, omegaz]
% Dim= [d1, D, d, z, conform_rat, cd]
% Bearing = struc array, joka sisältää vakioparametrit

ex=Disp(1);
ey=Disp(2);
ez=Disp(3);
gammax=Disp(4);
gammay=Disp(5);
thetaz=Disp(6);

vx=Vel(1);
vy=Vel(2);
vz=Vel(3);
omegax=Vel(4);
omegay=Vel(5);
omegaz=Vel(6);

% Lasketaan nominaali laakerin voima ulkoisesta kuormasta (tarvitaan
% öljykalvon pakususden laskennassa)
F_nominal=sqrt(Bearing.LoadVec(1)^2+Bearing.LoadVec(2)^2+Bearing.LoadVec(3)^2);

% Effective modulus of elasticity
E_effect=(0.5* ((1-Bearing.nu_a^2)/Bearing.E_a+(1-Bearing.nu_b^2)/Bearing.E_b))^-1;
fi=5/Bearing.z*F_nominal; % maximum ball load

%------------------------------------------------------------------
% radius calculation 
%r_bx_out=-(Bearing.dm/2+Bearing.d/2+Bearing.cd/4); %radii of outer race <0
% otetaan kontaktikulma huomioon
r_bx_out=-(Bearing.dm+(Bearing.d+(Bearing.cd/2))*cos(Bearing.alfa0) ) / (2*cos(Bearing.alfa0));

%r_bx_in=(Bearing.dm/2-Bearing.d/2-Bearing.cd/4); %radii of inner race >0
% otetaan kontaktikulma huomioon
r_bx_in=(Bearing.dm-(Bearing.d+(Bearing.cd/2))*cos(Bearing.alfa0) ) / (2*cos(Bearing.alfa0));

R_outer=abs(r_bx_out); % Näitä käytetään beta kulman laskennassa
R_inner=abs(r_bx_in); 

r_by_out=-(Bearing.conform_rat_out*Bearing.d); %outer groove radius <0
r_by_in=-(Bearing.conform_rat_in*Bearing.d);  %inner groove radius <0

%The effective radii of elliptical contact conjunctions
%outer race
R_x_out=(1/(Bearing.d/2)+1/r_bx_out)^-1;
R_y_out=(1/(Bearing.d/2)+1/r_by_out)^-1;
R_out=(1/R_x_out+1/R_y_out)^-1;
%inner race
R_x_in=(1/(Bearing.d/2)+1/r_bx_in)^-1;
R_y_in=(1/(Bearing.d/2)+1/r_by_in)^-1;
R_in=(1/R_x_in+1/R_y_in)^-1;

%ellicticity parameter
k_e_out=1.0339*(R_y_out/R_x_out)^0.6360;
k_e_in=1.0339*(R_y_in/R_x_in)^0.6360;

%----------------------------------------------------------------------------
%contact stiffness coefficient

Xi_in=1.0003+0.5968*(R_x_in/R_y_in);
zeta_in=1.5277+0.6023*log(R_y_in/R_x_in);

Xi_out=1.0003+0.5968*(R_x_out/R_y_out);
zeta_out=1.5277+0.6023*log(R_y_out/R_x_out);

Kc_in=pi*k_e_in*E_effect*sqrt(R_in*Xi_in /(4.5*zeta_in^3));
Kc_out=pi*k_e_out*E_effect*sqrt(R_out*Xi_out /(4.5*zeta_out^3));
% Total stiffness coefficient
Kc_tot=( ( ((1/Kc_in)^(2/3)) + ((1/Kc_out)^(2/3)) )^(3/2))^-1;


% FILM THICKNESS
%Dimensioless material parameter
G_viiva_out=Bearing.alfa*E_effect;
G_viiva_in=Bearing.alfa*E_effect;
%Dimensioless load parameter
W_viiva_out=fi/(E_effect*R_x_out^2);
W_viiva_in=fi/(E_effect*R_x_in^2);


%Sum of surface velocities 
U_in=(Bearing.dm^2-Bearing.d^2)/(4*Bearing.dm)*omegaz;
U_out=U_in;
U=U_in+U_out;

%Dimensioless speed parameter
U_viiva_out=(Bearing.eta0*U)/(2*E_effect*R_x_out);
U_viiva_in=(Bearing.eta0*U)/(2*E_effect*R_x_in);

%Nominal film thickness for an elliptical contact problem
h0_out=2.69*R_x_out*((U_viiva_out)^0.67)*(G_viiva_out^0.53)*(W_viiva_out^-0.067)*(1-0.61*exp(-0.73*k_e_out));
h0_in=2.69*R_x_in*((U_viiva_in)^0.67)*(G_viiva_in^0.53)*(W_viiva_in^-0.067)*(1-0.61*exp(-0.73*k_e_in));
%Total film thickness
h0=h0_out+h0_in;

% attitude angles of the balls stored in a vector 
beta0=0:2*pi/Bearing.z:1.99*pi;
beta0=beta0';
% beta0deg=beta0*180/pi;
theta_out=0;
theta_in=thetaz;
beta=beta0+(theta_in*R_inner+theta_out*R_outer) / (2*(Bearing.d/2 + R_in));

% radial eccentrity
er=ex*cos(beta) + ey*sin(beta);
% axial eccentrity
ea=ez - (-gammax*sin(beta) + gammay*cos(beta)) * (R_in+abs(r_by_in));
%contact angle
fii=atan(ea./(er+R_inner+abs(r_by_in)-R_outer+abs(r_by_out)));

% elastic deformation
di=abs(r_by_out)+abs(r_by_in)-( er+R_inner+abs(r_by_in)-R_outer+abs(r_by_out) ).*cos(fii).^(-1);
delta = Bearing.d+h0-di;
% contact force
Fi=zeros(Bearing.z,1);
for i=1:Bearing.z
    if delta(i)>0
        Fi(i)=Kc_tot*delta(i)^1.5;
    else
        Fi(i)=0;
    end
end

%Lasketaan akseliin vaikuttavat voimat
FX=-sum(Fi.*cos(fii).*cos(beta));
FY=-sum(Fi.*cos(fii).*sin(beta));
FZ=-sum(Fi.*sin(fii));
% Momentit
TX=-sum(Fi*(R_inner+Bearing.d/2).*sin(fii).*sin(beta));
TY=-sum(Fi*(R_inner+Bearing.d/2).*sin(fii).*(-cos(beta)));

% Kootaan voimavektori
Forcevec=[FX, FY, FZ, TX, TY]';
 