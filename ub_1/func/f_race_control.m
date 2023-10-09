function [lambda]=f_race_control(F_i,F_o,C_beta_i,C_beta_o,B)
% Function calculates the friction force both in the inner ring and outer
% ring, compares them and returns the race control parameter lambda
% lambda = 1 for outer-race control
% lambda = 0 for inner-race control

% Contact parameters
%The effective radii of elliptical contact conjunctions p. 55
%inner race
R_x_i=B.d*(B.d_e-B.d*C_beta_i)/(2*B.d_e); %Eq. (2.28)
R_y_i=B.f_i*B.d/(2*B.f_i-1);              %Eq. (2.29)
R_i=(1/R_x_i+1/R_y_i)^(-1);         %Eq. (2.24)
%outer race
R_x_o=B.d*(B.d_e+B.d*C_beta_o)/(2*B.d_e);
R_y_o=B.f_o*B.d/(2*B.f_o-1);
R_o=(1/R_x_o+1/R_y_o)^(-1);         %Eq. (2.24)

%ellipticity parameter p. 75
k_e_i=1.0339*(R_y_i/R_x_i)^0.6360;  %Eq. (3.28)
k_e_o=1.0339*(R_y_o/R_x_o)^0.6360;

% Elliptical integrals--
E2_i=1.0003+0.5968/(R_y_i/R_x_i);  % Eq. (3.29), elliptical integral of the second kind
% E1_i=1.5277+0.6023*log(R_y_i/R_x_i); % Eq. (3.30), elliptical integral of the first kind
E2_o=1.0003+0.5968/(R_y_o/R_x_o);
% E1_o=1.5277+0.6023*log(R_y_o/R_x_o);

% Effective modulus of elasticity p. 73
E_p=2 / ( (1-B.nu_a^2)/B.E_a + (1-B.nu_b^2)/B.E_b ); %Eq. (3.16)

% semimajor axis of the contact
a_i=(6*k_e_i^2*E2_i*F_i*R_i/(pi*E_p))^(1/3);
a_o=(6*k_e_o^2*E2_o*F_o*R_o/(pi*E_p))^(1/3);

mu=0.07; % coefficient of sliding friction

% Friction force of the contact ellipse
M_s_i=3/8*mu*F_i*a_i*E2_i;
M_s_o=3/8*mu*F_i*a_o*E2_o;

% Comparison
if M_s_o >= M_s_i 
    lambda=1; % outer race control
else
    lambda=0; % inner-race control
end