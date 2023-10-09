function [K_i,K_o]=f_Hertz_K(C_beta_i,C_beta_o,B)
% Function calculates the hertz contact stiffness coefficients for the
% inner and outer race contacts


% % the radii of curvature of the ball inner race contact p. 53
% r_ax_i=d/2;
% r_ay_i=d/2;
% r_bx_i=(d_e-d*cos(beta_i))/(2*cos(beta_i));
% r_by_i=-f_i*d;
% 
% % the radii of curvature of the ball outer race contact p. 54
% r_ax_o=d/2;
% r_ay_o=d/2;
% r_bx_o=-(d_e+d*cos(beta_o))/(2*cos(beta_o));
% r_by_o=-f_o*d;

% The effective radii of elliptical contact conjunctions p. 55
% inner race
R_x_i=B.d*(B.d_e-B.d*C_beta_i)/(2*B.d_e); % Eq. (2.28)
R_y_i=B.f_i*B.d/(2*B.f_i-1);              % Eq. (2.29)
R_i=(1/R_x_i+1/R_y_i)^(-1);         % Eq. (2.24)
% outer race
R_x_o=B.d*(B.d_e+B.d*C_beta_o)/(2*B.d_e);
R_y_o=B.f_o*B.d/(2*B.f_o-1);
R_o=(1/R_x_o+1/R_y_o)^(-1);         % Eq. (2.24)

% ellicticity parameter p. 75
k_e_i=1.0339*(R_y_i/R_x_i)^0.6360;  % Eq. (3.28)
k_e_o=1.0339*(R_y_o/R_x_o)^0.6360;

% Elliptical integrals--
E2_i=1.0003+0.5968/(R_y_i/R_x_i);  % Eq. (3.29), elliptical integral of the second kind
E1_i=1.5277+0.6023*log(R_y_i/R_x_i); % Eq. (3.30), elliptical integral of the first kind
E2_o=1.0003+0.5968/(R_y_o/R_x_o);
E1_o=1.5277+0.6023*log(R_y_o/R_x_o);

% Effective modulus of elasticity p. 73
E_p=2 / ( (1-B.nu_a^2)/B.E_a + (1-B.nu_b^2)/B.E_b ); % Eq. (3.16)


%contact stiffness coefficient
K_i=pi*k_e_i*E_p*sqrt(R_i*E2_i /(4.5*E1_i^3));
K_o=pi*k_e_o*E_p*sqrt(R_o*E2_o /(4.5*E1_o^3));






