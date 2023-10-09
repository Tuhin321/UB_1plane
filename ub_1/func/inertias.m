function [mpx,mp] = inertias(x1,x2,d1,d2,rho)
% mass points for cylinders
% mass,   Jp,       Jt,       L,         dia]
  dia=d2;
  L=(x2-x1);
  mpx=L/2+x1;
  m=pi*((0.5*d2)^2-(0.5*d1)^2)*L*rho; 
  Ip=0.5*m*((0.5*d1)^2+(0.5*d2)^2);
  It=1/12*m*(3*((0.5*d1)^2+(0.5*d2)^2)+L^2);
  mp=[m, Ip, It, L, dia]; 
end