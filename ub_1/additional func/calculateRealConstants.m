
function res = calculateRealConstants(Do, Di, nu)

A = pi*(Do^2-Di^2)/4;
Izz = pi*(Do^4-Di^4)/64;
Iyy = pi*(Do^4-Di^4)/64;
H = Do;
B = Di;
theta = 0;
Iv = 0;
Ixx = pi*(Do^4-Di^4)/32;
Rr=Di/Do;
shearz = (6*(1+nu)*(1+Rr^2)^2)/((7+6*nu)*(1+Rr^2)^2+(20+12*nu)*Rr^2);
sheary = (6*(1+nu)*(1+Rr^2)^2)/((7+6*nu)*(1+Rr^2)^2+(20+12*nu)*Rr^2);
res = [A, Izz,Iyy,H,B,theta,Iv,Ixx,shearz,sheary];

end
