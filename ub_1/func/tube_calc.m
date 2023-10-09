function res = tube_calc(Do, Di, nu_or_kx)
    A = pi*(Do^2-Di^2)/4;
    Izz = pi*(Do^4-Di^4)/64;
    Iyy = pi*(Do^4-Di^4)/64;
    H = Do;
    B = Di;
    theta = 0;
    Iv = 0;
    Ixx = pi*(Do^4-Di^4)/32;
    if Di==0 %barrel
        kx = nu_or_kx;
        res = [A, Izz,Iyy,H,B,theta,Iv,Ixx,kx,kx];
    else
        nu = nu_or_kx;
        Rr=Di/Do;
        shearz = (6*(1+nu)*(1+Rr^2)^2)/((7+6*nu)*(1+Rr^2)^2+(20+12*nu)*Rr^2);
        sheary = (6*(1+nu)*(1+Rr^2)^2)/((7+6*nu)*(1+Rr^2)^2+(20+12*nu)*Rr^2);
        res = [A, Izz,Iyy,H,B,theta,Iv,Ixx,shearz,sheary];
    end
end