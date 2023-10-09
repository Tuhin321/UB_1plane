function [Forcevec]=f_SRB_Forces(Disp, Vel,Inp, Bearing,Request,theta,phi)
% Ball Bearing force calculation routine
% Written by Jussi Sopanen

% Testing does default transformation matrix work on all files
Bearing.A=[0  0 -1
    0  1  0
    1  0  0];

Disp(1:3)=Bearing.A*Disp(1:3)';
Disp(4:6)=Bearing.A*Disp(4:6)';

% Effective modulus of elasticity
E_effect=(0.5* ((1-Bearing.nu_a^2)/Bearing.E_a+(1-Bearing.nu_b^2)/Bearing.E_b))^-1;
% fi=5/Bearing.z*F_nominal; % maximum ball load

%------------------------------------------------------------------
% radius calculation 
% taking contact angle into account
r_bx_out=-(Bearing.dm+(Bearing.d+(Bearing.cd/2))*cos(Bearing.alfa0) ) / (2*cos(Bearing.alfa0));

%%%%%%%%%%new to convert BB to SRB%%%%%%%%%%%%%%%
r_ax_in=Bearing.d/2; 
r_ay_in=Bearing.rr; 
r_ax_out=Bearing.d/2; 
r_ay_out=Bearing.rr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% r_bx_in=(Bearing.dm/2-Bearing.d/2-Bearing.cd/4); %radii of inner race >0
% taking contact angle into account
r_bx_in=(Bearing.dm-(Bearing.d+(Bearing.cd/2))*cos(Bearing.alfa0) ) / (2*cos(Bearing.alfa0));

R_outer=abs(r_bx_out); % these are used for calculating beta angle
R_inner=abs(r_bx_in);

%%%%%%%%%%%%%%%%%new to convert BB to SRB%%%%%%%%%%%%%%%

r_by_out=-Bearing.ro;
r_by_in=-Bearing.ri;
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%new to convert BB to SRB%%%%%%%%%%%%%%%
%The effective radii of elliptical contact conjunctions
%outer race
% R_x_out=(1/(Bearing.d/2)+1/r_bx_out)^-1;
% R_y_out=(1/(Bearing.d/2)+1/r_by_out)^-1;
R_x_out=(1/r_ax_out+1/r_bx_out)^-1;
R_y_out=(1/r_ay_out+1/r_by_out)^-1;
R_out=(1/R_x_out+1/R_y_out)^-1;
%inner race
% R_x_in=(1/(Bearing.d/2)+1/r_bx_in)^-1;
% R_y_in=(1/(Bearing.d/2)+1/r_by_in)^-1;
R_x_in=(1/1/r_ax_in+1/r_bx_in)^-1;
R_y_in=(1/1/r_ay_in+1/r_by_in)^-1;
R_in=(1/R_x_in+1/R_y_in)^-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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





    Forcevec=zeros(6,1);

% attitude angle increment
beta_inc=(2*pi)/Bearing.z; %(360/Bearing.z)*(pi/180);

if strcmp(Request.Defect.Type,'Outer_race_defect')
    
    
   for i=1:Bearing.nz

        if i==1
            fii=-Bearing.alfa0;
        else
            fii=Bearing.alfa0;
        end
        
        
          %loop across rollers
        for j=1:Bearing.z  
            
            %attitude angle of roller 'j' in row 'i'
            beta_ji=(j-1)*beta_inc+theta - 2*pi*floor(((j-1)*beta_inc+theta)/(2*pi));  %+0*(i-1)*10*pi/180; 
            
             if  Bearing.position-Inp.theta1<beta_ji && beta_ji<Bearing.position+Inp.theta1   %defect position
                SS=eval(Inp.y);
                delta_plus=SS;                
                if strcmp(Request.Defect.Form,'Simple')
                    delta_plus=Bearing.depth;
                end  

                
                
                A0=Bearing.ro+Bearing.ri-Bearing.d-Bearing.cd/2;
                er=Disp(1)*cos(beta_ji)+Disp(2)*sin(beta_ji);
                ez=0;
                A_beta=A0+er*cos(fii)+ez*sin(fii);
                dist=Bearing.ro+Bearing.ri-A_beta;
                delta_beta=Bearing.d-dist-delta_plus;
                
                if delta_beta >= 0
                    Q_ji=-Kc_tot*delta_beta^3./2; %SIGN???
                    %                               Q_ji=0; %SIGN???
                else
                    Q_ji=0;
                end
                Forcevec(1)=Forcevec(1)+Q_ji*cos(fii)*cos(beta_ji);
                Forcevec(2)=Forcevec(2)+Q_ji*cos(fii)*sin(beta_ji);
                Forcevec(3)=Forcevec(3)+Q_ji*sin(fii);
                
             else
                           
                 A0=Bearing.ro+Bearing.ri-Bearing.d-Bearing.cd/2;
                 er=Disp(1)*cos(beta_ji)+Disp(2)*sin(beta_ji);
                 ez=0;
                 A_beta=A0+er*cos(fii)+ez*sin(fii);
                 dist=Bearing.ro+Bearing.ri-A_beta;
                 delta_beta=Bearing.d-dist;
            
                 if delta_beta >= 0
                    Q_ji=-Kc_tot*delta_beta^3./2; %SIGN???
                 else
                    Q_ji=0;
                 end
                 Forcevec(1)=Forcevec(1)+Q_ji*cos(fii)*cos(beta_ji);
                 Forcevec(2)=Forcevec(2)+Q_ji*cos(fii)*sin(beta_ji);
                 Forcevec(3)=Forcevec(3)+Q_ji*sin(fii);
            
             end   
        end
   end
   
elseif strcmp(Request.Defect.Type,'Inner_race_defect')
    
    for i=1:Bearing.nz

        if i==1
            fii=-Bearing.alfa0;
        else
            fii=Bearing.alfa0;
        end
        
        
          %loop across rollers
        for j=1:Bearing.z  
            
            %attitude angle of roller 'j' in row 'i'
            beta_ji=(j-1)*beta_inc+theta - 2*pi*floor(((j-1)*beta_inc+theta)/(2*pi));
            
             
             if  phi-Inp.theta1<beta_ji && beta_ji<phi+Inp.theta1   
                
                SS=eval(Inp.y);
                delta_plus=SS; 
                
                if strcmp(Request.Defect.Form,'Simple')
                    delta_plus=Bearing.depth;
                end  

                
                
                A0=Bearing.ro+Bearing.ri-Bearing.d-Bearing.cd/2;
                er=Disp(1)*cos(beta_ji)+Disp(2)*sin(beta_ji);
                ez=0;
                A_beta=A0+er*cos(fii)+ez*sin(fii);
                dist=Bearing.ro+Bearing.ri-A_beta;
                delta_beta=Bearing.d-dist-delta_plus;
                
                if delta_beta >= 0
                    Q_ji=-Kc_tot*delta_beta^3./2; %SIGN???
                    %                               Q_ji=0; %SIGN???
                else
                    Q_ji=0;
                end
                Forcevec(1)=Forcevec(1)+Q_ji*cos(fii)*cos(beta_ji);
                Forcevec(2)=Forcevec(2)+Q_ji*cos(fii)*sin(beta_ji);
                Forcevec(3)=Forcevec(3)+Q_ji*sin(fii);
                
             else
                           
                 A0=Bearing.ro+Bearing.ri-Bearing.d-Bearing.cd/2;
                 er=Disp(1)*cos(beta_ji)+Disp(2)*sin(beta_ji);
                 ez=0;
                 A_beta=A0+er*cos(fii)+ez*sin(fii);
                 dist=Bearing.ro+Bearing.ri-A_beta;
                 delta_beta=Bearing.d-dist;
            
                 if delta_beta >= 0
                    Q_ji=-Kc_tot*delta_beta^3./2; %SIGN???
                 else
                    Q_ji=0;
                 end
                 Forcevec(1)=Forcevec(1)+Q_ji*cos(fii)*cos(beta_ji);
                 Forcevec(2)=Forcevec(2)+Q_ji*cos(fii)*sin(beta_ji);
                 Forcevec(3)=Forcevec(3)+Q_ji*sin(fii);
            
             end   
        end
   end   
    
else 
    
                 A0=Bearing.ro+Bearing.ri-Bearing.d-Bearing.cd/2;
                 er=Disp(1)*cos(beta_ji)+Disp(2)*sin(beta_ji);
                 ez=0;
                 A_beta=A0+er*cos(fii)+ez*sin(fii);
                 dist=Bearing.ro+Bearing.ri-A_beta;
                 delta_beta=Bearing.d-dist;
            
                 if delta_beta >= 0
                    Q_ji=-Kc_tot*delta_beta^3./2; %SIGN???
                 else
                    Q_ji=0;
                 end
                 Forcevec(1)=Forcevec(1)+Q_ji*cos(fii)*cos(beta_ji);
                 Forcevec(2)=Forcevec(2)+Q_ji*cos(fii)*sin(beta_ji);
                 Forcevec(3)=Forcevec(3)+Q_ji*sin(fii);   
    
end
Forcevec(1:3)=Bearing.A'*Forcevec(1:3);   
% % attitude angles of the balls stored in a vector 
%           beta0=(0:2*pi/Bearing.z:1.99*pi)';
% 
% %beta0deg=beta0*180/pi;
%           theta_out=0;
%           theta_in=Disp(6);
% %           beta=beta0+(theta_in*R_inner+theta_out*R_outer) / (2*(Bearing.d/2 + R_in));
%           beta=beta0+theta_in/2*(1-Bearing.d/Bearing.dm*cos(fii));
% 
%       
%             A0=Bearing.ro+Bearing.ri-Bearing.d-Bearing.cd/2;
% %             A0=ro+ri-dr-cd;
%             %axial deflection for roller 'j' in row 'i'
%             delta_z=A0*sin(fii)+Disp(3);
%             %radial displacement for roller 'j' in row 'i'
%             delta_r=A0*cos(fii)+Disp(1)*cos(beta)+Disp(2)*sin(beta);
%             %loaded distance  between curvature centers
%             A_beta=sqrt(delta_z.^2+delta_r.^2);
%             
%             %Distance between race surfaces along the common normal
%             dist=Bearing.ro+Bearing.ri-A_beta;
%             %Deflection of roller
%             delta=Bearing.d-dist;
%             fiii=atan(delta_z/delta_r)';
%             pd=delta>0; %vektori, jonka elementti on 1 jos delta >0
%             Fi=Kc_tot*pd.*delta.^1.5;
%             Force=zeros(6,1);
%             Force(1)=-sum(Fi.*cos(beta).*cos(fiii));%-sum(Fi.*cos(fii).*cos(beta));
%             Force(2)=-sum(Fi.*sin(beta).*cos(fiii));%-sum(Fi.*cos(fii).*sin(beta));
% 
%             
%             Forcevec(1:3)=Bearing.A'*Force(1:3);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % 


   

end





 