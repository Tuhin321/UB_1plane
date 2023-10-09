function f_MaximumCapacity_v3_4AMB(Fsat,wrangef,filename,gap,sys)
% FOR ABCD BEARING CONFIG


%Scale the input and output units to SI-units
% inch2meter=2.54e-2;     %inches [in] to meters [m]
% lbf2N=4.44822162;       %pounds [lbf] to newtons [N]

A_rotor=sys.A_1dir; 
B_rotor=sys.B_1dir; 
C_rotor=sys.C_1dir; % /lbf2N required if the model is in US-units
D_rotor=sys.D_1dir;
%Load the rotor state space matrices
% A_rotor=load(['A_' filename],'-ascii'); 
% B_rotor=load(['B_' filename],'-ascii'); 
% B_rotor=B_rotor/lbf2N; % /lbf2N required if the model is in US-units
% C_rotor=load(['C_' filename],'-ascii'); 
% C_rotor=C_rotor*inch2meter; %*inch2meter required if the model is in US-units
% G_rotor=load(['G_' filename],'-ascii'); %Gyroscopic matrix not used
Fsat1=Fsat(1);
Fsat2=Fsat(2);
Fsat3=Fsat(3);
Fsat4=Fsat(4);
gap1=gap(1);
gap2=gap(2);
gap3=gap(3);
gap4=gap(4);
%Organize inputs and outputs
% impeller_in=[1]; %Impeller in (m) (Impeller and Axial disks at end of shafts)
% journals_in=[3 4 11]; %AMB in (N)
% 
% sensors_out=[3 5 14]; %Sensor out (V)
% journals_out=[4 6 13]; %AMB out (N)
% impeller_out=[1 16]; %Impeller out (m) (Impeller and Axial disks)

%%Recompile A-matrix
modes=size(A_rotor,1)/2;
downleft=A_rotor(modes+1:2*modes,1:modes);
%wanted damping ratios (0.01 = 1% damping)
damps=[0.000 0.000 0.01 0.01 0.01 0.01];  % for flexible mode
omegas=diag(downleft);
%Calculate damping values -2*ksi*w
dampval=diag(sqrt(-omegas)*2*damps);
downright_damp=diag(-dampval);
A_rotor(modes+1:2*modes,modes+1:2*modes)=downright_damp;


%%%%%%% Input and outputs of the system:

%Input 1: Impeller Input Disturbance Force
%Input 2: AMB 1 Input Force
%Input 3: AMB 2 Input Force
%Input 4: AMB 3 Input Force

%Output 1: AMB 1 Output displacement
%Output 2: AMB 2 Output displacement
%Output 3: AMB 3 Output displacement
%Output 4: AMB 1 Output Force
%Output 5: AMB 2 Output Force
%Output 6: AMB 3 Output Force

A=A_rotor;
outscale=diag([1/gap1 1/gap2 1/gap3 1/gap4 1/Fsat1 1/Fsat2 1/Fsat3 1/Fsat4]); %normalize outputs to maximum values

Bu=B_rotor(:,2:end); %w(impeller)
Bw=B_rotor(:,1); %u1 and u2 (amb)

C_rotor=[C_rotor(1:end,:); zeros(1,size(A_rotor,1)); zeros(1,size(A_rotor,1)); zeros(1,size(A_rotor,1));zeros(1,size(A_rotor,1))]; %position (m) @ (journal %force (N) )  
C=outscale*C_rotor;
%            u1 u2 u3 u4
Du=outscale*[0  0  0 0 ;  %x1
             0  0  0 0;   %x2
             0  0  0 0 ;  %x3
             0  0  0 0 ;  %x4
             1  0  0 0 ;  %Famb1
             0  1  0 0 ;  %Famb2
             0  0  1 0 ;  %Famb3
             0  0  0 1 ;];%Famb4

Dw=outscale*[0 0 0 0 0 0 0 0]'; %w

I=eye(max(size(A)));
II=eye(8);
fm=wrangef;
rm=[fm;fm;fm;fm;fm;fm;fm;fm];
%Calculate the 2-norm of the system
%does not take into account sensor noise and due this load capacity migth
%be too conservative. 
for i=1:length(wrangef)
  omega=1i*wrangef(i); %calculate load capacity at each frequency
%   K=C*inv(omega*I-A);
  K=C/(omega*I-A);
  Gu=Du+K*Bu;
  Gw=Dw+K*Bw;
%   Gp=(II-Gu*inv(Gu'*Gu)*Gu')*Gw;
  Gp=(II-(Gu/(Gu'*Gu))*Gu')*Gw;
  p=max(abs(Gp));
  fm(i)=1/p; %maximum load is the inverse of the normalized value!
  rm(:,i)=abs(Gp*fm(i));
end

%Plot the maximum impeller force in the given frequency range
figure;  
plot(wrangef/2/pi,fm,'LineWidth',3)
xlabel('Frequency, Hz');
ylabel('Maximum impeller force, N');
title([filename,' RoBeDyn']);

% %%%%%SAVE HINF.dat the load capacity%%%%%%
% if sizee==1; %original bearing
% ei=fopen('Results\OL.dat','wt');
% elseif sizee==2; %40% longer bearing
% ei=fopen('Results\OL_BIG.dat','wt');
% end
% for ij=1:length(wrangef)
%     fprintf(ei,' %g',fm(ij));
% end
% fclose(ei);

end