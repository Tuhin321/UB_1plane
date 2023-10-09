% INPUT FILE FOR ROTOR-BEARING DYNAMICS CODE (RoBeDyn)
function Inp=IndataJeffcottGoodandBad(Request)
% 1= bad mass % 2= bad shaft dia % 3 = bad modulus of elasticity  % 4= bad bearing coefficients   % 5= bad support stiffness 
%-Model name--------------------------------------------%
Inp.model_title='Overhung Jeffcott Rotor';

badParameter=Request.badParameter;
if badParameter~=0
    percentOfTrue=1-Request.extentOfBadness;
    Request.percentOfTrue= percentOfTrue;% percentage of resemblance of bad parameter to the true model
    % for e.g. for 30 % variation, 0.3 is the extent of badness and the bad
    % parameter is 70% of the true parameter.
end
% akselin halkaisijat
if badParameter==2
    D1= 24.0*1e-3*percentOfTrue; 
else
    D1=24.0*1e-3;    % akselin halkaisija
end

% Materiaaliarvot
if badParameter==3
    E=2.1e11*percentOfTrue; % 30% error in modulus of elasticity
else
    E=2.1e11;
end
nuxy=0.3;

if Request.badParameter==1 % error in material density
        rho=7800*percentOfTrue;
else
    rho=7800;
end
ks=6*(1+nuxy)/(7+6*nuxy); %shear correction factor 

%-Real Constant-----------------------------------------% 
%--ID,A,Izz,Iyy,H,B,theta,istrn,Ixx,shearz,sheary
Inp.Real=[ 1, pi*D1^2/4, pi*D1^4/64, pi*D1^4/64, D1, D1, 0, 0, pi*D1^4/32, ks, ks %akseli D1
             ];

%keypointit
k1=[0.000,0,0]'; % takalaakeri
k2=[0.200,0,0]'; % etulaakeri 
k3=[0.240,0,0]'; % puhallin

%------------------------------------------------------------
% Solmujen ja elementtien generointi keypointtien välille.....
% [Node,Elem,MaxNodeNro,MaxElemNro]=
%                     Lmesh(startp,endp, Nelem,  MatID,RealID,   nodestart,elemstart,start_node)
[Node1,Elem1,MxN,MxE]=Lmesh(k1,k2         ,2,        1,1,                1,  1,1);   % disk
[Node2,Elem2,MxN,MxE]=Lmesh(k2,k3         ,1,        1,1,              MxN,MxE+1,0); % 


%-Nodes-------------------------------------------------%
%-----ID,   X,   Y,   Z
Inp.Node=[Node1; Node2]; % kootaan node matriisi
for ii=1:2; eval(strcat('clear Node', num2str(ii))); end % tuhotaan turhat muuttujat



Inp.Node = [Inp.Node
    MxN+1, 0, -0.03, -0.03
    MxN+2, 0.2, -0.03, -0.03]; % Y and Z  offset for better visualization




%-Elements----------------------------------------------%   
%------ID, I,  J, Mat, Real   
Inp.Elem=[Elem1; Elem2]; % kootaan Elem matriisi
for ii=1:2; eval(strcat('clear Elem', num2str(ii))); end % tuhotaan turhat muuttujat

%-Rajoitteet--------------------------------------------%
%-----Node, Dof, Value
Inp.Disp=[1, 1, 0
          1, 4, 0
          2, 1, 0
          2, 4, 0
          3, 1, 0
          3, 4, 0
          4, 1, 0
          4, 4, 0
%           1, 2, 0 %jäykästi tuettu
%           1, 3, 0
%           3, 2, 0
%           3, 3, 0
      ];
%



Inp.Disp = [Inp.Disp
    5, 1, 0
    5, 4, 0
    5, 5, 0
    5, 6, 0
    6, 1, 0
    6, 4, 0
    6, 5, 0
    6, 6, 0];



% Inp.Disp=[];     
%-PointMass Elements-------------------------------------%
% Pistemassojen inertiat määritelty globaalissa koordinaatistossa
% MassPoints=[ID, node, mass, Jxx,   Jyy,  Jzz,  ,h    L, dia, dia_in];
% if badParameter==1 % 30 % error in disc mass which is quite small so later on, 10x error will be introduced to the mass matrix as well
%     Inp.MassPoints=[ 1, 4,  6*percentOfTrue,  0.035*percentOfTrue, 0.055*percentOfTrue, 0.055*percentOfTrue,  0, 0.020, 0.150, 0 ];% puhallin
% else    
    Inp.MassPoints=[ 1, 4,  6,  0.035, 0.055, 0.055,  0, 0.020, 0.150, 0 ]; % puhallin                
% end
%

suppm1=10; % Support mass 1
suppm2=10; % Support mass 2

Inp.MassPoints=[Inp.MassPoints
    2, 5, suppm1, zeros(1,7)
    3, 6, suppm2, zeros(1,7)];



%-SpringDamper ------------------------------------------%                
%   SpringDamper=[ID, Inode, Jnode, Type, Dir,  Value];
Inp.SpringDamper=[ 
                 ];
%


    
if badParameter==5
    Inp.SpringDamper=[Inp.SpringDamper
                        1, 5, 0, 1, 2, 1e6*percentOfTrue
                        2, 5, 0, 1, 3, 1e6*percentOfTrue
                        3, 6, 0, 1, 2, 1e6*percentOfTrue
                        4, 6, 0, 1, 3, 1e6*percentOfTrue
                        5, 5, 0, 2, 2, 1e3*percentOfTrue
                        6, 5, 0, 2, 3, 1e3*percentOfTrue
                        7, 6, 0, 2, 2, 1e3*percentOfTrue
                        8, 6, 0, 2, 3, 1e3*percentOfTrue];
else
    Inp.SpringDamper=[Inp.SpringDamper
                        1, 5, 0, 1, 2, 1e6
                        2, 5, 0, 1, 3, 1e6
                        3, 6, 0, 1, 2, 1e6
                        4, 6, 0, 1, 3, 1e6
                        5, 5, 0, 2, 2, 1e3
                        6, 5, 0, 2, 3, 1e3
                        7, 6, 0, 2, 2, 1e3
                        8, 6, 0, 2, 3, 1e3]; 
end


% Type = 1 --> spring
% Type = 2 --> damper
% Dir 1=X, 2=Y, 3=Z
% Jos J node on 0, niin jousi on maassa kiinni


% Unbalance masses-------------------------------------------
% -- Node, value (kg*m), angle

if badParameter==6
    Inp.UB=[4, 0.8 *6*7.96e-5, pi/2];
else
     Inp.UB=[4, 6*7.96e-5, 0];
end
%-Materiaali--------------------------------------------%
%---- ID, E, nuxy, rho
Inp.Mat=[ 1,  E, nuxy, rho
          2,  E*100, nuxy, rho/10];
%-lumped massa, jos lumpm=1-----------------------------%
Inp.lumpm=0;

% modal damping ratios
%                      Mode#, damping
Inp.ModalDamping=[ 1,    10*1e-3   %1. taivutus
                   2,    10*1e-3 
                   3,    1.5e-3  %2. taivutus
                   4,    1.5e-3
                   5,    2e-3   %3. taivutus
                   6,    2e-3];



% Bearing parameters ------------------------------------%
if 1
     
    % Laakerin matriisien input-------------------------------------%
    jj=1;
    Inp.Bearing(jj).type='Bearing Matrix'; % Stringi, joka kuvaa laakeria
    Inp.Bearing(jj).Inode=1;     % Akselin solmu, johon laakeri on liitetty
    Inp.Bearing(jj).Jnode=5;     % Tuennan solmu, johon laakeri on liitetty

    % Stiffness matrix (directions in global coordinate system, x=axial dir. yz=radial)
    %     Kb=[kxx  kxy  kxz
    %         kyx  kyy  kyz
    %         kzx  kzy  kzz];
    Inp.Bearing(jj).Kb=[2e10    0       0
                         0    3e7       0
                         0      0     5e7];
    if badParameter==4 % 30% error in vertical stiffness and damping coefficient
        Inp.Bearing(jj).Kb=[2e10    0                    0
                            0    3e7*percentOfTrue       0
                            0       0     5e7*percentOfTrue];
    end
    % Damping matrix (directions in global coordinate system, x=axial dir. yz=radial)
    %     Cb=[cxx  cxy  cxz
    %         cyx  cyy  cyz
    %         czx  czy  czz];
    Inp.Bearing(jj).Cb=1e-5*Inp.Bearing(jj).Kb;
                        
    % Toinen laakeri
    Inp.Bearing(2)=Inp.Bearing(1);
    Inp.Bearing(jj).Inode=3;     % Akselin solmu, johon laakeri on liitetty
    Inp.Bearing(jj).Jnode=6;     % Tuennan solmu, johon laakeri on liitetty
                         
    % Laakeri matriisi input loppuu--------------------------------%
else
    Inp.Bearing=[];
end   






%-Lähtötietojen muokkaus parempaan muotoon---------------------%

% Lisätään massapisteet Elementti matriisiin indeksointia yms. varten
for jj=1:size(Inp.MassPoints,1)
    
    ElemMP(jj,1)=MxE+jj; % maksimielementtinumero +1
    ElemMP(jj,2)=Inp.MassPoints(jj,2); %Solmunumero
    ElemMP(jj,3:5)=0; % J node ja muut parametrit nollia
    
end
% Päivitetään 
if exist('ElemMP')
    Inp.Elem=[Inp.Elem
        ElemMP];
end

      
%-Voimat------------------------------------------------%
%-----Node, Dof, Value
% Force=[5,  1,   5000];
Inp.Force=[];

%-Cleanup------------------------------------------------------%
% Poistetaan väliaikaiset muuttujat
% Voi kommentoida pois debuggausta varte
% save ModelInp.mat Inp
% clear all
% load ModelInp.mat 
end
%---------------------------------------------------------------%


