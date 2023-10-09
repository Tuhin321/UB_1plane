
function Inp=INDATA_fan(Request,ub_magnphase)
%% Initialization
%-Model name--------------------------------------------%
Inp.model_title='Overhung Jeffcott Rotor';

%-Material--------------------------------------------%
E=2.1e11;
nuxy=0.3;
rho=7800;
Inp.Mat=[ 1,  E,     nuxy,   rho ];
     
% Masspoint and related calculations
Inp.MassPoints=[];
% MassPoints=[ID, node, mass, Jxx,   Jyy,  Jzz,   h,   L, dia_out, dia_in];


% akselin halkaisijat
D1 =75.*1e-3;    % juoksupyörän kiinnityksen halkaisija
D2=80.0*1e-3;    % laakerihalkaisija
D3=115.0*1e-3;   % keskihalkaisija
D4=70.0*1e-3;    % halkaisija kytkimen kohdalla

%% Adding sections
% [sectionLength Din Dout1 Mat1 Dout2 Mat2 ... ]; addSection 
Inp.Node = [1 0 0 0]; % adding the first node of particular shaft
[0.040, 0, D1, 1];addSection; 
Inp.MassPoints=[ Inp.MassPoints; size(Inp.MassPoints,1)+1, size(Inp.Node,1), 80,      5.9,    3.1,    3.1,  0, 0.120, 0.850,  0.050];% Jouksupyörä
[0.072, 0, D1, 1]; addSection;% puhallin 
[0.188, 0, D2, 1]; addSection; % olake
[1e8 1e4 0]; addBearingMatrix % dummy bearing parameters || ACTUAL BEARINGS ADDED LATER
Inp.MassPoints=[ Inp.MassPoints; size(Inp.MassPoints,1)+1, size(Inp.Node,1), 42,        0,      0,      0,  0, 0.100, 0.15,   0.080];% Housing
[0.150, 0, D2, 1]; addSection; % laakeri etu
[0.325, 0, D3, 1]; addSection; % olake 
[0.134, 0, D2, 1]; addSection; % olake 
[1e8 1e4 0]; addBearingMatrix % dummy bearing parameters || ACTUAL BEARINGS ADDED LATER
Inp.MassPoints=[ Inp.MassPoints; size(Inp.MassPoints,1)+1, size(Inp.Node,1), 42,        0,      0,      0,  0, 0.100, 0.15,   0.080];% Housing
[0.140, 0, D2, 1]; addSection; % laakeri taka
[0.080 0, D2, 1]; addSection; % kytkin
Inp.MassPoints=[ Inp.MassPoints; size(Inp.MassPoints,1)+1, size(Inp.Node,1), 30,      0.2,    0.3,    0.3,  0, 0.130, 0.100,  0.070];% impeller% Kytkin
[0.071, 0, D4, 1]; addSection; % loppu
                    
                    
%% -Rajoitteet--------------------------------------------%
%-----Node, Dof, Value
% estetään aksiaalisiirtymät ja vääntö
for k=1:size(Inp.Node,1)
    if k==1
        Inp.Disp=[Inp.Node(k,1), 1, 0
                  Inp.Node(k,1), 4, 0
                  ];
    else
        Inp.Disp=[Inp.Disp
                  Inp.Node(k,1), 1, 0
                  Inp.Node(k,1), 4, 0
                  ];
    end
end

% Pistemassat inertiat määritetly globaalissa koordinaatistossa

%Inp.MassPoints=[]; 
if ~isfield(Inp,'MassPoints')
    Inp.MassPoints = [];
end

%-SpringDamper ------------------------------------------%    
Inp.SpringDamper=[];
Inp.Support(1).type='No Support';
Inp.Support(2).type='No Support';
%   SpringDamper=[ID, Inode, Jnode, Type, Dir,  Value];
% Inp.SpringDamper=[ 1,    18,     0,    1,  2,    200e6 % Tuennnan jäykkyys 
%                    2,    18,     0,    1,  3,    600e6
%                    3,    19,     0,    1,  2,    200e6
%                    4,    19,     0,    1,  3,    600e6
%                    5,    18,     0,    2,  2,    100e3 % vaimennuskertoimet arvioitu
%                    6,    18,     0,    2,  3,    100e3
%                    7,    19,     0,    2,  2,    100e3
%                    8,    19,     0,    2,  3,    100e3
%                  ];
% Type = 1 --> spring
% Type = 2 --> damper
% Dir 1=X, 2=Y, 3=Z
% Jos J node on 0, niin jousi on maassa kiinni


% Unbalance masses-------------------------------------------
% -- Node, value (kg*m), angle
Inp.UB=[Inp.MassPoints(1,2), ub_magnphase(1), ub_magnphase(2)*pi/180
%         8, 5*0.01e-3,      pi
%        11, 5*0.01e-3,      pi
];


%-lumped massa, jos lumpm=1-----------------------------%
Inp.lumpm=0;


%-Lähtötietojen muokkaus parempaan muotoon---------------------%
% Lisätään massapisteet Elementti matriisiin indeksointia yms. varten
for jj=1:size(Inp.MassPoints,1)
    
    ElemMP(jj,1)=size(Inp.Elem,1)+jj; % maksimielementtinumero +1
    ElemMP(jj,2)=Inp.MassPoints(jj,2); %Solmunumero
    ElemMP(jj,3:5)=0; % J node ja muut parametrit nollia
    
end
% Päivitetään 
if exist('ElemMP')
    Inp.Elem=[Inp.Elem
        ElemMP];
end

Inp.grav = [0 -9.81];

%-Voimat------------------------------------------------%
%-----Node, Dof, Value
% Force=[5,  1,   5000];
Inp.Force=[];

%-Cleanup------------------------------------------------------%
% Poistetaan väliaikaiset muuttujat
% Poistetaan väliaikaiset muuttujat
clearvars -except Inp
%---------------------------------------------------------------%
end
