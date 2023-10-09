
function Inp=INDATA_fan(Request)
%% Initialization
%-Model name--------------------------------------------%
Inp.model_title='Volcar rotor v 1';

%-Material--------------------------------------------%
E=2.1e11;
nuxy=0.3;
rho=7800;
Inp.Mat=[ 1,  E,     nuxy,   rho ];
     
% Masspoint and related calculations
Inp.MassPoints=[];
% MassPoints=[ID, node, mass, Jxx,   Jyy,  Jzz,   h,   L, dia_out, dia_in];


%% Adding sections
% [sectionLength Din Dout1 Mat1 Dout2 Mat2 ... ]; addSection 
Inp.Node = [1 0 0 0]; % adding the first node of particular shaft
[0.098, 0.025, 0.035, 1];addSection; % section L1
[0.022, 0.025, 0.035, 1];addSection; % section L2
[0.0035, 0.025, 0.038, 1]; addSection;% Shoulder L3
[0.0025, 0.025, 0.040, 1]; addSection; % shoulder L4
[0.020,  0.025, 0.040, 1, 0.086 ,2]; addSection; % collar with magnets L5
[0.0005, 0.025, 0.040, 1, 0.086 ,2]; addSection; % collar disc L6
[0.098,  0.025, 0.040, 1, 0.086 ,2]; addSection; % EM
[0.0005, 0.025, 0.040, 1, 0.086 ,2]; addSection; % collar disc R1
[0.020,  0.025, 0.040, 1, 0.086 ,2]; addSection; % collar with magnets R2
[0.0025, 0.025, 0.040, 1]; addSection; % shoulder R3
[0.022,  0.025, 0.038, 1]; addSection; % shoulder R4
[1e8 1e4 0]; addBearingMatrix % dummy bearing parameters || ACTUAL BEARINGS ADDED LATER
[0.007,  0.025, 0.038, 1]; addSection; % shoulder R5
[0.016,  0.025, 0.038, 1, 0.045, 2]; addSection; % shoulder R5
[0.018,  0.025, 0.038, 1]; addSection; % shoulder R5
Inp.MassPoints=[ Inp.MassPoints; size(Inp.MassPoints,1)+1, size(Inp.Node,1), 42,        0,      0,      0,  0, 0.100, 0.15,   0.080];% Housing

                    
                    
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
Inp.UB=[2, 60*0.01e-3, 0
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
