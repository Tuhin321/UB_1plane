function []=arrow3D(startp,endp,rgb,segm)
%This function draws an 3D arrow
% []=arrow3D(startp,endp,rgb,segm)
% startp = arrow start point [x1, y1, z1]
% endp   = arrow start point [x2, y2, z2]
% rgb    = color of the arrow
%          E.g. [1 1 0]=yellow, [1 0 1]=magenta, [0 1 1]=cyan,
%               [1 0 0]=red,    [0 1 0]=green,   [0 0 1]=blue,
%               [1 1 1]=white,  [0 0 0]=black
% segm  = number of segments around circumference (Default=14)
% 
% Written by Jussi Sopanen LUT/IMVe 2002

%jos segmenttien lukumäärää ei ole annettu, niin default on 14
if nargin==3; segm=14; end

%luodaan sylinterin geometria
[xx1 yy1 zz1]=cylinder([0.1 0.1], segm);

%luodaan conen  geometria
[xx2 yy2 zz2]=cylinder([0.2 0.0], segm);

% nuolen pituus
L=sqrt((endp(1)-startp(1))^2+(endp(2)-startp(2))^2+(endp(3)-startp(3))^2);

%skaalataan pituus oikeaksi
zz1=zz1*0.8*L;
zz2=(zz2-1)*0.2*L;

%skaalataan leveys
sf=4;
xx1=xx1*L/sf; yy1=yy1*L/sf;
xx2=xx2*L/sf; yy2=yy2*L/sf;

%käännetään torvet globaalin koord suuntaiseksi
v=[0 1 0]';
theta=pi/2;
%Rodriquesis kaava
A=eye(3)+skew(v)*sin(theta)+2*skew(v)^2*sin(theta/2)^2;


% lasketaan kiertomatriisi (suuntakosinien avulla)
T=T3D(startp,endp);

% lasketaan kierrot
for ii=1:length(xx1)
    for jj=1:2
        %cylinterin coord
        rotv1=T'*A*[xx1(jj,ii) yy1(jj,ii) zz1(jj,ii)]';
        xr1(jj,ii)=rotv1(1);
        yr1(jj,ii)=rotv1(2);
        zr1(jj,ii)=rotv1(3);
        %conen coord
        rotv2=T'*A*[xx2(jj,ii) yy2(jj,ii) zz2(jj,ii)]';
        xr2(jj,ii)=rotv2(1);
        yr2(jj,ii)=rotv2(2);
        zr2(jj,ii)=rotv2(3);
    end
end

%siirretään oikeisiin paikkoihin
xr1=xr1+startp(1); yr1=yr1+startp(2); zr1=zr1+startp(3); 
xr2=xr2+endp(1); yr2=yr2+endp(2); zr2=zr2+endp(3); 

% sylinterin coord patch objekteiksi
[fac, vert, c]= surf2patch(xr1', yr1', zr1');
B=repmat(rgb, size(vert,1), 1);
patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp')
hold on
% conen coord patch objekteiksi
[fac, vert, c]= surf2patch(xr2', yr2', zr2');
B=repmat(rgb, size(vert,1), 1);
patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp')

%nuolen pohja
vert=[xr1(1,1:segm)', yr1(1,1:segm)', zr1(1,1:segm)'];
fac=1:1:segm;
% fac=[1 2 3 4 5 6 7 8 9 10 11 12 13 14];
B=repmat(rgb, size(vert,1), 1);
patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp')

%conen pohja
vert=[xr2(1,1:segm)', yr2(1,1:segm)', zr2(1,1:segm)'];
fac=1:1:segm;
% fac=[1 2 3 4 5 6 7 8 9 10 11 12 13 14];
B=repmat(rgb, size(vert,1), 1);
patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp')

% axis equal
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% grid on


% -----------------------------------------------------------------------------
function [T]=T3D(startp,endp)
% T3D palauttaa siirtomatriisin elementin lokaalista koordinaatistosta
% globaaliin koordinaatistoon.
% function [T]=T3D(start,end)
%
% startp = alkupiste
% endp   = loppupiste
% 
% Written by Jussi Sopanen, LUT/IMVe, 2002
% ======================================================================================

% Alustetaan 3x3 nollamatriisi
T=zeros(3);
theta=0;
% Elementin solmujen koordinaatit
X1=startp(1);
X2=endp(1);
Y1=startp(2);
Y2=endp(2);
Z1=startp(3);
Z2=endp(3);
% pituus
L=sqrt((X2-X1)^2+(Y2-Y1)^2+(Z2-Z1)^2);
% xy-tasoon projisoitu pituus
Lxy=sqrt(L^2-(Z2-Z1)^2);

% Trigonometristen funktioiden tarkkuusraja
d=0.0001*L;
% Siirtomatriisin termit
if Lxy>d
    S1=(Y2-Y1)/Lxy;
    C1=(X2-X1)/Lxy;
else
    S1=0;
    C1=1;
end

S2=(Z2-Z1)/L;
S3=sin(theta);
C2=Lxy/L;
C3=cos(theta);
% Osamatriisi
T=[  C1*C2             S1*C2            S2
   (-C1*S2*S3-S1*C3) (-S1*S2*S3+C1*C3)  S3*C2
   (-C1*S2*C3+S1*S3) (-S1*S2*C3-C1*S3)  C3*C2];
