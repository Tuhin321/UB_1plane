function f_RGeomPlotShaded(Inp,alpha)
% ===================================================================
% Function f_RGeomPlotShaded plots the Rotor Geometry as a Surface plot
%
% Inp   = Model data (Structure array)
% alpha = surface transparency (0..1) 1=not transparent, 0=fully transparent
%         default: alpha=1
% 
% Hint: for nice shaded graphics use commands after plotting: 
% light; lighting gouraud; shading interp;
% 
% Written by Jussi Sopanen, LUT, 2004
% ===================================================================

% Modifications:
% - 23.9.2014 Multiple colors for plot when using layered model by E.
%   Sikanen

if nargin==1; alpha=1; end;% no transparency

Node=Inp.Node; % copying node matrix

figure
set(gcf,'Units','normalized')		  % Units normalized (always)
set(gcf,'Position',[0.1 0.1 0.8 0.8]) % Position set to given
set(gcf,'Color',[1 1 1])			  % Background color white (always)
title(Inp.model_title,'Fontsize',14)

% Drawing elements----------------------------------------------------
for ii=1:size(Inp.Elem,1)
    if ~Inp.Elem(ii,3) == 0  % If not mass element

        %for jj=1:size(Inp.Node,1)
        %    % Getting row of I node from Node matrix
        %    if Inp.Elem(ii,2)==Inp.Node(jj,1)
        %        Iindx=jj; 
        %    end
        %    % Getting row of J node from Node matrix
        %    if Inp.Elem(ii,3)==Inp.Node(jj,1)
        %        Jindx=jj; 
        %    end
        %end
        % same than for-if loop
        Iindx= (Inp.Elem(ii,2)==Inp.Node(:,1));
        Jindx= (Inp.Elem(ii,3)==Inp.Node(:,1));
          
        % Element start and end point
        startp=[Inp.Node(Iindx,2), Inp.Node(Iindx,3), Inp.Node(Iindx,4)];
        endp=[Inp.Node(Jindx,2), Inp.Node(Jindx,3), Inp.Node(Jindx,4)];
        
        % element's rotation according to local x axis
        %thetaIJ=[Node(Iindx,5), Node(Jindx,5) Node(Iindx,6), Node(Jindx,6) Node(Iindx,7), Node(Jindx,7)];
        thetaIJ=[0 0 0 0 0 0]; % estimated to be zero, further on torsional vibration mode could be plotted using 3D model....
        
        %ind=find(Inp.Real(:,1)==Inp.Elem(ii,5)); % searching for corresponding row of real constant
        %dia=Inp.Real(ind,5);    % diameter
        % similar than find previously, but faster?
        d_out=Inp.Real(Inp.Real(:,1)==Inp.Elem(ii,5),5);    % outer diameter
        d_in=Inp.Real(Inp.Real(:,1)==Inp.Elem(ii,5),6);    % inner diameter
        % Modified 24.3.2010 6. column of Real matrix is inner diameter (previously 5. column -->H
        % and 6. column-->B)
        
        % Drawing element (cylinder as surfaces)
        %    shaft3D(startp,endp,thetaIJ,dia,dia_in,    rgb,alpha,segm)
        if isfield(Inp,'ElemColor')
            [TI]=shaft3D(startp,endp,thetaIJ,d_out,d_in,Inp.ElemColor(ii,:),alpha);
        else
            [TI]=shaft3D(startp,endp,thetaIJ,d_out,d_in,[0 1 1],alpha);
        end
        
    end
end

% Drawing disks-------------------------------------------------------
% MassPoints=[ID, node, mass, Jxx,   Jyy,  Jzz,   h,   L, dia_out, dia_in];
for i=1:size(Inp.MassPoints,1)
    
    dia=Inp.MassPoints(i,9); % disk diameter
    dia_in=Inp.MassPoints(i,10); % disk inner diameter
    h=Inp.MassPoints(i,7); % offset
    L=Inp.MassPoints(i,8); % width
    % searching row corresponding to node
    ind=find(Inp.Node(:,1)==Inp.MassPoints(i,2));
    startp=Inp.Node(ind,2:4)'+TI'*[-L/2+h 0 0 ]';
    endp=Inp.Node(ind,2:4)'+TI'*[L/2+h 0 0 ]';
    
    % Drawing cylinders
    %shaft3D(startp,endp,   thetaIJ,dia,dia_in,    rgb,alpha,segm)
    shaft3D( startp,endp,zeros(6,1),dia,dia_in,[1 0 0],alpha,30);

    % drawing lines of node inner diameter to start point
    gx=[Node(ind,2) Node(ind,2) Node(ind,2) Node(ind,2)
        startp(1)   startp(1)   startp(1)   startp(1)];
    gy=[Node(ind,3) Node(ind,3) Node(ind,3) Node(ind,3)
        startp(2)+dia_in/2   startp(2)-dia_in/2   startp(2)   startp(2)];
    gz=[Node(ind,4) Node(ind,4) Node(ind,4) Node(ind,4)
        startp(3)   startp(3) startp(3)+dia_in/2   startp(3)-dia_in/2   ];
    % drawing lines of node inner diameter to end point
    plot3(gx,gy,gz,'-r')
    gx=[Node(ind,2) Node(ind,2) Node(ind,2) Node(ind,2)
        endp(1)   endp(1)   endp(1)   endp(1)];
    gy=[Node(ind,3) Node(ind,3) Node(ind,3) Node(ind,3)
        endp(2)+dia_in/2   endp(2)-dia_in/2   endp(2)   endp(2)];
    gz=[Node(ind,4) Node(ind,4) Node(ind,4) Node(ind,4)
        endp(3)   endp(3) endp(3)+dia_in/2   endp(3)-dia_in/2   ];
    plot3(gx,gy,gz,'-r')
    % drawing horizontal lines
    gx=[startp(1) startp(1) startp(1) startp(1)
        endp(1)   endp(1)   endp(1)   endp(1)];
    gy=[startp(2)+dia_in/2 startp(2)-dia_in/2 startp(2) startp(2)
        endp(2)+dia_in/2   endp(2)-dia_in/2   endp(2)   endp(2)];
    gz=[startp(3) startp(3) startp(3)+dia_in/2 startp(3)-dia_in/2
        endp(3)   endp(3) endp(3)+dia_in/2   endp(3)-dia_in/2   ];
    plot3(gx,gy,gz,'-r')
end

xlabel('X','FontSize',13)
ylabel('Y','FontSize',13)
zlabel('Z','FontSize',13)
grid on

axis equal;  

end


function [TI]=shaft3D(startp,endp,thetaIJ,dia,dia_in,rgb,alpha,segm)
% This function draws an 3D shaft
% shaft3D(startp,endp,thetaIJ,dia,rgb,segm)
% startp = shaft start point [x1, y1, z1]
% endp   = shaft end point [x2, y2, z2]
% thetaIJ= rotations of nodes ([RotXI, RotXJ, RotYI, RotYJ, RotZI, RotZJ])
% dia_out= outer diameter of the shaft
% dia_in = inner diameter of the shaft
% rgb    = color of the shaft
%          E.g. [1 1 0]=yellow, [1 0 1]=magenta, [0 1 1]=cyan,
%               [1 0 0]=red,    [0 1 0]=green,   [0 0 1]=blue,
%               [1 1 1]=white,  [0 0 0]=black
% segm  = number of segments around circumference (Default=14)
% 
% Written by Jussi Sopanen LUT 2004

% if number of segments undefined, default of 14 is used
if nargin==7; segm=14; end

% creating cylinder geometry (r=0.5, L=1)
[xx yy zz]=cylinder([0.5 0.5], segm);

% arrow length
L=sqrt((endp(1)-startp(1))^2+(endp(2)-startp(2))^2+(endp(3)-startp(3))^2);

% scaling to proper length
zz1=zz*L;

% scaling outer diameter
xx1=xx*dia; yy1=yy*dia;

% scaling inner diameter, if dia_in > 0
if dia_in>0
    xx2=xx*dia_in; yy2=yy*dia_in;
    % updating coordinates
    xx1=[xx1 xx2 ]; yy1=[yy1 yy2 ]; zz1=[zz1 zz1 ];
end
% mapping according to global coordinate system
v=[0 1 0]';
theta=pi/2;
% Rodriques' equation
A=eye(3)+skew(v)*sin(theta)+2*skew(v)^2*sin(theta/2)^2;

% calculating transformation matrices (using directional cosine expressions)
TI=T3D(startp,endp); % for I node
TJ=T3D(startp,endp); % for J node

% Transformation matrices from rotations of I and J node
theta1=thetaIJ(1); theta2=thetaIJ(3); theta3=thetaIJ(5); 
A123_I=[cos(theta2)*cos(theta3)  sin(theta1)*sin(theta2)*cos(theta3)-cos(theta1)*sin(theta3)  cos(theta1)*sin(theta2)*cos(theta3)+sin(theta1)*sin(theta3)
        cos(theta2)*sin(theta3)  sin(theta1)*sin(theta2)*sin(theta3)+cos(theta1)*cos(theta3)  cos(theta1)*sin(theta2)*sin(theta3)-sin(theta1)*cos(theta3)
        -sin(theta2)             sin(theta1)*cos(theta2)                                      cos(theta1)*cos(theta2)];
theta1=thetaIJ(2); theta2=thetaIJ(4); theta3=thetaIJ(6); 
A123_J=[cos(theta2)*cos(theta3)  sin(theta1)*sin(theta2)*cos(theta3)-cos(theta1)*sin(theta3)  cos(theta1)*sin(theta2)*cos(theta3)+sin(theta1)*sin(theta3)
        cos(theta2)*sin(theta3)  sin(theta1)*sin(theta2)*sin(theta3)+cos(theta1)*cos(theta3)  cos(theta1)*sin(theta2)*sin(theta3)-sin(theta1)*cos(theta3)
        -sin(theta2)             sin(theta1)*cos(theta2)                                      cos(theta1)*cos(theta2)];

% mapping cylinder using transformation matrices
for ii=1:length(xx1)
    
    % calculating new coordinates for I node
    rotv1=A123_I*TI'*(A*[xx1(1,ii) yy1(1,ii) zz1(1,ii)]')+[startp(1) startp(2) startp(3)]';
    xr1(1,ii)=rotv1(1);
    yr1(1,ii)=rotv1(2);
    zr1(1,ii)=rotv1(3);

    % calculating new coordinates for J node
    rotv1=A123_J*TJ'*(A*[xx1(2,ii) yy1(2,ii) zz1(2,ii)-L]')+[endp(1) endp(2) endp(3)]';
    xr1(2,ii)=rotv1(1);
    yr1(2,ii)=rotv1(2);
    zr1(2,ii)=rotv1(3);
end

% coordinates of cylinder surfaces to patch objects
% cylinder outer surface
[fac, vert]= surf2patch(xr1(:,1:segm+1)', yr1(:,1:segm+1)', zr1(:,1:segm+1)');
B=repmat(rgb, size(vert,1), 1); % rgb colors
patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp','FaceAlpha', alpha)
hold on
if dia_in>0; % if inner D > 0, drawing also inner cylinder
    [fac, vert]= surf2patch(xr1(:,segm+2:2*segm+2)', yr1(:,segm+2:2*segm+2)', zr1(:,segm+2:2*segm+2)');
    B=repmat(rgb, size(vert,1), 1); %rgb colors
    patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp','FaceAlpha', alpha)
    
    % changing number of segments if inner diameter > 0 for drawing the end of cylinder
    segm=segm*2+2;
end;

% cylinder lower face
vert=[xr1(1,1:segm)', yr1(1,1:segm)', zr1(1,1:segm)']; %vertices
fac=1:1:segm; % faces
B=repmat(rgb, size(vert,1), 1); %rgb colors
patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp','FaceAlpha', alpha)

% cylinder upper face
vert=[xr1(2,1:segm)', yr1(2,1:segm)', zr1(2,1:segm)'];
% fac = same than previously
patch('Faces',fac,'Vertices', vert,'FaceVertexCData',B,'facecolor', 'interp','FaceAlpha', alpha)

end


% -----------------------------------------------------------------------------
function [T]=T3D(startp,endp)
% T3D returns transformation matrix of element's local coordinate system
% to global coordinate system.
% function [T]=T3D(start,end)
%
% startp = start point
% endp   = end point
% 
% Written by Jussi Sopanen, LUT/IMVe, 2002
% ======================================================================================

theta=0;
% Node coordinates of element
X1=startp(1);
X2=endp(1);
Y1=startp(2);
Y2=endp(2);
Z1=startp(3);
Z2=endp(3);
% length
L=sqrt((X2-X1)^2+(Y2-Y1)^2+(Z2-Z1)^2);
% projected length to xy plane
Lxy=sqrt(L^2-(Z2-Z1)^2);

% Error limits for trigonometric functions
d=0.0001*L;
% Transformation terms
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
% Sub matrix
T=[  C1*C2             S1*C2            S2
   (-C1*S2*S3-S1*C3) (-S1*S2*S3+C1*C3)  S3*C2
   (-C1*S2*C3+S1*S3) (-S1*S2*C3-C1*S3)  C3*C2];

end

