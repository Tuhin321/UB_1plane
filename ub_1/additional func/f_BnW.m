function f_BnW(Inp,xxplot,pnum,ltype)
% ===================================================================
% Function f_RGeomPlotWireFr plots Rotor Geometry as a WireFRame plot
% f_RGeomPlotWireFr(Inp,xxplot,pnum,ltype)
% Parameters:
% -----------------
% Inp     = Model input data (Structural array)
% xxplot  = component plot [E N B D M], where
%           E = plot elements (0=off,1=on)
%           N = plot nodes (0=off,1=on)
%           B = plot bearings (0=off,1=on)
%           D = plot constraints (0=off,1=on)
%           M = plot masspoints (0=off,1=on)
% pnum    = plot element and node numbers [E, N] (0=off,1=on)
% ltype   = cell array, which contains the line type and linewidth of the
%           element line e.g. {'-b',3}


% Written by Jussi Sopanen, LUT, 2004 (version 2.1 2006)
% ===================================================================

% default settings
if nargin<2
    xxplot=[1 1 1 1 1]; % plotting all components
end
if nargin<3 %
    pnum=[0 1]; % node numbers are plotted but no element numbers
end
if nargin<4;
    ltype={'-k',3}; % type of rotor center line
end

% Resaving Input variables
Node=Inp.Node;
Elem=Inp.Elem;
Disp=Inp.Disp;
Real=Inp.Real;

% new figure
figure   
% figure settings
set(gcf,'Units','normalized')		   % Units normalized (always)
set(gcf,'Position',[0.1 0.1 0.8 0.8],'renderer', 'Painter')  % Position set to given
% set(gcf,'Position',[100 100 5000 1000])  % Position set to given
set(gcf,'Color',[1 1 1])			   % Background color white (always)

% Calculating model dimensions--------------------------------------
model_size=[max(Node(:,2))-min(Node(:,2))  max(Node(:,3))-min(Node(:,3)) max(Node(:,4))-min(Node(:,4))];
% node number offset (position of node number=(xnode,ynode+nnos,znode+nnos)
nnos=max(model_size)*0.02; % 2 maximum model size

% counter for legend
NumLeg=0;

% Drawing elements----------------------------------------------------
if xxplot(1)==1
    
    for ii=1:size(Elem,1)
        
        if ~Elem(ii,3) == 0  % If not mass element
            %for jj=1:size(Node,1)
            %    % Getting row of I node from Node matrix
            %    if Elem(ii,2)==Node(jj,1)
            %        Iindx=jj; 
            %    end
            %    % Getting row of I node from Node matrix
            %    if Elem(ii,3)==Node(jj,1)
            %        Jindx=jj; 
            %    end
            %end
            % these lines replaces previous for-if loop
            Iindx= (Elem(ii,2)==Node(:,1));
            Jindx= (Elem(ii,3)==Node(:,1));
            
            % Element center line ltype according to input parameters
            h_leg(1)=plot3( [Node(Iindx,2) Node(Jindx,2)], [Node(Iindx,3) Node(Jindx,3)],...
                [Node(Iindx,4) Node(Jindx,4)],ltype{1},'linewidth',ltype{2});
            t_leg{1}='Centerline'; % legend text
            hold on
            
            %---------------------------------------------------------
            % Circle is drawn to node points to represent the thickness of axis
            % Should be working on general cases,
            % node points could be freely anywhere, not just along X axis
            
            %ind=find(Real(:,1)==Elem(ii,5)); % searching row corresponding real constant
            %dia=Real(ind,5); % Shaft diameter (Elem(ii,5)=Real constant number) 
            % old version
            %dia=Real(Real(:,1)==Elem(ii,5),5); % faster than previous...?
            %[xc,yc,zc]=cylinder(dia/2,20); % generating cylinder
            
            % Change 24.3.2010 Real constant column 5 is d_out and 6 d_in
            d_out=Real(Real(:,1)==Elem(ii,5),5); 
            d_in=Real(Real(:,1)==Elem(ii,5),6); 
            
            [xc,yc,zc]=cylinder(d_out/2,20); % generating cylinder
            if ~(d_in==d_out || d_in==0)
                [xci,yci,zci]=cylinder(d_in/2,20); % generating cylinder
            end
            %  mapping cylinder according to global coordinate system
            v=[0 1 0]'; % rotation vector
            theta=pi/2; % amount of rotation
            A=eye(3)+skew(v)*sin(theta)+2*skew(v)^2*sin(theta/2)^2;  % Transformation matrix (Rodriques' equation)      
            % calculation transformation matrix from element's node points (using cosines)
            startp=[Node(Iindx,2) Node(Iindx,3) Node(Iindx,4)]';
            endp=[Node(Jindx,2) Node(Jindx,3) Node(Jindx,4)]';
            TI=T3D(startp,endp); 
            % calculating circle point coordinates
            for kk=1:length(xc)
                % cylinder outside
                rotv1=(startp) + TI'*(A*[xc(1,kk) yc(1,kk) zc(1,kk)]');
                xr1(1,kk)=rotv1(1); yr1(1,kk)=rotv1(2); zr1(1,kk)=rotv1(3);
                rotv2=(endp) + TI'*(A*[xc(1,kk) yc(1,kk) zc(1,kk)]');
                xr2(1,kk)=rotv2(1); yr2(1,kk)=rotv2(2); zr2(1,kk)=rotv2(3);                
                
                % inside
                if ~(d_in==d_out || d_in==0)
                    rotv1=(startp) + TI'*(A*[xci(1,kk) yci(1,kk) zci(1,kk)]');
                    xr1i(1,kk)=rotv1(1); yr1i(1,kk)=rotv1(2); zr1i(1,kk)=rotv1(3);
                    rotv2=(endp) + TI'*(A*[xci(1,kk) yci(1,kk) zci(1,kk)]');
                    xr2i(1,kk)=rotv2(1); yr2i(1,kk)=rotv2(2); zr2i(1,kk)=rotv2(3);
                end
            end
            % Drawing circle outsides
            plot3(xr1,yr1,zr1,'-k','linewidth',1.5)
            plot3(xr2,yr2,zr2,'-k','linewidth',1.5)
            % inside
            if ~(d_in==d_out || d_in==0)
                plot3(xr1i,yr1i,zr1i,'--k','linewidth',1.5)
                plot3(xr2i,yr2i,zr2i,'--k','linewidth',1.5)
            end
            % drawing lines between circles at 90 degree spacing
            for kk=1:5:20
                h_leg(2)=plot3([xr1(kk) xr2(kk)],[yr1(kk) yr2(kk)],[zr1(kk) zr2(kk)],'-k','linewidth',1.5);
                if ~(d_in==d_out || d_in==0)
                plot3([xr1i(kk) xr2i(kk)],[yr1i(kk) yr2i(kk)],[zr1i(kk) zr2i(kk)],'--k','linewidth',1.5)
                end
                t_leg{2}='Shaft'; % setting legend text
            end
        end % drawing beam elements

        % ploting element numbers
        if pnum(1)==1
            % minor offset between text and node (same x coordinate)
%             text(mean([Node(Iindx,2) Node(Jindx,2)]),...
%                 mean([Node(Iindx,3) Node(Jindx,3)]+nnos),...
%                 mean([Node(Iindx,4) Node(Jindx,4)]+nnos),...
%                 int2str(Elem(ii,1) ),'color',[0 0 1])
        end
    end
    
    % counter for legend
    NumLeg=NumLeg+2;
    
    
end % if xxplot(1)==1

%Drawing nodes and node numbers---------------------------------- 
if xxplot(2)==1
    % counter for legend
    NumLeg=NumLeg+1;
    for ii=1:size(Node,1)
        % node plotted as sphere
        h_leg(NumLeg)=plot3( Node(ii,2), Node(ii,3), Node(ii,4), 'ok','linewidth',2);
        hold on
        t_leg{NumLeg}='Node'; % legend text
        % plotting node numbers
        if pnum(2)==1
%             text(Node(ii,2),Node(ii,3)+nnos,Node(ii,4)+nnos,int2str(Node(ii,1) ),...
%                 'HorizontalAlignment','center','color',[0 0 0])
        end
    end
end

% Drawing bearings
if xxplot(3)==1
    
    % counter for legend
    if size(Inp.Bearing,1)>0; NumLeg=NumLeg+1; end;
    
    for ii=1:size(Inp.Bearing,2)
        % bearing diameter and width
        if (~isfield(Inp.Bearing(ii),'D') || isempty(Inp.Bearing(ii).D)) % if diameter is not defined
            dia=max(model_size)*0.05; % if not defined using 5% of model size
        else
            dia=Inp.Bearing(ii).D*1e-3; % transformation mm-->m
        end
        if (~isfield(Inp.Bearing(ii),'B') || isempty(Inp.Bearing(ii).B)) % if diameter is not defined
            B=dia/4;
        else
            B=Inp.Bearing(ii).B*1e-3; % transformation mm-->m
        end
        [xc,yc,zc]=cylinder(dia/2,20); % generating cylinder
        % mapping cylinder according to global coordinate system
        v=[0 1 0]'; % rotation vector
        theta=pi/2; % amount of rotation
        A=eye(3)+skew(v)*sin(theta)+2*skew(v)^2*sin(theta/2)^2;  % Rotation matrix (Rodriques' equation)
        % searching row corresponding to node
        ind=find(Node(:,1)==Inp.Bearing(ii).Inode);
        % bearing start and end point
        startp=[Node(ind,2) Node(ind,3) Node(ind,4)]'+[-B/2 0 0]';
        endp=[Node(ind,2) Node(ind,3) Node(ind,4)]'+[B/2 0 0]';        
        
        % calculating circle points using transformation matrix
        for kk=1:length(xc)
            % cylinder coordinates
            rotv1=(startp) + (A*[xc(1,kk) yc(1,kk) zc(1,kk)]');
            xr1(1,kk)=rotv1(1);
            yr1(1,kk)=rotv1(2);
            zr1(1,kk)=rotv1(3);
            rotv2=(endp) + (A*[xc(1,kk) yc(1,kk) zc(1,kk)]');
            xr2(1,kk)=rotv2(1);
            yr2(1,kk)=rotv2(2);
            zr2(1,kk)=rotv2(3);
        end
        % Drawing circles
        plot3(xr1,yr1,zr1,'-k','linewidth',2); hold on;
        plot3(xr2,yr2,zr2,'-k','linewidth',2)
        % drawing lines between circles at 90 degree spacing
        for kk=1:5:20
            h_leg(NumLeg)=plot3([xr1(kk) xr2(kk)],[yr1(kk) yr2(kk)],[zr1(kk) zr2(kk)],'-k','linewidth',2);
            t_leg{NumLeg}='Bearing'; % legend text
        end
        
    end
end %xxplot(3)==1


%CONSTRAINTS******************************************************
% going through Disp matrix [Node, Dof, Value] (Dof: 1=UX, 2=UY, 3=UZ, 4=ROTX, 5=ROTY, 6=ROTZ)
if xxplot(4)==1
    
    % counter for legend
    if size(Disp,1)>0; NumLeg=NumLeg+1; end;
    
    % preallocating string variables
    %cnums{size(Disp,1)}='';
    cnums{max(Disp(:,1))}='';
    
    for ii=1:size(Disp,1)
        %for jj=1:size(Node,1)
        %    % Getting row of Disp node from Node matrix
        %    if Disp(ii,1)==Node(jj,1) %
        %        Dindx=jj; 
        %    end
        %end
        % same than previous for-if loop
        Dindx=Disp(ii,1)==Node(:,1); %
        
        % square side length
        a=nnos; % using same distance than node-number-offset
        % node point coordinates
        x1=Node(Dindx,2);
        y1=Node(Dindx,3);
        z1=Node(Dindx,4);
        
        % saving string variable from constraint number, e.g. 1256
        cnums{Disp(ii,1)}=strcat(cnums{Disp(ii,1)},num2str(Disp(ii,2)));
        % saving constraint coordinates for drawing the constraints
        numcoord(Disp(ii,1),:)=[x1 y1 z1];
       
        % drawing squares
        % xy plane
        h_leg(NumLeg)=plot3([x1-a/2 x1+a/2 x1+a/2 x1-a/2 x1-a/2], [y1-a/2 y1-a/2 y1+a/2 y1+a/2 y1-a/2],...
            [z1 z1 z1 z1 z1],'k','linewidth',2);
        hold on
        % xz plane
        plot3([x1-a/2 x1+a/2 x1+a/2 x1-a/2 x1-a/2], [y1 y1 y1 y1 y1],...
            [z1-a/2 z1-a/2 z1+a/2 z1+a/2 z1-a/2],'k','linewidth',2);
        
        t_leg{NumLeg}='Constraint'; % legend text
    end
    
    % plotting numbers
    for ii=1:size(Disp,1)
        text(numcoord(Disp(ii,1),1)+0*nnos,numcoord(Disp(ii,1),2)-nnos,numcoord(Disp(ii,1),3)-nnos,cnums{Disp(ii,1)},...
                'HorizontalAlignment','center','color','k')
    end
       
end

% Drawing amss points
% MassPoints=[ID, node, mass, Jxx,   Jyy,  Jzz,   h,   L, dia_out, dia_in];
if xxplot(5)==1

    % counter for legend
    if size(Inp.MassPoints,1)>0; NumLeg=NumLeg+1; end;
    
    for ii=1:size(Inp.MassPoints,1)
        
        ltm='-k'; % setting ltype for drawing mass points
        
        dia=Inp.MassPoints(ii,9); % disk diameter
        dia_in=Inp.MassPoints(ii,10); % disk inner diameter
        h=Inp.MassPoints(ii,7); % offset
        L=Inp.MassPoints(ii,8); % width
        
        [xc,yc,zc]=cylinder(dia/2,20); % generating sylinder
        % mapping cylinder according to global coordinate system
        v=[0 1 0]'; % Rotation vector
        theta=pi/2; % amount of rotation
        A=eye(3)+skew(v)*sin(theta)+2*skew(v)^2*sin(theta/2)^2;  % Rotation matrix (Rodriques' equation)
        % searching row corresponding to node
        ind=find(Node(:,1)==Inp.MassPoints(ii,2));
                
        % cylinder start and end point
        startp=[Node(ind,2) Node(ind,3) Node(ind,4)]'+...
               [-L/2+h 0 0]';
        endp=[Node(ind,2) Node(ind,3) Node(ind,4)]'+...
               [L/2+h 0 0]';        
        
        % calculating circle point coordinates using transformation matrix
        for kk=1:length(xc)
            % cylinder coordinates
            rotv1=(startp) + (A*[xc(1,kk) yc(1,kk) zc(1,kk)]');
            xr1(1,kk)=rotv1(1);
            yr1(1,kk)=rotv1(2);
            zr1(1,kk)=rotv1(3);
            rotv2=(endp) + (A*[xc(1,kk) yc(1,kk) zc(1,kk)]');
            xr2(1,kk)=rotv2(1);
            yr2(1,kk)=rotv2(2);
            zr2(1,kk)=rotv2(3);
        end
        % Drawing circles
        plot3(xr1,yr1,zr1,ltm,'color', [0 0 0], 'linewidth',1.5)
        plot3(xr2,yr2,zr2,ltm,'color', [0 0 0],'linewidth',1.5)
        % drawing circle inner diameter
        dr=dia_in/dia; % diameter ratio
        plot3(xr1,yr1*dr,zr1*dr,ltm,'color', [0 0 0],'linewidth',1.5)
        plot3(xr2,yr2*dr,zr2*dr,ltm,'color', [0 0 0],'linewidth',1.5)
        % drawing lines between circles at 90 degree spacing
        for kk=1:5:20
            h_leg(NumLeg)=plot3([xr1(kk) xr2(kk)],[yr1(kk) yr2(kk)],[zr1(kk) zr2(kk)],ltm,'color', [0 0 0],'linewidth',1.5);
            t_leg{NumLeg}='MassPoint';
        end
        % drawing lines from node to inner diameter
        gx=[Node(ind,2) Node(ind,2) Node(ind,2) Node(ind,2)
            startp(1)   startp(1)   startp(1)   startp(1)];
        gy=[Node(ind,3) Node(ind,3) Node(ind,3) Node(ind,3)
            startp(2)+dia_in/2   startp(2)-dia_in/2   startp(2)   startp(2)];
        gz=[Node(ind,4) Node(ind,4) Node(ind,4) Node(ind,4)
            startp(3)   startp(3) startp(3)+dia_in/2   startp(3)-dia_in/2   ];
        plot3(gx,gy,gz,ltm) % from start point to node
        gx=[Node(ind,2) Node(ind,2) Node(ind,2) Node(ind,2)
            endp(1)   endp(1)   endp(1)   endp(1)];
        gy=[Node(ind,3) Node(ind,3) Node(ind,3) Node(ind,3)
            endp(2)+dia_in/2   endp(2)-dia_in/2   endp(2)   endp(2)];
        gz=[Node(ind,4) Node(ind,4) Node(ind,4) Node(ind,4)
            endp(3)   endp(3) endp(3)+dia_in/2   endp(3)-dia_in/2   ];
        plot3(gx,gy,gz,ltm) % from end point to node
        % piirretään vaakaviivat sisähalkaisijan kohdalle
        gx=[startp(1) startp(1) startp(1) startp(1)
            endp(1)   endp(1)   endp(1)   endp(1)];
        gy=[startp(2)+dia_in/2 startp(2)-dia_in/2 startp(2) startp(2)
            endp(2)+dia_in/2   endp(2)-dia_in/2   endp(2)   endp(2)];
        gz=[startp(3) startp(3) startp(3)+dia_in/2 startp(3)-dia_in/2
            endp(3)   endp(3) endp(3)+dia_in/2   endp(3)-dia_in/2   ];
        plot3(gx,gy,gz,ltm)
    end % for masspoints
end %if

% figure general settings
grid on
xlabel('X','fontsize',13)
ylabel('Y','fontsize',13)
zlabel('Z','fontsize',13)
title(Inp.model_title,'fontsize',14)

% legend
% legend(h_leg,t_leg,'Location', 'South');
view(2);
axis equal; grid off; %view(2); axis vis3d;


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


