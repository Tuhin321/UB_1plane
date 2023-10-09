function f_RespShapePlot(Inp,Request,U_vector,MaxAmp,newfigure,sf)
% Mode plotting of rotor forced vibration
% INPUT:
% Plotting rotor deflection (requires calculated response) +++++++++++++++
% Request.RespPlot.Spin_speed=3400*2*pi/60;
%
% Written by Jussi Sopanen LTY/IMVE 2004

if nargin==4
    newfigure=1;
end
% searching speed index
[V I]=find(Request.RespPlot.Spin_speed==Request.UBResp.Spin_vec);
FiiR=U_vector(:,:,I);

if nargin<5
    % default values
    % calculating scale factor for mode (scaling um --> m)
    ModDim=max(max(Inp.Node(:,2:4)))-min(min(Inp.Node(:,2:4)));
    if ModDim/max(MaxAmp(:,I)) > 10 && ModDim/max(MaxAmp(:,I)) < 100
        ScaleFactor=10;
    elseif ModDim/max(MaxAmp(:,I)) > 100 && ModDim/max(MaxAmp(:,I)) < 1000
        ScaleFactor=100;
    elseif ModDim/max(MaxAmp(:,I)) > 1000 && ModDim/max(MaxAmp(:,I)) < 10000
        ScaleFactor=1000;
    elseif ModDim/max(MaxAmp(:,I)) > 10000 && ModDim/max(MaxAmp(:,I)) < 100000
        ScaleFactor=10000;
    else
        ScaleFactor=100000;
    end
else
    ScaleFactor=sf; % using given value
end

% Calculation of Rotor Matrices------------------------------------------
% These matrices are constant with respect to rotation speed
[Matrx]=f_Rotor_Matrix(Inp);
% Contents of structure Matrx:
% Matrx.Kc=Kc; % Stiffness Matrix
% Matrx.Cc=Cc; % Damping matrix
% Matrx.Gc=Gc; % Gyroscopic matrix (must be multiplied with Spin Speed)
% Matrx.Mc=Mc; % Mass matrix
% Matrx.FR=FR; % Vector of external applied forces
% Matrx.F2=F2; % Unbalance vector (see documentation)
% Matrx.F3=F3; % Unbalance vector (see documentation)
% Matrx.nnum2gdof=nnum2gdof; % Indexing matrix of unconstrained dofs
% Matrx.nnum2cgdof=nnum2cgdof; % Indexing matrix of constrained dofs

% copying index matrix
nnum2cgdof=Matrx.nnum2cgdof;

% Plotting modes--------------------------------------------------------

Node=Inp.Node;
Elem=Inp.Elem;

if newfigure==1
    figure
    set(gcf,'Units','normalized')		% Units normalized (always)
    set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
    set(gcf,'Color',[1 1 1])			% Background color white (always)
end
% mm=Request.RespPlot.ModeNro; % Plotting mode number

% t=(0:2*pi/20:1.8*pi); % time parameter, full revolution
% sf=5000; %scale factor
% mode plot
for ii=1:size(Node,1)
    %nnum=Node(ii,1); % node number
    % global DOF corresponding to node number
    UXindx=nnum2cgdof(Node(ii,1),1);
    UYindx=nnum2cgdof(Node(ii,1),2);
    UZindx=nnum2cgdof(Node(ii,1),3);
    
    % Node's x displacement
    if ~UXindx==0 % if not constrained gdof
        u1(ii,:)=Node(ii,2)*ones(1,size(FiiR,2))+...
                 0*ScaleFactor*FiiR(UXindx,:); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        u1(ii,:)=Node(ii,2)*ones(1,size(FiiR,2));
    end
    % Node's y displacement
    if ~UYindx==0 % if not constrained gdof
        v1(ii,:)=Node(ii,3)*ones(1,size(FiiR,2))+...
                 ScaleFactor*FiiR(UYindx,:);    
    else
        v1(ii,:)=Node(ii,3)*ones(1,size(FiiR,2));
    end
    % Node's z displacement
    if ~UZindx==0 % if not constrained gdof
        w1(ii,:)=Node(ii,4)*ones(1,size(FiiR,2))+...
                 ScaleFactor*FiiR(UZindx,:);        
    else
        w1(ii,:)=Node(ii,4)*ones(1,size(FiiR,2));
    end

end

% Transposes of displacement matrices for plo3 command
u1=u1'; v1=v1'; w1=w1';

% Plotting circles
H=plot3(u1, v1, w1,'b','linewidth',1);
h_leg(1)=H(1); % saving handle to legendi
t_leg{1}='Whirl (\Omegat=-1.8\pi...0)';
hold on
% dot to nodes
mu=size(u1,1); 
h_leg(2)=plot3(u1(mu,:), v1(mu,:), w1(mu,:),'ok','markersize',6,'linewidth',2);
t_leg{2}='Node';

% Plotting beam elements as lines
t=0;
for ii=1:size(Elem,1)
    
    if ~Elem(ii,3) == 0 % If not mass element
        
        for jj=1:size(Node,1)
            % Getting row of I node from Node matrix
            if Elem(ii,2)==Node(jj,1)
                Iindx=jj; 
            end
            % Getting row of J node from Node matrix
            if Elem(ii,3)==Node(jj,1)
                Jindx=jj; 
            end
        end
        
        % global DOF corresponding to node number    
        UXindxI=nnum2cgdof(Node(Iindx,1),1);
        UYindxI=nnum2cgdof(Node(Iindx,1),2);
        UZindxI=nnum2cgdof(Node(Iindx,1),3);
        UXindxJ=nnum2cgdof(Node(Jindx,1),1);
        UYindxJ=nnum2cgdof(Node(Jindx,1),2);
        UZindxJ=nnum2cgdof(Node(Jindx,1),3);        
        
        % Node's x displacement
        if ~UXindxI==0 % if not constrained gdof
            uI=Node(Iindx,2)+ 0*ScaleFactor*FiiR(UXindxI,mu) ; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            uI=Node(Iindx,2);
        end
        if ~UXindxJ==0 % if not constrained gdof
            uJ=Node(Jindx,2)+ 0*ScaleFactor*FiiR(UXindxJ,mu) ;  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            uJ=Node(Jindx,2);
        end
        
        % Node's y displacement
        if ~UYindxI==0 % if not constrained gdof
            vI=Node(Iindx,3)+ScaleFactor*FiiR(UYindxI,mu) ; 
        else
            vI=Node(Iindx,3);
        end
        if ~UYindxJ==0 % if not constrained gdof
            vJ=Node(Jindx,3)+ScaleFactor*FiiR(UYindxJ,mu) ; 
        else
            vJ=Node(Jindx,3);
        end
        % Node's z displacement
        if ~UZindxI==0 % if not constrained gdof
            wI=Node(Iindx,4)+ScaleFactor*FiiR(UZindxI,mu) ; 
        else
            wI=Node(Iindx,4);
        end
        if ~UZindxJ==0 % if not constrained gdof
            wJ=Node(Jindx,4)+ScaleFactor*FiiR(UZindxJ,mu) ; 
        else
            wJ=Node(Jindx,4);
        end
        
        % Element as line
        plot3( [Node(Iindx,2) Node(Jindx,2)], [Node(Iindx,3) Node(Jindx,3)],...
            [Node(Iindx,4) Node(Jindx,4)],'-k','linewidth',1) % center line
        h_leg(3)=plot3( [uI uJ], [vI vJ], [wI wJ ],'-k','linewidth',2); % line from node to node
        t_leg{3}='Centerline';
        if ii==1
            plot3( [Node(Iindx,2) uI ], [Node(Iindx,3) vI], [Node(Iindx,4) wI],'-k','linewidth',2)
        else
            plot3( [Node(Jindx,2) uJ ], [Node(Jindx,3) vJ], [Node(Jindx,4) wJ],'-k','linewidth',2)
        end
        
    else % if mass element
        
        for jj=1:size(Node,1)
            % Getting row of I node from Node matrix
            if Elem(ii,2)==Node(jj,1)
                Iindx=jj; 
            end
        end
        
        % global DOF corresponding to node number    
        UXindxI=nnum2cgdof(Node(Iindx,1),1);
        UYindxI=nnum2cgdof(Node(Iindx,1),2);
        UZindxI=nnum2cgdof(Node(Iindx,1),3);
        
        % Node's x displacement
        if ~UXindxI==0 % if not constrained gdof
            uI=Node(Iindx,2)+  0*ScaleFactor*FiiR(UXindxI,mu);  %¤!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
            uI=Node(Iindx,2);
        end
        
        % Node's y displacement
        if ~UYindxI==0 % if not constrained gdof
            vI=Node(Iindx,3)+ScaleFactor*FiiR(UYindxI,mu) ; 
        else
            vI=Node(Iindx,3);
        end
        
        % Node's z displacement
        if ~UZindxI==0 % if not constrained gdof
            wI=Node(Iindx,4)+ScaleFactor*FiiR(UZindxI,mu) ; 
        else
            wI=Node(Iindx,4);
        end
        
        % Mass element as dot
        h_leg(4)=plot3( [uI ], [vI], [wI ],'or','linewidth',3,'markersize',10);
        t_leg{4}='MassPoint';
    end

end

NumLeg=size(h_leg,2)+1; % checking if legend values are 3 or 4

% Plotting bearing
for kkk=1:size(Inp.Bearing,2)

    for jj=1:size(Node,1)
        % Getting row of I node from Node matrix
        if Inp.Bearing(kkk).Inode==Node(jj,1)
            Iindx=jj;
        end
    end

    % global DOF corresponding to node number
    UXindxI=nnum2cgdof(Node(Iindx,1),1);
    UYindxI=nnum2cgdof(Node(Iindx,1),2);
    UZindxI=nnum2cgdof(Node(Iindx,1),3);

    % Node's x displacement
    if ~UXindxI==0 % if not constrained gdof
        uI=Node(Iindx,2)+  0*ScaleFactor*FiiR(UXindxI,mu);  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        uI=Node(Iindx,2);
    end

    % Node's y displacement
    if ~UYindxI==0 % if not constrained gdof
        vI=Node(Iindx,3)+ScaleFactor*FiiR(UYindxI,mu) ;
    else
        vI=Node(Iindx,3);
    end

    % Node's z displacement
    if ~UZindxI==0 % if not constrained gdof
        wI=Node(Iindx,4)+ScaleFactor*FiiR(UZindxI,mu) ;
    else
        wI=Node(Iindx,4);
    end


    % Plotting bearing as green square
    h_leg(NumLeg)=plot3( [uI ], [vI], [wI ],'sg','linewidth',3,'markersize',12);
    t_leg{NumLeg}='Bearing';
end


% Setting figure title
str0=['Rotor Deflection @ ' num2str(Request.RespPlot.Spin_speed*60/(2*pi)) ' rpm , Maximum Amplitude =' num2str(max(MaxAmp(:,I))*1e6,3) ' \mum' ];
str2={ Inp.model_title; str0};
title(str2,'Fontsize',14)
xlabel('X','fontsize',13)
ylabel('Y','fontsize',13)
zlabel('Z','fontsize',13)

% legend
legend(h_leg,t_leg,'Location', 'SouthEastOutside');




% figure settings
view(3);  axis vis3d; axis equal; axis off;
view3d rot;


% Drawing XYZ coordinate system
%f_draw_caxis(origloc, sz)
f_draw_caxis;
