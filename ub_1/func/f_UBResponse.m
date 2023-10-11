function [maxAmp_Y_new, phaseAngle_Y_new, maxAmp_Z_new, phaseAngle_Z_new]=f_UBResponse(Inp,Request)

% funktio laskee ep‰tasapainon aiheuttaman steady state vasteen.
%         Vaatii Request-structuren
%         Request.UBResp.Spin_vec=[0:50:6000]*2*pi/60;
%         Request.UBResp.Type='Full'; % tai
%         Request.UBResp.Type='PseudoModal';
%         Request.UBResp.NroModes=30;
%         Request.UBResp.OutputNodes=[2, 6, 13, 16];

warning off MATLAB:nearlySingularMatrix

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

% kopioidaan indeksointimatriisi
nnum2cgdof=Matrx.nnum2cgdof;

% waitbar
% h = waitbar(0,'Computing Unbalance Response...');

% aikaparametri, h‰nt‰ l‰hes t‰ysi kierros
t=(-1.8*pi:2*pi/50:0); 
% alustetaan tulosmatriisit
MaxAmp=zeros(length(Request.UBResp.OutputNodes),length(Request.UBResp.Spin_vec));
U_vector=zeros(size(Matrx.Mc,1),size(t,2),length(Request.UBResp.Spin_vec));
% Phase angle for probe direction
ProbeAngle = Request.UBResp.ProxiProbeAngle; % Probe directions  

% Start of speed loop -----------------------------------------------------
for i=1:length(Request.UBResp.Spin_vec)
    
    % Roottorin matriisit (t‰ytyy m‰‰ritell‰ uudelleen joka kierroksella,
    % koska n‰ihin lis‰t‰‰n laakereiden nopeusriippuvaiset matriisit)
    Kc=Matrx.Kc; 
    Cc=Matrx.Cc;
    
    % Calculation of Bearing Matrices
    for ii=1:size(Inp.Bearing,2)
        
        
        % Calculation of bearing stiffness and damping matrix
        [Kb,Cb]=f_Bearing_Matrix(Inp.Bearing(ii),Request.UBResp.Spin_vec(i));
        
        % Laakerin j‰ykkyys- ja vaimennusmatriisin lis‰‰minen systeemin matriiseihin
        % Laakerinodet m‰‰ritetty muuttujilla INode ja JNode
        for iii=1:3
            for jjj=1:3
                
                % akselin solmu  (laakeri-laakeri)
                Igdof_r=nnum2cgdof(Inp.Bearing(ii).Inode,iii); %Haetaan globaali Dof riville
                Igdof_c=nnum2cgdof(Inp.Bearing(ii).Inode,jjj); %Haetaan globaali Dof sarakkeelle
                if ~ (Igdof_r==0 | Igdof_c==0)
                    Kc(Igdof_r,Igdof_c)=Kc(Igdof_r,Igdof_c) + Kb(iii,jjj); %lis‰t‰‰n laakerin termit j‰ykkyysmatriisiin
                    Cc(Igdof_r,Igdof_c)=Cc(Igdof_r,Igdof_c) + Cb(iii,jjj); %lis‰t‰‰n laakerin termit vaimennusmatriisiin
                end
                
                if ~Inp.Bearing(ii).Jnode==0
                    % tuennan solmu (tuenta-tuenta)
                    Jgdof_r=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); %Haetaan globaali Dof riville
                    Jgdof_c=nnum2cgdof(Inp.Bearing(ii).Jnode,jjj); %Haetaan globaali Dof sarakkeelle
                    if ~ (Jgdof_r==0 | Jgdof_c==0)
                        Kc(Jgdof_r,Jgdof_c)=Kc(Jgdof_r,Jgdof_c) + Kb(iii,jjj); %lis‰t‰‰n laakerin termit j‰ykkyysmatriisiin
                        Cc(Jgdof_r,Jgdof_c)=Cc(Jgdof_r,Jgdof_c) + Cb(iii,jjj); %lis‰t‰‰n laakerin termit vaimennusmatriisiin
                    end
                    
                    % laakeri-tuenta ristitermit (yl‰kolmio)
                    if ~ (Igdof_r==0 | Jgdof_c==0)
                        Kc(Igdof_r,Jgdof_c)=Kc(Igdof_r,Jgdof_c) - Kb(iii,jjj); %lis‰t‰‰n laakerin termit j‰ykkyysmatriisiin
                        Cc(Igdof_r,Jgdof_c)=Cc(Igdof_r,Jgdof_c) - Cb(iii,jjj); %lis‰t‰‰n laakerin termit vaimennusmatriisiin
                    end
                    % tuenta-laakeri ristitermit (alakolmio)
                    if ~ (Jgdof_r==0 | Igdof_c==0)
                        Kc(Jgdof_r,Igdof_c)=Kc(Jgdof_r,Igdof_c) - Kb(iii,jjj); %lis‰t‰‰n laakerin termit j‰ykkyysmatriisiin
                        Cc(Jgdof_r,Igdof_c)=Cc(Jgdof_r,Igdof_c) - Cb(iii,jjj); %lis‰t‰‰n laakerin termit vaimennusmatriisiin
                    end
                end
            end % jjj
        end %iii
    end % ii Bearing
    
    if strcmp(Inp.Support(ii).type,'RPM dependent') %
        % Calculation of Support Matrices
        for ii=1:size(Inp.Support,2)
            
            % Calculation of support stiffness and damping matrix
            [Ks,Cs]=f_Support_Matrix(Inp.Support(ii),Request.UBResp.Spin_vec(i));
            
            % Supportin j‰ykkyys- ja vaimennusmatriisin lis‰‰minen systeemin matriiseihin
            % Supportin node m‰‰ritetty laakerin JNode ja GROUND
            for iii=1:3
                for jjj=1:3
                    
                    % akselin solmu  (laakeri-laakeri)
                    Igdof_r=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); %Haetaan globaali Dof riville
                    Igdof_c=nnum2cgdof(Inp.Bearing(ii).Jnode,jjj); %Haetaan globaali Dof sarakkeelle
                    if ~ (Igdof_r==0 | Igdof_c==0)
                        Kc(Igdof_r,Igdof_c)=Kc(Igdof_r,Igdof_c) + Ks(iii,jjj); %lis‰t‰‰n supportin termit j‰ykkyysmatriisiin
                        Cc(Igdof_r,Igdof_c)=Cc(Igdof_r,Igdof_c) + Cs(iii,jjj); %lis‰t‰‰n supportin termit vaimennusmatriisiin
                    end
                end % jjj
            end %iii
        end % ii Bearing
    end

    

    % Calculation of steady state responce++++++++++++++++
    % laskenta t‰ysill‰ matriiseilla

    % pyˆrimisnopeus
    Spin=Request.UBResp.Spin_vec(i);

    % Gyromatriisi kerrotaan pyˆrimisnopeudella
    Gc=Spin*Matrx.Gc.*0;

    % Steady State tila Lalanne s. 102 kaikilla vapausasteilla
    AA=[Kc-Matrx.Mc*Spin^2     -Spin*(Cc+Gc)
    Spin*(Cc+Gc)           Kc-Matrx.Mc*Spin^2];
    
    % Voimavektori    
    fu=[Spin^2*Matrx.F2
        Spin^2*Matrx.F3];

    % ratkaistaan Steady State tila
    pp=AA\fu;


    % Muutetaan koordinaatit fyysisiksi siirtymisksi ja lasketaan
    % maksimiamplitudi
    % alustetaan matrisi
    U_vec_t=zeros(size(Matrx.Mc,1),size(t,2));

    for kkk=1:size(t,2)
        % Kaikki dof mukana
        U_vec_t(:,kkk)=( pp(1:size(Matrx.Mc,1))*sin(t(kkk)) + pp(size(Matrx.Mc,1)+1:2*size(Matrx.Mc,1))*cos(t(kkk)));
        %U_vec_t(:,kkk)=( pp(1:size(Matrx.Mc,1))*cos(t(kkk)) + pp(size(Matrx.Mc,1)+1:2*size(Matrx.Mc,1))*sin(t(kkk)));
        % palautettava 3-ulotteinen matriisi
        U_vector(:,kkk,i)=( pp(1:size(Matrx.Mc,1))*sin(t(kkk)) + pp(size(Matrx.Mc,1)+1:2*size(Matrx.Mc,1))*cos(t(kkk)));
        %U_vector(:,kkk,i)=( pp(1:size(Matrx.Mc,1))*cos(t(kkk)) + pp(size(Matrx.Mc,1)+1:2*size(Matrx.Mc,1))*sin(t(kkk)));     
    end
    


    for kkk=1:length(Request.UBResp.OutputNodes)
        %tulostettava solmu ja dof
        dofy=nnum2cgdof(Request.UBResp.OutputNodes(kkk),2);
        dofz=nnum2cgdof(Request.UBResp.OutputNodes(kkk),3);

        % Tallennetaan joka solmun max amplitudi vs. rpm matriisiin
%         MaxAmp(kkk,i)=max(sqrt(U_vec_t(dofy,:).^2+U_vec_t(dofz,:).^2));
        %Request.UBResp.OutputNodes=[2, 6, 13, 16];
        
        Psy_z(:,:) = pp(1:size(Matrx.Mc,1));
        Pcy_z(:,:) = pp(size(Matrx.Mc,1)+1:2*size(Matrx.Mc,1));
        
        Psy(kkk,:)= Psy_z(dofy,:);
        Pcy(kkk,:)= Pcy_z(dofy,:);
        Psz(kkk,:)= Psy_z(dofz,:);
        Pcz(kkk,:)= Pcy_z(dofz,:);
 
        if strcmp(Request.UBResp.OutputInProbeDirections,'Yes')
            theta = ProbeAngle; % Probe directions
        else
            theta = 0; % Global directions
        end
        
        Psy_probe(kkk,:)= Psy(kkk,:)*cos(theta)+Psz(kkk,:)*sin(theta);
        Pcy_probe(kkk,:)= Pcy(kkk,:)*cos(theta)+Pcz(kkk,:)*sin(theta);
        Psz_probe(kkk,:)= Psz(kkk,:)*cos(theta)-Psy(kkk,:)*sin(theta);
        Pcz_probe(kkk,:)= Pcz(kkk,:)*cos(theta)-Pcy(kkk,:)*sin(theta);
           
        maxAmp_Y_new(kkk,i)=sqrt((Psy_probe(kkk,:).^2)+(Pcy_probe(kkk,:).^2));
%         phaseAngle_Y_new(kkk,i) = atan(Psy_req_new(kkk,:)/Pcy_req_new(kkk,:));
        phaseAngle_Y_new(kkk,i) = atan2(Psy_probe(kkk,:),Pcy_probe(kkk,:));
            
        maxAmp_Z_new(kkk,i)=sqrt((Psz_probe(kkk,:).^2)+(Pcz_probe(kkk,:).^2));
%         phaseAngle_Z_new(kkk,i) = atan(Psz_req_new(kkk,:)/Pcz_req_new(kkk,:));
        phaseAngle_Z_new(kkk,i) = atan2(Psz_probe(kkk,:),Pcz_probe(kkk,:));

      
        if strcmp(Request.UBResp.PlotOrbitInYZplane,'Yes')
            % For plotting the orbits of the requested output nodes
            if Request.UBResp.Spin_vec(i) == Request.UBResp.SpinSpeed
                figure
                set(gcf,'Units','normalized')                   % Units normalized (always)
                set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
                set(gcf,'Color',[1 1 1])                        % Background color white (always)
                plot(U_vector(dofy,:,i),U_vector(dofz,:,i),'linewidth',2)
                axis equal
                xlabel('Y ( \mum)','fontsize',13)
                ylabel('Z ( \mum)','fontsize',13)
                grid on
                hold on
                title(['Orbits of Node- ' num2str(Request.UBResp.OutputNodes(kkk))],'fontsize',14)                
            end
        end
        
    end
    
    if strcmp(Request.UBResp.PlotRelativeAmp,'Yes')
        % For list of Rotor/Bearing/Sensor Nodes
        for bb=1:numel(Request.UBResp.RotorNodes)
            bdofy=nnum2cgdof(Request.UBResp.RotorNodes(bb),2);
            bdofz=nnum2cgdof(Request.UBResp.RotorNodes(bb),3);

            Psy_B(bb,:)= Psy_z(bdofy,:);
            Pcy_B(bb,:)= Pcy_z(bdofy,:);
            Psz_B(bb,:)= Psy_z(bdofz,:);
            Pcz_B(bb,:)= Pcy_z(bdofz,:);
            
            Psy_req_B(bb,i)= Psy_B(bb,:)*cos(ProbeAngle)+Psz_B(bb,:)*sin(ProbeAngle);
            Pcy_req_B(bb,i)= Pcy_B(bb,:)*cos(ProbeAngle)+Pcz_B(bb,:)*sin(ProbeAngle);
            Psz_req_B(bb,i)= Psz_B(bb,:)*cos(ProbeAngle)-Psy_B(bb,:)*sin(ProbeAngle);
            Pcz_req_B(bb,i)= Pcz_B(bb,:)*cos(ProbeAngle)-Pcy_B(bb,:)*sin(ProbeAngle);
        end
        
        % For list of Support Nodes
        for ss=1:numel(Request.UBResp.SupportNodes)
            sdofy=nnum2cgdof(Request.UBResp.SupportNodes(ss),2);
            sdofz=nnum2cgdof(Request.UBResp.SupportNodes(ss),3);

            Psy_S(ss,:)= Psy_z(sdofy,:);
            Pcy_S(ss,:)= Pcy_z(sdofy,:);
            Psz_S(ss,:)= Psy_z(sdofz,:);
            Pcz_S(ss,:)= Pcy_z(sdofz,:);
            
            Psy_req_S(ss,i)= Psy_S(ss,:)*cos(ProbeAngle)+Psz_S(ss,:)*sin(ProbeAngle);
            Pcy_req_S(ss,i)= Pcy_S(ss,:)*cos(ProbeAngle)+Pcz_S(ss,:)*sin(ProbeAngle);
            Psz_req_S(ss,i)= Psz_S(ss,:)*cos(ProbeAngle)-Psy_S(ss,:)*sin(ProbeAngle);
            Pcz_req_S(ss,i)= Pcz_S(ss,:)*cos(ProbeAngle)-Pcy_S(ss,:)*sin(ProbeAngle);  
        end
        
        % Relative displacement calculation
        Psy_rel(:,i)=Psy_req_B(:,i)-Psy_req_S(:,i);
        Pcy_rel(:,i)=Pcy_req_B(:,i)-Pcy_req_S(:,i);
        Psz_rel(:,i)=Psz_req_B(:,i)-Psz_req_S(:,i);
        Pcz_rel(:,i)=Pcz_req_B(:,i)-Pcz_req_S(:,i);

        % Relative maximum amplitude in probe directions
        maxAmp_Y_Rel(:,i)=sqrt((Psy_rel(:,i).^2)+(Pcy_rel(:,i).^2));
        maxAmp_Z_Rel(:,i)=sqrt((Psz_rel(:,i).^2)+(Pcz_rel(:,i).^2));
    end
    
    
    % waitbar
%     waitbar(i/length(Request.UBResp.Spin_vec))
    
end % Spin vec loop, indeksi i

% delete(h) 						% close waitbar

if 0
    % Tulostus
    figure 
    set(gcf,'Units','normalized')		% Units normalized (always)
    set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
    set(gcf,'Color',[1 1 1])			% Background color white (always)
    
    
    plot(Request.UBResp.Spin_vec*60/(2*pi),MaxAmp*1e6,'linewidth',2) %) %,'b-'
    
    % v=axis; v(3)=1e-8; axis(v); clear v
    xlabel('RPM','fontsize',13)
    ylabel('Maximum Amplitude (\mum)', 'fontsize',13)
    grid on
    
    % luodaan legendille string muuttuja
    for kkk=1:size(Request.UBResp.OutputNodes,2)
        str{kkk}=['Node ' num2str(Request.UBResp.OutputNodes(kkk)) ' Mag'];
    end
    
    legend(str)
    title({Inp.model_title; 'Unbalance Response'},'fontsize',14)
    
    
    if strcmp(Request.UBResp.OutputInGlobalYZDirections,'Yes') | strcmp(Request.UBResp.OutputInProbeDirections,'Yes')
        s=theta*(180/pi);
        % Plot for Amplitude vs. RPM
        figure
        subplot(2,1,1)
        set(gcf,'Units','normalized')                   % Units normalized (always)
        set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
        set(gcf,'Color',[1 1 1])                        % Background color white (always)
        plot(Request.UBResp.Spin_vec*60/(2*pi),maxAmp_Y_new*1e6,'linewidth',2) %) %,'b-'
        xlabel('RPM','fontsize',13)
        ylabel('Maximum Amplitude (\mum)', 'fontsize',13)
        grid on
        for kkk=1:size(Request.UBResp.OutputNodes,2)
            str{kkk}=['Node ' num2str(Request.UBResp.OutputNodes(kkk)) ' Mag'];
        end
        legend(str)
        if strcmp(Request.UBResp.OutputInProbeDirections,'Yes')
            title({Inp.model_title; ['Unbalance Response and Phase Angle in Probe1 Direction: ' num2str(s) ' Deg']},'fontsize',14)
        else
            title({Inp.model_title; 'Unbalance Response and Phase Angle in Global Y-Direction'},'fontsize',14)
        end
    
        % Plot for Phase Angle vs. RPM
        subplot(2,1,2)
        set(gcf,'Units','normalized')                   % Units normalized (always)
        set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
        set(gcf,'Color',[1 1 1])                        % Background color white (always)
        plot(Request.UBResp.Spin_vec*60/(2*pi),phaseAngle_Y_new*(180/pi),'linewidth',2)
        xlabel('RPM','fontsize',13)
        ylabel('Phase Angle (deg)', 'fontsize',13)
        grid on
        % String variable for the legend
        for kkk=1:size(Request.UBResp.OutputNodes,2)
            str{kkk}=['Node ' num2str(Request.UBResp.OutputNodes(kkk)) ' Mag'];
        end
        legend(str)
    
    
        % Plot for Amplitude vs. RPM
        figure
        subplot(2,1,1)
        set(gcf,'Units','normalized')                   % Units normalized (always)
        set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
        set(gcf,'Color',[1 1 1])                        % Background color white (always)
        plot(Request.UBResp.Spin_vec*60/(2*pi),maxAmp_Z_new*1e6,'linewidth',2) %) %,'b-'
        xlabel('RPM','fontsize',13)
        ylabel('Maximum Amplitude (\mum)', 'fontsize',13)
        grid on
        for kkk=1:size(Request.UBResp.OutputNodes,2)
            str{kkk}=['Node ' num2str(Request.UBResp.OutputNodes(kkk)) ' Mag'];
        end
        legend(str)
        s2=s-90;
        if strcmp(Request.UBResp.OutputInProbeDirections,'Yes')
            title({Inp.model_title; ['Unbalance Response and Phase Angle in Probe2 Direction: ' num2str(s2) ' Deg']},'fontsize',14)
        else
            title({Inp.model_title; 'Unbalance Response and Phase Angle in Global Z-Direction'},'fontsize',14)
        end
    
        % Plot for Phase Angle vs. RPM
        subplot(2,1,2)
        set(gcf,'Units','normalized')                   % Units normalized (always)
        set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
        set(gcf,'Color',[1 1 1])                        % Background color white (always)
        plot(Request.UBResp.Spin_vec*60/(2*pi),phaseAngle_Z_new*(180/pi),'linewidth',2)
        xlabel('RPM','fontsize',13)
        ylabel('Phase Angle (deg)', 'fontsize',13)
        grid on
        % String variable for the legend
        for kkk=1:size(Request.UBResp.OutputNodes,2)
            str{kkk}=['Node ' num2str(Request.UBResp.OutputNodes(kkk)) ' Mag'];
        end
        legend(str)
    end
    
    if strcmp(Request.UBResp.PlotRelativeAmp,'Yes')
        pa=ProbeAngle*(180/pi);
        % Plot for Relative Displacement w.r.t Support in Probe1 Direction vs. RPM
        figure
        set(gcf,'Units','normalized')                   % Units normalized (always)
        set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
        set(gcf,'Color',[1 1 1])                        % Background color white (always)
        plot(Request.UBResp.Spin_vec*60/(2*pi),maxAmp_Y_Rel*1e6,'linewidth',2) %) %,'b-'
        xlabel('RPM','fontsize',13)
        ylabel('Maximum Amplitude (\mum)', 'fontsize',13)
        grid on
        for kkk=1:numel(Request.UBResp.RotorNodes)
            str{kkk}=['Node ' num2str(Request.UBResp.RotorNodes(kkk)) ' Mag'];
        end
        legend(str)
        title({Inp.model_title; ['Relative Displacement w.r.t Support in Probe1 Direction: ' num2str(pa) ' Deg']},'fontsize',14)
        
        % Plot for Relative Displacement w.r.t Support in Probe2 Direction vs. RPM
        figure
        set(gcf,'Units','normalized')                   % Units normalized (always)
        set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
        set(gcf,'Color',[1 1 1])                        % Background color white (always)
        plot(Request.UBResp.Spin_vec*60/(2*pi),maxAmp_Z_Rel*1e6,'linewidth',2) %) %,'b-'
        xlabel('RPM','fontsize',13)
        ylabel('Maximum Amplitude (\mum)', 'fontsize',13)
        grid on
        for kkk=1:numel(Request.UBResp.RotorNodes)
            str{kkk}=['Node ' num2str(Request.UBResp.RotorNodes(kkk)) ' Mag'];
        end
        legend(str)
        pa2=pa-90;
        title({Inp.model_title; ['Relative Displacement w.r.t Support in Probe2 Direction: ' num2str(pa2) ' Deg']},'fontsize',14)
    end
end

warning on MATLAB:nearlySingularMatrix



