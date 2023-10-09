function [FrqR,DampR,Matrx,criticalspeeds]=f_CriticalSM(Inp,Request,Stiffness)
% Function calculates critical frequnecies according to stiffness of support
% Options: - full, using full Kc, Mc Cc, Gc matrices
%         - pseudomodal, calculating n number of modes and calculation is based on them
%         
%         Requires Request struct
%         Request.Campbell.NroFreqs=6;
%         Request.Campbell.Spin_vec=[0:1000:10000]*2*pi/60;
%         Request.Campbell.Type='Full'; % tai
%         Request.Campbell.Type='PseudoModal';
%         Request.Campbell.NroModes=30;
         %Requires vector which includes stiffnesses of support.
         %Esimerkiksi
         %z = 7:9; d=1:9; %EXAMPLE HERE: z=7:9 -> 10^7 -> 10^9 stiffnesses,
         %10 per power.
         %stiffness_vec=cell2mat( arrayfun(@(i) d.*10^z(i),1:length(z),'UniformOutput',false) );
% 
% Input:   Inp
%          Request
%                        
% Returns frequencies, rpm and dampings to matrices, the same ones that are plotted

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

% initializing speed and damping matrices
FrqR=zeros(Request.Campbell.NroFreqs+4,length(Request.Campbell.Spin_vec));
DampR=zeros(Request.Campbell.NroFreqs+4,length(Request.Campbell.Spin_vec));

%h = waitbar(0,'Computing Campbell diagram...');
% Start of stiffness loop
for k=1:length(Stiffness)
% Start of speed loop
for i=1:length(Request.Campbell.Spin_vec)
    
    % Rotor matrices (need to redefine on every cycle of loop
    % since speed dependent matrices of bearings are added)
    
    Kc=Matrx.Kc; 
    Cc=Matrx.Cc;
    Kc2=Matrx.Kc; % Using Lalanne markings K* matrix (does not include cross-coupled terms of bearings)
    
    
    % Calculation of Bearing Matrices
    for ii=1:size(Inp.Bearing,2)
        
        
        % Calculation of bearing stiffness and damping matrix
       % [Kb,Cb]=f_Bearing_Matrix(Inp.Bearing(ii),Request.Campbell.Spin_vec(i));
       Kb=[2.021e7    0       0
             0    2.8445e8       0
             0      0     Stiffness(k)];
        Cb=Kb*0;
        % Adding bearing's stiffness and damping matrix into system matrices
        % Bearing nodes are defined using variables INode and JNode
        for iii=1:3
            % Forming K* matrix, which cross-coupled terms of bearing are removed
            % (Lalanne)
            Igdof=nnum2cgdof(Inp.Bearing(ii).Inode,iii); % node of shaft, getting global DOF for row
            if ~ (Igdof==0) % if not constrained
                Kc2(Igdof,Igdof)=Kc2(Igdof,Igdof) + Kb(iii,iii); % adding bearing's terms to stiffness matrix
            end
            % node of support
            if ~Inp.Bearing(ii).Jnode==0
                Jgdof=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); % % node of support, getting global DOF for row
                if ~ (Jgdof==0) % if not constrained
                    Kc2(Jgdof,Jgdof)=Kc2(Jgdof,Jgdof) + Kb(iii,iii); % adding bearing's terms to stiffness matrix
                end     
            end
            
            for jjj=1:3
                
                % node of shaft  (bearing-bearing)
                Igdof_r=nnum2cgdof(Inp.Bearing(ii).Inode,iii); % Getting global DOF for row
                Igdof_c=nnum2cgdof(Inp.Bearing(ii).Inode,jjj); % Getting global DOF for column
                if ~ (Igdof_r==0 | Igdof_c==0) 
                    Kc(Igdof_r,Igdof_c)=Kc(Igdof_r,Igdof_c) + Kb(iii,jjj); % adding bearing's terms to stiffness matrix
                    Cc(Igdof_r,Igdof_c)=Cc(Igdof_r,Igdof_c) + Cb(iii,jjj); % adding bearing's terms to damping matrix
                end
                
                if ~Inp.Bearing(ii).Jnode==0
                    % node of support (support-support)
                    Jgdof_r=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); % Getting global DOF for row
                    Jgdof_c=nnum2cgdof(Inp.Bearing(ii).Jnode,jjj); % Getting global DOF for column
                    if ~ (Jgdof_r==0 | Jgdof_c==0)
                        Kc(Jgdof_r,Jgdof_c)=Kc(Jgdof_r,Jgdof_c) + Kb(iii,jjj); % adding bearing's terms to stiffness matrix
                        Cc(Jgdof_r,Jgdof_c)=Cc(Jgdof_r,Jgdof_c) + Cb(iii,jjj); % adding bearing's terms to damping matrix
                    end
                    
                    % bearing-support cross-coupled terms (upper triangle)
                    if ~ (Igdof_r==0 | Jgdof_c==0)
                        Kc(Igdof_r,Jgdof_c)=Kc(Igdof_r,Jgdof_c) - Kb(iii,jjj); % adding bearing's terms to stiffness matrix
                        Cc(Igdof_r,Jgdof_c)=Cc(Igdof_r,Jgdof_c) - Cb(iii,jjj); % adding bearing's terms to damping matrix
                    end
                    % support-bearing cross-coupled terms (lower triangle)
                    if ~ (Jgdof_r==0 | Igdof_c==0)
                        Kc(Jgdof_r,Igdof_c)=Kc(Jgdof_r,Igdof_c) - Kb(iii,jjj); % adding bearing's terms to stiffness matrix
                        Cc(Jgdof_r,Igdof_c)=Cc(Jgdof_r,Igdof_c) - Cb(iii,jjj); % adding bearing's terms to damping matrix
                    end   
                end
                

            end % jjj
            
        end %iii
    end % ii Bearing
    
    % Calculation of natural frequencies
    if strcmp(Request.Campbell.Type,'Full')
        % calculation using full matrices
        
        % Gyromatrix multiplied by rotational speed
        Gc=Request.Campbell.Spin_vec(i)*Matrx.Gc;
        
        % Eigenvalue problem A*X=lambda*X (State Space mode)
        zeroX=zeros(size(Matrx.Mc));  eyeX=eye(size(Matrx.Mc));
        AA=[  zeroX                           eyeX
            -(Matrx.Mc)^(-1)*Kc      -(Matrx.Mc)^(-1)*(Cc+Gc)];
        
        % Solution of eigenmode problem
        [Fii2 d2]=eig(AA); % Gyromatrix included
        clear AA
        
        % Sorting matrices according to frequencies
        [Freq2, idx2]=sort(diag(d2)/(2*pi));  %Freq=alfa +/- beta*i 
        
        % Freq2 include complex eigenfrequencies
        % Searching only modes corresponding positive eigenfrequency
        kk2=0; jj2=0;
        nfreq=size(Freq2,1);
        while kk2 < Request.Campbell.NroFreqs+2  % +2 purely empirical
            jj2=jj2+1;
            % Imaginary part is equal to frequency Hz
            if imag(Freq2(jj2)) > 0.001 %
                %  i = index of rotational speed
                %  kk2 = counter for found frequency
                kk2=kk2+1; 
                FrqR(kk2,i)=imag(Freq2(jj2)); % Saving frequencies
                DampR(kk2,i)=-real(Freq2(jj2))/(abs(Freq2(jj2))); % Saving damping ratios
                %FiiR(:,kk2,i)=Fii2(nfreq/2+1:nfreq,jj2)/(Freq2(jj2)*2*pi); % Removing other half
            end 
        end
        clear kk2 jj2
        
    else
        % Calculation using speudo modal method
        
        % Solution of eigenmode problem
        [FiiN d]=eig(Kc2,Matrx.Mc);
        % Sorting matrices according to frequencies
        [Freq, idx]=sort(sqrt(diag(d))/(2*pi));
        
        %for k=1:length(idx)
        %   Fiiapu(:,k)=FiiN(:,idx(k));
        %end
        %FiiN=Fiiapu;
        %clear Fiiapu idx k d        
        % same than above, sorting modes according to frequencies
        FiiN=FiiN(:,idx); clear idx d;
                
        % Selecting numver of modes (Number of NroModes first eigenmodes)
        FiiNn=FiiN(:,1:Request.Campbell.NroModes);
        
        % Mass normalization
        for k=1:Request.Campbell.NroModes
            FiiNormN(:,k)=FiiNn(:,k)/sqrt(FiiNn(:,k)'*Matrx.Mc*FiiNn(:,k));
        end
        
        % Modal matrices
        Mm=FiiNormN'*Matrx.Mc*FiiNormN;
        Km=FiiNormN'*Kc2*FiiNormN + FiiNormN'*(Kc-Kc2)*FiiNormN; % jälkimmäinen termi sisältää laakerin ristitermit
        
        % Gyromatrix multiplied by rotational speed
        Gc=Request.Campbell.Spin_vec(i)*Matrx.Gc;
        % damping matrix + gyromatrix
        Cm=FiiNormN'*(Cc+Gc)*FiiNormN;
        
        % Eigenvalue problem A*X=lambda*X
        zeroX=zeros(size(Mm));  eyeX=eye(size(Mm));
        AA=[  zeroX               eyeX
            -(Mm)^(-1)*Km      -(Mm)^(-1)*Cm];
        
        [Fii2 d2]=eig(AA); % Solution of state space eigenvalue problem
        clear AA
        
        % Sorting matrices according to frequencies
        [Freq2, idx2]=sort((diag(d2))/(2*pi));
        
        % Freq2 include complex eigenfrequencies
        % Fii2 include complex eigenmodes
        % Searching only modes that corresponds positive frequency
        kk2=0; jj2=0;
        nfreq=size(Freq2,1);
        while kk2 < Request.Campbell.NroFreqs+4 % +4 is purely empirical
            jj2=jj2+1;
            % Imaginary part is equal to frequency Hz
            if imag(Freq2(jj2)) > 0.001 %
                %  i = index of rotational speed
                %  kk2 = counter for found frequency
                kk2=kk2+1; 
                FrqR(kk2,i)=imag(Freq2(jj2)); % Saving frequencies
                DampR(kk2,i)=-real(Freq2(jj2))/(abs(Freq2(jj2))); % Saving damping ratios
                %FiiR(:,kk2,i)=Fii2(nfreq/2+1:nfreq,jj2)/(Freq2(jj2)*2*pi); % Removing other half
            end 
        end
        
        clear kk2 jj2
        
    end % if full / pseudo modal
    
    
    % waitbar
%     waitbar(i/length(Request.Campbell.Spin_vec))
end % Spin_vec loop, indeksi i
% delete(h) 						% close waitbar

% Frequencies cannot be sorted according to magnitude since modes of frequencies
% may turn out crossed or new frequencies may occur
% Frequencies are tried to be sorted using curve fitting either forward direction
% or backwards (rising rotational speed, diminishing rotational speed)
% Typically frequencies changes most at low rotational speeds

FrqR2=FrqR;
DampR2=DampR;

for kk=1:Request.Campbell.NroFreqs
    
        % Curve fitting from lower to highest speed
        for ii=2:size(FrqR2,2)-1 %looppi yli nopeuksien
            % calculating prediction for frequency ii+1 using previous points
            yi = interp1(Request.Campbell.Spin_vec(1:ii)',FrqR2(kk,1:ii)',Request.Campbell.Spin_vec(1:ii+1),'pchip') ;
            ero = FrqR2(:,ii+1) - yi(ii+1)*ones(size(FrqR2(:,ii+1)));
    
            % searching index, where is the minimum difference of frequencies
            [MinEro IndMin]=min(abs(ero));
    
            % Changing elements' places
            apu=FrqR2(kk,ii+1);
            FrqR2(kk,ii+1)=FrqR2(IndMin,ii+1);
            FrqR2(IndMin,ii+1)=apu;
            %
            apu2=DampR2(kk,ii+1);
            DampR2(kk,ii+1)=DampR2(IndMin,ii+1);
            DampR2(IndMin,ii+1)=apu2;
       % end

    % Curve fitting from lower to highest speed
%     nn=size(FrqR2,2);
%     for ii=1:1:nn-2 %looppi yli nopeuksien takaperin
%         % calculating prediction for frequency ii-1 using previous points
%         yi = interp1(Request.Campbell.Spin_vec(nn-ii:nn)',FrqR2(kk,nn-ii:nn)',Request.Campbell.Spin_vec(nn-ii-1:nn),'pchip');      
%         ero = FrqR2(:,nn-ii-1) - yi(1)*ones(size(FrqR2(:,nn-ii-1)));
% 
%         % searching index, where is the minimum difference of frequencies
%         [MinEro IndMin]=min(abs(ero));        
%         
%         % Changing elements' places
%         apu=FrqR2(kk,nn-ii-1); 
%         FrqR2(kk,nn-ii-1)=FrqR2(IndMin,nn-ii-1);
%         FrqR2(IndMin,nn-ii-1)=apu;
%         % 
%         apu2=DampR2(kk,nn-ii-1); 
%         DampR2(kk,nn-ii-1)=DampR2(IndMin,nn-ii-1);
%         DampR2(IndMin,nn-ii-1)=apu2;
%                 
    end
    FrqR(kk,:)=FrqR2(kk,:);
    DampR(kk,:)=DampR2(kk,:);
    clear yi apu apu2
end

FrqR111((1+(k-1)*8):(k*8),:)=FrqR;
end
[M,N]=size(FrqR111);

Herz=(Request.Campbell.Spin_vec/(2*pi/60))/60; %speeds in Hz
vec=[];
syms x;
for j = 1:M;
        a=FrqR111(j,:);
        b=Herz;
%         figure(3)
%         plot(a,b)
    X = 0:length(a)-1;
    Y = 0:length(b)-1;
    p = polyfit(X,a,3); %USE THIS, KEPT FOR USE
    p1 = poly2sym(polyfit(X,a,1),x);
    p2 = poly2sym(polyfit(Y,b,1),x);
    root = double(real(solve(p1-p2)));
    value1=polyval(p,root(1));
    criticalspeeds(j)=value1;
end

criticalspeeds; % saving critical speeds



% Plotting------------------------------------------------------------

% % figure   % Campbell
% % if ~isfield(Request.Campbell, 'PlotDamping') || ~strcmp(Request.Campbell.PlotDamping,'No')
% %     subplot(2,1,1)
% % end
% % 
% % set(gcf,'Units','normalized')		% Units normalized (always)
% % set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
% % set(gcf,'Color',[1 1 1])			% Background color white (always)
% % 
% % % plot(Request.Campbell.Spin_vec*60/(2*pi),FrqR(1:Request.Campbell.NroFreqs,:),'k','linewidth',2)
% % plot(Request.Campbell.Spin_vec*60/(2*pi),FrqR(1:Request.Campbell.NroFreqs,:),'linewidth',2)
% % hold on
% % % Pyörimistaajuus
% % plot(Request.Campbell.Spin_vec*60/(2*pi),Request.Campbell.Spin_vec/(2*pi),'--k','linewidth',2)
% % axis([0 max(Request.Campbell.Spin_vec*60/(2*pi)) 0 1.1*max(Request.Campbell.Spin_vec/(2*pi))])
% % xlabel('RPM','fontsize',13)
% % ylabel('Frequency [Hz]','fontsize',13)
% % grid on
% % title({Inp.model_title; 'Campbell Diagram'},'fontsize',14)

if ~isfield(Request.Campbell, 'PlotDamping') || ~strcmp(Request.Campbell.PlotDamping,'No')
    %figure  % Damping ratios
    subplot(2,1,2)

    set(gcf,'Units','normalized')		% Units normalized (always)
    set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
    set(gcf,'Color',[1 1 1])			% Background color white (always)

    %plot(Request.Campbell.Spin_vec*60/(2*pi),DampR(1:Request.Campbell.NroFreqs,:)*100,'k','linewidth',2)
    plot(Request.Campbell.Spin_vec*60/(2*pi),DampR(1:Request.Campbell.NroFreqs,:)*100,'linewidth',2)
    xlabel('RPM','fontsize',13)
    ylabel('Damping ratio [%]','fontsize',13)
    damp_limit=round(max(max(DampR(1:Request.Campbell.NroFreqs,:)*100))/5)*5+5;
    axis([0 max(Request.Campbell.Spin_vec*60/(2*pi)) 0 damp_limit])
    grid on
    title('Damping Ratios','fontsize',14)
end




