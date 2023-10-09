function [MaxAmp,U_vector]=f_UBResponse(Inp,Request)

% funktio laskee ep‰tasapainon aiheuttaman steady state vasteen.
%         
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
h = waitbar(0,'Computing Unbalance Response...');

% aika parametri, h‰nt‰ l‰hes t‰ysi kierros
t=(-1.8*pi:2*pi/20:0); 
% alustetaan tulosmatriisit
MaxAmp=zeros(length(Request.UBResp.OutputNodes),length(Request.UBResp.Spin_vec));
U_vector=zeros(size(Matrx.Mc,1),size(t,2),length(Request.UBResp.Spin_vec));

% Start of speed loop -----------------------------------------------------
for i=1:length(Request.UBResp.Spin_vec)
    
    % Roottorin matriisit (t‰ytyy m‰‰ritell‰ uudelleen joka kierroksella,
    % koska n‰ihin lis‰t‰‰n laakereiden nopeusriippuvaiset matriisit)
    Kc=Matrx.Kc; 
    Cc=Matrx.Cc;
    Kc2=Matrx.Kc; % Lalannen merkinnˆill‰ K* matriisi (ei sis‰ll‰ laakereiden ristitermej‰)
    
    
    % Calculation of Bearing Matrices
    for ii=1:size(Inp.Bearing,2)
        
        
        % Calculation of bearing stiffness and damping matrix
        [Kb,Cb]=f_Bearing_Matrix(Inp.Bearing(ii),Request.UBResp.Spin_vec(i));
        
        % Laakerin j‰ykkyys- ja vaimennusmatriisin lis‰‰minen systeemin matriiseihin
        % Laakerinodet m‰‰ritetty muuttujilla INode ja JNode
        for iii=1:3
            % Muodostetaan K* matriisi, jossa laakerin ristitermit on poistettu
            % (Lalanne)
            Igdof=nnum2cgdof(Inp.Bearing(ii).Inode,iii); % akselin solmu, haetaan globaali Dof riville
            if ~ (Igdof==0) % jos ei rajoitettu
                Kc2(Igdof,Igdof)=Kc2(Igdof,Igdof) + Kb(iii,iii); %lis‰t‰‰n laakerin termit j‰ykkyysmatriisiin
            end
            % tuennan solmu
            if ~Inp.Bearing(ii).Jnode==0
                Jgdof=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); % % tuennan solmu, haetaan globaali Dof riville
                if ~ (Jgdof==0) % jos ei rajoitettu
                    Kc2(Jgdof,Jgdof)=Kc2(Jgdof,Jgdof) + Kb(iii,iii); %lis‰t‰‰n laakerin termit j‰ykkyysmatriisiin
                end     
            end
            
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


    % Calculation of steady state responce++++++++++++++++
    % laskenta t‰ysill‰ matriiseilla

    % pyˆrimisnopeus
    Spin=Request.UBResp.Spin_vec(i);

    % Gyromatriisi kerrotaan pyˆrimisnopeudella
    Gc=Spin*Matrx.Gc;

    % Steady State tila Lalanne s. 102 kaikilla vapausasteilla
    AA=[Kc-Matrx.Mc*Spin^2     -Spin*(Cc+Gc)
    Spin*(Cc+Gc)           Kc-Matrx.Mc*Spin^2];
    
    % Voimavektori    
    fu=[Spin^2*Matrx.F2
        Spin^2*Matrx.F3];

    % ratkaistaan Steady State tila
    pp=AA\fu;


    % Muutetaan koordinaatit fyysisiksi siirtymisksi ja lasketaan
    % maksimi amplitudi
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
        MaxAmp(kkk,i)=max(sqrt(U_vec_t(dofy,:).^2+U_vec_t(dofz,:).^2));
        %Request.UBResp.OutputNodes=[2, 6, 13, 16];
    end
    
    % waitbar
    waitbar(i/length(Request.UBResp.Spin_vec))
    
end % Spin vec loop, indeksi i
    
delete(h) 						% close waitbar

% Tulostus
figure 
set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)


plot(Request.UBResp.Spin_vec*60/(2*pi),MaxAmp*1e6,'linewidth',2) %) %,'b-'

% v=axis; v(3)=1e-8; axis(v); clear v
xlabel('RPM','fontsize',13)
ylabel('Maximum Amplitude ( \mum)', 'fontsize',13)
grid on

% luodaan legendille string muuttuja
for kkk=1:size(Request.UBResp.OutputNodes,2)
    str{kkk}=['Node ' num2str(Request.UBResp.OutputNodes(kkk)) ' Mag'];
end
legend(str)
title({Inp.model_title; 'Unbalance Response'},'fontsize',14)



warning on MATLAB:nearlySingularMatrix

