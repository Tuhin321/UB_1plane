function f_ModeShapePlot(Inp,Request)
% Roottorin moodin tulostus
% INPUT:
% Inp
% Request
% Request.ModePlot.Spin_speed=5000*2*pi/60;
% Request.ModePlot.ModeNro=1; 
% Request.ModePlot.ScaleFactor=2;
%
% Written by Jussi Sopanen LTY/IMVE 2004



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

% Roottorin matriisit (n�ihin lis�t��n laakereiden nopeusriippuvaiset matriisit)
Kc=Matrx.Kc; 
Cc=Matrx.Cc;
Kc2=Matrx.Kc; % Lalannen merkinn�ill� K* matriisi (ei sis�ll� laakereiden ristitermej�)


% Calculation of Bearing Matrices
for ii=1:size(Inp.Bearing,2)
    
    
    % Calculation of bearing stiffness and damping matrix
    [Kb,Cb]=f_Bearing_Matrix(Inp.Bearing(ii),Request.ModePlot.Spin_speed);
    
    % Laakerin j�ykkyys- ja vaimennusmatriisin lis��minen systeemin matriiseihin
    % Laakerinodet m��ritetty muuttujilla INode ja JNode
    for iii=1:3
        % Muodostetaan K* matriisi, jossa laakerin ristitermit on poistettu
        % (Lalanne)
        Igdof=nnum2cgdof(Inp.Bearing(ii).Inode,iii); % akselin solmu, haetaan globaali Dof riville
        if ~ (Igdof==0) % jos ei rajoitettu
            Kc2(Igdof,Igdof)=Kc2(Igdof,Igdof) + Kb(iii,iii); %lis�t��n laakerin termit j�ykkyysmatriisiin
        end
        % tuennan solmu
        if ~Inp.Bearing(ii).Jnode==0
            Jgdof=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); % % tuennan solmu, haetaan globaali Dof riville
            if ~ (Jgdof==0) % jos ei rajoitettu
                Kc2(Jgdof,Jgdof)=Kc2(Jgdof,Jgdof) + Kb(iii,iii); %lis�t��n laakerin termit j�ykkyysmatriisiin
            end     
        end
        
        for jjj=1:3
            
            % akselin solmu  (laakeri-laakeri)
            Igdof_r=nnum2cgdof(Inp.Bearing(ii).Inode,iii); %Haetaan globaali Dof riville
            Igdof_c=nnum2cgdof(Inp.Bearing(ii).Inode,jjj); %Haetaan globaali Dof sarakkeelle
            if ~ (Igdof_r==0 || Igdof_c==0) 
                Kc(Igdof_r,Igdof_c)=Kc(Igdof_r,Igdof_c) + Kb(iii,jjj); %lis�t��n laakerin termit j�ykkyysmatriisiin
                Cc(Igdof_r,Igdof_c)=Cc(Igdof_r,Igdof_c) + Cb(iii,jjj); %lis�t��n laakerin termit vaimennusmatriisiin
            end
            
            if ~Inp.Bearing(ii).Jnode==0
                % tuennan solmu (tuenta-tuenta)
                Jgdof_r=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); %Haetaan globaali Dof riville
                Jgdof_c=nnum2cgdof(Inp.Bearing(ii).Jnode,jjj); %Haetaan globaali Dof sarakkeelle
                if ~ (Jgdof_r==0 || Jgdof_c==0)
                    Kc(Jgdof_r,Jgdof_c)=Kc(Jgdof_r,Jgdof_c) + Kb(iii,jjj); %lis�t��n laakerin termit j�ykkyysmatriisiin
                    Cc(Jgdof_r,Jgdof_c)=Cc(Jgdof_r,Jgdof_c) + Cb(iii,jjj); %lis�t��n laakerin termit vaimennusmatriisiin
                end
                
                % laakeri-tuenta ristitermit (yl�kolmio)
                if ~ (Igdof_r==0 || Jgdof_c==0)
                    Kc(Igdof_r,Jgdof_c)=Kc(Igdof_r,Jgdof_c) - Kb(iii,jjj); %lis�t��n laakerin termit j�ykkyysmatriisiin
                    Cc(Igdof_r,Jgdof_c)=Cc(Igdof_r,Jgdof_c) - Cb(iii,jjj); %lis�t��n laakerin termit vaimennusmatriisiin
                end
                % tuenta-laakeri ristitermit (alakolmio)
                if ~ (Jgdof_r==0 || Igdof_c==0)
                    Kc(Jgdof_r,Igdof_c)=Kc(Jgdof_r,Igdof_c) - Kb(iii,jjj); %lis�t��n laakerin termit j�ykkyysmatriisiin
                    Cc(Jgdof_r,Igdof_c)=Cc(Jgdof_r,Igdof_c) - Cb(iii,jjj); %lis�t��n laakerin termit vaimennusmatriisiin
                end   
            end
            
            
        end % jjj
        
    end %iii
end % ii Bearing

% Calculation of natural frequencies
% laskenta t�ysill� matriiseilla

% Gyromatriisi kerrotaan py�rimisnopeudella
Gc=Request.ModePlot.Spin_speed*Matrx.Gc;

% Eigenvalue problem A*X=lambda*X (State Space muoto)
zeroX=zeros(size(Matrx.Mc));  eyeX=eye(size(Matrx.Mc));
AA=[  zeroX                           eyeX
    -(Matrx.Mc)^(-1)*Kc      -(Matrx.Mc)^(-1)*(Cc+Gc)];

% Ominaisarvotehtavan ratkaisu 
[Fii2 d2]=eig(AA); % Gyromatriisi mukana
clear AA

% Matriisien jarjestely taajuuden mukaan.
[Freq2, idx2]=sort(diag(d2)/(2*pi));  %Freq=alfa +/- beta*i 
% Kompleksiset ominaisarvot
%for k=1:length(idx2)
%    Fiiapu(:,k)=Fii2(:,idx2(k));
%end
%Fii2=Fiiapu;  clear Fiiapu
% J�rjestell��n muodot paremmin kuin edell�
Fii2=Fii2(:,idx2); %muotojen j�rjestely

% Freq2 sis�lt�� kompleksiset ominaistaajuudet 
% Fii2 sis�lt�� kompleksiset ominaismuodot
% Etsit��n vain positiivista ominaistaajuutta vastaavat muodot 
kk2=0; 
nfreq=size(Freq2,1);
for jj2=1:nfreq
    % Imagin��riosa on yht� kuin taajuus Hz
    if imag(Freq2(jj2)) > 0.001 %
        %  kk2 = l�ydetyn taajuuden laskuri
        kk2=kk2+1; 
        FrqR(kk2)=imag(Freq2(jj2)); % Tallennetaan taajuudet
        DampR(kk2)=-real(Freq2(jj2))/(abs(Freq2(jj2))); % Tallennetaan vaimennussuhteet
        FiiR(:,kk2)=Fii2(nfreq/2+1:nfreq,jj2)/(Freq2(jj2)*2*pi); % Pudotetaan toinen puolikas pois
    end 
end
clear kk2 jj2

% Muotojen tulostus--------------------------------------------------------

% kopioidaan matriisit
Node=Inp.Node;
Elem=Inp.Elem;

% Piirret��n alkuper�inen rakenne ilman laakereita ja rajoitteita. My�s
% solmu ja elementtinumerot j�tet��n pois
if isfield(Request.ModePlot, 'PlotGeom') && strcmp(Request.ModePlot.PlotGeom,'Yes')
    %asetetaan default optiot
    origloc=[0 0 0]; %origon asema
    shade=1;   %shaded plot
    alpha=0.1; %pinnat osittain l�pin�kyvi�

    % tarkistetaan onko optiot m��ritelty
    if isfield(Request.ModePlot,'Options')
        for i=1:size(Request.ModePlot.Options,1)
            switch Request.ModePlot.Options{i,1}
                case 'originloc'
                    origloc=Request.ModePlot.Options{i,2}; %origon asema
                case 'shading'
                    switch Request.ModePlot.Options{i,2};
                        case 'off'; shade=0;
                        case 'Off'; shade=0;
                        case 'OFF'; shade=0;    
                    end
                case 'transparency'
                    alpha=Request.ModePlot.Options{i,2};
            end
        end
    end

    % Draw Model as a 3D WireFrame plot % Funktio avaa uuden figuren
    if shade==0; f_RGeomPlotWireFr(Inp, [1 1 0 0 1], [0 0],{'-k',1}); end

    % Draw Model as a Shaded 3D plot
    if shade==1; f_RGeomPlotShaded(Inp,alpha); end
       
else % jos ei plotata roottorin geometriaa, niin avataan uusi kuva
    figure  
    set(gcf,'Units','normalized')		% Units normalized (always)
    set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
    set(gcf,'Color',[1 1 1])			% Background color white (always)
end

mm=Request.ModePlot.ModeNro; % Tulostettavan muodon numero

t=(-1.8*pi:pi/20:0); % aika tai oikeammin omega*t
% sf=5000; %scale factor
% muodon tulostus eli piirto
for ii=1:size(Node,1)
    %nnum=Node(ii,1); %solmunumero
    %solmunumeroa vastaava globaali DOF    
    UXindx=nnum2cgdof(Node(ii,1),1);
    UYindx=nnum2cgdof(Node(ii,1),2);
    UZindx=nnum2cgdof(Node(ii,1),3);
    
    % Noden x-siirtym�
    if ~UXindx==0 % jos ei rajoitettu gdof
        u1(ii,:)=Node(ii,2)*ones(1,length(t))+...
                 Request.ModePlot.ScaleFactor*(real(FiiR(UXindx,mm))*cos(t) - imag(FiiR(UXindx,mm))*sin(t));
    else
        u1(ii,:)=Node(ii,2)*ones(1,length(t));
    end
    % Noden y-siirtym�
    if ~UYindx==0 % jos ei rajoitettu gdof
        v1(ii,:)=Node(ii,3)*ones(1,length(t))+...
                 Request.ModePlot.ScaleFactor*(real(FiiR(UYindx,mm))*cos(t) - imag(FiiR(UYindx,mm))*sin(t));    
    else
        v1(ii,:)=Node(ii,3)*ones(1,length(t));
    end
    % Noden z-siirtym�
    if ~UZindx==0 % jos ei rajoitettu gdof
        w1(ii,:)=Node(ii,4)*ones(1,length(t))+...
                 Request.ModePlot.ScaleFactor*(real(FiiR(UZindx,mm))*cos(t) - imag(FiiR(UZindx,mm))*sin(t));        
    else
        w1(ii,:)=Node(ii,4)*ones(1,length(t));
    end

end %node

% Transpoosit siirtym�matriiseista plot3 k�sky� varten
u1=u1'; v1=v1'; w1=w1';

% Tulostetaan ympyr�t
H=plot3(u1, v1, w1,'b','linewidth',1);
h_leg(1)=H(1); %tallennetaan handle legendiin
t_leg{1}='Whirl (\Omegat=-1.8\pi...0)';
hold on
% pallot solmuihin
% mu=length(u1); 
mu=size(u1,1);
h_leg(2)=plot3(u1(mu,:), v1(mu,:), w1(mu,:),'ok','markersize',6,'linewidth',2);
t_leg{2}='Node';

% Tulostetaan palkkielementit viivoina ja massaelementti pisteen�
t=0;
for ii=1:size(Elem,1)

    if ~Elem(ii,3) == 0  % Jos ei ole massaelementti
        %for jj=1:size(Node,1)
        %    %Haetaan I noden rivi Node matriisista
        %    if Elem(ii,2)==Node(jj,1)
        %        Iindx=jj;
        %    end
        %    %Haetaan I noden rivi Node matriisista
        %    if Elem(ii,3)==Node(jj,1)
        %        Jindx=jj;
        %    end
        %end
        % n�m� rivit tekev�t saman asian kuin edellinen for-if h�ss�kk�
        Iindx= (Elem(ii,2)==Node(:,1));
        Jindx= (Elem(ii,3)==Node(:,1));

        %solmunumeroa vastaava globaali DOF
        UXindxI=nnum2cgdof(Node(Iindx,1),1);
        UYindxI=nnum2cgdof(Node(Iindx,1),2);
        UZindxI=nnum2cgdof(Node(Iindx,1),3);
        UXindxJ=nnum2cgdof(Node(Jindx,1),1);
        UYindxJ=nnum2cgdof(Node(Jindx,1),2);
        UZindxJ=nnum2cgdof(Node(Jindx,1),3);

        % Noden x-siirtym�
        if ~UXindxI==0 % jos ei rajoitettu gdof
            uI=Node(Iindx,2)+Request.ModePlot.ScaleFactor*real(FiiR(UXindxI,mm))*cos(t) ;
        else
            uI=Node(Iindx,2);
        end
        if ~UXindxJ==0 % jos ei rajoitettu gdof
            uJ=Node(Jindx,2)+Request.ModePlot.ScaleFactor*real(FiiR(UXindxJ,mm))*cos(t) ;
        else
            uJ=Node(Jindx,2);
        end

        % Noden y-siirtym�
        if ~UYindxI==0 % jos ei rajoitettu gdof
            vI=Node(Iindx,3)+Request.ModePlot.ScaleFactor*real(FiiR(UYindxI,mm))*cos(t) ;
        else
            vI=Node(Iindx,3);
        end
        if ~UYindxJ==0 % jos ei rajoitettu gdof
            vJ=Node(Jindx,3)+Request.ModePlot.ScaleFactor*real(FiiR(UYindxJ,mm))*cos(t) ;
        else
            vJ=Node(Jindx,3);
        end
        % Noden z-siirtym�
        if ~UZindxI==0 % jos ei rajoitettu gdof
            wI=Node(Iindx,4)+Request.ModePlot.ScaleFactor*real(FiiR(UZindxI,mm))*cos(t) ;
        else
            wI=Node(Iindx,4);
        end
        if ~UZindxJ==0 % jos ei rajoitettu gdof
            wJ=Node(Jindx,4)+Request.ModePlot.ScaleFactor*real(FiiR(UZindxJ,mm))*cos(t) ;
        else
            wJ=Node(Jindx,4);
        end


        %Elementti viivana
        plot3( [Node(Iindx,2) Node(Jindx,2)], [Node(Iindx,3) Node(Jindx,3)],...
            [Node(Iindx,4) Node(Jindx,4)],'-k','linewidth',1) % keskiviva
        h_leg(3)=plot3( [uI uJ], [vI vJ], [wI wJ ],'-k','linewidth',2); % viiva solmusta solmuun
        t_leg{3}='Centerline';
        if ii==1 %viivat keskelt� solmuun
            plot3( [Node(Iindx,2) uI ], [Node(Iindx,3) vI], [Node(Iindx,4) wI],'-k','linewidth',2)
            plot3( [Node(Jindx,2) uJ ], [Node(Jindx,3) vJ], [Node(Jindx,4) wJ],'-k','linewidth',2)
        else
            plot3( [Node(Jindx,2) uJ ], [Node(Jindx,3) vJ], [Node(Jindx,4) wJ],'-k','linewidth',2)
        end

    else % jos on massaelementti

        %for jj=1:size(Node,1)
        %    %Haetaan I noden rivi Node matriisista
        %    if Elem(ii,2)==Node(jj,1)
        %       Iindx=jj;
        %    end
        %end
        % vastaava kuin edell�, mutta nopeampi
        Iindx=(Elem(ii,2)==Node(:,1));

        %solmunumeroa vastaava globaali DOF
        UXindxI=nnum2cgdof(Node(Iindx,1),1);
        UYindxI=nnum2cgdof(Node(Iindx,1),2);
        UZindxI=nnum2cgdof(Node(Iindx,1),3);

        % Noden x-siirtym�
        if ~UXindxI==0 % jos ei rajoitettu gdof
            uI=Node(Iindx,2)+Request.ModePlot.ScaleFactor*real(FiiR(UXindxI,mm))*cos(t) ;
        else
            uI=Node(Iindx,2);
        end

        % Noden y-siirtym�
        if ~UYindxI==0 % jos ei rajoitettu gdof
            vI=Node(Iindx,3)+Request.ModePlot.ScaleFactor*real(FiiR(UYindxI,mm))*cos(t) ;
        else
            vI=Node(Iindx,3);
        end

        % Noden z-siirtym�
        if ~UZindxI==0 % jos ei rajoitettu gdof
            wI=Node(Iindx,4)+Request.ModePlot.ScaleFactor*real(FiiR(UZindxI,mm))*cos(t) ;
        else
            wI=Node(Iindx,4);
        end

        %Massaelementttipisteen�
        h_leg(4)=plot3( uI, vI, wI ,'or','linewidth',3,'markersize',10);
        t_leg{4}='MassPoint';

    end %massaelementti
end % elem

NumLeg=size(h_leg,2)+1; %tarkistetaan onko legendin arvoja 3 vai 4

% Laakerin plottaus
for kkk=1:size(Inp.Bearing,2)

    for jj=1:size(Node,1)
        %Haetaan I noden rivi Node matriisista
        if Inp.Bearing(kkk).Inode==Node(jj,1)
            Iindx=jj;
        end
    end

    %solmunumeroa vastaava globaali DOF
    UXindxI=nnum2cgdof(Node(Iindx,1),1);
    UYindxI=nnum2cgdof(Node(Iindx,1),2);
    UZindxI=nnum2cgdof(Node(Iindx,1),3);
    
    % Noden x-siirtym�
    if ~UXindxI==0 % jos ei rajoitettu gdof
        uI=Node(Iindx,2)+Request.ModePlot.ScaleFactor*real(FiiR(UXindxI,mm))*cos(t) ;
    else
        uI=Node(Iindx,2);
    end

    % Noden y-siirtym�
    if ~UYindxI==0 % jos ei rajoitettu gdof
        vI=Node(Iindx,3)+Request.ModePlot.ScaleFactor*real(FiiR(UYindxI,mm))*cos(t) ;
    else
        vI=Node(Iindx,3);
    end

    % Noden z-siirtym�
    if ~UZindxI==0 % jos ei rajoitettu gdof
        wI=Node(Iindx,4)+Request.ModePlot.ScaleFactor*real(FiiR(UZindxI,mm))*cos(t) ;
    else
        wI=Node(Iindx,4);
    end

    %Piirret��n laakeri vihre�n� neli�n�
    h_leg(NumLeg)=plot3( uI, vI, wI,'sg','linewidth',3,'markersize',12);
    t_leg{NumLeg}='Bearing';
end
    
% BW/FW muodon tunnistus---------------------------------------------------
% fiidot=zeros(1,size(Node,1));
% sfiidot=zeros(1,size(Node,1)); %alustus
k=0;
for i=1:size(Node,1)
    % haetaan vapausasteet
    UYindx=nnum2cgdof(Node(i,1),2);
    UZindx=nnum2cgdof(Node(i,1),3);

    if ~UZindx==0 && ~UYindx==0 % jos ei rajoitettu gdof
        k=k+1;
        % % FW/BW laskenta by Chen&Gunter s. 44
        zc=real(FiiR(UZindx,mm));
        yc=real(FiiR(UYindx,mm));
        zs=-imag(FiiR(UZindx,mm)); % kokeilu! Oikea muoto onkin konjugaatti
        ys=-imag(FiiR(UYindx,mm));

        % kiertosuunta eri kuin Chen&Gunter
        fiidot(k)=-(yc*zs-ys*zc);
        if abs(fiidot(k)) < 1e-14; fiidot(k)=0; end; %py�ristet��n pienet luvut nollaan
        sfiidot(k)=sign(fiidot(k)); %kieppumisen suuntakulman nopeus
    end
end

% lasketaan solmujen arvoista yksi luku (=keskiarvo)
sfiidot=sum(sfiidot)/length(sfiidot);

% asetetaan kuvaa varten stringi
if sfiidot == 1
    modetype='BW';
elseif sfiidot == -1
    modetype='FW';
elseif sfiidot == 0
    modetype='Line Motion';
else 
    modetype='Mixed BW/FW';
end
%test+++++++++++++++
% if sfiidot > 0.95 
%     modetype='BW';
% elseif sfiidot < -0.95
%     modetype='FW';
% elseif sfiidot < -0.05 && sfiidot > 0.05
%     modetype='Line Motion';
% else 
%     modetype='Mixed BW/FW';
% end
% test end++++++++++++++
%--------------------------------------------------------------------------

% Kuvan otsikon asetus
str0=['Mode ' num2str(mm) ' @ ' num2str(Request.ModePlot.Spin_speed*60/(2*pi)) ' rpm (' modetype '), Frequency = ' num2str(FrqR(mm),4) ' Hz, Damping Ratio = ' num2str(DampR(mm)*100,3) ' %'];
str2={ Inp.model_title; str0};
title(str2,'Fontsize',14)
xlabel('X','fontsize',13)
ylabel('Y','fontsize',13)
zlabel('Z','fontsize',13)

% legend
% legend(h_leg,t_leg,'Location', 'SouthOutSide');

% kuvan asetukset
% Turn warning temporaly off
warning off MATLAB:hg:patch:RGBColorDataNotSupported
view(3);  axis equal; axis vis3d; axis off;
warning on MATLAB:hg:patch:RGBColorDataNotSupported
% view3d rot; 

% koordinaatiston piirto
%f_draw_caxis(origloc, sz)
f_draw_caxis;