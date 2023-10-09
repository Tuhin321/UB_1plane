function [FrqR,DampR,Matrx,cr_area]=f_Campbell(Inp,Request)
% Funktio laskee ominaistaajuudet pyˆrimisnopeuden funktiona
% Optiot: - full, jossa k‰ytet‰‰n t‰ysi‰ matriisiseja Kc, Mc Cc, Gc
%         - pseudomodal, jossa lasketaan n kpl muotoja ja tehd‰‰n laskenta niill‰
%         
%         Vaatii Request-structuren
%         Request.Campbell.NroFreqs=6;
%         Request.Campbell.Spin_vec=[0:1000:10000]*2*pi/60;
%         Request.Campbell.Type='Full'; % tai
%         Request.Campbell.Type='PseudoModal';
%         Request.Campbell.NroModes=30;
% 
% Input:   Inp
%          Request
%                        
% Palauttaa matriiseihin taajuudet, rpm ja vaimennukset eli sama kuin plotataan

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

% alustetaan nopeus ja vaimennus matriisit
FrqR=zeros(Request.Campbell.NroFreqs+4,length(Request.Campbell.Spin_vec));
DampR=zeros(Request.Campbell.NroFreqs+4,length(Request.Campbell.Spin_vec));

% h = waitbar(0,'Computing Campbell diagram...');

% Start of speed loop
for i=1:length(Request.Campbell.Spin_vec)
    
    % Roottorin matriisit (t‰ytyy m‰‰ritell‰ uudelleen joka kierroksella,
    % koska n‰ihin lis‰t‰‰n laakereiden nopeusriippuvaiset matriisit)
    Kc=Matrx.Kc; 
    Cc=Matrx.Cc;
    Kc2=Matrx.Kc; % Lalannen merkinnˆill‰ K* matriisi (ei sis‰ll‰ laakereiden ristitermej‰)
    
    
    % Calculation of Bearing Matrices
    for ii=1:size(Inp.Bearing,2)
        
        
        % Calculation of bearing stiffness and damping matrix
        [Kb,Cb]=f_Bearing_Matrix(Inp.Bearing(ii),Request.Campbell.Spin_vec(i));
        
        if ii==1
        Kbear(1,i)=Kb(2,2); %kyy
        Kbear(2,i)=Kb(3,3); %kzz
        Cbear(1,i)=Cb(2,2); %kyy
        Cbear(2,i)=Cb(3,3); %kzz
        elseif ii==2
        Kbear(3,i)=Kb(2,2); %kyy
        Kbear(4,i)=Kb(3,3); %kzz
        Cbear(3,i)=Cb(2,2); %kyy
        Cbear(4,i)=Cb(3,3); %kzz
        end
        
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
    
    % Calculation of natural frequencies
    if strcmp(Request.Campbell.Type,'Full')
        % laskenta t‰ysill‰ matriiseilla
        
        % Gyromatriisi kerrotaan pyˆrimisnopeudella
        Gc=Request.Campbell.Spin_vec(i)*Matrx.Gc;
        
        % Eigenvalue problem A*X=lambda*X (State Space muoto)
        zeroX=zeros(size(Matrx.Mc));  eyeX=eye(size(Matrx.Mc));
        AA=[  zeroX                           eyeX
            -(Matrx.Mc)^(-1)*Kc      -(Matrx.Mc)^(-1)*(Cc+Gc)];
        
        % Ominaisarvotehtavan ratkaisu 
        [Fii2 d2]=eig(AA); % Gyromatriisi mukana
        clear AA
        
        % Matriisien jarjestely taajuuden mukaan.
        [Freq2, idx2]=sort(diag(d2)/(2*pi));  %Freq=alfa +/- beta*i 
        
        % Freq2 sis‰lt‰‰ kompleksiset ominaistaajuudet 
        % Etsit‰‰n vain positiivista ominaistaajuutta vastaavat muodot 
        kk2=0; jj2=0;
        nfreq=size(Freq2,1);
        while kk2 < Request.Campbell.NroFreqs+2  % +2 on puhtaasti empiirinen 
            jj2=jj2+1;
            % Imagin‰‰riosa on yht‰ kuin taajuus Hz
            if imag(Freq2(jj2)) > 0.001 %
                %  i = pyˆrimisnopeuden indeksi
                %  kk2 = lˆydetyn taajuuden laskuri
                kk2=kk2+1; 
                FrqR(kk2,i)=imag(Freq2(jj2)); % Tallennetaan taajuudet
                DampR(kk2,i)=-real(Freq2(jj2))/(abs(Freq2(jj2))); % Tallennetaan vaimennussuhteet
                %FiiR(:,kk2,i)=Fii2(nfreq/2+1:nfreq,jj2)/(Freq2(jj2)*2*pi); % Pudotetaan toinen puolikas pois
            end 
        end
        clear kk2 jj2
        
    else
        % Laskenta speudo modal menetelm‰ll‰
        
        % Ominaisarvotehtavan ratkaisu
        [FiiN d]=eig(Kc2,Matrx.Mc);
        % Matriisien jarjestely taajuuden mukaan.
        [Freq, idx]=sort(sqrt(diag(d))/(2*pi));
        
        %for k=1:length(idx)
        %   Fiiapu(:,k)=FiiN(:,idx(k));
        %end
        %FiiN=Fiiapu;
        %clear Fiiapu idx k d        
        % tekee saman kuin edell‰, eli j‰rjest‰‰ muodot taajuuden mukaan
        FiiN=FiiN(:,idx); clear idx d;
                
        % Muotojen maaran valinta (Nmuoto alinta ominaismuotoa)
        FiiNn=FiiN(:,1:Request.Campbell.NroModes);
        
        % Massanormeeraus
        for k=1:Request.Campbell.NroModes
            FiiNormN(:,k)=FiiNn(:,k)/sqrt(FiiNn(:,k)'*Matrx.Mc*FiiNn(:,k));
        end
        
        % Moodaaliset matriisit
        Mm=FiiNormN'*Matrx.Mc*FiiNormN;
        Km=FiiNormN'*Kc2*FiiNormN + FiiNormN'*(Kc-Kc2)*FiiNormN; % j‰lkimm‰inen termi sis‰lt‰‰ laakerin ristitermit
        
        % Gyromatriisi kerrotaan pyˆrimisnopeudella
        Gc=Request.Campbell.Spin_vec(i)*Matrx.Gc;
        % vaimennusmatriisi + gyromatriisi
        Cm=FiiNormN'*(Cc+Gc)*FiiNormN;
        
        % Eigenvalue problem A*X=lambda*X
        zeroX=zeros(size(Mm));  eyeX=eye(size(Mm));
        AA=[  zeroX               eyeX
            -(Mm)^(-1)*Km      -(Mm)^(-1)*Cm];
        
        [Fii2 d2]=eig(AA); % State space ominaisarvoteht‰v‰n ratkaisu
        clear AA
        
        % Matriisien jarjestely taajuuden mukaan.
        [Freq2, idx2]=sort((diag(d2))/(2*pi));
        
        % Freq2 sis‰lt‰‰ kompleksiset ominaistaajuudet 
        % Fii2 sis‰lt‰‰ kompleksiset ominaismuodot
        % Etsit‰‰n vain positiivista ominaistaajuutta vastaavat muodot 
        kk2=0; jj2=0;
        nfreq=size(Freq2,1);
        while kk2 < Request.Campbell.NroFreqs+4 % +4 on puhtaasti empiirinen 
            jj2=jj2+1;
            % Imagin‰‰riosa on yht‰ kuin taajuus Hz
            if imag(Freq2(jj2)) > 0.001 %
                %  i = pyˆrimisnopeuden indeksi
                %  kk2 = lˆydetyn taajuuden laskuri
                kk2=kk2+1; 
                FrqR(kk2,i)=imag(Freq2(jj2)); % Tallennetaan taajuudet
                DampR(kk2,i)=-real(Freq2(jj2))/(abs(Freq2(jj2))); % Tallennetaan vaimennussuhteet
                %FiiR(:,kk2,i)=Fii2(nfreq/2+1:nfreq,jj2)/(Freq2(jj2)*2*pi); % Pudotetaan toinen puolikas pois
            end 
        end
        
        clear kk2 jj2
        
    end % if full / pseudo modal
    
    
    % waitbar
    waitbar(i/length(Request.Campbell.Spin_vec))
    
    
end % Spin_vec loop, indeksi i

% delete(h) 						% close waitbar


% Taajuuksia ei voi j‰rjest‰‰ suoraan suuruuden mukaan, koska muotojen
% taajuudet voivat menn‰ ristiin tai ratkaisussa tulee esille uusia
% taajuuksia
% Taajuudet pyrit‰‰n j‰rjest‰m‰‰n k‰yr‰n sovituksen avulla joko etuperin
% tai takaperin (kasvava pyˆrimisnopeus, pienenev‰ pyˆrimisnopeus)
% Yleens‰ taajuudet muuttuvat eniten pienill‰ pyˆrimisnopeuksilla
FrqR2=FrqR;
DampR2=DampR;
for kk=1:Request.Campbell.NroFreqs
    


    %     % K‰yr‰n sovitus pienimm‰st‰ suurimpaan nopeuteen
    %     for ii=2:size(FrqR2,2)-1 %looppi yli nopeuksien
    %         % lasketaan ennuste taajuudelle ii+1 edellisten pisteiden avulla
    %         yi = interp1(Request.Campbell.Spin_vec(1:ii)',FrqR2(kk,1:ii)',Request.Campbell.Spin_vec(1:ii+1),'pchip') ;
    %         ero = FrqR2(:,ii+1) - yi(ii+1)*ones(size(FrqR2(:,ii+1)));
    %
    %         %etsit‰‰n indeksi, jossa minimi taajuuksien ero
    %         [MinEro IndMin]=min(abs(ero));
    %
    %         % Vaihdetaan alkioiden paikkoja
    %         apu=FrqR2(kk,ii+1);
    %         FrqR2(kk,ii+1)=FrqR2(IndMin,ii+1);
    %         FrqR2(IndMin,ii+1)=apu;
    %         %
    %         apu2=DampR2(kk,ii+1);
    %         DampR2(kk,ii+1)=DampR2(IndMin,ii+1);
    %         DampR2(IndMin,ii+1)=apu2;
    %     end

    % K‰yr‰n sovitus suurimmassa nopeudesta alasp‰in
    nn=size(FrqR2,2);
    for ii=1:1:nn-2 %looppi yli nopeuksien takaperin
        % lasketaan ennuste taajuudelle ii-1 edellisten pisteiden avulla
        yi = interp1(Request.Campbell.Spin_vec(nn-ii:nn)',FrqR2(kk,nn-ii:nn)',Request.Campbell.Spin_vec(nn-ii-1:nn),'pchip');      
        ero = FrqR2(:,nn-ii-1) - yi(1)*ones(size(FrqR2(:,nn-ii-1)));

        %etsit‰‰n indeksi, jossa minimi taajuuksien ero
        [MinEro IndMin]=min(abs(ero));        
        
        % Vaihdetaan alkioiden paikkoja
        apu=FrqR2(kk,nn-ii-1); 
        FrqR2(kk,nn-ii-1)=FrqR2(IndMin,nn-ii-1);
        FrqR2(IndMin,nn-ii-1)=apu;
        % 
        apu2=DampR2(kk,nn-ii-1); 
        DampR2(kk,nn-ii-1)=DampR2(IndMin,nn-ii-1);
        DampR2(IndMin,nn-ii-1)=apu2;
                
    end
    FrqR(kk,:)=FrqR2(kk,:);
    DampR(kk,:)=DampR2(kk,:);
    clear yi apu apu2
end

% Operation area

n_cr_area=1;
cr_area=[];
for cr_r=Request.Campbell.plotset(1):Request.Campbell.plotset(1):Request.Campbell.NroFreqs
    for nX=Request.Campbell.nX(1:1:end)
        value=find(Request.Campbell.Spin_vec/(2*pi)*nX>FrqR(cr_r,:),1);
        if value>0
        cr_area(end+1,:)=[cr_r; nX; FrqR(cr_r,value); Request.Campbell.Spin_vec(value)*(1-Request.Campbell.margin);Request.Campbell.Spin_vec(value); Request.Campbell.Spin_vec(value)*(1+Request.Campbell.margin)];  
        end
    end
end


% Tulostus------------------------------------------------------------

figure   % Campbell
if ~isfield(Request.Campbell, 'PlotDamping') || ~strcmp(Request.Campbell.PlotDamping,'No')
    subplot(2,1,1)
end

set(gcf,'Units','normalized')		% Units normalized (always)
set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
set(gcf,'Color',[1 1 1])			% Background color white (always)

% plot(Request.Campbell.Spin_vec*60/(2*pi),FrqR(1:Request.Campbell.NroFreqs,:),'k','linewidth',2)
plot(Request.Campbell.Spin_vec*60/(2*pi),FrqR(1:Request.Campbell.NroFreqs,:),'linewidth',2)
hold on
% Pyˆrimistaajuus
for nX=Request.Campbell.nX(1:1:end)
    plot(Request.Campbell.Spin_vec*60/(2*pi),Request.Campbell.Spin_vec/(2*pi)*nX,'--k','linewidth',2)
end
if Request.Campbell.plotset(2)>0;  
    for ver_lin=1:1:length(cr_area(:,1))
        plot(ones(1,2).*cr_area(ver_lin,5)*60/(2*pi), [0 cr_area(ver_lin,3)],'k','linewidth',2)
    end
end
if Request.Campbell.plotset(3)>0;  
    for ver_lin2=1:1:length(cr_area(:,1))
        plot(ones(1,2).*cr_area(ver_lin2,4)*60/(2*pi), [0 cr_area(ver_lin2,3)],'r','linewidth',2)
        plot(ones(1,2).*cr_area(ver_lin2,6)*60/(2*pi), [0 cr_area(ver_lin2,3)],'r','linewidth',2)
    end
end
    
axis([0 max(Request.Campbell.Spin_vec*60/(2*pi)) 0 max(Request.Campbell.Spin_vec/(2*pi))])
xlabel('RPM','fontsize',13)
ylabel('Frequency [Hz]','fontsize',13)
grid on
title({Inp.model_title; 'Campbell Diagram'},'fontsize',14)

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

% plots the interpolated bearing stiffness parameters as functions of RPM
% in case of user defined speed dependent bearing stiffness coefficients
% 5.7.2017 // JEH
if ~isfield(Inp.Bearing(1).type, 'User Defined Speed Dependent B') 
figure
subplot(2,1,1)
plot(Request.Campbell.Spin_vec*60/(2*pi),Kbear(1,:),'k')
hold on
plot(Request.Campbell.Spin_vec*60/(2*pi),Kbear(2,:),'r')
plot(Request.Campbell.Spin_vec*60/(2*pi),Kbear(3,:),'b')
plot(Request.Campbell.Spin_vec*60/(2*pi),Kbear(4,:),'g')
legend('Drive - kyy','Drive - kzz','Service - kyy','Service - kzz')
xlabel('RPM','fontsize',13)
ylabel('Bearing stiffness','fontsize',13)
grid on
title('Interpolated bearing stiffness','fontsize',14)

subplot(2,1,2)
plot(Request.Campbell.Spin_vec*60/(2*pi),Cbear(1,:),'k')
hold on
plot(Request.Campbell.Spin_vec*60/(2*pi),Cbear(2,:),'r')
plot(Request.Campbell.Spin_vec*60/(2*pi),Cbear(3,:),'b')
plot(Request.Campbell.Spin_vec*60/(2*pi),Cbear(4,:),'g')
legend('Drive - cyy','Drive - czz','Service - cyy','Service - czz')
xlabel('RPM','fontsize',13)
ylabel('Bearing damping','fontsize',13)
grid on
title('Interpolated bearing damping','fontsize',14)
end
