function [FrqR,DampR]=f_Instability(Inp,Request)

% Funktio laskee epästabiilius rajan
% Optiot: - full, jossa käytetään täysiä matriisiseja Kc, Mc Cc, Gc
%         - pseudomodal, jossa lasketaan n kpl muotoja ja tehdään laskenta niillä
%         
%         --> Lähtötiedostoon Request-structure
%         Request.Instability.Spin_Range=[0 10000]*2*pi/60;
%         Request.Instability.Speed_Toler=1*2*pi/60;
%         Request.Instability.Type='Full'; % tai
%         Request.Instability.Type='PseudoModal';
%         Request.Instability.NroModes=30;
% 
%          
%         
% 
% Input:   Inp
%          Request
%          
% 
%                
% Palauttaa matriiseihin taajuudet, rpm ja vaimennukset eli samta kuin plotataan

% 7.1.2013 havaittu bugi ja korjattu muodon iteroinnissa J. Sopanen

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

FrqR=zeros(8,2);
DampR=zeros(8,2);
% FrqR=zeros(Request.Campbell.NroFreqs+4,length(Request.Campbell.Spin_vec));
% DampR=zeros(Request.Campbell.NroFreqs+4,length(Request.Campbell.Spin_vec));

h = waitbar(0,'Iterating threshold speed of instability...');

Spin(1)=Request.Instability.Spin_Range(1);
Spin(2)=mean(Request.Instability.Spin_Range);
Spin(3)=Request.Instability.Spin_Range(2);
lask=0;
while (Spin(3)-Spin(1) > Request.Instability.Speed_Toler) & lask < 3;
    
    % Start of speed loop
    for i=1:3% length(Request.Campbell.Spin_vec)
        
        % Roottorin matriisit (täytyy määritellä uudelleen joka kierroksella,
        % koska näihin lisätään laakereiden nopeusriippuvaiset matriisit)
        Kc=Matrx.Kc; 
        Cc=Matrx.Cc;
        Kc2=Matrx.Kc; % Lalannen merkinnöillä K* matriisi (ei sisällä laakereiden ristitermejä)
        
        
        % Calculation of Bearing Matrices
        for ii=1:size(Inp.Bearing,2)
            
            
            % Calculation of bearing stiffness and damping matrix
            [Kb,Cb]=f_Bearing_Matrix(Inp.Bearing(ii),Spin(i));
            
            % Laakerin jäykkyys- ja vaimennusmatriisin lisääminen systeemin matriiseihin
            % Laakerinodet määritetty muuttujilla INode ja JNode
            for iii=1:3
                % Muodostetaan K* matriisi, jossa laakerin ristitermit on poistettu
                % (Lalanne)
                Igdof=nnum2cgdof(Inp.Bearing(ii).Inode,iii); % akselin solmu, haetaan globaali Dof riville
                if ~ (Igdof==0) % jos ei rajoitettu
                    Kc2(Igdof,Igdof)=Kc2(Igdof,Igdof) + Kb(iii,iii); %lisätään laakerin termit jäykkyysmatriisiin
                end
                % tuennan solmu
                if ~Inp.Bearing(ii).Jnode==0
                    Jgdof=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); % % tuennan solmu, haetaan globaali Dof riville
                    if ~ (Jgdof==0) % jos ei rajoitettu
                        Kc2(Jgdof,Jgdof)=Kc2(Jgdof,Jgdof) + Kb(iii,iii); %lisätään laakerin termit jäykkyysmatriisiin
                    end
                end
                
                for jjj=1:3

                    % akselin solmu  (laakeri-laakeri)
                    Igdof_r=nnum2cgdof(Inp.Bearing(ii).Inode,iii); %Haetaan globaali Dof riville
                    Igdof_c=nnum2cgdof(Inp.Bearing(ii).Inode,jjj); %Haetaan globaali Dof sarakkeelle
                    if ~ (Igdof_r==0 | Igdof_c==0)
                        Kc(Igdof_r,Igdof_c)=Kc(Igdof_r,Igdof_c) + Kb(iii,jjj); %lisätään laakerin termit jäykkyysmatriisiin
                        Cc(Igdof_r,Igdof_c)=Cc(Igdof_r,Igdof_c) + Cb(iii,jjj); %lisätään laakerin termit vaimennusmatriisiin
                    end

                    if ~Inp.Bearing(ii).Jnode==0
                        % tuennan solmu (tuenta-tuenta)
                        Jgdof_r=nnum2cgdof(Inp.Bearing(ii).Jnode,iii); %Haetaan globaali Dof riville
                        Jgdof_c=nnum2cgdof(Inp.Bearing(ii).Jnode,jjj); %Haetaan globaali Dof sarakkeelle
                        if ~ (Jgdof_r==0 | Jgdof_c==0)
                            Kc(Jgdof_r,Jgdof_c)=Kc(Jgdof_r,Jgdof_c) + Kb(iii,jjj); %lisätään laakerin termit jäykkyysmatriisiin
                            Cc(Jgdof_r,Jgdof_c)=Cc(Jgdof_r,Jgdof_c) + Cb(iii,jjj); %lisätään laakerin termit vaimennusmatriisiin
                        end

                        % laakeri-tuenta ristitermit (yläkolmio)
                        if ~ (Igdof_r==0 | Jgdof_c==0)
                            Kc(Igdof_r,Jgdof_c)=Kc(Igdof_r,Jgdof_c) - Kb(iii,jjj); %lisätään laakerin termit jäykkyysmatriisiin
                            Cc(Igdof_r,Jgdof_c)=Cc(Igdof_r,Jgdof_c) - Cb(iii,jjj); %lisätään laakerin termit vaimennusmatriisiin
                        end
                        % tuenta-laakeri ristitermit (alakolmio)
                        if ~ (Jgdof_r==0 | Igdof_c==0)
                            Kc(Jgdof_r,Igdof_c)=Kc(Jgdof_r,Igdof_c) - Kb(iii,jjj); %lisätään laakerin termit jäykkyysmatriisiin
                            Cc(Jgdof_r,Igdof_c)=Cc(Jgdof_r,Igdof_c) - Cb(iii,jjj); %lisätään laakerin termit vaimennusmatriisiin
                        end
                    end
                end % jjj
                
            end %iii
        end % ii Bearing
        
        
        
        % Calculation of natural frequencies
        if strcmp(Request.Instability.Type,'Full')
            % laskenta täysillä matriiseilla
            
            % Gyromatriisi kerrotaan pyörimisnopeudella
            Gc=Spin(i)*Matrx.Gc;
            
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
            for k=1:length(idx2)
                Fiiapu(:,k)=Fii2(:,idx2(k));
            end
            Fii2=Fiiapu;  
            clear Fiiapu
            
            % Freq2 sisältää kompleksiset ominaistaajuudet 
            % Fii2 sisältää kompleksiset ominaismuodot
            % Etsitään vain positiivista ominaistaajuutta vastaavat muodot 
            kk2=0; jj2=0;
            nfreq=size(Freq2,1);
            while kk2 < 8  % +kahdeksan muotoa
                jj2=jj2+1;
                % Imaginääriosa on yhtä kuin taajuus Hz
                if imag(Freq2(jj2)) > 0.0 %
                    %  i = pyörimisnopeuden indeksi
                    %  kk2 = löydetyn taajuuden laskuri
                    kk2=kk2+1; 
                    FrqR(kk2,i)=imag(Freq2(jj2)); % Tallennetaan taajuudet
                    DampR(kk2,i)=-real(Freq2(jj2))/(abs(Freq2(jj2))); % Tallennetaan vaimennussuhteet
                    FiiR(:,kk2,i)=Fii2(nfreq/2+1:nfreq,jj2); % Pudotetaan toinen puolikas pois
                end 
            end
            clear kk2 jj2
            
        else
            % Laskenta speudo modal menetelmällä
            
            % Ominaisarvotehtavan ratkaisu
            [FiiN d]=eig(Kc2,Matrx.Mc);
            % Matriisien jarjestely taajuuden mukaan.
            [Freq, idx]=sort(sqrt(diag(d))/(2*pi));
            
            % Reaaliset ominaisarvot
            for k=1:length(idx)
                Fiiapu(:,k)=FiiN(:,idx(k));
            end
            FiiN=Fiiapu;
            clear Fiiapu idx k d
            
            % Muotojen maaran valinta (Nmuoto alinta ominaismuotoa)
            FiiNn=FiiN(:,1:Request.Instability.NroModes);
            
            % Massanormeeraus
            for k=1:Request.Instability.NroModes
                FiiNormN(:,k)=FiiNn(:,k)/sqrt(FiiNn(:,k)'*Matrx.Mc*FiiNn(:,k));
            end
            
            % Moodaaliset matriisit
            Mm=FiiNormN'*Matrx.Mc*FiiNormN;
            Km=FiiNormN'*Kc2*FiiNormN + FiiNormN'*(Kc-Kc2)*FiiNormN; % jälkimmäinen termi sisältää laakerin ristitermit
            
            % Gyromatriisi kerrotaan pyörimisnopeudella
            Gc=Spin(i)*Matrx.Gc;
            % vaimennusmatriisi + gyromatriisi
            Cm=FiiNormN'*(Cc+Gc)*FiiNormN;
            
            % Eigenvalue problem A*X=lambda*X
            zeroX=zeros(size(Mm));  eyeX=eye(size(Mm));
            AA=[  zeroX               eyeX
                -(Mm)^(-1)*Km      -(Mm)^(-1)*Cm];
            
            [Fii2 d2]=eig(AA); % State space ominaisarvotehtävän ratkaisu
            clear AA
            
            % Matriisien jarjestely taajuuden mukaan.
            [Freq2, idx2]=sort((diag(d2))/(2*pi));
            
            % Freq2 sisältää kompleksiset ominaistaajuudet 
            % Fii2 sisältää kompleksiset ominaismuodot
            % Etsitään vain positiivista ominaistaajuutta vastaavat muodot 
            kk2=0; jj2=0;
            nfreq=size(Freq2,1);
            while kk2 < 8 % 8 muotoa
                jj2=jj2+1;
                % Imaginääriosa on yhtä kuin taajuus Hz
                if imag(Freq2(jj2)) > 0.001 %
                    %  i = pyörimisnopeuden indeksi
                    %  kk2 = löydetyn taajuuden laskuri
                    kk2=kk2+1; 
                    FrqR(kk2,i)=imag(Freq2(jj2)); % Tallennetaan taajuudet
                    DampR(kk2,i)=-real(Freq2(jj2))/(abs(Freq2(jj2))); % Tallennetaan vaimennussuhteet
                    FiiR(:,kk2,i)=Fii2(nfreq/2+1:nfreq,jj2); % Pudotetaan toinen puolikas pois
                end 
            end
            
            clear kk2 jj2
            
        end % if full / pseudo modal
        
        
        
        
        
        
        
        % waitbar
        waitbar(Request.Instability.Speed_Toler/(Spin(3)-Spin(1)));
        
        %disp(['Speed ' num2str(i) ' of ' num2str(length(Request.Campbell.Spin_vec)) ' Processed...'])
        
    end % Spin_vec loop, indeksi i
    
    
    %DampR
    % Etsitään negatiivisten taajuuksien indeksit
    [I J]=find(DampR < 0); % 
    
    % Jos ei negatiivisia taajuuksia, niin kasvatetaan laskuria ja positutaan
    % loopista
    if isempty(J)
        lask=lask+1;
        break
    end
    
    % negatiivisia vaimennuksia voi olla useilla muodoilla. etsitään alin
    % taajuus, jonka pitäisi olla 1. FW
    minF=1e6; % alkuarvo alimmalle taajuudelle
    if exist('Jmin'); clear Jmin; end; %Tuhotaan vanha Jmin, jos olemassa
    for a=1:size(I,1)
        if FrqR(I(a),J(a)) < minF % jos taajuus pienempi kuin minF
            minF=FrqR(I(a),J(a));
            Imin=I(a);
            %Jmin=J(a); % 7.1.2013 havaittu bugi toiminnassa, korjaus
            %seuraavilla riveillä
            if exist('Jmin'); 
                Jmin=min([Jmin J(a)]); 
            else
                Jmin=J(a);
            end;
        end
    end
    
    % Puolitetaan nopeusväli
    if Jmin==2
        aa(1)=Spin(1);
        aa(2)=0.5*(Spin(1)+Spin(2));
        aa(3)=Spin(2);
        Spin=aa;
    elseif Jmin==3
        aa(1)=Spin(2);
        aa(2)=0.5*(Spin(2)+Spin(3));
        aa(3)=Spin(3);
        Spin=aa;
    else % jos vaimenuus ensimmäisellä nopeudella negatiivinen
        lask=lask+1;
    end
    
    
    
end %while

delete(h) 						% close waitbar

% Tulostetaan 
disp(' ')
disp('Stability Analysis Results:')
disp('---------------------------')

if lask==0 % Jos iteraatio löysi alimman rajanopeuden
    
    %rajanopeus rpm:nä
    RajaRPM=round(mean(Spin*60/(2*pi)));
        
    disp(['   + Threshold Speed: ' num2str(RajaRPM) ' rpm +/- ' num2str(Request.Instability.Speed_Toler*60/(2*pi)) ' rpm'])
    disp(['   + Whirl Frequency: ' num2str(round(FrqR(Imin,Jmin)*60)) ' cpm'])
    disp(['   + Frequency Ratio: ' num2str(FrqR(Imin,Jmin)*60/RajaRPM,4) ])
    
%     mutoojen tulostus rajanopeudella
Request.ModePlot.Spin_speed=RajaRPM*2*pi/60;
Request.ModePlot.ModeNro=Imin; 
Request.ModePlot.ScaleFactor=1e3; %7.1.2013 scale factor päivitetty 2--> 1e3. Johtuu aikaisemmasta muutoksesta f_ModeShapePlot funktiossa
Request.ModePlot.Type='Full'; % tai
% Request.ModePlot.Type='PseudoModal';  % PseudoModal ei vielä toimi
% Request.ModePlot.NroModes=30;

% for ii=1:4
%     Request.ModePlot.ModeNro=ii; 
    f_ModeShapePlot(Inp,Request);
% end
else 
    % Jos ei löydetty rajanopeutta
    disp(['No Threshold Speed Found in the Range of ' num2str(Spin(1)*60/(2*pi)) ' - ' num2str(Spin(3)*60/(2*pi)) ' rpm'])
end
disp(' ')
disp(' ')

