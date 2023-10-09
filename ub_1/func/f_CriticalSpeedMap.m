%% this file plots so called 'critical speed map' using the bearing
% locations that are defined in INDATA. The bearings are assumed to be
% isotropic and the bearing stiffness at both ends are simultaneously
% changing as a function of 'pp' in the main for-loop
% 5.7.2017 // JEH

count=0;

for pp=2:0.1:11
    
    count=count+1;
    
    [Matrx]=f_Rotor_Matrix(Inp);
    nnum2cgdof=Matrx.nnum2cgdof;
    
    Kc=Matrx.Kc; 
    Cc=Matrx.Cc;
    Kc2=Matrx.Kc; % Lalannen merkinnˆill‰ K* matriisi (ei sis‰ll‰ laakereiden ristitermej‰)
    Mc=Matrx.Mc;

    Kbb(count)=1*10^(pp);
   
    
    Kb=[0 0 0
        0 Kbb(count) 0
        0 0 Kbb(count)];
    
    Cb=0*Kb;
    
    % Calculation of Bearing Matrices
    for ii=1:size(Inp.Bearing,2)
        
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
        end
    end
            
            
            
    % this is the new FORMAT FROM CAMPBELL
    % Eigenvalue problem A*X=lambda*X (State Space muoto)
    zeroX=zeros(size(Mc));  eyeX=eye(size(Mc));
    AA=[  zeroX                           eyeX
        -(Mc)^(-1)*Kc      -(Mc)^(-1)*(Cc)];
    
    % Ominaisarvotehtavan ratkaisu
    [Fii2 d2]=eig(AA); % Gyromatriisi mukana
    clear AA
    
    % Matriisien jarjestely taajuuden mukaan.
    [Freq2, idx2]=sort(diag(d2)/(2*pi));  %Freq=alfa +/- beta*i
    %END OF TEST
    
    % Kompleksiset ominaisarvot
    Fii2=Fii2(:,idx2); %muotojen j‰rjestely
    kk=0; ii=0;
    nfreq=size(Freq2,1);
    
    while kk < 10%Request.FreeModes.NroFreqs  % etsit‰‰n NroFreqs alinta taajuutta
        ii=ii+1;
        % Imagin‰‰riosa on yht‰ kuin taajuus Hz
        if imag(Freq2(ii)) > .01 % positiiviset taajuudet, ei j‰yk‰n kappaleen liikkeit‰
            %kk = lˆydetyn taajuuden laskuri
            kk=kk+1;
            FrqR(kk,count)=imag(Freq2(ii)); % Tallennetaan taajuudet
            DampR(kk,count)=-real(Freq2(ii))/(abs(Freq2(ii))); % Tallennetaan vaimennussuhteet
            %Test: this was the original
            %         FiiR(:,kk)=Fii2(nfreq/2+1:nfreq,ii); % Pudotetaan toinen puolikas pois
            % This should be used with new formulation
            FiiR(:,kk)=Fii2(nfreq/2+1:nfreq,ii)/(Freq2(ii)*2*pi); % Pudotetaan toinen puolikas pois
        end
    end
    clear kk Fii2
    
end

figure
semilogx(Kbb,FrqR(:,:))
grid on
