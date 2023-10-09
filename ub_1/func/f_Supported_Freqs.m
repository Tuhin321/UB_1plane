function [FrqS,DampS,FiiS,Matrx]=f_Supported_Freqs(Inp,Request)
%% function calculates the supported natural frequencies of the system.
% The support system requires the usage of linearized bearing matrices.
% Note that the supported modes are at zero (0) rpm
% 5.7.2017 // JEH

[Matrx]=f_Rotor_Matrix(Inp);
nnum2cgdof=Matrx.nnum2cgdof;

%% Piirretään alkuperäinen rakenne ilman laakereita ja rajoitteita. Myös
% solmu ja elementtinumerot jätetään pois
if strcmp(Request.SuppModes.PlotModes,'Yes')
    %asetetaan default optiot
    az=0; el=90; %view(a,el) --> XY-näkymä
    origloc=[0 0 0]; %origon asema
    shade=0; %wireframe plot
    alpha=1; %pinnat ei läpinäkyviä

    % tarkistetaan onko optiot määritelty
    if isfield(Request.SuppModes,'Options')
        for i=1:size(Request.SuppModes.Options,1)
            switch Request.SuppModes.Options{i,1}
                case 'view'
                    switch Request.SuppModes.Options{i,2}
                        %case 'XY'
                        %    az=0; el=90; %view(a,el) --> XY-näkymä
                        case 'XZ'
                            az=0; el=0; %view(a,el) --> XZ-näkymä
                        case 'xy'
                            az=0; el=0; %view(a,el) --> XZ-näkymä
                        case '3D'
                            az=-37.5; el=30; %view(a,el) --> 3D-näkymä
                        case '3d'
                            az=-37.5; el=30; %view(a,el) --> 3D-näkymä
                    end
                case 'originloc'
                    origloc=Request.SuppModes.Options{i,2}; %origon asema
                case 'shading'
                    switch Request.SuppModes.Options{i,2};
                        case 'on'; shade=1;
                        case 'On'; shade=1;
                        case 'ON'; shade=1;    
                    end
                case 'transparency'
                    alpha=Request.SuppModes.Options{i,2};
            end
        end
    end

    % Draw Model as a 3D WireFrame plot % Funktio avaa uuden figuren
    if shade==0; f_RGeomPlotWireFr(Inp, [1 1 0 0 1], [0 0],{'-b',2}); end

    % Draw Model as a Shaded 3D plot
    if shade==1; f_RGeomPlotShaded(Inp,alpha); end
    

    
    % koordinaatiston piirto
    %f_draw_caxis(origloc, sz)
    f_draw_caxis(origloc, 0.1);
    
     % Turn warning temporaly off
    warning off MATLAB:hg:patch:RGBColorDataNotSupported   
    %figure settings
    view(az,el); axis equal; axis vis3d; axis off;   
    view3d rot;
    warning on MATLAB:hg:patch:RGBColorDataNotSupported   
end


% bearing parameters
Mc=Matrx.Mc;
Kc=Matrx.Kc;
Cc=Matrx.Cc;


% This loop using iii is taken from the f_Campbell. It includes the
% bearing stiffness properties into the system matrices in the defined
% locations
for ii=1:size(Inp.Bearing,2)
    
         % Calculation of bearing stiffness and damping matrix at zero rpm
         % // JEH
        [Kb,Cb]=f_Bearing_Matrix(Inp.Bearing(ii),0);

    
    % Calculation of bearing stiffness and damping matrix
    if ii==1
        Request.Bearing.position='drive side';
    else
        Request.Bearing.position='servise side';
    end
    
    for iii=1:3
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

% Eigenvalue problem A*X=lambda*X (State Space form)
% rotational speed zero --> no gyroscopic effect Gc = 0
zeroX=zeros(size(Mc));  
eyeX=eye(size(Mc));

AA=[  zeroX                           eyeX
    -(Mc)^(-1)*Kc      -(Mc)^(-1)*(Cc)];

% Ominaisarvotehtavan ratkaisu
[Fii2 d2]=eig(AA); % Gyromatriisi mukana
clear AA

% Matriisien jarjestely taajuuden mukaan.
[Freq2, idx2]=sort(diag(d2)/(2*pi));  %Freq=alfa +/- beta*i
%END OF TEST

% Kompleksiset ominaisarvot
Fii2=Fii2(:,idx2); %muotojen järjestely
kk=0; ii=0;
nfreq=size(Freq2,1);

while kk < Request.SuppModes.NroFreqs  % etsitään NroFreqs alinta taajuutta
    ii=ii+1;
    % Imaginääriosa on yhtä kuin taajuus Hz
    if imag(Freq2(ii)) > .01 % positiiviset taajuudet, ei jäykän kappaleen liikkeitä
        %kk = löydetyn taajuuden laskuri
        kk=kk+1;
        FrqS(kk)=imag(Freq2(ii)); % Tallennetaan taajuudet
        DampS(kk)=-real(Freq2(ii))/(abs(Freq2(ii))); % Tallennetaan vaimennussuhteet
        %Test: this was the original
        %         FiiS(:,kk)=Fii2(nfreq/2+1:nfreq,ii); % Pudotetaan toinen puolikas pois
        % This should be used with new formulation
        FiiS(:,kk)=Fii2(nfreq/2+1:nfreq,ii)/(Freq2(ii)*2*pi); % Pudotetaan toinen puolikas pois
    end
end
clear kk Fii2

% Tulostetaan taajuudet komentoikkunaan
disp(' ')
disp('Supported Frequencies and Damping Ratios:')
disp('-------------------------------------------------------')

for jj=1:Request.SuppModes.NroFreqs
    str=sprintf('%5.3f', DampS(jj)*100);
    disp(['Mode ' num2str(jj) ',  Frequency: ' num2str(FrqS(jj),8) ' Hz, Damping Ratio: ' str ' %'])
end
disp('-------------------------------------------------------')

%% Drawing the modes -------------------------------------------------------
if strcmp(Request.SuppModes.PlotModes,'Yes')
    
    Node=Inp.Node;
    Elem=Inp.Elem;
    % Optiot .ModePlot structuresta
    %                           Mode#, ScaleFactor, linetype, linewidth
    % Request.SuppModes.ModePlot={1, 1000, '-k', 2
    %                             2, 1000, '-b', 2  };
    for a=1:size(Request.SuppModes.ModePlot,1)

        % kopioidaan muodon tulostuksen parametrit
        mm=Request.SuppModes.ModePlot{a,1}; %mode number
        ScaleFactor=Request.SuppModes.ModePlot{a,2};
        ltype=Request.SuppModes.ModePlot{a,3};
        lwidth=Request.SuppModes.ModePlot{a,4};

        % Tulostetaan palkkielementit viivoina
        t=pi/2; % koska muodossa yleensää imaginääriosa merkittävä --> sin(t)=1
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
                % nämä rivit tekevät saman asian kuin edellinen for-if hässäkkä
                Iindx= (Elem(ii,2)==Node(:,1));
                Jindx= (Elem(ii,3)==Node(:,1));

                %solmunumeroa vastaava globaali DOF
                UXindxI=nnum2cgdof(Node(Iindx,1),1);
                UYindxI=nnum2cgdof(Node(Iindx,1),2);
                UZindxI=nnum2cgdof(Node(Iindx,1),3);
                RotXindxI=nnum2cgdof(Node(Iindx,1),4); % torsio vapausaste
                UXindxJ=nnum2cgdof(Node(Jindx,1),1);
                UYindxJ=nnum2cgdof(Node(Jindx,1),2);
                UZindxJ=nnum2cgdof(Node(Jindx,1),3);
                RotXindxJ=nnum2cgdof(Node(Jindx,1),4); % torsio vapausaste

                % Noden x-siirtymä.....................................
                if ~UXindxI==0 % jos ei rajoitettu gdof
                    uI=Node(Iindx,2)+...
                        ScaleFactor*(real(FiiS(UXindxI,mm))*cos(t) - (imag(FiiS(UXindxI,mm)))*sin(t)) ;
                else
                    uI=Node(Iindx,2);
                end
                if ~UXindxJ==0 % jos ei rajoitettu gdof
                    uJ=Node(Jindx,2)+...
                        ScaleFactor*(real(FiiS(UXindxJ,mm))*cos(t) - (imag(FiiS(UXindxJ,mm)))*sin(t)) ;
                else
                    uJ=Node(Jindx,2);
                end
                % Noden y-siirtymä.....................................
                if ~UYindxI==0 % jos ei rajoitettu gdof
                    vI=Node(Iindx,3)+...
                        ScaleFactor*(real(FiiS(UYindxI,mm))*cos(t) - (imag(FiiS(UYindxI,mm)))*sin(t)) ;
                else
                    vI=Node(Iindx,3);
                end
                if ~UYindxJ==0 % jos ei rajoitettu gdof
                    vJ=Node(Jindx,3)+...
                        ScaleFactor*(real(FiiS(UYindxJ,mm))*cos(t) - (imag(FiiS(UYindxJ,mm)))*sin(t)) ;
                else
                    vJ=Node(Jindx,3);
                end
                % Noden z-siirtymä.....................................
                if ~UZindxI==0 % jos ei rajoitettu gdof
                    wI=Node(Iindx,4)+...
                        ScaleFactor*(real(FiiS(UZindxI,mm))*cos(t) - (imag(FiiS(UZindxI,mm)))*sin(t)) ;
                else
                    wI=Node(Iindx,4);
                end
                if ~UZindxJ==0 % jos ei rajoitettu gdof
                    wJ=Node(Jindx,4)+...
                        ScaleFactor*(real(FiiS(UZindxJ,mm))*cos(t) - (imag(FiiS(UZindxJ,mm)))*sin(t)) ;
                else
                    wJ=Node(Jindx,4);
                end
                % Noden x-rotaatio (vääntö)...........................
                if ~RotXindxI==0 % jos ei rajoitettu gdof
                    thxI=Node(Iindx,4)+...
                        ScaleFactor*(real(FiiS(RotXindxI,mm))*cos(t) - (imag(FiiS(RotXindxI,mm)))*sin(t));
                else
                    thxI=Node(Iindx,4);
                end
                if ~RotXindxJ==0 % jos ei rajoitettu gdof
                    thxJ=Node(Jindx,4)+...
                        ScaleFactor*(real(FiiS(RotXindxJ,mm))*cos(t) - (imag(FiiS(RotXindxJ,mm)))*sin(t)) ;
                else
                    thxJ=Node(Jindx,4);
                end

                % muotojen tulostus-------------------------------
                if abs(thxI) > 1e-9 & abs(thxJ) > 1e-9 % onko kyseessä vääntö
                    % vääntömuodon tulostus
                    h1(a)=plot3( [uI uJ], [thxI thxJ ], [0 0],ltype,'linewidth',lwidth);
                    % otsikkotiedot
                    str=sprintf('%5.3f', DampS(mm)*100);
                    t1{a}=['Mode ' num2str(mm) ' (Torsion): \itf\rm = ' num2str(FrqS(mm),4) ' Hz, \it\zeta\rm = ' str ' %'];
                else
                    % X, Y ja Z siirtymien tulostus
                    h1(a)=plot3( [uI uJ], [vI vJ], [wI wJ ],ltype,'linewidth',lwidth);
                    % otsikkotiedot
                    str=sprintf('%5.3f', DampS(mm)*100);
                    t1{a}=['Mode ' num2str(mm) ': \itf\rm = ' num2str(FrqS(mm),8) ' Hz, \it\zeta\rm = ' str ' %'];
                end
            end %if
        end %for elements
    end %for modes
    % Turn warning temporaly off
    warning off MATLAB:hg:patch:RGBColorDataNotSupported
    title({ Inp.model_title; 'Supported Modes'},'Fontsize',14)
    if exist('h1')
        legend(h1,t1, 'Location', 'NorthEast');
    end
    warning on MATLAB:hg:patch:RGBColorDataNotSupported
end %if PlotModes

% this plots the support nodes in case of of each studied mode. Please
% observe the relative behavior of each support DOF
if strcmp(Request.SuppModes.PlotSupport,'Yes')
    Supp1=nnum2cgdof(Inp.SpringDamper(1,2),2);
    figure
    plot(1:length(FiiS(1,:)),imag(FiiS(Supp1:Supp1+3,:)))
    grid on
    legend('1','2','3','4')
end

