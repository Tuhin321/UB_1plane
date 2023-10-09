function [FrqR,DampR,FiiR,Matrx]=f_FreeModes(Inp,Request)
% Function f_FreeModes ignores bearing definitions in the Inp structure and
% calculates free-free modes of the system (bending, axial and torsion).
%
% Input:   Inp structure containing model data
%          Request structure array with the following fields:
%          Request.FreeModes.NroFreqs  = Number of frequencies to be solved and plotted
%          Request.FreeModes.PlotModes = Option for mode shape plot ('Yes'
%             or 'No'
%          Request.FreeModes.ModePlot = Structural array containing mode
%             plot parameters. Syntax is: {Mode#, ScaleFactor, linetype,
%             linewidth}. Example: {1, 1000, '-k',  2
%                                   6, 10000, '--g', 2 }
%          Request.FreeModes.Options = options for model plot. Format is:
%             {'property1',value}. Valid properties and values are:
%             'view': 'XY'(view XY-plane,default), 'XZ'(view XZ-plane), '3D'(3D-view)
%             'originloc': vector of origintriad location, default [0,0,0]
%             'shading': 'on'(plot surfaces), otherwise wireframe plot (default)
%             'transparency': face transpacency if 'shading','on'. Value
%             between 0..1 (default 1 --> not transparent)
%
% Output:  FrqR = vector of frequencies
%          DampR =vector of damping ratios
%          FiiR = modeshape matrix
%          Matrx = structure of rotor matrices
%
% Written by Jussi Sopanen, LUT/IMVE, 2004-2006

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

% Piirret‰‰n alkuper‰inen rakenne ilman laakereita ja rajoitteita. Myˆs
% solmu ja elementtinumerot j‰tet‰‰n pois
if strcmp(Request.FreeModes.PlotModes,'Yes')
    %asetetaan default optiot
    az=0; el=90; %view(a,el) --> XY-n‰kym‰
    origloc=[0 0 0]; %origon asema
    shade=0; %wireframe plot
    alpha=1; %pinnat ei l‰pin‰kyvi‰

    % tarkistetaan onko optiot m‰‰ritelty
    if isfield(Request.FreeModes,'Options')
        for i=1:size(Request.FreeModes.Options,1)
            switch Request.FreeModes.Options{i,1}
                case 'view'
                    switch Request.FreeModes.Options{i,2}
                        %case 'XY'
                        %    az=0; el=90; %view(a,el) --> XY-n‰kym‰
                        case 'XZ'
                            az=0; el=0; %view(a,el) --> XZ-n‰kym‰
                        case 'xy'
                            az=0; el=0; %view(a,el) --> XZ-n‰kym‰
                        case '3D'
                            az=-37.5; el=30; %view(a,el) --> 3D-n‰kym‰
                        case '3d'
                            az=-37.5; el=30; %view(a,el) --> 3D-n‰kym‰
                    end
                case 'originloc'
                    origloc=Request.FreeModes.Options{i,2}; %origon asema
                case 'shading'
                    switch Request.FreeModes.Options{i,2};
                        case 'on'; shade=1;
                        case 'On'; shade=1;
                        case 'ON'; shade=1;    
                    end
                case 'transparency'
                    alpha=Request.FreeModes.Options{i,2};
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

% Ratkaistaan muodot ja taajuudet-----------------------------------------

% TEST: THIS WAS THE ORIGINAL FORMAT-------------
% % Generalized eigenvalue problem A*X=lambda*B*X
% zeroX=zeros(size(Matrx.Mc));
% AA=[zeroX   Matrx.Mc
%     Matrx.Mc Matrx.Cc];
% BB=[Matrx.Mc   zeroX
%     zeroX  -Matrx.Kc];
% % Ominaisarvotehtavan ratkaisu
% [Fii2 d2]=eig(AA,BB); % Gyromatriisi mukana
% clear AA BB
% % Matriisien jarjestely taajuuden mukaan.
% [Freq2, idx2]=sort(1./diag(d2)/(2*pi));

% this is the new FORMAT FROM CAMPBELL
        % Eigenvalue problem A*X=lambda*X (State Space muoto)
        zeroX=zeros(size(Matrx.Mc));  eyeX=eye(size(Matrx.Mc));
        AA=[  zeroX                           eyeX
            -(Matrx.Mc)^(-1)*Matrx.Kc      -(Matrx.Mc)^(-1)*(Matrx.Cc)];
        
        % Ominaisarvotehtavan ratkaisu 
        [Fii2 d2]=eig(AA); % Gyromatriisi mukana
        clear AA
        
        % Matriisien jarjestely taajuuden mukaan.
        [Freq2, idx2]=sort(diag(d2)/(2*pi));  %Freq=alfa +/- beta*i 
%END OF TEST



% Kompleksiset ominaisarvot
Fii2=Fii2(:,idx2); %muotojen j‰rjestely
% for i=1:length(idx2)
%     Fiiapu(:,i)=Fii2(:,idx2(i));
% end
% Fii2=Fiiapu;
% clear Fiiapu

% Freq2 sis‰lt‰‰ kompleksiset ominaistaajuudet
% Fii2 sis‰lt‰‰ kompleksiset ominaismuodot
% Etsit‰‰n vain positiivista ominaistaajuutta vastaavat muodot
kk=0; ii=0;
nfreq=size(Freq2,1);
while kk < Request.FreeModes.NroFreqs  % etsit‰‰n NroFreqs alinta taajuutta
    ii=ii+1;
    % Imagin‰‰riosa on yht‰ kuin taajuus Hz
    if imag(Freq2(ii)) > 0.005 % positiiviset taajuudet, ei j‰yk‰n kappaleen liikkeit‰ 
        %kk = lˆydetyn taajuuden laskuri
        kk=kk+1;
        FrqR(kk,1)=imag(Freq2(ii)); % Tallennetaan taajuudet
        DampR(kk,1)=-real(Freq2(ii))/(abs(Freq2(ii))); % Tallennetaan vaimennussuhteet
        %Test: this was the original
%         FiiR(:,kk)=Fii2(nfreq/2+1:nfreq,ii); % Pudotetaan toinen puolikas pois
        % This should be used with new formulation
        FiiR(:,kk)=Fii2(nfreq/2+1:nfreq,ii)/(Freq2(ii)*2*pi); % Pudotetaan toinen puolikas pois
    end
end
clear kk Fii2

% Tulostetaan taajuudet komentoikkunaan
disp(' ')
disp('Free-Free Frequencies and Damping Ratios:')
disp('-------------------------------------------------------')

for jj=1:Request.FreeModes.NroFreqs
    str=sprintf('%5.3f', DampR(jj)*100);
    disp(['Mode ' num2str(jj) ',  Frequency: ' num2str(FrqR(jj),4) ' Hz, Damping Ratio: ' str ' %'])
end
disp('-------------------------------------------------------')

% Muotojen piirto--------------------------------------------------------
if strcmp(Request.FreeModes.PlotModes,'Yes')
    Node=Inp.Node;
    Elem=Inp.Elem;
    % Optiot .ModePlot structuresta
    %                           Mode#, ScaleFactor, linetype, linewidth
    % Request.FreeModes.ModePlot={1, 1000, '-k', 2
    %                             2, 1000, '-b', 2  };
    for a=1:size(Request.FreeModes.ModePlot,1)

        % kopioidaan muodon tulostuksen parametrit
        mm=Request.FreeModes.ModePlot{a,1}; %mode number
        ScaleFactor=Request.FreeModes.ModePlot{a,2};
        ltype=Request.FreeModes.ModePlot{a,3};
        lwidth=Request.FreeModes.ModePlot{a,4};

        % Tulostetaan palkkielementit viivoina
        t=pi/2; % koska muodossa yleens‰‰ imagin‰‰riosa merkitt‰v‰ --> sin(t)=1
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
                % n‰m‰ rivit tekev‰t saman asian kuin edellinen for-if h‰ss‰kk‰
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

                % Noden x-siirtym‰.....................................
                if ~UXindxI==0 % jos ei rajoitettu gdof
                    uI=Node(Iindx,2)+...
                        ScaleFactor*(real(FiiR(UXindxI,mm))*cos(t) - (imag(FiiR(UXindxI,mm)))*sin(t)) ;
                else
                    uI=Node(Iindx,2);
                end
                if ~UXindxJ==0 % jos ei rajoitettu gdof
                    uJ=Node(Jindx,2)+...
                        ScaleFactor*(real(FiiR(UXindxJ,mm))*cos(t) - (imag(FiiR(UXindxJ,mm)))*sin(t)) ;
                else
                    uJ=Node(Jindx,2);
                end
                % Noden y-siirtym‰.....................................
                if ~UYindxI==0 % jos ei rajoitettu gdof
                    vI=Node(Iindx,3)+...
                        ScaleFactor*(real(FiiR(UYindxI,mm))*cos(t) - (imag(FiiR(UYindxI,mm)))*sin(t)) ;
                else
                    vI=Node(Iindx,3);
                end
                if ~UYindxJ==0 % jos ei rajoitettu gdof
                    vJ=Node(Jindx,3)+...
                        ScaleFactor*(real(FiiR(UYindxJ,mm))*cos(t) - (imag(FiiR(UYindxJ,mm)))*sin(t)) ;
                else
                    vJ=Node(Jindx,3);
                end
                % Noden z-siirtym‰.....................................
                if ~UZindxI==0 % jos ei rajoitettu gdof
                    wI=Node(Iindx,4)+...
                        ScaleFactor*(real(FiiR(UZindxI,mm))*cos(t) - (imag(FiiR(UZindxI,mm)))*sin(t)) ;
                else
                    wI=Node(Iindx,4);
                end
                if ~UZindxJ==0 % jos ei rajoitettu gdof
                    wJ=Node(Jindx,4)+...
                        ScaleFactor*(real(FiiR(UZindxJ,mm))*cos(t) - (imag(FiiR(UZindxJ,mm)))*sin(t)) ;
                else
                    wJ=Node(Jindx,4);
                end
                % Noden x-rotaatio (v‰‰ntˆ)...........................
                % EDIT 16.1.2014 by JSo: removed Node(Iindx,4)+... from
                % thxI and thxJ calculation
                if ~RotXindxI==0 % jos ei rajoitettu gdof
                    thxI=ScaleFactor*(real(FiiR(RotXindxI,mm))*cos(t) - (imag(FiiR(RotXindxI,mm)))*sin(t));
                else
                    thxI=Node(Iindx,3);
                end
                if ~RotXindxJ==0 % jos ei rajoitettu gdof
                    thxJ=ScaleFactor*(real(FiiR(RotXindxJ,mm))*cos(t) - (imag(FiiR(RotXindxJ,mm)))*sin(t)) ;
                else
                    thxJ=Node(Jindx,3);
                end

                % muotojen tulostus-------------------------------
                if abs(thxI) > 1e-9 && abs(thxJ) > 1e-9 % onko kyseess‰ v‰‰ntˆ
                    % v‰‰ntˆmuodon tulostus
                    % EDIT 16.1.2014 by JSo: added [Node(Iindx,3)+thxI
                    % Node(Jindx,3)+thxJ ]. Now branched torsional system can be plotted correctly 
                    h1(a)=plot3( [uI uJ], [Node(Iindx,3)+thxI Node(Jindx,3)+thxJ ], [0 0],ltype,'linewidth',lwidth);
                    % otsikkotiedot
                    str=sprintf('%5.3f', DampR(mm)*100);
                    t1{a}=['Mode ' num2str(mm) ' (Torsion): \itf\rm = ' num2str(FrqR(mm),4) ' Hz, \it\zeta\rm = ' str ' %'];
                else
                    % X, Y ja Z siirtymien tulostus
                    h1(a)=plot3( [uI uJ], [vI vJ], [wI wJ ],ltype,'linewidth',lwidth);
                    % otsikkotiedot
                    str=sprintf('%5.3f', DampR(mm)*100);
                    t1{a}=['Mode ' num2str(mm) ': \itf\rm = ' num2str(FrqR(mm),4) ' Hz, \it\zeta\rm = ' str ' %'];
                end
             end %if
        end %for elements
    end %for modes
    % Turn warning temporaly off
    warning off MATLAB:hg:patch:RGBColorDataNotSupported
    title({ Inp.model_title; 'Free-Free Modes'},'Fontsize',14)
    if exist('h1')
        legend(h1,t1, 'Location', 'South');
    end
    warning on MATLAB:hg:patch:RGBColorDataNotSupported
end %if PlotModes
% axis tight;
















