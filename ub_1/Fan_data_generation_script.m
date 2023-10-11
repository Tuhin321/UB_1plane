% Tuhin, sept 2023|| using Template for using kalman filter
% ==================================================================
% The robedyn main file is directly called as a function, which in turn uses the indata
% as a function to determine the system matrices using a function f_systemMatrices.
% 
% A variable 'label' is added to the request which will alter the ub in the model.
% based on that, the response amplitude (max amp) at certain speeds will be
% saved for creating dataset for deep learning model

%% the rest of the text is INVALID for this package
% for the true model Request.badParameter= 0;
% 1= bad mass % 2= bad shaft dia % 3 = bad modulus of elasticity  % 4= bad bearing coefficients  % 5= bad support stiffness
% IMP: when using Jeffcott without support
% (IndataJeffcottGoodandBadWoSupp), option 5 as bad parameter is invalid
% 3. Note that even when both good and bad models have same value as 'badParameter', the
%  response would still differ for both because of different initial values
%===============================================================
close all;clearvars;clc
restoredefaultpath
addpath ('additional func')
addpath func

% parallel pooling
% feature('numcores') %PARALLEL 
% mypool=parpool(4);
% profile on


%% Provide inputs for indata
omega_range=0:20:30000*2*pi/60;  % sped range

% for ub magnitude: generate 1000 labels
a = 0;
mag_max = 1e-2;
mag_1p = (mag_max-a).*rand(100,1) + a;
% mag_range = [min(mag_1p) max(mag_1p)]

% for ub phase: generate 1000 labels
phase_max = 360;
phase_1p = (phase_max-a).*rand(100,1) + a;
% phase_range = [min(phase_1p) max(phase_1p)]

labels_1p =[mag_1p,phase_1p];

%% Define request varaibles
Request.ubplanes= '1p'; % always keep it zero to load good model here
Request.show_plots='No';
Request.UnbalanceResponse='Yes';
% Note: Request.UBResp.OutputInGlobalYZDirections and Request.UBResp.OutputInProbeDirections, both cannot be 'Yes'
Request.UBResp.OutputInGlobalYZDirections='Yes';   % In the global Y, Z directions 'Yes' or 'No'
Request.UBResp.OutputInProbeDirections='No';   % In the directions of proximity probes 'Yes' or 'No'
Request.UBResp.ProxiProbeAngle=pi/4; % Phase angle (in radian) for the direction of the proximity probe in the direction of rotation (probes are 90 degrees apart).
Request.UBResp.PlotOrbitInYZplane='Yes';   % For plotting orbits of the defined nodes in the Y-Z plane
Request.UBResp.SpinSpeed=10000*2*pi/60; % For Plotting the Orbits: Rotational Speed must be multiple of Speed Range's Time Step and within its limit.
Request.UBResp.PlotRelativeAmp='No';  % Relative displacement is calculated in the direction of proximity probes. It is the displacement of each rotor node in the location of proximity probe relative to the corresponding support node.
Request.UBResp.Spin_vec=omega_range;    % Speed range
Request.UBResp.OutputNodes=[2 9 22 28]; % List of Output Nodes

% Note: Equal number of Rotor/Bearing/Sensor Nodes and Support Nodes should be defined
% Request.UBResp.RotorNodes=[5, 17];         % List of Rotor/Bearing/Sensor Nodes
% Request.UBResp.SupportNodes=[21, 22];      % List of corresponding/respective Support Nodes

%% Introduce the random generated labels to MAIN code to generate response (features)
% waitbar
h = waitbar(0,'Generating dataset...');

tic
features_matrix=cell (size(labels_1p,1),1);
for i=1: size(labels_1p,1)
features_matrix{i}= robedynMain(Request,labels_1p(i,:));
waitbar(i/size(labels_1p,1))
end
delete(h)
% Features_cells= mat2cell(features_matrix,[4],[1 1 1 1].*158)
toc

metadata= [repelem(labels_1p(:,1),size(omega_range,2)),repelem(labels_1p(:,2), size(omega_range,2)), (repmat(omega_range, 1, size(labels_1p,1)))',cat(2,features_matrix{:,:})'];
% end

% T= table (metadata(:1),metadata(:2), metadata(:end))
metadata=metadata(randperm(size(metadata, 1)), :);

writematrix(metadata,'Dataset_v1_1.csv') 

return

% % %% use plots for later 
% % % plots 
% % figure
% % set(gcf,'Units','normalized')		% Units normalized (always)
% % set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
% % set(gcf,'Color',[1 1 1])			% Background color white (always)
% % 
% % Inp.outDof= [1 2 5 6 9 10 13 14 17 18];
% % 
% % for jj= 1:size(Inp.outDof,2)
% %     i= Inp.outDof(jj);
% %     subplot(5,4,i)
% %     plot(time, Estim(:,i),'-.b', 'linewidth', 1);hold on
% %     plot(time, True(:,i),'--g','linewidth', 1);hold on
% %     plot(time,Badmodel(:,i).*1e6,'--r', 'linewidth', 1) % brg2, z
% %     xlabel('time')
% %     ylabel('displacment ( \mum)')
% %     title(['dof ',num2str(i)]);
% %     xlim([0.5e-4,1])
% %     
% %     subplot(5,4,i+2)
% %     plot(time, Estim(:,i+20),'-.b', 'linewidth', 1);hold on
% %     plot(time, True(:,i+20),'--g','linewidth', 1);hold on
% %     plot(time,Badmodel(:,i+20).*1e6,'--r', 'linewidth', 1) % brg2, z
% %     grid minor
% %     xlabel('time')
% %     ylabel('velocity ( \mum/s)')
% %     title(['dof ',num2str(i)]);
% %     if i==1
% %         legend('kalman', 'true','bad model')
% %     else
% %     end
% %     xlim([0.5e-4,1])
% % end
% % 
% % [str] = f_alteredPara(Request);
% % sgtitle(['Displacement at all dofs with measurement at dofs: ', num2str(Inp.obs1),' and ', num2str(Inp.obs2),'. Bad parameter = ', str])
% % % for i= 1:Inp.ndof
% % %     subplot(5,4,i)
% % %     plot(time, Estim(:,i),'-.b', 'linewidth', 1);hold on
% % %     plot(time, True(:,i),'--g','linewidth', 1);hold on
% % %     plot(time,Badmodel(:,i).*1e6,'--r', 'linewidth', 1) % brg2, z
% % %    grid minor
% % %     title(['dof ',num2str(i)]);
% % %     if i==1
% % %         legend('kalman', 'true','bad model')
% % %     else
% % %     end
% % %     xlabel('time')
% % %     ylabel('displacment ( \mum)')
% % % end
% % % 
% % % [str] = f_alteredPara(Request);
% % % sgtitle(['Displacement at all dofs with measurement at dofs: ', num2str(Inp.obs1),' and ', num2str(Inp.obs2),'. Bad parameter = ', str])
% % 
% % cyan = [0, 0.75, 0.75];
% % % plot errors at all dofs
% % figure
% % set(gcf,'Units','normalized')		% Units normalized (always)
% % set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
% % set(gcf,'Color',[1 1 1])			% Background color white (always)
% % for jj= 1:size(Inp.outDof,2)
% %     i= Inp.outDof(jj);
% %     subplot(5,4,i)
% %     hold on
% %     pos_coVar = 1.96.*sqrt(coVar(:,i).');
% %     inte      = fill([time.', fliplr(time.')], [-pos_coVar, fliplr(+pos_coVar)],cyan,'LineStyle','-','FaceAlpha',0.5,'edgecolor',cyan); % everything must be a row vector for fill
% %     plot(time, Estim(:,i)-True(:,i),'-k', 'linewidth', 1);
% %     xlabel('time')
% %     ylabel('displacment ( \mum)')
% %     title(['dof ',num2str(i)]);
% %     xlim([5e-4,1])
% %     
% %     subplot(5,4,i+2)
% %     hold on
% %     vel_coVar = 1.96.*sqrt(coVar(:,i+20).');
% %     inte      = fill([time.', fliplr(time.')], [-vel_coVar, fliplr(+vel_coVar)],cyan,'LineStyle','-','FaceAlpha',0.5,'edgecolor',cyan); % everything must be a row vector for fill
% %     plot(time, Estim(:,i+20)-True(:,i+20),'-k', 'linewidth', 1);
% %     grid minor
% %     xlabel('time')
% %     ylabel('velocity ( \mum/s)')
% %     title(['dof ',num2str(i)]);
% %     xlim([5e-4,1])
% %     if i==1
% %         legend('kalman', 'true','bad model')
% %     else
% %     end
% %     
% % end
% % [str] = f_alteredPara(Request);
% % sgtitle(['Error in prediciton (Estimator - Good Model). Measured dofs: ', num2str(Inp.obs1),' and ', num2str(Inp.obs2),'. Bad parameter = ', str])
% % toc
% % 
% % %delete(gcp);
% % %exit;
% % %connector off
% % % delete(mypool); %PARALLEL 
% % profsave
% % profile off
% % 
% % % figure
% % % set(gcf,'Units','normalized')		% Units normalized (always)
% % % set(gcf,'Position',[0.1 0.1 0.8 0.8])			% Position set to given
% % % set(gcf,'Color',[1 1 1])			% Background color white (always)
% % % 
% % % for i= 1:Inp.ndof
% % %     subplot(5,4,i)
% % %     
% % %     plot(time, Estim(:,i)-True(:,i),'-k', 'linewidth', 1);hold on
% % %     grid minor
% % %     title(['dof ',num2str(i)]);
% % %     if i==1
% % %         legend('error')
% % %     else
% % %     end
% % %     xlabel('time')
% % %     ylabel('displacment ( \mum)')
% % % end
% % % 
% % % [str] = f_alteredPara(Request);
% % % sgtitle(['Error in prediciton (Estimator - Good Model). Measured dofs: ', num2str(Inp.obs1),' and ', num2str(Inp.obs2),'. Bad parameter = ', str])
% %  