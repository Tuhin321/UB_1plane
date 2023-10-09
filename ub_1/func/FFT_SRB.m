function FFT=FFT_SRB(Inp,tt,FF)
% FFT for global data
% close all
% clear all

% Signal parameters
Fs=1/Inp.tstep;    % sampling frequency [Hz] 
tend=Inp.tend;      % sample duration [s]
tstart=0.02;  % start time for analyzing data

ps=Fs*tend; % amount of data points in the sample

% Number of FFT points (2^n powers. Needs to be smaller than ps
% 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144 etc.  
NFFT=Inp.NFFT;


% resolution of frequency
Fres=(1/(NFFT/Fs));

% generating test function (amplitude and frequency is known)
% t=0:1/Fs:tend; % time vector
% y=100*sin(10*2*pi*t)+200*sin(20*2*pi*t)+150*sin(30*2*pi*t);

% data into matrix (will be replaced with simulated data)
%D=[t' y'];

% load Drop_v8.mat 
% node=4;
% % simulated time
%  D=[tout 1e3*Xdot((node-1)*4+2,:)'];
D = [ tt FF];
% global D
% printing frequency Hz
%# freqmin=5; freqmax=300; % max frequency is Fs/2
freqmin=0; freqmax=300; % max frequency is Fs/2
% indeces of frequency vector (for printing)
sfrq_ind=max([ 1 floor(freqmin/Fres)]);
efrq_ind=floor(freqmax/Fres);
% Cutting sample to be same length than NFFT, the start time tstart


X1=D(floor(tstart*Fs)+1 : floor(tstart*Fs)+NFFT,2);
% X1=D(0:100,2);
%  Calculating FFT in function using Hanning windowing
[MX,freq]=fft_hanning(X1,NFFT,Fs);




% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figure (6)  
% plot settings
set(gcf,'Units','normalized')		   % Units normalized (always)
set(gcf,'Position',[0.1 0.1 0.8 0.8])  % Position set to given
set(gcf,'Color',[1 1 1])			   % Background color white (always)

% Time plane plotting. Sample window is drawn as red box
% subplot(2,1,1)
% plot(D(:,1), D(:,2))
% hold on
% v=axis; % getting axis info of plot
% % x and y coordinates of the sample
% xb=[D(floor(tstart*Fs)+1,1) D(floor(tstart*Fs)+1,1) ...
%     D(floor(tstart*Fs)+NFFT,1) D(floor(tstart*Fs)+NFFT,1) ...
%     D(floor(tstart*Fs)+1,1)];
% yb=[v(3) v(4) v(4) v(3) v(3)];
% plot(xb,yb,'r','linewidth',2)
% xlabel('Time [s]')
% ylabel('Force Amplitude [N]')
% grid on

% subplot(2,1,2)
% Spectrum plotting

% semilogy(freq(sfrq_ind:efrq_ind),MX(sfrq_ind:efrq_ind),'-k', 'linewidth',2)
plot(freq(sfrq_ind:efrq_ind),MX(sfrq_ind:efrq_ind),'-k', 'linewidth',2)


% hold on 
% for jj=1:Inp.z
%    
%     text( Inp.Omega_cage/(2*pi)*jj, 20,' X ');
% end


xlabel('Frequency [Hz]')
ylabel('Amplitude')
grid on
end
