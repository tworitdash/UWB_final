clear;clc;
close all;

Sfactork = [1];

%fck = 5e9:1e9:21e9;
% 
% dely = zeros(size(S));
% delc = zeros(size(Bk));



for k = 1:size(Sfactork, 2)


%% Load data
 %load Sphere_LinearSAR.mat;
load Scissor_LinearSAR.mat;
% load Scissor_LinearSAR.mat;
%% Data Parameters
% Number of Transmit/Receive antennas
NTRx = length(rawdata(:,1));
% Number of freNfrequency points
Nfre = length(rawdata(1,:));
% Start and end frequency of the VNA data set
Fstart = 1e9;
Fstop = 21e9;
%% Select frequency band from the original data set
%fc = fck(k);
%fc = fck(k);
fc = 16e9;
B = 10e9;
%% MIMO Topology configuration
TRx = zeros(NTRx,2); 
TRx(:,1) = [-70e-2:1e-2:70e-2].';
%% imaging parameters
% Set Image area
% X-cross range; Y-Range;
focX = [-0.25 0.25];
focY = [0.2 0.6];

% focX = [-0.6 0.6];
% focY = [0.2 0.7];


% Image resolution
detas = 1e-2;
% Reduce sampling by the factor of
Sfactor = Sfactork(k);
% Dynamic range of image display
dynRng2D = 20;
%%
c = (21 - k)/25;
color = [c, c, 1, 1];
Final_image;


end
% %print(['Q2_cross_range_vs_GHz'], '-depsc');
% figure(150);
% plot(Bk*10^(-9), dely*10^(3), 'LineWidth', 3, 'color', [0.6350, 0.0780, 0.1840]);
% hold on;
% plot(Bk*10^(-9), c./(2 .* Bk).*10^(3), 'LineWidth', 3, 'color', [0.25, 0.25, 0.25]);
% legend({'Observation', 'Theoretical'}, 'FontSize', 12, 'FontWeight', 'bold');
% %xlabel('Bandwidth (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Bandwidth (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Down Range resolution [mm]', 'FontSize', 12, 'FontWeight', 'bold');
% box on;
% %title(strcat('Center frequency:', num2str(fc./1e9),'GHz, Bandwidth:', num2str(B./1e9),'GHz, Sampling step:', num2str(Sfactor),'cm'));
% %title('Down Range resolution vs central frequency f_c');
% title('Down Range resolution vs Bandwidth', 'FontSize', 12, 'FontWeight', 'bold');
% grid on;
% print('Q2_dr_2', '-depsc');
% 
% figure(151);
% plot(Bk*10^(-9), delc*10^(3), 'LineWidth', 3, 'color', [0.6350, 0.0780, 0.1840]);
% hold on;
% plot(Bk*10^(-9), ((0.4/1.4).*c./fc).*10^(3).*Bk./Bk, 'LineWidth', 3, 'color', [0.25, 0.25, 0.25]);
% legend({'Observation', 'Theoretical'}, 'FontSize', 12, 'FontWeight', 'bold');
% xlabel('Bandwidth (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% %xlabel('Bandwidth (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Cross Range resolution [mm]' , 'FontSize', 12, 'FontWeight', 'bold');
% box on;
% %title(strcat('Center frequency:', num2str(fc./1e9),'GHz, Bandwidth:', num2str(B./1e9),'GHz, Sampling step:', num2str(Sfactor),'cm'));
% %title('Cross Range resolution vs central frequency f_c', 'FontSize', 12, 'FontWeight', 'bold');
% title('Cross Range resolution vs Bandwidth', 'FontSize', 12, 'FontWeight', 'bold');
% grid on;
% print('Q2_cr_2', '-depsc');