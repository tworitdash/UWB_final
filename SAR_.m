clear;clc;
close all;

Sfactork = 1;
BW = 10e9;
%BW = 12e9;
%% 
for k = 1:size(Sfactork, 2)
    for l = 1:size(BW, 2)

load Sphere_LinearSAR.mat;
 %load Scissor_LinearSAR.mat;
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
fc = 11e9;
B = BW(l);
%% MIMO Topology configuration
TRx = zeros(NTRx,2); 
TRx(:,1) = [-70e-2:1e-2:70e-2].';
%% imaging parameters
% Set Image area
% X-cross range; Y-Range;
% focX = [-0.25 0.25];
% focY = [0.2 0.6];
% 
focX = [-0.1 0.1];
focY = [0.35 0.5];

%focX = [-0.6 0.6];
%focY = [0.2 0.7];


% Image resolution
detas = 0.1e-2;
% Reduce sampling by the factor of
Sfactor = Sfactork(k);
% Dynamic range of image display
dynRng2D = 20;
%%
c = (5 - l)/7;
color = [c, c, 1, 1];
Image_;
%Image_2;
%Image_3;
    end

end
% figure(3);
% legend(legendInfo, 'FontSize', 12, 'FontWeight', 'bold');
% %print('Singular_Spectra_BW_changing', '-depsc');
% print('Singular_Spectra_Spacing_changing', '-depsc');
% figure(4);
% legend(legendInfo, 'FontSize', 12, 'FontWeight', 'bold');
% %print('Singular_Spectra_BW_changing_dB', '-depsc');
% print('Singular_Spectra_Spacing_changing_dB', '-depsc');

