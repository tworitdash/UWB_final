%%



%load Sphere_LinearSAR.mat;
    
c = 2.998e8;
qj = sqrt(-1);
%% 
TRx = TRx(1:Sfactor:end,:);
rawdata = rawdata(1:Sfactor:end,:);
NTRx = length(TRx(:,1));

figure;
scatter(TRx(:,1),TRx(:,2),40,'k*','LineWidth',2);
grid on;axis image;
axis([1.2*(-0.7) 1.2*(0.7) 1.2*(-0.2) 1.2*(0.2)]);
title(strcat('Array Topology,',' NTRx=',num2str(NTRx)));xlabel('X [m]');ylabel('Z [m]');box on;
% %% 
f = linspace(Fstart,Fstop,Nfre);
detaf = f(2)-f(1);
fmin = fc-B/2; 
fmax = fc+B/2; 
lambda = c/fc;
fminn = find(abs(f-fmin)==min(abs(f-fmin)));
fmaxn = find(abs(f-fmax)==min(abs(f-fmax)));
fmin = f(fminn);
fmax = f(fmaxn);
disp(strcat('Centre freq', num2str(fc), 'Operational Bandwidth ', num2str(B./1e9), ' GHz ', 'Start Frequency: ',num2str(fmin./1e9),'GHz, End frequency: ', num2str(fmax./1e9),'GHz'));
f = f(fminn:fmaxn);
Nfre = length(f);
%rawdata = rawdata(:,fminn:fmaxn);
rawdata = rawdata(:,fminn:fmaxn);

%%
winsize = round(Nfre*1); 
windows = kaiser(winsize,3).';
Windowh = meshgrid(windows,[1:NTRx]);
rawdata = rawdata.*Windowh; clear Windowh;
rawdata = permute(rawdata,[2 1]); 
%% Transform to Time-domain
detat = 5e-12;
NF = round((1/detat)/detaf +1);
t = [0:detat:(NF-1)*detat].';
Temp = zeros(NF,NTRx);
nstart = round(fmin/detaf+1);
nstartc = NF - round(fmin/detaf)+1;
Temp(nstart:nstart+Nfre-1,:) = rawdata;
Temp(nstartc:-1:nstartc-Nfre+1,:) = conj(rawdata);
rawdatat = ifft(Temp,[],1);
clear temp;
Tmax = 10e-9;
Nt = round(Tmax/detat+1);
rawdatat = rawdatat(1:Nt,:);
t = t(1:Nt);

range = c .* t/2;

figure;
imagesc(TRx(:,1),t./1e-9,rawdatat);
xlabel('X [m]'); ylabel('Time [ns]');
title('Data in time-doman');
colorbar;grid on;

figure;
imagesc(TRx(:,1),range,rawdatat);
xlabel('X [m]'); ylabel('Range [m]');
title('Data in time-doman');
colorbar;grid on;
%%
NfocXY = [round(abs(focX(2)-focX(1))/detas)+1 round(abs(focY(2)-focY(1))/detas)+1];
X = linspace(focX(1),focX(2),NfocXY(1));
Y = linspace(focY(1),focY(2),NfocXY(2));
Image = zeros(NfocXY);
Green_func = zeros(NfocXY);
Observation_matrix = ones(NfocXY(2), 1);
Antenna_data = zeros(NfocXY);
tic
sigma = ones(1001, 1);
for trxi = 1:NTRx
    txi = TRx(trxi,:);
     sig = squeeze(rawdatat(:,trxi));
%     sig_freq = squeeze(rawdata(:, trxi));
    for xi = 1:NfocXY(1)
        for yi = 1:NfocXY(2)
            D = 2*sqrt((X(xi)-txi(1))^2 + (Y(yi)-txi(2))^2);
            %G = exp(-1j.*(4*pi/lambda).*d)/(4 * pi * d);
            %freq = round(d/c*detaf) + 1;
            tpr = round(D/c/detat)+1;
            %Green_func(xi, yi) = G;
            Image(xi,yi) = Image(xi,yi) + sig(tpr);
        end
    end
%     Observation_matrix = Green_func.' * sigma;
%     P = rawdata(1:1001, trxi) * Observation_matrix.';
%     Antenna_data = Antenna_data + P;
%     
    detat = 5e-12;
%     NF = round((1/detat)/detaf +1);
%     t = [0:detat:(NF-1)*detat].';
%     Temp = zeros(NF,size(Y, 2));
%     nstart = round(fmin/detaf+1);
%     nstartc = NF - round(fmin/detaf)+1;
%     Temp(nstart:nstart+Nfre-1,:) = Antenna_data;
%     Temp(nstartc:-1:nstartc-Nfre+1,:) = conj(Antenna_data);
%     Antenna_data_ifft = ifft(Temp, [], 1);
end

% figure;
% imagesc(X, Y, db(abs(Antenna_data_ifft)), [-dynRng2D 0]); shading flat;
%colormap('gray');colorbar;axis image;axis xy;grid on;
% figure;
% imagesc(X, Y, abs(Green_func)); shading flat;
toc
%% Plot figures
Imageab = abs(Image);
Image_nor = Imageab/max(max(Imageab));
Image_dB = 20*log10(Image_nor);
figure;
imagesc(X,Y,Image_dB.',[-dynRng2D 0]);
colormap('gray');colorbar;axis image;axis xy;grid on; 
%caxis([-30 0]);
xlabel('X [m]');ylabel('Y [m]');box on;
title(strcat('Center frequency:', num2str(fc./1e9),'GHz, Bandwidth:', num2str(B./1e9),'GHz, Sampling step:', num2str(Sfactor),'cm'));
%print(['Q4_Image_Sf_', num2str(Sfactor), '_cm'], '-depsc');

Imageab = abs(Image);
Projection = max(Imageab,[],2);
Projection_nor = Projection./max(Projection);
Projection_dB = 20*log10(Projection_nor);

% ind2 = X(find(Projection_dB >= max(Projection_dB) - 5));
% 
% delc(k) = max(ind2) - min(ind2);

%half_power1 = find(Projection_dB)

% figure(100);
% plot(X,Projection_dB,'Linewidth',2);
% %xlim(focX);
% grid on;
% hold on;
% xlabel('X [m]');ylabel('Beam pattern [dB]');
% box on;
% legendInfo{k} = ['Spacing = ',  num2str(Sfactor), ' cm'];
% 
% title(strcat('Center frequency:', num2str(fc./1e9),'GHz, Bandwidth:', num2str(B./1e9),'GHz, Sampling step:', num2str(Sfactor),'cm'));

%% Down range

Imageab_2 = abs(Image);
Projection_2 = max(Imageab_2,[],1);
Projection_nor_2 = Projection_2./max(Projection_2);
Projection_dB_2 = 20*log10(Projection_nor_2);
%half_power1 = find(Projection_dB)


% ind = interp1(Projection_dB_2, Y, max(Projection_dB_2) - 3,'nearest');
% index = find(Y == ind);
% 
% temp = Projection_dB_2;
%  
% temp(index-2:index+2) = eps.*[-1, -2, -3, -4, -5];
% % 
% ind2 = interp1(temp,Y,-3, 'linear');

% ind = Y(find(Projection_dB_2 >= max(Projection_dB_2) - 5));
% dely(k) = max(ind) - min(ind);

% ind3 = ind(find(ind > 0.3960 & ind < 0.4020));
% 
% dely(k) = max(ind3) - min(ind3);

% figure(200);
% hold on;
% plot(Y,Projection_dB_2,'k','Linewidth',2, 'color', rand(1, 3));
% xlim(focY);
% grid on;
% xlabel('Y [m]');ylabel('Down Range vs Intensity [dB]');box on;
% title(strcat('Center frequency:', num2str(fc./1e9),'GHz, Bandwidth:', num2str(B./1e9),'GHz, Sampling step:', num2str(Sfactor),'cm'));
% figure(100);
%legend(legendInfo, 'FontSize', 12, 'FontWeight', 'bold');
