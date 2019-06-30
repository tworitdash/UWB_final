%%



%load Scissor_LinearSAR.mat;
    
c = 2.998e8;
qj = sqrt(-1);
%% 
TRx = TRx(1:Sfactor:end,:);
rawdata = rawdata(1:Sfactor:end,:);
NTRx = length(TRx(:,1));

% figure;
% scatter(TRx(:,1),TRx(:,2),40,'k*','LineWidth',2);
% grid on;axis image;
% axis([1.2*(-0.7) 1.2*(0.7) 1.2*(-0.2) 1.2*(0.2)]);
% title(strcat('Array Topology,',' NTRx=',num2str(NTRx)));xlabel('X [m]');ylabel('Z [m]');box on;
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
%Nfre = 1000;
%f = f(fminn:fmaxn);
%f = linspace(fmin, fmax, 100);
detaf = 100e6;
f = fmin:detaf:fmax;

Nfre = length(f);

rawdata_2 = rawdata(:,fminn:fmaxn);

%% 

NfocXY = [round(abs(focX(2)-focX(1))/detas)+1 round(abs(focY(2)-focY(1))/detas)+1];
X = linspace(focX(1),focX(2),NfocXY(1));
Y = linspace(focY(1),focY(2),NfocXY(2));
a = size(X, 2); b = size(Y, 2); ab = a * b;

Observ = zeros(NTRx*Nfre, ab);
%Data = zeros(NTRx*Nfre, 1);

Data = rawdata(:);

[X_i, Y_i] = meshgrid(X, Y);

    for m = 1:NTRx
        txi = TRx(m,:);
        d_i = sqrt((X_i-txi(1))^2 + (Y_i-txi(2))^2);
        d = d_i(:);
        for p = 1:Nfre
            f0 = f(p);
            k0 = (2 * pi * f0)/c;
            Observ((m-1)*Nfre+p, :) = (exp(-1j*k0*d)/(4 * pi * d)).^2;
        end
    end

% Image

Image_vec = pinv(Observ) * Data;
Image_mat = zeros(a, b);

 for yj = 1:NfocXY(2)
        for xj = 1:NfocXY(1)
            Image_mat(xj, yj) = Image_vec((yj - 1)*a+xj);
        end
 end
 
figure(l + 100);
imagesc(X, Y, (abs(Image_mat))); shading flat;colormap('gray');
xlabel('X [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y [m]', 'FontSize', 12, 'FontWeight', 'bold');
title(['Sphere Image, BW = ', num2str(B*10^(-9)), 'GHz, Spacing = ', num2str(Sfactor), 'cm, f_c = ', num2str(fc*10^(-9)), ' GHz'],...
    'FontSize', 12, 'FontWeight', 'bold');
% %print(['Image_BW_', num2str(B*10^(-9)), '_GHz'], '-depsc');
% print(['Spacing_', num2str(Sfactor), '_cm'], '-depsc');


%% SVD 

[U, S, V] = svd(Observ);
figure(3);
hold on;
plot((diag(S)/max(diag(S))), 'LineWidth', 3, 'Color', color);
grid on;
hold on;
%legendInfo{k} = ['Spacing = ',  num2str(Sfactor/(lambdac*100)), '\lambda cm'];
%legendInfo{l} = ['BW = ',  num2str(B*10^(-9)), ' GHz'];
xlabel('Serial Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Singular Value linear scale Normalized', 'FontSize', 12, 'FontWeight', 'bold');
title(['Singular value at BW = ', num2str(B*10^(-9)), 'GHz f_c = ', num2str(fc*10^(-9)), ' GHz'], ...
    'FontSize', 12, 'FontWeight', 'bold');
%title(['Spacing = ', num2str(Sfactor), 'cm f_c = ', num2str(fc*10^(-9)), ' GHz'], ...
    %'FontSize', 12, 'FontWeight', 'bold');
% 
% figure(4);
% hold on;
% plot(db(diag(S)/max(diag(S))), 'LineWidth', 3, 'Color', color);
% grid on;
% hold on;
% legendInfo{k} = ['Spacing = ',  num2str(Sfactor/(lambdac*100)), '\lambda cm'];
% %legendInfo{l} = ['BW = ',  num2str(B*10^(-9)), ' GHz'];
% xlabel('Serial Number', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Singular Value dB scale Normalized', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Singular value at BW = ', num2str(B*10^(-9)), 'GHz f_c = ', num2str(fc*10^(-9)), ' GHz'], ...
%     'FontSize', 12, 'FontWeight', 'bold');
% %title(['Spacing = ', num2str(Sfactor), 'cm f_c = ', num2str(fc*10^(-9)), ' GHz'], ...
%     %'FontSize', 12, 'FontWeight', 'bold');
% 
% 
% %% Reconstructing image
% 
% % A_new = U(1:1500, 1:1500)*S(1:1500, 1:1500)*V(1:1500, 1:1500)';
% % 
% %Image_vec = pinv(A_new) * Data(1:1500, :);
% % Sum_i = zeros(ab, 1);
% % for i = 1:2091
% %     sigma_i = S(i, i);
% %     U_i = U(:, i)';
% %     V_i = V(:, i);
% %     Sum_i = Sum_i + (1/sigma_i) * U_i * Data * V_i;
% % end
% % % S_dig = diag(S);
% % % Image_vec = (diag(S_dig(1:1400))^(-1)\U(:, 1:1400)' * Data) * V(:, 1:1400);
% % Image_vec = Sum_i;
% % 
% % Image_mat = zeros(a, b);
% % 
% %  for xj = 1:NfocXY(1)
% %         for yj = 1:NfocXY(2)
% %             Image_mat(xj, yj) = Image_vec((xj - 1)*b+yj);
% %         end
% %  end
% %  figure;
% %  imagesc(X, Y, db(abs((Image_mat)))); shading flat;colormap('gray');
% %  
% %  
% %  %X_ls = V * (S(1:ab, :))\U' * Data;

%% With dominant singular vlaues:
m = 323;

X_ls = V(:, 1:m) * inv(S(1:m, 1:m)) * U(:, 1:m)' * Data;
Image_mat = reshape(X_ls, a, b);
% for xj = 1:NfocXY(1)
%         for yj = 1:NfocXY(2)
%             Image_mat(xj, yj) = X_ls((xj - 1)*b+yj);
%         end
% end
 figure;
 imagesc(X, Y, (abs(Image_mat))); shading flat;colormap('gray');
