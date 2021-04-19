%%
clear
close all
clc
%%

%% Modifications HC pour tests opérateurs
% Test pour images couleur

nmpdef;

% HC Opération Convolution
H = ones(9); H=H/sum(H(:));

% ------- Start Local modification Val
hand_K = @(x)conv2col(x,H,'full');
hand_KT = @(x)conv2col(x,fliplr(flip(H)),'valid');
% hand_K = @(x)conv2col(x,H,'full');
% hand_KT = @(x)conv2col(x,flip(H),'valid'); %%% A VERIFIER s'il ne faut pas un fliplr et un flipud

%I = double( imread('gray_imgs/lena_gray_512.png') ) / 255;
I = double( imread('color_imgs/lena_color_512.png') ) / 255;

n = size(I,3);
for i = 1:n
    I(:,:,n+i) = I(:,:,(n+1)-i);
end

% I(:,:,4) = I(:,:,3);
% I(:,:,5) = I(:,:,2);
% I(:,:,6) = I(:,:,1);
% ------- End Local Modification Val

% In = I + 0.05*randn(size(I));
% HC Opération Convolution
In = hand_K(I);
In = In + 0.05*randn(size(In));

%%
figure(1)
subplot(3,2,1)
imagesc(I(:,:,1)); colorbar; colormap gray
axis square;
axis equal;
axis image;
subplot(3,2,2)
imagesc(In(:,:,1)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,3)
imagesc(I(:,:,2)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,4)
imagesc(In(:,:,2)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,5)
imagesc(I(:,:,3)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,6)
imagesc(In(:,:,3)); colorbar
axis square;
axis equal;
axis image;

figure(2)
subplot(2,2,1)
imagesc(I(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(2,2,2)
imagesc(In(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(2,2,3)
imagesc(I(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;
subplot(2,2,4)
imagesc(In(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;

%% --------------------------
% --- Adapt cutoff value ---

pars_irn = irntvInputPars('l2tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
%%%% pars_irn.adapt_epsF   = 1;
%%%% pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.pcgtol_ini = 1e-4;

%pars_irn.U0         = In;
% HC Opération Convolution
pars_irn.U0         = hand_KT(In);


%pars_irn.pcgtol_ini     = 1e-1;
pars_irn.adaptPCGtol    = 0;

pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;

t = tic;
% I_Threshold = irntv(In, {}, 0.075, pars_irn);
% HC Opération Convolution
pars_irn.rrs=1e1;
pars_irn.loops=10;
%lambda = 0.075;
% lambda = 0.008;
lambda = 0.0049; % obtained from the trade-off lambda
I_Threshold = irntv_HC(In, {hand_K,hand_KT}, lambda, pars_irn);
toc(t)
% ------- Start Local Modification Val

% Peut-on relancer en mettant 
n = 10;
for i = 1:n
    pars_irn.U0         = I_Threshold;
    I_Threshold = irntv_HC(In, {hand_K,hand_KT}, lambda, pars_irn);
end

% for i = 1 : n
%     lambda = 0.0005 + i * 0.0001;
%     I_Threshold = irntv_HC(In, {hand_K,hand_KT}, lambda, pars_irn);
%     %pars_irn.U0         = I_Threshold;
%     totSNR(i) = snr(I,I_Threshold);
%     totRMSE(i) = sqrt(mean(mean((sum(I,3) - sum(I_Threshold,3)).^2)));
% end
    
SNR1 = snr(I, I_Threshold);
% % SAM = jmsam(I, I_Threshold);
% % 
% [val1,lbdOpti1] = max(totSNR);  % opti SNR = 15.2556 for lbd = 0.0076
% [val2,lbdOpti2] = min(totRMSE); % opti RMSE = 0.1970 for lbd = 0.0034
% 
% [val3,lbdOpti3] = sort(totSNR,'descend');  
% [val4,lbdOpti4] = sort(totRMSE,'ascend');
% 
% for i = 1:n
%     maxi = 2*n;
%     sumIndex = find(lbdOpti3 == i) + find(lbdOpti4 == i);
%     if sumIndex < maxi
%         maxi = sumIndex;
%         lbdTradeOff = find(lbdOpti3 == i); % trade-off lambda = 0.0049
%     end
% end

% figure(5)
% 
% subplot 211
% plot(totSNR)
% title('Evolution of the SNR')
% xlabel('Iterations')
% ylabel('SNR')
% 
% subplot 212
% plot(totRMSE)
% title('Evolution of the RMSE')
% xlabel('Iterations')
% ylabel('RMSE')

% ------ End Local Modifications Val

%%
figure(3)
subplot(3,2,1)
imagesc(I(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,2)
imagesc(I(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,3)
imagesc(In(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,4)
imagesc(In(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,5)
imagesc(I_Threshold(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,6)
imagesc(I_Threshold(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;

%% --------------------------
% --- Fixed cutoff value ---


pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;



% ------------------------
% ---   Add epsilon    ---


pars_irn.weight_scheme = NMP_WEIGHTS_EPSILON;

% ------------------------------
% --- Matrix inversion lemma ---


pars_irn = irntvInputPars('l2tv');

pars_irn.pcgtol_ini = 1e-4;
pars_irn.loops      = 5;
%pars_irn.U0         = In;
% HC Opération Convolution
pars_irn.U0         = hand_KT(In);

pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;

pars_irn.pcgtol_ini     = 1e-1;
pars_irn.adaptPCGtol    = 0;

t = tic;

pars_irn.rrs=1e1;
pars_irn.loops=10;

lambda = 0.16; % new reference value
I_mil = irntv_HC(In, {hand_K,hand_KT}, lambda, pars_irn);

% -------Start Local Modifications Val
% n = 100;
% for i = 1:n
%     lambda = 0.15 + i*0.001; % ------- 0.008 initiallement
%     I_mil = irntv_HC(In, {hand_K,hand_KT}, lambda, pars_irn);
%     toc(t)
%     totSNR(i) = snr(I,I_mil);
%     totRMSE(i) = sqrt(mean(mean((sum(I,3) - sum(I_mil,3)).^2)));
% end

% figure(6)
% 
% subplot 211
% plot(totSNR)
% title('Evolution of the SNR')
% xlabel('Iterations')
% ylabel('SNR')
% 
% subplot 212
% plot(totRMSE)
% title('Evolution of the RMSE')
% xlabel('Iterations')
% ylabel('RMSE')

SNR = snr(I, I_mil);
sam = SAM_map(I,I_mil);
ssim = SSIM_map(I,I_mil);

% ------- End Local Modifications Val

%% Add HC affichage des résultats

figure(4)
subplot(3,2,1)
imagesc(I(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,2)
imagesc(I(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,3)
imagesc(In(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,4)
imagesc(In(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,5)
imagesc(I_mil(:,:,1:3)); colorbar
axis square;
axis equal;
axis image;
subplot(3,2,6)
imagesc(I_mil(:,:,4:6)); colorbar
axis square;
axis equal;
axis image;

%% ------- Start Local Modifications Val
% ------------------------------
% --- Non-Negative Quadratic Programming ---


pars_irn = irntvInputPars('l2tv');

pars_irn.pcgtol_ini = 1e-4;
pars_irn.loops      = 5;
%pars_irn.U0         = In;
% HC Opération Convolution
% pars_irn.U0         = hand_KT(In);
pars_irn.U0         = [];

pars_irn.variant       = NMP_TV_NQP;
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;

pars_irn.pcgtol_ini     = 1e-1;
pars_irn.adaptPCGtol    = 0;

pars_irn.vmax_NQP    = 20;      % has to be tuned
pars_irn.gamma_NQP   = 1e-2;    % has to be between 1e-3 and 5e-1
pars_irn.alpha_NQP     = 0.75;  % has to be between 0.5 and 1
pars_irn.loops_NQP    = 20;     % has to be tuned

t = tic;
pars_irn.rrs=1e1;
pars_irn.loops=10;

n = 100;
for i = 1:n
    lambda = 0.15 + i*0.001; % ------- 0.008 initiallement
    I_NQP = irntv_HC(In, {hand_K,hand_KT}, lambda, pars_irn);
    toc(t)
    totSNR(i) = snr(I,I_NQP);
    totRMSE(i) = sqrt(mean(mean((sum(I,3) - sum(I_NQP,3)).^2)));
end

figure(7)

subplot 211
plot(totSNR)
title('Evolution of the SNR')
xlabel('Iterations')
ylabel('SNR')

subplot 212
plot(totRMSE)
title('Evolution of the RMSE')
xlabel('Iterations')
ylabel('RMSE')

SNR2 = snr(I, I_NQP);
sam = SAM_map(I,I_NQP);
ssim = SSIM_map(I,I_NQP);

% ------ End Local Modifications Val
