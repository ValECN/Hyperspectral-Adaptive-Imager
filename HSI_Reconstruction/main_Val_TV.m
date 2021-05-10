%% Reconstruction using Total Variation regularization
% Original author: Valentin NOËL based on Dr.Rodriguez and Dr.Wohlberg's
% work (NUMIPAD library)

addpath 'irntv'
addpath 'Functions'
addpath Visualisation/

param_REC.data = 'A_PSF_10';

switch param_REC.data
    case 'K-poon'
        path_directory = [pwd, '/acquisition_data_1'];       
    case 'A'
        path_directory = [pwd, '/acquisition_data'];
    case 'A_PSF'
        path_directory = [pwd, '/acquisition_data_PSF'];      
    case 'A_PSF_10'
        path_directory = [pwd, '/acquisition_data_PSF10'];  
end

nmpdef;

[FCS, I, panchro, DMD_conf, IC] = rebuild_aqucube(path_directory);
% [I,FCS] = data_processing();

%% Script TV

[R,C,W,S] = size(FCS);
Filter = [];

for s = 1:S
    for k = 1:W
        slice = sparse(FCS(:,:,k,s));
        [i,j,~] = find(slice);
        M = length(i); N = length(j);
        cdd = mat2cell(slice,R,C);
        Filter = [Filter;cdd]; 
    end
end

F = reshape(Filter,W,S);
        
hand_K = @(o)direct(o,F);
hand_KT = @(z)adj(z,F);

%% TV

pars_irn = irntvInputPars('l2tv');

pars_irn.epsR         = 10;
pars_irn.epsF         = 0;
pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 1;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 0;
pars_irn.epsF_cutoff  = 0;   % This is the percentage cutoff
pars_irn.pcgtol_ini   = 0.1;
pars_irn.adaptPCGtol  = 1;

pars_irn.problem = NMP_L2TV;
pars_irn.variant = NMP_TV_STANDARD;

NMP_WEIGHTS_THRESHOLD   = 30; % Fichier nmpdef.m
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;

% HC Opération Convolution
pars_irn.rrs=1.5e1;
pars_irn.loops=20;
pars_irn.lambda_ini = 0.1;
lambda = 20;

fprintf('---------------------------------\n')
fprintf('Lancement irntv_HC lambda = %.2e (%d it)\n',lambda,pars_irn.loops)
fprintf('---------------------------------\n')

% Initialisation
pars_irn.U0 = hand_KT(I);

% Minimization
t = tic;
I_Threshold = irntv_HC(I, {hand_K,hand_KT}, lambda, pars_irn);
toc(t)
e_time=toc(t);

%% --- Plots
figure(1)
imagesc(sum(I_Threshold,3))
title('Panchro reconstructed')
xlabel('X\_cam')
ylabel('Y\_cam')

figure(2)
x1 = 8; y1 = 13; x2 = 8; y2 = 24;
plot(squeeze(I_Threshold(x1,y1,:)))
hold on 
plot(squeeze(I_Threshold(x2,y2,:)))
title('Spectra of the reconstructed cube using TV')
xlabel('Bandwidth')
ylabel('Amplitude')
legend(sprintf('Spectra Reconstruction for x = %d and y = %d',x1,y1),sprintf('Spectra Reconstruction for x = %d and y = %d',x2,y1))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

%% --- Visualization
param_Visualization.run = 1;
if param_Visualization.run == 1
    % -- Saving the reconstructed cube as the .mat to use for Visualization
    Object = I_Threshold;
    save('Object.mat','Object')

    % -- Visualization
    Visual_HyperSpectral_2019(I_Threshold)
end

%% --- Metrics
IC_bin = spectral_binning(IC,4);
IC_bin = IC_bin(:,:,1:W);
[ssim_val, ssim_map] = SSIM_map(IC_bin, I_Threshold);
[RMSE_map] = RMSE_map(IC_bin, I_Threshold);
[SAM_map] = SAM_map(IC_bin, I_Threshold);
