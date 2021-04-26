%% Reconstruction using Total Variation regularization
% Original author: Valentin NOËL based on Dr.Rodriguez and Dr.Wohlberg's
% work (NUMIPAD library)

addpath 'irntv'
addpath 'Functions'
addpath Visualisation/
path_directory = [pwd, '/acquisition_data'];

[FCS, I, panchro, DMD_conf, IC] = rebuild_aqucube(path_directory);

%% Script TV

[R,C,W,S] = size(FCS);
Filter = [];

for s = 1:S
    for k = 1:W
        slice = sparse(FCS(:,:,k,s));
        [i,j,~] = find(slice);
        M = length(i); N = length(j);
        cdd = mat2cell(slice,R,C);
        Filter = [Filter;cdd]; % Quelles dimensions pour le FC?
    end
end

F = reshape(Filter,W,S);
        
hand_K = @(o)direct(o,F);
hand_KT = @(z)adj(z,F);

%% TV
% addpath ~/Matlab/irntv
nmpdef;

pars_irn = irntvInputPars('l2tv');

pars_irn.epsR         = 1;
pars_irn.epsF         = 0;
pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 100;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 0;
pars_irn.epsF_cutoff  = 0;   % This is the percentage cutoff
pars_irn.pcgtol_ini   = 0.1;
pars_irn.adaptPCGtol  = 1;

pars_irn.problem = NMP_L2TV;
pars_irn.variant = NMP_TV_LamdaAdapt;

NMP_WEIGHTS_THRESHOLD   = 30; % Fichier nmpdef.m
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;

% HC Opération Convolution
pars_irn.rrs=1e1;
pars_irn.loops=10;
pars_irn.lambda_ini = 1;
lambda = 100;

fprintf('---------------------------------\n')
fprintf('Lancement irntv_HC lambda = %.2e (%d it)\n',lambda,pars_irn.loops)
fprintf('---------------------------------\n')

% Initialisation
pars_irn.U0 = hand_KT(I);
% pars_irn.U0 = [];

% Minimization
t = tic;
I_Threshold = irntv_HC(I, {hand_K,hand_KT}, lambda, pars_irn);
% I_Threshold = irntv_HC(I, {}, lambda, pars_irn);
toc(t)
e_time=toc(t);

%%
figure(1)
imagesc(sum(I_Threshold,3))
title('Panchro reconstructed')
xlabel('X\_cam')
ylabel('Y\_cam')

figure(2)
plot(squeeze(I_Threshold(12,12,:)))
title('Spectra of the reconstructed HSI for x = y = 10')
xlabel('Bandwidth')
ylabel('Amplitude')

param_Visualization.run = 0;
if param_Visualization.run == 1
    % -- Saving the reconstructed cube as the .mat to use for Visualization
    Object = I_Threshold;
    save('Object.mat','Object')

    % -- Visualization
    Visual_HyperSpectral_2019(I_Threshold)
end
