path_directory = '.../acquisition_data';

[FCS, I, panchro, DMD_conf, IC] = rebuild_aqucube(path_directory);


%% Script TV
%% Parameters and Constants
x0=0;
ml=-5;my=-2;mx=my;
mu_y=10^(my);mu_x=10^(mx);mu_lambda=10^(ml);

r_l=1200; c_l=r_l;
binsize=3;
spec_binsize=2;
BandStart=425;
BandNo=110;
regionx=400; %r,c, coords for region A
regiony=400; %r,c, coords for region A
region_label='full';
N=10;

CCD_ADU_noise = I;
%%
F = FCS;
[R,C,W,S] = size(F);
Filter = [];

for s = 1:S
    for k = 1:W
        slice = F(:,:,k,s);
        [i,j,~] = find(slice);
        M = length(i); N = length(j);
        cdd = mat2cell(slice,R,C);
        Filter = [Filter;cdd]; % Quelles dimensions pour le FC?
    end
end

Filter = reshape(Filter,W,S);


F = Filter;

%%
% On dispose des données aqucube et du cube de filtrage FiltCube_Sparse
% Attention à voir le bruit...

% Modèle adjoint :
% Entrée : z de taille (R,C,S)
% Sortie : ATz de taille (R,C,W)
aqucube = I;
[R,C,S] = size(aqucube);
W = size(F,1);
z = aqucube;
ATz = zeros(R,C,W);
for w=1:W
    Tmp = zeros(R,C);
    for s = 1:S
        Tmp = Tmp + F{w,s}.*z(:,:,s);
    end
    ATz(:,:,w)= Tmp;
end

% Modèle direct
% Entrée : o de taille (R,C,W)
% Sortie : z de taille (R,C,S)
o = ATz;
z = zeros(R,C,S);
for s=1:S
    Tmp = zeros(R,C);
    for w=1:W
        Tmp = Tmp + F{w,s}.*o(:,:,w);
    end
    z(:,:,s) = Tmp;
end
        
sig_noise = sqrt(aqucube)/binsize;
data = zeros(size(aqucube));
for s=1:S
    data(:,:,s) = aqucube(:,:,s)./squeeze(sig_noise(:,:,s));
end
for w=1:W
    for s=1:S
        Filter{w,s}=F{w,s}./squeeze(sig_noise(:,:,s));
    end
end
hand_K = @(o)direct(o,F);
hand_KT = @(z)adj(z,F);

%% TV
% addpath ~/Matlab/irntv
nmpdef;

pars_irn = irntvInputPars('l2tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.pcgtol_ini = 1e-4;

pars_irn.pcgtol_ini     = 1e-1;
pars_irn.adaptPCGtol    = 0;

pars_irn.problem = NMP_L2TV;
pars_irn.variant = NMP_TV_STANDARD;

NMP_WEIGHTS_THRESHOLD   = 30; % Fichier nmpdef.m
pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;

% I_Threshold = irntv(In, {}, 0.075, pars_irn);
% HC Opération Convolution
pars_irn.rrs=1e1;
pars_irn.loops=10;
%lambda = 0.075;
pars_irn.lambda = 1;
lambda = 1;

% for lambda = % %1e-3 1e-2 1e-1 1 10 100]
fprintf('---------------------------------\n')
fprintf('Lancement irntv_HC lambda = %.2e (%d it)\n',lambda,pars_irn.loops)
fprintf('---------------------------------\n')

% Initialisation
pars_irn.U0 = hand_KT(data);
% Minimization
t = tic;
I_Threshold = irntv_HC(data, {hand_K,hand_KT}, lambda, pars_irn);
toc(t)
e_time=toc(t);
fprintf('Sauvegarde\n')
eval(sprintf('save Result_TV_lambda%.2e_it%d.mat I_Threshold e_time \n', lambda,pars_irn.loops))

close all

% Peut-on relancer en prenant 
% for k = 1:9
%     pars_irn.U0 = I_Threshold ;
%     t = tic;
%     I_Threshold = irntv_HC(data, {hand_K,hand_KT}, lambda, pars_irn);
%     toc(t)
%     e_time=toc(t);
%     fprintf('Sauvegarde\n')
%     eval(sprintf('save Result_TV_lambda%.2e_it%d.mat I_Threshold e_time \n', lambda,pars_irn.loops))
% end

% [ssim_val,mapSSIM]= SSIM_map(GroundTruth,I_Threshold);
% [mapSam]= SAM_map(GroundTruth,I_Threshold);
%%
figure(1)
imagesc(sum(I_Threshold,3))
title('Panchro reconstructed')

figure(2)
plot(squeeze(I_Threshold(10,10,:)))
