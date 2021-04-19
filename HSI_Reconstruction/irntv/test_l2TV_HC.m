%% Modifications HC pour tests opérateurs

nmpdef;

% HC Opération Convolution
H = ones(10); H = H/sum(H(:));
hand_K = @(x)conv2(x,H,'full');
hand_KT = @(x)conv2(x,flip(H),'valid');



I = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

% In = I + 0.05*randn(size(I));
% HC Opération Convolution
In = hand_K(I);
In = In + 0.05*randn(size(In));


%% --------------------------
% --- Adapt cutoff value ---

pars_irn = irntvInputPars('l2tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 5;
%pars_irn.U0         = In;
% HC Opération Convolution
pars_irn.U0         = hand_KT(In);


pars_irn.pcgtol_ini     = 1e-1;
pars_irn.adaptPCGtol    = 0;

pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;

t = tic;
% I_Threshold = irntv(In, {}, 0.075, pars_irn);
% HC Opération Convolution
pars_irn.rrs=1e1;
pars_irn.loops=10;
lambda = 0.075;
lambda = 0.0081;
I_Threshold = irntv_HC(In, {hand_K,hand_KT}, lambda, pars_irn);
toc(t)
subplot(2,2,3)
imagesc(I_Threshold); title('Threshold')


snr(I, I_Threshold)

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
I_mil = irntv_HC(In, {hand_K,hand_KT}, 0.075, pars_irn);
toc(t)

snr(I, I_mil)

%% Add HC affichage des résultats

subplot(2,2,1)
imagesc(I); title('Original'); colormap gray
subplot(2,2,2)
imagesc(In); title('Noisy')
subplot(2,2,3)
imagesc(I_Threshold); title('Threshold')
subplot(2,2,4)
imagesc(I_mil); title('Matrix Inversion Lemma')


