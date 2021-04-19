

nmpdef;


I = double( imread('gray_imgs/lena_gray_512.png') ) / 255;

In = I + 0.05*randn(size(I));


% --------------------------
% --- Adapt cutoff value ---

pars_irn = irntvInputPars('l2tv');

pars_irn.adapt_epsR   = 1;
pars_irn.epsR_cutoff  = 0.01;   % This is the percentage cutoff
pars_irn.adapt_epsF   = 1;
pars_irn.epsF_cutoff  = 0.05;   % This is the percentage cutoff

pars_irn.pcgtol_ini = 1e-4;

pars_irn.loops      = 5;
pars_irn.U0             = In;


pars_irn.pcgtol_ini     = 1e-1;
pars_irn.adaptPCGtol    = 0;

pars_irn.weight_scheme = NMP_WEIGHTS_THRESHOLD;

t = tic;
I_Threshold = irntv(In, {}, 0.075, pars_irn);
toc(t)

snr(I, I_Threshold)

% --------------------------
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
pars_irn.U0         = In;

pars_irn.variant       = NMP_TV_MIL;
pars_irn.weight_scheme = NMP_WEIGHTS_MIL;

pars_irn.pcgtol_ini     = 1e-1;
pars_irn.adaptPCGtol    = 0;

t = tic;
I_mil = irntv(In, {}, 0.075, pars_irn);
toc(t)

snr(I, I_mil)

%% Add HC affichagedes résultats

subplot(2,2,1)
imagesc(I); title('Original')
subplot(2,2,2)
imagesc(In); title('Noisy')
subplot(2,2,3)
imagesc(I_Threshold); title('Threshold')
subplot(2,2,4)
imagesc(I_mil); title('Matrix Inversion Lemma')


