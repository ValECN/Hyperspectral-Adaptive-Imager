param_REC.method = 1;
param_REC.data = 'K-poon';
                                        
param_ED.run = 1;
param_ED.Nit = 1;
param_ED.Meth = 6;
param_ED.Thresh = 0;
param_ED.Gap = 1;
param_ED.Mode = 't';

param_OBJ.G = 2;
param_OBJ.rho = .5;
param_OBJ.MAX_ADU = 3800;
param_OBJ.bits = 12;

param_CGNE.lambda =  1;
param_CGNE.mu_x = 30e-1;%15e-1 %logspace(-7,2,11);
param_CGNE.mu_y = param_CGNE.mu_x;
param_CGNE.mu_lambda = 3e-1; %12e-1; % param_CGNE.mu_x_vect * param_CGNE.lambda;
param_CGNE.Stopcriter = 'Eps';
param_CGNE.epsilonx   = 10^(-4);
param_CGNE.noise_type='i' ;                                                 % = 'i';'ii';'iis'
param_CGNE.Sig_noise = sqrt(I);                                             % std of the Poissoian noise
param_CGNE.x0 = 0;
param_CGNE.diff_l2 = zeros(length(param_CGNE.mu_x),1);
param_CGNE.Error_l2 = zeros(length(param_CGNE.mu_x),1);
param_CGNE.Error_l1 = zeros(length(param_CGNE.mu_x),1);
param_CGNE.stop_crit = 500;

param_Visualization.run = 1;