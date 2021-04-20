%% Reconstruction using Tikhonov regularization
% Original author: Valentin NOËL based on Dr.Ardi's work
% path_directory = '/home2/vnoel/Documents/Simulateur/simulator_playground/acquisition_data';

addpath Functions/
addpath Visualisation/
path_directory = [pwd, '/acquisition_data'];

% [FCS, I, panchro, DMD_conf, IC] = rebuild_aqucube_Val(path_directory);
[FCS, I, panchro, DMD_conf, IC] = rebuild_aqucube(path_directory);
panchro = panchro(:,:,1);
%%
[liste_val, nb_pixels_in_spectra, panchro_modified] = reconfiguration_2D(panchro, 100, 1000);  % reconfigured panchro
l = 5;
RA_def;

%% RA Ibrahim

[R, C, W, S] = size(FCS);

%% Init object

OC = FCS(:,:,:,1); % Output Cube
coef_NPS = param_OBJ.MAX_ADU / (max(max(sum(OC,3))) * param_OBJ.G * param_OBJ.rho);
Cube_PS = OC * coef_NPS;

% figure(1)
% histogram(I,[-.5 .5 1:100 101:50:4000])

%% *** RLS && CGNE: Parameters

if param_ED.run == 1
    [Gr, Gr_x, Gr_y] = contour_from_isolines(panchro_modified,0.92);
end

nGr_x = [Gr_x, zeros(R,1)];
nGr_y = [Gr_y; zeros(1,R)];

Gr_x = 1 - nGr_x;
Gr_y = 1 - nGr_y;

[L1, C1] = find(nGr_x);
[L2, C2] = find(nGr_y);

%%
Num_iter = zeros(length(param_CGNE.mu_x_vect),1);

%% The iterative method + direct inverse method

switch param_REC.method
    
    case 2
    
        H_C = sparse(R*C*S,R*C*W);  % compensated H
        H_NC = sparse(R*C*S,R*C*W);

        I_in = (1:C*W*R); 
        J_in = ceil(I_in/W);

        for s = 1:S

            FCS_C = permute(FCS(:,:,:,s),[3 2 1]); % to get WxCxR cube
            FCS_C = FCS_C(:);

            FCS_NC = permute(FCS(:,:,:,s),[3 2 1]); % to get WxCxR cube
            FCS_NC = FCS_NC(:);

            A_C  = (sparse(I_in,J_in,FCS_C)');
            A_NC = (sparse(I_in,J_in,FCS_NC)');

            H_C([(s-1)*R*C+1:s*R*C],:) = A_C;
            H_NC([(s-1)*R*C+1:s*R*C],:)= A_NC;  

        end

        [D_x, D_y, D_lambda] = derivative_operators(Gr_x,Gr_y,1,1,[R C W]);

        % *** Construction of CCD images for Direct Inversion
        if(S == 1)
            I_V = I';
            I_V = I_V(:);
        else
            I_V = permute(I,[2 1 3]);
            I_V = I_V(:);
        end

        % *** Construction of Gamma Matrix (variance & covariance Matrix) for Direct Inversion

        gamma = spdiags((1./I_V),0,R*C*S,R*C*S);
        Temp_vect = H_C'*gamma*I_V;

        for l = 1:length(param_CGNE.mu_x_vect)

        param_CGNE.mu_x = param_CGNE.mu_x_vect(l);
        param_CGNE.mu_y = param_CGNE.mu_x_vect(l);
        param_CGNE.mu_lambda = param_CGNE.mu_lambda_vect(l);
        
        tic
        [Cube_REC,epsilon_dx,epsilon_gradx,Num_iter(l)] = CGNE_Val(I,FCS,Gr_x,Gr_y, param_CGNE);
        toc
        
        I_vect = (H_C'*gamma*H_C+ param_CGNE.mu_x*(D_x)'*D_x+param_CGNE.mu_y*(D_y)'*D_y + param_CGNE.mu_lambda*(D_lambda')*D_lambda)\Temp_vect;
        I_vect = reshape(I_vect,W,C,R);
            
        Cube_REC_direct = permute(I_vect,[3 2 1]);

        param_CGNE.Error_l2(l)= error_HSI(Cube_PS, Cube_REC, 'l2');
        param_CGNE.Error_l1(l)= error_HSI(Cube_PS, Cube_REC, 'l1');
        
        figure
        imagesc(sum(Cube_REC,3))
        title(sprintf('panchro reconstructed (CGNE) for mu_x = %d and lambda = %d ', param_CGNE.mu_x_vect(l),param_CGNE.lambda))
        xlabel('X\_cam')
        ylabel('Y\_cam')
        
        figure
        imagesc(sum(Cube_REC_direct,3))
        title(sprintf('panchro reconstructed (inversion) for mu_x = %d and lambda = %d ', param_CGNE.mu_x_vect(l),param_CGNE.lambda))
        xlabel('X\_cam')
        ylabel('Y\_cam')

        end
    case 1
        
        for l = 1:length(param_CGNE.mu_x_vect)
            
            param_CGNE.mu_x = param_CGNE.mu_x_vect(l);
            param_CGNE.mu_y = param_CGNE.mu_y_vect(l);
            param_CGNE.mu_lambda = param_CGNE.mu_lambda_vect(l);

            tic
            [Cube_REC(:,:,:,l),epsilon_dx,epsilon_gradx,Num_iter(l)] = CGNE_Val(I,FCS,Gr_x,Gr_y, param_CGNE);
            toc

            param_CGNE.Error_l2(l)= error_HSI(Cube_PS, Cube_REC(:,:,:,l), 'l2');
            param_CGNE.Error_l1(l)= error_HSI(Cube_PS, Cube_REC(:,:,:,l), 'l1');
            
%             figure
%             plot(squeeze(Cube_REC(10,10,:,l)))
%             title(sprintf('Spectra of the reconstructed cube for x = y = 10 and mu_x = %d', param_CGNE.mu_x_vect(l)))
%             xlabel('Bandwidth')
%             ylabel('Amplitude')
            
        end
        % --- Plots
        figure
        imagesc(sum(Cube_REC(:,:,:,l),3))
        title(sprintf('panchro reconstructed for mu_x = %d and lambda = %d ', param_CGNE.mu_x_vect(l),param_CGNE.lambda))
        xlabel('X\_cam')
        ylabel('Y\_cam')
        
    case 0
        
        iter = 1;
        stop_iter = 10;
        mini = zeros(stop_iter,1);
        lbd = zeros(stop_iter,1);
        mini(1) = 1e9;
        step = 0.01;
        
        while iter < stop_iter
            
            N = length(param_CGNE.mu_x_vect);
            Cube_REC = zeros([R, C, W, N]);

            for l = 1:N

            param_CGNE.mu_x = param_CGNE.mu_x_vect(l);
            param_CGNE.mu_y = param_CGNE.mu_x_vect(l);
            param_CGNE.mu_lambda = param_CGNE.mu_lambda_vect(l);

            tic
            [Cube_REC(:,:,:,l),epsilon_dx,epsilon_gradx,Num_iter(l)] = CGNE_Val(I,FCS,Gr_x,Gr_y, param_CGNE);
            toc

            param_CGNE.Error_l2(l)= error_HSI(Cube_PS, Cube_REC(:,:,:,l), 'l2');
            param_CGNE.Error_l1(l)= error_HSI(Cube_PS, Cube_REC(:,:,:,l), 'l1');
            
            ecart = (sum(Cube_REC(:,:,3,l)) - sum(IC,3))^2;
            diff_lbd(l) = sum(sum(ecart));
            
            %**************************************************************************
            end
            if min(diff_lbd) < mini(iter)
                mini(iter+1) = min(diff_lbd);
                lbd(iter+1) = param_CGNE.lambda;
                param_CGNE.lambda = param_CGNE.lambda + step;
            else 
                mini(iter+1) = mini(iter);
                lbd(iter+1) = lbd(iter);
                param_CGNE.lambda = param_CGNE.lambda + step;
            end
            iter = iter +1;
        end

    otherwise
        disp('Enter another value for param_REC.method')
    
    figure
    imagesc(sum(Cube_REC,3))
    title(sprintf('panchro reconstructed for mu_x = %d and lambda = %d ', param_CGNE.mu_x_vect(1),param_CGNE.lambda))
    xlabel('X\_cam')
    ylabel('Y\_cam')
    
    figure
    for k = 2:S
        ecart = (sum(Cube_REC(:,:,3,k)) - sum(IC,3))^2;
        diff(k) = sum(sum(sum(ecart)));
        subplot(3,3,k-1)
        imagesc(sum(Cube_REC(:,:,:,k),3))
        title(sprintf('Panchro rec. for mu\_x = %d and lambda = %d', param_CGNE.mu_x_vect(k),param_CGNE.lambda))
        hold on
    end
    xlabel('X\_cam')
    ylabel('Y\_cam')
end
%% --- Spectral binning
% Cube_REC = Cube_REC(:,:,:,1); % binning to do on the FC
% Cube_REC_bin = spectral_binning(Cube_REC, 2);
% IC_bin = spectral_binning(IC, 8);

%% --- Plots
% l = 6;
% figure
% plot(squeeze(Cube_REC(10,10,:)))
% hold on
% plot(squeeze(IC(10,10,:)))
% title('Spectrally binned input and reconstructed cube bandwidth')
% xlabel('Bandwidth')
% ylabel('Amplitude')
% legend('Reconstructed spectra','Original spectra')   

figure
plot(squeeze(Cube_REC(2,10,:)))

figure 
imagesc(sum(I,3))
title('panchro input')
xlabel('X\_cam')
ylabel('Y\_cam')

figure
if length(size(Cube_REC)) > 3
    imagesc(sum(Cube_REC(:,:,:,l),3))
else
    imagesc(sum(Cube_REC,3))
end

hold on
plot([C1+ 0.5, C1+0.5]',[L1+0.5, L1-0.5]','r','LineWidth',2,'PickableParts','none');
plot([C2-0.5, C2+0.5]', [L2+0.5, L2+0.5]','r','LineWidth',2,'PickableParts','none');
hold off
title('Contours detected plotted on the reconstruted HSI')
xlabel('X\_cam')
ylabel('Y\_cam')

%% --- Visualization

% addpath Visualisation/
% 
% Object = Cube_REC;
% save('Object.mat','Object')
% 
% if length(size(Cube_REC)) > 3
%     Visual_HyperSpectral_2019(Cube_REC(:,:,:,l))
% else
%     Visual_HyperSpectral_2019(Cube_REC,3)
% end


% plot_spectra_neighbors(Cube_REC_bin,C1,C2,L1,L2)

