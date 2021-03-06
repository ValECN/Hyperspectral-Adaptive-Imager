%% *** Reconstruction using Tikhonov regularization
% *** Original author: Valentin NOËL based on Dr.Ardi's work

%% --- Initialization:

close all

%% --- Directories management:

addpath Functions/
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

%% --- Acquiring the wanted information from the simulated HSI:

[FCS, I, panchro, DMD_conf, IC] = rebuild_aqucube(path_directory);
% [I,FCS] = data_processing();
% FCS = spectral_binning(FCS, 2);
[R, C, W, S] = size(FCS);
RA_def; % file containing all the parameters

%% --- Re-configuration of the panchro for better contour detection:

[liste_val, nb_pixels_in_spectra, panchro_modified] = reconfiguration_2D(panchro(:,:,1), 100, 100000);  % reconfigured panchro

%% --- Contour detection:
param_ED.run = 1;
if param_ED.run == 1
    [Gr, Gr_x, Gr_y] = contour_from_isolines(panchro_modified,1);
    % -- Slight manipulations since contours are to be set to 0 for the CGNE:

    nGr_x = [Gr_x, zeros(R,1)];
    nGr_y = [Gr_y; zeros(1,R)];

    Gr_x = 1 - nGr_x;
    Gr_y = 1 - nGr_y;

    % -- Lx and Cx for following contour plot:

    [L1, C1] = find(nGr_x);
    [L2, C2] = find(nGr_y);
    
elseif param_ED.run == 2
    
    [nGr_x,nGr_y]=EdgeDetect_Function2(panchro_modified);
    % -- Slight manipulations since contours are to be set to 0 for the CGNE:

    Gr_x = 1 - nGr_x;
    Gr_y = 1 - nGr_y;

    % -- Lx and Cx for following contour plot:

    [L1, C1] = find(nGr_x);
    [L2, C2] = find(nGr_y);
end

%% --- The iterative method + direct inverse method

switch param_REC.method
    
    case 2 % Both methods
    
        H_C = sparse(R*C*S,R*C*W);  % compensated H
        H_NC = sparse(R*C*S,R*C*W); % non-compensated H

        I_in = (1:C*W*R); 
        J_in = ceil(I_in/W);

        for s = 1:S

            FCS_C = permute(FCS(:,:,:,s),[3 2 1]); % to get WxCxR cube
            FCS_C = FCS_C(:);

            FCS_NC = permute(FCS(:,:,:,s),[3 2 1]); % to get WxCxR cube
            FCS_NC = FCS_NC(:);

            A_C  = (sparse(I_in,J_in,FCS_C)');
            A_NC = (sparse(I_in,J_in,FCS_NC)');

            H_C((s-1)*R*C+1:s*R*C,:) = A_C;
            H_NC((s-1)*R*C+1:s*R*C,:)= A_NC;  

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

        for l = 1:length(param_CGNE.mu_x)
        
        tic
        [Cube_REC,epsilon_dx,epsilon_gradx,~] = CGNE_Val(I,FCS,Gr_x,Gr_y, param_CGNE);
        toc
        
        I_vect = (H_C'*gamma*H_C+ param_CGNE.mu_x*(D_x)'*D_x+param_CGNE.mu_y*(D_y)'*D_y + param_CGNE.mu_lambda*(D_lambda')*D_lambda)\Temp_vect;
        I_vect = reshape(I_vect,W,C,R);
            
        Cube_REC_direct = permute(I_vect,[3 2 1]);

        %param_CGNE.Error_l2(l)= error_HSI(Cube_PS, Cube_REC, 'l2');
        %param_CGNE.Error_l1(l)= error_HSI(Cube_PS, Cube_REC, 'l1');
    
        % -- Specific plot for the direct inversion reconstruction plot:
        figure
        imagesc(sum(Cube_REC_direct,3))
        title(sprintf('panchro reconstructed (inversion) for mu_x = %d and lambda = %d ', param_CGNE.mu_x,param_CGNE.mu_lambda))
        xlabel('X\_cam')
        ylabel('Y\_cam')

        end
        
    case 1 % Iterative case
        
        for l = 1:length(param_CGNE.mu_x)

            tic
            [Cube_REC,epsilon_dx,epsilon_gradx,~] = CGNE_Val(I,FCS,Gr_x,Gr_y, param_CGNE);
            toc

           % param_CGNE.Error_l2(l)= error_HSI(Cube_PS, Cube_REC, 'l2');
           % param_CGNE.Error_l1(l)= error_HSI(Cube_PS, Cube_REC, 'l1');
                      
        end

    otherwise    
        disp('Enter another value for param_REC.method')      
end

%% --- Spectral binning

IC_bin = spectral_binning(IC, 4);

%% --- Plots

figure(1)
imagesc(sum(Cube_REC,3))
hold on
plot([C1+ 0.5, C1+0.5]',[L1+0.5, L1-0.5]','r','LineWidth',2,'PickableParts','none');
plot([C2-0.5, C2+0.5]', [L2+0.5, L2+0.5]','r','LineWidth',2,'PickableParts','none');
hold off
title('Contours detected plotted on the reconstruted HSI')
xlabel('X\_cam')
ylabel('Y\_cam')

figure(2)
x1 = 8; x2 = 8; y1 = 13; y2 = 24;
plot(squeeze(Cube_REC(x1,y1,:)))
hold on 
plot(squeeze(Cube_REC(x2,y2,:)))
title('Spectra of the reconstructed cube using RA')
xlabel('Bandwidth')
ylabel('Amplitude')
legend(sprintf('Spectra Reconstruction for x = %d and y = %d',x1,y1),sprintf('Spectra Reconstruction for x = %d and y = %d',x2,y2))
grid on
ax = gca;
ax.GridColor = [0 .5 .5];
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';

%% --- Visualization 

if param_Visualization.run == 1
    % -- Saving the reconstructed cube as the .mat to use for Visualization
    Object = Cube_REC;
    save('Object.mat','Object')

    % -- Visualization
    Visual_HyperSpectral_2019(Cube_REC)
end

%% --- Metrics

IC_bin = spectral_binning(IC,1);
if size(IC_bin,3) > W
    IC_bin = IC_bin(:,:,1:W);
end
[ssim_val, ssim_map] = SSIM_map(IC_bin, Cube_REC);
[RMSE_map] = RMSE_map(IC_bin, Cube_REC);
[SAM_map] = SAM_map(IC_bin, Cube_REC);
RMSE = error_HSI(IC_bin, Cube_REC,'l2');