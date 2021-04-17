%% Re-spectralization and re-configuration of Hyperspectral Images and Panchromatic Images
% Original author: Valentin NO�L

% --- Initialization
load('A.mat');                    % loading a HSI of size (W,R,C)
IC      = permute(IC2, [2 3 1]);  % permuting to format (R,C,W)

IC_modified      = respectralized_IC(IC,200,800);          % respectralized HSI (moving Gaussian)

panchro = sum(IC_modified,3);     % panchromatic image
panchro_modified = reconfiguration_2D(panchro, 100, 600);  % reconfigured panchro

% --- Plots
figure(1)
plot(squeeze(IC(10,10,:)))
hold on
plot(squeeze(IC_modified(10,10,:)))
hold on
plot(squeeze(IC(2,10,:)))
hold on 
plot(squeeze(IC_modified(2,10,:)))
legend('IC spectra 2','IC modified spectra 2','IC spectra 1','IC modified spectra 1')
title('Spectra comparison')

figure(2)
subplot 211
imagesc(panchro)
title('Original panchromatic image')
colorbar

subplot 212
imagesc(panchro_modified)
title('Modified panchromatic image')
colorbar

% --- Function respectralized_IC
function IC_modified = respectralized_IC(IC, lba_min, lba_max)

[R,C,W] = size(IC);

if isempty(lba_min) || lba_min > lba_max || lba_min < 0 % generalizing configuration
    lba_min = 400; 
end
if isempty(lba_max) || lba_max < lba_min || lba_max > W
    lba_max = 800; 
end

liste_spectra = [];

for x = 1:R
    for y = 1:C
        s = squeeze(IC(x,y,:));                         % obtaining the spectrum (Gaussian shape)
        [~, index] = max(s);                            % finding its peak
        if isempty(liste_spectra)
            liste_spectra = index;                      % initializing the list with the first spectra
        else
            if ismember(index,liste_spectra) ~= 1
                liste_spectra(length(liste_spectra)+1) = index;
            end
        end
    end
end

N = length(liste_spectra);
liste_spectra = sort(liste_spectra, 'ascend');          % ascending order for easier attribution
pts = linspace(lba_min, lba_max, N);                    % equally separated points
% interval = 0.5 * (pts(2) - pts(1));

IC_modified = zeros(size(IC));
for x = 1:R
    for y = 1:C
        s = squeeze(IC(x,y,:));
        [~, index] = max(s);
        for n = 1:N
            if index == liste_spectra(n)
                gap = pts(n)-liste_spectra(n);
                abs_gap = abs(pts(n)-liste_spectra(n));
                for w = 1:W
                    if gap > 0
                        if mod(w+abs_gap,W) == 0        % trick to avoid mod(x,x) = 0
                            modulo = 1;
                        else 
                            modulo = mod(w+abs_gap,W);
                        end
                        IC_modified(x,y,modulo) = IC(x,y,w);
                    else                                % case where the spectrum is < or > than pts(n)
                        if mod(w-abs_gap,W) == 0
                            modulo = 1;
                        else 
                            modulo = mod(w-abs_gap,W);
                        end
                        IC_modified(x,y,modulo) = IC(x,y,w);
                    end
                end
            end
        end
    end
end

end

% --- Function reconfiguration_2D (same principle but in 2D)
function I_modified = reconfiguration_2D(I, val_min, val_max) 

[R, C] = size(I);
liste_val = [];
for x = 1:R
    for y = 1:C
        s = I(x,y);
         if isempty(liste_val)
            liste_val = s;
        else
            if ismember(s,liste_val) ~= 1
                liste_val(length(liste_val)+1) = s;
            end
         end
    end
end

N = length(liste_val);
liste_val = sort(liste_val, 'ascend');
pts = linspace(val_min, val_max, N);

I_modified = zeros(size(I));
for x = 1:R
    for y = 1:C
        s = I(x,y);
        for n = 1:N
            if s == liste_val(n)
                I_modified(x,y) = pts(n);
            end
        end
    end
end
end       







