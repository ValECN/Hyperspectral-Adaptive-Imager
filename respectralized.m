load('A.mat') % example input cube
IC = permute(IC2, [2 3 1]); % input cube in format (R,C,W)
 
IC_modified = respectralized_IC(IC,300,800);

figure
plot(squeeze(IC(10,10,:)))
hold on
plot(squeeze(IC_modified(10,10,:)))
hold on
plot(squeeze(IC(2,10,:)))
hold on 
plot(squeeze(IC_modified(2,10,:)))
legend('IC spectra 2','IC modified spectra 2','IC spectra 1','IC modified spectra 1')
title('Spectra comparison')

function IC_modified = respectralized_IC(IC, lba_min, lba_max) % lba : lambda

if isempty(lba_min)
    lba_min = 400; 
end
if isempty(lba_max)
    lba_max = 800; 
end

[R,C,W] = size(IC);
liste_spectra = [];


for x = 1:R
    for y = 1:C
        s = squeeze(IC(x,y,:));
        [~, index] = max(s);
        if isempty(liste_spectra)
            liste_spectra = index;
        else
            if ismember(index,liste_spectra) ~= 1
                liste_spectra(length(liste_spectra)+1) = index;  % storing all the different spectrums
            end
        end
    end
end

N = length(liste_spectra);
liste_spectra = sort(liste_spectra, 'ascend');                   % sorting in ascending order the spectrums
pts = linspace(lba_min, lba_max, N);                             % creating equally separed points in the range(lba_min:lba_max)
% interval = 0.5 * (pts(2) - pts(1));

IC_modified = zeros(size(IC));
for x = 1:R
    for y = 1:C
        s = squeeze(IC(x,y,:));
        [~, index] = max(s);
        for n = 1:N
            if index == liste_spectra(n)
                gap = pts(n)-liste_spectra(n);                   % quantify the difference between the real spectrum value and its attributed value
                abs_gap = abs(pts(n)-liste_spectra(n));
                for w = 1:W
                    if gap > 0
                        if mod(w+abs_gap,W) == 0                 % trick to skip the 0 value obtained for mod(x,x)
                            modulo = 1;
                        else 
                            modulo = mod(w+abs_gap,W);
                        end
                        IC_modified(x,y,modulo) = IC(x,y,w);
                    else
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