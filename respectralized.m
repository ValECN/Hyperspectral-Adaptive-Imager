load('A.mat')
IC = permute(IC2, [2 3 1]);
 
IC_modified = respectralized_IC(IC,200,800);

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

function IC_modified = respectralized_IC(IC, lba_min, lba_max)
[R,C,W] = size(IC);

if isempty(lba_min) || lba_min > lba_max || lba_min < 0
    lba_min = 400; 
end
if isempty(lba_max) || lba_max < lba_min || lba_max > W
    lba_max = 800; 
end

liste_spectra = [];

for x = 1:R
    for y = 1:C
        s = squeeze(IC(x,y,:));
        [~, index] = max(s);
        if isempty(liste_spectra)
            liste_spectra = index;
        else
            if ismember(index,liste_spectra) ~= 1
                liste_spectra(length(liste_spectra)+1) = index;
            end
        end
    end
end

N = length(liste_spectra);
liste_spectra = sort(liste_spectra, 'ascend');
pts = linspace(lba_min, lba_max, N);
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
                        if mod(w+abs_gap,W) == 0
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