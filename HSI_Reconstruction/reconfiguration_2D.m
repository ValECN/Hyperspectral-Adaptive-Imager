function [liste_val, nb_pixels_in_spectra, I_modified] = reconfiguration_2D(I, val_min, val_max) 

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


factor = 4;

while length(liste_val) > R/4
    uniq = unique(round(liste_val,factor));
    D = diff(uniq);
    for k = 1:length(D)
        if D(k) <= 1^(-factor)
        uniq(k) = uniq(k+1);
        end
    end
    factor = factor-1;
    liste_val = uniq;
end

N = length(liste_val);
liste_val = sort(liste_val, 'ascend');
pts = linspace(val_min, val_max, N);
nb_pixels_in_spectra = zeros(size(liste_val));
I_modified = zeros(size(I));
    

for x = 1:R
    for y = 1:C
        s = round(I(x,y),factor+1);
        for n = 1:N
            if (s >= liste_val(n)-1) && (s <= liste_val(n)+1)
                I_modified(x,y) = pts(n);
                nb_pixels_in_spectra(n) = nb_pixels_in_spectra(n) + 1;
            end
        end
    end
end
end       