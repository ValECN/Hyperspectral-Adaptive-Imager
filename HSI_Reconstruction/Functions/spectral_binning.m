%% --- Spectral binning function
function I_bin = spectral_binning(I, factor)

if length(size(I)) > 3
    I = I(:,:,:,1);
end

[R,C,W] = size(I);

K = floor(W/factor);

I_bin = zeros(R,C,K);

for r = 1:R
    for c = 1:C
        for k = 1:K
                I_bin(r,c,k) = sum((I(r,c,(factor*(k-1))+1 : factor*k)))/factor;
        end
    end
end
end