%% Computing the contours of the image using its isolines

I = contour_from_isolines(panchro);

% --- Plot
figure 
imagesc(I)
title('Contours of the input image')
xlabel('X\_cam')
ylabel('Y\_cam')

function Gr = contour_from_isolines(I)

Gr = zeros(size(I));
[R, ~] = size(I);
c = contour(I);
c(c>R) = -1;                                % To prevent incoherent values given by contour

for k = 1:length(c)
    if c(1,k) >= 0 && c(2,k) >=0
        Gr(round(c(2,k)),round(c(1,k))) = 1;
    end
end
end

    