%% Computing the contours of the image using its isolines

[I, Gr_x, Gr_y] = contour_from_isolines(panchro);

% --- Plot
figure 
imagesc(I)
title('Contours of the input image')
xlabel('X\_cam')
ylabel('Y\_cam')

function [Gr, Gr_x, Gr_y] = contour_from_isolines(I)

Gr = zeros(size(I));
[R, ~] = size(I);
isolines = contour(I);
isolines(isolines>R) = -1; % To prevent incoherent values given by contour()

for k = 1:length(isolines)
    if isolines(1,k) >= 0 && isolines(2,k) >=0
        Gr(round(isolines(2,k)),round(isolines(1,k))) = 1;
    end
end

Gr_x = diff(I,1,2);
Gr_y = diff(I,1,1);

end

    