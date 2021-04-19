function [Gr, Gr_x, Gr_y] = contour_from_isolines(I,LOD) % LOD: Level Of Detail

Gr = zeros(size(I));
[R, ~] = size(I);
isolines = contour(I);
isolines(isolines>R) = -1; % To prevent incoherent values given by contour()

for k = 1:length(isolines)
    if isolines(1,k) >= 0 && isolines(2,k) >=0
        Gr(round(isolines(2,k)),round(isolines(1,k))) = 1;
    end
end

% Classical gradient calculation on the x and y axis
Gr_x = diff(I,1,2);
Gr_y = diff(I,1,1);

% Apply the LOD to keep the X% higher values
thresh_x = sort(unique(max(abs(Gr_x))),'descend');
crit_x = round(LOD * length(thresh_x));
Gr_x(abs(Gr_x) < thresh_x(crit_x)) = 0;

thresh_y = sort(unique(max(Gr_y)),'descend');
crit_y = round(LOD * length(thresh_y));
Gr_y(abs(Gr_y) < thresh_y(crit_y)) = 0;

% Binarize Gr_x and Gr_y
Gr_x = logical(Gr_x); % or 1 - logical(Gr_x)
Gr_y = logical(Gr_y); % or 1 - logical(Gr_y)

end