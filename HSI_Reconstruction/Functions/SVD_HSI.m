[R,C] = size(H_C);

% val = zeros(R*C,R*C*W);

val = svds(H_C,1e3);


% val = round(val,3);
% l = length(find(val(:,:,3)));
% size_2D = size(val(:,:,1));
% l_tot = size_2D(1) * size_2D(2);

% ratio_SVD = (l/l_tot);
