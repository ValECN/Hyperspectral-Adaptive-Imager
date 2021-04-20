function RGB = hypercube_to_RGB(reflectances)

[R,C,L] = size(reflectances);

% To transform from reflectances to radiances
% load illum_6500.mat; illum = illum_6500;
% radiances = zeros(size(reflectances)); % initialize array
% for k = 1:L,
%   radiances(:,:,k) = reflectances(:,:,k)*illum(k);
% end
% 
% % To transform into RGB
% rad = reshape(radiances, R*C,L);
% load xyzbar.mat;
% XYZ = (xyzbar(1:L,:)'*rad')';
% XYZ = reshape(XYZ, R, C, 3);
% XYZ = max(XYZ, 0);
% XYZ = XYZ/max(XYZ(:));
% 
% % Image RGB
% 
% RGB = XYZ2sRGB_exgamma(XYZ);
% RGB = max(RGB, 0);
% RGB = min(RGB, 1);

RGB = sum(reflectances,3);

% R = reflectances(:,:,20);
% G = reflectances(:,:,14);
% B = reflectances(:,:,5);
% 
% RGB = zeros(size(reflectances,1), size(reflectances,2), 3);
% RGB(:,:,1) = R;
% RGB(:,:,2) = G;
% RGB(:,:,3) = B;
% 
