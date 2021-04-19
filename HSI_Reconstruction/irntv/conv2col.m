function y = conv2col(x,H,opt)
y1 = conv2(x(:,:,1),H,opt);
if size(x,3)>1,
    [N,M] = size(y1);
   y = zeros(N,N,size(x,3));
   y(:,:,1) = y1;
   for k = 2: size(x,3)
       y(:,:,k) = conv2(x(:,:,k),H,opt);
   end
else
    y = y1;
end
end

