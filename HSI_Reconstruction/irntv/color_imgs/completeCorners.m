function f = completeCorners(f,it)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   completeCorners function will fill the gaps (in pixels)
%                   to create a corner when needed to prevent leaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 1;
[R,C] = size(f);
for k = 1:it
    for i = (a+1):(R-a)
        for j = (a+1):(C-a)
            if ((f(i,j-a) ~= 0) && (f(i+a,j) ~= 0) && (f(i-a,j) == 0) && (f(i-a,j+a) == 0) && (f(i,j+a) == 0) && (f(i+a,j-a) == 0)) || (((f(i-a,j) ~= 0) && (f(i,j+a) ~= 0)) && (f(i,j-a) == 0) && (f(i+a,j-a) == 0) && (f(i+a,j) == 0) && (f(i-a,j+a) == 0)) || (((f(i,j-a) ~= 0) && (f(i-a,j) ~= 0)) && (f(i,j+a) == 0) && (f(i+a,j+a) == 0) && (f(i+a,j) == 0) && (f(i-a,j-a) == 0)) || (((f(i,j+a) ~= 0) && (f(i+a,j) ~= 0)) && (f(i,j-a) == 0) && (f(i-a,j-a) == 0) && (f(i-a,j) == 0) && (f(i+a,j+a) == 0))
                f(i,j) = 255;
            end
        end
    end
end

