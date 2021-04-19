function ATz = adj(z,Filter)
% Entrée : z de taille (R,C,S)
%          Filter (W,S) cellules de taille (R,C) sparse
% Sortie : ATz de taille (R,C,W)

[R,C,S] = size(z);
W = size(Filter,1);

ATz = zeros(R,C,W);
for w=1:W
    Tmp = zeros(R,C);
    for s = 1:S
        Tmp = Tmp + Filter{w,s}.*z(:,:,s);
    end
    ATz(:,:,w)= Tmp;
end
end


