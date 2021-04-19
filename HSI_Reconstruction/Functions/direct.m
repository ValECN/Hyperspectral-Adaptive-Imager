function Ao = direct(o,Filter)
% Entrée : o de taille (R,C,W)
%          Filter (W,S) cellules de taille (R,C) sparse
% Sortie : d de taille (R,C,S)
[R,C,W] = size(o);
S = size(Filter,2);

Ao = zeros(R,C,S);
for s=1:S
    Tmp = zeros(R,C);
    for w=1:W
        Tmp = Tmp + Filter{w,s}.*o(:,:,w);
    end
    Ao(:,:,s) = Tmp;
end

end

